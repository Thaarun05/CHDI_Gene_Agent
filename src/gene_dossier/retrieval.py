"""Evidence retrieval: keyword/metadata first, optional Chroma semantic second.

Chroma is an **index only** — never the source of truth. Structured truth stays
in the provenance DB (Postgres/SQLite); raw material stays on disk via
``raw_store``.

MVP scope (no hybrid RAG / reranking yet):

1. In-memory keyword + metadata filters over :class:`EvidenceRecord`
2. Optional Chroma collection keyed by ``{dossier_run_id}:{evidence_record_id}``
   with evidence text plus bounded structured context, and filterable metadata
   (gene, run, section, source, grade, assertion_type)

Missing chromadb, embedding model downloads, or disk errors soft-fail so the
dossier pipeline still runs without an LLM or vector index.
"""

from __future__ import annotations

import hashlib
import logging
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import EvidenceGrade, EvidenceRecord

logger = logging.getLogger(__name__)

DEFAULT_COLLECTION = "evidence_records"
_TOKEN_RE = re.compile(r"[a-z0-9_]+", re.IGNORECASE)
INDEX_SCHEMA_VERSION = 1


# --------------------------------------------------------------------------------------
# Keyword / metadata retrieval (always available)
# --------------------------------------------------------------------------------------
def _norm(value: str | None) -> str:
    return " ".join((value or "").strip().lower().split())


def _enum_value(value: Any) -> str:
    if value is None:
        return ""
    return str(getattr(value, "value", value))


def vector_id_for_record(record: EvidenceRecord) -> str:
    """Return the Chroma primary key for one evidence record.

    ``source_id`` can repeat across historical dossier runs, so vector IDs must
    include both the provenance run and the persisted EvidenceRecord primary key.
    """
    return f"{record.dossier_run_id}:{record.id}"


def tokenize(text: str) -> list[str]:
    """Lowercase alphanumeric tokens for keyword matching."""
    return [t.lower() for t in _TOKEN_RE.findall(text or "")]


def evidence_matches_filters(
    record: EvidenceRecord,
    *,
    gene_symbol: str | None = None,
    section: str | None = None,
    source_name: str | None = None,
    evidence_grade: EvidenceGrade | str | None = None,
    assertion_type: str | None = None,
    source_id: str | None = None,
) -> bool:
    """Return True when ``record`` satisfies all provided metadata filters.

    Keyword-path ``section`` matching is substring-based (normalized contains).
    Other fields are exact matches after normalization. Chroma ``where`` filters
    used by :class:`ChromaEvidenceIndex` are exact-match only.
    """
    if gene_symbol and _norm(record.gene_symbol) != _norm(gene_symbol):
        return False
    if section and _norm(section) not in _norm(record.section):
        return False
    if source_name and _norm(record.source_name) != _norm(source_name):
        return False
    if evidence_grade is not None:
        want = _enum_value(evidence_grade).upper()
        have = _enum_value(record.evidence_grade).upper()
        if want and have != want:
            return False
    if assertion_type is not None:
        want = _enum_value(assertion_type).lower()
        have = _enum_value(record.assertion_type).lower()
        if want and have != want:
            return False
    if source_id and record.source_id != source_id:
        return False
    return True


def keyword_score(query: str, record: EvidenceRecord) -> float:
    """Simple token-overlap score against display_text / section / source / ids."""
    q_tokens = set(tokenize(query))
    if not q_tokens:
        return 0.0
    haystack = " ".join(
        [
            record.display_text or "",
            record.section or "",
            record.subsection or "",
            record.source_name or "",
            record.fact_type or "",
            record.source_id or "",
            _enum_value(record.assertion_type),
        ]
    )
    doc_tokens = set(tokenize(haystack))
    if not doc_tokens:
        return 0.0
    overlap = q_tokens & doc_tokens
    if not overlap:
        return 0.0
    return len(overlap) / float(len(q_tokens))


@dataclass(frozen=True)
class RetrievalHit:
    """One ranked evidence hit from keyword and/or semantic retrieval."""

    record: EvidenceRecord
    score: float
    method: str  # "keyword" | "semantic" | "hybrid"
    source_id: str


@dataclass
class KeywordEvidenceIndex:
    """In-memory keyword + metadata index over EvidenceRecords."""

    records: list[EvidenceRecord] = field(default_factory=list)

    def add(self, records: Iterable[EvidenceRecord]) -> int:
        batch = list(records)
        self.records.extend(batch)
        return len(batch)

    def clear(self) -> None:
        self.records.clear()

    def search(
        self,
        query: str = "",
        *,
        gene_symbol: str | None = None,
        section: str | None = None,
        source_name: str | None = None,
        evidence_grade: EvidenceGrade | str | None = None,
        assertion_type: str | None = None,
        source_id: str | None = None,
        limit: int = 20,
    ) -> list[RetrievalHit]:
        """Filter by metadata, then rank by keyword overlap (if query given)."""
        filtered = [
            r
            for r in self.records
            if evidence_matches_filters(
                r,
                gene_symbol=gene_symbol,
                section=section,
                source_name=source_name,
                evidence_grade=evidence_grade,
                assertion_type=assertion_type,
                source_id=source_id,
            )
        ]
        hits: list[RetrievalHit] = []
        q = (query or "").strip()
        for record in filtered:
            score = keyword_score(q, record) if q else 1.0
            if q and score <= 0.0:
                continue
            hits.append(
                RetrievalHit(
                    record=record,
                    score=score,
                    method="keyword",
                    source_id=record.source_id,
                )
            )
        hits.sort(key=lambda h: (-h.score, h.source_id))
        return hits[: max(0, limit)]


def search_evidence_keyword(
    records: Iterable[EvidenceRecord],
    query: str = "",
    *,
    gene_symbol: str | None = None,
    section: str | None = None,
    source_name: str | None = None,
    evidence_grade: EvidenceGrade | str | None = None,
    assertion_type: str | None = None,
    source_id: str | None = None,
    limit: int = 20,
) -> list[RetrievalHit]:
    """Convenience wrapper: build a one-shot keyword index and search it."""
    index = KeywordEvidenceIndex()
    index.add(records)
    return index.search(
        query,
        gene_symbol=gene_symbol,
        section=section,
        source_name=source_name,
        evidence_grade=evidence_grade,
        assertion_type=assertion_type,
        source_id=source_id,
        limit=limit,
    )


# --------------------------------------------------------------------------------------
# Chroma semantic skeleton (optional)
# --------------------------------------------------------------------------------------
class HashEmbeddingFunction:
    """Deterministic local embedding — no model download / no network.

    Used as a safe default so Chroma indexing works offline in CI and when the
    ONNX MiniLM model cannot be downloaded. Semantic quality is limited; swap
    in a real embedding function when configured.
    """

    def __init__(self, dimensions: int = 64) -> None:
        self.dimensions = dimensions

    def name(self) -> str:
        return "gene_dossier_hash_embedding"

    def __call__(self, input: list[str]) -> list[list[float]]:  # noqa: A003
        return [self._embed_one(text) for text in input]

    def embed_query(self, input: list[str]) -> list[list[float]]:  # noqa: A003
        """Chroma 1.x query path; mirrors ``__call__`` for this hash embedder."""
        return self(input)

    def _embed_one(self, text: str) -> list[float]:
        digest = hashlib.sha256((text or "").encode("utf-8")).digest()
        values: list[float] = []
        seed = digest
        while len(values) < self.dimensions:
            seed = hashlib.sha256(seed).digest()
            for byte in seed:
                if len(values) >= self.dimensions:
                    break
                values.append((byte / 127.5) - 1.0)
        norm = math.sqrt(sum(v * v for v in values)) or 1.0
        return [v / norm for v in values]


class LangChainEmbeddingFunction:
    """Chroma-compatible wrapper around a real LangChain embedding provider."""

    _gene_dossier_embedding_backend = "real"

    def __init__(self, embeddings: Any, *, provider: str, model: str) -> None:
        self.embeddings = embeddings
        self.provider = provider
        self.model = model

    def name(self) -> str:
        return f"gene_dossier_{self.provider}_{self.model}".replace("/", "_").replace(":", "_")

    def __call__(self, input: list[str]) -> list[list[float]]:  # noqa: A003
        return self.embeddings.embed_documents(list(input))

    def embed_query(self, input: list[str]) -> list[list[float]]:  # noqa: A003
        return self(input)


def build_local_minilm_embedding_function() -> Any | None:
    """Build Chroma's local ONNX all-MiniLM-L6-v2 embedding function.

    Chroma downloads public model weights on first use, then performs document
    and query embedding locally. Evidence text is never sent to a provider.
    """
    try:
        from chromadb.utils.embedding_functions import DefaultEmbeddingFunction

        embedding_function = DefaultEmbeddingFunction()
        embedding_function._gene_dossier_embedding_backend = "local_minilm"  # type: ignore[attr-defined]
        embedding_function._gene_dossier_embedding_model = "all-MiniLM-L6-v2"  # type: ignore[attr-defined]
        return embedding_function
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to build local MiniLM embedding model: %s", exc)
        return None


def build_real_embedding_function(settings: Settings | None = None) -> Any | None:
    """Build an optional external semantic provider for non-demo callers."""
    cfg = settings or get_settings()
    if cfg.has_key("openai_api_key"):
        try:
            from langchain_openai import OpenAIEmbeddings

            model_name = "text-embedding-3-small"
            kwargs: dict[str, Any] = {
                "model": model_name,
                "api_key": cfg.openai_api_key,
                "tiktoken_enabled": False,
                "check_embedding_ctx_length": False,
            }
            base_url = (cfg.openai_base_url or "").strip()
            if base_url:
                kwargs["base_url"] = base_url
            return LangChainEmbeddingFunction(
                OpenAIEmbeddings(**kwargs),
                provider="openai",
                model=model_name,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to build OpenAI embedding model: %s", exc)
    return None


def _record_metadata(record: EvidenceRecord) -> dict[str, str]:
    """Chroma metadata values must be scalar (str/int/float/bool)."""
    return {
        "source_id": record.source_id or "",
        "evidence_record_id": record.id or "",
        "gene_symbol": record.gene_symbol or "",
        "section": record.section or "",
        "subsection": record.subsection or "",
        "source_name": record.source_name or "",
        "evidence_grade": _enum_value(record.evidence_grade),
        "assertion_type": _enum_value(record.assertion_type),
        "fact_type": record.fact_type or "",
        "dossier_run_id": record.dossier_run_id or "",
    }


def _record_document(record: EvidenceRecord) -> str:
    """Bound semantic text to the supplied EvidenceRecord fields only."""
    parts = [
        record.display_text or "",
        f"gene: {record.gene_symbol or ''}",
        f"section: {record.section or ''}",
        f"subsection: {record.subsection or ''}",
        f"source: {record.source_name or ''}",
        f"fact_type: {record.fact_type or ''}",
        f"assertion_type: {_enum_value(record.assertion_type)}",
        f"source_id: {record.source_id or ''}",
        f"evidence_record_id: {record.id or ''}",
        f"dossier_run_id: {record.dossier_run_id or ''}",
    ]
    return "\n".join(part for part in parts if part.strip())


def _where_from_filters(
    *,
    gene_symbol: str | None = None,
    section: str | None = None,
    source_name: str | None = None,
    evidence_grade: EvidenceGrade | str | None = None,
    assertion_type: str | None = None,
    source_id: str | None = None,
    dossier_run_ids: list[str] | tuple[str, ...] | set[str] | None = None,
) -> dict[str, Any] | None:
    """Build a Chroma ``where`` clause (exact-match only; no substring section)."""
    clauses: list[dict[str, Any]] = []
    if gene_symbol:
        clauses.append({"gene_symbol": gene_symbol})
    if section:
        clauses.append({"section": section})
    if source_name:
        clauses.append({"source_name": source_name})
    if evidence_grade is not None:
        clauses.append({"evidence_grade": _enum_value(evidence_grade)})
    if assertion_type is not None:
        clauses.append({"assertion_type": _enum_value(assertion_type)})
    if source_id:
        clauses.append({"source_id": source_id})
    run_ids = [rid for rid in dict.fromkeys(dossier_run_ids or []) if rid]
    if len(run_ids) == 1:
        clauses.append({"dossier_run_id": run_ids[0]})
    elif len(run_ids) > 1:
        clauses.append({"$or": [{"dossier_run_id": rid} for rid in run_ids]})
    if not clauses:
        return None
    if len(clauses) == 1:
        return clauses[0]
    return {"$and": clauses}


@dataclass
class ChromaIndexStatus:
    """Soft status for Chroma availability / last operation."""

    available: bool
    backend: str  # "persistent" | "ephemeral" | "unavailable"
    collection: str | None = None
    path: str | None = None
    error: str | None = None
    indexed_count: int = 0
    embedding_backend: str = "unavailable"  # local_minilm | real | hash_test_fallback | unavailable
    embedding_model: str = "unknown"
    embedding_dimension: int | None = None
    index_identity: str | None = None


@dataclass(frozen=True)
class EmbeddingIndexIdentity:
    """Stable identity for a Chroma collection's embedding configuration."""

    provider: str
    model: str
    dimension: int
    schema_version: int = INDEX_SCHEMA_VERSION

    @property
    def key(self) -> str:
        return (
            f"schema{self.schema_version}"
            f"__{self.provider}"
            f"__{self.model}"
            f"__dim{self.dimension}"
        )

    def metadata(self) -> dict[str, str | int]:
        return {
            "gene_dossier_index_schema": self.schema_version,
            "gene_dossier_embedding_provider": self.provider,
            "gene_dossier_embedding_model": self.model,
            "gene_dossier_embedding_dimension": self.dimension,
            "gene_dossier_index_identity": self.key,
        }


def _safe_collection_component(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    text = re.sub(r"_+", "_", text).strip("._-")
    return text[:40] or "unknown"


def _dimension_of_vector(vector: Any) -> int | None:
    try:
        return len(vector)
    except TypeError:
        return None


def _infer_embedding_dimension(embedding_function: Any) -> int | None:
    explicit = getattr(embedding_function, "dimensions", None)
    if isinstance(explicit, int) and explicit > 0:
        return explicit
    try:
        sample = embedding_function(["gene dossier embedding dimension probe"])
    except Exception as exc:  # noqa: BLE001
        logger.warning("Embedding dimension probe failed: %s", exc)
        return None
    if not sample:
        return None
    return _dimension_of_vector(sample[0])


def _embedding_identity(
    embedding_function: Any,
    *,
    embedding_backend: str,
) -> EmbeddingIndexIdentity | None:
    dimension = _infer_embedding_dimension(embedding_function)
    if not dimension:
        return None
    provider = (
        getattr(embedding_function, "provider", None)
        or getattr(embedding_function, "_gene_dossier_embedding_provider", None)
        or embedding_backend
        or "unknown"
    )
    model = (
        getattr(embedding_function, "model", None)
        or getattr(embedding_function, "_gene_dossier_embedding_model", None)
        or getattr(embedding_function, "name", lambda: "unknown")()
    )
    return EmbeddingIndexIdentity(
        provider=_safe_collection_component(str(provider)),
        model=_safe_collection_component(str(model)),
        dimension=dimension,
    )


def collection_name_for_embedding(
    base_collection_name: str,
    identity: EmbeddingIndexIdentity,
) -> str:
    """Return a deterministic Chroma collection name scoped to embedding identity."""
    base = _safe_collection_component(base_collection_name)
    digest = hashlib.sha1(identity.key.encode("utf-8")).hexdigest()[:12]
    suffix = f"{identity.provider}_{identity.model}_d{identity.dimension}_{digest}"
    name = f"{base}__{suffix}"
    # Chroma collection names must be 3-63 characters.
    return name[:63].rstrip("._-")


class ChromaEvidenceIndex:
    """Optional Chroma collection over EvidenceRecord display_text + metadata.

    Soft-fails construction and operations when chromadb is missing or broken.
    """

    def __init__(
        self,
        *,
        settings: Settings | None = None,
        persist_directory: str | Path | None = None,
        collection_name: str = DEFAULT_COLLECTION,
        embedding_function: Any | None = None,
        ephemeral: bool = False,
        allow_hash_fallback: bool = True,
        allow_external_embedding_provider: bool = True,
        read_only: bool = False,
        read_only_client_factory: Callable[[Path], Any] | None = None,
    ) -> None:
        self.settings = settings or get_settings()
        self.requested_collection_name = collection_name
        self.collection_name = collection_name
        self.embedding_backend = "unavailable"
        self.read_only = read_only
        if embedding_function is not None:
            self.embedding_function = embedding_function
            self.embedding_backend = getattr(
                embedding_function,
                "_gene_dossier_embedding_backend",
                "real",
            )
        elif allow_external_embedding_provider:
            real_embedding = build_real_embedding_function(self.settings)
            if real_embedding is not None:
                self.embedding_function = real_embedding
                self.embedding_backend = "real"
            elif allow_hash_fallback:
                self.embedding_function = HashEmbeddingFunction()
                self.embedding_function._gene_dossier_embedding_backend = "hash_test_fallback"  # type: ignore[attr-defined]
                self.embedding_backend = "hash_test_fallback"
            else:
                self.embedding_function = None
        elif allow_hash_fallback:
            self.embedding_function = HashEmbeddingFunction()
            self.embedding_function._gene_dossier_embedding_backend = "hash_test_fallback"  # type: ignore[attr-defined]
            self.embedding_backend = "hash_test_fallback"
        else:
            self.embedding_function = None
        self._client: Any | None = None
        self._collection: Any | None = None
        self.status = ChromaIndexStatus(available=False, backend="unavailable")

        if self.embedding_function is None:
            self.status = ChromaIndexStatus(
                available=False,
                backend="unavailable",
                error="semantic embedding function unavailable",
                embedding_backend="unavailable",
            )
            return

        self.index_identity = _embedding_identity(
            self.embedding_function,
            embedding_backend=self.embedding_backend,
        )
        if self.index_identity is None:
            self.status = ChromaIndexStatus(
                available=False,
                backend="unavailable",
                error="embedding dimension could not be determined",
                embedding_backend=self.embedding_backend,
            )
            return
        self.collection_name = collection_name_for_embedding(
            self.requested_collection_name,
            self.index_identity,
        )

        try:
            import chromadb
        except Exception as exc:  # noqa: BLE001
            self.status = ChromaIndexStatus(
                available=False,
                backend="unavailable",
                error=f"chromadb import failed: {exc}",
                embedding_backend=self.embedding_backend,
            )
            logger.warning("Chroma unavailable: %s", exc)
            return

        try:
            if ephemeral:
                if read_only:
                    raise ValueError("read-only Chroma requires a persistent collection")
                self._client = chromadb.EphemeralClient()
                backend = "ephemeral"
                path_str = None
            else:
                path = Path(
                    persist_directory
                    if persist_directory is not None
                        else self.settings.index_path
                )
                if read_only and not path.is_dir():
                    raise FileNotFoundError("read-only Chroma directory does not exist")
                if not read_only:
                    path.mkdir(parents=True, exist_ok=True)
                if read_only:
                    if read_only_client_factory is None:
                        raise RuntimeError(
                            "installed Chroma PersistentClient cannot guarantee "
                            "non-mutating access; read-only semantic retrieval disabled"
                        )
                    self._client = read_only_client_factory(path)
                else:
                    self._client = chromadb.PersistentClient(path=str(path))
                backend = "persistent"
                path_str = str(path)
            if read_only:
                self._collection = self._client.get_collection(
                    name=self.collection_name,
                    embedding_function=self.embedding_function,
                )
            else:
                self._collection = self._client.get_or_create_collection(
                    name=self.collection_name,
                    embedding_function=self.embedding_function,
                    metadata={
                        "hnsw:space": "cosine",
                        **self.index_identity.metadata(),
                    },
                )
            metadata = getattr(self._collection, "metadata", None) or {}
            metadata_identity = metadata.get("gene_dossier_index_identity")
            if metadata_identity and metadata_identity != self.index_identity.key:
                raise ValueError(
                    "Chroma collection embedding identity mismatch: "
                    f"{metadata_identity!r} != {self.index_identity.key!r}"
                )
            self.status = ChromaIndexStatus(
                available=True,
                backend=backend,
                collection=self.collection_name,
                path=path_str,
                indexed_count=int(self._collection.count()),
                embedding_backend=self.embedding_backend,
                embedding_model=self.index_identity.model,
                embedding_dimension=self.index_identity.dimension,
                index_identity=self.index_identity.key,
            )
        except Exception as exc:  # noqa: BLE001
            self.status = ChromaIndexStatus(
                available=False,
                backend="unavailable",
                error=str(exc),
                embedding_backend=self.embedding_backend,
                embedding_model=getattr(self, "index_identity", None).model
                if getattr(self, "index_identity", None)
                else "unknown",
                embedding_dimension=getattr(self, "index_identity", None).dimension
                if getattr(self, "index_identity", None)
                else None,
                index_identity=getattr(self, "index_identity", None).key
                if getattr(self, "index_identity", None)
                else None,
            )
            logger.warning("Chroma init failed; semantic search disabled: %s", exc)
            self._client = None
            self._collection = None

    @property
    def available(self) -> bool:
        return bool(self.status.available and self._collection is not None)

    def upsert_evidence(self, records: Iterable[EvidenceRecord]) -> int:
        """Index / update records by run-qualified EvidenceRecord ID. Returns count upserted."""
        if self.read_only:
            return 0
        if not self.available or self._collection is None:
            return 0
        ids: list[str] = []
        documents: list[str] = []
        metadatas: list[dict[str, str]] = []
        for record in records:
            sid = vector_id_for_record(record)
            text = _record_document(record).strip()
            if not text:
                continue
            ids.append(sid)
            documents.append(text)
            metadatas.append(_record_metadata(record))
        if not ids:
            return 0
        try:
            metadata = getattr(self._collection, "metadata", None) or {}
            expected = metadata.get("gene_dossier_embedding_dimension")
            actual = getattr(getattr(self, "index_identity", None), "dimension", None)
            if expected is not None and actual is not None and int(expected) != int(actual):
                raise ValueError(
                    "Chroma collection embedding dimension mismatch before upsert: "
                    f"{expected} != {actual}"
                )
            self._collection.upsert(
                ids=ids, documents=documents, metadatas=metadatas
            )
            self.status.indexed_count = int(self._collection.count())
            return len(ids)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Chroma upsert failed: %s", exc)
            self.status.error = str(exc)
            return 0

    def query(
        self,
        query: str,
        *,
        gene_symbol: str | None = None,
        section: str | None = None,
        source_name: str | None = None,
        evidence_grade: EvidenceGrade | str | None = None,
        assertion_type: str | None = None,
        source_id: str | None = None,
        dossier_run_ids: list[str] | tuple[str, ...] | set[str] | None = None,
        limit: int = 10,
        record_lookup: dict[str, EvidenceRecord] | None = None,
    ) -> list[RetrievalHit]:
        """Semantic query. Returns hits (with records when ``record_lookup`` given)."""
        if not self.available or self._collection is None:
            return []
        q = (query or "").strip()
        if not q:
            return []
        where = _where_from_filters(
            gene_symbol=gene_symbol,
            section=section,
            source_name=source_name,
            evidence_grade=evidence_grade,
            assertion_type=assertion_type,
            source_id=source_id,
            dossier_run_ids=dossier_run_ids,
        )
        try:
            kwargs: dict[str, Any] = {
                "query_texts": [q],
                "n_results": max(1, limit),
                "include": ["documents", "metadatas", "distances"],
            }
            if where is not None:
                kwargs["where"] = where
            raw = self._collection.query(**kwargs)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Chroma query failed: %s", exc)
            self.status.error = str(exc)
            return []

        ids = (raw.get("ids") or [[]])[0]
        distances = (raw.get("distances") or [[]])[0]
        hits: list[RetrievalHit] = []
        lookup = record_lookup or {}
        for i, sid in enumerate(ids):
            distance = distances[i] if i < len(distances) else None
            # Cosine distances can exceed 1.0; use 1/(1+d) so ranking is preserved
            # instead of clamping many hits to score 0.0.
            if distance is None:
                score = 0.0
            else:
                distance_f = max(0.0, float(distance))
                score = 1.0 / (1.0 + distance_f)
            record = lookup.get(sid)
            if record is None:
                continue
            hits.append(
                RetrievalHit(
                    record=record,
                    score=score,
                    method="semantic",
                    source_id=sid,
                )
            )
        return hits


def index_evidence_in_chroma(
    records: Iterable[EvidenceRecord],
    *,
    settings: Settings | None = None,
    persist_directory: str | Path | None = None,
    collection_name: str = DEFAULT_COLLECTION,
    ephemeral: bool = False,
) -> ChromaIndexStatus:
    """Upsert evidence into Chroma; soft-fails to ``available=False`` status."""
    index = ChromaEvidenceIndex(
        settings=settings,
        persist_directory=persist_directory,
        collection_name=collection_name,
        ephemeral=ephemeral,
    )
    count = index.upsert_evidence(records)
    index.status.indexed_count = count if index.available else 0
    return index.status


def search_evidence(
    records: Iterable[EvidenceRecord],
    query: str,
    *,
    gene_symbol: str | None = None,
    section: str | None = None,
    source_name: str | None = None,
    evidence_grade: EvidenceGrade | str | None = None,
    assertion_type: str | None = None,
    source_id: str | None = None,
    dossier_run_ids: list[str] | tuple[str, ...] | set[str] | None = None,
    limit: int = 20,
    use_chroma: bool = True,
    chroma_index: ChromaEvidenceIndex | None = None,
    settings: Settings | None = None,
) -> list[RetrievalHit]:
    """Keyword/metadata first; optionally merge Chroma semantic hits.

    Dedupes by ``source_id``, keeping the higher score. Semantic hits are only
    attempted when ``use_chroma`` is True and an index is available.
    """
    record_list = list(records)
    keyword_hits = search_evidence_keyword(
        record_list,
        query,
        gene_symbol=gene_symbol,
        section=section,
        source_name=source_name,
        evidence_grade=evidence_grade,
        assertion_type=assertion_type,
        source_id=source_id,
        limit=limit,
    )
    by_id: dict[str, RetrievalHit] = {
        vector_id_for_record(h.record): RetrievalHit(
            record=h.record,
            score=h.score,
            method=h.method,
            source_id=vector_id_for_record(h.record),
        )
        for h in keyword_hits
    }

    if use_chroma and (query or "").strip():
        index = chroma_index
        if index is None:
            index = ChromaEvidenceIndex(
                settings=settings or get_settings(),
                ephemeral=True,
            )
            if index.available:
                index.upsert_evidence(record_list)
        if index is not None and index.available:
            lookup = {vector_id_for_record(r): r for r in record_list}
            semantic_hits = index.query(
                query,
                gene_symbol=gene_symbol,
                section=section,
                source_name=source_name,
                evidence_grade=evidence_grade,
                assertion_type=assertion_type,
                source_id=source_id,
                dossier_run_ids=dossier_run_ids,
                limit=limit,
                record_lookup=lookup,
            )
            for hit in semantic_hits:
                existing = by_id.get(hit.source_id)
                if existing is None:
                    by_id[hit.source_id] = hit
                elif hit.score > existing.score:
                    by_id[hit.source_id] = RetrievalHit(
                        record=hit.record,
                        score=hit.score,
                        method="hybrid",
                        source_id=hit.source_id,
                    )
                elif existing.method == "keyword":
                    by_id[hit.source_id] = RetrievalHit(
                        record=existing.record,
                        score=existing.score,
                        method="hybrid",
                        source_id=existing.source_id,
                    )

    merged = sorted(by_id.values(), key=lambda h: (-h.score, h.source_id))
    return merged[: max(0, limit)]


__all__ = [
    "DEFAULT_COLLECTION",
    "RetrievalHit",
    "vector_id_for_record",
    "KeywordEvidenceIndex",
    "HashEmbeddingFunction",
    "build_local_minilm_embedding_function",
    "EmbeddingIndexIdentity",
    "collection_name_for_embedding",
    "ChromaIndexStatus",
    "ChromaEvidenceIndex",
    "tokenize",
    "evidence_matches_filters",
    "keyword_score",
    "search_evidence_keyword",
    "index_evidence_in_chroma",
    "search_evidence",
]
