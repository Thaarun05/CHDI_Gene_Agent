"""Evidence retrieval: keyword/metadata first, optional Chroma semantic second.

Chroma is an **index only** — never the source of truth. Structured truth stays
in the provenance DB (Postgres/SQLite); raw material stays on disk via
``raw_store``.

MVP scope (no hybrid RAG / reranking yet):

1. In-memory keyword + metadata filters over :class:`EvidenceRecord`
2. Optional Chroma collection keyed by ``source_id`` with ``display_text``
   documents and filterable metadata (gene, section, source, grade,
   assertion_type)

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
from typing import Any, Iterable

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import EvidenceGrade, EvidenceRecord

logger = logging.getLogger(__name__)

DEFAULT_COLLECTION = "evidence_records"
_TOKEN_RE = re.compile(r"[a-z0-9_]+", re.IGNORECASE)


# --------------------------------------------------------------------------------------
# Keyword / metadata retrieval (always available)
# --------------------------------------------------------------------------------------
def _norm(value: str | None) -> str:
    return " ".join((value or "").strip().lower().split())


def _enum_value(value: Any) -> str:
    if value is None:
        return ""
    return str(getattr(value, "value", value))


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


def _where_from_filters(
    *,
    gene_symbol: str | None = None,
    section: str | None = None,
    source_name: str | None = None,
    evidence_grade: EvidenceGrade | str | None = None,
    assertion_type: str | None = None,
    source_id: str | None = None,
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
    ) -> None:
        self.settings = settings or get_settings()
        self.collection_name = collection_name
        self.embedding_function = embedding_function or HashEmbeddingFunction()
        self._client: Any | None = None
        self._collection: Any | None = None
        self.status = ChromaIndexStatus(available=False, backend="unavailable")

        try:
            import chromadb
        except Exception as exc:  # noqa: BLE001
            self.status = ChromaIndexStatus(
                available=False,
                backend="unavailable",
                error=f"chromadb import failed: {exc}",
            )
            logger.warning("Chroma unavailable: %s", exc)
            return

        try:
            if ephemeral:
                self._client = chromadb.EphemeralClient()
                backend = "ephemeral"
                path_str = None
            else:
                path = Path(
                    persist_directory
                    if persist_directory is not None
                    else self.settings.index_path
                )
                path.mkdir(parents=True, exist_ok=True)
                self._client = chromadb.PersistentClient(path=str(path))
                backend = "persistent"
                path_str = str(path)
            self._collection = self._client.get_or_create_collection(
                name=self.collection_name,
                embedding_function=self.embedding_function,
                metadata={"hnsw:space": "cosine"},
            )
            self.status = ChromaIndexStatus(
                available=True,
                backend=backend,
                collection=self.collection_name,
                path=path_str,
                indexed_count=int(self._collection.count()),
            )
        except Exception as exc:  # noqa: BLE001
            self.status = ChromaIndexStatus(
                available=False,
                backend="unavailable",
                error=str(exc),
            )
            logger.warning("Chroma init failed; semantic search disabled: %s", exc)
            self._client = None
            self._collection = None

    @property
    def available(self) -> bool:
        return bool(self.status.available and self._collection is not None)

    def upsert_evidence(self, records: Iterable[EvidenceRecord]) -> int:
        """Index / update records by ``source_id``. Returns count upserted."""
        if not self.available or self._collection is None:
            return 0
        ids: list[str] = []
        documents: list[str] = []
        metadatas: list[dict[str, str]] = []
        for record in records:
            sid = (record.source_id or "").strip()
            text = (record.display_text or "").strip()
            if not sid or not text:
                continue
            ids.append(sid)
            documents.append(text)
            metadatas.append(_record_metadata(record))
        if not ids:
            return 0
        try:
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
    by_id: dict[str, RetrievalHit] = {h.source_id: h for h in keyword_hits}

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
            lookup = {r.source_id: r for r in record_list if r.source_id}
            semantic_hits = index.query(
                query,
                gene_symbol=gene_symbol,
                section=section,
                source_name=source_name,
                evidence_grade=evidence_grade,
                assertion_type=assertion_type,
                source_id=source_id,
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
    "KeywordEvidenceIndex",
    "HashEmbeddingFunction",
    "ChromaIndexStatus",
    "ChromaEvidenceIndex",
    "tokenize",
    "evidence_matches_filters",
    "keyword_score",
    "search_evidence_keyword",
    "index_evidence_in_chroma",
    "search_evidence",
]
