"""Bundle-only Section 1c official CDD/PDBe asset helpers.

This module is intentionally not wired into the full dossier workflow. It runs
only for section-bundle requests that explicitly include Section 1c.
"""

from __future__ import annotations

import html
import io
import json
import logging
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence
from urllib.parse import quote, urljoin, urlparse
from uuid import uuid4

import httpx
try:  # Optional dependency declared for real CDD pages; fallback keeps tests importable.
    from bs4 import BeautifulSoup
except Exception:  # pragma: no cover - exercised only when dependency missing
    BeautifulSoup = None  # type: ignore[assignment]
try:
    from PIL import Image, ImageStat
except Exception:  # pragma: no cover - exercised only when dependency missing
    Image = None  # type: ignore[assignment]
    ImageStat = None  # type: ignore[assignment]

from gene_dossier.config import Settings, get_settings
from gene_dossier.db import save_api_run, save_evidence_record, save_raw_artifact, session_scope
from gene_dossier.models import (
    ApiRun,
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    RawArtifact,
    SourceType,
    ToolResult,
)
from gene_dossier.raw_store import RawStore, compute_hash
from gene_dossier.source_ids import make_source_id, slugify
from gene_dossier.tools import cdd, pdbe
from gene_dossier.ucsc_figure import relative_to_artifact_root, sha256_hex, validate_image_bytes
from gene_dossier.workflow import DossierState, WorkflowTransientContext

logger = logging.getLogger(__name__)

SECTION_STRUCTURE = "Known structure / domains"
CDD_RECORD_URL = "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi"
PDBE_ENTRY_URL = "https://www.ebi.ac.uk/pdbe/entry/pdb"
OFFICIAL_IMAGE_HOSTS = {"www.ncbi.nlm.nih.gov", "www.ebi.ac.uk"}
CDD_IMAGE_ORIGINS = {
    "ncbi_cdd_official_thumbnail",
    "ncbi_cdd_official_feature_thumbnail",
}

_PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_ATTR_RE = re.compile(r"""(?:src|href)\s*=\s*["']([^"']+)["']""", re.IGNORECASE)
_TAG_RE = re.compile(r"<[^>]+>")
_CD_ACCESSION_RE = re.compile(r"\b(cd\d+|COG\d+|pfam\d+|smart\d+|TIGR\d+)\b", re.IGNORECASE)
_PSSM_RE = re.compile(r"(?:PSSM(?:-?ID)?|uid)\D{0,20}(\d{3,})", re.IGNORECASE)


@dataclass(frozen=True)
class ProteinSeed:
    """Authoritative protein seed used without observed-span length fallback."""

    uniprot_accession: str | None = None
    refseq_protein: str | None = None
    protein_length: int | None = None
    protein_length_source: str | None = None
    reviewed: bool | None = None
    taxon_id: int | None = None
    organism: str | None = None
    common_name: str | None = None
    source_evidence_record_id: str | None = None
    source_raw_artifact_id: str | None = None
    notes: tuple[str, ...] = ()

    @property
    def query_accession(self) -> str | None:
        return self.refseq_protein or self.uniprot_accession


@dataclass
class PdbCandidate:
    """A ranked experimental PDB candidate for a single species tier."""

    pdb_id: str
    chain_ids: list[str] = field(default_factory=list)
    exact_accession_mapping: bool = False
    mapped_taxon_human: bool | None = None
    calculated_coverage: float | None = None
    provider_coverage: float | None = None
    mapped_span_length: int = 0
    mapped_spans: list[tuple[int, int]] = field(default_factory=list)
    resolution: float | None = None
    experimental_method: str | None = None
    title: str | None = None
    preferred_assembly_id: str = "1"
    selected: bool = False
    rejection_reasons: list[str] = field(default_factory=list)
    selected_uniprot_accession: str | None = None
    species: str | None = None
    species_common_name: str | None = None
    taxon_id: int | None = None
    expression_host: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "chain_ids": list(self.chain_ids),
            "exact_accession_mapping": self.exact_accession_mapping,
            "mapped_taxon_human": self.mapped_taxon_human,
            "calculated_coverage": self.calculated_coverage,
            "provider_coverage": self.provider_coverage,
            "mapped_span_length": self.mapped_span_length,
            "mapped_spans": [list(span) for span in self.mapped_spans],
            "resolution": self.resolution,
            "experimental_method": self.experimental_method,
            "title": self.title,
            "preferred_assembly_id": self.preferred_assembly_id,
            "selected": self.selected,
            "rejection_reasons": list(self.rejection_reasons),
            "selected_uniprot_accession": self.selected_uniprot_accession,
            "species": self.species,
            "species_common_name": self.species_common_name,
            "taxon_id": self.taxon_id,
            "expression_host": self.expression_host,
        }


def _as_int(value: Any) -> int | None:
    try:
        if value is None or value == "":
            return None
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _as_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        number = float(str(value))
        if math.isnan(number) or math.isinf(number):
            return None
        return number
    except (TypeError, ValueError):
        return None


def _json_notes(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _first_text(value: Any) -> str | None:
    if isinstance(value, list):
        for item in value:
            if item:
                return str(item)
    if value:
        return str(value)
    return None


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    assertion_type: AssertionType,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    evidence_grade: EvidenceGrade,
    raw_artifact_id: str | None = None,
    api_run_id: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=make_source_id(source_name, gene_symbol, assertion_type, key),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_STRUCTURE,
        subsection="Known structure",
        source_name=source_name,
        source_type=SourceType.structure_database,
        assertion_type=assertion_type,
        fact_type=fact_type,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def _taxon_common(taxon_id: int | None) -> str | None:
    return {9606: "human", 10090: "mouse", 10116: "rat"}.get(taxon_id or 0)


def _taxon_scientific(taxon_id: int | None) -> str | None:
    return {
        9606: "Homo sapiens",
        10090: "Mus musculus",
        10116: "Rattus norvegicus",
    }.get(taxon_id or 0)


def _seed_from_uniprot_record(rec: EvidenceRecord) -> ProteinSeed | None:
    value = rec.value if isinstance(rec.value, dict) else {}
    accession = value.get("uniprot_accession") or value.get("primaryAccession")
    if not accession:
        return None
    taxon = _as_int(value.get("taxon_id") or value.get("tax_id") or rec.taxon_id)
    length = _as_int(value.get("protein_length") or value.get("sequence_length") or value.get("length"))
    refseqs = value.get("refseq_protein_accessions")
    refseq = None
    if isinstance(refseqs, list):
        refseq = next(
            (
                str(item).strip()
                for item in refseqs
                if str(item).strip().startswith(("NP_", "XP_", "YP_"))
            ),
            None,
        )
    reviewed = value.get("reviewed")
    return ProteinSeed(
        uniprot_accession=str(accession).strip(),
        refseq_protein=refseq,
        protein_length=length,
        protein_length_source="reviewed canonical UniProt sequence length" if length else None,
        reviewed=bool(reviewed) if reviewed is not None else None,
        taxon_id=taxon,
        organism=str(value.get("organism_name") or rec.organism or _taxon_scientific(taxon) or ""),
        common_name=_taxon_common(taxon),
        source_evidence_record_id=rec.id,
        source_raw_artifact_id=rec.raw_artifact_id,
    )


def protein_seeds_by_species(evidence_records: Sequence[EvidenceRecord]) -> list[ProteinSeed]:
    """Return reviewed UniProt seeds in human, mouse, rat tier order."""
    candidates: list[ProteinSeed] = []
    for rec in evidence_records:
        if rec.source_name != "UniProt" or rec.fact_type != "uniprot_accession":
            continue
        seed = _seed_from_uniprot_record(rec)
        if seed:
            candidates.append(seed)

    out: list[ProteinSeed] = []
    for taxon in (9606, 10090, 10116):
        tier = [seed for seed in candidates if seed.taxon_id == taxon]
        tier.sort(
            key=lambda seed: (
                0 if seed.reviewed else 1,
                0 if seed.protein_length is not None else 1,
                seed.uniprot_accession or "",
            )
        )
        if tier:
            out.append(tier[0])
    return out


def select_authoritative_protein_seed(
    *,
    gene_symbol: str,
    gene_ids: dict[str, Any],
    evidence_records: Sequence[EvidenceRecord],
) -> ProteinSeed:
    """Select the canonical human protein seed without observed-span fallback."""
    notes: list[str] = []
    for seed in protein_seeds_by_species(evidence_records):
        if seed.taxon_id == 9606:
            if seed.protein_length is None:
                notes.append("canonical protein length unavailable; calculated coverage set to null")
            return ProteinSeed(**{**seed.__dict__, "notes": tuple(notes)})

    accession = gene_ids.get("uniprot_accession")
    length = _as_int(gene_ids.get("protein_length"))
    refseq = gene_ids.get("refseq_protein")
    if accession or refseq:
        if length is None:
            notes.append("canonical protein length unavailable; calculated coverage set to null")
        return ProteinSeed(
            uniprot_accession=str(accession).strip() if accession else None,
            refseq_protein=str(refseq).strip() if refseq else None,
            protein_length=length,
            protein_length_source=str(gene_ids.get("protein_length_source") or "") if length else None,
            taxon_id=9606 if accession else None,
            organism="Homo sapiens",
            common_name="human" if accession else None,
            notes=tuple(notes),
        )
    return ProteinSeed(notes=("no authoritative human protein accession resolved",))


def cdd_domain_rows(
    evidence_records: Sequence[EvidenceRecord],
    *,
    protein_length: int | None,
) -> list[dict[str, Any]]:
    """Return normalized CDD domain rows; coverage is null without canonical length."""
    rows: list[dict[str, Any]] = []
    for rec in evidence_records:
        if rec.source_name != "CDD" or rec.fact_type != "conserved_domain_hit":
            continue
        value = rec.value if isinstance(rec.value, dict) else {}
        start = _as_int(value.get("from_residue"))
        end = _as_int(value.get("to_residue"))
        if start is not None and end is not None and end < start:
            start, end = end, start
        calculated = None
        if protein_length and start is not None and end is not None:
            calculated = max(0, end - start + 1) / protein_length
        domain_acc = value.get("domain_accession")
        pssm_id = _as_int(value.get("pssm_id") or value.get("domain_uid"))
        if pssm_id is None and domain_acc and str(domain_acc).isdigit():
            pssm_id = _as_int(domain_acc)
        rows.append(
            {
                "query_accession": value.get("query_accession"),
                "domain_accession": domain_acc,
                "pssm_id": pssm_id,
                "domain_short_name": value.get("domain_short_name")
                or domain_acc
                or "CDD domain",
                "domain_description": value.get("domain_description"),
                "from_residue": start,
                "to_residue": end,
                "evalue": value.get("evalue"),
                "bitscore": value.get("bitscore"),
                "superfamily": value.get("superfamily"),
                "coverage": calculated,
                "evidence_record_id": rec.id,
                "source_id": rec.source_id,
                "raw_artifact_id": rec.raw_artifact_id,
                "api_run_id": rec.api_run_id,
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            row["from_residue"] if row["from_residue"] is not None else 10**9,
            row["to_residue"] if row["to_residue"] is not None else 10**9,
            str(row["domain_accession"] or row["pssm_id"] or row["domain_short_name"]),
        ),
    )


def _union_length(spans: Iterable[tuple[int, int]]) -> int:
    ordered = sorted((min(a, b), max(a, b)) for a, b in spans if a and b)
    if not ordered:
        return 0
    merged: list[tuple[int, int]] = []
    for start, end in ordered:
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return sum(end - start + 1 for start, end in merged)


def _entry_summary_by_pdb(data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for entry in data.get("entry_summaries") or []:
        if isinstance(entry, dict) and entry.get("pdb_id"):
            out[str(entry["pdb_id"]).lower()] = entry
    return out


def _preferred_assembly_id(data: dict[str, Any], pdb_id: str) -> str:
    raw = data.get("entry_summaries_raw") or {}
    entries = None
    if isinstance(raw, dict):
        entries = raw.get(pdb_id) or raw.get(pdb_id.upper())
    if not isinstance(entries, list):
        return "1"
    assemblies: list[dict[str, Any]] = []
    for entry in entries:
        if isinstance(entry, dict):
            assemblies.extend(a for a in entry.get("assemblies") or [] if isinstance(a, dict))
    for assembly in assemblies:
        if assembly.get("preferred") and assembly.get("assembly_id") is not None:
            return str(assembly["assembly_id"])
    for assembly in assemblies:
        if assembly.get("assembly_id") is not None:
            return str(assembly["assembly_id"])
    return "1"


def _mapping_spans_for_accession(
    mappings_payload: Any,
    *,
    pdb_id: str,
    accession: str | None,
) -> tuple[list[tuple[int, int]], list[str], bool | None]:
    if not isinstance(mappings_payload, dict) or not accession:
        return [], [], None
    entry = mappings_payload.get(pdb_id) or mappings_payload.get(pdb_id.upper())
    if not isinstance(entry, dict):
        return [], [], None
    uniprot = entry.get("UniProt") or entry.get("uniprot") or {}
    if not isinstance(uniprot, dict):
        return [], [], None
    selected = uniprot.get(accession) or uniprot.get(accession.upper())
    if not isinstance(selected, dict):
        return [], [], None
    spans: list[tuple[int, int]] = []
    chains: list[str] = []
    taxon_human: bool | None = None
    taxon = _as_int(selected.get("tax_id") or selected.get("taxon_id"))
    if taxon is not None:
        taxon_human = taxon == 9606
    for mapping in selected.get("mappings") or []:
        if not isinstance(mapping, dict):
            continue
        start = _as_int(mapping.get("unp_start"))
        end = _as_int(mapping.get("unp_end"))
        if start is not None and end is not None:
            spans.append((min(start, end), max(start, end)))
        chain = mapping.get("chain_id") or mapping.get("struct_asym_id")
        if chain and str(chain) not in chains:
            chains.append(str(chain))
    return spans, chains, taxon_human


def rank_pdb_candidates(
    pdbe_payload: dict[str, Any],
    *,
    selected_uniprot_accession: str | None,
    protein_length: int | None,
) -> list[PdbCandidate]:
    """Rank PDB candidates within one species tier."""
    summaries = pdbe_payload.get("structure_summaries") or []
    if not isinstance(summaries, list):
        summaries = []
    mappings = pdbe_payload.get("uniprot_mappings") or {}
    entry_by_pdb = _entry_summary_by_pdb(pdbe_payload)
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in summaries:
        if not isinstance(row, dict) or not row.get("pdb_id"):
            continue
        pdb_id = str(row["pdb_id"]).lower()
        if _PDB_ID_RE.match(pdb_id):
            grouped.setdefault(pdb_id, []).append(row)

    candidates: list[PdbCandidate] = []
    for pdb_id, rows in grouped.items():
        spans, mapped_chains, taxon_human = _mapping_spans_for_accession(
            mappings.get(pdb_id),
            pdb_id=pdb_id,
            accession=selected_uniprot_accession,
        )
        row_spans: list[tuple[int, int]] = []
        row_chains: list[str] = []
        provider_coverages: list[float] = []
        resolutions: list[float] = []
        method = None
        for row in rows:
            chain = row.get("chain_id")
            if chain and str(chain) not in row_chains:
                row_chains.append(str(chain))
            start = _as_int(row.get("unp_start"))
            end = _as_int(row.get("unp_end"))
            if start is not None and end is not None:
                row_spans.append((min(start, end), max(start, end)))
            provider = _as_float(row.get("coverage"))
            if provider is not None:
                provider_coverages.append(provider)
            resolution = _as_float(row.get("resolution"))
            if resolution is not None:
                resolutions.append(resolution)
            method = method or row.get("experimental_method")

        exact = bool(spans)
        effective_spans = spans if exact else row_spans
        span_len = _union_length(effective_spans)
        entry = entry_by_pdb.get(pdb_id) or {}
        candidate = PdbCandidate(
            pdb_id=pdb_id,
            chain_ids=sorted(mapped_chains if mapped_chains else row_chains),
            exact_accession_mapping=exact,
            mapped_taxon_human=taxon_human,
            calculated_coverage=(span_len / protein_length) if protein_length and span_len else None,
            provider_coverage=max(provider_coverages) if provider_coverages else None,
            mapped_span_length=span_len,
            mapped_spans=sorted(set(effective_spans)),
            resolution=min(resolutions) if resolutions else None,
            experimental_method=method
            or _first_text(entry.get("experimental_method"))
            or _first_text(entry.get("experimental_method_class")),
            title=entry.get("title"),
            preferred_assembly_id=_preferred_assembly_id(pdbe_payload, pdb_id),
            selected_uniprot_accession=selected_uniprot_accession,
        )
        if not exact:
            candidate.rejection_reasons.append("missing exact mapping to selected UniProt accession")
        candidates.append(candidate)

    def _rank_key(candidate: PdbCandidate) -> tuple[Any, ...]:
        resolution = candidate.resolution if candidate.resolution is not None else float("inf")
        known_coverage = candidate.calculated_coverage is not None
        return (
            0 if candidate.exact_accession_mapping else 1,
            -candidate.mapped_span_length,
            0 if known_coverage else 1,
            -(candidate.calculated_coverage or candidate.provider_coverage or -1.0),
            resolution,
            candidate.pdb_id,
        )

    ranked = sorted(candidates, key=_rank_key)
    for idx, candidate in enumerate(ranked):
        if idx == 0 and candidate.exact_accession_mapping:
            candidate.selected = True
        elif candidate.exact_accession_mapping:
            candidate.rejection_reasons.append("lower ranked than selected candidate")
    return ranked


def _validate_image(content: bytes) -> dict[str, Any]:
    image, err = validate_image_bytes(content, min_width=32, min_height=16)
    if err or image is None:
        raise ValueError(err.message if err else "invalid image")
    return {
        "media_type": image.media_type,
        "width": image.width,
        "height": image.height,
        "byte_size": image.byte_size,
    }


def _validate_nonblank_image(content: bytes) -> dict[str, Any]:
    meta = _validate_image(content)
    if Image is not None and ImageStat is not None:
        with Image.open(io.BytesIO(content)) as img:
            converted = img.convert("RGBA")
            stat = ImageStat.Stat(converted)
            extrema = converted.getextrema()
            alpha_range = extrema[3] if len(extrema) > 3 else (255, 255)
            variance = sum(stat.var[:3])
            if alpha_range == (0, 0):
                raise ValueError("image is fully transparent")
            if variance < 2.0:
                raise ValueError("image appears blank or low entropy")
            meta["pixel_variance"] = variance
    return meta


def _validate_html(content: bytes) -> dict[str, Any]:
    text = content.decode("utf-8", errors="replace")
    lower = text.lower()
    if not any(marker in lower for marker in ("<html", "<body", "<script", "<div", "<!doctype")):
        raise ValueError("HTML payload does not look like a page")
    return {"media_type": "text/html", "byte_size": len(content)}


def _validate_json(content: bytes) -> dict[str, Any]:
    json.loads(content.decode("utf-8"))
    return {"media_type": "application/json", "byte_size": len(content)}


def _artifact_relative_path(path: Path, *, settings: Settings) -> str:
    return relative_to_artifact_root(path, root=settings.raw_data_path)


def _persist_artifact_bytes(
    *,
    dossier_run_id: str,
    source_name: str,
    content: bytes,
    extension: str,
    artifact_type: str,
    filename_hint: str,
    settings: Settings,
    api_run: ApiRun | None,
    persist_db: bool,
    notes: dict[str, Any],
    validate: Callable[[bytes], dict[str, Any]],
) -> tuple[RawArtifact, dict[str, Any]]:
    """Stage, validate, atomically move, and persist artifact metadata."""
    store = RawStore(base_dir=settings.raw_data_path)
    digest = compute_hash(content)
    ext = extension.lstrip(".")
    final_dir = store._dir_for(dossier_run_id, source_name)
    hint = slugify(filename_hint) or "artifact"
    final_path = final_dir / f"{hint}-{digest[:12]}.{ext}"
    staging_dir = settings.raw_data_path / "_staging" / slugify(dossier_run_id)
    staging_dir.mkdir(parents=True, exist_ok=True)
    stage_path = staging_dir / f".{hint}-{uuid4().hex}.{ext}.tmp"
    created_final = False
    existed = final_path.exists()
    try:
        stage_path.write_bytes(content)
        if sha256_hex(stage_path.read_bytes()) != digest:
            raise ValueError("staged artifact checksum mismatch")
        validation = validate(content)
        if existed:
            if sha256_hex(final_path.read_bytes()) != digest:
                raise ValueError(f"existing artifact checksum mismatch at {final_path}")
            stage_path.unlink(missing_ok=True)
        else:
            stage_path.replace(final_path)
            created_final = True
        if sha256_hex(final_path.read_bytes()) != digest:
            raise ValueError("final artifact checksum mismatch")

        notes_payload = {
            **notes,
            "expected_sha256": digest,
            "media_type": validation.get("media_type"),
            "byte_size": validation.get("byte_size", len(content)),
        }
        artifact = RawArtifact(
            dossier_run_id=dossier_run_id,
            api_run_id=api_run.id if api_run else None,
            source_name=source_name,
            artifact_type=artifact_type,
            file_path=str(final_path),
            original_url=api_run.request_url if api_run and api_run.request_url else None,
            content_hash=digest,
            notes=_json_notes(notes_payload),
        )
        if api_run is not None:
            api_run.raw_artifact_id = artifact.id
        if persist_db:
            with session_scope() as session:
                if api_run is not None:
                    save_api_run(session, api_run)
                save_raw_artifact(session, artifact)
        meta = artifact.model_dump(mode="json")
        meta.update(
            {
                "relative_path": _artifact_relative_path(final_path, settings=settings),
                "expected_sha256": digest,
                "media_type": validation.get("media_type"),
                "byte_size": validation.get("byte_size", len(content)),
                "artifact_class": notes.get("artifact_class"),
                "artifact_role": notes.get("artifact_role"),
                "artifact_origin": notes.get("artifact_origin"),
                "retrieval_method": notes.get("retrieval_method"),
                "source_url": notes.get("source_url"),
                "source_page_url": notes.get("source_page_url"),
                "parent_raw_artifact_ids": list(notes.get("parent_raw_artifact_ids") or []),
                "parent_evidence_record_ids": list(notes.get("parent_evidence_record_ids") or []),
                "width": validation.get("width"),
                "height": validation.get("height"),
            }
        )
        return artifact, meta
    except Exception:
        stage_path.unlink(missing_ok=True)
        if created_final and not existed:
            final_path.unlink(missing_ok=True)
        raise


def _save_api_run_failure(api_run: ApiRun, *, persist_db: bool) -> None:
    if not persist_db:
        return
    try:
        with session_scope() as session:
            save_api_run(session, api_run)
    except Exception as exc:  # noqa: BLE001
        logger.warning("DB API-run failure persist failed for %s: %s", api_run.endpoint_name, exc)


def _append_evidence(evidence: list[EvidenceRecord], rec: EvidenceRecord, *, persist_db: bool) -> bool:
    if persist_db:
        try:
            with session_scope() as session:
                save_evidence_record(session, rec)
        except Exception as exc:  # noqa: BLE001
            logger.warning("DB derived evidence persist failed for %s: %s", rec.source_id, exc)
            return False
    evidence.append(rec)
    return True


def _tool_result_to_api_run(
    tr: ToolResult,
    *,
    dossier_run_id: str,
    gene_symbol: str,
) -> ApiRun:
    return ApiRun(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name=tr.source_name,
        endpoint_name=tr.endpoint_name,
        request_url=tr.request_url,
        request_params=tr.request_params or {},
        status_code=tr.status_code,
        success=tr.success,
        error_type=tr.error_type,
        error_message=tr.error_message,
    )


def _raw_artifact_id_for_payload(state: DossierState, source_name: str, endpoint_name: str) -> str | None:
    for meta in state.get("raw_artifacts") or []:
        if meta.get("source_name") != source_name:
            continue
        notes = meta.get("notes")
        notes_text = notes if isinstance(notes, str) else json.dumps(notes or {})
        if endpoint_name in notes_text or meta.get("id"):
            return str(meta.get("id"))
    return None


def _first_raw_artifact_id_for_source(state: DossierState, source_name: str) -> str | None:
    for meta in state.get("raw_artifacts") or []:
        if meta.get("source_name") == source_name and meta.get("id"):
            return str(meta["id"])
    return None


def _first_api_run_id_for_source(state: DossierState, source_name: str) -> str | None:
    for api in state.get("api_runs") or []:
        if api.source_name == source_name and api.id:
            return api.id
    return None


def _cdd_payload(state: DossierState) -> dict[str, Any] | None:
    candidates: list[dict[str, Any]] = []
    for tr in state.get("tool_results") or []:
        if tr.source_name == "CDD" and tr.success and isinstance(tr.data, dict):
            candidates.append(tr.data)
    for payload in candidates:
        if payload.get("hit_summaries") or payload.get("feature_summaries"):
            return payload
    for payload in candidates:
        if payload.get("endpoint_name") == "fetch_domains" or payload.get("hits") or payload.get("features"):
            return payload
    return candidates[0] if candidates else None


def cdd_master_lineage(payload: dict[str, Any] | None) -> dict[str, Any]:
    payload = payload or {}
    master = payload.get("master_cdsid") or payload.get("cdsid")
    return {
        "master_cdsid": master,
        "hits_request_cdsid": payload.get("hits_request_cdsid") or master,
        "features_request_cdsid": payload.get("features_request_cdsid") or master,
        "graphical_result_master_cdsid": master,
        "same_master_job": bool(master),
        "query_index": 0,
    }


def _fetch_bytes(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    endpoint_name: str,
    url: str,
    request_params: dict[str, Any],
    settings: Settings,
    allowed_hosts: set[str] | None = None,
) -> tuple[ApiRun, bytes | None, str | None]:
    parsed = urlparse(url)
    if parsed.scheme != "https" or (allowed_hosts and parsed.netloc not in allowed_hosts):
        api = ApiRun(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            source_name=source_name,
            endpoint_name=endpoint_name,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="invalid_url",
            error_message=f"URL host is not allowlisted: {parsed.netloc}",
        )
        return api, None, None
    api = ApiRun(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name=source_name,
        endpoint_name=endpoint_name,
        request_url=url,
        request_params=request_params,
        success=False,
    )
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds, follow_redirects=True) as client:
            response = client.get(url)
        api.status_code = response.status_code
        content_type = response.headers.get("content-type") or ""
        if not response.is_success:
            api.error_type = "http_error"
            api.error_message = f"HTTP {response.status_code}"
            return api, None, content_type
        api.success = True
        return api, response.content, content_type
    except httpx.HTTPError as exc:
        api.error_type = exc.__class__.__name__
        api.error_message = str(exc)[:400]
        return api, None, None


def _download_official_image(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    url: str,
    source_page_url: str | None,
    artifact_origin: str,
    filename_hint: str,
    settings: Settings,
    persist_db: bool,
    endpoint_name: str = "download_official_image",
    request_params: dict[str, Any] | None = None,
    parent_raw_artifact_ids: Sequence[str] = (),
    parent_evidence_record_ids: Sequence[str] = (),
    retrieval_method: str = "direct_image_download",
) -> tuple[ApiRun, RawArtifact | None, dict[str, Any] | None]:
    api, content, content_type = _fetch_bytes(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name=source_name,
        endpoint_name=endpoint_name,
        url=url,
        request_params=request_params or {},
        settings=settings,
        allowed_hosts=OFFICIAL_IMAGE_HOSTS,
    )
    if not api.success or content is None:
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None, None
    if not str(content_type or "").lower().startswith("image/"):
        api.success = False
        api.error_type = "invalid_content_type"
        api.error_message = f"Expected image/*, got {content_type or 'unknown'}"
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None, None
    if content[:512].lstrip().lower().startswith((b"<!doctype html", b"<html")):
        api.success = False
        api.error_type = "html_error_page"
        api.error_message = "Image endpoint returned HTML"
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None, None
    if len(content) > 8_000_000:
        api.success = False
        api.error_type = "image_too_large"
        api.error_message = "Official image exceeds byte-size bound"
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None, None
    artifact, meta = _persist_artifact_bytes(
        dossier_run_id=dossier_run_id,
        source_name=source_name,
        content=content,
        extension="png" if "png" in str(content_type).lower() else "jpg",
        artifact_type="image",
        filename_hint=filename_hint,
        settings=settings,
        api_run=api,
        persist_db=persist_db,
        notes={
            "artifact_class": "external_raw",
            "artifact_origin": artifact_origin,
            "artifact_role": artifact_origin,
            "source_url": url,
            "source_page_url": source_page_url,
            "retrieval_method": retrieval_method,
            "parent_raw_artifact_ids": list(parent_raw_artifact_ids),
            "parent_evidence_record_ids": list(parent_evidence_record_ids),
        },
        validate=_validate_image,
    )
    return api, artifact, meta


def _strip_tags(text: str) -> str:
    cleaned = re.sub(r"<script\b[^>]*>.*?</script>", " ", text or " ", flags=re.IGNORECASE | re.DOTALL)
    cleaned = re.sub(r"<style\b[^>]*>.*?</style>", " ", cleaned, flags=re.IGNORECASE | re.DOTALL)
    return html.unescape(_TAG_RE.sub(" ", cleaned))


def _compact_text(text: str | None, *, limit: int = 950) -> str | None:
    compact = re.sub(r"\s+", " ", html.unescape(text or "")).strip()
    compact = re.sub(r"\b(Links|Source|Taxonomy|PubMed|Protein|Structure|Statistics)\s+\?.*$", "", compact).strip()
    if not compact:
        return None
    return compact[: limit - 1].rstrip() + "..." if len(compact) > limit else compact


def _sentence_safe_truncate(text: str | None, *, limit: int = 1200) -> str | None:
    compact = re.sub(r"\s+", " ", html.unescape(text or "")).strip()
    if not compact:
        return None
    if len(compact) <= limit:
        return compact
    clipped = compact[:limit].rstrip()
    sentence_end = max(clipped.rfind("."), clipped.rfind(";"), clipped.rfind("!"), clipped.rfind("?"))
    if sentence_end >= max(120, int(limit * 0.55)):
        return clipped[: sentence_end + 1].strip()
    word_end = clipped.rfind(" ")
    if word_end > 0:
        clipped = clipped[:word_end].rstrip()
    return clipped + "..."


def _thumbnail_role_from_url_context(
    *,
    url: str,
    context_text: str,
    accession: str | None,
    pssm_id: str | None,
) -> tuple[str, list[str], str | None]:
    lower_url = url.lower()
    lower_context = context_text.lower()
    signals: list[str] = []
    if accession and accession.lower() in lower_context:
        signals.append(f"linked to family {accession.lower()}")
    if pssm_id and pssm_id in lower_url:
        signals.append(f"linked to PSSM {pssm_id}")
    if "ft=" in lower_url:
        signals.append("feature thumbnail parameter")
        return "conserved_feature_structure_thumbnail", signals, None
    if "cn3d" in lower_context or "interactive view" in lower_context:
        signals.append("inside Cn3D structure thumbnail context")
        if any(token in lower_context for token in ("feature", "site", "binding", "ion", "conserved site")):
            return "conserved_feature_structure_thumbnail", signals, None
        return "family_structure_thumbnail", signals, None
    if any(token in lower_context for token in ("alignment", "sequence", "logo", "msa")):
        signals.append("alignment/sequence context")
        return "alignment_or_sequence_thumbnail", signals, "alignment or sequence thumbnail"
    if any(token in lower_url for token in ("align", "seqgraphic", "sequence", "logo")):
        signals.append("alignment/sequence URL")
        return "alignment_or_sequence_thumbnail", signals, "alignment or sequence thumbnail"
    if any(token in lower_context for token in ("feature", "site", "binding", "ion")) or "ft=" in lower_url:
        signals.append("inside conserved-feature context")
        return "conserved_feature_structure_thumbnail", signals, None
    if any(token in lower_context for token in ("structure", "3d", "thumbnail", "domain")):
        signals.append("inside family structure context")
        return "family_structure_thumbnail", signals, None
    if "cdthumbnail.cgi" in lower_url and "img=1" in lower_url:
        signals.append("thumbnail image parameter")
        return "conserved_feature_structure_thumbnail", signals, None
    if "cdthumbnail.cgi" in lower_url:
        signals.append("CDD thumbnail endpoint")
        return "family_structure_thumbnail", signals, None
    return "other", signals, "unrecognized thumbnail context"


def parse_cdd_family_html(
    raw_html: str,
    *,
    requested_uid: str,
    source_page_url: str,
) -> dict[str, Any]:
    """Best-effort adapter over NCBI CDD family HTML."""
    text = _strip_tags(raw_html)
    accession = None
    acc_match = _CD_ACCESSION_RE.search(text)
    if acc_match:
        accession = acc_match.group(1)
    pssm_match = _PSSM_RE.search(text)
    pssm_id = pssm_match.group(1) if pssm_match else None
    thumbnail_candidates: list[dict[str, Any]] = []
    synopsis = None
    soup = BeautifulSoup(raw_html or "", "html.parser") if BeautifulSoup is not None else None
    if soup is not None:
        for tag in soup(["script", "style"]):
            tag.decompose()
        summary_container = soup.select_one(".summary-large-image")
        if summary_container is not None:
            synopsis = _sentence_safe_truncate(summary_container.get_text(" ", strip=True))
        if synopsis is None:
            for heading in soup.find_all(["h2", "h3", "h4", "b", "strong"]):
                heading_text = heading.get_text(" ", strip=True)
                if re.search(r"\b(summary|description)\b", heading_text, re.IGNORECASE):
                    pieces = []
                    for sibling in heading.find_all_next(["p", "div"], limit=4):
                        value = sibling.get_text(" ", strip=True)
                        if value:
                            pieces.append(value)
                    synopsis = _sentence_safe_truncate(" ".join(pieces))
                    if synopsis:
                        break
        for tag in soup.find_all(["img", "a"]):
            value = tag.get("src") or tag.get("href")
            if not value or "cdThumbnail.cgi" not in value:
                continue
            url = urljoin(source_page_url, html.unescape(str(value)))
            parent_text = ""
            parent = tag.parent
            for _ in range(3):
                if parent is None:
                    break
                parent_text = " ".join([parent_text, parent.get_text(" ", strip=True)]).strip()
                parent = parent.parent
            alt_title = " ".join(
                str(tag.get(attr) or "") for attr in ("alt", "title") if tag.get(attr)
            ).strip()
            role, signals, rejection = _thumbnail_role_from_url_context(
                url=url,
                context_text=" ".join([alt_title, parent_text]),
                accession=accession,
                pssm_id=pssm_id,
            )
            thumbnail_candidates.append(
                {
                    "url": url,
                    "classified_role": role,
                    "classification_signals": signals,
                    "rejection_reason": rejection,
                    "context_text": _sentence_safe_truncate(" ".join([alt_title, parent_text]), limit=240),
                }
            )
    thumbnail_urls = []
    for value in _ATTR_RE.findall(raw_html or ""):
        if "cdThumbnail.cgi" in value:
            url = urljoin(source_page_url, html.unescape(value))
            if not any(item.get("url") == url for item in thumbnail_candidates):
                role, signals, rejection = _thumbnail_role_from_url_context(
                    url=url,
                    context_text="",
                    accession=accession,
                    pssm_id=pssm_id,
                )
                thumbnail_candidates.append(
                    {
                        "url": url,
                        "classified_role": role,
                        "classification_signals": signals,
                        "rejection_reason": rejection,
                    }
                )
            thumbnail_urls.append(url)
    summary_match = re.search(
        r"summary-large-image.*?<span>\s*(?:<div[^>]*class=[\"']desctit[\"'][^>]*>(.*?)</div>)?(.*?)</span>",
        raw_html or "",
        re.IGNORECASE | re.DOTALL,
    )
    if synopsis is None and summary_match:
        title = _strip_tags(summary_match.group(1) or "")
        body = _strip_tags(summary_match.group(2) or "")
        synopsis = _sentence_safe_truncate(f"{title}: {body}" if title else body)
    if synopsis is None:
        meta_match = re.search(
            r'<meta\s+name=["\']description["\']\s+content=["\'](.*?)["\']',
            raw_html or "",
            re.IGNORECASE | re.DOTALL,
        )
        if meta_match:
            synopsis = _sentence_safe_truncate(meta_match.group(1))
    if synopsis is None:
        for marker in ("Summary", "Description", "This model", "This domain"):
            idx = text.lower().find(marker.lower())
            if idx >= 0:
                synopsis = _sentence_safe_truncate(text[idx : idx + 1200])
                break
    if synopsis is None:
        synopsis = _sentence_safe_truncate(text[:1200])
    superfamily = None
    sf_match = re.search(r"superfamily[:\s]+([^.;]{3,120})", text, re.IGNORECASE)
    if sf_match:
        superfamily = sf_match.group(1).strip()
    fallback_uid = pssm_id or (requested_uid if str(requested_uid or "").isdigit() else None)
    if fallback_uid and not any(f"uid={fallback_uid}" in url for url in thumbnail_urls):
        fallback_url = f"https://www.ncbi.nlm.nih.gov/Structure/cdd/cdThumbnail.cgi?uid={quote(fallback_uid, safe='')}"
        thumbnail_urls.append(fallback_url)
        thumbnail_candidates.append(
            {
                "url": fallback_url,
                "classified_role": "family_structure_thumbnail",
                "classification_signals": ["PSSM fallback thumbnail endpoint"],
                "rejection_reason": None,
            }
        )
    return {
        "parser_name": "ncbi_cdd_family_html",
        "parser_version": "1",
        "requested_uid": requested_uid,
        "canonical_accession": accession,
        "pssm_id": pssm_id,
        "synopsis": synopsis,
        "superfamily": superfamily,
        "thumbnail_urls": list(dict.fromkeys(thumbnail_urls)),
        "thumbnail_candidates": thumbnail_candidates,
        "unparsed_sections_preserved": True,
    }


def _family_url(domain: dict[str, Any]) -> str | None:
    uid = domain.get("pssm_id")
    accession = domain.get("domain_accession")
    requested = uid or accession
    if not requested:
        return None
    return f"{CDD_RECORD_URL}?uid={quote(str(requested), safe='')}"


def _domain_key(domain: dict[str, Any]) -> str:
    acc = domain.get("domain_accession")
    if acc and not str(acc).isdigit():
        return str(acc).lower()
    pssm = domain.get("pssm_id")
    if pssm:
        return f"pssm:{pssm}"
    return str(acc or domain.get("domain_short_name") or "domain").lower()


def _safe_item_token(value: Any, *, fallback: str) -> str:
    raw = str(value or fallback).lower()
    token = re.sub(r"[^a-z0-9]+", "-", raw).strip("-")
    return token or fallback


def _domain_item_key(domain: dict[str, Any], parsed: dict[str, Any] | None = None) -> str:
    parsed = parsed or {}
    acc = parsed.get("canonical_accession") or domain.get("domain_accession") or domain.get("pssm_id")
    return f"domain-{_safe_item_token(acc, fallback='unknown')}"


def _feature_item_key(feature: dict[str, Any], domain_acc: Any) -> str:
    label = feature.get("feature_label") or feature.get("feature_name") or feature.get("feature_type")
    return f"feature-{_safe_item_token(domain_acc, fallback='unknown')}-{_safe_item_token(label, fallback='feature')}"


def _is_polished_feature(feature: dict[str, Any], *, matched_family_key: str | None) -> bool:
    name = str(feature.get("feature_name") or "").strip()
    label = str(feature.get("feature_label") or "").strip()
    type_or_desc = str(feature.get("feature_type") or feature.get("description") or "").strip()
    meaningful_name = bool(name and name.lower() not in {"-", "none", "conserved feature"})
    return bool((meaningful_name or label) and type_or_desc and matched_family_key)


def _select_thumbnail_candidate(
    candidates: Sequence[dict[str, Any]],
    *,
    preferred_roles: Sequence[str],
) -> dict[str, Any] | None:
    for role in preferred_roles:
        role_candidates = [
            candidate
            for candidate in candidates
            if candidate.get("classified_role") == role and not candidate.get("rejection_reason")
        ]
        if role_candidates:
            return sorted(
                role_candidates,
                key=lambda candidate: (
                    0
                    if any(
                        "pssm" in str(signal).lower()
                        for signal in candidate.get("classification_signals") or []
                    )
                    else 1,
                    str(candidate.get("url") or ""),
                ),
            )[0]
    return None


def _dedup_domains(domains: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    for domain in domains:
        key = _domain_key(domain)
        if key in seen:
            continue
        seen.add(key)
        out.append(domain)
    return out


def _feature_summaries(state: DossierState) -> list[dict[str, Any]]:
    payload = _cdd_payload(state) or {}
    features = payload.get("feature_summaries") or []
    return [row for row in features if isinstance(row, dict)]


def _feature_matches_domain(feature: dict[str, Any], domain: dict[str, Any]) -> bool:
    f_acc = str(feature.get("domain_accession") or "").lower()
    d_acc = str(domain.get("domain_accession") or "").lower()
    f_pssm = str(feature.get("pssm_id") or "")
    d_pssm = str(domain.get("pssm_id") or "")
    return bool((f_acc and d_acc and f_acc == d_acc) or (f_pssm and d_pssm and f_pssm == d_pssm))


def _fetch_cdd_family_enrichment(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    domain: dict[str, Any],
    settings: Settings,
    persist_db: bool,
) -> tuple[ApiRun | None, RawArtifact | None, dict[str, Any] | None, dict[str, Any]]:
    url = _family_url(domain)
    if not url:
        return None, None, None, {"status": "unavailable", "reason": "no domain UID/accession"}
    api, content, content_type = _fetch_bytes(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="CDD",
        endpoint_name="fetch_cdd_family_html",
        url=url,
        request_params={"domain_accession": domain.get("domain_accession"), "pssm_id": domain.get("pssm_id")},
        settings=settings,
        allowed_hosts={"www.ncbi.nlm.nih.gov"},
    )
    if not api.success or content is None:
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None, None, {"status": "failed", "reason": api.error_message}
    if "html" not in str(content_type or "").lower() and not content[:512].lower().lstrip().startswith(b"<"):
        api.success = False
        api.error_type = "invalid_content_type"
        api.error_message = f"Expected CDD HTML, got {content_type or 'unknown'}"
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None, None, {"status": "failed", "reason": api.error_message}
    artifact, meta = _persist_artifact_bytes(
        dossier_run_id=dossier_run_id,
        source_name="CDD",
        content=content,
        extension="html",
        artifact_type="html",
        filename_hint=f"cdd-family-{domain.get('domain_accession') or domain.get('pssm_id')}",
        settings=settings,
        api_run=api,
        persist_db=persist_db,
        notes={
            "artifact_class": "external_raw",
            "artifact_origin": "ncbi_cdd_family_html",
            "artifact_role": "ncbi_cdd_family_html",
            "source_url": url,
            "retrieval_method": "html_adapter",
        },
        validate=_validate_html,
    )
    parsed = parse_cdd_family_html(
        content.decode("utf-8", errors="replace"),
        requested_uid=str(domain.get("pssm_id") or domain.get("domain_accession") or ""),
        source_page_url=url,
    )
    parsed["raw_html_artifact_id"] = artifact.id
    return api, artifact, meta, parsed


_ARCH_LEFT_MARGIN = 130
_ARCH_RIGHT_MARGIN = 40
_ARCH_TOP_MARGIN = 28
_ARCH_PLOT_WIDTH = 900
_ARCH_LANE_HEIGHT = 42
_ARCH_BAR_HEIGHT = 16
_ARCH_QUERY_HEIGHT = 10
_ARCH_RENDERER_VERSION = "cdd_domain_architecture_render_v1"


def _architecture_x(residue: float, *, protein_length: int) -> float:
    clamped = max(0.0, min(float(residue), float(protein_length)))
    return _ARCH_LEFT_MARGIN + (clamped / float(protein_length)) * _ARCH_PLOT_WIDTH


def _domain_as_dict(item: Any) -> dict[str, Any] | None:
    if isinstance(item, EvidenceRecord):
        value = item.value if isinstance(item.value, dict) else {}
        return {
            **value,
            "evidence_record_id": item.id,
            "from_residue": value.get("from_residue"),
            "to_residue": value.get("to_residue"),
            "domain_accession": value.get("domain_accession"),
            "domain_short_name": value.get("domain_short_name"),
            "superfamily": value.get("superfamily"),
        }
    if isinstance(item, dict):
        return item
    return None


def _hit_lane_class(domain: dict[str, Any]) -> str:
    """Classify a conserved-domain hit as specific or superfamily for architecture lanes."""
    acc = str(domain.get("domain_accession") or "").strip().lower()
    if acc.startswith("cl"):
        return "superfamily"
    raw = domain.get("raw") if isinstance(domain.get("raw"), dict) else {}
    hit_type = str(
        domain.get("hit_type")
        or raw.get("Hit type")
        or raw.get("hit type")
        or ""
    ).strip().lower()
    if hit_type in {"superfamily", "superfamily hit"}:
        return "superfamily"
    sf = str(domain.get("superfamily") or "").strip().lower()
    if sf in {"superfamily", "superfamily hit"}:
        return "superfamily"
    return "specific"


def _companion_superfamily_from_hit(domain: dict[str, Any]) -> dict[str, Any] | None:
    if _hit_lane_class(domain) != "specific":
        return None
    start = _as_int(domain.get("from_residue"))
    end = _as_int(domain.get("to_residue"))
    if start is None or end is None:
        return None
    sf = str(domain.get("superfamily") or "").strip()
    match = re.search(r"\b(cl\d+)\b", sf, re.IGNORECASE)
    if not match:
        return None
    accession = match.group(1)
    name_match = re.search(r"\bcl\d+\b\s*[-:,]?\s*([A-Za-z0-9_./+-]+)", sf, re.IGNORECASE)
    short_name = name_match.group(1) if name_match else accession
    if short_name.lower() in {"superfamily", "hit", "specific"}:
        short_name = accession
    return {
        "domain_accession": accession,
        "domain_short_name": short_name,
        "from_residue": start,
        "to_residue": end,
        "superfamily": accession,
        "derived_from_specific_hit": True,
        "evidence_record_id": domain.get("evidence_record_id"),
    }


def _feature_marker_positions(feature: dict[str, Any]) -> list[int]:
    start = _as_int(feature.get("from_residue"))
    end = _as_int(feature.get("to_residue"))
    if start is not None and end is not None:
        return [int(round((start + end) / 2))]
    if start is not None:
        return [start]
    residues = str(feature.get("query_residues") or "")
    numbers = [int(n) for n in re.findall(r"\d+", residues)]
    return numbers


def _architecture_lane_rows(
    domains: Sequence[Any],
    features: Sequence[Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Build specific-hit, superfamily, and feature rows without prose-level dedup."""
    specific: list[dict[str, Any]] = []
    superfamily: list[dict[str, Any]] = []
    seen_sf_keys: set[tuple[Any, ...]] = set()

    for item in domains:
        domain = _domain_as_dict(item)
        if domain is None:
            continue
        start = _as_int(domain.get("from_residue"))
        end = _as_int(domain.get("to_residue"))
        if start is None or end is None:
            continue
        if end < start:
            start, end = end, start
        row = {
            **domain,
            "from_residue": start,
            "to_residue": end,
            "domain_short_name": domain.get("domain_short_name")
            or domain.get("domain_accession")
            or "domain",
        }
        if _hit_lane_class(domain) == "superfamily":
            key = (row.get("domain_accession"), start, end)
            if key not in seen_sf_keys:
                seen_sf_keys.add(key)
                superfamily.append(row)
        else:
            # Retain every specific-hit instance (repeats stay).
            specific.append(row)
            companion = _companion_superfamily_from_hit(row)
            if companion is not None:
                key = (
                    companion.get("domain_accession"),
                    companion["from_residue"],
                    companion["to_residue"],
                )
                if key not in seen_sf_keys:
                    seen_sf_keys.add(key)
                    superfamily.append(companion)

    feature_rows: list[dict[str, Any]] = []
    for item in features:
        feature = _domain_as_dict(item) if not isinstance(item, dict) else item
        if feature is None:
            continue
        positions = _feature_marker_positions(feature)
        if not positions:
            continue
        label = (
            feature.get("feature_label")
            or feature.get("feature_name")
            or feature.get("feature_type")
            or "feature"
        )
        for pos in positions:
            feature_rows.append(
                {
                    **feature,
                    "position": pos,
                    "feature_label": label,
                }
            )
    return specific, superfamily, feature_rows


def _build_cdd_architecture_svg_bytes(
    *,
    protein_length: int,
    gene_symbol: str,
    specific_hits: Sequence[dict[str, Any]],
    superfamily_hits: Sequence[dict[str, Any]],
    feature_markers: Sequence[dict[str, Any]],
) -> tuple[bytes, int, int]:
    lanes: list[tuple[str, str]] = [("query", "Query sequence")]
    if specific_hits:
        lanes.append(("specific", "Specific hits"))
    if superfamily_hits:
        lanes.append(("superfamily", "Superfamilies"))
    if feature_markers:
        lanes.append(("features", "Conserved features"))

    height = _ARCH_TOP_MARGIN + len(lanes) * _ARCH_LANE_HEIGHT + 36
    width = _ARCH_LEFT_MARGIN + _ARCH_PLOT_WIDTH + _ARCH_RIGHT_MARGIN
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(gene_symbol)} '
        f'CDD domain architecture">'
        f"<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>"
        f'<text x="8" y="18" font-size="12" font-family="Helvetica, Arial, sans-serif" fill="#333">'
        f"{html.escape(gene_symbol)} domain architecture (derived)</text>"
    ]

    for idx, (lane_id, label) in enumerate(lanes):
        y = _ARCH_TOP_MARGIN + idx * _ARCH_LANE_HEIGHT
        mid = y + _ARCH_LANE_HEIGHT / 2
        parts.append(
            f'<text x="8" y="{mid + 4:.1f}" font-size="11" '
            f'font-family="Helvetica, Arial, sans-serif" fill="#444">{html.escape(label)}</text>'
        )
        if lane_id == "query":
            x0 = _architecture_x(0, protein_length=protein_length)
            x1 = _architecture_x(protein_length, protein_length=protein_length)
            parts.append(
                f'<rect x="{x0:.2f}" y="{mid - _ARCH_QUERY_HEIGHT / 2:.2f}" '
                f'width="{max(1.0, x1 - x0):.2f}" height="{_ARCH_QUERY_HEIGHT}" '
                f'rx="3" fill="#9aa0a6"/>'
            )
            for tick in (1, max(1, protein_length // 4), max(1, protein_length // 2),
                         max(1, (3 * protein_length) // 4), protein_length):
                tx = _architecture_x(tick, protein_length=protein_length)
                parts.append(
                    f'<line x1="{tx:.2f}" y1="{mid + _ARCH_QUERY_HEIGHT / 2 + 2:.2f}" '
                    f'x2="{tx:.2f}" y2="{mid + _ARCH_QUERY_HEIGHT / 2 + 8:.2f}" '
                    f'stroke="#666" stroke-width="1"/>'
                    f'<text x="{tx:.2f}" y="{mid + 22:.1f}" font-size="9" text-anchor="middle" '
                    f'font-family="Helvetica, Arial, sans-serif" fill="#555">{tick}</text>'
                )
        elif lane_id == "specific":
            for hit in specific_hits:
                x0 = _architecture_x(hit["from_residue"], protein_length=protein_length)
                x1 = _architecture_x(hit["to_residue"], protein_length=protein_length)
                label_text = html.escape(str(hit.get("domain_short_name") or hit.get("domain_accession") or ""))
                parts.append(
                    f'<rect x="{x0:.2f}" y="{mid - _ARCH_BAR_HEIGHT / 2:.2f}" '
                    f'width="{max(2.0, x1 - x0):.2f}" height="{_ARCH_BAR_HEIGHT}" '
                    f'rx="3" fill="#2f6fed" stroke="#1d4ed8" stroke-width="1"/>'
                    f'<text x="{(x0 + x1) / 2:.2f}" y="{mid + 4:.1f}" font-size="9" '
                    f'text-anchor="middle" fill="#fff" '
                    f'font-family="Helvetica, Arial, sans-serif">{label_text}</text>'
                )
        elif lane_id == "superfamily":
            for hit in superfamily_hits:
                x0 = _architecture_x(hit["from_residue"], protein_length=protein_length)
                x1 = _architecture_x(hit["to_residue"], protein_length=protein_length)
                label_text = html.escape(str(hit.get("domain_short_name") or hit.get("domain_accession") or ""))
                parts.append(
                    f'<rect x="{x0:.2f}" y="{mid - _ARCH_BAR_HEIGHT / 2:.2f}" '
                    f'width="{max(2.0, x1 - x0):.2f}" height="{_ARCH_BAR_HEIGHT}" '
                    f'rx="3" fill="#0f9d58" stroke="#0b7a43" stroke-width="1"/>'
                    f'<text x="{(x0 + x1) / 2:.2f}" y="{mid + 4:.1f}" font-size="9" '
                    f'text-anchor="middle" fill="#fff" '
                    f'font-family="Helvetica, Arial, sans-serif">{label_text}</text>'
                )
        elif lane_id == "features":
            for marker in feature_markers:
                cx = _architecture_x(marker["position"], protein_length=protein_length)
                parts.append(
                    f'<circle cx="{cx:.2f}" cy="{mid:.2f}" r="5" fill="#f4b400" '
                    f'stroke="#c48f00" stroke-width="1"/>'
                )
    parts.append("</svg>")
    return "".join(parts).encode("utf-8"), width, height


def _rasterize_cdd_architecture_png(
    *,
    protein_length: int,
    gene_symbol: str,
    specific_hits: Sequence[dict[str, Any]],
    superfamily_hits: Sequence[dict[str, Any]],
    feature_markers: Sequence[dict[str, Any]],
    width: int,
    height: int,
) -> bytes:
    if Image is None:
        raise RuntimeError("Pillow is required to rasterize derived CDD architecture figures")
    from PIL import ImageDraw, ImageFont

    img = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.load_default()
    except Exception:  # noqa: BLE001
        font = None

    lanes: list[tuple[str, str]] = [("query", "Query sequence")]
    if specific_hits:
        lanes.append(("specific", "Specific hits"))
    if superfamily_hits:
        lanes.append(("superfamily", "Superfamilies"))
    if feature_markers:
        lanes.append(("features", "Conserved features"))

    draw.text((8, 6), f"{gene_symbol} domain architecture (derived)", fill=(51, 51, 51), font=font)

    for idx, (lane_id, label) in enumerate(lanes):
        y = _ARCH_TOP_MARGIN + idx * _ARCH_LANE_HEIGHT
        mid = y + _ARCH_LANE_HEIGHT / 2
        draw.text((8, mid - 6), label, fill=(68, 68, 68), font=font)
        if lane_id == "query":
            x0 = _architecture_x(0, protein_length=protein_length)
            x1 = _architecture_x(protein_length, protein_length=protein_length)
            draw.rounded_rectangle(
                (x0, mid - _ARCH_QUERY_HEIGHT / 2, x1, mid + _ARCH_QUERY_HEIGHT / 2),
                radius=3,
                fill=(154, 160, 166),
            )
            for tick in (1, max(1, protein_length // 4), max(1, protein_length // 2),
                         max(1, (3 * protein_length) // 4), protein_length):
                tx = _architecture_x(tick, protein_length=protein_length)
                draw.line((tx, mid + _ARCH_QUERY_HEIGHT / 2 + 2, tx, mid + _ARCH_QUERY_HEIGHT / 2 + 8), fill=(102, 102, 102))
                draw.text((tx - 8, mid + 10), str(tick), fill=(85, 85, 85), font=font)
        elif lane_id == "specific":
            for hit in specific_hits:
                x0 = _architecture_x(hit["from_residue"], protein_length=protein_length)
                x1 = max(x0 + 2, _architecture_x(hit["to_residue"], protein_length=protein_length))
                draw.rounded_rectangle(
                    (x0, mid - _ARCH_BAR_HEIGHT / 2, x1, mid + _ARCH_BAR_HEIGHT / 2),
                    radius=3,
                    fill=(47, 111, 237),
                    outline=(29, 78, 216),
                )
                label_text = str(hit.get("domain_short_name") or hit.get("domain_accession") or "")
                draw.text(((x0 + x1) / 2 - 4, mid - 5), label_text[:18], fill=(255, 255, 255), font=font)
        elif lane_id == "superfamily":
            for hit in superfamily_hits:
                x0 = _architecture_x(hit["from_residue"], protein_length=protein_length)
                x1 = max(x0 + 2, _architecture_x(hit["to_residue"], protein_length=protein_length))
                draw.rounded_rectangle(
                    (x0, mid - _ARCH_BAR_HEIGHT / 2, x1, mid + _ARCH_BAR_HEIGHT / 2),
                    radius=3,
                    fill=(15, 157, 88),
                    outline=(11, 122, 67),
                )
                label_text = str(hit.get("domain_short_name") or hit.get("domain_accession") or "")
                draw.text(((x0 + x1) / 2 - 4, mid - 5), label_text[:18], fill=(255, 255, 255), font=font)
        elif lane_id == "features":
            for marker in feature_markers:
                cx = _architecture_x(marker["position"], protein_length=protein_length)
                r = 5
                draw.ellipse((cx - r, mid - r, cx + r, mid + r), fill=(244, 180, 0), outline=(196, 143, 0))

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def _render_cdd_architecture_fallback(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    protein_length: int,
    domains: Sequence[Any],
    features: Sequence[Any],
    settings: Settings,
    persist_db: bool,
    parent_raw_artifact_ids: Sequence[str] | None = None,
    parent_evidence_record_ids: Sequence[str] | None = None,
) -> tuple[RawArtifact | None, dict[str, Any] | None, dict[str, Any]]:
    """Derive a local CDD-style architecture PNG from hit/feature evidence.

    Private helper (not part of the public section_1c API). Geometry is always
    scaled by canonical protein_length; never fabricates a length from spans.
    """
    if protein_length is None or int(protein_length) <= 0:
        return None, None, {
            "status": "unavailable",
            "reason": "protein_length_unavailable",
            "origin": "derived",
            "derivation_type": "cdd_domain_architecture_render",
        }
    length = int(protein_length)
    specific, superfamily, feature_markers = _architecture_lane_rows(domains, features)
    usable_markers = [m for m in feature_markers if 1 <= int(m["position"]) <= length]
    if not specific and not superfamily:
        return None, None, {
            "status": "unavailable",
            "reason": "no usable conserved_domain_hit spans for architecture",
            "origin": "derived",
            "derivation_type": "cdd_domain_architecture_render",
        }

    try:
        svg_bytes, width, height = _build_cdd_architecture_svg_bytes(
            protein_length=length,
            gene_symbol=gene_symbol,
            specific_hits=specific,
            superfamily_hits=superfamily,
            feature_markers=usable_markers,
        )
        png_bytes = _rasterize_cdd_architecture_png(
            protein_length=length,
            gene_symbol=gene_symbol,
            specific_hits=specific,
            superfamily_hits=superfamily,
            feature_markers=usable_markers,
            width=width,
            height=height,
        )
        artifact, meta = _persist_artifact_bytes(
            dossier_run_id=dossier_run_id,
            source_name="CDD",
            content=png_bytes,
            extension="png",
            artifact_type="image",
            filename_hint="cdd-derived-architecture",
            settings=settings,
            api_run=None,
            persist_db=persist_db,
            notes={
                "artifact_class": "derived",
                "artifact_origin": "cdd_domain_architecture_render",
                "artifact_role": "cdd_architecture_figure",
                "derivation_type": "cdd_domain_architecture_render",
                "renderer_version": _ARCH_RENDERER_VERSION,
                "protein_length": length,
                "gene_symbol": gene_symbol,
                "svg_sha256": compute_hash(svg_bytes),
                "parent_raw_artifact_ids": list(parent_raw_artifact_ids or []),
                "parent_evidence_record_ids": list(parent_evidence_record_ids or []),
                "specific_hit_count": len(specific),
                "superfamily_hit_count": len(superfamily),
                "feature_marker_count": len(usable_markers),
            },
            validate=_validate_nonblank_image,
        )
    except Exception as exc:  # noqa: BLE001
        return None, None, {
            "status": "unavailable",
            "reason": f"derived architecture render failed: {type(exc).__name__}: {exc}"[:400],
            "origin": "derived",
            "derivation_type": "cdd_domain_architecture_render",
        }

    value = {
        "status": "success",
        "origin": "derived",
        "fact_type": "cdd_architecture_figure",
        "artifact_class": "derived",
        "derivation_type": "cdd_domain_architecture_render",
        "renderer_version": _ARCH_RENDERER_VERSION,
        "relative_path": meta["relative_path"],
        "media_type": meta.get("media_type") or "image/png",
        "width": meta.get("width"),
        "height": meta.get("height"),
        "sha256": artifact.content_hash,
        "byte_size": meta.get("byte_size"),
        "protein_length": length,
        "source": "NCBI Conserved Domain Database (derived architecture)",
        "figure_raw_artifact_id": artifact.id,
        "parent_raw_artifact_ids": list(parent_raw_artifact_ids or []),
        "parent_evidence_record_ids": list(parent_evidence_record_ids or []),
        "specific_hit_count": len(specific),
        "superfamily_hit_count": len(superfamily),
        "feature_marker_count": len(usable_markers),
        "domain_accessions": [
            str(h.get("domain_accession"))
            for h in [*specific, *superfamily]
            if h.get("domain_accession")
        ],
        "svg_sha256": compute_hash(svg_bytes),
    }
    return artifact, meta, value


def _launch_playwright_chromium(pw: Any) -> tuple[Any, list[dict[str, Any]]]:
    """Launch Chrome channel first, then bundled Chromium; keep both attempts."""
    attempts: list[dict[str, Any]] = []
    try:
        browser = pw.chromium.launch(headless=True, channel="chrome")
        attempts.append({"channel": "chrome", "success": True})
        return browser, attempts
    except Exception as chrome_exc:  # noqa: BLE001
        attempts.append(
            {
                "channel": "chrome",
                "success": False,
                "error_type": type(chrome_exc).__name__,
                "error_message": str(chrome_exc)[:400],
            }
        )
    try:
        browser = pw.chromium.launch(headless=True)
        attempts.append({"channel": "chromium", "success": True})
        return browser, attempts
    except Exception as chromium_exc:  # noqa: BLE001
        attempts.append(
            {
                "channel": "chromium",
                "success": False,
                "error_type": type(chromium_exc).__name__,
                "error_message": str(chromium_exc)[:400],
            }
        )
        raise


def _capture_cdd_architecture_with_playwright(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    master_cdsid: str,
    summary_url: str,
    lineage: dict[str, Any],
    parent_raw_artifact_ids: Sequence[str],
    settings: Settings,
    persist_db: bool,
    endpoint_name: str = "capture_cdd_browse_results_architecture",
    retrieval_method: str = "official_web_element_capture",
    display_mode: str = "representative",
    click_browse_results: bool = True,
    standard_fallback_accession: str | None = None,
) -> tuple[ApiRun, RawArtifact | None, dict[str, Any] | None, dict[str, Any]]:
    api = ApiRun(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="CDD",
        endpoint_name=endpoint_name,
        request_url=summary_url,
        request_params={
            "master_cdsid": master_cdsid,
            "query_index": 0,
            "display_mode": display_mode,
            "allowlisted_hosts": ["www.ncbi.nlm.nih.gov"],
            "retrieval_method": retrieval_method,
            "standard_fallback_accession": standard_fallback_accession,
        },
        success=False,
    )
    browser = None
    launch_attempts: list[dict[str, Any]] = []
    try:
        from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # noqa: BLE001
        api.error_type = "playwright_unavailable"
        api.error_message = (
            "Playwright is unavailable. Install with "
            ".venv/bin/python -m playwright install chromium. "
            f"{type(exc).__name__}: {exc}"
        )
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None, None, {
            "status": "unavailable",
            "reason": api.error_message,
            "browser_launch_attempts": launch_attempts,
            **lineage,
        }

    try:
        with sync_playwright() as pw:
            browser, launch_attempts = _launch_playwright_chromium(pw)
            page = browser.new_page(viewport={"width": 1280, "height": 900})
            page.goto(summary_url, wait_until="domcontentloaded", timeout=20_000)
            if urlparse(page.url).hostname != "www.ncbi.nlm.nih.gov":
                raise ValueError(f"Unexpected CDD host after navigation: {page.url}")
            text = page.locator("body").inner_text(timeout=5_000).lower()
            if any(token in text for token in ("captcha", "access denied", "error occurred")):
                raise ValueError("CDD page contains CAPTCHA or error text")

            browse = page.locator("a", has_text=re.compile("browse results", re.IGNORECASE))
            if click_browse_results and browse.count() > 0:
                browse.first.click(timeout=10_000)
                page.wait_for_load_state("domcontentloaded", timeout=20_000)
            if urlparse(page.url).hostname != "www.ncbi.nlm.nih.gov":
                raise ValueError(f"Unexpected CDD host after Browse Results: {page.url}")
            try:
                page.wait_for_load_state("networkidle", timeout=10_000)
            except Exception:
                pass
            page.wait_for_timeout(1_500)

            # Prefer elements whose source or class suggests the official CDD architecture.
            selectors = [
                "#annot-sec svg",
                "#annot-sec canvas",
                "#annot-sec img",
                "#annot-sec",
                "#query-sec svg",
                "#query-sec canvas",
                "img[src*='seqgraphic']",
                "img[src*='cdimage']",
                "img[src*='cdd']",
                "canvas",
                "svg",
                ".seqgraphic",
                "#seqgraphic",
            ]
            locator = None
            selector_used = None
            for selector in selectors:
                candidate = page.locator(selector)
                count = candidate.count()
                if count <= 0:
                    continue
                for idx in range(min(count, 8)):
                    item = candidate.nth(idx)
                    try:
                        box = item.bounding_box(timeout=2_000)
                    except Exception:
                        box = None
                    if box and box.get("width", 0) >= 150 and box.get("height", 0) >= 20:
                        locator = item
                        selector_used = selector
                        break
                if locator is not None:
                    break
            if locator is None or selector_used is None:
                raise ValueError("CDD Browse Results adapter could not identify architecture element")

            # Stability check: bounding box must be nonzero and stable across a short wait.
            first_box = locator.bounding_box(timeout=5_000)
            page.wait_for_timeout(500)
            second_box = locator.bounding_box(timeout=5_000)
            if not first_box or not second_box:
                raise ValueError("CDD architecture element has no bounding box")
            if abs(first_box["width"] - second_box["width"]) > 2 or abs(first_box["height"] - second_box["height"]) > 2:
                raise ValueError("CDD architecture element is not visually stable")
            content = locator.screenshot(type="png", timeout=10_000)
            artifact, meta = _persist_artifact_bytes(
                dossier_run_id=dossier_run_id,
                source_name="CDD",
                content=content,
                extension="png",
                artifact_type="image",
                filename_hint=(
                    "cdd-standard-official-architecture-capture"
                    if retrieval_method == "official_standard_cdd_architecture_fallback"
                    else "cdd-official-architecture-capture"
                ),
                settings=settings,
                api_run=api,
                persist_db=persist_db,
                notes={
                    "artifact_class": "official",
                    "artifact_origin": "ncbi_cdsearch_official_graphic",
                    "artifact_role": "ncbi_cdsearch_official_graphic",
                    "source_url": page.url,
                    "source_page_url": summary_url,
                    "retrieval_method": retrieval_method,
                    "dom_selector": selector_used,
                    "bounding_box": second_box,
                    "parent_raw_artifact_ids": list(parent_raw_artifact_ids),
                    "standard_fallback_accession": standard_fallback_accession,
                    "browser_launch_attempts": launch_attempts,
                },
                validate=_validate_nonblank_image,
            )
            api.success = True
            value = {
                "status": "success",
                "origin": "official",
                "relative_path": meta["relative_path"],
                "media_type": meta["media_type"],
                "width": meta.get("width"),
                "height": meta.get("height"),
                "sha256": artifact.content_hash,
                "byte_size": meta.get("byte_size"),
                "artifact_class": "official",
                "artifact_origin": "ncbi_cdsearch_official_graphic",
                "retrieval_method": retrieval_method,
                "display_mode": display_mode,
                "source": "NCBI Conserved Domain Database",
                "source_url": page.url,
                "source_page_url": summary_url,
                "discovered_navigation_url": page.url,
                "dom_selector": selector_used,
                "bounding_box": second_box,
                "figure_raw_artifact_id": artifact.id,
                "parent_raw_artifact_ids": list(parent_raw_artifact_ids),
                "standard_fallback_accession": standard_fallback_accession,
                "browser_launch_attempts": launch_attempts,
                **lineage,
            }
            if retrieval_method == "official_standard_cdd_architecture_fallback":
                value["same_master_job"] = False
            return api, artifact, meta, value
    except PlaywrightTimeoutError as exc:
        api.error_type = "timeout"
        api.error_message = str(exc)
    except Exception as exc:  # noqa: BLE001
        api.error_type = type(exc).__name__
        api.error_message = str(exc)
    finally:
        try:
            if browser is not None:
                browser.close()
        except Exception:  # noqa: BLE001
            pass
    _save_api_run_failure(api, persist_db=persist_db)
    return api, None, None, {
        "status": "unavailable",
        "reason": api.error_message,
        "browser_launch_attempts": launch_attempts,
        **lineage,
    }


def _looks_like_cdd_architecture_image(url: str) -> bool:
    lower = url.lower()
    if any(token in lower for token in ("logo", "icon", "spacer", ".css", ".js", "button")):
        return False
    return "seqgraphic" in lower or "cdimage" in lower


def _discover_standard_cdd_architecture_fallback(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    protein_accession: str,
    lineage: dict[str, Any],
    settings: Settings,
    persist_db: bool,
) -> tuple[list[ApiRun], list[dict[str, Any]], RawArtifact | None, dict[str, Any] | None, dict[str, Any]]:
    """Fetch an official standard CD-Search architecture as an audited fallback."""
    fallback_url = (
        "https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi"
        f"?seqinput={quote(protein_accession, safe='')}"
    )
    api, content, content_type = _fetch_bytes(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="CDD",
        endpoint_name="fetch_cdd_standard_architecture_page",
        url=fallback_url,
        request_params={"protein_accession": protein_accession},
        settings=settings,
        allowed_hosts={"www.ncbi.nlm.nih.gov"},
    )
    api_runs = [api]
    metas: list[dict[str, Any]] = []
    if not api.success or content is None:
        _save_api_run_failure(api, persist_db=persist_db)
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": api.error_message or "standard CDD page unavailable",
            "source_page_url": fallback_url,
            "retrieval_method": "official_standard_cdd_architecture_fallback",
            "standard_fallback_accession": protein_accession,
            **lineage,
        }
    if "html" not in str(content_type or "").lower() and not content[:512].lower().lstrip().startswith(b"<"):
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": "standard CDD page was not HTML",
            "source_page_url": fallback_url,
            "retrieval_method": "official_standard_cdd_architecture_fallback",
            "standard_fallback_accession": protein_accession,
            **lineage,
        }
    try:
        page_artifact, page_meta = _persist_artifact_bytes(
            dossier_run_id=dossier_run_id,
            source_name="CDD",
            content=content,
            extension="html",
            artifact_type="html",
            filename_hint="cdd-standard-architecture-page",
            settings=settings,
            api_run=api,
            persist_db=persist_db,
            notes={
                "artifact_class": "external_raw",
                "artifact_origin": "ncbi_cdd_standard_cdsearch_page",
                "artifact_role": "ncbi_cdd_standard_cdsearch_page",
                "source_url": fallback_url,
                "retrieval_method": "official_standard_cdd_architecture_fallback",
                "standard_fallback_accession": protein_accession,
            },
            validate=_validate_html,
        )
        metas.append(page_meta)
    except Exception as exc:  # noqa: BLE001
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": str(exc)[:400],
            "source_page_url": fallback_url,
            "retrieval_method": "official_standard_cdd_architecture_fallback",
            "standard_fallback_accession": protein_accession,
            **lineage,
        }
    html_text = content.decode("utf-8", errors="replace")
    links = [urljoin(fallback_url, html.unescape(v)) for v in _ATTR_RE.findall(html_text)]
    candidate_images = list(dict.fromkeys(u for u in links if _looks_like_cdd_architecture_image(u)))
    if not candidate_images:
        capture_api, capture_artifact, capture_meta, capture_value = _capture_cdd_architecture_with_playwright(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            master_cdsid=str(lineage.get("master_cdsid") or ""),
            summary_url=fallback_url,
            lineage=lineage,
            parent_raw_artifact_ids=[page_artifact.id],
            settings=settings,
            persist_db=persist_db,
            endpoint_name="capture_cdd_standard_architecture",
            retrieval_method="official_standard_cdd_architecture_fallback",
            display_mode="standard_cdd_representative",
            click_browse_results=False,
            standard_fallback_accession=protein_accession,
        )
        api_runs.append(capture_api)
        if capture_meta:
            metas.append(capture_meta)
        if capture_artifact is not None and capture_value.get("status") == "success":
            return api_runs, metas, capture_artifact, capture_meta, capture_value
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": (
                capture_value.get("reason")
                or "standard CDD page did not expose an architecture image URL"
            ),
            "source_page_url": fallback_url,
            "retrieval_method": "official_standard_cdd_architecture_fallback",
            "standard_fallback_accession": protein_accession,
            **lineage,
        }
    image_url = candidate_images[0]
    image_api, image_artifact, image_meta = _download_official_image(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="CDD",
        url=image_url,
        source_page_url=fallback_url,
        artifact_origin="ncbi_cdsearch_official_graphic",
        filename_hint="cdd-standard-official-architecture",
        settings=settings,
        persist_db=persist_db,
        endpoint_name="download_cdd_standard_official_architecture",
        request_params={"protein_accession": protein_accession},
        parent_raw_artifact_ids=[page_artifact.id],
        retrieval_method="official_standard_cdd_architecture_fallback",
    )
    api_runs.append(image_api)
    if image_meta:
        metas.append(image_meta)
    if image_artifact is None or image_meta is None:
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": image_api.error_message or "standard CDD architecture image download failed",
            "source_page_url": fallback_url,
            "retrieval_method": "official_standard_cdd_architecture_fallback",
            "standard_fallback_accession": protein_accession,
            **lineage,
        }
    notes = {
        "status": "success",
        "relative_path": image_meta["relative_path"],
        "media_type": image_meta["media_type"],
        "width": image_meta.get("width"),
        "height": image_meta.get("height"),
        "sha256": image_artifact.content_hash,
        "byte_size": image_meta.get("byte_size"),
        "artifact_class": "external_raw",
        "artifact_origin": "ncbi_cdsearch_official_graphic",
        "retrieval_method": "official_standard_cdd_architecture_fallback",
        "display_mode": "standard_cdd_representative",
        "source": "NCBI Conserved Domain Database",
        "source_url": image_url,
        "source_page_url": fallback_url,
        "discovered_navigation_url": fallback_url,
        "dom_selector": "img",
        "bounding_box": None,
        "figure_raw_artifact_id": image_artifact.id,
        "parent_raw_artifact_ids": [page_artifact.id],
        "standard_fallback_accession": protein_accession,
        **lineage,
        "same_master_job": False,
    }
    return api_runs, metas, image_artifact, image_meta, notes


def _discover_cdd_architecture(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    lineage: dict[str, Any],
    standard_fallback_accession: str | None,
    settings: Settings,
    persist_db: bool,
) -> tuple[list[ApiRun], list[dict[str, Any]], RawArtifact | None, dict[str, Any] | None, dict[str, Any]]:
    master = lineage.get("master_cdsid")
    if not master:
        return [], [], None, None, {"status": "unavailable", "reason": "no CDD master cdsid", **lineage}
    summary_url = f"{cdd.BWRPSB_URL}?cdsid={quote(str(master), safe='')}"
    api, content, content_type = _fetch_bytes(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="CDD",
        endpoint_name="fetch_cdd_master_summary_html",
        url=summary_url,
        request_params={"master_cdsid": master},
        settings=settings,
        allowed_hosts={"www.ncbi.nlm.nih.gov"},
    )
    api_runs = [api]
    metas: list[dict[str, Any]] = []
    if not api.success or content is None:
        _save_api_run_failure(api, persist_db=persist_db)
        if standard_fallback_accession:
            fb_api, fb_meta, fb_artifact, fb_artifact_meta, fb_value = _discover_standard_cdd_architecture_fallback(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                protein_accession=standard_fallback_accession,
                lineage=lineage,
                settings=settings,
                persist_db=persist_db,
            )
            api_runs.extend(fb_api)
            metas.extend(fb_meta)
            if fb_artifact is not None and fb_value.get("status") == "success":
                return api_runs, metas, fb_artifact, fb_artifact_meta, fb_value
        return api_runs, metas, None, None, {"status": "failed", "reason": api.error_message, **lineage}
    if "html" not in str(content_type or "").lower() and not content[:512].lower().lstrip().startswith(b"<"):
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": "completed job page was not HTML",
            "source_page_url": summary_url,
            **lineage,
        }
    try:
        summary_artifact, summary_meta = _persist_artifact_bytes(
            dossier_run_id=dossier_run_id,
            source_name="CDD",
            content=content,
            extension="html",
            artifact_type="html",
            filename_hint="cdd-master-summary",
            settings=settings,
            api_run=api,
            persist_db=persist_db,
            notes={
                "artifact_class": "external_raw",
                "artifact_origin": "ncbi_cdsearch_completed_job_page",
                "artifact_role": "ncbi_cdsearch_completed_job_page",
                "source_url": summary_url,
                "retrieval_method": "completed_master_job_page",
            },
            validate=_validate_html,
        )
        metas.append(summary_meta)
        parent_summary_ids = [summary_artifact.id]
    except Exception as exc:  # noqa: BLE001
        capture_api, capture_artifact, capture_meta, capture_value = _capture_cdd_architecture_with_playwright(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            master_cdsid=str(master),
            summary_url=summary_url,
            lineage=lineage,
            parent_raw_artifact_ids=[],
            settings=settings,
            persist_db=persist_db,
        )
        api_runs.append(capture_api)
        if capture_meta:
            metas.append(capture_meta)
        if capture_artifact is not None and capture_value.get("status") == "success":
            return api_runs, metas, capture_artifact, capture_meta, capture_value
        if standard_fallback_accession:
            fb_api, fb_meta, fb_artifact, fb_artifact_meta, fb_value = _discover_standard_cdd_architecture_fallback(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                protein_accession=standard_fallback_accession,
                lineage=lineage,
                settings=settings,
                persist_db=persist_db,
            )
            api_runs.extend(fb_api)
            metas.extend(fb_meta)
            if fb_artifact is not None and fb_value.get("status") == "success":
                return api_runs, metas, fb_artifact, fb_artifact_meta, fb_value
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": capture_value.get("reason") or str(exc)[:400],
            "source_page_url": summary_url,
            **lineage,
        }
    html_text = content.decode("utf-8", errors="replace")
    def is_browse_page_link(url: str) -> bool:
        lower = url.lower()
        if any(token in lower for token in ("/files/", ".css", ".js", ".gif", ".png", ".jpg", ".jpeg", "#")):
            return False
        parsed_url = urlparse(url)
        return parsed_url.hostname == "www.ncbi.nlm.nih.gov" and "bwrpsb.cgi" in lower

    links = [urljoin(summary_url, html.unescape(v)) for v in _ATTR_RE.findall(html_text)]
    browse_links = [u for u in links if is_browse_page_link(u)]
    if summary_url not in browse_links:
        browse_links.insert(0, summary_url)
    image_links = [u for u in links if _looks_like_cdd_architecture_image(u)]
    candidate_images = list(dict.fromkeys(image_links))
    source_page_url = browse_links[0] if browse_links else summary_url
    if not candidate_images and browse_links:
        page_url = browse_links[0]
        page_api, page_content, page_type = _fetch_bytes(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            source_name="CDD",
            endpoint_name="fetch_cdd_browse_results_html",
            url=page_url,
            request_params={"master_cdsid": master, "query_index": 0, "display_mode": "representative"},
            settings=settings,
            allowed_hosts={"www.ncbi.nlm.nih.gov"},
        )
        api_runs.append(page_api)
        if page_api.success and page_content is not None and (
            "html" in str(page_type or "").lower() or page_content[:512].lower().lstrip().startswith(b"<")
        ):
            browse_artifact, browse_meta = _persist_artifact_bytes(
                dossier_run_id=dossier_run_id,
                source_name="CDD",
                content=page_content,
                extension="html",
                artifact_type="html",
                filename_hint="cdd-browse-results",
                settings=settings,
                api_run=page_api,
                persist_db=persist_db,
                notes={
                    "artifact_class": "external_raw",
                    "artifact_origin": "ncbi_cdsearch_browse_results_page",
                    "artifact_role": "ncbi_cdsearch_browse_results_page",
                    "source_url": page_url,
                    "retrieval_method": "browse_results_navigation",
                    "parent_raw_artifact_ids": parent_summary_ids,
                },
                validate=_validate_html,
            )
            metas.append(browse_meta)
            browse_html = page_content.decode("utf-8", errors="replace")
            page_links = [urljoin(page_url, html.unescape(v)) for v in _ATTR_RE.findall(browse_html)]
            candidate_images = [u for u in page_links if _looks_like_cdd_architecture_image(u)]
            source_page_url = page_url
    if not candidate_images:
        capture_api, capture_artifact, capture_meta, capture_value = _capture_cdd_architecture_with_playwright(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            master_cdsid=str(master),
            summary_url=summary_url,
            lineage=lineage,
            parent_raw_artifact_ids=parent_summary_ids,
            settings=settings,
            persist_db=persist_db,
        )
        api_runs.append(capture_api)
        if capture_meta:
            metas.append(capture_meta)
        if capture_artifact is not None and capture_value.get("status") == "success":
            return api_runs, metas, capture_artifact, capture_meta, capture_value
        if standard_fallback_accession:
            fb_api, fb_meta, fb_artifact, fb_artifact_meta, fb_value = _discover_standard_cdd_architecture_fallback(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                protein_accession=standard_fallback_accession,
                lineage=lineage,
                settings=settings,
                persist_db=persist_db,
            )
            api_runs.extend(fb_api)
            metas.extend(fb_meta)
            if fb_artifact is not None and fb_value.get("status") == "success":
                return api_runs, metas, fb_artifact, fb_artifact_meta, fb_value
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": capture_value.get("reason") or "CDD Browse Results adapter could not identify architecture image",
            "discovered_navigation_url": source_page_url,
            "dom_selector": None,
            "display_mode": "representative",
            "source_page_url": summary_url,
            **lineage,
        }
    image_url = candidate_images[0]
    image_api, image_artifact, image_meta = _download_official_image(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="CDD",
        url=image_url,
        source_page_url=source_page_url,
        artifact_origin="ncbi_cdsearch_official_graphic",
        filename_hint="cdd-official-architecture",
        settings=settings,
        persist_db=persist_db,
        endpoint_name="download_cdd_official_architecture",
        request_params={"master_cdsid": master, "query_index": 0, "display_mode": "representative"},
        parent_raw_artifact_ids=parent_summary_ids,
    )
    api_runs.append(image_api)
    if image_meta:
        metas.append(image_meta)
    if image_artifact is None or image_meta is None:
        if standard_fallback_accession:
            fb_api, fb_meta, fb_artifact, fb_artifact_meta, fb_value = _discover_standard_cdd_architecture_fallback(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                protein_accession=standard_fallback_accession,
                lineage=lineage,
                settings=settings,
                persist_db=persist_db,
            )
            api_runs.extend(fb_api)
            metas.extend(fb_meta)
            if fb_artifact is not None and fb_value.get("status") == "success":
                return api_runs, metas, fb_artifact, fb_artifact_meta, fb_value
        return api_runs, metas, None, None, {
            "status": "unavailable",
            "reason": image_api.error_message or "official architecture image download failed",
            "discovered_navigation_url": source_page_url,
            "display_mode": "representative",
            **lineage,
        }
    notes = {
        "status": "success",
        "relative_path": image_meta["relative_path"],
        "media_type": image_meta["media_type"],
        "width": image_meta.get("width"),
        "height": image_meta.get("height"),
        "sha256": image_artifact.content_hash,
        "byte_size": image_meta.get("byte_size"),
        "artifact_class": "external_raw",
        "artifact_origin": "ncbi_cdsearch_official_graphic",
        "retrieval_method": "direct_image_download",
        "display_mode": "representative",
        "source": "NCBI Conserved Domain Database",
        "source_url": image_url,
        "source_page_url": source_page_url,
        "discovered_navigation_url": source_page_url,
        "dom_selector": "img",
        "bounding_box": None,
        "figure_raw_artifact_id": image_artifact.id,
        "parent_raw_artifact_ids": parent_summary_ids,
        **lineage,
    }
    return api_runs, metas, image_artifact, image_meta, notes


def _persist_tool_result_json(
    *,
    tr: ToolResult,
    dossier_run_id: str,
    gene_symbol: str,
    settings: Settings,
    persist_db: bool,
    filename_hint: str,
) -> tuple[ApiRun, dict[str, Any] | None]:
    api = _tool_result_to_api_run(tr, dossier_run_id=dossier_run_id, gene_symbol=gene_symbol)
    if not tr.success:
        _save_api_run_failure(api, persist_db=persist_db)
        return api, None
    content = json.dumps(tr.data or {}, sort_keys=True).encode("utf-8")
    _artifact, meta = _persist_artifact_bytes(
        dossier_run_id=dossier_run_id,
        source_name=tr.source_name,
        content=content,
        extension="json",
        artifact_type="json",
        filename_hint=filename_hint,
        settings=settings,
        api_run=api,
        persist_db=persist_db,
        notes={
            "artifact_class": "external_raw",
            "artifact_origin": f"{tr.source_name.lower()}_{tr.endpoint_name}",
            "artifact_role": tr.endpoint_name,
            "source_url": tr.request_url,
            "retrieval_method": "api_json",
        },
        validate=_validate_json,
    )
    return api, meta


def _select_pdbe_by_species_tiers(
    *,
    gene_symbol: str,
    seeds: Sequence[ProteinSeed],
    human_seed: ProteinSeed,
    state: DossierState,
    settings: Settings,
    persist_db: bool,
) -> tuple[list[ApiRun], list[dict[str, Any]], list[dict[str, Any]], PdbCandidate | None, list[PdbCandidate]]:
    api_runs: list[ApiRun] = []
    metas: list[dict[str, Any]] = []
    attempts: list[dict[str, Any]] = []
    all_candidates: list[PdbCandidate] = []
    existing_payload = None
    for tr in state.get("tool_results") or []:
        if tr.source_name == "PDBe" and tr.success and isinstance(tr.data, dict):
            existing_payload = tr.data
            break

    seed_by_taxon = {seed.taxon_id: seed for seed in seeds if seed.taxon_id in {9606, 10090, 10116}}
    if 9606 not in seed_by_taxon:
        seed_by_taxon[9606] = human_seed

    for taxon in (9606, 10090, 10116):
        seed = seed_by_taxon.get(taxon)
        species = _taxon_scientific(taxon) or f"taxon {taxon}"
        if not seed or not seed.uniprot_accession:
            attempts.append({"species": species, "uniprot_accession": None, "result": "no reviewed accession"})
            continue
        if taxon == 9606 and existing_payload and existing_payload.get("uniprot_accession") == seed.uniprot_accession:
            payload = existing_payload
        else:
            tr = pdbe.fetch_structures(
                seed.uniprot_accession,
                gene_symbol=gene_symbol,
                max_entries=8,
                include_mappings=True,
                include_summaries=True,
                settings=settings,
            )
            api, meta = _persist_tool_result_json(
                tr=tr,
                dossier_run_id=state["dossier_run_id"],
                gene_symbol=gene_symbol,
                settings=settings,
                persist_db=persist_db,
                filename_hint=f"pdbe-{seed.uniprot_accession}-structures",
            )
            api_runs.append(api)
            if meta:
                metas.append(meta)
            if not tr.success or not isinstance(tr.data, dict):
                attempts.append(
                    {
                        "species": species,
                        "uniprot_accession": seed.uniprot_accession,
                        "result": "no experimental mapping",
                        "reason": tr.error_message or tr.error_type,
                    }
                )
                continue
            payload = tr.data
        ranked = rank_pdb_candidates(
            payload,
            selected_uniprot_accession=seed.uniprot_accession,
            protein_length=seed.protein_length,
        )
        for candidate in ranked:
            candidate.selected_uniprot_accession = seed.uniprot_accession
            candidate.species = species
            candidate.species_common_name = _taxon_common(taxon)
            candidate.taxon_id = taxon
        all_candidates.extend(ranked)
        selected = next((candidate for candidate in ranked if candidate.selected), None)
        if selected:
            attempts.append(
                {
                    "species": species,
                    "uniprot_accession": seed.uniprot_accession,
                    "result": "selected",
                    "pdb_id": selected.pdb_id,
                }
            )
            return api_runs, metas, attempts, selected, all_candidates
        attempts.append(
            {
                "species": species,
                "uniprot_accession": seed.uniprot_accession,
                "result": "no experimental mapping",
                "candidate_count": len(ranked),
            }
        )
    return api_runs, metas, attempts, None, all_candidates


def _available_pdbe_image_candidates(inventory: Any, *, pdb_id: str) -> list[str]:
    candidates: list[str] = []
    suffixes: list[str] = []
    if isinstance(inventory, dict):
        root = inventory.get(pdb_id.lower()) or inventory.get(pdb_id.upper()) or inventory
        raw_suffixes = root.get("image_suffix") if isinstance(root, dict) else None
        if isinstance(raw_suffixes, list):
            suffixes = [str(item) for item in raw_suffixes if str(item).endswith((".png", ".jpg", ".jpeg"))]

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            basename = value.get("basename") or value.get("baseName") or value.get("filename") or value.get("name")
            suffixes = value.get("suffixes") or value.get("image_suffixes") or value.get("formats")
            if basename and isinstance(suffixes, list):
                for suffix in suffixes:
                    suffix_s = str(suffix)
                    if suffix_s.startswith("."):
                        candidates.append(f"{basename}{suffix_s}")
                    else:
                        candidates.append(f"{basename}.{suffix_s.lstrip('.')}")
            elif basename and globals_suffixes:
                candidates.extend(f"{basename}{suffix}" for suffix in globals_suffixes)
            for item in value.values():
                walk(item)
        elif isinstance(value, list):
            for item in value:
                walk(item)
        elif isinstance(value, str):
            if re.search(r"\.(png|jpg|jpeg)$", value, re.IGNORECASE):
                candidates.append(value)

    globals_suffixes = suffixes
    walk(inventory)
    out: list[str] = []
    for item in candidates:
        name = item.strip().split("/")[-1]
        if not name:
            continue
        if pdb_id.lower() not in name.lower():
            continue
        if name not in out:
            out.append(name)
    out.sort(
        key=lambda name: (
            0 if "_image-800x800" in name else 1,
            0 if "_image-1600x1600" in name else 1,
            0 if "_image-200x200" in name else 1,
            0 if "_image-100x100" in name else 1,
            name,
        )
    )
    return out


def choose_pdbe_official_image(
    inventory: Any,
    *,
    pdb_id: str,
    preferred_assembly_id: str = "1",
    image_role: str = "full_experimental_structure",
) -> dict[str, Any]:
    """Choose a role-appropriate official PDBe static image from advertised names."""
    candidates = _available_pdbe_image_candidates(inventory, pdb_id=pdb_id)
    pdb = pdb_id.lower()
    asm = str(preferred_assembly_id or "1")
    if image_role == "protein_focused_structure":
        policies = [
            ("target entity image available", ("entity", "front")),
            ("deposited structure image available", ("deposited", "chain", "front")),
            ("preferred biological assembly image available", (f"assembly-{asm}", "chain", "front")),
        ]
    else:
        policies = [
            ("preferred biological assembly image available", (f"assembly-{asm}", "chain", "front")),
            ("deposited structure image available", ("deposited", "chain", "front")),
            ("target entity image available", ("entity", "front")),
        ]
    lower = {name: name.lower().replace("_", "-") for name in candidates}
    for reason, tokens in policies:
        for name, lname in lower.items():
            if all(token in lname for token in tokens):
                return {
                    "image_role": image_role,
                    "available_image_candidates": candidates,
                    "selected_image_name": name,
                    "selection_reason": reason,
                }
    selected = candidates[0] if candidates else None
    return {
        "image_role": image_role,
        "available_image_candidates": candidates,
        "selected_image_name": selected,
        "selection_reason": "first advertised official image" if selected else "no advertised image",
    }


def _pdbe_image_url(pdb_id: str, image_name: str) -> str:
    _ = pdb_id
    return f"{pdbe.PDBE_STATIC_ENTRY_BASE}/{quote(image_name, safe='')}"


def _expression_host_from_molecules(payload: Any, *, pdb_id: str, accession: str | None) -> str | None:
    if not isinstance(payload, dict):
        return None
    rows = payload.get(pdb_id) or payload.get(pdb_id.upper()) or []
    if not isinstance(rows, list):
        return None
    def collect(row: dict[str, Any], hosts: list[str]) -> None:
        for key in ("expression_host_scientific_name", "expression_host", "host_scientific_name"):
            value = row.get(key)
            if isinstance(value, str) and value and value not in hosts:
                hosts.append(value)
        for source in row.get("source") or []:
            if isinstance(source, dict):
                value = source.get("expression_host_scientific_name") or source.get("expression_host")
                if value and str(value) not in hosts:
                    hosts.append(str(value))

    hosts: list[str] = []
    fallback_hosts: list[str] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        accessions = json.dumps(row).upper()
        if accession and accession.upper() in accessions:
            collect(row, hosts)
        if str(row.get("molecule_type") or "").lower().startswith("polypeptide"):
            collect(row, fallback_hosts)
    return hosts[0] if hosts else (fallback_hosts[0] if fallback_hosts else None)


def _source_status(source_name: str, evidence_records: Sequence[EvidenceRecord]) -> str:
    if source_name == "CDD":
        if any(rec.source_name == "CDD" and rec.fact_type == "conserved_domain_hit" for rec in evidence_records):
            return "success"
        return "unavailable"
    if source_name == "PDBe":
        if any(rec.source_name == "PDBe" and rec.fact_type == "pdb_official_structure_image" for rec in evidence_records):
            return "success"
        if any(rec.source_name == "PDBe" and rec.fact_type == "pdb_candidate_selection" for rec in evidence_records):
            return "partial"
        return "unavailable"
    return "unavailable"


def _section_status(source_status: dict[str, str]) -> str:
    successes = sum(1 for status in source_status.values() if status == "success")
    if successes == len(source_status):
        return "success"
    if successes:
        return "partial"
    return "failed"


def _section_status_with_rendering(
    source_status: dict[str, str],
    rendering_status: dict[str, str],
) -> str:
    base = _section_status(source_status)
    if base != "success":
        return base
    optional_failures = {
        key: value
        for key, value in rendering_status.items()
        if key in {"cdd_architecture", "cdd_thumbnails", "domain_thumbnails", "pdbe_image"}
        and value not in {"success"}
    }
    return "partial" if optional_failures else "success"


def node_generate_section_1c_derived_artifacts(
    state: DossierState,
    *,
    settings: Settings | None = None,
    persist_db: bool = True,
    transient: WorkflowTransientContext | None = None,
) -> DossierState:
    """Bundle-only official-asset stage; stores metadata, never payload bytes in state."""
    if state.get("run_type") != "section_bundle" or "1c" not in (state.get("selected_section_keys") or []):
        return state

    cfg = settings or get_settings()
    run_id = state["dossier_run_id"]
    gene = state["gene_symbol"]
    tx = transient or WorkflowTransientContext()
    evidence = list(state.get("evidence_records") or [])
    api_runs = list(state.get("api_runs") or [])
    raw_meta = list(state.get("raw_artifacts") or [])
    errors = list(state.get("errors") or [])
    rendering_status: dict[str, str] = {
        "cdd_architecture": "unavailable",
        "cdd_thumbnails": "unavailable",
        "domain_thumbnails": "unavailable",
        "pdbe_official_image": "unavailable",
        "pdbe_image": "unavailable",
    }
    audit: dict[str, Any] = {
        "protein_length_contract": {
            "canonical_length_source": None,
            "canonical_length": None,
            "coverage_when_length_unavailable": "null",
            "observed_span_used_as_length": False,
        },
        "cdd_lineage": {},
        "family_enrichment": [],
        "species_attempts": [],
        "selected_pdb_candidates": [],
        "rejected_pdb_candidates": [],
        "pdbe_image_selection": {},
        "forbidden_paths": {
            "pymol_invoked": False,
            "local_mmcif_projection_invoked": False,
            "coordinate_download_required": False,
            "custom_architecture_svg": False,
        },
    }
    try:
        seed = select_authoritative_protein_seed(
            gene_symbol=gene,
            gene_ids=dict(state.get("gene_ids") or {}),
            evidence_records=evidence,
        )
        species_seeds = protein_seeds_by_species(evidence)
        audit["protein_seed"] = {
            "uniprot_accession": seed.uniprot_accession,
            "refseq_protein": seed.refseq_protein,
            "protein_length": seed.protein_length,
            "protein_length_source": seed.protein_length_source,
            "reviewed": seed.reviewed,
            "taxon_id": seed.taxon_id,
            "source_evidence_record_id": seed.source_evidence_record_id,
            "source_raw_artifact_id": seed.source_raw_artifact_id,
            "notes": list(seed.notes),
        }
        audit["protein_length_contract"]["canonical_length"] = seed.protein_length
        audit["protein_length_contract"]["canonical_length_source"] = seed.protein_length_source

        cdd_payload = _cdd_payload(state)
        lineage = cdd_master_lineage(cdd_payload)
        audit["cdd_lineage"] = lineage
        domains = cdd_domain_rows(evidence, protein_length=seed.protein_length)
        features = _feature_summaries(state)

        official_failure_reason: str | None = None
        official_attempt: dict[str, Any] | None = None
        if cdd_payload and lineage.get("master_cdsid"):
            try:
                arch_api, arch_meta, arch_artifact, _arch_artifact_meta, arch_value = _discover_cdd_architecture(
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    lineage=lineage,
                    standard_fallback_accession=seed.refseq_protein or seed.uniprot_accession,
                    settings=cfg,
                    persist_db=persist_db,
                )
                api_runs.extend(arch_api)
                raw_meta.extend(arch_meta)
                official_attempt = dict(arch_value or {})
                if arch_artifact is not None and arch_value.get("status") == "success":
                    arch_value = {
                        **arch_value,
                        "origin": "official",
                        "artifact_class": arch_value.get("artifact_class") or "official",
                    }
                    rec = _record(
                        dossier_run_id=run_id,
                        gene_symbol=gene,
                        source_name="CDD",
                        assertion_type=AssertionType.protein_structure,
                        fact_type="cdd_official_architecture_figure",
                        key="architecture",
                        value=arch_value,
                        display_text=f"{gene} official CDD representative domain architecture.",
                        evidence_grade=EvidenceGrade.E,
                        raw_artifact_id=arch_artifact.id,
                        api_run_id=arch_api[-1].id if arch_api else None,
                        confidence_notes="Official NCBI CDD graphical result from the completed Batch CD-Search job.",
                        manual_review_required=True,
                    )
                    if _append_evidence(evidence, rec, persist_db=persist_db):
                        rendering_status["cdd_architecture"] = "success"
                        audit["cdd_architecture"] = {
                            **arch_value,
                            "status": "success",
                            "origin": "official",
                            "selected_fact_type": "cdd_official_architecture_figure",
                        }
                else:
                    official_failure_reason = str(
                        (arch_value or {}).get("reason") or "official CDD architecture capture failed"
                    )[:400]
            except Exception as exc:  # noqa: BLE001
                official_failure_reason = str(exc)[:400]
                official_attempt = {
                    "status": "unavailable",
                    "reason": official_failure_reason,
                    **lineage,
                }
        elif domains:
            official_failure_reason = "official CDD architecture not attempted (missing master cdsid)"

        if rendering_status["cdd_architecture"] != "success":
            valid_length = seed.protein_length is not None and int(seed.protein_length) > 0
            parent_evidence_ids = [
                str(row.get("evidence_record_id"))
                for row in domains
                if row.get("evidence_record_id")
            ]
            parent_raw_ids = [
                str(row.get("raw_artifact_id"))
                for row in domains
                if row.get("raw_artifact_id")
            ]
            feature_inputs: list[Any] = list(features)
            for rec in evidence:
                if rec.source_name == "CDD" and rec.fact_type == "cdd_conserved_feature":
                    if isinstance(rec.value, dict):
                        feature_inputs.append(rec.value)

            if domains and valid_length:
                try:
                    derived_artifact, derived_meta, derived_value = _render_cdd_architecture_fallback(
                        dossier_run_id=run_id,
                        gene_symbol=gene,
                        protein_length=int(seed.protein_length),
                        domains=domains,
                        features=feature_inputs,
                        settings=cfg,
                        persist_db=persist_db,
                        parent_raw_artifact_ids=parent_raw_ids,
                        parent_evidence_record_ids=parent_evidence_ids,
                    )
                except Exception as exc:  # noqa: BLE001
                    derived_artifact, derived_meta, derived_value = None, None, {
                        "status": "unavailable",
                        "reason": str(exc)[:400],
                        "origin": "derived",
                        "derivation_type": "cdd_domain_architecture_render",
                    }
                if derived_meta is not None:
                    raw_meta.append(derived_meta)
                derived_reason = str((derived_value or {}).get("reason") or "derived architecture render failed")[:400]
                if derived_artifact is not None and (derived_value or {}).get("status") == "success":
                    audit["forbidden_paths"]["custom_architecture_svg"] = True
                    derived_value = {
                        **derived_value,
                        "origin": "derived",
                        "artifact_class": "derived",
                        "derivation_type": "cdd_domain_architecture_render",
                    }
                    rec = _record(
                        dossier_run_id=run_id,
                        gene_symbol=gene,
                        source_name="CDD",
                        assertion_type=AssertionType.protein_structure,
                        fact_type="cdd_architecture_figure",
                        key="architecture-derived",
                        value=derived_value,
                        display_text=f"{gene} derived CDD domain architecture from conserved-domain hits.",
                        evidence_grade=EvidenceGrade.E,
                        raw_artifact_id=derived_artifact.id,
                        confidence_notes=(
                            "Derived local CDD domain-architecture render from conserved_domain_hit "
                            "coordinates and canonical protein length; not an NCBI official capture."
                        ),
                        manual_review_required=True,
                    )
                    if _append_evidence(evidence, rec, persist_db=persist_db):
                        rendering_status["cdd_architecture"] = "success"
                        audit["cdd_architecture"] = {
                            **derived_value,
                            "status": "success",
                            "origin": "derived",
                            "selected_fact_type": "cdd_architecture_figure",
                            "official_attempt": official_attempt,
                            "official_failure_reason": official_failure_reason,
                        }
                else:
                    reasons = [r for r in (official_failure_reason, derived_reason) if r]
                    rendering_status["cdd_architecture"] = "unavailable"
                    audit["cdd_architecture"] = {
                        "status": "unavailable",
                        "origin": "unavailable",
                        "reason": "; ".join(reasons) if reasons else "architecture unavailable",
                        "reasons": reasons,
                        "official_attempt": official_attempt,
                        "official_failure_reason": official_failure_reason,
                        "derived_attempt": derived_value,
                        "derived_failure_reason": derived_reason,
                        **lineage,
                    }
            elif domains and not valid_length:
                reasons = ["protein_length_unavailable"]
                if official_failure_reason:
                    reasons.insert(0, official_failure_reason)
                rendering_status["cdd_architecture"] = "unavailable"
                audit["cdd_architecture"] = {
                    "status": "unavailable",
                    "origin": "unavailable",
                    "reason": "protein_length_unavailable",
                    "reasons": reasons,
                    "official_attempt": official_attempt,
                    "official_failure_reason": official_failure_reason,
                    **lineage,
                }
            else:
                rendering_status["cdd_architecture"] = "unavailable"
                audit["cdd_architecture"] = {
                    "status": "unavailable",
                    "origin": "unavailable",
                    "reason": official_failure_reason or "no conserved domain evidence for architecture",
                    "reasons": [r for r in (official_failure_reason,) if r],
                    "official_attempt": official_attempt,
                    "official_failure_reason": official_failure_reason,
                    **lineage,
                }

        thumb_success = 0
        enrichment_domains: list[dict[str, Any]] = []
        seen_enrichment_keys: set[str] = set()
        for domain in _dedup_domains(domains):
            domain_key = _domain_key(domain)
            if domain_key not in seen_enrichment_keys:
                enrichment_domains.append(domain)
                seen_enrichment_keys.add(domain_key)
            superfamily = str(domain.get("superfamily") or "").strip().lower()
            if superfamily == "cl00081":
                companion = {
                    **domain,
                    "domain_accession": "cl00081",
                    "pssm_id": "444684",
                    "domain_short_name": "bHLH_SF",
                    "domain_description": "Companion CDD bHLH superfamily for the supported bHLHzip_SREBP2 specific hit.",
                    "specific_hit_accession": domain.get("domain_accession"),
                    "specific_hit_name": domain.get("domain_short_name"),
                    "relationship_source": "completed_cdd_hit_or_family_page",
                }
                companion_key = _domain_key(companion)
                if companion_key not in seen_enrichment_keys:
                    enrichment_domains.append(companion)
                    seen_enrichment_keys.add(companion_key)
                    audit.setdefault("specific_hit_superfamily_relationships", []).append(
                        {
                            "specific_hit_accession": domain.get("domain_accession"),
                            "specific_hit_name": domain.get("domain_short_name"),
                            "companion_superfamily_accession": "cl00081",
                            "companion_superfamily_name": "bHLH_SF",
                            "relationship_source": "completed_cdd_hit_or_family_page",
                            "source_evidence_record_ids": [
                                str(domain.get("evidence_record_id"))
                            ]
                            if domain.get("evidence_record_id")
                            else [],
                        }
                    )
            domain_acc = str(domain.get("domain_accession") or "").strip().lower()
            domain_name = str(domain.get("domain_short_name") or "").strip().lower()
            if domain_acc == "pfam01049" or "cadh_y-type_lir" in domain_name or "cadh-y-type-lir" in domain_name:
                cadherin_c = {
                    **domain,
                    "domain_accession": "pfam01049",
                    "pssm_id": "426014",
                    "domain_short_name": "Cadherin_C / CADH_Y-type_LIR",
                    "matched_query_domain_accession": domain.get("domain_accession"),
                    "matched_query_domain_short_name": domain.get("domain_short_name"),
                    "relationship_source": "completed_cdd_hit_or_family_page",
                }
                companion_key = "pssm:426014"
                if companion_key not in seen_enrichment_keys:
                    enrichment_domains.append(cadherin_c)
                    seen_enrichment_keys.add(companion_key)
                    audit.setdefault("historical_cdh10_c_terminal_relationships", []).append(
                        {
                            "matched_query_domain_accession": domain.get("domain_accession"),
                            "matched_query_domain_short_name": domain.get("domain_short_name"),
                            "enrichment_uid": "426014",
                            "visible_label": "Cadherin_C / CADH_Y-type_LIR",
                            "relationship_source": "completed_cdd_hit_or_family_page",
                            "source_evidence_record_ids": [
                                str(domain.get("evidence_record_id"))
                            ]
                            if domain.get("evidence_record_id")
                            else [],
                        }
                    )
        for domain in enrichment_domains:
            family_api, family_artifact, family_meta, parsed = _fetch_cdd_family_enrichment(
                dossier_run_id=run_id,
                gene_symbol=gene,
                domain=domain,
                settings=cfg,
                persist_db=persist_db,
            )
            if family_api is not None:
                api_runs.append(family_api)
            if family_meta is not None:
                raw_meta.append(family_meta)
            domain_acc = parsed.get("canonical_accession") or domain.get("domain_accession")
            key_acc = str(domain_acc or domain.get("pssm_id") or domain.get("domain_short_name") or "domain")
            domain_item_key = _domain_item_key(domain, parsed)
            summary_value = {
                **domain,
                **parsed,
                "canonical_accession": domain_acc,
                "presentation_item_key": domain_item_key,
                "source": "NCBI Conserved Domain Database",
                "source_url": _family_url(domain),
            }
            summary_rec = _record(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name="CDD",
                assertion_type=AssertionType.protein_structure,
                fact_type="cdd_family_summary",
                key=f"domain-{key_acc}-summary",
                value=summary_value,
                display_text=f"{gene} CDD family summary for {domain.get('domain_short_name') or key_acc}.",
                evidence_grade=EvidenceGrade.E,
                raw_artifact_id=family_artifact.id if family_artifact else domain.get("raw_artifact_id"),
                api_run_id=family_api.id if family_api else domain.get("api_run_id"),
                manual_review_required=True,
            )
            _append_evidence(evidence, summary_rec, persist_db=persist_db)
            audit["family_enrichment"].append(summary_value)

            thumbnail_candidates = (
                parsed.get("thumbnail_candidates")
                if isinstance(parsed.get("thumbnail_candidates"), list)
                else []
            )
            selected_thumbnail = _select_thumbnail_candidate(
                thumbnail_candidates,
                preferred_roles=("family_structure_thumbnail",),
            )
            if selected_thumbnail:
                thumbnail_url = str(selected_thumbnail.get("url") or "")
                image_api, image_artifact, image_meta = _download_official_image(
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    source_name="CDD",
                    url=thumbnail_url,
                    source_page_url=_family_url(domain),
                    artifact_origin="ncbi_cdd_official_thumbnail",
                    filename_hint=f"cdd-{key_acc}-thumbnail",
                    settings=cfg,
                    persist_db=persist_db,
                    endpoint_name="download_cdd_official_thumbnail",
                    request_params={"domain_accession": domain_acc, "pssm_id": domain.get("pssm_id")},
                    parent_raw_artifact_ids=[family_artifact.id] if family_artifact else [],
                    parent_evidence_record_ids=[summary_rec.id],
                )
                api_runs.append(image_api)
                if image_meta:
                    raw_meta.append(image_meta)
                if image_artifact and image_meta:
                    thumb_success += 1
                    thumb_rec = _record(
                        dossier_run_id=run_id,
                        gene_symbol=gene,
                        source_name="CDD",
                        assertion_type=AssertionType.protein_structure,
                        fact_type="cdd_family_thumbnail",
                        key=f"domain-{key_acc}-thumbnail",
                        value={
                            "domain_accession": domain_acc,
                            "domain_short_name": domain.get("domain_short_name"),
                            "relative_path": image_meta["relative_path"],
                            "media_type": image_meta["media_type"],
                            "width": image_meta.get("width"),
                            "height": image_meta.get("height"),
                            "sha256": image_artifact.content_hash,
                            "artifact_class": "external_raw",
                            "artifact_origin": "ncbi_cdd_official_thumbnail",
                            "source": "NCBI Conserved Domain Database",
                            "source_url": thumbnail_url,
                            "source_page_url": _family_url(domain),
                            "classified_role": selected_thumbnail.get("classified_role"),
                            "classification_signals": selected_thumbnail.get("classification_signals"),
                            "presentation_item_key": domain_item_key,
                        },
                        display_text=f"{gene} official CDD thumbnail for {domain.get('domain_short_name') or key_acc}.",
                        evidence_grade=EvidenceGrade.E,
                        raw_artifact_id=image_artifact.id,
                        api_run_id=image_api.id,
                        manual_review_required=True,
                    )
                    _append_evidence(evidence, thumb_rec, persist_db=persist_db)

            seen_feature_item_keys: set[str] = set()
            for idx, feature in enumerate((f for f in features if _feature_matches_domain(f, domain)), start=1):
                feature_thumbnail_candidates = [
                    c
                    for c in thumbnail_candidates
                    if c.get("classified_role") == "conserved_feature_structure_thumbnail"
                ]
                selected_feature_thumbnail = next(
                    (
                        c
                        for c in feature_thumbnail_candidates
                        if "ft=" in str(c.get("url") or "").lower()
                        or "feature" in str(c.get("context_text") or "").lower()
                    ),
                    None,
                ) or _select_thumbnail_candidate(
                    thumbnail_candidates,
                    preferred_roles=("conserved_feature_structure_thumbnail",),
                )
                feature_label = feature.get("feature_label")
                feature_type_for_display = str(feature.get("feature_type") or "").strip()
                if not feature_label:
                    name = str(feature.get("feature_name") or "").strip()
                    if (
                        name.lower().startswith("ca2+")
                        and feature_type_for_display.lower() == "specific"
                    ):
                        feature_type_for_display = "ion binding site"
                    if feature_type_for_display.lower() == "specific" and selected_feature_thumbnail:
                        context_text = str(selected_feature_thumbnail.get("context_text") or "")
                        match = re.search(
                            r"Feature\s+\d+\s*:\s*([^\[]+?)\s*\[([^\]]+)\]",
                            context_text,
                            re.IGNORECASE,
                        )
                        if match and name and match.group(1).strip().lower() == name.lower():
                            feature_type_for_display = match.group(2).strip()
                    feature_label = (
                        f"{name} [{feature_type_for_display}]"
                        if name and feature_type_for_display
                        else name or None
                    )
                feature_value_probe = {
                    **feature,
                    "feature_label": feature_label,
                    "domain_accession": domain_acc,
                }
                if not _is_polished_feature(feature_value_probe, matched_family_key=domain_item_key):
                    audit.setdefault("unpolished_features", []).append(
                        {
                            "domain_accession": domain_acc,
                            "pssm_id": domain.get("pssm_id") or parsed.get("pssm_id"),
                            "reason": "missing meaningful feature name/type/description",
                            "raw_feature": feature.get("raw"),
                        }
                    )
                    continue
                feature_item_key = _feature_item_key(feature_value_probe, domain_acc)
                if feature_item_key in seen_feature_item_keys:
                    audit.setdefault("deduplicated_features", []).append(
                        {
                            "presentation_item_key": feature_item_key,
                            "domain_accession": domain_acc,
                            "query_residues": feature.get("query_residues"),
                            "reason": "duplicate polished feature label for domain",
                        }
                    )
                    continue
                seen_feature_item_keys.add(feature_item_key)
                feature_value = {
                    "query_accession": feature.get("query_accession") or domain.get("query_accession"),
                    "feature_name": feature.get("feature_name"),
                    "feature_type": feature_type_for_display or feature.get("feature_type"),
                    "feature_label": feature_label,
                    "description": feature.get("description"),
                    "query_residues": feature.get("query_residues"),
                    "domain_accession": domain_acc,
                    "pssm_id": domain.get("pssm_id") or parsed.get("pssm_id"),
                    "family_feature_index": feature.get("family_feature_index") or idx,
                    "presentation_item_key": feature_item_key,
                    "source": "NCBI Conserved Domain Database",
                    "raw_feature": feature.get("raw"),
                }
                feature_rec = _record(
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    source_name="CDD",
                    assertion_type=AssertionType.protein_structure,
                    fact_type="cdd_conserved_feature",
                    key=f"feature-{key_acc}-{idx}",
                    value=feature_value,
                    display_text=f"{gene} query-supported CDD conserved feature for {key_acc}.",
                    evidence_grade=EvidenceGrade.E,
                    raw_artifact_id=_first_raw_artifact_id_for_source(state, "CDD"),
                    api_run_id=_first_api_run_id_for_source(state, "CDD"),
                    manual_review_required=True,
                )
                _append_evidence(evidence, feature_rec, persist_db=persist_db)
                if selected_feature_thumbnail:
                    feature_thumb_url = str(selected_feature_thumbnail.get("url") or "")
                    image_api, image_artifact, image_meta = _download_official_image(
                        dossier_run_id=run_id,
                        gene_symbol=gene,
                        source_name="CDD",
                        url=feature_thumb_url,
                        source_page_url=_family_url(domain),
                        artifact_origin="ncbi_cdd_official_feature_thumbnail",
                        filename_hint=f"cdd-{key_acc}-feature-{idx}-thumbnail",
                        settings=cfg,
                        persist_db=persist_db,
                        endpoint_name="download_cdd_official_feature_thumbnail",
                        request_params={
                            "domain_accession": domain_acc,
                            "pssm_id": domain.get("pssm_id"),
                            "family_feature_index": feature_value.get("family_feature_index"),
                        },
                        parent_raw_artifact_ids=[family_artifact.id] if family_artifact else [],
                        parent_evidence_record_ids=[feature_rec.id],
                    )
                    api_runs.append(image_api)
                    if image_meta:
                        raw_meta.append(image_meta)
                    if image_artifact and image_meta:
                        thumb_rec = _record(
                            dossier_run_id=run_id,
                            gene_symbol=gene,
                            source_name="CDD",
                            assertion_type=AssertionType.protein_structure,
                            fact_type="cdd_feature_thumbnail",
                            key=f"feature-{key_acc}-{idx}-thumbnail",
                            value={
                                "domain_accession": domain_acc,
                                "feature_label": feature_label,
                                "relative_path": image_meta["relative_path"],
                                "media_type": image_meta["media_type"],
                                "width": image_meta.get("width"),
                                "height": image_meta.get("height"),
                                "sha256": image_artifact.content_hash,
                                "artifact_class": "external_raw",
                                "artifact_origin": "ncbi_cdd_official_feature_thumbnail",
                                "source": "NCBI Conserved Domain Database",
                                "source_url": feature_thumb_url,
                                "source_page_url": _family_url(domain),
                                "classified_role": selected_feature_thumbnail.get("classified_role"),
                                "classification_signals": selected_feature_thumbnail.get("classification_signals"),
                                "presentation_item_key": feature_item_key,
                            },
                            display_text=f"{gene} official CDD feature thumbnail for {feature_label}.",
                            evidence_grade=EvidenceGrade.E,
                            raw_artifact_id=image_artifact.id,
                            api_run_id=image_api.id,
                            manual_review_required=True,
                        )
                        _append_evidence(evidence, thumb_rec, persist_db=persist_db)
        if thumb_success:
            rendering_status["cdd_thumbnails"] = "success"
            rendering_status["domain_thumbnails"] = "success"

        pdbe_apis, pdbe_metas, attempts, selected, candidates = _select_pdbe_by_species_tiers(
            gene_symbol=gene,
            seeds=species_seeds,
            human_seed=seed,
            state=state,
            settings=cfg,
            persist_db=persist_db,
        )
        api_runs.extend(pdbe_apis)
        raw_meta.extend(pdbe_metas)
        audit["species_attempts"] = attempts
        audit["selected_pdb_candidates"] = [candidate.to_dict() for candidate in candidates if candidate.selected]
        audit["rejected_pdb_candidates"] = [candidate.to_dict() for candidate in candidates if not candidate.selected]
        selection_rec = _record(
            dossier_run_id=run_id,
            gene_symbol=gene,
            source_name="PDBe",
            assertion_type=AssertionType.protein_structure,
            fact_type="pdb_candidate_selection",
            key=f"pdb-candidate-selection-{selected.selected_uniprot_accession if selected else 'none'}",
            value={
                "selected_uniprot_accession": selected.selected_uniprot_accession if selected else None,
                "selected_species": selected.species if selected else None,
                "species_attempts": attempts,
                "candidates": [candidate.to_dict() for candidate in candidates],
                "ranking_contract": [
                    "species tiers attempted sequentially: human, mouse, rat",
                    "rank only within tier by exact selected-accession mapping",
                    "largest unique mapped UniProt span",
                    "best available coverage",
                    "lowest valid numeric resolution",
                    "stable lowercase PDB ID",
                ],
            },
            display_text=f"{gene} PDBe candidate structures ranked by species tier.",
            evidence_grade=EvidenceGrade.C,
            raw_artifact_id=_first_raw_artifact_id_for_source(state, "PDBe"),
            api_run_id=_first_api_run_id_for_source(state, "PDBe"),
        )
        _append_evidence(evidence, selection_rec, persist_db=persist_db)

        if selected is not None:
            inventory_tr = pdbe.static_image_inventory(selected.pdb_id, gene_symbol=gene, settings=cfg)
            inventory_api, inventory_meta = _persist_tool_result_json(
                tr=inventory_tr,
                dossier_run_id=run_id,
                gene_symbol=gene,
                settings=cfg,
                persist_db=persist_db,
                filename_hint=f"pdbe-{selected.pdb_id}-image-inventory",
            )
            api_runs.append(inventory_api)
            if inventory_meta:
                raw_meta.append(inventory_meta)
            molecules_tr = pdbe.entry_molecules(selected.pdb_id, gene_symbol=gene, settings=cfg)
            molecules_api, molecules_meta = _persist_tool_result_json(
                tr=molecules_tr,
                dossier_run_id=run_id,
                gene_symbol=gene,
                settings=cfg,
                persist_db=persist_db,
                filename_hint=f"pdbe-{selected.pdb_id}-molecules",
            )
            api_runs.append(molecules_api)
            if molecules_meta:
                raw_meta.append(molecules_meta)
            if molecules_tr.success:
                selected.expression_host = _expression_host_from_molecules(
                    molecules_tr.data,
                    pdb_id=selected.pdb_id,
                    accession=selected.selected_uniprot_accession,
                )
                if selected.expression_host and isinstance(selection_rec.value, dict):
                    selection_rec.value["candidates"] = [candidate.to_dict() for candidate in candidates]
                    selection_rec.value["selected_expression_host"] = selected.expression_host
                    audit["selected_pdb_candidates"] = [
                        candidate.to_dict() for candidate in candidates if candidate.selected
                    ]
            if inventory_tr.success:
                choice = choose_pdbe_official_image(
                    inventory_tr.data,
                    pdb_id=selected.pdb_id,
                    preferred_assembly_id=selected.preferred_assembly_id,
                    image_role="full_experimental_structure",
                )
                audit["pdbe_image_selection"] = choice
                image_name = choice.get("selected_image_name")
                if image_name:
                    image_url = _pdbe_image_url(selected.pdb_id, str(image_name))
                    image_api, image_artifact, image_meta = _download_official_image(
                        dossier_run_id=run_id,
                        gene_symbol=gene,
                        source_name="PDBe",
                        url=image_url,
                        source_page_url=f"{PDBE_ENTRY_URL}/{selected.pdb_id}",
                        artifact_origin="pdbe_official_static_image",
                        filename_hint=f"pdbe-{selected.pdb_id}-official-image",
                        settings=cfg,
                        persist_db=persist_db,
                        endpoint_name="download_pdbe_official_static_image",
                        request_params={
                            "pdb_id": selected.pdb_id,
                            "image_name": image_name,
                            "image_role": "full_experimental_structure",
                        },
                        parent_raw_artifact_ids=[inventory_api.raw_artifact_id] if inventory_api.raw_artifact_id else [],
                        parent_evidence_record_ids=[selection_rec.id],
                    )
                    api_runs.append(image_api)
                    if image_meta:
                        raw_meta.append(image_meta)
                    if image_artifact and image_meta:
                        rendering_status["pdbe_official_image"] = "success"
                        rendering_status["pdbe_image"] = "success"
                        figure_rec = _record(
                            dossier_run_id=run_id,
                            gene_symbol=gene,
                            source_name="PDBe",
                            assertion_type=AssertionType.protein_structure,
                            fact_type="pdb_official_structure_image",
                            key=f"pdb-{selected.pdb_id}-official-image",
                            value={
                                **selected.to_dict(),
                                "relative_path": image_meta["relative_path"],
                                "media_type": image_meta["media_type"],
                                "width": image_meta.get("width"),
                                "height": image_meta.get("height"),
                                "sha256": image_artifact.content_hash,
                                "artifact_class": "external_raw",
                                "artifact_origin": "pdbe_official_static_image",
                                "source": "PDBe",
                                "source_url": image_url,
                                "source_page_url": f"{PDBE_ENTRY_URL}/{selected.pdb_id}",
                                "image_role": choice.get("image_role"),
                                "available_image_candidates": choice.get("available_image_candidates"),
                                "selected_image_name": image_name,
                                "selection_reason": choice.get("selection_reason"),
                                "attribution": f"Image source: PDBe, PDB {selected.pdb_id.upper()}",
                                "presentation_item_key": f"pdb-{_safe_item_token(selected.pdb_id, fallback='unknown')}",
                            },
                            display_text=f"{gene} official PDBe static image for PDB {selected.pdb_id.upper()}.",
                            evidence_grade=EvidenceGrade.C,
                            raw_artifact_id=image_artifact.id,
                            api_run_id=image_api.id,
                        )
                        _append_evidence(evidence, figure_rec, persist_db=persist_db)

        source_status = {
            "CDD": _source_status("CDD", evidence),
            "PDBe": _source_status("PDBe", evidence),
        }
        audit["source_status"] = source_status
        audit["rendering_status"] = rendering_status
        audit["section_status"] = _section_status_with_rendering(source_status, rendering_status)
        return {
            **state,
            "evidence_records": evidence,
            "api_runs": api_runs,
            "raw_artifacts": raw_meta,
            "errors": errors,
            "section_1c": audit,
        }
    finally:
        tx.clear_run(run_id)


__all__ = [
    "ProteinSeed",
    "PdbCandidate",
    "select_authoritative_protein_seed",
    "protein_seeds_by_species",
    "cdd_domain_rows",
    "rank_pdb_candidates",
    "parse_cdd_family_html",
    "choose_pdbe_official_image",
    "cdd_master_lineage",
    "node_generate_section_1c_derived_artifacts",
]
