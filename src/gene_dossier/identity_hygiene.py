"""Deduplicate stale species-identity EvidenceRecords after query-symbol migration.

When identity records began retaining the dossier query symbol (e.g. ``SREBF2``)
instead of the species symbol (e.g. ``Srebf2``), re-normalization could leave
duplicate logical facts. This module collapses those duplicates without
touching RawArtifacts or ApiRuns.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence

from gene_dossier.models import AssertionType, EvidenceRecord

IDENTITY_SOURCE_NAMES = frozenset({"NCBI Gene", "Ensembl", "UniProt"})


@dataclass(frozen=True)
class IdentityDedupResult:
    """Outcome of an identity-evidence deduplication pass."""

    retained: tuple[EvidenceRecord, ...]
    removed: tuple[EvidenceRecord, ...]


def is_species_identity_record(rec: EvidenceRecord) -> bool:
    """True when ``rec`` is a normalized gene-identity fact from identity sources."""
    if rec.source_name not in IDENTITY_SOURCE_NAMES:
        return False
    assertion = rec.assertion_type
    if assertion == AssertionType.gene_identity:
        return True
    return str(getattr(assertion, "value", assertion)) == "gene_identity"


def normalize_entrez_id(raw: Any) -> str | None:
    """Normalize an Entrez gene ID (strip non-digits / leading zeros)."""
    text = str(raw or "").strip()
    if not text:
        return None
    digits = "".join(ch for ch in text if ch.isdigit())
    if not digits:
        return None
    return str(int(digits))


def normalize_ensembl_gene_id(raw: Any) -> str | None:
    """Normalize an Ensembl gene ID (uppercase, drop version suffix)."""
    text = str(raw or "").strip()
    if not text:
        return None
    base = text.split(".", 1)[0].upper()
    if base.startswith(("ENST", "ENSMUST", "ENSRNOT", "ENSP", "ENSMUSP", "ENSRNOP")):
        return None
    if base.startswith(("ENSG", "ENSMUSG", "ENSRNOG")):
        return base
    return None


def normalize_uniprot_accession(raw: Any) -> str | None:
    """Normalize a UniProt accession to uppercase."""
    text = str(raw or "").strip()
    if not text:
        return None
    return text.upper()


def resolve_taxon_id(rec: EvidenceRecord) -> int | None:
    """Resolve taxon ID from record metadata, then structured value fields."""
    if rec.taxon_id is not None:
        try:
            return int(rec.taxon_id)
        except (TypeError, ValueError):
            pass
    value = rec.value if isinstance(rec.value, dict) else {}
    for key in ("tax_id", "taxon_id", "organism_id"):
        raw = value.get(key)
        if raw is None:
            continue
        try:
            return int(raw)
        except (TypeError, ValueError):
            continue
    return None


def canonical_primary_identifier(rec: EvidenceRecord) -> str | None:
    """Return the canonical primary ID used for identity-record deduplication."""
    value = rec.value if isinstance(rec.value, dict) else {}
    fact = rec.fact_type or ""
    if fact == "entrez_gene_id":
        return normalize_entrez_id(value.get("entrez_gene_id"))
    if fact == "ensembl_gene_id":
        return normalize_ensembl_gene_id(
            value.get("ensembl_gene_id") or value.get("id")
        )
    if fact == "uniprot_accession":
        return normalize_uniprot_accession(
            value.get("uniprot_accession") or value.get("primaryAccession")
        )
    if fact == "ensembl_xref":
        accession = normalize_uniprot_accession(value.get("uniprot_accession")) or ""
        gene = normalize_ensembl_gene_id(value.get("ensembl_gene_id")) or ""
        if not accession and not gene:
            return None
        return f"{accession}:{gene}"
    if fact == "canonical_transcript":
        text = str(value.get("canonical_transcript") or "").strip()
        if not text:
            return None
        # Transcript IDs keep version when present; compare case-insensitively.
        return text.upper()
    if fact == "genomic_location":
        chrom = value.get("seq_region_name")
        start = value.get("start")
        end = value.get("end")
        if chrom is None or start is None or end is None:
            return None
        assembly = (
            value.get("assembly_name")
            or value.get("assembly")
            or value.get("genome")
            or "unknown"
        )
        return f"{assembly}:{chrom}:{start}-{end}"
    return None


def identity_dedup_key(rec: EvidenceRecord) -> tuple[Any, ...] | None:
    """Build the dedup key, or ``None`` when the record is not identity-dedupable."""
    if not is_species_identity_record(rec):
        return None
    primary = canonical_primary_identifier(rec)
    if not primary:
        return None
    return (
        rec.dossier_run_id,
        rec.source_name,
        rec.fact_type,
        resolve_taxon_id(rec),
        primary,
    )


def _score_identity_record(rec: EvidenceRecord, *, query_symbol: str) -> tuple:
    """Higher score wins. Prefers corrected query-symbol records with provenance."""
    value = rec.value if isinstance(rec.value, dict) else {}
    query = (query_symbol or "").strip()
    gene_matches = 1 if (rec.gene_symbol or "").strip() == query else 0
    has_official = 1 if str(rec.official_symbol or "").strip() else 0
    has_query_field = 1 if str(value.get("query_gene_symbol") or "").strip() else 0
    has_species_field = 1 if str(value.get("species_gene_symbol") or "").strip() else 0
    has_raw = 1 if rec.raw_artifact_id else 0
    has_api = 1 if rec.api_run_id else 0
    created = rec.created_at.timestamp() if rec.created_at is not None else 0.0
    return (
        gene_matches,
        has_official,
        has_query_field,
        has_species_field,
        has_raw,
        has_api,
        created,
        rec.id or "",
    )


def assert_dossier_gene_matches(run_gene_symbol: str | None, query_symbol: str) -> None:
    """Raise ``ValueError`` when ``--gene`` does not match the dossier run gene."""
    run_gene = (run_gene_symbol or "").strip()
    query = (query_symbol or "").strip()
    if run_gene.casefold() != query.casefold():
        raise ValueError(
            f"Gene mismatch: --gene={query!r} but dossier run has "
            f"gene_symbol={run_gene!r}."
        )


def select_preferred_identity_record(
    records: Sequence[EvidenceRecord],
    *,
    query_symbol: str,
) -> EvidenceRecord:
    """Return the preferred record from a duplicate identity group."""
    if not records:
        raise ValueError("cannot select preferred identity record from empty group")
    return max(records, key=lambda r: _score_identity_record(r, query_symbol=query_symbol))


def dedupe_species_identity_records(
    records: Iterable[EvidenceRecord],
    *,
    query_symbol: str,
) -> IdentityDedupResult:
    """Partition identity records into retained vs removed duplicates.

    Non-identity records are always retained. Within each
    ``(run, source, fact_type, taxon, primary_id)`` group, only the preferred
    corrected record is retained.
    """
    identity: list[EvidenceRecord] = []
    passthrough: list[EvidenceRecord] = []
    for rec in records:
        if is_species_identity_record(rec) and identity_dedup_key(rec) is not None:
            identity.append(rec)
        else:
            passthrough.append(rec)

    groups: dict[tuple[Any, ...], list[EvidenceRecord]] = {}
    for rec in identity:
        key = identity_dedup_key(rec)
        assert key is not None
        groups.setdefault(key, []).append(rec)

    retained: list[EvidenceRecord] = list(passthrough)
    removed: list[EvidenceRecord] = []
    for group in groups.values():
        preferred = select_preferred_identity_record(group, query_symbol=query_symbol)
        retained.append(preferred)
        for rec in group:
            if rec.id != preferred.id:
                removed.append(rec)

    return IdentityDedupResult(
        retained=tuple(retained),
        removed=tuple(removed),
    )


__all__ = [
    "IDENTITY_SOURCE_NAMES",
    "IdentityDedupResult",
    "is_species_identity_record",
    "normalize_entrez_id",
    "normalize_ensembl_gene_id",
    "normalize_uniprot_accession",
    "resolve_taxon_id",
    "canonical_primary_identifier",
    "identity_dedup_key",
    "assert_dossier_gene_matches",
    "select_preferred_identity_record",
    "dedupe_species_identity_records",
]
