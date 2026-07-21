"""Normalize PPI ToolResults into EvidenceRecords.

Consumes successful client payloads from STRING and BioGRID. Does **not** call
the network.

Rules:
- One EvidenceRecord per interaction / partner row present in the payload
- Do not invent partners, scores, or experimental systems
- STRING combined scores are computational (grade E)
- BioGRID experimental interactions are curated PPI (grade C)
"""

from __future__ import annotations

from typing import Any

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
    ToolResult,
)
from gene_dossier.source_ids import make_source_id

SECTION_PPI = "Protein-protein interactions"


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    evidence_grade: EvidenceGrade,
    subsection: str | None = None,
    organism: str | None = None,
    taxon_id: int | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one PPI EvidenceRecord."""
    source_id = make_source_id(source_name, gene_symbol, AssertionType.ppi, key)
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_PPI,
        subsection=subsection,
        source_name=source_name,
        source_type=SourceType.interaction_database,
        assertion_type=AssertionType.ppi,
        fact_type=fact_type,
        organism=organism,
        taxon_id=taxon_id,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def _partner_from_string_row(
    row: dict[str, Any],
    *,
    gene_symbol: str,
    string_id: str | None,
) -> tuple[str | None, str | None, str | None, str | None]:
    """Return (partner_name, partner_string_id, name_a, name_b) for a STRING row."""
    name_a = row.get("preferredName_A")
    name_b = row.get("preferredName_B")
    id_a = row.get("stringId_A")
    id_b = row.get("stringId_B")
    target = gene_symbol.strip().upper()

    a_match = str(name_a or "").strip().upper() == target
    b_match = str(name_b or "").strip().upper() == target
    if string_id:
        if str(id_a) == str(string_id):
            a_match = True
        if str(id_b) == str(string_id):
            b_match = True

    if a_match and not b_match:
        return (
            str(name_b) if name_b else None,
            str(id_b) if id_b else None,
            str(name_a) if name_a else None,
            str(name_b) if name_b else None,
        )
    if b_match and not a_match:
        return (
            str(name_a) if name_a else None,
            str(id_a) if id_a else None,
            str(name_a) if name_a else None,
            str(name_b) if name_b else None,
        )
    # Ambiguous / neither side matches: keep both names, no guessed partner.
    return (
        None,
        None,
        str(name_a) if name_a else None,
        str(name_b) if name_b else None,
    )


def _partner_from_biogrid_row(
    row: dict[str, Any],
    *,
    gene_symbol: str,
) -> tuple[str | None, str | None, Any, Any]:
    """Return (partner_symbol, query_side, symbol_a, symbol_b)."""
    symbol_a = row.get("official_symbol_a")
    symbol_b = row.get("official_symbol_b")
    target = gene_symbol.strip().upper()
    a_match = str(symbol_a or "").strip().upper() == target
    b_match = str(symbol_b or "").strip().upper() == target
    if a_match and not b_match:
        return (
            str(symbol_b) if symbol_b else None,
            "A",
            symbol_a,
            symbol_b,
        )
    if b_match and not a_match:
        return (
            str(symbol_a) if symbol_a else None,
            "B",
            symbol_a,
            symbol_b,
        )
    return (None, None, symbol_a, symbol_b)


def normalize_string(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize STRING ``fetch_interaction_partners`` payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    string_id = data.get("string_id")
    species = data.get("species")
    try:
        taxon_id = int(species) if species is not None else None
    except (TypeError, ValueError):
        taxon_id = None
    organism = "Homo sapiens" if taxon_id == 9606 else None

    partners = data.get("partners") or []
    if not isinstance(partners, list):
        partners = []

    records: list[EvidenceRecord] = []
    for idx, row in enumerate(partners, start=1):
        if not isinstance(row, dict):
            continue
        partner_name, partner_sid, name_a, name_b = _partner_from_string_row(
            row, gene_symbol=gene_symbol, string_id=str(string_id) if string_id else None
        )
        score = row.get("score")
        key_partner = partner_name or partner_sid or f"row-{idx}"
        key = f"{string_id or gene_symbol}-{key_partner}-{idx}"
        if partner_name:
            display = (
                f"{gene_symbol} STRING partner {partner_name}"
                + (f" (score={score})." if score is not None else ".")
            )
        else:
            display = (
                f"{gene_symbol} STRING interaction"
                + (f" {name_a}–{name_b}" if name_a or name_b else "")
                + (f" (score={score})." if score is not None else ".")
            )
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="STRING",
                fact_type="string_partner",
                key=key,
                value={
                    "string_id": string_id,
                    "partner_preferred_name": partner_name,
                    "partner_string_id": partner_sid,
                    "preferred_name_a": name_a,
                    "preferred_name_b": name_b,
                    "string_id_a": row.get("stringId_A"),
                    "string_id_b": row.get("stringId_B"),
                    "score": score,
                    "nscore": row.get("nscore"),
                    "fscore": row.get("fscore"),
                    "pscore": row.get("pscore"),
                    "ascore": row.get("ascore"),
                    "escore": row.get("escore"),
                    "dscore": row.get("dscore"),
                    "tscore": row.get("tscore"),
                    "species": species,
                },
                display_text=display,
                evidence_grade=EvidenceGrade.E,
                subsection="STRING interaction partners",
                organism=organism,
                taxon_id=taxon_id,
                confidence_notes=(
                    "STRING combined score is computational/integrated evidence, "
                    "not a single experimental assay."
                ),
                manual_review_required=partner_name is None,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_biogrid(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize BioGRID ``fetch_interactions`` summary payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    tax_id = data.get("tax_id")
    try:
        taxon_id = int(tax_id) if tax_id is not None else None
    except (TypeError, ValueError):
        taxon_id = None
    organism = "Homo sapiens" if taxon_id == 9606 else None

    summaries = data.get("interaction_summaries") or []
    if not isinstance(summaries, list):
        summaries = []

    records: list[EvidenceRecord] = []
    for idx, row in enumerate(summaries, start=1):
        if not isinstance(row, dict):
            continue
        interaction_id = row.get("biogrid_interaction_id")
        partner, query_side, symbol_a, symbol_b = _partner_from_biogrid_row(
            row, gene_symbol=gene_symbol
        )
        system = row.get("experimental_system")
        pubmed_id = row.get("pubmed_id")
        key = str(interaction_id) if interaction_id is not None else f"row-{idx}"
        if partner:
            display = (
                f"{gene_symbol} BioGRID interactor {partner}"
                + (f" via {system}" if system else "")
                + (f" (PMID {pubmed_id})." if pubmed_id else ".")
            )
        else:
            display = (
                f"{gene_symbol} BioGRID interaction"
                + (f" {symbol_a}–{symbol_b}" if symbol_a or symbol_b else "")
                + (f" via {system}" if system else "")
                + (f" (PMID {pubmed_id})." if pubmed_id else ".")
            )
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="BioGRID",
                fact_type="biogrid_interaction",
                key=key,
                value={
                    "biogrid_interaction_id": interaction_id,
                    "partner_symbol": partner,
                    "query_side": query_side,
                    "official_symbol_a": symbol_a,
                    "official_symbol_b": symbol_b,
                    "entrez_gene_a": row.get("entrez_gene_a"),
                    "entrez_gene_b": row.get("entrez_gene_b"),
                    "experimental_system": system,
                    "experimental_system_type": row.get("experimental_system_type"),
                    "author": row.get("author"),
                    "pubmed_id": pubmed_id,
                    "organism_a": row.get("organism_a"),
                    "organism_b": row.get("organism_b"),
                    "throughput": row.get("throughput"),
                    "source_database": row.get("source_database"),
                    "tax_id": tax_id,
                },
                display_text=display,
                evidence_grade=EvidenceGrade.C,
                subsection="BioGRID interactions",
                organism=organism,
                taxon_id=taxon_id,
                manual_review_required=partner is None,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_ppi(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch PPI normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    if source == "STRING":
        return normalize_string(tool_result, **kwargs)
    if source == "BioGRID":
        return normalize_biogrid(tool_result, **kwargs)
    return []


__all__ = [
    "normalize_string",
    "normalize_biogrid",
    "normalize_ppi",
]
