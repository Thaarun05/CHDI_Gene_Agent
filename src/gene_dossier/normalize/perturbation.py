"""Normalize perturbation ToolResults into EvidenceRecords.

Consumes successful GEO ``fetch_perturbations`` payloads. Does **not** call
the network.

Rules:
- One EvidenceRecord per GDS summary
- GEO search hits are discovery/mention evidence, not proof the dataset
  causally perturbs the gene
- Do not invent expression values, sample phenotypes, or treatment effects
- Preserve search organism/context/term metadata from the client payload
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

SECTION_GEO = "GEO perturbations"

DEFAULT_MAX_GDS = 100


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    subsection: str | None = None,
    organism: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = True,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one GEO perturbation EvidenceRecord."""
    source_id = make_source_id(
        "GEO", gene_symbol, AssertionType.perturbation, key
    )
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_GEO,
        subsection=subsection,
        source_name="GEO",
        source_type=SourceType.expression_database,
        assertion_type=AssertionType.perturbation,
        fact_type=fact_type,
        organism=organism,
        evidence_grade=EvidenceGrade.F,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def normalize_geo(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_records: int = DEFAULT_MAX_GDS,
) -> list[EvidenceRecord]:
    """Normalize GEO ``fetch_perturbations`` GDS summaries."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    summaries = data.get("gds_summaries") or []
    if not isinstance(summaries, list):
        summaries = []

    organism = data.get("organism")
    context = data.get("context")
    search_term = data.get("search_term")
    caveat = (
        "GEO Profiles→GDS search hit; dataset membership is not proof that the "
        "experiment causally perturbs this gene."
    )

    records: list[EvidenceRecord] = []
    for row in summaries[: max(0, max_records)]:
        if not isinstance(row, dict):
            continue
        gds_uid = row.get("gds_uid")
        accession = row.get("accession")
        if not gds_uid and not accession:
            continue
        title = row.get("title")
        key = str(accession or gds_uid)
        taxon = row.get("taxon") or organism
        display = f"{gene_symbol} GEO dataset {accession or gds_uid}"
        if title:
            display += f": {title}"
        display += "."

        # Keep sample list compact in value; do not invent phenotype labels.
        samples = row.get("samples") or []
        if not isinstance(samples, list):
            samples = []
        sample_accessions = [
            s.get("accession")
            for s in samples
            if isinstance(s, dict) and s.get("accession")
        ]

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                fact_type="geo_gds_dataset",
                key=key,
                value={
                    "gds_uid": gds_uid,
                    "accession": accession,
                    "title": title,
                    "summary": row.get("summary"),
                    "gpl": row.get("gpl"),
                    "gse": row.get("gse"),
                    "taxon": row.get("taxon"),
                    "gdstype": row.get("gdstype"),
                    "valtype": row.get("valtype"),
                    "ssinfo": row.get("ssinfo"),
                    "subsetinfo": row.get("subsetinfo"),
                    "n_samples": row.get("n_samples"),
                    "sample_accessions": sample_accessions,
                    "sample_count": len(sample_accessions) or row.get("n_samples"),
                    "pubmedids": row.get("pubmedids"),
                    "ftplink": row.get("ftplink"),
                    "search_organism": organism,
                    "search_context": context,
                    "search_term": search_term,
                    "profile_id_count": len(data.get("profile_ids") or [])
                    if isinstance(data.get("profile_ids"), list)
                    else None,
                    "caveat": caveat,
                },
                display_text=display,
                subsection="GEO DataSets",
                organism=str(taxon) if taxon else None,
                confidence_notes=caveat,
                manual_review_required=True,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_perturbation(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch perturbation normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    if source == "GEO":
        return normalize_geo(
            tool_result,
            dossier_run_id=dossier_run_id,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    return []


__all__ = [
    "DEFAULT_MAX_GDS",
    "normalize_geo",
    "normalize_perturbation",
]
