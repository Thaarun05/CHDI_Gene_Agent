"""Normalize pathway ToolResults into EvidenceRecords.

Consumes successful client payloads from Reactome and WikiPathways. Does **not**
call the network.

Rules:
- One EvidenceRecord per pathway summary (or pathway-detail payload)
- Do not invent pathway membership or biology beyond payload fields
- Reactome UniProt mappings are curated (grade C)
- WikiPathways identifier-filtered matches are curated community pathways (grade C)
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

SECTION_PATHWAYS = "Pathways"


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
    species: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one pathway EvidenceRecord."""
    source_id = make_source_id(
        source_name, gene_symbol, AssertionType.pathway_membership, key
    )
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_PATHWAYS,
        subsection=subsection,
        source_name=source_name,
        source_type=SourceType.pathway_database,
        assertion_type=AssertionType.pathway_membership,
        fact_type=fact_type,
        organism=organism,
        species=species,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def normalize_reactome(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize Reactome ``fetch_pathways`` or ``pathway_detail`` payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    accession = data.get("uniprot_accession")
    records: list[EvidenceRecord] = []

    summaries = data.get("pathway_summaries")
    if isinstance(summaries, list):
        for row in summaries:
            if not isinstance(row, dict):
                continue
            st_id = row.get("st_id")
            display_name = row.get("display_name")
            if not st_id and not display_name:
                continue
            key = str(st_id or display_name)
            species_name = row.get("species_name")
            records.append(
                _record(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    source_name="Reactome",
                    fact_type="reactome_pathway",
                    key=key,
                    value={
                        "uniprot_accession": accession,
                        "db_id": row.get("db_id"),
                        "st_id": st_id,
                        "st_id_version": row.get("st_id_version"),
                        "display_name": display_name,
                        "species_name": species_name,
                        "doi": row.get("doi"),
                        "has_diagram": row.get("has_diagram"),
                        "release_date": row.get("release_date"),
                        "last_updated_date": row.get("last_updated_date"),
                        "schema_class": row.get("schema_class"),
                        "detail_url": row.get("detail_url"),
                        "browser_url": row.get("browser_url"),
                    },
                    display_text=(
                        f"{gene_symbol} Reactome pathway"
                        + (f" {display_name}" if display_name else "")
                        + (f" ({st_id})." if st_id else ".")
                    ),
                    evidence_grade=EvidenceGrade.C,
                    subsection="Reactome pathways",
                    organism=str(species_name) if species_name else None,
                    species=str(species_name) if species_name else None,
                    api_run_id=api_run_id,
                    raw_artifact_id=raw_artifact_id,
                )
            )
        return records

    # Optional single pathway-detail payload (summarize_pathway_detail shape or raw).
    detail = data
    if isinstance(data.get("summary"), dict):
        detail = data["summary"]
    st_id = detail.get("st_id") or detail.get("stId")
    display_name = detail.get("display_name") or detail.get("displayName")
    if st_id or display_name:
        summation_texts = detail.get("summation_texts") or []
        if not isinstance(summation_texts, list):
            summation_texts = []
        pubmed_ids = detail.get("pubmed_ids") or []
        if not isinstance(pubmed_ids, list):
            pubmed_ids = []
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="Reactome",
                fact_type="reactome_pathway_detail",
                key=str(st_id or display_name),
                value={
                    "st_id": st_id,
                    "display_name": display_name,
                    "summation_texts": summation_texts,
                    "pubmed_ids": pubmed_ids,
                    "detail_url": detail.get("detail_url"),
                },
                display_text=(
                    f"{gene_symbol} Reactome pathway detail"
                    + (f" {display_name}" if display_name else "")
                    + (f" ({st_id})." if st_id else ".")
                ),
                evidence_grade=EvidenceGrade.C,
                subsection="Reactome pathway detail",
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def _wikipathways_records_from_summaries(
    summaries: list[Any],
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_bulk: str,
    filter_meta: dict[str, Any],
    seen_ids: set[str],
    api_run_id: str | None,
    raw_artifact_id: str | None,
) -> list[EvidenceRecord]:
    """Build WikiPathways EvidenceRecords from pathway summary dicts."""
    records: list[EvidenceRecord] = []
    for idx, row in enumerate(summaries, start=1):
        if not isinstance(row, dict):
            continue
        pathway_id = row.get("id")
        name = row.get("name")
        if not pathway_id and not name:
            continue
        id_key = str(pathway_id) if pathway_id is not None else f"name-{name}-{idx}"
        if pathway_id is not None and str(pathway_id) in seen_ids:
            continue
        if pathway_id is not None:
            seen_ids.add(str(pathway_id))

        species = row.get("species")
        if source_bulk == "text":
            evidence_grade = EvidenceGrade.E
            manual_review_required = True
            fact_type = "wikipathways_text_match"
            confidence_notes = (
                "Matched by WikiPathways text-symbol search; may be a false "
                "positive and needs review."
            )
        else:
            evidence_grade = EvidenceGrade.C
            manual_review_required = False
            fact_type = "wikipathways_pathway"
            confidence_notes = (
                "Matched via WikiPathways bulk JSON structured identifier filter "
                f"({source_bulk})."
            )

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="WikiPathways",
                fact_type=fact_type,
                key=id_key,
                value={
                    "pathway_id": pathway_id,
                    "name": name,
                    "url": row.get("url"),
                    "species": species,
                    "revision": row.get("revision"),
                    "description": row.get("description"),
                    "ncbigene": row.get("ncbigene"),
                    "uniprot": row.get("uniprot"),
                    "ensembl": row.get("ensembl"),
                    "source_bulk": source_bulk,
                    **filter_meta,
                },
                display_text=(
                    f"{gene_symbol} WikiPathways"
                    + (f" {name}" if name else "")
                    + (f" ({pathway_id})." if pathway_id else ".")
                ),
                evidence_grade=evidence_grade,
                subsection="WikiPathways",
                organism=str(species) if species else None,
                species=str(species) if species else None,
                confidence_notes=confidence_notes,
                manual_review_required=manual_review_required,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_wikipathways(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize WikiPathways ``fetch_pathways`` matched summaries."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    filter_meta = {
        "entrez_ids": data.get("entrez_ids"),
        "uniprot_ids": data.get("uniprot_ids"),
        "ensembl_ids": data.get("ensembl_ids"),
        "allow_text_symbol_match": data.get("allow_text_symbol_match"),
    }
    seen: set[str] = set()
    records: list[EvidenceRecord] = []

    xref_summaries = data.get("pathway_summaries") or []
    if isinstance(xref_summaries, list):
        records.extend(
            _wikipathways_records_from_summaries(
                xref_summaries,
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_bulk="xref",
                filter_meta=filter_meta,
                seen_ids=seen,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )

    text_summaries = data.get("text_pathway_summaries") or []
    if isinstance(text_summaries, list) and text_summaries:
        records.extend(
            _wikipathways_records_from_summaries(
                text_summaries,
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_bulk="text",
                filter_meta=filter_meta,
                seen_ids=seen,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_pathways(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch pathway normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    if source == "Reactome":
        return normalize_reactome(tool_result, **kwargs)
    if source == "WikiPathways":
        return normalize_wikipathways(tool_result, **kwargs)
    return []


__all__ = [
    "normalize_reactome",
    "normalize_wikipathways",
    "normalize_pathways",
]
