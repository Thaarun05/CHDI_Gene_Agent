"""Normalize grant ToolResults into EvidenceRecords.

Consumes successful NIH RePORTER ``fetch_grants`` payloads. Does **not** call
the network.

Rules:
- One EvidenceRecord per project summary
- Exact search hits are weak mention evidence (grade F), not proof the grant
  focuses on the gene
- Broader pathway search hits are weaker still and require manual review
- Do not invent award amounts, PIs, or relevance beyond payload fields
- Prefer exact-search records when the same project appears in both result sets
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

SECTION_GRANTS = "NIH/ERC grants"

DEFAULT_MAX_EXACT = 100
DEFAULT_MAX_BROADER = 50


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _project_key(row: dict[str, Any]) -> str | None:
    """Stable project identity from project_num / core_project_num."""
    for key in ("project_num", "core_project_num"):
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = True,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one NIH RePORTER grant EvidenceRecord."""
    source_id = make_source_id(
        "NIH RePORTER", gene_symbol, AssertionType.grant_project, key
    )
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_GRANTS,
        subsection=subsection,
        source_name="NIH RePORTER",
        source_type=SourceType.grant_database,
        assertion_type=AssertionType.grant_project,
        fact_type=fact_type,
        evidence_grade=EvidenceGrade.F,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def _project_records(
    summaries: list[Any],
    *,
    dossier_run_id: str,
    gene_symbol: str,
    search_kind: str,
    search_text: str | None,
    seen_keys: set[str],
    max_records: int,
    api_run_id: str | None,
    raw_artifact_id: str | None,
) -> list[EvidenceRecord]:
    """Emit grant records for one search result set."""
    records: list[EvidenceRecord] = []
    if search_kind == "exact":
        fact_type = "nih_reporter_project"
        subsection = "NIH RePORTER exact search"
        notes = (
            "NIH RePORTER exact text search hit; project mention is not proof "
            "the grant focuses on this gene."
        )
    else:
        fact_type = "nih_reporter_broader_project"
        subsection = "NIH RePORTER broader pathway search"
        notes = (
            "NIH RePORTER broader pathway search hit; may be a false positive "
            "and requires review before citing as gene-focused funding."
        )

    for idx, row in enumerate(summaries, start=1):
        if len(records) >= max(0, max_records):
            break
        if not isinstance(row, dict):
            continue
        project_key = _project_key(row)
        dedupe = project_key or f"{search_kind}-{idx}"
        if dedupe in seen_keys:
            continue
        seen_keys.add(dedupe)

        title = row.get("project_title")
        project_num = row.get("project_num") or row.get("core_project_num")
        org = row.get("organization_name")
        display = f"{gene_symbol} NIH project"
        if project_num:
            display += f" {project_num}"
        if title:
            display += f": {title}"
        display += "."

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                fact_type=fact_type,
                key=dedupe,
                value={
                    "search_kind": search_kind,
                    "search_text": search_text,
                    "project_title": title,
                    "project_num": row.get("project_num"),
                    "core_project_num": row.get("core_project_num"),
                    "organization_name": org,
                    "budget_start": row.get("budget_start"),
                    "budget_end": row.get("budget_end"),
                    "agency_ic_admin": row.get("agency_ic_admin"),
                    "agency_ic_fundings": row.get("agency_ic_fundings"),
                    "fiscal_year": row.get("fiscal_year"),
                    "award_amount": row.get("award_amount"),
                    "principal_investigators": row.get("principal_investigators"),
                    "project_detail_url": row.get("project_detail_url"),
                    "abstract_text": row.get("abstract_text"),
                    "terms": row.get("terms"),
                    "caveat": notes,
                },
                display_text=display,
                subsection=subsection,
                confidence_notes=notes,
                manual_review_required=True,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_nih_reporter(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_exact: int = DEFAULT_MAX_EXACT,
    max_broader: int = DEFAULT_MAX_BROADER,
) -> list[EvidenceRecord]:
    """Normalize NIH RePORTER ``fetch_grants`` exact/broader summaries."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    seen: set[str] = set()
    records: list[EvidenceRecord] = []

    exact = data.get("exact_summaries") or []
    if isinstance(exact, list):
        records.extend(
            _project_records(
                exact,
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                search_kind="exact",
                search_text=data.get("exact_search_text"),
                seen_keys=seen,
                max_records=max_exact,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )

    broader = data.get("broader_summaries") or []
    if isinstance(broader, list) and broader:
        records.extend(
            _project_records(
                broader,
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                search_kind="broader",
                search_text=data.get("broader_search_text"),
                seen_keys=seen,
                max_records=max_broader,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_grants(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch grant normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    if source == "NIH RePORTER":
        return normalize_nih_reporter(
            tool_result,
            dossier_run_id=dossier_run_id,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    return []


__all__ = [
    "DEFAULT_MAX_EXACT",
    "DEFAULT_MAX_BROADER",
    "normalize_nih_reporter",
    "normalize_grants",
]
