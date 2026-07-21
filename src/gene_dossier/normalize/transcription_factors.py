"""Normalize transcription-factor ToolResults into EvidenceRecords.

Consumes successful Harmonizome payloads, preferring ``tf_summaries`` (or
TF-dataset rows within ``association_summaries``). Does **not** call the
network.

Rules:
- Emit only validated TF / regulator dataset associations
- Do not invent TF–target relationships beyond payload fields
- Computational / curated-resource hits are grade E and need review
- Non-TF Harmonizome associations belong in ``normalize/expression.py``
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

SECTION_TF = "Transcription factors"

# Must match tools/harmonizome.py TF_DATASET_NAMES.
TF_DATASET_NAMES = {
    "ENCODE Transcription Factor Binding Site Profiles",
    "ENCODE Transcription Factor Targets",
    "ChEA Transcription Factor Binding Site Profiles",
    "ChEA Transcription Factor Targets",
    "JASPAR Predicted Transcription Factor Targets",
    "MotifMap Predicted Transcription Factor Targets",
}

PREDICTED_TF_DATASETS = {
    "JASPAR Predicted Transcription Factor Targets",
    "MotifMap Predicted Transcription Factor Targets",
}

DEFAULT_MAX_TF = 200


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
    confidence_notes: str | None = None,
    manual_review_required: bool = True,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one transcription-factor EvidenceRecord."""
    source_id = make_source_id(
        "Harmonizome",
        gene_symbol,
        AssertionType.transcription_factor_association,
        key,
    )
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_TF,
        subsection=subsection,
        source_name="Harmonizome",
        source_type=SourceType.expression_database,
        assertion_type=AssertionType.transcription_factor_association,
        fact_type=fact_type,
        evidence_grade=EvidenceGrade.E,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def _tf_summary_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Prefer ``tf_summaries``; else TF-filtered ``association_summaries``."""
    tf_summaries = data.get("tf_summaries")
    if isinstance(tf_summaries, list) and tf_summaries:
        return [row for row in tf_summaries if isinstance(row, dict)]

    associations = data.get("association_summaries") or []
    if not isinstance(associations, list):
        return []
    out: list[dict[str, Any]] = []
    for row in associations:
        if not isinstance(row, dict):
            continue
        if row.get("dataset_name") in TF_DATASET_NAMES:
            out.append(row)
    return out


def normalize_harmonizome_tf(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_records: int = DEFAULT_MAX_TF,
) -> list[EvidenceRecord]:
    """Normalize Harmonizome TF / regulator association summaries only."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    rows = _tf_summary_rows(data)
    records: list[EvidenceRecord] = []
    for idx, row in enumerate(rows[: max(0, max_records)], start=1):
        dataset = row.get("dataset_name")
        attribute = row.get("attribute_name")
        if not dataset and not attribute:
            continue
        # Defense in depth if tf_summaries ever contains non-TF rows.
        if dataset and dataset not in TF_DATASET_NAMES:
            continue

        predicted = dataset in PREDICTED_TF_DATASETS
        fact_type = (
            "tf_predicted_association" if predicted else "tf_association"
        )
        notes = (
            "Predicted Harmonizome TF association (motif/resource prediction); "
            "not primary experimental evidence."
            if predicted
            else (
                "Harmonizome association from a TF/regulator dataset; "
                "computational / curated-resource hit, not primary evidence."
            )
        )
        key = "-".join(
            str(x)
            for x in (dataset or "dataset", attribute or "attribute", idx)
            if x is not None and str(x) != ""
        )
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                fact_type=fact_type,
                key=key,
                value={
                    "dataset_name": dataset,
                    "attribute_name": attribute,
                    "attribute_href": row.get("attribute_href"),
                    "threshold_value": row.get("threshold_value"),
                    "standardized_value": row.get("standardized_value"),
                    "associated_gene_symbol": row.get("associated_gene_symbol"),
                    "predicted_dataset": predicted,
                    "caveat": notes,
                },
                display_text=(
                    f"{gene_symbol} Harmonizome TF association: "
                    f"{attribute or 'attribute'} ({dataset})."
                ),
                subsection=(
                    "Harmonizome predicted TF datasets"
                    if predicted
                    else "Harmonizome TF datasets"
                ),
                confidence_notes=notes,
                manual_review_required=True,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_transcription_factors(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch TF normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    if source == "Harmonizome":
        return normalize_harmonizome_tf(
            tool_result,
            dossier_run_id=dossier_run_id,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    return []


__all__ = [
    "TF_DATASET_NAMES",
    "PREDICTED_TF_DATASETS",
    "DEFAULT_MAX_TF",
    "normalize_harmonizome_tf",
    "normalize_transcription_factors",
]
