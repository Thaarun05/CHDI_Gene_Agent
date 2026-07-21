"""Normalize expression ToolResults into EvidenceRecords.

Consumes successful client payloads from GTEx, Harmonizome, Allen Brain Atlas,
and BrainRNASeq. Does **not** call the network.

Rules:
- Emit tissue / cell-type / eQTL facts only from payload fields
- Do not invent expression values or regional Allen summaries
- Cap GTEx eQTL EvidenceRecords (full count preserved on a summary record)
- Route Harmonizome TF datasets to Transcription factors; other associations
  to Tissue and cell expression
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

SECTION_EXPRESSION = "Tissue and cell expression"
SECTION_EQTL = "eQTLs"
SECTION_TF = "Transcription factors"

# Harmonizome TF / regulator datasets (must match tools/harmonizome.py).
TF_DATASET_NAMES = {
    "ENCODE Transcription Factor Binding Site Profiles",
    "ENCODE Transcription Factor Targets",
    "ChEA Transcription Factor Binding Site Profiles",
    "ChEA Transcription Factor Targets",
    "JASPAR Predicted Transcription Factor Targets",
    "MotifMap Predicted Transcription Factor Targets",
}

# Unbounded GTEx eQTL pages can be huge; keep dossier records bounded.
DEFAULT_MAX_EQTL_RECORDS = 100


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _data_list(payload: Any) -> list[Any]:
    """Return GTEx-style ``data`` list from a payload dict or bare list."""
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        data = payload.get("data")
        if isinstance(data, list):
            return data
    return []


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    source_type: SourceType,
    assertion_type: AssertionType,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    evidence_grade: EvidenceGrade,
    section: str,
    organism: str | None = None,
    species: str | None = None,
    taxon_id: int | None = None,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one EvidenceRecord with a deterministic source_id."""
    source_id = make_source_id(source_name, gene_symbol, assertion_type, key)
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=section,
        subsection=subsection,
        source_name=source_name,
        source_type=source_type,
        assertion_type=assertion_type,
        fact_type=fact_type,
        organism=organism,
        species=species,
        taxon_id=taxon_id,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def normalize_gtex(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_eqtl_records: int = DEFAULT_MAX_EQTL_RECORDS,
) -> list[EvidenceRecord]:
    """Normalize GTEx ``fetch_expression_and_eqtl`` (or compatible) payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    gencode_id = data.get("gencode_id")
    dataset_id = data.get("dataset_id")
    records: list[EvidenceRecord] = []
    common = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "source_name": "GTEx",
        "source_type": SourceType.expression_database,
        "organism": "Homo sapiens",
        "taxon_id": 9606,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }

    if data.get("median_expression_success", True) and data.get("median_expression") is not None:
        for row in _data_list(data.get("median_expression")):
            if not isinstance(row, dict):
                continue
            tissue = row.get("tissueSiteDetailId") or row.get("tissueSiteDetail")
            median = row.get("median")
            if tissue is None or median is None:
                continue
            unit = row.get("unit")
            records.append(
                _record(
                    **common,
                    assertion_type=AssertionType.expression,
                    fact_type="tissue_median_expression",
                    key=f"{gencode_id or gene_symbol}-{tissue}",
                    value={
                        "gencode_id": gencode_id,
                        "dataset_id": dataset_id,
                        "tissue_site_detail_id": tissue,
                        "median": median,
                        "unit": unit,
                    },
                    display_text=(
                        f"{gene_symbol} GTEx median expression in {tissue} is "
                        f"{median}"
                        + (f" {unit}" if unit else "")
                        + "."
                    ),
                    evidence_grade=EvidenceGrade.B,
                    section=SECTION_EXPRESSION,
                    subsection="GTEx median expression",
                    confidence_notes="GTEx is human-only.",
                )
            )

    if data.get("eqtl_success", True) and data.get("eqtl") is not None:
        eqtl_rows = [
            r for r in _data_list(data.get("eqtl")) if isinstance(r, dict)
        ]
        total = len(eqtl_rows)
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.expression,
                fact_type="eqtl_summary",
                key=f"{gencode_id or gene_symbol}-eqtl-summary",
                value={
                    "gencode_id": gencode_id,
                    "dataset_id": dataset_id,
                    "eqtl_count": total,
                    "records_emitted_cap": max_eqtl_records,
                },
                display_text=(
                    f"{gene_symbol} has {total} GTEx single-tissue eQTL row(s)"
                    + (f" for {gencode_id}." if gencode_id else ".")
                ),
                evidence_grade=EvidenceGrade.B,
                section=SECTION_EQTL,
                subsection="GTEx eQTL summary",
                confidence_notes="GTEx is human-only.",
            )
        )

        # Prefer stronger associations when capping.
        def _pval(row: dict[str, Any]) -> float:
            try:
                return float(row.get("pValue"))
            except (TypeError, ValueError):
                return float("inf")

        ordered = sorted(eqtl_rows, key=_pval)
        for row in ordered[: max(0, max_eqtl_records)]:
            snp = row.get("snpId") or row.get("variantId")
            tissue = row.get("tissueSiteDetailId")
            if not snp:
                continue
            key = f"{snp}-{tissue or 'na'}"
            records.append(
                _record(
                    **common,
                    assertion_type=AssertionType.expression,
                    fact_type="eqtl",
                    key=key,
                    value={
                        "gencode_id": row.get("gencodeId") or gencode_id,
                        "dataset_id": row.get("datasetId") or dataset_id,
                        "snp_id": row.get("snpId"),
                        "variant_id": row.get("variantId"),
                        "chromosome": row.get("chromosome"),
                        "pos": row.get("pos"),
                        "tissue_site_detail_id": tissue,
                        "nes": row.get("nes"),
                        "p_value": row.get("pValue"),
                    },
                    display_text=(
                        f"{gene_symbol} GTEx eQTL {snp}"
                        + (f" in {tissue}" if tissue else "")
                        + (
                            f" (NES={row.get('nes')}, p={row.get('pValue')})."
                            if row.get("nes") is not None or row.get("pValue") is not None
                            else "."
                        )
                    ),
                    evidence_grade=EvidenceGrade.B,
                    section=SECTION_EQTL,
                    subsection="GTEx single-tissue eQTL",
                    confidence_notes="GTEx is human-only.",
                )
            )

    return records


def normalize_harmonizome(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize Harmonizome association summaries (expression + TF routing)."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    summaries = data.get("association_summaries") or []
    if not isinstance(summaries, list):
        summaries = []

    records: list[EvidenceRecord] = []
    for idx, row in enumerate(summaries, start=1):
        if not isinstance(row, dict):
            continue
        dataset = row.get("dataset_name")
        attribute = row.get("attribute_name")
        if not dataset and not attribute:
            continue
        is_tf = dataset in TF_DATASET_NAMES
        key = "-".join(
            str(x)
            for x in (dataset or "dataset", attribute or "attribute", idx)
            if x is not None and str(x) != ""
        )
        if is_tf:
            records.append(
                _record(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    source_name="Harmonizome",
                    source_type=SourceType.expression_database,
                    assertion_type=AssertionType.transcription_factor_association,
                    fact_type="tf_association",
                    key=key,
                    value={
                        "dataset_name": dataset,
                        "attribute_name": attribute,
                        "attribute_href": row.get("attribute_href"),
                        "threshold_value": row.get("threshold_value"),
                        "standardized_value": row.get("standardized_value"),
                        "associated_gene_symbol": row.get("associated_gene_symbol"),
                    },
                    display_text=(
                        f"{gene_symbol} Harmonizome TF association: "
                        f"{attribute or 'attribute'} ({dataset})."
                    ),
                    evidence_grade=EvidenceGrade.E,
                    section=SECTION_TF,
                    subsection="Harmonizome TF datasets",
                    confidence_notes=(
                        "Harmonizome association from a TF/regulator dataset; "
                        "computational / curated-resource hit, not primary evidence."
                    ),
                    manual_review_required=True,
                    api_run_id=api_run_id,
                    raw_artifact_id=raw_artifact_id,
                )
            )
        else:
            records.append(
                _record(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    source_name="Harmonizome",
                    source_type=SourceType.expression_database,
                    assertion_type=AssertionType.expression,
                    fact_type="gene_association",
                    key=key,
                    value={
                        "dataset_name": dataset,
                        "attribute_name": attribute,
                        "attribute_href": row.get("attribute_href"),
                        "threshold_value": row.get("threshold_value"),
                        "standardized_value": row.get("standardized_value"),
                        "associated_gene_symbol": row.get("associated_gene_symbol"),
                    },
                    display_text=(
                        f"{gene_symbol} Harmonizome association: "
                        f"{attribute or 'attribute'} ({dataset})."
                    ),
                    evidence_grade=EvidenceGrade.E,
                    section=SECTION_EXPRESSION,
                    subsection="Harmonizome associations",
                    confidence_notes=(
                        "Harmonizome gene–attribute association; verify before "
                        "treating as primary expression evidence."
                    ),
                    manual_review_required=True,
                    api_run_id=api_run_id,
                    raw_artifact_id=raw_artifact_id,
                )
            )
    return records


def normalize_allen_brain(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize Allen HBA probe + expression fetch payloads.

    Does not invent per-structure expression values from the raw microarray
    service payload; records probe identity and that expression was retrieved.
    """
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    records: list[EvidenceRecord] = []
    common = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "source_name": "Allen Brain Atlas",
        "source_type": SourceType.expression_database,
        "organism": "Homo sapiens",
        "taxon_id": 9606,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }

    summaries = data.get("probe_summaries") or []
    if not isinstance(summaries, list):
        summaries = []
    for probe in summaries:
        if not isinstance(probe, dict):
            continue
        probe_id = probe.get("id")
        if probe_id is None:
            continue
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.expression,
                fact_type="hba_probe",
                key=str(probe_id),
                value={
                    "probe_id": probe_id,
                    "probe_name": probe.get("name"),
                    "gene_acronym": probe.get("gene_acronym"),
                    "gene_name": probe.get("gene_name"),
                    "entrez_id": probe.get("entrez_id"),
                },
                display_text=(
                    f"{gene_symbol} Allen HBA microarray probe {probe_id}"
                    + (f" ({probe.get('name')})." if probe.get("name") else ".")
                ),
                evidence_grade=EvidenceGrade.B,
                section=SECTION_EXPRESSION,
                subsection="Allen HBA probes",
            )
        )

    expressions = data.get("expressions") or {}
    if not isinstance(expressions, dict):
        expressions = {}
    for probe_id, payload in expressions.items():
        msg = _as_dict(payload).get("msg")
        msg_dict = msg if isinstance(msg, dict) else {}
        probes = msg_dict.get("probes")
        samples = msg_dict.get("samples")
        expression = msg_dict.get("expression")
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.expression,
                fact_type="hba_expression_payload",
                key=f"expr-{probe_id}",
                value={
                    "probe_id": probe_id,
                    "probe_count": len(probes) if isinstance(probes, list) else None,
                    "sample_count": len(samples) if isinstance(samples, list) else None,
                    "expression_matrix_present": expression is not None,
                    "caveat": (
                        "Raw Allen microarray expression payload preserved; "
                        "regional values are not invented here."
                    ),
                },
                display_text=(
                    f"{gene_symbol} Allen HBA expression retrieved for probe "
                    f"{probe_id}."
                ),
                evidence_grade=EvidenceGrade.B,
                section=SECTION_EXPRESSION,
                subsection="Allen HBA expression",
                confidence_notes=(
                    "Expression payload retrieved; do not invent structure-level "
                    "values without explicit parsing."
                ),
            )
        )
    return records


def normalize_brainrnaseq(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize BrainRNASeq matched gene expression summaries."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    species = str(data.get("species") or "human")
    summaries = data.get("expression_summaries") or []
    if not isinstance(summaries, list):
        summaries = []

    organism = "Mus musculus" if species == "mouse" else "Homo sapiens"
    taxon_id = 10090 if species == "mouse" else 9606
    grade = EvidenceGrade.D if species == "mouse" else EvidenceGrade.B

    records: list[EvidenceRecord] = []
    for idx, row in enumerate(summaries, start=1):
        if not isinstance(row, dict):
            continue
        cell_types = row.get("cell_type_values") or {}
        if not isinstance(cell_types, dict) or not cell_types:
            continue
        gene_id = row.get("gene_id") or row.get("id") or gene_symbol
        key = f"{species}-{gene_id}-{idx}"
        # Compact display: list a few cell-type keys, full map in value.
        preview_keys = list(cell_types.keys())[:5]
        preview = ", ".join(preview_keys)
        if len(cell_types) > 5:
            preview += f", … ({len(cell_types)} columns)"
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="BrainRNASeq",
                source_type=SourceType.expression_database,
                assertion_type=AssertionType.cell_type_expression,
                fact_type="cell_type_expression",
                key=key,
                value={
                    "species": species,
                    "gene_id": row.get("gene_id"),
                    "id": row.get("id"),
                    "cell_type_values": cell_types,
                    "cell_type_count": row.get("cell_type_count") or len(cell_types),
                },
                display_text=(
                    f"{gene_symbol} BrainRNASeq ({species}) cell-type expression "
                    f"columns include {preview}."
                ),
                evidence_grade=grade,
                section=SECTION_EXPRESSION,
                subsection=f"BrainRNASeq {species}",
                organism=organism,
                species=species,
                taxon_id=taxon_id,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_expression(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch expression normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    if source == "GTEx":
        return normalize_gtex(tool_result, **kwargs)
    if source == "Harmonizome":
        return normalize_harmonizome(tool_result, **kwargs)
    if source == "Allen Brain Atlas":
        return normalize_allen_brain(tool_result, **kwargs)
    if source == "BrainRNASeq":
        return normalize_brainrnaseq(tool_result, **kwargs)
    return []


__all__ = [
    "DEFAULT_MAX_EQTL_RECORDS",
    "normalize_gtex",
    "normalize_harmonizome",
    "normalize_allen_brain",
    "normalize_brainrnaseq",
    "normalize_expression",
]
