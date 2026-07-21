"""Normalize variant / disease-association ToolResults into EvidenceRecords.

Consumes successful client payloads from ClinVar, Open Targets, and OMIM.
Does **not** call the network.

Rules:
- Do not invent pathogenicity, disease causality, or OMIM phenotype links
- Skip OMIM when ``selection_method`` is ``ambiguous`` / no selected MIM
- Open Targets association scores are computational (grade E)
- ClinVar classifications are curated reporting (grade C), not automatic causal proof
- Preserve client caveats on OMIM phenotype absence
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

SECTION_VARIANTS = "ClinVar / OMIM / Open Targets / SNPs"
SECTION_CHEMICAL = "Chemical tools"

DEFAULT_MAX_CLINVAR = 200
DEFAULT_MAX_OT_ASSOCIATIONS = 100


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


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
    section: str = SECTION_VARIANTS,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one variant/disease EvidenceRecord."""
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
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def normalize_clinvar(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_records: int = DEFAULT_MAX_CLINVAR,
) -> list[EvidenceRecord]:
    """Normalize ClinVar ``fetch_clinvar_variants`` summaries."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    summaries = data.get("variant_summaries") or []
    if not isinstance(summaries, list):
        summaries = []

    records: list[EvidenceRecord] = []
    for row in summaries[: max(0, max_records)]:
        if not isinstance(row, dict):
            continue
        uid = row.get("uid")
        accession = row.get("accession")
        if not uid and not accession:
            continue
        classification = row.get("germline_classification")
        title = row.get("title") or row.get("variation_name")
        traits = row.get("trait_names") or []
        if not isinstance(traits, list):
            traits = []
        key = str(accession or uid)
        review = False
        notes = (
            "ClinVar classification is curated reporting; not automatic proof of "
            "disease causality for this gene."
        )
        if classification:
            clf = str(classification).lower()
            if "uncertain" in clf or "conflicting" in clf or "not provided" in clf:
                review = True

        display = f"{gene_symbol} ClinVar {accession or uid}"
        if title:
            display += f": {title}"
        if classification:
            display += f" [{classification}]"
        display += "."

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="ClinVar",
                source_type=SourceType.genetic_database,
                assertion_type=AssertionType.variant_association,
                fact_type="clinvar_variant",
                key=key,
                value={
                    "uid": uid,
                    "accession": accession,
                    "title": row.get("title"),
                    "obj_type": row.get("obj_type"),
                    "variation_name": row.get("variation_name"),
                    "cdna_change": row.get("cdna_change"),
                    "protein_change": row.get("protein_change"),
                    "canonical_spdi": row.get("canonical_spdi"),
                    "variation_locs": row.get("variation_locs"),
                    "gene_symbol_reported": row.get("gene_symbol"),
                    "gene_id": row.get("gene_id"),
                    "germline_classification": classification,
                    "review_status": row.get("review_status"),
                    "last_evaluated": row.get("last_evaluated"),
                    "trait_names": traits,
                    "molecular_consequence_list": row.get(
                        "molecular_consequence_list"
                    ),
                    "search_term": data.get("search_term"),
                    "caveat": notes,
                },
                display_text=display,
                evidence_grade=EvidenceGrade.C,
                subsection="ClinVar variants",
                confidence_notes=notes,
                manual_review_required=review,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_opentargets(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_associations: int = DEFAULT_MAX_OT_ASSOCIATIONS,
) -> list[EvidenceRecord]:
    """Normalize Open Targets disease-association or tractability payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = (
        tool_result.gene_symbol
        or str(data.get("gene_symbol") or data.get("approved_symbol") or "")
    )
    if not gene_symbol:
        return []

    ensembl_id = data.get("ensembl_id") or data.get("target_id")
    records: list[EvidenceRecord] = []
    common = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "source_name": "Open Targets",
        "source_type": SourceType.genetic_database,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }

    associations = data.get("association_summaries") or []
    if isinstance(associations, list) and associations:
        # Prefer higher scores first when capping.
        def _score(row: dict[str, Any]) -> float:
            try:
                return float(row.get("score") or 0.0)
            except (TypeError, ValueError):
                return 0.0

        ordered = sorted(
            [r for r in associations if isinstance(r, dict)],
            key=_score,
            reverse=True,
        )
        for row in ordered[: max(0, max_associations)]:
            disease_id = row.get("disease_id")
            disease_name = row.get("disease_name")
            if not disease_id and not disease_name:
                continue
            score = row.get("score")
            key = str(disease_id or disease_name)
            records.append(
                _record(
                    **common,
                    assertion_type=AssertionType.disease_association,
                    fact_type="opentargets_disease_association",
                    key=key,
                    value={
                        "ensembl_id": ensembl_id,
                        "approved_symbol": data.get("approved_symbol"),
                        "approved_name": data.get("approved_name"),
                        "disease_id": disease_id,
                        "disease_name": disease_name,
                        "score": score,
                        "datatype_scores": row.get("datatype_scores"),
                        "datasource_scores": row.get("datasource_scores"),
                        "caveat": (
                            "Open Targets association score is computational "
                            "evidence of association, not proof of causality."
                        ),
                    },
                    display_text=(
                        f"{gene_symbol} Open Targets association with "
                        f"{disease_name or disease_id}"
                        + (f" (score={score})." if score is not None else ".")
                    ),
                    evidence_grade=EvidenceGrade.E,
                    subsection="Open Targets disease associations",
                    confidence_notes=(
                        "Open Targets association score is computational evidence "
                        "of association, not proof of causality."
                    ),
                    manual_review_required=True,
                )
            )

    for idx, item in enumerate(data.get("tractability_summaries") or [], start=1):
        if not isinstance(item, dict):
            continue
        modality = item.get("modality")
        label = item.get("label")
        if modality is None and label is None:
            continue
        key = f"tract-{modality or 'na'}-{label or idx}"
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.disease_association,
                fact_type="opentargets_tractability",
                key=key,
                value={
                    "ensembl_id": ensembl_id,
                    "modality": modality,
                    "label": label,
                    "value": item.get("value"),
                },
                display_text=(
                    f"{gene_symbol} Open Targets tractability"
                    + (f" {modality}" if modality else "")
                    + (f" / {label}" if label else "")
                    + (f"={item.get('value')}." if item.get("value") is not None else ".")
                ),
                evidence_grade=EvidenceGrade.E,
                subsection="Open Targets tractability",
                manual_review_required=True,
                confidence_notes=(
                    "Tractability annotation is platform metadata, not disease "
                    "causality evidence."
                ),
            )
        )

    probes = data.get("chemical_probes") or []
    if isinstance(probes, list):
        for idx, probe in enumerate(probes, start=1):
            if not isinstance(probe, dict):
                continue
            # Preserve reported probe fields without inventing activity claims.
            probe_id = (
                probe.get("id")
                or probe.get("chemicalProbeId")
                or probe.get("drugId")
                or probe.get("name")
            )
            key = f"probe-{probe_id or idx}"
            records.append(
                _record(
                    **common,
                    assertion_type=AssertionType.chemical_interaction,
                    fact_type="opentargets_chemical_probe",
                    key=key,
                    value={"ensembl_id": ensembl_id, "probe": probe},
                    display_text=(
                        f"{gene_symbol} Open Targets chemical probe"
                        + (f" {probe_id}." if probe_id else f" #{idx}.")
                    ),
                    evidence_grade=EvidenceGrade.E,
                    section=SECTION_CHEMICAL,
                    subsection="Open Targets chemical probes",
                    manual_review_required=True,
                    confidence_notes=(
                        "Chemical probe metadata from Open Targets; verify before "
                        "treating as a validated chemical tool."
                    ),
                )
            )
    return records


def _phenotype_map_rows(phenotype_map_list: Any) -> list[dict[str, Any]]:
    """Normalize OMIM phenotypeMapList items to dicts with a phenotype label."""
    if not isinstance(phenotype_map_list, list):
        return []
    out: list[dict[str, Any]] = []
    for item in phenotype_map_list:
        if not isinstance(item, dict):
            continue
        # Common shapes: {phenotypeMap: {...}} or flat phenotype fields.
        pmap = item.get("phenotypeMap") if isinstance(item.get("phenotypeMap"), dict) else item
        if not isinstance(pmap, dict):
            continue
        out.append(pmap)
    return out


def normalize_omim(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize OMIM ``fetch_gene_entry`` payloads (safe MIM selection only)."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    method = data.get("selection_method")
    selected_mim = data.get("selected_mim")
    if method == "ambiguous" or not selected_mim:
        return []

    summary = data.get("entry_summary")
    if not isinstance(summary, dict):
        return []

    caveat = str(
        data.get("caveat")
        or (
            "OMIM gene entries may lack disease/phenotype maps; do not invent "
            "phenotype relationships from absence of data."
        )
    )
    preferred_title = summary.get("preferred_title")
    records: list[EvidenceRecord] = [
        _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            source_name="OMIM",
            source_type=SourceType.genetic_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="omim_gene_entry",
            key=str(selected_mim),
            value={
                "mim_number": selected_mim,
                "selection_method": method,
                "preferred_title": preferred_title,
                "alternative_titles": summary.get("alternative_titles"),
                "chromosome": summary.get("chromosome"),
                "cyto_location": summary.get("cyto_location"),
                "gene_symbols": summary.get("gene_symbols"),
                "gene_name": summary.get("gene_name"),
                "phenotype_map_count": summary.get("phenotype_map_count"),
                "caveat": caveat,
            },
            display_text=(
                f"{gene_symbol} OMIM entry {selected_mim}"
                + (f" ({preferred_title})." if preferred_title else ".")
            ),
            evidence_grade=EvidenceGrade.C,
            subsection="OMIM gene entry",
            confidence_notes=caveat,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    ]

    for idx, pmap in enumerate(
        _phenotype_map_rows(summary.get("phenotype_map_list")), start=1
    ):
        phenotype = (
            pmap.get("phenotype")
            or pmap.get("phenotypeName")
            or pmap.get("phenotypeMapPhenotype")
        )
        phenotype_mim = (
            pmap.get("phenotypeMimNumber")
            or pmap.get("mimNumber")
            or pmap.get("phenotypeMapMimNumber")
        )
        if not phenotype and phenotype_mim is None:
            continue
        key = f"{selected_mim}-pheno-{phenotype_mim or idx}"
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="OMIM",
                source_type=SourceType.genetic_database,
                assertion_type=AssertionType.disease_association,
                fact_type="omim_phenotype_map",
                key=key,
                value={
                    "gene_mim_number": selected_mim,
                    "phenotype": phenotype,
                    "phenotype_mim_number": phenotype_mim,
                    "phenotype_map": pmap,
                    "caveat": (
                        "OMIM phenotype map entry as reported; do not invent "
                        "additional disease relationships."
                    ),
                },
                display_text=(
                    f"{gene_symbol} OMIM phenotype map"
                    + (f": {phenotype}" if phenotype else "")
                    + (f" (MIM {phenotype_mim})." if phenotype_mim else ".")
                ),
                evidence_grade=EvidenceGrade.C,
                subsection="OMIM phenotype map",
                confidence_notes=(
                    "OMIM phenotype map entry as reported; do not invent "
                    "additional disease relationships."
                ),
                manual_review_required=True,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_variants(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch variant/disease normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    if source == "ClinVar":
        return normalize_clinvar(tool_result, **kwargs)
    if source == "Open Targets":
        return normalize_opentargets(tool_result, **kwargs)
    if source == "OMIM":
        return normalize_omim(tool_result, **kwargs)
    return []


__all__ = [
    "DEFAULT_MAX_CLINVAR",
    "DEFAULT_MAX_OT_ASSOCIATIONS",
    "normalize_clinvar",
    "normalize_opentargets",
    "normalize_omim",
    "normalize_variants",
]
