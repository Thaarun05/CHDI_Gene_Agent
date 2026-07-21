"""Normalize chemical / bioassay ToolResults into EvidenceRecords.

Consumes successful client payloads from CTD, ChEMBL, and PubChem. Does **not**
call the network.

Rules:
- Do not invent chemicals, activities, targets, or assay claims
- ChEMBL target identity only when ``target_selection_method == "matched"``
- Ambiguous / unmatched ChEMBL targets do not invent a preferred target
- CTD rows are curated chemical–gene interactions (grade C)
- PubChem AID listings are gene-linked assay discovery (grade E); activity rows
  preserve reported outcomes only (also grade E, not curated chemical-tool proof)
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

SECTION_CTD = "CTD perturbations"
SECTION_CHEMICAL = "Chemical tools"

DEFAULT_MAX_CTD_RECORDS = 500
DEFAULT_MAX_CHEMBL_ACTIVITIES = 200
DEFAULT_MAX_PUBCHEM_AIDS = 100


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
    section: str,
    organism: str | None = None,
    taxon_id: int | None = None,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one chemical-related EvidenceRecord."""
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
        taxon_id=taxon_id,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def normalize_ctd(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_records: int = DEFAULT_MAX_CTD_RECORDS,
) -> list[EvidenceRecord]:
    """Normalize CTD ``fetch_chemical_gene_interactions`` summaries."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    summaries = data.get("interaction_summaries") or []
    if not isinstance(summaries, list):
        summaries = []

    records: list[EvidenceRecord] = []
    for idx, row in enumerate(summaries[: max(0, max_records)], start=1):
        if not isinstance(row, dict):
            continue
        chemical_name = row.get("chemical_name")
        chemical_id = row.get("chemical_id")
        if not chemical_name and not chemical_id:
            continue
        interaction = row.get("interaction")
        organism = row.get("organism")
        organism_id = row.get("organism_id")
        try:
            taxon_id = int(organism_id) if organism_id not in (None, "") else None
        except (TypeError, ValueError):
            taxon_id = None
        key = "-".join(
            str(x)
            for x in (chemical_id or chemical_name, row.get("gene_id"), idx)
            if x is not None and str(x) != ""
        )
        display = (
            f"{gene_symbol} CTD interaction with "
            f"{chemical_name or chemical_id}"
        )
        if interaction:
            display += f": {interaction}"
        else:
            display += "."

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="CTD",
                source_type=SourceType.chemical_database,
                assertion_type=AssertionType.chemical_interaction,
                fact_type="ctd_chemical_gene_interaction",
                key=key,
                value={
                    "chemical_name": chemical_name,
                    "chemical_id": chemical_id,
                    "cas_rn": row.get("cas_rn"),
                    "gene_symbol_reported": row.get("gene_symbol"),
                    "gene_id": row.get("gene_id"),
                    "organism": organism,
                    "organism_id": organism_id,
                    "interaction": interaction,
                    "interaction_actions": row.get("interaction_actions"),
                    "pubmed_ids": row.get("pubmed_ids"),
                },
                display_text=display if display.endswith(".") else display + ".",
                evidence_grade=EvidenceGrade.C,
                section=SECTION_CTD,
                subsection="CTD chemical–gene interactions",
                organism=str(organism) if organism else None,
                taxon_id=taxon_id,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_chembl(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_activities: int = DEFAULT_MAX_CHEMBL_ACTIVITIES,
) -> list[EvidenceRecord]:
    """Normalize ChEMBL ``fetch_chemical_tools`` payloads.

    Emits a target identity record only when selection is ``matched``.
    Assay/activity rows are still emitted when present, with weaker grading
    when the target was not safely matched.
    """
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    method = str(data.get("target_selection_method") or "")
    target_id = data.get("target_chembl_id")
    target_matched = method == "matched" and bool(target_id)
    weaker = not target_matched
    grade = EvidenceGrade.C if target_matched else EvidenceGrade.E
    review = weaker

    records: list[EvidenceRecord] = []
    common = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "source_name": "ChEMBL",
        "source_type": SourceType.chemical_database,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }

    if target_matched:
        # Prefer summary for the selected target when available.
        pref_name = None
        organism = None
        for row in data.get("target_summaries") or []:
            if isinstance(row, dict) and str(row.get("target_chembl_id")) == str(
                target_id
            ):
                pref_name = row.get("pref_name")
                organism = row.get("organism")
                break
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.chemical_interaction,
                fact_type="chembl_target",
                key=str(target_id),
                value={
                    "target_chembl_id": target_id,
                    "pref_name": pref_name,
                    "organism": organism,
                    "target_selection_method": method,
                },
                display_text=(
                    f"{gene_symbol} ChEMBL target {target_id}"
                    + (f" ({pref_name})." if pref_name else ".")
                ),
                evidence_grade=EvidenceGrade.C,
                section=SECTION_CHEMICAL,
                subsection="ChEMBL target",
                organism=str(organism) if organism else None,
            )
        )

    for assay in data.get("assay_summaries") or []:
        if not isinstance(assay, dict):
            continue
        assay_id = assay.get("assay_chembl_id")
        if not assay_id:
            continue
        desc = assay.get("description")
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.chemical_interaction,
                fact_type="chembl_assay",
                key=str(assay_id),
                value={
                    "assay_chembl_id": assay_id,
                    "description": desc,
                    "assay_type": assay.get("assay_type"),
                    "assay_organism": assay.get("assay_organism"),
                    "assay_cell_type": assay.get("assay_cell_type"),
                    "target_chembl_id": assay.get("target_chembl_id") or target_id,
                    "document_chembl_id": assay.get("document_chembl_id"),
                    "target_selection_method": method,
                    "assay_terms": data.get("assay_terms"),
                },
                display_text=(
                    f"{gene_symbol} ChEMBL assay {assay_id}"
                    + (f": {desc}" if desc else ".")
                ),
                evidence_grade=grade,
                section=SECTION_CHEMICAL,
                subsection="ChEMBL assays",
                confidence_notes=(
                    None
                    if target_matched
                    else (
                        f"ChEMBL target selection was {method or 'unresolved'}; "
                        "assay linked via gene-specific search terms, not a "
                        "confirmed matched target."
                    )
                ),
                manual_review_required=review,
            )
        )

    activities = data.get("activity_summaries") or []
    if not isinstance(activities, list):
        activities = []
    for idx, act in enumerate(activities[: max(0, max_activities)], start=1):
        if not isinstance(act, dict):
            continue
        mol = act.get("molecule_chembl_id")
        assay_id = act.get("assay_chembl_id")
        if not mol and not assay_id:
            continue
        key = "-".join(
            str(x)
            for x in (mol or "mol", assay_id or "assay", act.get("standard_type"), idx)
            if x is not None and str(x) != ""
        )
        std_type = act.get("standard_type")
        std_value = act.get("standard_value")
        std_units = act.get("standard_units")
        display = f"{gene_symbol} ChEMBL activity"
        if mol:
            display += f" for {mol}"
        if std_type is not None:
            display += f" {std_type}"
            if act.get("standard_relation"):
                display += f" {act['standard_relation']}"
            if std_value is not None:
                display += f" {std_value}"
            if std_units:
                display += f" {std_units}"
        display += "."
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.chemical_interaction,
                fact_type="chembl_activity",
                key=key,
                value={
                    "molecule_chembl_id": mol,
                    "canonical_smiles": act.get("canonical_smiles"),
                    "standard_type": std_type,
                    "standard_relation": act.get("standard_relation"),
                    "standard_value": std_value,
                    "standard_units": std_units,
                    "pchembl_value": act.get("pchembl_value"),
                    "assay_chembl_id": assay_id,
                    "document_chembl_id": act.get("document_chembl_id"),
                    "target_chembl_id": target_id,
                    "target_selection_method": method,
                },
                display_text=display,
                evidence_grade=grade,
                section=SECTION_CHEMICAL,
                subsection="ChEMBL activities",
                confidence_notes=(
                    None
                    if target_matched
                    else (
                        f"ChEMBL target selection was {method or 'unresolved'}; "
                        "review before treating as on-target chemical evidence."
                    )
                ),
                manual_review_required=review,
            )
        )
    return records


def normalize_pubchem(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_aids: int = DEFAULT_MAX_PUBCHEM_AIDS,
) -> list[EvidenceRecord]:
    """Normalize PubChem ``fetch_bioassays`` payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    gene_id = data.get("gene_id")
    aids = data.get("aids") or []
    if not isinstance(aids, list):
        aids = []

    records: list[EvidenceRecord] = []
    common = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "source_name": "PubChem",
        "source_type": SourceType.chemical_database,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }

    if aids:
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.chemical_interaction,
                fact_type="pubchem_aid_summary",
                key=f"geneid-{gene_id or gene_symbol}-aids",
                value={
                    "gene_id": gene_id,
                    "aid_count": data.get("aid_count") or len(aids),
                    "aids_sample": [str(a) for a in aids[:max_aids]],
                },
                display_text=(
                    f"{gene_symbol} has {data.get('aid_count') or len(aids)} "
                    f"PubChem BioAssay AID(s) for Gene ID {gene_id}."
                ),
                evidence_grade=EvidenceGrade.E,
                section=SECTION_CHEMICAL,
                subsection="PubChem BioAssay AIDs",
                confidence_notes=(
                    "PubChem AID list is gene-linked assay discovery; not proof "
                    "of chemical tool quality or on-target activity."
                ),
                manual_review_required=True,
            )
        )

    for summary in data.get("description_summaries") or []:
        if not isinstance(summary, dict):
            continue
        aid = summary.get("aid")
        if aid is None:
            continue
        name = summary.get("name")
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.chemical_interaction,
                fact_type="pubchem_assay_description",
                key=f"aid-{aid}",
                value={
                    "gene_id": gene_id,
                    "aid": aid,
                    "name": name,
                    "comment": summary.get("comment"),
                },
                display_text=(
                    f"{gene_symbol} PubChem assay AID {aid}"
                    + (f": {name}." if name else ".")
                ),
                evidence_grade=EvidenceGrade.E,
                section=SECTION_CHEMICAL,
                subsection="PubChem assay descriptions",
                manual_review_required=True,
                confidence_notes=(
                    "Assay description metadata only; does not invent activity "
                    "outcomes."
                ),
            )
        )

    assay_csv = data.get("assay_csv") or {}
    if isinstance(assay_csv, dict):
        for aid_key, payload in assay_csv.items():
            if not isinstance(payload, dict):
                continue
            for idx, act in enumerate(payload.get("activity_summaries") or [], start=1):
                if not isinstance(act, dict):
                    continue
                cid = act.get("pubchem_cid")
                outcome = act.get("activity_outcome")
                if cid is None and outcome is None and act.get("standard_value") is None:
                    continue
                key = f"aid-{aid_key}-cid-{cid or 'na'}-{idx}"
                display = f"{gene_symbol} PubChem AID {aid_key} activity"
                if cid is not None:
                    display += f" CID {cid}"
                if outcome:
                    display += f" outcome={outcome}"
                if act.get("standard_type") is not None:
                    display += f" {act['standard_type']}"
                    if act.get("standard_value") is not None:
                        display += f"={act['standard_value']}"
                    if act.get("standard_units"):
                        display += f" {act['standard_units']}"
                display += "."
                records.append(
                    _record(
                        **common,
                        assertion_type=AssertionType.chemical_interaction,
                        fact_type="pubchem_activity",
                        key=key,
                        value={
                            "gene_id": gene_id,
                            "aid": aid_key,
                            "pubchem_cid": cid,
                            "activity_outcome": outcome,
                            "standard_type": act.get("standard_type"),
                            "standard_value": act.get("standard_value"),
                            "standard_units": act.get("standard_units"),
                            "caveat": (
                                "PubChem activity rows preserve reported assay "
                                "outcomes only and are not proof of chemical tool "
                                "quality or on-target gene activity."
                            ),
                        },
                        display_text=display,
                        evidence_grade=EvidenceGrade.E,
                        section=SECTION_CHEMICAL,
                        subsection="PubChem assay activities",
                        confidence_notes=(
                            "PubChem activity rows preserve reported assay outcomes "
                            "only and are not proof of chemical tool quality or "
                            "on-target gene activity."
                        ),
                        manual_review_required=True,
                    )
                )
    return records


def normalize_chemicals(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch chemical normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    if source == "CTD":
        return normalize_ctd(tool_result, **kwargs)
    if source == "ChEMBL":
        return normalize_chembl(tool_result, **kwargs)
    if source == "PubChem":
        return normalize_pubchem(tool_result, **kwargs)
    return []


__all__ = [
    "DEFAULT_MAX_CTD_RECORDS",
    "DEFAULT_MAX_CHEMBL_ACTIVITIES",
    "DEFAULT_MAX_PUBCHEM_AIDS",
    "normalize_ctd",
    "normalize_chembl",
    "normalize_pubchem",
    "normalize_chemicals",
]
