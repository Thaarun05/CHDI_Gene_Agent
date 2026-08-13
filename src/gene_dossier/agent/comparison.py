"""Deterministic HD modifier analysis lens over qualifying EvidenceRecords."""

from __future__ import annotations

from gene_dossier.models import EvidenceRecord

from .capabilities import NEED_CONTRIBUTORS, qualifying_records, record_matches_need
from .models import (
    ComparisonCell,
    ComparisonRow,
    ComparisonStrength,
    EvidenceNeed,
    EvidenceRequirement,
    EvidenceRequirementAssessment,
    RequirementStatus,
)


HD_MODIFIER_DIMENSIONS: tuple[tuple[str, EvidenceNeed, bool, int], ...] = (
    ("Human Genetic Modifier Evidence", EvidenceNeed.human_genetic_association, True, 1),
    ("HD Literature", EvidenceNeed.hd_literature, True, 2),
    ("Repeat Instability / Mechanism", EvidenceNeed.repeat_instability_mechanism, True, 2),
    ("Experimental / Model", EvidenceNeed.experimental_evidence, False, 1),
    ("HD Expression", EvidenceNeed.brain_expression, False, 1),
    ("Pathway / PPI", EvidenceNeed.pathway_membership, False, 1),
    ("Therapeutic / Perturbability", EvidenceNeed.therapeutic_perturbability, False, 1),
)


def hd_modifier_requirements(genes: list[str]) -> list[EvidenceRequirement]:
    requirements: list[EvidenceRequirement] = []
    for index, (label, need, required, minimum) in enumerate(HD_MODIFIER_DIMENSIONS, start=1):
        requirements.append(
            EvidenceRequirement(
                id=f"hd_modifier_{index}",
                label=label,
                description=f"Provenance-backed evidence for the {label.lower()} comparison dimension.",
                genes=genes,
                evidence_need=need,
                capability_ids=list(NEED_CONTRIBUTORS[need]),
                required=required,
                minimum_support=minimum,
                rationale="Required by the deterministic HD modifier comparison lens.",
            )
        )
    return requirements


def _grade_value(record: EvidenceRecord) -> str:
    return str(getattr(record.evidence_grade, "value", record.evidence_grade)).upper()


def _record_text(record: EvidenceRecord) -> str:
    value = " ".join(f"{key} {val}" for key, val in (record.value or {}).items())
    return " ".join((record.fact_type, value, record.display_text)).lower()


def _assertion(record: EvidenceRecord) -> str:
    return str(getattr(record.assertion_type, "value", record.assertion_type)).lower()


def _source_type(record: EvidenceRecord) -> str:
    return str(getattr(record.source_type, "value", record.source_type)).lower()


def _is_sufficient(assessment: EvidenceRequirementAssessment | None) -> bool:
    return assessment is not None and assessment.status is RequirementStatus.sufficient


def grade_hd_modifier_cell(
    evidence_need: EvidenceNeed,
    matched_records: list[EvidenceRecord],
    assessment: EvidenceRequirementAssessment | None,
) -> ComparisonStrength:
    """Grade one HD rubric cell with dimension-specific, non-numeric rules."""
    if assessment and assessment.status in {
        RequirementStatus.missing,
        RequirementStatus.unsupported_capability,
    }:
        return "Missing"
    if not matched_records:
        return "Missing"

    sources = {record.source_name.strip().lower() for record in matched_records}
    if evidence_need is EvidenceNeed.human_genetic_association:
        direct = [
            record
            for record in matched_records
            if record_matches_need(record, EvidenceNeed.human_genetic_association)
            and _source_type(record) == "genetic_database"
            and any(
                term in _record_text(record)
                for term in ("modifier", "gwas", "age at onset", "somatic expansion")
            )
        ]
        direct_sources = {record.source_name.strip().lower() for record in direct}
        if (
            _is_sufficient(assessment)
            and len(direct) >= 2
            and len(direct_sources) >= 2
            and any(_grade_value(record) == "A" for record in direct)
        ):
            return "Strong"
        if direct:
            return "Moderate"
        return "Weak"

    if evidence_need is EvidenceNeed.hd_literature:
        distinct_papers = {record.source_id for record in matched_records}
        if _is_sufficient(assessment) and len(distinct_papers) >= 3:
            return "Strong"
        if len(distinct_papers) >= 2:
            return "Moderate"
        return "Limited"

    if evidence_need is EvidenceNeed.repeat_instability_mechanism:
        mechanistic = [
            record
            for record in matched_records
            if any(
                term in _record_text(record)
                for term in (
                    "repeat instability",
                    "repeat expansion",
                    "somatic expansion",
                    "mismatch repair",
                    "dna repair",
                    "repair mechanism",
                )
            )
        ]
        mechanistic_sources = {record.source_name.strip().lower() for record in mechanistic}
        if _is_sufficient(assessment) and len(mechanistic) >= 3 and len(mechanistic_sources) >= 2:
            return "Strong"
        if len(mechanistic) >= 2:
            return "Moderate"
        if mechanistic:
            return "Limited"
        return "Weak"

    if evidence_need is EvidenceNeed.experimental_evidence:
        experimental = [
            record
            for record in matched_records
            if _assertion(record) in {"perturbation", "knockout_phenotype"}
        ]
        experimental_sources = {record.source_name.strip().lower() for record in experimental}
        if _is_sufficient(assessment) and len(experimental) >= 3 and len(experimental_sources) >= 2:
            return "Strong"
        if len(experimental) >= 2:
            return "Moderate"
        if experimental:
            return "Limited"
        return "Weak"

    if evidence_need in {
        EvidenceNeed.brain_expression,
        EvidenceNeed.pathway_membership,
        EvidenceNeed.therapeutic_perturbability,
    }:
        # These dimensions provide context or tractability, not direct modifier validity.
        if len(matched_records) >= 2 or len(sources) >= 2:
            return "Moderate"
        return "Limited"

    if len(matched_records) >= 2:
        return "Moderate"
    return "Limited"


def build_hd_modifier_matrix(
    *,
    genes: list[str],
    requirements: list[EvidenceRequirement],
    assessments: list[EvidenceRequirementAssessment],
    records: list[EvidenceRecord],
) -> tuple[list[str], list[ComparisonRow]]:
    assessment_map = {
        (assessment.requirement_id, assessment.gene_symbol): assessment
        for assessment in assessments
    }
    rows: list[ComparisonRow] = []
    for requirement in requirements:
        cells: dict[str, ComparisonCell] = {}
        for gene in genes:
            gene_requirement = requirement.model_copy(update={"genes": [gene]})
            matched = qualifying_records(records, gene_requirement)
            assessment = assessment_map.get((requirement.id, gene))
            status = grade_hd_modifier_cell(requirement.evidence_need, matched, assessment)
            source_count = len({record.source_name for record in matched})
            summary = (
                f"{len(matched)} qualifying record(s) from {source_count} source(s)."
                if matched
                else "No qualifying provenance-backed evidence in the selected universe."
            )
            cells[gene] = ComparisonCell(
                status=status,
                summary=summary,
                evidence_count=len(matched),
                evidence_record_ids=[record.id for record in matched],
            )
        rows.append(ComparisonRow(dimension=requirement.label, cells=cells))
    return [requirement.label for requirement in requirements], rows


__all__ = [
    "HD_MODIFIER_DIMENSIONS",
    "build_hd_modifier_matrix",
    "grade_hd_modifier_cell",
    "hd_modifier_requirements",
]
