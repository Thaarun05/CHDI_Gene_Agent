"""Deterministic HD modifier analysis lens over qualifying EvidenceRecords."""

from __future__ import annotations

from gene_dossier.models import EvidenceRecord

from .capabilities import NEED_CONTRIBUTORS
from .evidence import canonicalize_requirement_evidence, public_evidence_reference
from .models import (
    ComparisonCell,
    ComparisonDecision,
    ComparisonDecisionOutcome,
    ComparisonRow,
    ComparisonStrength,
    EvidenceNeed,
    EvidenceRequirement,
    EvidenceRequirementAssessment,
    RequirementStatus,
    ScientificQuestionPlan,
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
            if _source_type(record) == "genetic_database"
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
        distinct_papers = {public_evidence_reference(record) for record in matched_records}
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
    plan: ScientificQuestionPlan | None = None,
    ordinal_by_id: dict[str, int] | None = None,
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
            canonical = canonicalize_requirement_evidence(
                records,
                gene_requirement,
                gene=gene,
                query_policy=plan.query_policy if plan else None,
                disease_contexts=plan.entities.diseases if plan else (),
            )
            groups = list(canonical.qualifying)
            matched = [group.canonical_record for group in groups]
            assessment = assessment_map.get((requirement.id, gene))
            status = grade_hd_modifier_cell(requirement.evidence_need, matched, assessment)
            source_count = len({group.source_namespace for group in groups})
            directions = {
                group.eligibility.direction for group in groups if group.eligibility.direction
            }
            summary = (
                f"{len(groups)} unique qualifying item(s) from {source_count} source(s)."
                if matched
                else "No qualifying provenance-backed evidence in the selected universe."
            )
            ordinals = [
                ordinal_by_id[group.canonical_record.id]
                for group in groups
                if ordinal_by_id and group.canonical_record.id in ordinal_by_id
            ][:3]
            cells[gene] = ComparisonCell(
                status=status,
                summary=summary,
                evidence_count=len(groups),
                evidence_record_ids=[record.id for record in matched],
                public_evidence_refs=[group.public_reference for group in groups],
                citation_ordinals=ordinals,
                distinct_source_count=source_count,
                direct_count=sum(
                    group.eligibility.designation.value == "direct" for group in groups
                ),
                supporting_count=sum(
                    group.eligibility.designation.value == "supporting" for group in groups
                ),
                excluded_count=len(canonical.contextual) + len(canonical.excluded),
                directionality_known=bool(directions),
                has_conflict="increase" in directions and "decrease" in directions,
            )
        rows.append(ComparisonRow(dimension=requirement.label, cells=cells))
    return [requirement.label for requirement in requirements], rows


def build_comparison_decision(
    *,
    plan: ScientificQuestionPlan,
    matrix: list[ComparisonRow],
    assessments: list[EvidenceRequirementAssessment],
) -> ComparisonDecision:
    """Apply a conservative, deterministic decision policy to a comparison."""
    if len(plan.entities.genes) < 2 or not matrix:
        return ComparisonDecision(
            outcome=ComparisonDecisionOutcome.not_rankable,
            summary="No multi-gene comparison decision is available.",
            limitations=[
                "A decision requires at least two genes and qualifying comparison dimensions."
            ],
        )

    requirement_by_need = {
        requirement.evidence_need: requirement for requirement in plan.evidence_requirements
    }
    row_by_label = {row.dimension: row for row in matrix}
    criteria = list(plan.query_policy.comparison_criteria)
    if not criteria:
        return ComparisonDecision(
            outcome=ComparisonDecisionOutcome.not_rankable,
            summary=(
                "The evidence supports a dimension-by-dimension comparison, but no explicit "
                "decision criterion supports an overall preference."
            ),
            limitations=[
                "No overall winner or numeric ranking was inferred from mixed evidence dimensions."
            ],
        )

    strength = {"Missing": 0, "Weak": 1, "Limited": 2, "Moderate": 3, "Strong": 4}
    criterion_winners: list[tuple[EvidenceNeed, str, ComparisonCell]] = []
    unsupported: list[str] = []
    for criterion in criteria:
        requirement = requirement_by_need.get(criterion)
        row = row_by_label.get(requirement.label) if requirement else None
        if row is None:
            unsupported.append(criterion.value)
            continue
        ranked = sorted(
            row.cells.items(),
            key=lambda item: (strength[item[1].status], item[1].evidence_count),
            reverse=True,
        )
        if len(ranked) < 2 or strength[ranked[0][1].status] == strength[ranked[1][1].status]:
            continue
        if ranked[0][1].status in {"Missing", "Weak"} or ranked[0][1].has_conflict:
            continue
        criterion_winners.append((criterion, ranked[0][0], ranked[0][1]))

    if not criterion_winners:
        return ComparisonDecision(
            outcome=ComparisonDecisionOutcome.not_rankable,
            summary="The requested comparison criteria do not support a defensible preference.",
            limitations=[
                *(
                    [f"Unsupported criteria: {', '.join(sorted(unsupported))}."]
                    if unsupported
                    else []
                ),
                "Missing, tied, conflicting, or indirect evidence was not converted into a winner.",
            ],
        )
    winners = {winner for _criterion, winner, _cell in criterion_winners}
    required_complete = all(
        assessment.status is RequirementStatus.sufficient
        for assessment in assessments
        if assessment.required
    )
    if not required_complete:
        if len(winners) == 1:
            gene = next(iter(winners))
            criteria_label = ", ".join(
                criterion.value for criterion, _winner, _cell in criterion_winners
            )
            return ComparisonDecision(
                outcome=ComparisonDecisionOutcome.dimension_specific_difference,
                summary=(
                    f"{gene} has stronger qualifying support in {criteria_label}, but missing "
                    "required evidence prevents an overall preference."
                ),
                criterion=criteria_label,
                limitations=[
                    "No unconditional winner is supported while required dimensions remain incomplete."
                ],
            )
        return ComparisonDecision(
            outcome=ComparisonDecisionOutcome.not_rankable,
            summary="Required evidence gaps prevent a defensible overall preference.",
            limitations=["No winner was inferred from incomplete required evidence."],
        )
    if len(winners) > 1:
        return ComparisonDecision(
            outcome=ComparisonDecisionOutcome.dimension_specific_difference,
            summary="Different genes are better supported on different requested dimensions.",
            limitations=["No overall winner is supported across the requested criteria."],
        )

    winner = next(iter(winners))
    complete = len(criterion_winners) == len(criteria) and not unsupported
    direction_required = {
        EvidenceNeed.experimental_evidence,
        EvidenceNeed.chemical_perturbation,
        EvidenceNeed.therapeutic_perturbability,
    }
    strongest = all(
        cell.status == "Strong"
        and cell.direct_count >= 2
        and cell.distinct_source_count >= 2
        and not cell.has_conflict
        and (criterion not in direction_required or cell.directionality_known)
        for criterion, _gene, cell in criterion_winners
    )
    no_required_conflict = True
    for criterion in criteria:
        requirement = requirement_by_need.get(criterion)
        row = row_by_label.get(requirement.label) if requirement else None
        if row is not None and any(cell.has_conflict for cell in row.cells.values()):
            no_required_conflict = False
            break
    outcome = (
        ComparisonDecisionOutcome.supported_preference
        if complete and strongest and no_required_conflict
        else ComparisonDecisionOutcome.conditional_preference
    )
    return ComparisonDecision(
        outcome=outcome,
        preferred_gene=winner,
        criterion=", ".join(criterion.value for criterion, _gene, _cell in criterion_winners),
        summary=(
            f"{winner} has stronger qualifying support for the explicit criterion set."
            if outcome is ComparisonDecisionOutcome.supported_preference
            else f"{winner} is conditionally better supported for the available requested criterion evidence."
        ),
        limitations=(
            []
            if outcome is ComparisonDecisionOutcome.supported_preference
            else ["This is criterion-specific and does not establish an overall biological winner."]
        ),
    )


__all__ = [
    "HD_MODIFIER_DIMENSIONS",
    "build_hd_modifier_matrix",
    "build_comparison_decision",
    "grade_hd_modifier_cell",
    "hd_modifier_requirements",
]
