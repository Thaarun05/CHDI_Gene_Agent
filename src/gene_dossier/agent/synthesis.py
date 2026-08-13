"""Server-owned grounding slots with optional LLM prose fragments."""

from __future__ import annotations

import logging
import re
from collections import Counter
from decimal import ROUND_HALF_UP, Decimal
from enum import Enum
from typing import Iterable

from pydantic import BaseModel, ConfigDict, Field

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import EvidenceRecord

from .models import (
    AnswerStatus,
    CitationReference,
    ComparisonRow,
    EvidenceGap,
    EvidenceNeed,
    EvidenceRequirementAssessment,
    ExperimentRecommendation,
    RequirementStatus,
    ScientificQuestionPlan,
)

logger = logging.getLogger(__name__)

_CITATION_RE = re.compile(r"\[\[(\d+)\]\]")
_ANY_CITATION_MARKER_RE = re.compile(r"\[\[|\]\]|(?<!\[)\[\s*-?\d+\s*\](?!\])")
_MALFORMED_CITATION_RE = re.compile(r"(?<!\[)\[\d+\](?!\])|\[\[[^\]\d][^\]]*\]\]|\[\[\d+\](?!\])")
_RANKING_RE = re.compile(
    r"\b(winner|ranked first|best overall|score|outperforms|superior)\b", re.IGNORECASE
)
_CAUSAL_RE = re.compile(
    r"\b(cause[sd]?|causal(?:ly)?|drive[sn]?|led to|leads to|result(?:s|ed)? in|prove[sd]?|validate[sd]?)\b",
    re.IGNORECASE,
)
_DIRECTIONAL_RE = re.compile(
    r"\b(increase[sd]?|decrease[sd]?|reduce[sd]?|enhance[sd]?|suppress(?:es|ed)?|up-?regulat(?:e[sd]?|ion)|down-?regulat(?:e[sd]?|ion)|gain|loss)\b",
    re.IGNORECASE,
)
_CONFLICT_RE = re.compile(
    r"\b(conflict(?:s|ing)?|contradict(?:s|ory)?|disagree(?:s|ment)?|argues? against)\b",
    re.IGNORECASE,
)
_PREDICTED_RESULT_RE = re.compile(
    r"\b(will|would)\s+(validate|prove|confirm|demonstrate|show|establish)\b",
    re.IGNORECASE,
)

_GRADE_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}
_STATUS_ORDER = {
    RequirementStatus.sufficient: 0,
    RequirementStatus.limited: 1,
    RequirementStatus.missing: 2,
    RequirementStatus.unsupported_capability: 3,
}


class EpistemicRole(str, Enum):
    direct_evidence = "direct_evidence"
    supporting_evidence = "supporting_evidence"
    mechanistic_inference = "mechanistic_inference"
    hypothesis = "hypothesis"


class ProseSection(str, Enum):
    direct_answer = "direct_answer"
    evidence = "evidence"
    comparison = "comparison"
    recommendation = "recommendation"


class ClaimLanguagePolicy(str, Enum):
    descriptive_only = "descriptive_only"
    association_only = "association_only"
    mechanistic_support = "mechanistic_support"
    directional_supported = "directional_supported"
    causal_supported = "causal_supported"
    hypothesis_only = "hypothesis_only"
    recommendation = "recommendation"


class GroundedProseSlot(BaseModel):
    """Private server-owned grounding contract for one prose fragment."""

    model_config = ConfigDict(extra="forbid")

    slot_id: str
    section: ProseSection
    gene_symbols: tuple[str, ...]
    evidence_category: EvidenceNeed | None = None
    epistemic_role: EpistemicRole
    citation_ordinals: tuple[int, ...] = ()
    language_policy: ClaimLanguagePolicy
    evidence_text: tuple[str, ...] = ()
    public_sources: tuple[str, ...] = ()
    public_source_ids: tuple[str, ...] = ()
    record_ids: tuple[str, ...] = ()
    fallback_text: str
    stable_order: int
    heading_label: str
    allow_directional_language: bool = False
    allow_causal_language: bool = False
    allow_conflict_language: bool = False
    recommendation_index: int | None = None


class GroundedProseFragment(BaseModel):
    """The complete model-controlled surface for one prose slot."""

    model_config = ConfigDict(extra="forbid")

    slot_id: str
    text: str


class GroundedProseDraft(BaseModel):
    model_config = ConfigDict(extra="forbid")

    fragments: list[GroundedProseFragment]


class GroundingValidationIssue(BaseModel):
    code: str
    answer_section: str
    rule: str
    explanation: str
    slot_id: str | None = None
    blocking: bool = True


class GroundedSynthesisResult(BaseModel):
    summary: str
    generation_method: str
    recommendations: list[ExperimentRecommendation] = Field(default_factory=list)
    requested_slot_count: int = 0
    accepted_slot_count: int = 0
    fallback_slot_count: int = 0
    diagnostic_counts: dict[str, int] = Field(default_factory=dict)
    failure_type: str | None = None
    failure_message: str | None = None
    validation_issues: list[GroundingValidationIssue] = Field(default_factory=list)


def _assertion(record: EvidenceRecord) -> str:
    return str(getattr(record.assertion_type, "value", record.assertion_type)).lower()


def _grade(record: EvidenceRecord) -> str:
    return str(getattr(record.evidence_grade, "value", record.evidence_grade)).upper()


def _record_has_direction_metadata(record: EvidenceRecord) -> bool:
    for key, value in (record.value or {}).items():
        key_text = str(key).lower()
        if any(
            term in key_text
            for term in ("direction", "effect", "change", "modulation", "regulation")
        ) and value not in (None, "", []):
            return True
    return False


def _exact_private_values(records: Iterable[EvidenceRecord]) -> set[str]:
    private: set[str] = set()
    for record in records:
        for value in (
            record.id,
            record.dossier_run_id,
            record.api_run_id,
            record.raw_artifact_id,
            record.raw_response_pointer,
        ):
            if value:
                private.add(str(value))
    return private


def _remove_private_values(text: str, private_values: Iterable[str]) -> str:
    cleaned = text
    for value in sorted({item for item in private_values if item}, key=len, reverse=True):
        cleaned = cleaned.replace(value, "stored evidence")
    return cleaned.strip()


def _contains_actual_record_id(text: str, record_ids: Iterable[str]) -> bool:
    return any(record_id and record_id in text for record_id in record_ids)


def _policy_for_records(records: list[EvidenceRecord]) -> ClaimLanguagePolicy:
    assertions = {_assertion(record) for record in records}
    if assertions & {"perturbation", "knockout_phenotype"}:
        return ClaimLanguagePolicy.causal_supported
    if any(_record_has_direction_metadata(record) for record in records):
        return ClaimLanguagePolicy.directional_supported
    if assertions & {
        "variant_association",
        "disease_association",
        "ppi",
        "expression",
        "cell_type_expression",
        "transcription_factor_association",
        "chemical_interaction",
        "chemical_tool",
    }:
        return ClaimLanguagePolicy.association_only
    if assertions & {"protein_function", "pathway_membership", "literature_summary"}:
        return ClaimLanguagePolicy.mechanistic_support
    return ClaimLanguagePolicy.descriptive_only


def _role_for_records(
    records: list[EvidenceRecord], category: EvidenceNeed | None
) -> EpistemicRole:
    assertions = {_assertion(record) for record in records}
    if assertions & {
        "variant_association",
        "disease_association",
        "literature_summary",
        "perturbation",
        "knockout_phenotype",
        "ppi",
        "pathway_membership",
        "expression",
        "cell_type_expression",
    }:
        return EpistemicRole.direct_evidence
    if category is EvidenceNeed.repeat_instability_mechanism:
        return EpistemicRole.mechanistic_inference
    return EpistemicRole.supporting_evidence


def _record_sort_key(
    record: EvidenceRecord, input_order: dict[str, int]
) -> tuple[int, int, str, str]:
    return (
        _GRADE_ORDER.get(_grade(record), 99),
        input_order.get(record.id, 9999),
        record.source_name.casefold(),
        record.source_id,
    )


def _fallback_for_evidence_slot(
    *,
    genes: tuple[str, ...],
    category: EvidenceNeed,
    evidence_text: tuple[str, ...],
) -> str:
    gene_label = " and ".join(genes)
    category_label = category.value.replace("_", " ")
    if not evidence_text:
        return f"No qualifying {category_label} evidence was available for {gene_label}."
    joined = "; ".join(text.rstrip(". ") for text in evidence_text[:3])
    return f"For {gene_label}, the selected {category_label} evidence reports: {joined}."


def _slot_from_records(
    *,
    slot_id: str,
    section: ProseSection,
    genes: tuple[str, ...],
    category: EvidenceNeed,
    records: list[EvidenceRecord],
    ordinal_by_id: dict[str, int],
    private_values: set[str],
    stable_order: int,
    heading_label: str,
) -> GroundedProseSlot:
    evidence_text = tuple(
        text
        for record in records
        if (text := _remove_private_values(record.display_text, private_values))
    )
    policy = _policy_for_records(records)
    return GroundedProseSlot(
        slot_id=slot_id,
        section=section,
        gene_symbols=genes,
        evidence_category=category,
        epistemic_role=_role_for_records(records, category),
        citation_ordinals=tuple(
            ordinal_by_id[record.id] for record in records if record.id in ordinal_by_id
        ),
        language_policy=policy,
        evidence_text=evidence_text,
        public_sources=tuple(dict.fromkeys(record.source_name for record in records)),
        public_source_ids=tuple(
            dict.fromkeys(
                record.source_id
                for record in records
                if record.source_id and record.source_id not in private_values
            )
        ),
        record_ids=tuple(record.id for record in records),
        fallback_text=_fallback_for_evidence_slot(
            genes=genes,
            category=category,
            evidence_text=evidence_text,
        ),
        stable_order=stable_order,
        heading_label=heading_label,
        allow_directional_language=any(
            _record_has_direction_metadata(record) for record in records
        ),
        allow_causal_language=any(
            _assertion(record) in {"perturbation", "knockout_phenotype"} for record in records
        ),
    )


def build_grounded_prose_slots(
    *,
    plan: ScientificQuestionPlan,
    records: list[EvidenceRecord],
    assessments: list[EvidenceRequirementAssessment],
    recommendations: list[ExperimentRecommendation] | None = None,
    comparison_matrix: list[ComparisonRow] | None = None,
) -> list[GroundedProseSlot]:
    """Build stable, server-owned slots from deterministic evidence structures."""
    ordinal_by_id = {record.id: ordinal for ordinal, record in enumerate(records, start=1)}
    record_by_id = {record.id: record for record in records}
    input_order = {record.id: index for index, record in enumerate(records)}
    private_values = _exact_private_values(records)
    plan_gene_order = {gene: index for index, gene in enumerate(plan.entities.genes)}
    requirement_order = {item.id: index for index, item in enumerate(plan.evidence_requirements)}
    ordered_assessments = sorted(
        assessments,
        key=lambda item: (
            not item.required,
            plan_gene_order.get(item.gene_symbol, 999),
            requirement_order.get(item.requirement_id, 999),
            _STATUS_ORDER.get(item.status, 999),
            item.evidence_need.value,
        ),
    )

    assessment_records: list[tuple[EvidenceRequirementAssessment, list[EvidenceRecord]]] = []
    for assessment in ordered_assessments:
        matched = [
            record_by_id[item] for item in assessment.evidence_record_ids if item in record_by_id
        ]
        matched.sort(key=lambda record: _record_sort_key(record, input_order))
        if matched:
            assessment_records.append((assessment, matched))

    slots: list[GroundedProseSlot] = []
    stable_order = 0
    for gene in plan.entities.genes:
        candidates = [item for item in assessment_records if item[0].gene_symbol == gene]
        if not candidates:
            continue
        assessment, matched = min(
            candidates,
            key=lambda item: (
                not item[0].required,
                _STATUS_ORDER.get(item[0].status, 999),
                min((_GRADE_ORDER.get(_grade(record), 99) for record in item[1]), default=99),
                requirement_order.get(item[0].requirement_id, 999),
            ),
        )
        stable_order += 1
        slots.append(
            _slot_from_records(
                slot_id=f"slot_{stable_order:03d}",
                section=ProseSection.direct_answer,
                genes=(gene,),
                category=assessment.evidence_need,
                records=matched[:3],
                ordinal_by_id=ordinal_by_id,
                private_values=private_values,
                stable_order=stable_order,
                heading_label="Direct answer",
            )
        )

    matrix = comparison_matrix or []
    if matrix:
        requirement_by_label = {item.label: item for item in plan.evidence_requirements}
        for row in matrix:
            requirement = requirement_by_label.get(row.dimension)
            if requirement is None:
                continue
            for gene in plan.entities.genes:
                cell = row.cells.get(gene)
                if cell is None:
                    continue
                matched = [
                    record_by_id[item] for item in cell.evidence_record_ids if item in record_by_id
                ]
                matched.sort(key=lambda record: _record_sort_key(record, input_order))
                if not matched:
                    continue
                stable_order += 1
                slots.append(
                    _slot_from_records(
                        slot_id=f"slot_{stable_order:03d}",
                        section=ProseSection.comparison,
                        genes=(gene,),
                        category=requirement.evidence_need,
                        records=matched[:3],
                        ordinal_by_id=ordinal_by_id,
                        private_values=private_values,
                        stable_order=stable_order,
                        heading_label=f"{row.dimension} · {gene} · {cell.status}",
                    )
                )
    else:
        for assessment, matched in assessment_records:
            stable_order += 1
            slots.append(
                _slot_from_records(
                    slot_id=f"slot_{stable_order:03d}",
                    section=ProseSection.evidence,
                    genes=(assessment.gene_symbol,),
                    category=assessment.evidence_need,
                    records=matched[:3],
                    ordinal_by_id=ordinal_by_id,
                    private_values=private_values,
                    stable_order=stable_order,
                    heading_label=f"{assessment.gene_symbol} · {assessment.evidence_need.value}",
                )
            )

    for index, recommendation in enumerate(recommendations or []):
        gap_category = next(
            (
                assessment.evidence_need
                for assessment in ordered_assessments
                if recommendation.gap_ids
                and any(
                    gap_id.startswith(
                        f"{assessment.gene_symbol.lower()}:{assessment.evidence_need.value}:"
                    )
                    for gap_id in recommendation.gap_ids
                )
            ),
            None,
        )
        rationale_records = [
            records[ordinal - 1]
            for ordinal in recommendation.rationale_citation_ordinals
            if 1 <= ordinal <= len(records)
        ]
        genes = tuple(
            dict.fromkeys(
                assessment.gene_symbol
                for assessment in ordered_assessments
                if recommendation.gap_ids
                and any(
                    gap_id.startswith(f"{assessment.gene_symbol.lower()}:")
                    for gap_id in recommendation.gap_ids
                )
            )
        ) or tuple(plan.entities.genes[:1])
        stable_order += 1
        fallback = recommendation.description.removeprefix("Recommendation:").strip()
        slots.append(
            GroundedProseSlot(
                slot_id=f"slot_{stable_order:03d}",
                section=ProseSection.recommendation,
                gene_symbols=genes,
                evidence_category=gap_category,
                epistemic_role=EpistemicRole.hypothesis,
                citation_ordinals=tuple(recommendation.rationale_citation_ordinals),
                language_policy=ClaimLanguagePolicy.recommendation,
                evidence_text=tuple(
                    _remove_private_values(record.display_text, private_values)
                    for record in rationale_records
                ),
                public_sources=tuple(
                    dict.fromkeys(record.source_name for record in rationale_records)
                ),
                public_source_ids=tuple(
                    dict.fromkeys(
                        record.source_id
                        for record in rationale_records
                        if record.source_id and record.source_id not in private_values
                    )
                ),
                record_ids=tuple(record.id for record in rationale_records),
                fallback_text=fallback,
                stable_order=stable_order,
                heading_label="Recommendation",
                recommendation_index=index,
            )
        )
    return slots


def _policy_instruction(slot: GroundedProseSlot) -> str:
    instructions = {
        ClaimLanguagePolicy.descriptive_only: "Describe only what the supplied evidence reports.",
        ClaimLanguagePolicy.association_only: "Use association language only; do not imply mechanism or causality.",
        ClaimLanguagePolicy.mechanistic_support: "Describe mechanistic support without asserting causal proof or direction.",
        ClaimLanguagePolicy.directional_supported: "Directional wording is allowed only for the supplied structured effect.",
        ClaimLanguagePolicy.causal_supported: "Causal wording is allowed only to the extent explicitly supported by the supplied perturbational evidence.",
        ClaimLanguagePolicy.hypothesis_only: "State only a bounded hypothesis or uncertainty.",
        ClaimLanguagePolicy.recommendation: "Phrase a bounded experiment; do not predict its result or claim it will validate the target.",
    }
    return instructions[slot.language_policy]


def build_grounded_prose_prompt(
    *,
    question: str,
    status: AnswerStatus,
    slots: list[GroundedProseSlot],
    private_values: Iterable[str] = (),
) -> str:
    """Render a prompt containing no private provenance identifiers."""
    slot_blocks: list[str] = []
    for slot in sorted(slots, key=lambda item: item.stable_order):
        sources = ", ".join(slot.public_sources) or "none"
        source_ids = ", ".join(slot.public_source_ids) or "none"
        evidence = (
            "\n".join(f"- {text}" for text in slot.evidence_text)
            or "- No compatible rationale record; this is a gap-driven recommendation."
        )
        slot_blocks.append(
            "\n".join(
                (
                    f"slot_id: {slot.slot_id}",
                    f"section: {slot.section.value}",
                    f"genes: {', '.join(slot.gene_symbols)}",
                    f"evidence_category: {slot.evidence_category.value if slot.evidence_category else 'none'}",
                    f"epistemic_role: {slot.epistemic_role.value}",
                    f"language_policy: {slot.language_policy.value}",
                    f"public_sources: {sources}",
                    f"public_source_ids: {source_ids}",
                    f"instruction: {_policy_instruction(slot)}",
                    "evidence_text:",
                    evidence,
                )
            )
        )
    safe_question = _remove_private_values(question, private_values)
    return f"""You are filling server-defined prose slots.

For every supplied slot, return only:
- slot_id
- one concise prose fragment

Do not create headings, citations, citation markers, categories, genes, evidence labels,
gap identifiers, recommendation labels, comparison scores, or provenance.
Do not use outside knowledge. Do not include internal identifiers.
Follow each slot's language policy exactly.
If the supplied evidence does not support a useful statement, return a brief uncertainty
statement for that slot without adding external facts.

Question: {safe_question}
Scientific status: {status.value}

Slots:
{chr(10).join(chr(10) + block for block in slot_blocks)}
"""


def _fragment_issue(
    slot: GroundedProseSlot, code: str, rule: str, explanation: str
) -> GroundingValidationIssue:
    return GroundingValidationIssue(
        code=code,
        answer_section=slot.section.value,
        slot_id=slot.slot_id,
        rule=rule,
        explanation=explanation,
    )


def validate_prose_fragment(
    fragment: GroundedProseFragment,
    *,
    slot: GroundedProseSlot,
    all_plan_genes: tuple[str, ...],
    actual_record_ids: set[str],
) -> list[GroundingValidationIssue]:
    """Validate scientific wording only within its assigned prose slot."""
    text = fragment.text.strip()
    if not text:
        return [
            _fragment_issue(
                slot,
                "empty_fragment",
                "slot_fragment_must_be_nonempty",
                "The prose fragment was empty.",
            )
        ]
    issues: list[GroundingValidationIssue] = []
    if _ANY_CITATION_MARKER_RE.search(text):
        issues.append(
            _fragment_issue(
                slot,
                "model_generated_citation_marker",
                "citations_are_server_rendered",
                "The model attempted to generate citation syntax.",
            )
        )
    if _contains_actual_record_id(text, actual_record_ids):
        issues.append(
            _fragment_issue(
                slot,
                "private_evidence_record_id_exposed",
                "private_record_ids_must_not_be_rendered",
                "The fragment exposed an actual internal EvidenceRecord identifier.",
            )
        )
    assigned = set(slot.gene_symbols)
    mentioned = {
        gene
        for gene in all_plan_genes
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])", text, re.IGNORECASE)
    }
    if mentioned and not mentioned & assigned:
        issues.append(
            _fragment_issue(
                slot,
                "wrong_assigned_gene",
                "fragment_must_address_its_server_assigned_gene",
                "The fragment addressed a different resolved gene than its server-assigned slot.",
            )
        )
    if _RANKING_RE.search(text):
        issues.append(
            _fragment_issue(
                slot,
                "comparison_ranking_or_winner_language",
                "comparisons_must_not_rank_or_declare_winners",
                "The fragment introduced ranking or winner language.",
            )
        )
    if _CAUSAL_RE.search(text) and not slot.allow_causal_language:
        issues.append(
            _fragment_issue(
                slot,
                "unsupported_causal_language",
                "causal_language_requires_server_authorized_perturbational_evidence",
                "The fragment used causal language outside its assigned policy.",
            )
        )
    if _DIRECTIONAL_RE.search(text) and not slot.allow_directional_language:
        issues.append(
            _fragment_issue(
                slot,
                "unsupported_directional_language",
                "directional_language_requires_server_authorized_direction_metadata",
                "The fragment used directional language outside its assigned policy.",
            )
        )
    if _CONFLICT_RE.search(text) and not slot.allow_conflict_language:
        issues.append(
            _fragment_issue(
                slot,
                "unauthorized_conflict_language",
                "conflict_language_requires_deterministic_conflict_status",
                "The fragment asserted a conflict not authorized by deterministic analysis.",
            )
        )
    if slot.language_policy is ClaimLanguagePolicy.recommendation and _PREDICTED_RESULT_RE.search(
        text
    ):
        issues.append(
            _fragment_issue(
                slot,
                "recommendation_presented_as_established_fact",
                "recommendations_must_not_predict_results",
                "The recommendation predicted a result or validation outcome.",
            )
        )
    return issues


def _render_slot(slot: GroundedProseSlot, text: str) -> str:
    citations = " ".join(f"[[{ordinal}]]" for ordinal in slot.citation_ordinals)
    return f"{text.strip()} {citations}".strip()


def validate_rendered_answer(
    answer: str,
    *,
    citation_registry: list[CitationReference],
    actual_record_ids: set[str],
    rendered_slot_ids: list[str],
    expected_slot_ids: list[str],
) -> list[GroundingValidationIssue]:
    """Validate only final citation and structural integrity, never scientific wording."""
    issues: list[GroundingValidationIssue] = []
    if not answer.strip():
        issues.append(
            GroundingValidationIssue(
                code="empty_rendered_answer",
                answer_section="answer",
                rule="rendered_answer_must_be_nonempty",
                explanation="The server renderer produced an empty answer.",
            )
        )
    if _MALFORMED_CITATION_RE.search(answer):
        issues.append(
            GroundingValidationIssue(
                code="malformed_citation_syntax",
                answer_section="answer",
                rule="server_citations_must_use_double_bracket_ordinals",
                explanation="The rendered answer contained malformed citation syntax.",
            )
        )
    registry_ordinals = {item.ordinal for item in citation_registry}
    for ordinal in (int(match.group(1)) for match in _CITATION_RE.finditer(answer)):
        if ordinal < 1 or ordinal not in registry_ordinals:
            issues.append(
                GroundingValidationIssue(
                    code="unresolved_citation_ordinal",
                    answer_section="answer",
                    rule="rendered_citations_must_resolve_to_registry",
                    explanation="A rendered citation did not resolve to the citation registry.",
                )
            )
    if _contains_actual_record_id(answer, actual_record_ids):
        issues.append(
            GroundingValidationIssue(
                code="private_evidence_record_id_exposed",
                answer_section="answer",
                rule="rendered_answer_must_not_expose_private_record_ids",
                explanation="The rendered answer exposed an actual internal EvidenceRecord identifier.",
            )
        )
    if rendered_slot_ids != expected_slot_ids:
        issues.append(
            GroundingValidationIssue(
                code="unstable_section_order",
                answer_section="answer",
                rule="server_slot_order_must_remain_stable",
                explanation="Rendered prose slots were missing or out of deterministic order.",
            )
        )
    return issues


def _render_with_fragments(
    *,
    slots: list[GroundedProseSlot],
    fragments: list[GroundedProseFragment] | None,
    status: AnswerStatus,
    plan: ScientificQuestionPlan,
    records: list[EvidenceRecord],
    assessments: list[EvidenceRequirementAssessment],
    citation_registry: list[CitationReference],
    recommendations: list[ExperimentRecommendation],
) -> GroundedSynthesisResult:
    slot_by_id = {slot.slot_id: slot for slot in slots}
    supplied_fragments = fragments or []
    counts = Counter(
        fragment.slot_id for fragment in supplied_fragments if fragment.slot_id in slot_by_id
    )
    fragments_by_id = {
        fragment.slot_id: fragment
        for fragment in supplied_fragments
        if fragment.slot_id in slot_by_id and counts[fragment.slot_id] == 1
    }
    actual_record_ids = {record.id for record in records if record.id}
    issues: list[GroundingValidationIssue] = []
    accepted: dict[str, str] = {}
    rendered_slot_ids: list[str] = []
    rendered_lines: list[str] = []
    updated_recommendations = [item.model_copy(deep=True) for item in recommendations]

    for slot in sorted(slots, key=lambda item: item.stable_order):
        fragment = fragments_by_id.get(slot.slot_id)
        if counts[slot.slot_id] > 1:
            issues.append(
                _fragment_issue(
                    slot,
                    "duplicate_slot_id",
                    "each_expected_slot_may_appear_once",
                    "The model returned the expected slot more than once.",
                )
            )
        elif fragment is None:
            if fragments is not None:
                issues.append(
                    _fragment_issue(
                        slot,
                        "missing_slot",
                        "every_requested_slot_requires_a_fragment",
                        "The model omitted a requested prose slot.",
                    )
                )
        else:
            fragment_issues = validate_prose_fragment(
                fragment,
                slot=slot,
                all_plan_genes=tuple(plan.entities.genes),
                actual_record_ids=actual_record_ids,
            )
            issues.extend(fragment_issues)
            if not fragment_issues:
                accepted[slot.slot_id] = fragment.text.strip()

        text = accepted.get(slot.slot_id, slot.fallback_text)
        rendered_slot_ids.append(slot.slot_id)
        if slot.section is ProseSection.recommendation and slot.recommendation_index is not None:
            recommendation = updated_recommendations[slot.recommendation_index]
            if recommendation.rationale_citation_ordinals:
                description = text
            else:
                description = (
                    text
                    if text.casefold().startswith("gap-driven recommendation:")
                    else f"Gap-driven recommendation: {text}"
                )
            updated_recommendations[slot.recommendation_index] = recommendation.model_copy(
                update={"description": description}
            )
            continue
        rendered_lines.append(_render_slot(slot, text))

    opening = (
        "The selected provenance-backed evidence is insufficient to support every required evidence dimension."
        if status is AnswerStatus.insufficient_evidence
        else "The selected provenance-backed evidence supports the following grounded synthesis."
    )
    gap_lines = [
        f"Evidence gap for {item.gene_symbol} {item.evidence_need.value.replace('_', ' ')}: {item.detail}"
        for item in sorted(
            assessments,
            key=lambda item: (item.gene_symbol, item.requirement_id, item.evidence_need.value),
        )
        if item.status is not RequirementStatus.sufficient
    ]
    answer = " ".join([opening, *rendered_lines, *gap_lines]).strip()
    expected_summary_slot_ids = [
        slot.slot_id
        for slot in sorted(slots, key=lambda item: item.stable_order)
        if slot.section is not ProseSection.recommendation
    ]
    rendered_summary_slot_ids = [
        slot_id
        for slot_id in rendered_slot_ids
        if slot_by_id[slot_id].section is not ProseSection.recommendation
    ]
    structural_issues = validate_rendered_answer(
        answer,
        citation_registry=citation_registry,
        actual_record_ids=actual_record_ids,
        rendered_slot_ids=rendered_summary_slot_ids,
        expected_slot_ids=expected_summary_slot_ids,
    )
    issues.extend(structural_issues)
    accepted_count = len(accepted)
    fallback_count = len(slots) - accepted_count
    method = (
        "grounded_llm" if fallback_count == 0 else "hybrid" if accepted_count else "deterministic"
    )
    return GroundedSynthesisResult(
        summary=answer,
        generation_method=method,
        recommendations=updated_recommendations,
        requested_slot_count=len(slots),
        accepted_slot_count=accepted_count,
        fallback_slot_count=fallback_count,
        diagnostic_counts=dict(sorted(Counter(issue.code for issue in issues).items())),
        validation_issues=issues,
    )


def citations_are_valid(text: str, evidence_ids: list[str]) -> bool:
    """Compatibility helper for ordinal syntax/range and exact private-ID checks."""
    if _MALFORMED_CITATION_RE.search(text) or _contains_actual_record_id(text, evidence_ids):
        return False
    return all(
        1 <= int(match.group(1)) <= len(evidence_ids) for match in _CITATION_RE.finditer(text)
    )


def deterministic_summary(
    *,
    status: AnswerStatus,
    plan: ScientificQuestionPlan,
    records: list[EvidenceRecord],
    assessments: list[EvidenceRequirementAssessment],
) -> str:
    """Compatibility deterministic summary using only selected record text."""
    citation_registry = [
        CitationReference(
            ordinal=index,
            evidence_record_id=record.id,
            source_id=record.source_id,
            source_name=record.source_name,
        )
        for index, record in enumerate(records, start=1)
    ]
    slots = build_grounded_prose_slots(plan=plan, records=records, assessments=assessments)
    result = _render_with_fragments(
        slots=slots,
        fragments=None,
        status=status,
        plan=plan,
        records=records,
        assessments=assessments,
        citation_registry=citation_registry,
        recommendations=[],
    )
    return result.summary


def try_grounded_synthesis(
    *,
    question: str,
    status: AnswerStatus,
    plan: ScientificQuestionPlan,
    records: list[EvidenceRecord],
    assessments: list[EvidenceRequirementAssessment],
    gaps: list[EvidenceGap] | None = None,
    recommendations: list[ExperimentRecommendation] | None = None,
    comparison_matrix: list[ComparisonRow] | None = None,
    citation_registry: list[CitationReference] | None = None,
    settings: Settings | None = None,
) -> GroundedSynthesisResult:
    """Fill deterministic prose slots with optional provider-generated fragments."""
    del gaps  # Gap identity remains server-owned through recommendation objects.
    recs = [item.model_copy(deep=True) for item in recommendations or []]
    registry = citation_registry or [
        CitationReference(
            ordinal=index,
            evidence_record_id=record.id,
            source_id=record.source_id,
            source_name=record.source_name,
        )
        for index, record in enumerate(records, start=1)
    ]
    slots = build_grounded_prose_slots(
        plan=plan,
        records=records,
        assessments=assessments,
        recommendations=recs,
        comparison_matrix=comparison_matrix,
    )
    deterministic = _render_with_fragments(
        slots=slots,
        fragments=None,
        status=status,
        plan=plan,
        records=records,
        assessments=assessments,
        citation_registry=registry,
        recommendations=recs,
    )
    cfg = settings or get_settings()
    if not cfg.has_llm() or not records or not slots:
        return deterministic
    try:
        from gene_dossier.synthesis import build_chat_model_candidates
    except Exception as exc:  # noqa: BLE001
        logger.warning("grounded scientific synthesis unavailable: %s", type(exc).__name__)
        return deterministic.model_copy(
            update={"failure_type": "provider_failure", "failure_message": type(exc).__name__}
        )

    prompt = build_grounded_prose_prompt(
        question=question,
        status=status,
        slots=slots,
        private_values=_exact_private_values(records),
    )
    last_failure: tuple[str, str, list[GroundingValidationIssue]] | None = None
    for candidate in build_chat_model_candidates(cfg, purpose="answer"):
        try:
            if candidate.provider == "openai":
                structured = candidate.model.with_structured_output(
                    GroundedProseDraft,
                    method="json_schema",
                    strict=True,
                )
            else:
                structured = candidate.model.with_structured_output(GroundedProseDraft)
            raw = structured.invoke(prompt)
            if isinstance(raw, dict) and raw.get("refusal"):
                last_failure = (
                    "refusal",
                    "Provider refused grounded-answer synthesis.",
                    [
                        GroundingValidationIssue(
                            code="refusal",
                            answer_section="provider",
                            rule="provider_must_return_structured_fragments",
                            explanation="The provider returned a refusal instead of prose fragments.",
                        )
                    ],
                )
                continue
            draft = (
                raw
                if isinstance(raw, GroundedProseDraft)
                else GroundedProseDraft.model_validate(raw)
            )
            result = _render_with_fragments(
                slots=slots,
                fragments=draft.fragments,
                status=status,
                plan=plan,
                records=records,
                assessments=assessments,
                citation_registry=registry,
                recommendations=recs,
            )
            return result
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "grounded synthesis failed via %s: %s", candidate.provider, type(exc).__name__
            )
            message = str(exc).lower()
            failure_type = "provider_failure"
            issue_code = "provider_failure"
            if "refusal" in message:
                failure_type = "refusal"
                issue_code = "refusal"
            elif "length" in message or "incomplete" in message or "max" in message:
                failure_type = "incomplete_response"
                issue_code = "incomplete_output_cap_response"
            elif "validation" in message or "schema" in message:
                failure_type = "schema_failure"
                issue_code = "malformed_structured_answer"
            last_failure = (
                failure_type,
                type(exc).__name__,
                [
                    GroundingValidationIssue(
                        code=issue_code,
                        answer_section="provider",
                        rule="provider_output_must_validate_against_prose_fragment_schema",
                        explanation=f"Provider output failed before slot validation: {type(exc).__name__}.",
                    )
                ],
            )
    if last_failure is None:
        return deterministic
    combined_issues = [*deterministic.validation_issues, *last_failure[2]]
    return deterministic.model_copy(
        update={
            "failure_type": last_failure[0],
            "failure_message": last_failure[1],
            "validation_issues": combined_issues,
            "diagnostic_counts": dict(
                sorted(Counter(issue.code for issue in combined_issues).items())
            ),
        }
    )


def terra_usage_cost_usd(
    *,
    input_tokens: int,
    output_tokens: int,
    cached_input_tokens: int = 0,
) -> Decimal:
    """Return usage-based Terra standard-processing cost in USD."""
    uncached_input = max(input_tokens - cached_input_tokens, 0)
    cost = (
        Decimal(uncached_input) * Decimal("2.00")
        + Decimal(max(cached_input_tokens, 0)) * Decimal("0.20")
        + Decimal(max(output_tokens, 0)) * Decimal("12.00")
    ) / Decimal(1_000_000)
    return cost.quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP)


__all__ = [
    "ClaimLanguagePolicy",
    "EpistemicRole",
    "GroundedProseDraft",
    "GroundedProseFragment",
    "GroundedProseSlot",
    "GroundedSynthesisResult",
    "GroundingValidationIssue",
    "ProseSection",
    "build_grounded_prose_prompt",
    "build_grounded_prose_slots",
    "citations_are_valid",
    "deterministic_summary",
    "terra_usage_cost_usd",
    "try_grounded_synthesis",
    "validate_prose_fragment",
    "validate_rendered_answer",
]
