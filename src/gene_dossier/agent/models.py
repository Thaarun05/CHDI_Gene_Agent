"""Strict domain models for scientific question planning and assessment."""

from __future__ import annotations

from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

from gene_dossier.models import EvidenceRecord


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ScientificIntent(str, Enum):
    single_gene_question = "single_gene_question"
    multi_gene_question = "multi_gene_question"
    comparison = "comparison"
    mechanistic_question = "mechanistic_question"
    evidence_gap_question = "evidence_gap_question"
    general_biomedical_question = "general_biomedical_question"
    out_of_scope = "out_of_scope"


class AnswerMode(str, Enum):
    fact = "fact"
    synthesis = "synthesis"
    comparison = "comparison"
    mechanistic = "mechanistic"
    gap_analysis = "gap_analysis"


class PlannerMethod(str, Enum):
    llm_structured = "llm_structured"
    deterministic_fallback = "deterministic_fallback"


class AnswerStatus(str, Enum):
    answered = "answered"
    insufficient_evidence = "insufficient_evidence"
    clarification_required = "clarification_required"
    out_of_scope = "out_of_scope"


class RequirementStatus(str, Enum):
    sufficient = "sufficient"
    limited = "limited"
    missing = "missing"
    unsupported_capability = "unsupported_capability"


class EvidenceSelection(str, Enum):
    accepted_only = "accepted_only"
    accepted_or_latest_generated = "accepted_or_latest_generated"
    explicit_only = "explicit_only"


class CapabilityId(str, Enum):
    identity_function = "identity_function"
    orthology_conservation = "orthology_conservation"
    structure_domain = "structure_domain"
    expression_context = "expression_context"
    brain_expression = "brain_expression"
    experimental_expression = "experimental_expression"
    transcriptional_regulation = "transcriptional_regulation"
    ppi = "ppi"
    pathway = "pathway"
    hd_literature = "hd_literature"
    disease_association = "disease_association"
    human_genetic_association = "human_genetic_association"
    model_organism = "model_organism"
    chemical_perturbation = "chemical_perturbation"
    chemical_tools = "chemical_tools"


class EvidenceNeed(str, Enum):
    identity_function = "identity_function"
    orthology_conservation = "orthology_conservation"
    structure_domain = "structure_domain"
    expression_context = "expression_context"
    brain_expression = "brain_expression"
    experimental_evidence = "experimental_evidence"
    transcriptional_regulation = "transcriptional_regulation"
    protein_interaction = "protein_interaction"
    pathway_membership = "pathway_membership"
    hd_literature = "hd_literature"
    disease_association = "disease_association"
    human_genetic_association = "human_genetic_association"
    repeat_instability_mechanism = "repeat_instability_mechanism"
    model_organism = "model_organism"
    chemical_perturbation = "chemical_perturbation"
    therapeutic_perturbability = "therapeutic_perturbability"


class ScientificEntities(StrictModel):
    genes: list[str] = Field(default_factory=list, max_length=6)
    diseases: list[str] = Field(default_factory=list)
    biological_processes: list[str] = Field(default_factory=list)
    pathways: list[str] = Field(default_factory=list)
    chemicals: list[str] = Field(default_factory=list)

    @field_validator("genes")
    @classmethod
    def normalize_genes(cls, genes: list[str]) -> list[str]:
        return list(dict.fromkeys(gene.strip().upper() for gene in genes if gene.strip()))


class EvidenceRequirement(StrictModel):
    id: str = Field(min_length=1, max_length=80, pattern=r"^[a-z0-9][a-z0-9_-]*$")
    label: str = Field(min_length=1, max_length=160)
    description: str = Field(min_length=1, max_length=600)
    genes: list[str] = Field(min_length=1, max_length=6)
    evidence_need: EvidenceNeed
    capability_ids: list[CapabilityId] = Field(default_factory=list, max_length=8)
    required: bool
    minimum_support: int = Field(ge=1, le=10)
    rationale: str = Field(min_length=1, max_length=600)

    @field_validator("genes")
    @classmethod
    def normalize_genes(cls, genes: list[str]) -> list[str]:
        return list(dict.fromkeys(gene.strip().upper() for gene in genes if gene.strip()))


class ScientificQuestionPlanDraft(StrictModel):
    intent: ScientificIntent
    entities: ScientificEntities
    primary_gene: str | None = None
    objective: str = Field(min_length=1, max_length=500)
    analysis_lens: str = Field(default="general", max_length=80)
    answer_mode: AnswerMode
    evidence_requirements: list[EvidenceRequirement] = Field(max_length=10)
    requires_multi_gene: bool


class ScientificQuestionPlan(ScientificQuestionPlanDraft):
    planner_method: PlannerMethod


class EvidenceRequirementAssessment(StrictModel):
    requirement_id: str
    gene_symbol: str
    evidence_need: EvidenceNeed
    required: bool
    minimum_support: int
    status: RequirementStatus
    qualifying_count: int
    evidence_record_ids: list[str] = Field(default_factory=list)
    contributing_capability_ids: list[CapabilityId] = Field(default_factory=list)
    detail: str


class AgentEvidenceUniverse(StrictModel):
    gene_symbol: str
    base_evidence_run_id: str | None = None
    explicit_run_ids: list[str] = Field(default_factory=list)
    reused_tool_run_ids: list[str] = Field(default_factory=list)
    created_tool_run_ids: list[str] = Field(default_factory=list)
    tool_run_ids: list[str] = Field(default_factory=list)
    dossier_run_ids: list[str] = Field(default_factory=list)
    evidence_universe: str


class ToolActivity(StrictModel):
    gene_symbol: str
    capability_ids: list[CapabilityId]
    executor_kind: str
    status: str
    dossier_run_id: str | None = None
    section_keys: list[str] = Field(default_factory=list)
    sources: list[str] = Field(default_factory=list)
    evidence_records_persisted: int = 0
    indexed_records: int = 0
    reused: bool = False
    errors: list[str] = Field(default_factory=list)


class ScientificFailure(StrictModel):
    failure_type: str
    message: str
    provider: str | None = None


class CitationReference(StrictModel):
    ordinal: int
    evidence_record_id: str
    source_id: str
    source_name: str


class EvidenceCategoryBlock(StrictModel):
    gene_symbol: str
    category: str
    evidence_need: EvidenceNeed
    evidence_system: str
    claim_type: str
    evidence_record_ids: list[str] = Field(default_factory=list)
    citation_ordinals: list[int] = Field(default_factory=list)
    source_names: list[str] = Field(default_factory=list)
    retrieval_timestamps: list[str] = Field(default_factory=list)
    status: str
    summary: str


class EvidenceGap(StrictModel):
    id: str
    gene_symbol: str
    requirement_id: str
    evidence_need: EvidenceNeed
    status: RequirementStatus
    required: bool
    detail: str


class ExperimentRecommendation(StrictModel):
    label: Literal["Recommendation"] = "Recommendation"
    description: str
    gap_ids: list[str]
    decision_uncertainty: str
    rationale_citation_ordinals: list[int] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class SourceAttempt(StrictModel):
    gene_symbol: str
    dossier_run_id: str
    source_name: str
    status: str
    retrieved_at: str | None = None
    error_message: str | None = None


ComparisonStrength = Literal["Strong", "Moderate", "Limited", "Weak", "Missing"]


class ComparisonCell(StrictModel):
    status: ComparisonStrength
    summary: str
    evidence_count: int
    evidence_record_ids: list[str] = Field(default_factory=list)


class ComparisonRow(StrictModel):
    dimension: str
    cells: dict[str, ComparisonCell]


class ScientificAgentResult(BaseModel):
    status: AnswerStatus
    question: str
    context_gene: str | None = None
    plan: ScientificQuestionPlan | None = None
    evidence_universes: dict[str, AgentEvidenceUniverse] = Field(default_factory=dict)
    assessments: list[EvidenceRequirementAssessment] = Field(default_factory=list)
    selected_records: list[EvidenceRecord] = Field(default_factory=list, exclude=True)
    summary: str
    retrieval_method: str = "abstain"
    generation_method: str = "abstain"
    embedding_backend: str = "unavailable"
    limitations: list[str] = Field(default_factory=list)
    evidence_gaps: list[str] = Field(default_factory=list)
    tool_activity: list[ToolActivity] = Field(default_factory=list)
    agent_activity: list[str] = Field(default_factory=list)
    comparison_dimensions: list[str] = Field(default_factory=list)
    comparison_matrix: list[ComparisonRow] = Field(default_factory=list)
    evidence_categories: list[EvidenceCategoryBlock] = Field(default_factory=list)
    structured_gaps: list[EvidenceGap] = Field(default_factory=list)
    recommendations: list[ExperimentRecommendation] = Field(default_factory=list)
    citation_registry: list[CitationReference] = Field(default_factory=list)
    source_attempts: list[SourceAttempt] = Field(default_factory=list)
    retrieval_timestamps: list[str] = Field(default_factory=list)
    failures: list[ScientificFailure] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)
