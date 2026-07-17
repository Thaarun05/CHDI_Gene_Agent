"""Core domain models and enums for the Gene Dossier Platform.

These are pure Pydantic v2 models describing the provenance chain:

    ApiRun (a call) -> RawArtifact (stored response) -> EvidenceRecord (normalized fact)
    -> ReportSection / Claim -> VerificationResult

Persistence (SQLModel tables) is defined separately in ``db.py`` so the domain models
stay decoupled from storage.
"""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum
from typing import Any, Literal
from uuid import uuid4

from pydantic import BaseModel, Field


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def new_id() -> str:
    """Return a new random hex identifier."""
    return uuid4().hex


def utcnow() -> datetime:
    """Return the current timezone-aware UTC timestamp."""
    return datetime.now(timezone.utc)


# Three-way and four-way verdict aliases used by verification.
Verdict3 = Literal["pass", "warning", "fail"]
Verdict4 = Literal["pass", "warning", "fail", "human_review"]


# --------------------------------------------------------------------------------------
# Enums
# --------------------------------------------------------------------------------------
class EvidenceGrade(str, Enum):
    """Strength of an evidence record.

    A = direct human genetic or curated causal evidence
    B = human expression, eQTL, or disease association evidence
    C = curated protein, pathway, or PPI evidence
    D = mouse or cell model evidence
    E = predicted, computational, or text-mining evidence
    F = weak mention, broad list, or needs review
    """

    A = "A"
    B = "B"
    C = "C"
    D = "D"
    E = "E"
    F = "F"


class SourceType(str, Enum):
    """Category of a data source."""

    curated_database = "curated_database"
    literature = "literature"
    expression_database = "expression_database"
    genetic_database = "genetic_database"
    pathway_database = "pathway_database"
    interaction_database = "interaction_database"
    chemical_database = "chemical_database"
    structure_database = "structure_database"
    model_organism_database = "model_organism_database"
    patent_database = "patent_database"
    grant_database = "grant_database"
    commercial_source = "commercial_source"
    manual_visual_reference = "manual_visual_reference"
    hd_specific_database = "hd_specific_database"
    semi_structured_source = "semi_structured_source"


class AssertionType(str, Enum):
    """The kind of factual assertion an evidence record makes."""

    gene_identity = "gene_identity"
    orthology = "orthology"
    expression = "expression"
    cell_type_expression = "cell_type_expression"
    protein_function = "protein_function"
    protein_structure = "protein_structure"
    ppi = "ppi"
    pathway_membership = "pathway_membership"
    variant_association = "variant_association"
    disease_association = "disease_association"
    chemical_interaction = "chemical_interaction"
    perturbation = "perturbation"
    transcription_factor_association = "transcription_factor_association"
    knockout_phenotype = "knockout_phenotype"
    patent_claim = "patent_claim"
    grant_project = "grant_project"
    visual_observation = "visual_observation"
    literature_summary = "literature_summary"


class SourceStatus(str, Enum):
    """Per-source outcome in a dossier run. A source is never silently omitted."""

    success = "success"
    failed = "failed"
    deferred = "deferred"
    manual = "manual"
    requires_key = "requires_key"
    partial = "partial"
    skipped = "skipped"
    not_implemented = "not_implemented"


# --------------------------------------------------------------------------------------
# Provenance / run models
# --------------------------------------------------------------------------------------
class DossierRun(BaseModel):
    """A single dossier generation run for one gene."""

    id: str = Field(default_factory=new_id)
    gene_symbol: str
    official_symbol: str | None = None
    run_type: str = "full_dossier"
    focus: str | None = None
    status: str = "created"
    started_at: datetime = Field(default_factory=utcnow)
    completed_at: datetime | None = None
    config: dict[str, Any] | None = None
    notes: str | None = None


class ApiRun(BaseModel):
    """Metadata for a single API call made during a run."""

    id: str = Field(default_factory=new_id)
    dossier_run_id: str
    gene_symbol: str
    source_name: str
    endpoint_name: str
    method: str = "GET"
    request_url: str
    request_params: dict[str, Any] = Field(default_factory=dict)
    status_code: int | None = None
    success: bool = False
    error_type: str | None = None
    error_message: str | None = None
    retrieved_at: datetime = Field(default_factory=utcnow)
    raw_artifact_id: str | None = None


class RawArtifact(BaseModel):
    """A raw source response or artifact stored on disk with a content hash."""

    id: str = Field(default_factory=new_id)
    dossier_run_id: str
    api_run_id: str | None = None
    source_name: str
    artifact_type: str
    file_path: str
    original_url: str | None = None
    content_hash: str
    captured_at: datetime = Field(default_factory=utcnow)
    notes: str | None = None


class EvidenceRecord(BaseModel):
    """A normalized, source-level factual unit extracted from a raw response.

    This is the unit of truth for reports: every claim must cite one or more
    ``source_id`` values that resolve to evidence records.
    """

    id: str = Field(default_factory=new_id)
    source_id: str
    dossier_run_id: str
    gene_symbol: str
    official_symbol: str | None = None
    section: str
    subsection: str | None = None
    source_name: str
    source_type: SourceType
    assertion_type: AssertionType
    fact_type: str
    organism: str | None = None
    species: str | None = None
    taxon_id: int | None = None
    evidence_grade: EvidenceGrade
    manual_review_required: bool = False
    confidence_notes: str | None = None
    value: dict[str, Any] = Field(default_factory=dict)
    display_text: str
    api_run_id: str | None = None
    raw_artifact_id: str | None = None
    raw_response_pointer: str | None = None
    created_at: datetime = Field(default_factory=utcnow)


# --------------------------------------------------------------------------------------
# Report / claim / verification models
# --------------------------------------------------------------------------------------
class ReportSection(BaseModel):
    """A rendered section of the dossier, with cited sources/figures/tables."""

    id: str = Field(default_factory=new_id)
    dossier_run_id: str
    section_name: str
    subsection_name: str | None = None
    content_markdown: str = ""
    source_ids: list[str] = Field(default_factory=list)
    figure_ids: list[str] = Field(default_factory=list)
    table_ids: list[str] = Field(default_factory=list)
    status: str = "draft"
    created_at: datetime = Field(default_factory=utcnow)
    updated_at: datetime = Field(default_factory=utcnow)


class Claim(BaseModel):
    """A discrete factual claim extracted from a report section; must cite sources."""

    id: str = Field(default_factory=new_id)
    dossier_run_id: str
    section_id: str | None = None
    claim_text: str
    source_ids: list[str] = Field(default_factory=list)
    evidence_grade: EvidenceGrade | None = None
    claim_type: str | None = None
    created_at: datetime = Field(default_factory=utcnow)


class VerificationResult(BaseModel):
    """The outcome of verifying a single claim against evidence records."""

    id: str = Field(default_factory=new_id)
    claim_id: str
    source_id_presence_passed: bool
    source_exists_passed: bool
    semantic_support: Verdict3 = "pass"
    causal_language_check: Verdict3 = "pass"
    evidence_strength_check: Verdict3 = "pass"
    verdict: Verdict4 = "pass"
    reason: str | None = None
    needs_human_review: bool = False
    created_at: datetime = Field(default_factory=utcnow)


# --------------------------------------------------------------------------------------
# Transport / reporting models (not persisted as provenance)
# --------------------------------------------------------------------------------------
class ToolResult(BaseModel):
    """Uniform return type for every API client in ``tools/``.

    Clients never raise: transport/HTTP errors are captured here so the workflow can
    record coverage and continue with other sources.
    """

    source_name: str
    endpoint_name: str
    success: bool
    gene_symbol: str
    request_url: str
    request_params: dict[str, Any] = Field(default_factory=dict)
    status_code: int | None = None
    data: Any | None = None
    raw_artifact_id: str | None = None
    error_type: str | None = None
    error_message: str | None = None


class SourceCoverageResult(BaseModel):
    """Per-source status line for the source coverage report."""

    dossier_run_id: str
    source_name: str
    status: SourceStatus
    raw_artifact_path: str | None = None
    evidence_record_count: int = 0
    error_message: str | None = None
    report_sections_supported: list[str] = Field(default_factory=list)
    notes: str | None = None


__all__ = [
    "new_id",
    "utcnow",
    "Verdict3",
    "Verdict4",
    "EvidenceGrade",
    "SourceType",
    "AssertionType",
    "SourceStatus",
    "DossierRun",
    "ApiRun",
    "RawArtifact",
    "EvidenceRecord",
    "ReportSection",
    "Claim",
    "VerificationResult",
    "ToolResult",
    "SourceCoverageResult",
]
