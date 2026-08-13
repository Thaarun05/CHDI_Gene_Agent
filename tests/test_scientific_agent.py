"""Offline regression coverage for the bounded general scientific agent."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from gene_dossier.agent.capabilities import (
    record_matches_capability,
    record_matches_need,
    validate_source_capabilities,
)
from gene_dossier.agent.audit import audit_evidence_overlap, scientific_fingerprint
from gene_dossier.agent.comparison import grade_hd_modifier_cell
from gene_dossier.agent.models import (
    AnswerMode,
    AnswerStatus,
    CapabilityId,
    EvidenceNeed,
    EvidenceRequirement,
    EvidenceRequirementAssessment,
    EvidenceSelection,
    PlannerMethod,
    RequirementStatus,
    ScientificEntities,
    ScientificIntent,
    ScientificQuestionPlan,
    ScientificQuestionPlanDraft,
)
from gene_dossier.agent.orchestrator import (
    ScientificAgentRequest,
    ScientificAgentService,
    summarize_agent_result_runs,
)
from gene_dossier.agent.planner import (
    PlanResult,
    _try_structured_llm_plan,
    _validate_plan,
    build_gemini_planner_schema,
    deterministic_fallback_plan,
    plan_scientific_question,
)
from gene_dossier.agent.synthesis import (
    ClaimLanguagePolicy,
    EpistemicRole,
    GroundedProseDraft,
    GroundedProseFragment,
    GroundedProseSlot,
    ProseSection,
    build_grounded_prose_slots,
    citations_are_valid,
    terra_usage_cost_usd,
    try_grounded_synthesis,
    validate_prose_fragment,
    validate_rendered_answer,
)
from gene_dossier.config import Settings
from gene_dossier.db import (
    GeneratedReportRow,
    canonical_generated_report_id,
    init_db,
    list_evidence_for_run,
    save_api_run,
    save_dossier_run,
    save_evidence_record,
    save_generated_report,
    save_raw_artifact,
    session_scope,
)
from gene_dossier.models import (
    ApiRun,
    AssertionType,
    DossierRun,
    EvidenceGrade,
    EvidenceRecord,
    RawArtifact,
    SourceType,
    ToolResult,
    new_id,
)
from gene_dossier.retrieval import (
    ChromaEvidenceIndex,
    HashEmbeddingFunction,
    collection_name_for_embedding,
    EmbeddingIndexIdentity,
)
from gene_dossier.workflow import normalize_tool_result


class _NoSemanticIndex:
    available = False
    status = SimpleNamespace(embedding_backend="unavailable")

    def upsert_evidence(self, _records: list[EvidenceRecord]) -> int:
        return 0


def _settings(tmp_path: Path) -> Settings:
    return Settings(
        database_url=f"sqlite:///{tmp_path / 'agent.db'}",
        raw_data_dir=tmp_path / "raw",
        output_dir=tmp_path / "outputs",
        index_dir=tmp_path / "indexes",
        openai_api_key=None,
        anthropic_api_key=None,
        nvidia_nim_api_key=None,
        google_api_key=None,
    )


def _record(
    *,
    gene: str,
    run_id: str,
    source: str,
    assertion: AssertionType,
    text: str,
    source_type: SourceType = SourceType.curated_database,
    grade: EvidenceGrade = EvidenceGrade.C,
    fact_type: str = "test_fact",
) -> EvidenceRecord:
    return EvidenceRecord(
        id=new_id(),
        source_id=f"{source.lower()}:{gene.lower()}:{new_id()}",
        dossier_run_id=run_id,
        gene_symbol=gene,
        official_symbol=gene,
        section="Scientific evidence",
        source_name=source,
        source_type=source_type,
        assertion_type=assertion,
        fact_type=fact_type,
        evidence_grade=grade,
        value={"text": text},
        display_text=text,
    )


def _persist_run(service: ScientificAgentService, gene: str, run_id: str, records: list[EvidenceRecord]) -> None:
    now = datetime.now(timezone.utc)
    with session_scope(service.engine) as session:
        save_dossier_run(
            session,
            DossierRun(
                id=run_id,
                gene_symbol=gene,
                official_symbol=gene,
                run_type="scientific_agent_test",
                status="completed",
                started_at=now,
                completed_at=now,
                config={},
            ),
        )
        for index, record in enumerate(records, start=1):
            api_id = f"api-{run_id}-{index}"
            raw_id = f"raw-{run_id}-{index}"
            save_api_run(
                session,
                ApiRun(
                    id=api_id,
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    source_name=record.source_name,
                    endpoint_name="offline_fixture",
                    request_url=f"https://example.invalid/{record.source_name}",
                    success=True,
                    status_code=200,
                    raw_artifact_id=raw_id,
                ),
            )
            save_raw_artifact(
                session,
                RawArtifact(
                    id=raw_id,
                    dossier_run_id=run_id,
                    api_run_id=api_id,
                    source_name=record.source_name,
                    artifact_type="json",
                    file_path=f"raw/{run_id}/{raw_id}.json",
                    content_hash=(f"{index:064x}"[-64:]),
                ),
            )
            save_evidence_record(
                session,
                record.model_copy(update={"api_run_id": api_id, "raw_artifact_id": raw_id}),
            )


def _requirement(
    need: EvidenceNeed,
    genes: list[str],
    capability: CapabilityId,
    *,
    required: bool = True,
    minimum: int = 2,
    requirement_id: str | None = None,
) -> EvidenceRequirement:
    return EvidenceRequirement(
        id=requirement_id or need.value,
        label=need.value.replace("_", " ").title(),
        description=f"Qualifying {need.value} evidence.",
        genes=genes,
        evidence_need=need,
        capability_ids=[capability],
        required=required,
        minimum_support=minimum,
        rationale="Offline scientific-agent test requirement.",
    )


def _plan(genes: list[str], requirements: list[EvidenceRequirement], *, lens: str = "general") -> ScientificQuestionPlan:
    return ScientificQuestionPlan(
        intent=ScientificIntent.comparison if len(genes) > 1 else ScientificIntent.single_gene_question,
        entities=ScientificEntities(genes=genes),
        primary_gene=genes[0] if len(genes) == 1 else None,
        objective="hd_modifier_relevance" if lens == "hd_modifier_relevance" else "offline evidence objective",
        analysis_lens=lens,
        answer_mode=AnswerMode.comparison if len(genes) > 1 else AnswerMode.synthesis,
        evidence_requirements=requirements,
        requires_multi_gene=len(genes) > 1,
        planner_method=PlannerMethod.llm_structured,
    )


def _planner(plan: ScientificQuestionPlan):
    def fake_planner(_question: str, **_kwargs: Any) -> PlanResult:
        return PlanResult(plan)

    return fake_planner


def _assessment_for(
    need: EvidenceNeed,
    *,
    status: RequirementStatus = RequirementStatus.sufficient,
    count: int = 3,
) -> EvidenceRequirementAssessment:
    return EvidenceRequirementAssessment(
        requirement_id=f"assessment-{need.value}",
        gene_symbol="GENE1",
        evidence_need=need,
        required=True,
        minimum_support=2,
        status=status,
        qualifying_count=count,
        detail="Offline HD rubric test assessment.",
    )


def _service(tmp_path: Path, *, planner=None, baselines: dict[str, str] | None = None, section_executor=None, source_executor=None) -> ScientificAgentService:
    kwargs: dict[str, Any] = {
        "accepted_baselines": baselines or {},
        "settings": _settings(tmp_path),
        "index_factory": lambda: _NoSemanticIndex(),
    }
    if planner is not None:
        kwargs["planner"] = planner
    if section_executor is not None:
        kwargs["section_executor"] = section_executor
    if source_executor is not None:
        kwargs["source_executor"] = source_executor
    service = ScientificAgentService(**kwargs)
    init_db(service.engine)
    return service


def test_cdh10_ppi_uses_sufficient_baseline_without_acquisition(tmp_path: Path) -> None:
    service = _service(tmp_path, baselines={"CDH10": "cdh10-base"})
    records = [
        _record(gene="CDH10", run_id="cdh10-base", source="STRING", assertion=AssertionType.ppi, text="CDH10 protein interaction partner A.", source_type=SourceType.interaction_database),
        _record(gene="CDH10", run_id="cdh10-base", source="BioGRID", assertion=AssertionType.ppi, text="CDH10 protein interaction partner B.", source_type=SourceType.interaction_database),
    ]
    _persist_run(service, "CDH10", "cdh10-base", records)

    result = service.execute(ScientificAgentRequest(question="What proteins interact with CDH10?", context_gene="SREBF2"))

    assert result.status is AnswerStatus.answered
    assert result.plan is not None and result.plan.entities.genes == ["CDH10"]
    assert result.evidence_universes["CDH10"].dossier_run_ids == ["cdh10-base"]
    assert not result.tool_activity


def test_optional_gap_does_not_make_answer_insufficient(tmp_path: Path) -> None:
    required = _requirement(EvidenceNeed.identity_function, ["GENE1"], CapabilityId.identity_function, minimum=1)
    optional = _requirement(EvidenceNeed.protein_interaction, ["GENE1"], CapabilityId.ppi, required=False, minimum=2)
    service = _service(tmp_path, planner=_planner(_plan(["GENE1"], [required, optional])), baselines={"GENE1": "base"})
    _persist_run(
        service,
        "GENE1",
        "base",
        [_record(gene="GENE1", run_id="base", source="NCBI Gene", assertion=AssertionType.gene_identity, text="GENE1 identity record.")],
    )

    result = service.execute(
        ScientificAgentRequest(question="Summarize GENE1 with optional interaction context.", context_gene="GENE1", allow_tool_acquisition=False)
    )

    assert result.status is AnswerStatus.answered
    assert any("protein_interaction" in gap for gap in result.evidence_gaps)


def test_open_ended_srebf2_question_combines_cross_category_evidence(tmp_path: Path) -> None:
    requirements = [
        _requirement(EvidenceNeed.identity_function, ["SREBF2"], CapabilityId.identity_function, minimum=1),
        _requirement(EvidenceNeed.pathway_membership, ["SREBF2"], CapabilityId.pathway, minimum=1),
        _requirement(EvidenceNeed.hd_literature, ["SREBF2"], CapabilityId.hd_literature, minimum=2),
        _requirement(EvidenceNeed.brain_expression, ["SREBF2"], CapabilityId.brain_expression, required=False, minimum=1),
    ]
    service = _service(
        tmp_path,
        planner=_planner(_plan(["SREBF2"], requirements)),
        baselines={"SREBF2": "srebf2-base"},
    )
    records = [
        _record(gene="SREBF2", run_id="srebf2-base", source="UniProt", assertion=AssertionType.protein_function, text="SREBF2 protein function includes cholesterol biology."),
        _record(gene="SREBF2", run_id="srebf2-base", source="Reactome", assertion=AssertionType.pathway_membership, text="SREBF2 Reactome cholesterol biosynthesis pathway.", source_type=SourceType.pathway_database),
        _record(gene="SREBF2", run_id="srebf2-base", source="PubMed", assertion=AssertionType.literature_summary, text="SREBF2 cholesterol biology in a Huntington disease study.", source_type=SourceType.literature, grade=EvidenceGrade.F),
        _record(gene="SREBF2", run_id="srebf2-base", source="PubMed", assertion=AssertionType.literature_summary, text="SREBF2 pathway context in Huntington disease literature.", source_type=SourceType.literature, grade=EvidenceGrade.F),
    ]
    _persist_run(service, "SREBF2", "srebf2-base", records)

    result = service.execute(
        ScientificAgentRequest(
            question="What evidence connects SREBF2 cholesterol biology to Huntington's disease?",
            context_gene="SREBF2",
            allow_tool_acquisition=False,
        )
    )

    assert result.status is AnswerStatus.answered
    assert {record.source_name for record in result.selected_records} == {"UniProt", "Reactome", "PubMed"}
    assert len(result.assessments) == 4
    assert any(item.evidence_need is EvidenceNeed.brain_expression and not item.required for item in result.assessments)
    assert "[Evidence gap]" in result.summary


def test_disease_association_never_satisfies_human_genetic_modifier_need(tmp_path: Path) -> None:
    requirement = _requirement(
        EvidenceNeed.human_genetic_association,
        ["GENE2"],
        CapabilityId.human_genetic_association,
        minimum=1,
    )
    service = _service(tmp_path, planner=_planner(_plan(["GENE2"], [requirement])), baselines={"GENE2": "base"})
    disease_record = _record(
        gene="GENE2",
        run_id="base",
        source="Open Targets",
        assertion=AssertionType.disease_association,
        text="GENE2 has an Open Targets Huntington disease association score.",
        source_type=SourceType.genetic_database,
        grade=EvidenceGrade.E,
    )
    _persist_run(service, "GENE2", "base", [disease_record])

    result = service.execute(
        ScientificAgentRequest(question="Is GENE2 a human genetic HD modifier?", context_gene="GENE2", allow_tool_acquisition=False)
    )

    assert result.status is AnswerStatus.insufficient_evidence
    assert result.assessments[0].status.value == "unsupported_capability"
    assert not record_matches_need(disease_record, EvidenceNeed.human_genetic_association)


def test_complex_multi_gene_fallback_never_leaks_context_gene() -> None:
    outcome = deterministic_fallback_plan(
        "Compare MSH3, FAN1, PMS2, and HTT for HD modifier relevance.",
        "SREBF2",
        known_genes=["SREBF2"],
    )

    assert outcome.plan is None
    assert outcome.status is AnswerStatus.clarification_required
    assert "SREBF2" not in (outcome.message or "")


def test_planner_failure_preserves_safe_single_gene_fallback(tmp_path: Path) -> None:
    settings = _settings(tmp_path)
    outcome = plan_scientific_question(
        "What proteins interact with CDH10?",
        context_gene="SREBF2",
        settings=settings,
        known_genes=["SREBF2", "CDH10"],
    )

    assert outcome.plan is not None
    assert outcome.plan.entities.genes == ["CDH10"]
    assert outcome.plan.planner_method is PlannerMethod.deterministic_fallback


@pytest.mark.parametrize("symbol", ["FAN1", "fan1", "Fan1"])
def test_unknown_single_gene_fallback_never_substitutes_context_gene(
    symbol: str,
) -> None:
    unresolved = deterministic_fallback_plan(
        f"What proteins interact with {symbol}?",
        "SREBF2",
        known_genes=["SREBF2"],
    )
    contextual = deterministic_fallback_plan(
        "What proteins interact?",
        "CDH10",
        known_genes=["CDH10"],
    )
    known = deterministic_fallback_plan(
        f"What proteins interact with {symbol}?",
        "SREBF2",
        known_genes=["SREBF2", "FAN1"],
    )

    assert unresolved.plan is None
    assert unresolved.status is AnswerStatus.clarification_required
    assert "could not be safely resolved" in (unresolved.message or "")
    assert contextual.plan is not None
    assert contextual.plan.entities.genes == ["CDH10"]
    assert known.plan is not None
    assert known.plan.entities.genes == ["FAN1"]


@pytest.mark.parametrize("symbol", ["FAN1", "fan1", "Fan1"])
def test_structured_plan_rejects_context_gene_when_question_names_other_target(
    symbol: str,
) -> None:
    requirement = _requirement(
        EvidenceNeed.protein_interaction,
        ["SREBF2"],
        CapabilityId.ppi,
    )
    draft = ScientificQuestionPlanDraft(
        intent=ScientificIntent.single_gene_question,
        entities=ScientificEntities(genes=["SREBF2"]),
        primary_gene="SREBF2",
        objective="Identify interaction evidence.",
        analysis_lens="general",
        answer_mode=AnswerMode.fact,
        evidence_requirements=[requirement],
        requires_multi_gene=False,
    )

    with pytest.raises(ValueError, match="conflicts with the structured plan"):
        _validate_plan(
            draft,
            question=f"What proteins interact with {symbol}?",
            context_gene="SREBF2",
        )


def test_structured_plan_replaces_invalid_capability_with_server_mapping() -> None:
    requirement = _requirement(
        EvidenceNeed.human_genetic_association,
        ["SREBF2"],
        CapabilityId.hd_literature,
        minimum=1,
    )
    draft = ScientificQuestionPlanDraft(
        intent=ScientificIntent.mechanistic_question,
        entities=ScientificEntities(genes=["SREBF2"]),
        primary_gene="SREBF2",
        objective="Connect SREBF2 biology to Huntington disease.",
        analysis_lens="general",
        answer_mode=AnswerMode.synthesis,
        evidence_requirements=[requirement],
        requires_multi_gene=False,
    )

    plan = _validate_plan(
        draft,
        question="What connects SREBF2 biology to Huntington disease?",
        context_gene="SREBF2",
    )

    assert plan.evidence_requirements[0].capability_ids == [
        CapabilityId.human_genetic_association
    ]


def test_structured_plan_canonicalizes_expression_requirement_metadata() -> None:
    requirement = _requirement(
        EvidenceNeed.expression_context,
        ["SREBF2"],
        CapabilityId.brain_expression,
        minimum=2,
    ).model_copy(
        update={
            "label": "Human Genetic Association",
            "description": "Planner-supplied category name that conflicts with evidence_need.",
            "rationale": "Planner rationale should remain intact.",
        }
    )
    draft = ScientificQuestionPlanDraft(
        intent=ScientificIntent.mechanistic_question,
        entities=ScientificEntities(genes=["SREBF2"]),
        primary_gene="SREBF2",
        objective="Connect SREBF2 biology to Huntington disease.",
        analysis_lens="general",
        answer_mode=AnswerMode.synthesis,
        evidence_requirements=[requirement],
        requires_multi_gene=False,
    )

    plan = _validate_plan(
        draft,
        question="What connects SREBF2 biology to Huntington disease?",
        context_gene="SREBF2",
    )
    validated = plan.evidence_requirements[0]

    assert validated.evidence_need is EvidenceNeed.expression_context
    assert validated.label == "Expression Context"
    assert validated.description == "Tissue, cell-type, or biological-context expression evidence."
    assert validated.required == requirement.required
    assert validated.minimum_support == 2
    assert validated.genes == ["SREBF2"]
    assert validated.rationale == "Planner rationale should remain intact."


def test_structured_plan_canonicalizes_human_genetic_requirement_metadata() -> None:
    requirement = _requirement(
        EvidenceNeed.human_genetic_association,
        ["SREBF2"],
        CapabilityId.human_genetic_association,
        minimum=1,
    ).model_copy(
        update={
            "label": "Expression Context",
            "description": "Planner-supplied category name that conflicts with evidence_need.",
        }
    )
    draft = ScientificQuestionPlanDraft(
        intent=ScientificIntent.mechanistic_question,
        entities=ScientificEntities(genes=["SREBF2"]),
        primary_gene="SREBF2",
        objective="Connect SREBF2 biology to Huntington disease.",
        analysis_lens="general",
        answer_mode=AnswerMode.synthesis,
        evidence_requirements=[requirement],
        requires_multi_gene=False,
    )

    plan = _validate_plan(
        draft,
        question="What connects SREBF2 biology to Huntington disease?",
        context_gene="SREBF2",
    )
    validated = plan.evidence_requirements[0]

    assert validated.evidence_need is EvidenceNeed.human_genetic_association
    assert validated.label == "Human Genetic Association"
    assert validated.description == "Direct human genetic association or modifier evidence."


def test_gemini_native_json_schema_output_passes_server_validation(
    monkeypatch,
) -> None:
    captured: dict[str, Any] = {}
    question = "Compare MSH3, FAN1, PMS2, and HTT for HD modifier relevance."
    draft = ScientificQuestionPlanDraft(
        intent=ScientificIntent.comparison,
        entities=ScientificEntities(genes=["MSH3", "FAN1", "PMS2", "HTT"]),
        primary_gene=None,
        objective="Compare direct human genetic modifier evidence.",
        analysis_lens="hd_modifier_relevance",
        answer_mode=AnswerMode.comparison,
        evidence_requirements=[
            _requirement(
                EvidenceNeed.human_genetic_association,
                ["MSH3", "FAN1", "PMS2", "HTT"],
                CapabilityId.hd_literature,
                minimum=1,
            ).model_copy(
                update={
                    "label": "Planner-controlled label",
                    "description": "Planner-controlled description.",
                }
            )
        ],
        requires_multi_gene=True,
    )

    class FakeStructured:
        def invoke(self, prompt: str) -> dict[str, Any]:
            captured["prompt"] = prompt
            return draft.model_dump(mode="json")

    class FakeGemini:
        def with_structured_output(self, **kwargs: Any) -> FakeStructured:
            captured.update(kwargs)
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **_kwargs: [
            SimpleNamespace(provider="google_gemini", model=FakeGemini())
        ],
    )
    settings = _settings(Path("/tmp")).model_copy(
        update={"google_api_key": "test-key", "default_llm_provider": "google_gemini"}
    )

    plan = _try_structured_llm_plan(question, "SREBF2", settings)

    assert captured["schema"] == build_gemini_planner_schema()
    assert captured["method"] == "json_schema"
    assert plan is not None
    assert plan.entities.genes == ["MSH3", "FAN1", "PMS2", "HTT"]
    assert plan.planner_method is PlannerMethod.llm_structured
    requirement = plan.evidence_requirements[0]
    assert requirement.genes == ["MSH3", "FAN1", "PMS2", "HTT"]
    assert requirement.required is True
    assert requirement.minimum_support == 1
    assert requirement.evidence_need is EvidenceNeed.human_genetic_association
    assert requirement.label == "Human Genetic Association"
    assert requirement.capability_ids == [CapabilityId.human_genetic_association]


def test_gemini_planner_schema_preserves_canonical_model_contract() -> None:
    canonical = ScientificQuestionPlanDraft.model_json_schema()
    gemini_schema = build_gemini_planner_schema()

    assert build_gemini_planner_schema() == gemini_schema
    assert canonical["properties"]["evidence_requirements"]["maxItems"] == 10
    assert "maxItems" not in gemini_schema["properties"]["evidence_requirements"]
    assert ScientificQuestionPlanDraft.model_json_schema() == canonical
    assert (
        gemini_schema["$defs"]["EvidenceRequirement"]["properties"]["genes"][
            "maxItems"
        ]
        == 6
    )
    assert (
        gemini_schema["$defs"]["EvidenceRequirement"]["properties"]["capability_ids"][
            "maxItems"
        ]
        == 8
    )


def test_nvidia_nim_uses_canonical_structured_planner_model(monkeypatch) -> None:
    captured: dict[str, Any] = {}
    question = "Compare MSH3 and FAN1 for HD modifier relevance."
    draft = ScientificQuestionPlanDraft(
        intent=ScientificIntent.comparison,
        entities=ScientificEntities(genes=["MSH3", "FAN1"]),
        primary_gene=None,
        objective="Compare direct human genetic modifier evidence.",
        analysis_lens="hd_modifier_relevance",
        answer_mode=AnswerMode.comparison,
        evidence_requirements=[
            _requirement(
                EvidenceNeed.human_genetic_association,
                ["MSH3", "FAN1"],
                CapabilityId.hd_literature,
                minimum=1,
            )
        ],
        requires_multi_gene=True,
    )

    class FakeStructured:
        def invoke(self, prompt: str) -> ScientificQuestionPlanDraft:
            captured["prompt"] = prompt
            return draft

    class FakeNim:
        def with_structured_output(self, schema: Any) -> FakeStructured:
            captured["schema"] = schema
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **_kwargs: [SimpleNamespace(provider="nvidia_nim", model=FakeNim())],
    )
    settings = _settings(Path("/tmp")).model_copy(
        update={"nvidia_nim_api_key": "test-key", "default_llm_provider": "nvidia_nim"}
    )

    plan = _try_structured_llm_plan(question, "SREBF2", settings)

    assert captured["schema"] is ScientificQuestionPlanDraft
    assert captured["schema"] is not build_gemini_planner_schema()
    assert plan is not None
    assert plan.entities.genes == ["MSH3", "FAN1"]
    assert plan.intent is ScientificIntent.comparison
    assert plan.answer_mode is AnswerMode.comparison
    assert plan.analysis_lens == "hd_modifier_relevance"
    assert plan.requires_multi_gene is True
    requirement = plan.evidence_requirements[0]
    assert requirement.evidence_need is EvidenceNeed.human_genetic_association
    assert requirement.required is True
    assert requirement.genes == ["MSH3", "FAN1"]
    assert requirement.capability_ids == [CapabilityId.human_genetic_association]


def test_openai_planner_uses_canonical_strict_json_schema(monkeypatch) -> None:
    captured: dict[str, Any] = {}
    question = "Compare MSH3 and FAN1 as Huntington disease modifiers."
    draft = ScientificQuestionPlanDraft(
        intent=ScientificIntent.comparison,
        entities=ScientificEntities(genes=["MSH3", "FAN1"]),
        primary_gene=None,
        objective="Compare HD modifier evidence.",
        analysis_lens="hd_modifier_relevance",
        answer_mode=AnswerMode.comparison,
        evidence_requirements=[
            _requirement(
                EvidenceNeed.human_genetic_association,
                ["MSH3", "FAN1"],
                CapabilityId.hd_literature,
                minimum=1,
            )
        ],
        requires_multi_gene=True,
    )

    class FakeStructured:
        def invoke(self, _prompt: str) -> ScientificQuestionPlanDraft:
            return draft

    class FakeOpenAI:
        def with_structured_output(self, schema: Any, **kwargs: Any) -> FakeStructured:
            captured["schema"] = schema
            captured.update(kwargs)
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **kwargs: [
            SimpleNamespace(provider="openai", model=FakeOpenAI(), purpose=kwargs.get("purpose"))
        ],
    )
    settings = _settings(Path("/tmp")).model_copy(
        update={"openai_api_key": "test-key", "default_llm_provider": "openai"}
    )

    plan = _try_structured_llm_plan(question, None, settings)

    assert captured["schema"] is ScientificQuestionPlanDraft
    assert captured["method"] == "json_schema"
    assert captured["strict"] is True
    assert plan is not None
    assert plan.entities.genes == ["MSH3", "FAN1"]
    assert plan.primary_gene is None
    assert plan.analysis_lens == "hd_modifier_relevance"


def test_invalid_nvidia_nim_structured_output_fails_safely(monkeypatch) -> None:
    class FakeStructured:
        def invoke(self, _prompt: str) -> dict[str, Any]:
            return {"intent": "comparison", "entities": {"genes": ["MSH3"]}}

    class FakeNim:
        def with_structured_output(self, _schema: Any) -> FakeStructured:
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **_kwargs: [SimpleNamespace(provider="nvidia_nim", model=FakeNim())],
    )
    settings = _settings(Path("/tmp")).model_copy(
        update={"nvidia_nim_api_key": "test-key", "default_llm_provider": "nvidia_nim"}
    )

    assert (
        _try_structured_llm_plan(
            "Compare MSH3 and FAN1 for HD modifier relevance.",
            "SREBF2",
            settings,
        )
        is None
    )


def test_invalid_gemini_structured_output_fails_safely(monkeypatch) -> None:
    class FakeStructured:
        def invoke(self, _prompt: str) -> dict[str, Any]:
            return {"intent": "comparison", "entities": {"genes": ["MSH3"]}}

    class FakeGemini:
        def with_structured_output(self, **_kwargs: Any) -> FakeStructured:
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **_kwargs: [
            SimpleNamespace(provider="google_gemini", model=FakeGemini())
        ],
    )
    settings = _settings(Path("/tmp")).model_copy(
        update={"google_api_key": "test-key", "default_llm_provider": "google_gemini"}
    )

    assert (
        _try_structured_llm_plan(
            "Compare MSH3, FAN1, PMS2, and HTT for HD modifier relevance.",
            "SREBF2",
            settings,
        )
        is None
    )


def test_overlong_gemini_requirement_list_still_fails_canonical_validation(
    monkeypatch,
) -> None:
    requirement = {
        "id": "human_genetic_association",
        "label": "Planner-controlled label",
        "description": "Planner-controlled description.",
        "genes": ["MSH3", "FAN1", "PMS2", "HTT"],
        "evidence_need": "human_genetic_association",
        "capability_ids": ["hd_literature"],
        "required": True,
        "minimum_support": 1,
        "rationale": "Needed for direct HD modifier evidence.",
    }

    class FakeStructured:
        def invoke(self, _prompt: str) -> dict[str, Any]:
            return {
                "intent": "comparison",
                "entities": {"genes": ["MSH3", "FAN1", "PMS2", "HTT"]},
                "primary_gene": None,
                "objective": "Compare direct human genetic modifier evidence.",
                "analysis_lens": "hd_modifier_relevance",
                "answer_mode": "comparison",
                "evidence_requirements": [
                    {**requirement, "id": f"human_genetic_association_{idx}"}
                    for idx in range(11)
                ],
                "requires_multi_gene": True,
            }

    class FakeGemini:
        def with_structured_output(self, **_kwargs: Any) -> FakeStructured:
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **_kwargs: [
            SimpleNamespace(provider="google_gemini", model=FakeGemini())
        ],
    )
    settings = _settings(Path("/tmp")).model_copy(
        update={"google_api_key": "test-key", "default_llm_provider": "google_gemini"}
    )

    assert (
        _try_structured_llm_plan(
            "Compare MSH3, FAN1, PMS2, and HTT for HD modifier relevance.",
            "SREBF2",
            settings,
        )
        is None
    )


def test_non_biomedical_and_ambiguous_statuses_are_distinct() -> None:
    outside = deterministic_fallback_plan("What is the weather tomorrow?", "SREBF2")
    ambiguous = deterministic_fallback_plan("Explain this biomedical mechanism.", "SREBF2")

    assert outside.status is AnswerStatus.out_of_scope
    assert ambiguous.status is AnswerStatus.clarification_required


def test_ordinal_citation_validation_rejects_raw_and_out_of_range_ids() -> None:
    evidence_ids = ["abc123", "def456"]

    assert citations_are_valid("Supported evidence [[1]] and [[2]].", evidence_ids)
    assert not citations_are_valid("Unsupported marker [[3]].", evidence_ids)
    assert not citations_are_valid("EvidenceRecord ID: abc123 [[1]].", evidence_ids)
    assert citations_are_valid("Public identifier abcdef123456abcdef123456abcdef12 [[1]].", evidence_ids)


def _grounding_fixture() -> tuple[EvidenceRecord, ScientificQuestionPlan, EvidenceRequirementAssessment]:
    record = _record(
        gene="FAN1",
        run_id="private-run-id",
        source="PubMed",
        assertion=AssertionType.literature_summary,
        text="FAN1 is discussed in Huntington disease CAG-repeat literature.",
        source_type=SourceType.literature,
    )
    record = record.model_copy(
        update={
            "id": "private-record-id",
            "api_run_id": "private-api-id",
            "raw_artifact_id": "private-raw-id",
        }
    )
    requirement = _requirement(EvidenceNeed.hd_literature, ["FAN1"], CapabilityId.hd_literature)
    plan = _plan(["FAN1"], [requirement])
    assessment = EvidenceRequirementAssessment(
        requirement_id=requirement.id,
        gene_symbol="FAN1",
        evidence_need=EvidenceNeed.hd_literature,
        required=True,
        minimum_support=1,
        status=RequirementStatus.sufficient,
        qualifying_count=1,
        evidence_record_ids=[record.id],
        detail="Threshold met.",
    )
    return record, plan, assessment


def _fake_grounded_result(
    monkeypatch,
    *,
    fragments: list[GroundedProseFragment],
) -> tuple[Any, dict[str, Any]]:
    record, plan, assessment = _grounding_fixture()
    captured: dict[str, Any] = {}

    class FakeStructured:
        def invoke(self, prompt: str) -> GroundedProseDraft:
            captured["prompt"] = prompt
            return GroundedProseDraft(fragments=fragments)

    class FakeModel:
        def with_structured_output(self, schema: Any, **_kwargs: Any) -> FakeStructured:
            captured["schema"] = schema
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **_kwargs: [SimpleNamespace(provider="openai", model=FakeModel())],
    )
    result = try_grounded_synthesis(
        question="What links FAN1 to HD?",
        status=AnswerStatus.answered,
        plan=plan,
        records=[record],
        assessments=[assessment],
        settings=_settings(Path("/tmp")).model_copy(update={"openai_api_key": "test-key"}),
    )
    captured.update({"record": record, "plan": plan, "assessment": assessment})
    return result, captured


def test_grounded_model_schema_is_fragments_only_and_prompt_excludes_private_ids(monkeypatch) -> None:
    result, captured = _fake_grounded_result(
        monkeypatch,
        fragments=[
            GroundedProseFragment(slot_id="slot_001", text="FAN1 has supplied HD literature evidence."),
            GroundedProseFragment(slot_id="slot_002", text="The supplied literature discusses FAN1 in HD."),
        ],
    )

    fragment_properties = GroundedProseFragment.model_json_schema()["properties"]
    assert set(fragment_properties) == {"slot_id", "text"}
    assert captured["schema"] is GroundedProseDraft
    prompt = captured["prompt"]
    record = captured["record"]
    assert record.id not in prompt
    assert record.dossier_run_id not in prompt
    assert record.api_run_id not in prompt
    assert "gap_id" not in prompt
    assert result.generation_method == "grounded_llm"
    assert result.accepted_slot_count == result.requested_slot_count == 2
    assert "[[1]]" in result.summary


def test_server_assigns_slot_science_and_citations() -> None:
    record, plan, assessment = _grounding_fixture()
    slots = build_grounded_prose_slots(plan=plan, records=[record], assessments=[assessment])

    assert len(slots) == 2
    assert all(slot.gene_symbols == ("FAN1",) for slot in slots)
    assert all(slot.evidence_category is EvidenceNeed.hd_literature for slot in slots)
    assert all(slot.epistemic_role is EpistemicRole.direct_evidence for slot in slots)
    assert all(slot.citation_ordinals == (1,) for slot in slots)
    assert all(slot.language_policy is ClaimLanguagePolicy.mechanistic_support for slot in slots)


@pytest.mark.parametrize(
    ("unsafe_text", "expected_code"),
    [
        ("FAN1 has evidence [[-1]].", "model_generated_citation_marker"),
        ("FAN1 causes Huntington disease.", "unsupported_causal_language"),
        ("FAN1 decreases somatic expansion.", "unsupported_directional_language"),
        ("FAN1 evidence conflicts with MSH3 evidence.", "unauthorized_conflict_language"),
        ("FAN1 is the best overall winner.", "comparison_ranking_or_winner_language"),
    ],
)
def test_unsafe_fragment_falls_back_without_rejecting_complete_answer(
    monkeypatch,
    unsafe_text: str,
    expected_code: str,
) -> None:
    result, _captured = _fake_grounded_result(
        monkeypatch,
        fragments=[
            GroundedProseFragment(slot_id="slot_001", text=unsafe_text),
            GroundedProseFragment(slot_id="slot_002", text="The supplied literature discusses FAN1 in HD."),
        ],
    )

    assert result.generation_method == "hybrid"
    assert result.accepted_slot_count == 1
    assert result.fallback_slot_count == 1
    assert expected_code in result.diagnostic_counts
    assert unsafe_text not in result.summary
    assert "[[1]]" in result.summary
    assert result.failure_type is None


def test_actual_private_id_falls_back_but_public_hex_identifier_is_allowed(monkeypatch) -> None:
    _record_value, _plan_value, _assessment = _grounding_fixture()
    private_result, _ = _fake_grounded_result(
        monkeypatch,
        fragments=[
            GroundedProseFragment(slot_id="slot_001", text="FAN1 record private-record-id supports HD."),
            GroundedProseFragment(slot_id="slot_002", text="The supplied literature discusses FAN1 in HD."),
        ],
    )
    public_hex = "abcdef123456abcdef123456abcdef12"
    public_result, _ = _fake_grounded_result(
        monkeypatch,
        fragments=[
            GroundedProseFragment(slot_id="slot_001", text=f"FAN1 public identifier {public_hex} is reported."),
            GroundedProseFragment(slot_id="slot_002", text="The supplied literature discusses FAN1 in HD."),
        ],
    )

    assert private_result.generation_method == "hybrid"
    assert "private_evidence_record_id_exposed" in private_result.diagnostic_counts
    assert public_result.generation_method == "grounded_llm"
    assert public_hex in public_result.summary


def test_missing_duplicate_and_unknown_slots_are_handled_independently(monkeypatch) -> None:
    result, _captured = _fake_grounded_result(
        monkeypatch,
        fragments=[
            GroundedProseFragment(slot_id="slot_001", text="First duplicate."),
            GroundedProseFragment(slot_id="slot_001", text="Second duplicate."),
            GroundedProseFragment(slot_id="slot_unknown", text="Ignore this fragment."),
        ],
    )

    assert result.generation_method == "deterministic"
    assert result.fallback_slot_count == result.requested_slot_count == 2
    assert result.diagnostic_counts == {"duplicate_slot_id": 1, "missing_slot": 1}
    assert "Ignore this fragment" not in result.summary
    assert result.failure_type is None


def test_final_validation_is_structural_only() -> None:
    record, _plan_value, _assessment = _grounding_fixture()
    registry = [
        SimpleNamespace(ordinal=1, evidence_record_id=record.id, source_id=record.source_id, source_name=record.source_name)
    ]

    assert validate_rendered_answer(
        "A causal directional conflict claim is structurally cited [[1]].",
        citation_registry=registry,
        actual_record_ids={record.id},
        rendered_slot_ids=["slot_001"],
        expected_slot_ids=["slot_001"],
    ) == []


def test_fragment_validator_uses_server_policy_not_model_metadata() -> None:
    slot = GroundedProseSlot(
        slot_id="slot_001",
        section=ProseSection.evidence,
        gene_symbols=("FAN1",),
        evidence_category=EvidenceNeed.hd_literature,
        epistemic_role=EpistemicRole.supporting_evidence,
        citation_ordinals=(1,),
        language_policy=ClaimLanguagePolicy.association_only,
        evidence_text=("FAN1 is discussed in HD literature.",),
        record_ids=("private-id",),
        fallback_text="FAN1 is discussed in HD literature.",
        stable_order=1,
        heading_label="Evidence",
    )
    fragment = GroundedProseFragment(slot_id="slot_001", text="FAN1 causes HD.")

    issues = validate_prose_fragment(
        fragment,
        slot=slot,
        all_plan_genes=("FAN1",),
        actual_record_ids={"private-id"},
    )

    assert {issue.code for issue in issues} == {"unsupported_causal_language"}


def test_q1_usage_cost_calculation_matches_terra_rates() -> None:
    assert str(terra_usage_cost_usd(input_tokens=3451, output_tokens=1952)) == "0.030326"


def test_malformed_structured_grounded_output_falls_back_with_schema_detail(monkeypatch) -> None:
    record = _record(
        gene="MSH3",
        run_id="run1",
        source="PubMed",
        assertion=AssertionType.literature_summary,
        text="MSH3 evidence is discussed in Huntington disease literature.",
        source_type=SourceType.literature,
    )
    requirement = _requirement(EvidenceNeed.hd_literature, ["MSH3"], CapabilityId.hd_literature)
    plan = _plan(["MSH3"], [requirement])
    assessment = EvidenceRequirementAssessment(
        requirement_id=requirement.id,
        gene_symbol="MSH3",
        evidence_need=EvidenceNeed.hd_literature,
        required=True,
        minimum_support=1,
        status=RequirementStatus.sufficient,
        qualifying_count=1,
        evidence_record_ids=[record.id],
        detail="Threshold met.",
    )

    class FakeStructured:
        def invoke(self, _prompt: str) -> dict[str, str]:
            return {"not_direct_answer": "bad"}

    class FakeModel:
        def with_structured_output(self, *_args: Any, **_kwargs: Any) -> FakeStructured:
            return FakeStructured()

    monkeypatch.setattr(
        "gene_dossier.synthesis.build_chat_model_candidates",
        lambda _settings, **_kwargs: [SimpleNamespace(provider="openai", model=FakeModel())],
    )

    result = try_grounded_synthesis(
        question="What links MSH3 to HD?",
        status=AnswerStatus.answered,
        plan=plan,
        records=[record],
        assessments=[assessment],
        settings=_settings(Path("/tmp")).model_copy(update={"openai_api_key": "test-key"}),
    )

    assert result.generation_method == "deterministic"
    assert result.failure_type == "schema_failure"
    assert result.validation_issues[0].code == "malformed_structured_answer"


def test_gap_ids_are_stable_and_recommendations_reference_existing_gaps(
    tmp_path: Path,
) -> None:
    requirement = _requirement(
        EvidenceNeed.human_genetic_association,
        ["MSH3"],
        CapabilityId.human_genetic_association,
    )
    service = _service(tmp_path, planner=_planner(_plan(["MSH3"], [requirement])), baselines={})
    result = service.execute(
        ScientificAgentRequest(
            question="What evidence is missing for MSH3 as an HD modifier?",
            context_gene=None,
            allow_tool_acquisition=False,
        )
    )
    gap_ids = [gap.id for gap in result.structured_gaps]

    assert gap_ids == sorted(gap_ids)
    assert gap_ids
    assert all(gap_id.startswith("msh3:") for gap_id in gap_ids)
    assert {
        gap_id
        for recommendation in result.recommendations
        for gap_id in recommendation.gap_ids
    } <= set(gap_ids)
    assert all(recommendation.label == "Recommendation" for recommendation in result.recommendations)


def test_zero_base_bootstrap_then_exact_run_reuse(tmp_path: Path) -> None:
    requirement = _requirement(EvidenceNeed.protein_interaction, ["NOVEL1"], CapabilityId.ppi)
    plan = _plan(["NOVEL1"], [requirement])
    calls: list[list[str]] = []
    holder: dict[str, ScientificAgentService] = {}

    def section_executor(gene: str, *, section_keys: list[str], **_kwargs: Any):
        calls.append(list(section_keys))
        run_id = f"novel-run-{len(calls)}"
        records = [
            _record(gene=gene, run_id=run_id, source="NCBI Gene", assertion=AssertionType.gene_identity, text=f"{gene} identity record."),
            _record(gene=gene, run_id=run_id, source="STRING", assertion=AssertionType.ppi, text=f"{gene} interaction partner A.", source_type=SourceType.interaction_database),
            _record(gene=gene, run_id=run_id, source="BioGRID", assertion=AssertionType.ppi, text=f"{gene} interaction partner B.", source_type=SourceType.interaction_database),
        ]
        _persist_run(holder["service"], gene, run_id, records)
        return SimpleNamespace(status="completed", dossier_run_id=run_id, errors=[])

    service = _service(tmp_path, planner=_planner(plan), section_executor=section_executor)
    holder["service"] = service
    first = service.execute(ScientificAgentRequest(question="What proteins interact with NOVEL1?", context_gene="SREBF2"))

    assert first.status is AnswerStatus.answered
    assert calls == [["1a", "5a", "5b"]]
    assert first.evidence_universes["NOVEL1"].created_tool_run_ids == ["novel-run-1"]
    assert first.evidence_universes["NOVEL1"].evidence_universe == "tool_overlay_only"

    def must_not_execute(*_args: Any, **_kwargs: Any):
        raise AssertionError("qualifying tagged evidence should have been reused")

    second_service = _service(tmp_path, planner=_planner(plan), section_executor=must_not_execute)
    second = second_service.execute(ScientificAgentRequest(question="Which proteins interact with NOVEL1?", context_gene="CDH10"))

    assert second.status is AnswerStatus.answered
    assert second.evidence_universes["NOVEL1"].reused_tool_run_ids == ["novel-run-1"]
    assert second.evidence_universes["NOVEL1"].created_tool_run_ids == []
    assert second.evidence_universes["NOVEL1"].tool_run_ids == ["novel-run-1"]
    assert all(item.reused for item in second.tool_activity)

    refresh_service = _service(tmp_path, planner=_planner(plan), section_executor=section_executor)
    holder["service"] = refresh_service
    refreshed = refresh_service.execute(
        ScientificAgentRequest(
            question="Refresh the proteins interacting with NOVEL1.",
            context_gene="CDH10",
            refresh_if_available=True,
        )
    )

    assert refreshed.evidence_universes["NOVEL1"].reused_tool_run_ids == []
    assert refreshed.evidence_universes["NOVEL1"].created_tool_run_ids == ["novel-run-2"]
    assert calls[-1] == ["1a", "5a", "5b"]


def test_capability_verification_requires_capability_owned_records() -> None:
    pubmed = _record(
        gene="FAN1",
        run_id="run",
        source="PubMed",
        assertion=AssertionType.literature_summary,
        text="FAN1 Huntington disease somatic CAG repeat expansion study.",
        source_type=SourceType.literature,
    )
    pathway = _record(
        gene="FAN1",
        run_id="run",
        source="Reactome",
        assertion=AssertionType.pathway_membership,
        text="FAN1 DNA repair pathway in repeat expansion biology.",
        source_type=SourceType.pathway_database,
    )
    expression = _record(
        gene="FAN1",
        run_id="run",
        source="GTEx",
        assertion=AssertionType.expression,
        text="FAN1 expression in cortex.",
        source_type=SourceType.expression_database,
    )

    assert record_matches_capability(pubmed, CapabilityId.hd_literature)
    assert not record_matches_capability(pubmed, CapabilityId.pathway)
    assert record_matches_capability(pathway, CapabilityId.pathway)
    assert not record_matches_capability(pathway, CapabilityId.hd_literature)
    assert not record_matches_capability(expression, CapabilityId.ppi)
    assert ScientificAgentService._verified_capabilities(
        [CapabilityId.hd_literature, CapabilityId.pathway],
        "FAN1",
        [pubmed],
        include_identity=False,
    ) == [CapabilityId.hd_literature]


def test_reusable_run_rejects_falsely_tagged_unrelated_capability(tmp_path: Path) -> None:
    requirement = _requirement(
        EvidenceNeed.repeat_instability_mechanism,
        ["FAN1"],
        CapabilityId.pathway,
        minimum=1,
    )
    service = _service(tmp_path)
    run_id = "false-pathway-manifest"
    pubmed = _record(
        gene="FAN1",
        run_id=run_id,
        source="PubMed",
        assertion=AssertionType.literature_summary,
        text="FAN1 Huntington disease somatic CAG repeat expansion study.",
        source_type=SourceType.literature,
    )
    _persist_run(service, "FAN1", run_id, [pubmed])
    service._tag_acquisition_run(
        gene="FAN1",
        run_id=run_id,
        executor_kind="source_workflow",
        capabilities=[CapabilityId.pathway],
        section_keys=[],
        sources=["Reactome"],
        status="completed",
        successful=True,
    )

    reusable = service._find_reusable_run(
        gene="FAN1",
        requirement=requirement,
        excluded_run_ids=set(),
    )

    assert reusable is None


def test_manifest_verifies_executed_capability_not_triggering_requirement(tmp_path: Path) -> None:
    requirement = EvidenceRequirement(
        id="repeat_mechanism",
        label="Repeat instability mechanism",
        description="Mechanistic repeat-instability evidence.",
        genes=["SREBF2"],
        evidence_need=EvidenceNeed.repeat_instability_mechanism,
        capability_ids=[CapabilityId.hd_literature],
        required=True,
        minimum_support=1,
        rationale="HD literature may contribute to mechanistic evidence.",
    )
    plan = _plan(["SREBF2"], [requirement])
    holder: dict[str, ScientificAgentService] = {}

    def source_executor(gene: str, *, sources: list[str], **_kwargs: Any):
        run_id = "srebf2-hd-literature-run"
        records = [
            _record(
                gene=gene,
                run_id=run_id,
                source="PubMed",
                assertion=AssertionType.literature_summary,
                text="SREBP2 gene therapy targeting striatal astrocytes ameliorates Huntington disease phenotypes.",
                source_type=SourceType.literature,
                grade=EvidenceGrade.F,
            )
        ]
        _persist_run(holder["service"], gene, run_id, records)
        return SimpleNamespace(
            status="completed",
            dossier_run_id=run_id,
            errors=["Ensembl/lookup_symbol_human: HTTP 503"],
        )

    service = _service(tmp_path, planner=_planner(plan), source_executor=source_executor)
    holder["service"] = service

    result = service.execute(
        ScientificAgentRequest(
            question="How does SREBF2 connect to repeat instability?",
            context_gene="SREBF2",
            evidence_selection=EvidenceSelection.accepted_or_latest_generated,
        )
    )

    assert result.evidence_universes["SREBF2"].created_tool_run_ids == ["srebf2-hd-literature-run"]
    assert result.assessments[0].status is RequirementStatus.missing
    assert result.tool_activity[0].errors == ["Ensembl/lookup_symbol_human: HTTP 503"]
    run = service._run_for_gene("srebf2-hd-literature-run", "SREBF2")
    assert run is not None
    manifest = service._manifest(run)
    assert manifest is not None
    assert manifest["capability_ids"] == ["hd_literature"]
    assert manifest["status"] == "completed_with_errors"
    assert manifest["successful"] is True
    assert manifest["had_errors"] is True

    hd_requirement = _requirement(
        EvidenceNeed.hd_literature,
        ["SREBF2"],
        CapabilityId.hd_literature,
        minimum=1,
    )
    assert service._find_reusable_run(
        gene="SREBF2",
        requirement=hd_requirement,
        excluded_run_ids=set(),
    ) == ("srebf2-hd-literature-run", [CapabilityId.hd_literature])

    unrelated_requirement = _requirement(
        EvidenceNeed.pathway_membership,
        ["SREBF2"],
        CapabilityId.pathway,
        minimum=1,
    )
    assert service._find_reusable_run(
        gene="SREBF2",
        requirement=unrelated_requirement,
        excluded_run_ids=set(),
    ) is None


def test_partial_generated_dossier_does_not_imply_missing_capability_exists(tmp_path: Path) -> None:
    requirement = _requirement(EvidenceNeed.protein_interaction, ["PARTIAL1"], CapabilityId.ppi)
    service = _service(tmp_path, planner=_planner(_plan(["PARTIAL1"], [requirement])))
    run_id = "partial-generated-run"
    _persist_run(
        service,
        "PARTIAL1",
        run_id,
        [_record(gene="PARTIAL1", run_id=run_id, source="NCBI Gene", assertion=AssertionType.gene_identity, text="PARTIAL1 identity evidence.")],
    )
    with session_scope(service.engine) as session:
        save_generated_report(
            session,
            GeneratedReportRow(
                id=canonical_generated_report_id(run_id),
                dossier_run_id=run_id,
                gene_symbol="PARTIAL1",
                title="Partial PARTIAL1 dossier",
                status="Completed",
                created_at=datetime.now(timezone.utc),
                sections=["Section 1a"],
                html_path=f"section_validation/PARTIAL1/{run_id}/section_1.html",
            ),
        )

    result = service.execute(
        ScientificAgentRequest(
            question="What proteins interact with PARTIAL1?",
            context_gene="PARTIAL1",
            evidence_selection=EvidenceSelection.accepted_or_latest_generated,
            allow_tool_acquisition=False,
        )
    )

    assert result.evidence_universes["PARTIAL1"].evidence_universe == "latest_generated"
    assert result.status is AnswerStatus.insufficient_evidence
    assert result.assessments[0].qualifying_count == 0


def test_explicit_run_ownership_mismatch_is_rejected(tmp_path: Path) -> None:
    requirement = _requirement(EvidenceNeed.identity_function, ["RIGHT1"], CapabilityId.identity_function, minimum=1)
    service = _service(tmp_path, planner=_planner(_plan(["RIGHT1"], [requirement])))
    _persist_run(
        service,
        "WRONG1",
        "wrong-run",
        [_record(gene="WRONG1", run_id="wrong-run", source="NCBI Gene", assertion=AssertionType.gene_identity, text="WRONG1 identity.")],
    )

    result = service.execute(
        ScientificAgentRequest(
            question="What is RIGHT1?",
            context_gene="RIGHT1",
            explicit_run_ids={"RIGHT1": ["wrong-run"]},
            evidence_selection=EvidenceSelection.explicit_only,
        )
    )

    assert result.status is AnswerStatus.clarification_required
    assert "does not belong" in result.summary


def test_hd_modifier_lens_keeps_four_independent_zero_base_universes(tmp_path: Path) -> None:
    genes = ["MSH3", "FAN1", "PMS2", "HTT"]
    seed = _requirement(EvidenceNeed.hd_literature, genes, CapabilityId.hd_literature)
    plan = _plan(genes, [seed], lens="hd_modifier_relevance")
    holder: dict[str, ScientificAgentService] = {}
    section_calls: list[tuple[str, list[str]]] = []
    source_calls: list[tuple[str, list[str]]] = []

    def section_executor(gene: str, *, section_keys: list[str], **_kwargs: Any):
        section_calls.append((gene, list(section_keys)))
        run_id = f"section-{gene}"
        records = [
            _record(gene=gene, run_id=run_id, source="NCBI Gene", assertion=AssertionType.gene_identity, text=f"{gene} identity record."),
            _record(gene=gene, run_id=run_id, source="GEO", assertion=AssertionType.perturbation, text=f"{gene} GEO experimental perturbation evidence.", source_type=SourceType.expression_database, grade=EvidenceGrade.D),
            _record(gene=gene, run_id=run_id, source="Allen Brain Atlas", assertion=AssertionType.expression, text=f"{gene} expression in human brain cortex.", source_type=SourceType.expression_database, grade=EvidenceGrade.B),
        ]
        _persist_run(holder["service"], gene, run_id, records)
        return SimpleNamespace(status="completed", dossier_run_id=run_id, errors=[])

    def source_executor(gene: str, *, sources: list[str], **_kwargs: Any):
        source_calls.append((gene, list(sources)))
        run_id = f"source-{gene}"
        records = [
            _record(gene=gene, run_id=run_id, source="PubMed", assertion=AssertionType.literature_summary, text=f"{gene} Huntington disease CAG repeat expansion study A.", source_type=SourceType.literature, grade=EvidenceGrade.F),
            _record(gene=gene, run_id=run_id, source="PubMed", assertion=AssertionType.literature_summary, text=f"{gene} somatic CAG repeat instability study B.", source_type=SourceType.literature, grade=EvidenceGrade.F),
            _record(gene=gene, run_id=run_id, source="Reactome", assertion=AssertionType.pathway_membership, text=f"{gene} Reactome DNA repair pathway.", source_type=SourceType.pathway_database),
        ]
        _persist_run(holder["service"], gene, run_id, records)
        return SimpleNamespace(status="completed", dossier_run_id=run_id, errors=[])

    service = _service(
        tmp_path,
        planner=_planner(plan),
        section_executor=section_executor,
        source_executor=source_executor,
    )
    holder["service"] = service
    result = service.execute(
        ScientificAgentRequest(
            question="Compare MSH3, FAN1, PMS2, and HTT for HD modifier relevance.",
            context_gene="SREBF2",
            evidence_selection=EvidenceSelection.accepted_or_latest_generated,
        )
    )

    assert set(result.evidence_universes) == set(genes)
    assert "SREBF2" not in result.evidence_universes
    assert len([item for item in result.tool_activity if not item.reused]) == 4
    assert len({capability for item in result.tool_activity for capability in item.capability_ids}) <= 4
    assert section_calls == []
    assert {gene for gene, _ in source_calls} == set(genes)
    readiness = result.metadata["readiness"]
    assert {
        row["operational_state"]
        for row in readiness
        if row["required"] is False
    } == {"not_required"}
    assert result.status is AnswerStatus.insufficient_evidence
    assert result.comparison_dimensions[0] == "Human Genetic Modifier Evidence"
    assert len(result.comparison_matrix) == 7
    assert "winner" not in result.summary.lower()


def test_hd_modifier_lens_drops_general_extras_already_covered_by_rubric(tmp_path: Path) -> None:
    genes = ["MSH3", "FAN1"]
    duplicate_expression = _requirement(
        EvidenceNeed.expression_context,
        genes,
        CapabilityId.expression_context,
        minimum=1,
        requirement_id="planner_expression_context",
    )
    duplicate_ppi = _requirement(
        EvidenceNeed.protein_interaction,
        genes,
        CapabilityId.ppi,
        minimum=1,
        requirement_id="planner_ppi_context",
    )
    distinct_extra = _requirement(
        EvidenceNeed.structure_domain,
        genes,
        CapabilityId.structure_domain,
        minimum=1,
        requirement_id="planner_structure_context",
    )
    plan = _plan(
        genes,
        [duplicate_expression, duplicate_ppi, distinct_extra],
        lens="hd_modifier_relevance",
    )
    service = _service(tmp_path)

    augmented = service._augment_hd_plan(plan)

    requirement_ids = [requirement.id for requirement in augmented.evidence_requirements]
    assert "planner_expression_context" not in requirement_ids
    assert "planner_ppi_context" not in requirement_ids
    assert "planner_structure_context" in requirement_ids
    assert len([req for req in augmented.evidence_requirements if req.id.startswith("hd_modifier_")]) == 7


def test_hd_matrix_uses_full_final_evidence_not_synthesis_context(tmp_path: Path) -> None:
    genes = ["GENE1", "GENE2", "GENE3", "GENE4"]
    seed = _requirement(EvidenceNeed.pathway_membership, genes, CapabilityId.pathway, minimum=1)
    service = _service(
        tmp_path,
        planner=_planner(_plan(genes, [seed], lens="hd_modifier_relevance")),
        baselines={gene: f"{gene.lower()}-base" for gene in genes},
    )
    for gene in genes:
        run_id = f"{gene.lower()}-base"
        records = [
            _record(
                gene=gene,
                run_id=run_id,
                source=f"Reactome-{index}",
                assertion=AssertionType.pathway_membership,
                text=f"{gene} DNA repair pathway record {index}.",
                source_type=SourceType.pathway_database,
            )
            for index in range(35 if gene == "GENE1" else 6)
        ]
        _persist_run(service, gene, run_id, records)

    result = service.execute(
        ScientificAgentRequest(
            question="Compare GENE1, GENE2, GENE3, and GENE4 for HD modifier relevance.",
            context_gene="SREBF2",
            allow_tool_acquisition=False,
        )
    )

    pathway_row = next(
        row for row in result.comparison_matrix if row.dimension == "Pathway / PPI"
    )
    assert len(result.selected_records) == 20
    assert pathway_row.cells["GENE1"].evidence_count == 35
    assert len(pathway_row.cells["GENE1"].evidence_record_ids) == 35
    assert sum(cell.evidence_count for cell in pathway_row.cells.values()) == 53
    assert all(cell.status == "Moderate" for cell in pathway_row.cells.values())


def test_hd_modifier_grading_is_dimension_aware_and_conservative() -> None:
    human_records = [
        _record(
            gene="GENE1",
            run_id="run",
            source=source,
            assertion=AssertionType.variant_association,
            text=f"GENE1 direct Huntington modifier GWAS age at onset evidence from {source}.",
            source_type=SourceType.genetic_database,
            grade=EvidenceGrade.A,
        )
        for source in ("GWAS Catalog", "GeM-HD")
    ]
    hd_literature = [
        _record(
            gene="GENE1",
            run_id="run",
            source="PubMed",
            assertion=AssertionType.literature_summary,
            text=f"GENE1 Huntington disease paper {index}.",
            source_type=SourceType.literature,
            grade=EvidenceGrade.F,
        )
        for index in range(3)
    ]
    generic_cag = [
        _record(
            gene="GENE1",
            run_id="run",
            source="PubMed",
            assertion=AssertionType.literature_summary,
            text="GENE1 was mentioned with CAG.",
            source_type=SourceType.literature,
        )
    ]
    mechanistic = [
        _record(
            gene="GENE1",
            run_id="run",
            source=source,
            assertion=assertion,
            text=f"GENE1 mismatch repair mechanism supports somatic repeat expansion {index}.",
            source_type=source_type,
        )
        for index, (source, assertion, source_type) in enumerate(
            (
                ("PubMed", AssertionType.literature_summary, SourceType.literature),
                ("Reactome", AssertionType.pathway_membership, SourceType.pathway_database),
                ("MouseMine", AssertionType.knockout_phenotype, SourceType.model_organism_database),
            )
        )
    ]
    experimental = [
        _record(
            gene="GENE1",
            run_id="run",
            source=source,
            assertion=assertion,
            text=f"GENE1 experimental perturbation {index}.",
            source_type=source_type,
        )
        for index, (source, assertion, source_type) in enumerate(
            (
                ("GEO", AssertionType.perturbation, SourceType.expression_database),
                ("MouseMine", AssertionType.knockout_phenotype, SourceType.model_organism_database),
                ("CTD", AssertionType.perturbation, SourceType.chemical_database),
            )
        )
    ]
    expression = [
        _record(
            gene="GENE1",
            run_id="run",
            source=f"Expression-{index}",
            assertion=AssertionType.expression,
            text=f"GENE1 expression in human brain cortex {index}.",
            source_type=SourceType.expression_database,
            grade=EvidenceGrade.A,
        )
        for index in range(4)
    ]
    pathways = [
        _record(
            gene="GENE1",
            run_id="run",
            source=f"Pathway-{index}",
            assertion=AssertionType.pathway_membership,
            text=f"GENE1 pathway context {index}.",
            source_type=SourceType.pathway_database,
            grade=EvidenceGrade.A,
        )
        for index in range(4)
    ]
    therapeutic = [
        _record(
            gene="GENE1",
            run_id="run",
            source=f"Chemical-{index}",
            assertion=AssertionType.chemical_tool,
            text=f"GENE1 chemical tool {index}.",
            source_type=SourceType.chemical_database,
            grade=EvidenceGrade.A,
        )
        for index in range(4)
    ]

    assert grade_hd_modifier_cell(
        EvidenceNeed.human_genetic_association,
        human_records,
        _assessment_for(EvidenceNeed.human_genetic_association, count=2),
    ) == "Strong"
    assert grade_hd_modifier_cell(
        EvidenceNeed.hd_literature,
        hd_literature,
        _assessment_for(EvidenceNeed.hd_literature),
    ) == "Strong"
    assert grade_hd_modifier_cell(
        EvidenceNeed.repeat_instability_mechanism,
        generic_cag,
        _assessment_for(EvidenceNeed.repeat_instability_mechanism, count=1),
    ) == "Weak"
    assert grade_hd_modifier_cell(
        EvidenceNeed.repeat_instability_mechanism,
        mechanistic,
        _assessment_for(EvidenceNeed.repeat_instability_mechanism),
    ) == "Strong"
    assert grade_hd_modifier_cell(
        EvidenceNeed.experimental_evidence,
        experimental,
        _assessment_for(EvidenceNeed.experimental_evidence),
    ) == "Strong"
    assert grade_hd_modifier_cell(
        EvidenceNeed.brain_expression,
        expression,
        _assessment_for(EvidenceNeed.brain_expression, count=4),
    ) == "Moderate"
    assert grade_hd_modifier_cell(
        EvidenceNeed.pathway_membership,
        pathways,
        _assessment_for(EvidenceNeed.pathway_membership, count=4),
    ) == "Moderate"
    assert grade_hd_modifier_cell(
        EvidenceNeed.therapeutic_perturbability,
        therapeutic,
        _assessment_for(EvidenceNeed.therapeutic_perturbability, count=4),
    ) == "Moderate"


def test_enabled_source_capabilities_normalize_persist_and_index_offline(tmp_path: Path) -> None:
    assert validate_source_capabilities() == []
    settings = _settings(tmp_path)
    service = ScientificAgentService(accepted_baselines={}, settings=settings, index_factory=lambda: _NoSemanticIndex())
    init_db(service.engine)
    run_id = "source-contract-run"
    tool_results = [
        ToolResult(
            source_name="PubMed",
            endpoint_name="search_hd_literature",
            success=True,
            gene_symbol="FAN1",
            request_url="https://example.invalid/pubmed",
            data={"pmids": ["123"], "search_term": "FAN1 Huntington", "esummary": {"result": {"123": {"uid": "123", "title": "FAN1 and CAG expansion"}}}},
        ),
        ToolResult(
            source_name="Reactome",
            endpoint_name="fetch_pathways",
            success=True,
            gene_symbol="FAN1",
            request_url="https://example.invalid/reactome",
            data={"uniprot_accession": "Q9Y2M0", "pathway_summaries": [{"st_id": "R-HSA-1", "display_name": "DNA repair", "species_name": "Homo sapiens"}]},
        ),
        ToolResult(
            source_name="MouseMine",
            endpoint_name="fetch_mouse_annotations",
            success=True,
            gene_symbol="FAN1",
            request_url="https://example.invalid/mousemine",
            data={"mgi_id": "MGI:1", "gene_summaries": [{"mgi_id": "MGI:1", "symbol": "Fan1"}], "allele_rows": [], "phenotype_rows": [], "stock_rows": []},
        ),
    ]
    records = [record for result in tool_results for record in normalize_tool_result(result, dossier_run_id=run_id)]
    _persist_run(service, "FAN1", run_id, records)

    with session_scope(service.engine) as session:
        persisted = list_evidence_for_run(session, run_id)
    index = ChromaEvidenceIndex(
        settings=settings,
        collection_name="source_contract_test",
        embedding_function=HashEmbeddingFunction(),
        ephemeral=True,
        allow_external_embedding_provider=False,
    )

    assert {record.source_name for record in persisted} == {"PubMed", "Reactome", "MouseMine"}
    assert all(record.api_run_id and record.raw_artifact_id for record in persisted)
    assert index.available
    assert index.upsert_evidence(persisted) == len(persisted)


def test_pubmed_controlled_source_uses_validated_hd_request_contract(monkeypatch, tmp_path: Path) -> None:
    from gene_dossier.tools import pubmed

    captured: dict[str, Any] = {}

    def fake_request(**kwargs: Any) -> ToolResult:
        captured.update(kwargs)
        return ToolResult(
            source_name="PubMed",
            endpoint_name=kwargs["endpoint_name"],
            success=True,
            gene_symbol=kwargs["gene_symbol"],
            request_url="https://example.invalid/esearch.fcgi",
            request_params=kwargs["params"],
            status_code=200,
            data={"esearchresult": {"idlist": []}},
        )

    monkeypatch.setattr(pubmed, "_request", fake_request)
    result = pubmed.esearch("FAN1", retmax=25, settings=_settings(tmp_path))

    assert result.success
    assert captured["path"] == "esearch.fcgi"
    assert captured["params"] == {
        "db": "pubmed",
        "term": '"FAN1"[Title/Abstract] AND ("Huntington Disease"[MeSH Terms] OR "Huntington disease"[Title/Abstract] OR "Huntington\'s disease"[Title/Abstract] OR "huntingtin"[Title/Abstract])',
        "retmode": "json",
        "retmax": "25",
        "sort": "relevance",
    }


def test_pubmed_gene_title_abstract_clause_uses_validated_aliases() -> None:
    from gene_dossier.tools.pubmed import build_gene_title_abstract_clause

    clause, terms = build_gene_title_abstract_clause(
        "SREBF2",
        ["SREBP2", "SREBP-2"],
    )

    assert clause == '("SREBF2"[Title/Abstract] OR "SREBP2"[Title/Abstract] OR "SREBP-2"[Title/Abstract])'
    assert terms.canonical_symbol == "SREBF2"
    assert terms.aliases_used == ("SREBP2", "SREBP-2")


def test_pubmed_gene_title_abstract_clause_deduplicates_aliases() -> None:
    from gene_dossier.tools.pubmed import build_gene_title_abstract_clause

    clause, terms = build_gene_title_abstract_clause(
        "SREBF2",
        ["srebp2", "SREBP2", "SREBF2", "HD", "SREBP-2"],
    )

    assert clause == '("SREBF2"[Title/Abstract] OR "srebp2"[Title/Abstract] OR "SREBP-2"[Title/Abstract])'
    assert terms.aliases_used == ("srebp2", "SREBP-2")


def test_pubmed_gene_title_abstract_clause_falls_back_to_canonical_only() -> None:
    from gene_dossier.tools.pubmed import build_gene_title_abstract_clause

    clause, terms = build_gene_title_abstract_clause("FAN1", [])

    assert clause == '"FAN1"[Title/Abstract]'
    assert terms.canonical_symbol == "FAN1"
    assert terms.aliases_used == ()


def test_hd_pubmed_query_includes_alias_and_hd_context_clause() -> None:
    from gene_dossier.tools.pubmed import build_hd_search_term

    query, terms = build_hd_search_term(
        "SREBF2",
        aliases=["SREBP2", "SREBP-2"],
        full_name="sterol regulatory element binding transcription factor 2",
    )

    assert '"SREBF2"[Title/Abstract]' in query
    assert '"SREBP2"[Title/Abstract]' in query
    assert '"SREBP-2"[Title/Abstract]' in query
    assert '"sterol regulatory element binding transcription factor 2"[Title/Abstract]' in query
    assert '"Huntington Disease"[MeSH Terms]' in query
    assert '"Huntington disease"[Title/Abstract]' in query
    assert '"Huntington\'s disease"[Title/Abstract]' in query
    assert '"huntingtin"[Title/Abstract]' in query
    assert terms.aliases_used == ("SREBP2", "SREBP-2")


def test_hd_pubmed_search_preserves_final_query_provenance(monkeypatch, tmp_path: Path) -> None:
    from gene_dossier.tools import pubmed

    captured: list[dict[str, Any]] = []

    def fake_request(**kwargs: Any) -> ToolResult:
        captured.append(kwargs)
        return ToolResult(
            source_name="PubMed",
            endpoint_name=kwargs["endpoint_name"],
            success=True,
            gene_symbol=kwargs["gene_symbol"],
            request_url="https://example.invalid/esearch.fcgi",
            request_params=kwargs["params"],
            status_code=200,
            data={"esearchresult": {"idlist": [], "count": "0"}},
        )

    monkeypatch.setattr(pubmed, "_request", fake_request)
    result = pubmed.search_hd_literature(
        "SREBF2",
        aliases=["SREBP2", "SREBP-2"],
        fetch_abstracts=False,
        settings=_settings(tmp_path),
    )

    assert result.success
    assert result.request_params["canonical_symbol"] == "SREBF2"
    assert result.request_params["aliases_used"] == ["SREBP2", "SREBP-2"]
    assert result.request_params["final_search_term"] == result.data["final_search_term"]
    assert '"SREBP2"[Title/Abstract]' in result.request_params["final_search_term"]
    assert '"Huntington Disease"[MeSH Terms]' in result.request_params["final_search_term"]


def test_chroma_collection_name_is_scoped_to_embedding_identity() -> None:
    identity_64 = EmbeddingIndexIdentity(
        provider="hash",
        model="fixture",
        dimension=64,
    )
    identity_1536 = EmbeddingIndexIdentity(
        provider="hash",
        model="fixture",
        dimension=1536,
    )

    assert collection_name_for_embedding("scientific_agent", identity_64) != collection_name_for_embedding(
        "scientific_agent",
        identity_1536,
    )
    assert "d64" in collection_name_for_embedding("scientific_agent", identity_64)
    assert "d1536" in collection_name_for_embedding("scientific_agent", identity_1536)


def test_chroma_incompatible_embedding_dimension_selects_new_collection(tmp_path: Path) -> None:
    run_id = "run-chroma"
    record = _record(
        gene="MSH3",
        run_id=run_id,
        source="PubMed",
        assertion=AssertionType.literature_summary,
        source_type=SourceType.literature,
        text="MSH3 Huntington disease CAG repeat instability evidence.",
    )
    index_64 = ChromaEvidenceIndex(
        persist_directory=tmp_path,
        collection_name="scientific_agent",
        embedding_function=HashEmbeddingFunction(64),
        allow_external_embedding_provider=False,
    )
    index_1536 = ChromaEvidenceIndex(
        persist_directory=tmp_path,
        collection_name="scientific_agent",
        embedding_function=HashEmbeddingFunction(1536),
        allow_external_embedding_provider=False,
    )

    assert index_64.available
    assert index_1536.available
    assert index_64.status.collection != index_1536.status.collection
    assert index_64.status.embedding_dimension == 64
    assert index_1536.status.embedding_dimension == 1536
    assert index_64.upsert_evidence([record]) == 1
    assert index_1536.upsert_evidence([record]) == 1
    assert "dimension" not in (index_1536.status.error or "").lower()


def test_warm_reuse_uses_tagged_runs_without_external_source_execution(tmp_path: Path) -> None:
    settings = _settings(tmp_path)
    plan = _plan(
        ["MSH3"],
        [
            _requirement(
                EvidenceNeed.repeat_instability_mechanism,
                ["MSH3"],
                CapabilityId.hd_literature,
                minimum=1,
                requirement_id="msh3_repeat",
            )
        ],
    )

    def fail_source_executor(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("source acquisition should not run for warm tagged evidence")

    service = ScientificAgentService(
        accepted_baselines={},
        settings=settings,
        planner=_planner(plan),
        index_factory=lambda: _NoSemanticIndex(),
        source_executor=fail_source_executor,
    )
    init_db(service.engine)
    record = _record(
        gene="MSH3",
        run_id="msh3-hd-lit-run",
        source="PubMed",
        assertion=AssertionType.literature_summary,
        source_type=SourceType.literature,
        text="MSH3 Huntington disease CAG repeat instability and somatic expansion evidence.",
    )
    _persist_run(service, "MSH3", "msh3-hd-lit-run", [record])
    service._tag_acquisition_run(
        gene="MSH3",
        run_id="msh3-hd-lit-run",
        executor_kind="source_workflow",
        capabilities=[CapabilityId.hd_literature],
        section_keys=[],
        sources=["PubMed"],
        status="completed",
        successful=True,
    )

    before = len(service._load_run_records("MSH3", "msh3-hd-lit-run"))
    result = service.execute(
        ScientificAgentRequest(
            question="Through what evidence-supported mechanism could MSH3 influence somatic CAG-repeat expansion in Huntington disease?",
            context_gene=None,
            evidence_selection=EvidenceSelection.accepted_or_latest_generated,
            allow_tool_acquisition=True,
        )
    )
    after = len(service._load_run_records("MSH3", "msh3-hd-lit-run"))
    summary = summarize_agent_result_runs(result)

    assert before == after == 1
    assert result.status is AnswerStatus.answered
    assert result.context_gene is None
    assert result.evidence_universes["MSH3"].reused_tool_run_ids == ["msh3-hd-lit-run"]
    assert result.tool_activity[0].reused is True
    assert summary["reusedToolRunIds"] == ["msh3-hd-lit-run"]
    assert summary["createdToolRunIds"] == []
    timings = result.metadata["timings"]
    assert timings
    assert all(value >= 0 for value in timings.values())


def test_readiness_table_prevents_unsupported_optional_and_sufficient_acquisition(tmp_path: Path) -> None:
    genes = ["FAN1"]
    required_sufficient = _requirement(
        EvidenceNeed.hd_literature,
        genes,
        CapabilityId.hd_literature,
        minimum=1,
        requirement_id="required_hd_lit",
    )
    unsupported = _requirement(
        EvidenceNeed.human_genetic_association,
        genes,
        CapabilityId.human_genetic_association,
        minimum=1,
        requirement_id="required_human_genetic",
    )
    optional = _requirement(
        EvidenceNeed.brain_expression,
        genes,
        CapabilityId.brain_expression,
        required=False,
        minimum=1,
        requirement_id="optional_brain_expression",
    )
    plan = _plan(genes, [required_sufficient, unsupported, optional])

    def fail_section_executor(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("optional supporting evidence must not launch acquisition")

    service = _service(tmp_path, planner=_planner(plan), section_executor=fail_section_executor)
    record = _record(
        gene="FAN1",
        run_id="fan1-hd-lit-run",
        source="PubMed",
        assertion=AssertionType.literature_summary,
        source_type=SourceType.literature,
        text="FAN1 Huntington disease CAG repeat instability literature evidence.",
    )
    _persist_run(service, "FAN1", "fan1-hd-lit-run", [record])
    service._tag_acquisition_run(
        gene="FAN1",
        run_id="fan1-hd-lit-run",
        executor_kind="source_workflow",
        capabilities=[CapabilityId.hd_literature],
        section_keys=[],
        sources=["PubMed"],
        status="completed",
        successful=True,
    )

    result = service.execute(
        ScientificAgentRequest(
            question="Compare FAN1 as an HD modifier.",
            context_gene=None,
            evidence_selection=EvidenceSelection.accepted_or_latest_generated,
            allow_tool_acquisition=True,
        )
    )
    readiness = {
        (row["requirement_id"], row["gene"]): row
        for row in result.metadata["readiness"]
    }

    assert readiness[("required_hd_lit", "FAN1")]["operational_state"] == "sufficient_persisted_evidence"
    assert readiness[("required_hd_lit", "FAN1")]["acquisition_needed"] is False
    assert readiness[("required_human_genetic", "FAN1")]["operational_state"] == "unsupported_capability"
    assert readiness[("required_human_genetic", "FAN1")]["acquisition_needed"] is False
    assert readiness[("optional_brain_expression", "FAN1")]["operational_state"] == "not_required"
    assert readiness[("optional_brain_expression", "FAN1")]["acquisition_needed"] is False
    assert result.evidence_universes["FAN1"].reused_tool_run_ids == ["fan1-hd-lit-run"]
    assert all(item.reused for item in result.tool_activity)


def test_previous_source_failure_is_not_silently_treated_as_missing(tmp_path: Path) -> None:
    requirement = _requirement(
        EvidenceNeed.hd_literature,
        ["FAN1"],
        CapabilityId.hd_literature,
        minimum=1,
        requirement_id="fan1_hd_lit",
    )
    service = _service(tmp_path, planner=_planner(_plan(["FAN1"], [requirement])))
    now = datetime.now(timezone.utc)
    with session_scope(service.engine) as session:
        save_dossier_run(
            session,
            DossierRun(
                id="fan1-failed-pubmed-run",
                gene_symbol="FAN1",
                official_symbol="FAN1",
                run_type="scientific_agent_test",
                status="completed",
                started_at=now,
                completed_at=now,
                config={},
            ),
        )
        save_api_run(
            session,
            ApiRun(
                id="api-fan1-pubmed-failed",
                dossier_run_id="fan1-failed-pubmed-run",
                gene_symbol="FAN1",
                source_name="PubMed",
                endpoint_name="offline_fixture",
                request_url="https://example.invalid/pubmed",
                success=False,
                status_code=504,
                error_type="timeout",
                error_message="offline timeout",
            ),
        )

    def fail_source_executor(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("previous source failure should block silent reacquisition")

    service.source_executor = fail_source_executor
    result = service.execute(
        ScientificAgentRequest(
            question="What evidence links FAN1 to Huntington disease?",
            context_gene=None,
            evidence_selection=EvidenceSelection.accepted_or_latest_generated,
            allow_tool_acquisition=True,
        )
    )

    row = result.metadata["readiness"][0]
    assert row["operational_state"] == "previous_source_failure"
    assert row["acquisition_possible"] is True
    assert row["acquisition_needed"] is False
    assert "PubMed" in row["source_freshness_eligibility"]
    assert result.tool_activity == []


def test_source_workflow_activity_reports_all_persisted_adapters(tmp_path: Path) -> None:
    requirement = _requirement(
        EvidenceNeed.hd_literature,
        ["FAN1"],
        CapabilityId.hd_literature,
        minimum=1,
        requirement_id="fan1_hd_lit",
    )
    plan = _plan(["FAN1"], [requirement])
    holder: dict[str, ScientificAgentService] = {}

    def source_executor(gene: str, *, sources: list[str], **_kwargs: Any):
        run_id = "fan1-composite-source-run"
        records = [
            _record(
                gene=gene,
                run_id=run_id,
                source="NCBI Gene",
                assertion=AssertionType.gene_identity,
                text="FAN1 identity evidence.",
            ),
            _record(
                gene=gene,
                run_id=run_id,
                source="UniProt",
                assertion=AssertionType.protein_function,
                text="FAN1 protein function evidence.",
            ),
            _record(
                gene=gene,
                run_id=run_id,
                source="PubMed",
                assertion=AssertionType.literature_summary,
                source_type=SourceType.literature,
                text="FAN1 Huntington disease CAG repeat literature evidence.",
            ),
        ]
        _persist_run(holder["service"], gene, run_id, records)
        return SimpleNamespace(status="completed", dossier_run_id=run_id, errors=[])

    service = _service(tmp_path, planner=_planner(plan), source_executor=source_executor)
    holder["service"] = service

    result = service.execute(
        ScientificAgentRequest(
            question="What evidence links FAN1 to Huntington disease?",
            context_gene=None,
            evidence_selection=EvidenceSelection.accepted_or_latest_generated,
            allow_tool_acquisition=True,
        )
    )

    activity = [item for item in result.tool_activity if not item.reused][0]
    assert activity.capability_ids == [CapabilityId.hd_literature]
    assert activity.sources == ["NCBI Gene", "PubMed", "UniProt"]


def test_scientific_fingerprint_ignores_provenance_and_audit_detects_duplicates() -> None:
    previous = _record(
        gene="FAN1",
        run_id="old-run",
        source="PubMed",
        assertion=AssertionType.literature_summary,
        text="FAN1 Huntington disease CAG repeat evidence.",
        source_type=SourceType.literature,
    )
    exact_new = previous.model_copy(update={"id": new_id(), "dossier_run_id": "new-run"})
    distinct_new = previous.model_copy(
        update={
            "id": new_id(),
            "dossier_run_id": "new-run",
            "display_text": "FAN1 unrelated distinct literature statement.",
            "value": {"text": "FAN1 unrelated distinct literature statement."},
        }
    )

    assert scientific_fingerprint(previous) == scientific_fingerprint(exact_new)
    audit = audit_evidence_overlap([previous], [exact_new, distinct_new])

    assert audit["exactDuplicateCount"] == 1
    assert previous.source_id in audit["exactDuplicateSourceIds"]
    assert previous.source_id in audit["semanticallyDistinctSharedSourceIds"]
    assert audit["uniqueNewCount"] == 1


def test_direct_answer_slot_prefers_stronger_evidence_grade() -> None:
    weak_requirement = _requirement(
        EvidenceNeed.identity_function,
        ["GENE1"],
        CapabilityId.identity_function,
        requirement_id="first_requirement",
    )
    strong_requirement = _requirement(
        EvidenceNeed.hd_literature,
        ["GENE1"],
        CapabilityId.hd_literature,
        requirement_id="second_requirement",
    )
    weak = _record(
        gene="GENE1",
        run_id="run",
        source="Source E",
        assertion=AssertionType.gene_identity,
        text="Lower-grade identity evidence.",
        grade=EvidenceGrade.E,
    )
    strong = _record(
        gene="GENE1",
        run_id="run",
        source="Source A",
        assertion=AssertionType.literature_summary,
        text="Higher-grade HD literature evidence.",
        source_type=SourceType.literature,
        grade=EvidenceGrade.A,
    )
    assessments = [
        EvidenceRequirementAssessment(
            requirement_id=weak_requirement.id,
            gene_symbol="GENE1",
            evidence_need=weak_requirement.evidence_need,
            required=True,
            minimum_support=1,
            status=RequirementStatus.sufficient,
            qualifying_count=1,
            evidence_record_ids=[weak.id],
            detail="Threshold met.",
        ),
        EvidenceRequirementAssessment(
            requirement_id=strong_requirement.id,
            gene_symbol="GENE1",
            evidence_need=strong_requirement.evidence_need,
            required=True,
            minimum_support=1,
            status=RequirementStatus.sufficient,
            qualifying_count=1,
            evidence_record_ids=[strong.id],
            detail="Threshold met.",
        ),
    ]

    slots = build_grounded_prose_slots(
        plan=_plan(["GENE1"], [weak_requirement, strong_requirement]),
        records=[weak, strong],
        assessments=assessments,
    )

    assert slots[0].section is ProseSection.direct_answer
    assert slots[0].record_ids == (strong.id,)


def test_gap_driven_recommendation_has_no_unrelated_citations(tmp_path: Path) -> None:
    requirement = _requirement(
        EvidenceNeed.human_genetic_association,
        ["FAN1"],
        CapabilityId.human_genetic_association,
    )
    service = _service(tmp_path, planner=_planner(_plan(["FAN1"], [requirement])))

    result = service.execute(
        ScientificAgentRequest(
            question="What human genetic evidence exists for FAN1?",
            context_gene=None,
            allow_tool_acquisition=False,
        )
    )

    assert result.recommendations
    recommendation = result.recommendations[0]
    assert recommendation.description.startswith("Gap-driven recommendation:")
    assert recommendation.rationale_citation_ordinals == []


def test_read_only_service_skips_init_and_never_upserts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    settings = _settings(tmp_path)
    requirement = _requirement(EvidenceNeed.hd_literature, ["FAN1"], CapabilityId.hd_literature)
    plan = _plan(["FAN1"], [requirement])
    writer = ScientificAgentService(
        accepted_baselines={"FAN1": "fan1-base"},
        settings=settings,
        planner=_planner(plan),
        index_factory=lambda: _NoSemanticIndex(),
    )
    init_db(writer.engine)
    _persist_run(
        writer,
        "FAN1",
        "fan1-base",
        [
            _record(
                gene="FAN1",
                run_id="fan1-base",
                source="PubMed",
                assertion=AssertionType.literature_summary,
                text="FAN1 Huntington disease literature evidence.",
                source_type=SourceType.literature,
            )
        ],
    )

    class ReadOnlyIndex:
        available = True
        status = SimpleNamespace(embedding_backend="local_minilm", index_identity="test-index")

        def upsert_evidence(self, _records: list[EvidenceRecord]) -> int:
            raise AssertionError("read-only retrieval must not upsert")

        def query(self, *_args: Any, **_kwargs: Any) -> list[Any]:
            return []

    monkeypatch.setattr(
        "gene_dossier.agent.orchestrator.init_db",
        lambda _engine: (_ for _ in ()).throw(AssertionError("read-only execution must not initialize the database")),
    )
    reader = ScientificAgentService(
        accepted_baselines={"FAN1": "fan1-base"},
        settings=settings,
        planner=_planner(plan),
        index_factory=lambda: ReadOnlyIndex(),
        read_only=True,
    )

    result = reader.execute(
        ScientificAgentRequest(
            question="What evidence links FAN1 to Huntington disease?",
            context_gene=None,
            allow_tool_acquisition=False,
        )
    )

    assert result.selected_records
    assert result.retrieval_method == "keyword"
    with pytest.raises(ValueError, match="allow_tool_acquisition=False"):
        reader.execute(
            ScientificAgentRequest(
                question="Refresh FAN1 evidence.",
                context_gene="FAN1",
                allow_tool_acquisition=True,
            )
        )


def test_chroma_read_only_opens_existing_collection_without_creation(
    tmp_path: Path,
) -> None:
    calls: list[str] = []
    client_paths: list[Path] = []

    class FakeCollection:
        metadata: dict[str, str] = {}

        @staticmethod
        def count() -> int:
            return 7

    class FakeClient:
        def get_collection(self, **_kwargs: Any) -> FakeCollection:
            calls.append("get_collection")
            return FakeCollection()

        def get_or_create_collection(self, **_kwargs: Any) -> FakeCollection:
            calls.append("get_or_create_collection")
            raise AssertionError("read-only Chroma must not create collections")

    index_dir = tmp_path / "existing-index"
    index_dir.mkdir()

    def read_only_client_factory(path: Path) -> FakeClient:
        client_paths.append(path)
        return FakeClient()

    index = ChromaEvidenceIndex(
        persist_directory=index_dir,
        collection_name="read-only-agent-index",
        embedding_function=HashEmbeddingFunction(dimensions=16),
        allow_external_embedding_provider=False,
        allow_hash_fallback=False,
        read_only=True,
        read_only_client_factory=read_only_client_factory,
    )

    assert index.available
    assert index.status.indexed_count == 7
    assert calls == ["get_collection"]
    assert client_paths == [index_dir]
    assert index.upsert_evidence([]) == 0


def test_chroma_read_only_fails_closed_without_verified_client(
    tmp_path: Path,
    monkeypatch,
) -> None:
    index_dir = tmp_path / "existing-index"
    index_dir.mkdir()
    monkeypatch.setattr(
        "chromadb.PersistentClient",
        lambda **_kwargs: pytest.fail("read-only mode must not open PersistentClient"),
    )

    index = ChromaEvidenceIndex(
        persist_directory=index_dir,
        collection_name="read-only-agent-index",
        embedding_function=HashEmbeddingFunction(dimensions=16),
        allow_external_embedding_provider=False,
        allow_hash_fallback=False,
        read_only=True,
    )

    assert not index.available
    assert "cannot guarantee non-mutating access" in (index.status.error or "")


def test_chroma_read_only_does_not_create_missing_directory(tmp_path: Path) -> None:
    missing = tmp_path / "missing-index"

    index = ChromaEvidenceIndex(
        persist_directory=missing,
        collection_name="read-only-agent-index",
        embedding_function=HashEmbeddingFunction(dimensions=16),
        allow_external_embedding_provider=False,
        allow_hash_fallback=False,
        read_only=True,
    )

    assert not index.available
    assert not missing.exists()
