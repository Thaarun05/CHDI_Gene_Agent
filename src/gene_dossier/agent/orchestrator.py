"""Bounded scientific-agent orchestration over provenance-backed evidence."""

from __future__ import annotations

import logging
import time
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable

from sqlmodel import select

from gene_dossier.config import Settings, get_settings
from gene_dossier.db import (
    ApiRunRow,
    DossierRunRow,
    EvidenceRecordRow,
    SourceCoverageResultRow,
    evidence_from_row,
    get_dossier_run,
    get_engine,
    get_read_only_engine,
    init_db,
    list_generated_reports,
    save_dossier_run,
    session_scope,
)
from gene_dossier.models import DossierRun, EvidenceRecord
from gene_dossier.retrieval import (
    ChromaEvidenceIndex,
    RetrievalHit,
    build_local_minilm_embedding_function,
    search_evidence_keyword,
    vector_id_for_record,
)
from gene_dossier.section_bundle import run_section_bundle, sanitize_credentials
from gene_dossier.workflow import run_gene_dossier_full_api_pass

from .capabilities import (
    CAPABILITY_REGISTRY,
    NEED_CONTRIBUTORS,
    acquisition_capabilities,
    record_matches_capability,
    record_matches_need,
    validated_capability_ids,
    validate_source_capabilities,
)
from .comparison import (
    build_comparison_decision,
    build_hd_modifier_matrix,
    hd_modifier_requirements,
)
from .evidence import (
    CanonicalEvidenceGroup,
    canonicalize_requirement_evidence,
    public_evidence_reference,
    public_identifier,
    public_item_from_exclusion,
    public_item_from_group,
    public_run_reference,
    record_source_url,
    record_title,
)
from .models import (
    ActivitySummary,
    AgentEvidenceUniverse,
    AnswerSection,
    AnswerStatus,
    CapabilityId,
    CitationReference,
    CostSummary,
    EvidenceCategoryBlock,
    EvidenceGap,
    EvidenceNeed,
    EvidenceRequirement,
    EvidenceRequirementAssessment,
    EvidenceSelection,
    ExperimentRecommendation,
    PublicEvidenceItem,
    ResearchMode,
    RequirementStatus,
    ScientificAgentResult,
    ScientificFailure,
    SourceAttempt,
    ScientificQuestionPlan,
    ToolActivity,
)
from .planner import PlanResult, plan_scientific_question
from .synthesis import try_grounded_synthesis

logger = logging.getLogger(__name__)

ACQUISITION_MANIFEST_KEY = "scientific_agent_acquisition"
ACQUISITION_MANIFEST_SCHEMA = 1
MAX_GENES = 6
MAX_REQUIREMENTS = 10
MAX_DISTINCT_ACQUISITIONS = 4
MAX_TOTAL_EXECUTIONS = 8
MAX_EXECUTIONS_PER_GENE = 4
MAX_FINAL_RECORDS = 20
_AGENT_COLLECTION = "scientific_agent_minilm_l6_v2_v1"
_SEMANTIC_BACKENDS = {"local_minilm", "real"}


@dataclass
class ScientificAgentRequest:
    question: str
    context_gene: str | None = None
    evidence_selection: EvidenceSelection = EvidenceSelection.accepted_only
    explicit_run_ids: dict[str, list[str]] = field(default_factory=dict)
    explicit_tool_run_ids: dict[str, list[str]] = field(default_factory=dict)
    refresh_if_available: bool = False
    allow_tool_acquisition: bool = True
    research_mode: ResearchMode = ResearchMode.auto


@dataclass
class RequirementRetrieval:
    requirement: EvidenceRequirement
    gene: str
    hits: list[RetrievalHit]
    qualifying_records: list[EvidenceRecord]
    canonical_groups: list[CanonicalEvidenceGroup]
    contextual: list[tuple[EvidenceRecord, Any]]
    excluded: list[tuple[EvidenceRecord, Any]]
    method: str
    embedding_backend: str


@dataclass
class AcquisitionGroup:
    gene: str
    executor_kind: str
    capabilities: list[CapabilityId] = field(default_factory=list)
    requirements: list[EvidenceRequirement] = field(default_factory=list)


@dataclass(frozen=True)
class RequirementReadinessRow:
    gene: str
    requirement_id: str
    evidence_need: str
    required: bool
    minimum_support: int
    qualifying_persisted_record_count: int
    qualifying_dossier_run_refs: list[str]
    qualifying_tool_run_refs: list[str]
    source_freshness_eligibility: str
    current_assessment: str
    registered_acquisition_capability: list[str]
    registered_source_adapter: list[str]
    acquisition_possible: bool
    acquisition_needed: bool
    operational_state: str
    reason: str


def _unique(values: list[str]) -> list[str]:
    return [value for value in dict.fromkeys(item.strip() for item in values if item.strip())]


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


class ScientificAgentService:
    """Execute one non-recursive evidence-planning request."""

    def __init__(
        self,
        *,
        accepted_baselines: dict[str, str],
        settings: Settings | None = None,
        planner: Callable[..., PlanResult] = plan_scientific_question,
        section_executor: Callable[..., Any] = run_section_bundle,
        source_executor: Callable[..., Any] = run_gene_dossier_full_api_pass,
        index_factory: Callable[[], ChromaEvidenceIndex] | None = None,
        read_only: bool = False,
    ) -> None:
        self.accepted_baselines = {
            gene.upper(): run_id for gene, run_id in accepted_baselines.items()
        }
        self.settings = settings or get_settings()
        self.planner = planner
        self.section_executor = section_executor
        self.source_executor = source_executor
        self.read_only = read_only
        self.index_factory = index_factory or self._default_index
        self.engine = (
            get_read_only_engine(self.settings.database_url)
            if read_only
            else get_engine(self.settings.database_url)
        )
        self._request_indexed_universes: set[tuple[str, tuple[str, ...], str | None]] = set()

    def _default_index(self) -> ChromaEvidenceIndex:
        return ChromaEvidenceIndex(
            collection_name=_AGENT_COLLECTION,
            embedding_function=build_local_minilm_embedding_function(),
            ephemeral=False,
            allow_hash_fallback=False,
            allow_external_embedding_provider=False,
            read_only=self.read_only,
        )

    @staticmethod
    def _acquisition_specs(
        requirement: EvidenceRequirement,
        plan: ScientificQuestionPlan,
    ) -> list[Any]:
        specs = acquisition_capabilities(requirement)
        restrictions = {item.casefold() for item in plan.query_policy.source_restrictions}
        if not restrictions:
            return specs
        return [
            spec
            for spec in specs
            if spec.executor_kind == "source_workflow"
            and spec.sources
            and all(source.casefold() in restrictions for source in spec.sources)
        ]

    def _run_for_gene(self, run_id: str, gene: str) -> DossierRun | None:
        with session_scope(self.engine) as session:
            run = get_dossier_run(session, run_id)
        if run is None or run.gene_symbol.strip().upper() != gene:
            return None
        return run

    def _load_records(self, gene: str, run_ids: list[str]) -> list[EvidenceRecord]:
        if not run_ids:
            return []
        with session_scope(self.engine) as session:
            rows = session.exec(
                select(EvidenceRecordRow).where(
                    EvidenceRecordRow.gene_symbol == gene,
                    EvidenceRecordRow.dossier_run_id.in_(run_ids),
                )
            ).all()
            return [evidence_from_row(row) for row in rows]

    def _load_run_records(self, gene: str, run_id: str) -> list[EvidenceRecord]:
        return self._load_records(gene, [run_id])

    def _actual_run_sources(self, run_id: str) -> list[str]:
        with session_scope(self.engine) as session:
            rows = session.exec(
                select(ApiRunRow.source_name)
                .where(ApiRunRow.dossier_run_id == run_id)
                .order_by(ApiRunRow.source_name)
            ).all()
        return _unique([str(row) for row in rows])

    @staticmethod
    def _manifest(run: DossierRun) -> dict[str, Any] | None:
        config = run.config if isinstance(run.config, dict) else {}
        manifest = config.get(ACQUISITION_MANIFEST_KEY)
        if not isinstance(manifest, dict):
            return None
        if manifest.get("schema_version") != ACQUISITION_MANIFEST_SCHEMA:
            return None
        if (
            manifest.get("status") not in {"completed", "completed_with_errors"}
            or manifest.get("successful") is not True
        ):
            return None
        if manifest.get("executor_kind") not in {"section_bundle", "source_workflow"}:
            return None
        raw_capabilities = manifest.get("capability_ids")
        if not isinstance(raw_capabilities, list) or not raw_capabilities:
            return None
        try:
            capabilities = [CapabilityId(item).value for item in raw_capabilities]
        except ValueError:
            return None
        return {**manifest, "capability_ids": capabilities}

    def _validate_tool_overlay(self, gene: str, run_id: str) -> None:
        run = self._run_for_gene(run_id, gene)
        if run is None or run.status.lower() != "completed" or self._manifest(run) is None:
            raise ValueError(f"Tool run {run_id!r} is not a reusable acquisition run for {gene}.")

    def _latest_generated_run(self, gene: str) -> str | None:
        with session_scope(self.engine) as session:
            reports = list_generated_reports(session)
        for report in reports:
            if (
                report.gene_symbol.strip().upper() != gene
                or report.status.strip().lower() != "completed"
            ):
                continue
            run = self._run_for_gene(report.dossier_run_id, gene)
            if run is not None and run.status.lower() == "completed":
                return run.id
        return None

    def _known_genes(self) -> list[str]:
        with session_scope(self.engine) as session:
            stored = session.exec(select(DossierRunRow.gene_symbol).distinct()).all()
        return sorted(
            set(self.accepted_baselines)
            | {str(gene).strip().upper() for gene in stored if str(gene).strip()}
        )

    @staticmethod
    def _refresh_universe(universe: AgentEvidenceUniverse) -> None:
        universe.reused_tool_run_ids = _unique(universe.reused_tool_run_ids)
        universe.created_tool_run_ids = _unique(universe.created_tool_run_ids)
        universe.tool_run_ids = _unique(
            [*universe.reused_tool_run_ids, *universe.created_tool_run_ids]
        )
        universe.dossier_run_ids = _unique([*universe.explicit_run_ids, *universe.tool_run_ids])
        universe.base_evidence_ref = (
            public_run_reference(
                universe.base_evidence_run_id,
                gene_symbol=universe.gene_symbol,
            )
            if universe.base_evidence_run_id
            else None
        )
        universe.explicit_run_refs = [
            public_run_reference(run_id, gene_symbol=universe.gene_symbol)
            for run_id in universe.explicit_run_ids
        ]
        universe.reused_tool_run_refs = [
            public_run_reference(run_id, gene_symbol=universe.gene_symbol)
            for run_id in universe.reused_tool_run_ids
        ]
        universe.created_tool_run_refs = [
            public_run_reference(run_id, gene_symbol=universe.gene_symbol)
            for run_id in universe.created_tool_run_ids
        ]
        universe.tool_run_refs = [
            public_run_reference(run_id, gene_symbol=universe.gene_symbol)
            for run_id in universe.tool_run_ids
        ]
        universe.dossier_run_refs = [
            public_run_reference(run_id, gene_symbol=universe.gene_symbol)
            for run_id in universe.dossier_run_ids
        ]
        has_tools = bool(universe.tool_run_ids)
        if universe.base_evidence_run_id and universe.evidence_universe.startswith("accepted_demo"):
            universe.evidence_universe = (
                "accepted_demo_with_tool_overlay" if has_tools else "accepted_demo"
            )
        elif universe.base_evidence_run_id and universe.evidence_universe.startswith(
            "latest_generated"
        ):
            universe.evidence_universe = (
                "latest_generated_with_tool_overlay" if has_tools else "latest_generated"
            )
        elif universe.base_evidence_run_id:
            universe.evidence_universe = (
                "explicit_run_with_tool_overlay" if has_tools else "explicit_run"
            )
        else:
            universe.evidence_universe = "tool_overlay_only" if has_tools else "no_base_evidence"

    def _resolve_initial_universe(
        self,
        gene: str,
        request: ScientificAgentRequest,
    ) -> AgentEvidenceUniverse:
        explicit = _unique(request.explicit_run_ids.get(gene, []))
        requested_tools = _unique(request.explicit_tool_run_ids.get(gene, []))
        base_id: str | None = None
        universe_name = "no_base_evidence"
        selected: list[str] = []
        if explicit:
            for run_id in explicit:
                if self._run_for_gene(run_id, gene) is None:
                    raise ValueError(f"Explicit run {run_id!r} does not belong to {gene}.")
            base_id = explicit[0]
            selected = list(explicit)
            universe_name = "explicit_run"
        elif (
            request.evidence_selection is not EvidenceSelection.explicit_only
            and gene in self.accepted_baselines
        ):
            base_id = self.accepted_baselines[gene]
            selected = [base_id]
            universe_name = "accepted_demo"
        elif request.evidence_selection is EvidenceSelection.accepted_or_latest_generated:
            base_id = self._latest_generated_run(gene)
            if base_id:
                selected = [base_id]
                universe_name = "latest_generated"

        for run_id in requested_tools:
            self._validate_tool_overlay(gene, run_id)
        universe = AgentEvidenceUniverse(
            gene_symbol=gene,
            base_evidence_run_id=base_id,
            explicit_run_ids=selected,
            reused_tool_run_ids=requested_tools,
            evidence_universe=universe_name,
        )
        self._refresh_universe(universe)
        return universe

    def _retrieve_requirement(
        self,
        *,
        question: str,
        gene: str,
        requirement: EvidenceRequirement,
        universe: AgentEvidenceUniverse,
        plan: ScientificQuestionPlan,
    ) -> RequirementRetrieval:
        records = self._load_records(gene, universe.dossier_run_ids)
        scoped_requirement = requirement.model_copy(update={"genes": [gene]})
        canonical = canonicalize_requirement_evidence(
            records,
            scoped_requirement,
            gene=gene,
            query_policy=plan.query_policy,
            disease_contexts=plan.entities.diseases,
        )
        eligible = [group.canonical_record for group in canonical.qualifying]
        query = f"{question} {requirement.label} {requirement.description}"
        hits: list[RetrievalHit] = []
        method = "abstain"
        backend = "unavailable"
        try:
            index = self.index_factory()
            backend = getattr(index.status, "embedding_backend", "unavailable")
            if eligible and index.available and backend in _SEMANTIC_BACKENDS:
                index_key = (
                    gene,
                    tuple(universe.dossier_run_ids),
                    getattr(index.status, "index_identity", None),
                )
                if not self.read_only and index_key not in self._request_indexed_universes:
                    index.upsert_evidence(records)
                    self._request_indexed_universes.add(index_key)
                lookup = {vector_id_for_record(record): record for record in eligible}
                semantic_hits = index.query(
                    query,
                    gene_symbol=gene,
                    dossier_run_ids=universe.dossier_run_ids,
                    limit=max(20, requirement.minimum_support * 5),
                    record_lookup=lookup,
                )
                if semantic_hits:
                    hits.extend(semantic_hits)
                    method = "semantic"
        except Exception as exc:  # noqa: BLE001
            logger.warning("scientific-agent semantic retrieval failed: %s", exc)
            backend = "unavailable"

        if eligible:
            keyword_hits = search_evidence_keyword(eligible, query, gene_symbol=gene, limit=20)
            seen = {hit.record.id for hit in hits}
            for hit in keyword_hits:
                if hit.record.id not in seen:
                    hits.append(hit)
                    seen.add(hit.record.id)
            if hits and method == "abstain":
                method = "keyword"
            for record in eligible:
                if record.id in seen:
                    continue
                hits.append(
                    RetrievalHit(
                        record=record,
                        score=0.0,
                        method="metadata",
                        source_id=vector_id_for_record(record),
                    )
                )
                seen.add(record.id)
            if hits and method == "abstain":
                method = "metadata"
        return RequirementRetrieval(
            requirement=requirement,
            gene=gene,
            hits=hits[:20],
            qualifying_records=eligible,
            canonical_groups=list(canonical.qualifying),
            contextual=list(canonical.contextual),
            excluded=list(canonical.excluded),
            method=method,
            embedding_backend=backend,
        )

    @staticmethod
    def _assessment(
        retrieval: RequirementRetrieval,
        *,
        acquisition_supported: bool | None = None,
    ) -> EvidenceRequirementAssessment:
        requirement = retrieval.requirement
        count = len(retrieval.canonical_groups)
        if acquisition_supported is None:
            acquisition_supported = bool(acquisition_capabilities(requirement))
        if count >= requirement.minimum_support:
            status = RequirementStatus.sufficient
            detail = f"{count} qualifying relevant record(s) meet the threshold of {requirement.minimum_support}."
        elif count > 0:
            status = RequirementStatus.limited
            detail = f"{count} qualifying record(s) are below the threshold of {requirement.minimum_support}."
        elif acquisition_supported:
            status = RequirementStatus.missing
            detail = "No qualifying record is present; an approved acquisition capability exists."
        else:
            status = RequirementStatus.unsupported_capability
            detail = "No qualifying record is present and no end-to-end validated acquisition is enabled."
        return EvidenceRequirementAssessment(
            requirement_id=requirement.id,
            gene_symbol=retrieval.gene,
            evidence_need=requirement.evidence_need,
            required=requirement.required,
            minimum_support=requirement.minimum_support,
            status=status,
            qualifying_count=count,
            evidence_record_ids=[hit.record.id for hit in retrieval.hits],
            public_evidence_refs=[group.public_reference for group in retrieval.canonical_groups],
            distinct_source_count=len(
                {group.source_namespace for group in retrieval.canonical_groups}
            ),
            direct_count=sum(
                group.eligibility.designation.value == "direct"
                for group in retrieval.canonical_groups
            ),
            supporting_count=sum(
                group.eligibility.designation.value == "supporting"
                for group in retrieval.canonical_groups
            ),
            contextual_count=len(retrieval.contextual),
            excluded_count=len(retrieval.excluded),
            contributing_capability_ids=validated_capability_ids(requirement),
            detail=detail,
        )

    def _retrieve_all(
        self,
        *,
        question: str,
        plan: ScientificQuestionPlan,
        universes: dict[str, AgentEvidenceUniverse],
    ) -> list[RequirementRetrieval]:
        return [
            self._retrieve_requirement(
                question=question,
                gene=gene,
                requirement=requirement,
                universe=universes[gene],
                plan=plan,
            )
            for requirement in plan.evidence_requirements
            for gene in requirement.genes
        ]

    def _find_reusable_run(
        self,
        *,
        gene: str,
        requirement: EvidenceRequirement,
        excluded_run_ids: set[str],
        plan: ScientificQuestionPlan | None = None,
    ) -> tuple[str, list[CapabilityId]] | None:
        runs = self._find_reusable_runs(
            gene=gene,
            requirement=requirement,
            excluded_run_ids=excluded_run_ids,
            require_threshold=True,
            plan=plan,
        )
        return runs[0] if runs else None

    def _find_reusable_runs(
        self,
        *,
        gene: str,
        requirement: EvidenceRequirement,
        excluded_run_ids: set[str],
        require_threshold: bool = False,
        plan: ScientificQuestionPlan | None = None,
    ) -> list[tuple[str, list[CapabilityId]]]:
        needed_capabilities = set(validated_capability_ids(requirement))
        with session_scope(self.engine) as session:
            runs = [
                DossierRun.model_validate(row.model_dump())
                for row in session.exec(
                    select(DossierRunRow)
                    .where(DossierRunRow.gene_symbol == gene, DossierRunRow.status == "completed")
                    .order_by(DossierRunRow.completed_at.desc(), DossierRunRow.started_at.desc())
                ).all()
            ]
        scoped = requirement.model_copy(update={"genes": [gene]})
        selected: list[tuple[str, list[CapabilityId]]] = []
        selected_identities: set[str] = set()
        selected_capabilities: set[CapabilityId] = set()
        for run in runs:
            if run.id in excluded_run_ids:
                continue
            manifest = self._manifest(run)
            if manifest is None:
                continue
            manifest_capabilities = {CapabilityId(item) for item in manifest["capability_ids"]}
            matching = needed_capabilities & manifest_capabilities - selected_capabilities
            if not matching:
                continue
            canonical = canonicalize_requirement_evidence(
                self._load_run_records(gene, run.id),
                scoped,
                gene=gene,
                query_policy=plan.query_policy if plan else None,
                disease_contexts=plan.entities.diseases if plan else (),
            )
            records = [group.canonical_record for group in canonical.qualifying]
            verified_matching = {
                capability
                for capability in matching
                if any(record_matches_capability(record, capability) for record in records)
            }
            capability_identities = {
                group.identity
                for group in canonical.qualifying
                if any(
                    record_matches_capability(group.canonical_record, capability)
                    for capability in verified_matching
                )
            }
            if not capability_identities:
                continue
            selected.append((run.id, sorted(verified_matching, key=lambda item: item.value)))
            selected_capabilities.update(verified_matching)
            selected_identities.update(capability_identities)
            excluded_run_ids.add(run.id)
            if (
                len(selected_identities) >= requirement.minimum_support
                or selected_capabilities >= needed_capabilities
            ):
                break
        if require_threshold and len(selected_identities) < requirement.minimum_support:
            return []
        return selected

    def _tag_acquisition_run(
        self,
        *,
        gene: str,
        run_id: str,
        executor_kind: str,
        capabilities: list[CapabilityId],
        section_keys: list[str],
        sources: list[str],
        status: str,
        successful: bool,
        execution_status: str | None = None,
        had_errors: bool = False,
    ) -> None:
        with session_scope(self.engine) as session:
            run = get_dossier_run(session, run_id)
            if run is None or run.gene_symbol.strip().upper() != gene:
                raise ValueError(
                    "Persisted acquisition run could not be resolved for manifest tagging."
                )
            config = dict(run.config or {})
            config[ACQUISITION_MANIFEST_KEY] = {
                "schema_version": ACQUISITION_MANIFEST_SCHEMA,
                "capability_ids": [capability.value for capability in capabilities],
                "executor_kind": executor_kind,
                "section_keys": list(section_keys),
                "sources": list(sources),
                "status": status,
                "successful": successful,
                "execution_status": execution_status or status,
                "had_errors": bool(had_errors),
                "tagged_at": _utcnow().isoformat(),
            }
            save_dossier_run(session, run.model_copy(update={"config": config}))

    @staticmethod
    def _verified_capabilities(
        capabilities: list[CapabilityId],
        gene: str,
        records: list[EvidenceRecord],
        *,
        include_identity: bool,
    ) -> list[CapabilityId]:
        verified: list[CapabilityId] = []
        for capability in capabilities:
            if any(
                record.gene_symbol.upper() == gene and record_matches_capability(record, capability)
                for record in records
            ):
                verified.append(capability)
        if include_identity and any(
            record_matches_capability(record, CapabilityId.identity_function) for record in records
        ):
            verified.append(CapabilityId.identity_function)
        return list(dict.fromkeys(verified))

    def _execute_group(
        self,
        group: AcquisitionGroup,
        *,
        bootstrap_identity: bool,
    ) -> ToolActivity:
        if self.read_only:
            raise RuntimeError(
                "Controlled acquisition is disabled in read-only scientific-agent mode."
            )
        specs = [CAPABILITY_REGISTRY[capability] for capability in group.capabilities]
        section_keys = _unique([key for spec in specs for key in spec.section_keys])
        sources = _unique([source for spec in specs for source in spec.sources])
        if (
            bootstrap_identity
            and group.executor_kind == "section_bundle"
            and "1a" not in section_keys
        ):
            section_keys = ["1a", *section_keys]
        activity = ToolActivity(
            gene_symbol=group.gene,
            capability_ids=group.capabilities,
            executor_kind=group.executor_kind,
            status="started",
            section_keys=section_keys,
            sources=sources,
        )
        try:
            if group.executor_kind == "section_bundle":
                result = self.section_executor(
                    group.gene,
                    section_keys=section_keys,
                    settings=self.settings,
                    persist_db=True,
                    write_pdf=False,
                )
            else:
                result = self.source_executor(
                    group.gene,
                    settings=self.settings,
                    sources=sources,
                    call_network=True,
                    force_deterministic=True,
                    write_rancho=False,
                    write_pdf=False,
                    persist_db=True,
                )
            activity.status = str(result.status)
            activity.dossier_run_id = str(result.dossier_run_id)
            activity.public_run_ref = public_run_reference(activity.dossier_run_id)
            activity.execution_succeeded = activity.status.lower() == "completed"
            activity.errors = [
                sanitize_credentials(str(error)) for error in getattr(result, "errors", [])
            ]
            records = self._load_run_records(group.gene, activity.dossier_run_id)
            activity.evidence_records_persisted = len(records)
            actual_sources = self._actual_run_sources(activity.dossier_run_id)
            if actual_sources:
                activity.sources = actual_sources
            try:
                index = self.index_factory()
                activity.indexed_records = (
                    index.upsert_evidence(records) if index.available and records else 0
                )
            except Exception as exc:  # noqa: BLE001
                activity.errors.append(f"Evidence indexing failed: {type(exc).__name__}")
            verified = self._verified_capabilities(
                group.capabilities,
                group.gene,
                records,
                include_identity=bootstrap_identity,
            )
            completed = activity.status.lower() == "completed"
            manifest_status = (
                "completed_with_errors"
                if completed and activity.errors
                else "completed"
                if completed
                else activity.status.lower()
            )
            try:
                self._tag_acquisition_run(
                    gene=group.gene,
                    run_id=activity.dossier_run_id,
                    executor_kind=group.executor_kind,
                    capabilities=verified,
                    section_keys=section_keys,
                    sources=sources,
                    status=manifest_status,
                    successful=completed and bool(verified),
                    execution_status=activity.status.lower(),
                    had_errors=bool(activity.errors),
                )
            except Exception as exc:  # noqa: BLE001
                logger.exception(
                    "acquisition manifest tagging failed for %s", activity.dossier_run_id
                )
                activity.errors.append(
                    f"Acquisition reuse manifest was not persisted: {type(exc).__name__}"
                )
        except Exception as exc:  # noqa: BLE001
            logger.exception("controlled scientific acquisition failed for %s", group.gene)
            activity.status = "failed"
            activity.errors.append(sanitize_credentials(str(exc)))
        return activity

    def _augment_hd_plan(self, plan: ScientificQuestionPlan) -> ScientificQuestionPlan:
        requirements = list(plan.evidence_requirements)
        if plan.analysis_lens == "hd_modifier_relevance":
            rubric = hd_modifier_requirements(plan.entities.genes)
            rubric_needs = {requirement.evidence_need for requirement in rubric}
            hd_rubric_covered_needs = {
                *rubric_needs,
                # General planner needs already represented by the HD modifier lens.
                EvidenceNeed.expression_context,
                EvidenceNeed.protein_interaction,
                EvidenceNeed.model_organism,
                EvidenceNeed.chemical_perturbation,
            }
            extra = [
                requirement
                for requirement in requirements
                if requirement.evidence_need not in hd_rubric_covered_needs
            ][: max(0, MAX_REQUIREMENTS - len(rubric))]
            requirements = [*rubric, *extra]

        if plan.query_policy.ranking_requested:
            decision_needs = list(plan.query_policy.comparison_criteria)
            if EvidenceNeed.therapeutic_perturbability in decision_needs:
                decision_needs.extend(
                    [
                        EvidenceNeed.safety_tolerability,
                        EvidenceNeed.clinical_translational,
                    ]
                )
            decision_needs = list(dict.fromkeys(decision_needs))
            by_need = {requirement.evidence_need: index for index, requirement in enumerate(requirements)}
            for need in decision_needs:
                existing_index = by_need.get(need)
                if existing_index is not None:
                    requirements[existing_index] = requirements[existing_index].model_copy(
                        update={"required": True}
                    )
                    continue
                if len(requirements) >= MAX_REQUIREMENTS:
                    continue
                label = need.value.replace("_", " ").title()
                requirements.append(
                    EvidenceRequirement(
                        id=f"decision_{need.value}",
                        label=label,
                        description=f"Qualifying {label.lower()} evidence for comparison.",
                        genes=plan.entities.genes,
                        evidence_need=need,
                        capability_ids=list(NEED_CONTRIBUTORS[need]),
                        required=True,
                        minimum_support=1,
                        rationale=(
                            "Required by the deterministic comparison decision policy for the "
                            "scientist's explicit criterion."
                        ),
                    )
                )
                by_need[need] = len(requirements) - 1
        return plan.model_copy(update={"evidence_requirements": requirements})

    @staticmethod
    def _select_final_records(retrievals: list[RequirementRetrieval]) -> list[EvidenceRecord]:
        buckets = [list(retrieval.hits) for retrieval in retrievals]
        selected: list[EvidenceRecord] = []
        seen: set[str] = set()
        while buckets and len(selected) < MAX_FINAL_RECORDS:
            next_buckets: list[list[RetrievalHit]] = []
            for bucket in buckets:
                while bucket and public_evidence_reference(bucket[0].record) in seen:
                    bucket.pop(0)
                if bucket and len(selected) < MAX_FINAL_RECORDS:
                    hit = bucket.pop(0)
                    selected.append(hit.record)
                    seen.add(public_evidence_reference(hit.record))
                if bucket:
                    next_buckets.append(bucket)
            buckets = next_buckets
        return selected

    @staticmethod
    def _comparison_records(
        retrievals: list[RequirementRetrieval],
        requirements: list[EvidenceRequirement],
    ) -> list[EvidenceRecord]:
        """Return the full deduplicated final evidence union for HD rubric cells."""
        requirement_ids = {requirement.id for requirement in requirements}
        records: list[EvidenceRecord] = []
        seen: set[str] = set()
        for retrieval in retrievals:
            if retrieval.requirement.id not in requirement_ids:
                continue
            for group in retrieval.canonical_groups:
                record = group.canonical_record
                if group.public_reference in seen:
                    continue
                records.append(record)
                seen.add(group.public_reference)
        return records

    @staticmethod
    def _overall_retrieval_metadata(retrievals: list[RequirementRetrieval]) -> tuple[str, str]:
        methods = {item.method for item in retrievals}
        backends = {item.embedding_backend for item in retrievals}
        method = (
            "semantic"
            if "semantic" in methods
            else "keyword"
            if "keyword" in methods
            else "metadata"
            if "metadata" in methods
            else "abstain"
        )
        backend = (
            "local_minilm"
            if "local_minilm" in backends
            else "real"
            if "real" in backends
            else "unavailable"
        )
        return method, backend

    def _update_tool_scientific_counts(
        self,
        activities: list[ToolActivity],
        retrievals: list[RequirementRetrieval],
    ) -> None:
        """Keep execution success separate from qualifying scientific retrieval."""
        for activity in activities:
            run_id = activity.dossier_run_id
            if not run_id:
                continue
            all_records = self._load_run_records(activity.gene_symbol, run_id)
            all_refs = {public_evidence_reference(record) for record in all_records}
            qualifying_refs = {
                group.public_reference
                for retrieval in retrievals
                if retrieval.gene == activity.gene_symbol
                for group in retrieval.canonical_groups
                if any(record.dossier_run_id == run_id for record in group.backing_records)
            }
            activity.qualifying_evidence_count = len(qualifying_refs)
            activity.rejected_evidence_count = len(all_refs - qualifying_refs)
            activity.scientific_retrieval_succeeded = bool(qualifying_refs)

    @staticmethod
    def _public_evidence_items(
        retrievals: list[RequirementRetrieval],
    ) -> tuple[list[PublicEvidenceItem], list[PublicEvidenceItem]]:
        accepted: dict[str, PublicEvidenceItem] = {}
        contextual: dict[tuple[str, str], PublicEvidenceItem] = {}
        for retrieval in retrievals:
            for group in retrieval.canonical_groups:
                accepted.setdefault(group.public_reference, public_item_from_group(group))
            for record, eligibility in [*retrieval.contextual, *retrieval.excluded]:
                item = public_item_from_exclusion(
                    record,
                    eligibility,
                    evidence_need=retrieval.requirement.evidence_need,
                )
                contextual.setdefault(
                    (item.public_evidence_ref, item.exclusion_reason or "excluded"),
                    item,
                )
        accepted_items = sorted(
            accepted.values(),
            key=lambda item: (item.gene_symbol, item.evidence_need.value, item.public_evidence_ref),
        )
        contextual_items = sorted(
            contextual.values(),
            key=lambda item: (
                item.gene_symbol,
                item.evidence_need.value,
                item.exclusion_reason or "",
                item.public_evidence_ref,
            ),
        )
        return accepted_items, contextual_items

    @staticmethod
    def _activity_summary(
        *,
        plan: ScientificQuestionPlan,
        retrievals: list[RequirementRetrieval],
        activities: list[ToolActivity],
        readiness_rows: list[RequirementReadinessRow],
    ) -> ActivitySummary:
        accepted_refs = {
            group.public_reference
            for retrieval in retrievals
            for group in retrieval.canonical_groups
        }
        rejected_keys = {
            (public_evidence_reference(record), eligibility.reason_code)
            for retrieval in retrievals
            for record, eligibility in [*retrieval.contextual, *retrieval.excluded]
        }
        rejection_reasons: dict[str, int] = defaultdict(int)
        for _reference, reason in rejected_keys:
            rejection_reasons[reason] += 1
        executed = [item for item in activities if not item.reused]
        reused = [item for item in activities if item.reused]
        skip_reasons: dict[str, int] = defaultdict(int)
        for row in readiness_rows:
            if not row.acquisition_needed:
                skip_reasons[row.operational_state] += 1
        return ActivitySummary(
            requirements_planned=len(plan.evidence_requirements),
            persisted_retrieval_completed=True,
            tools_executed=len(executed),
            tools_failed=sum(not item.execution_succeeded for item in executed),
            runs_reused=len(reused),
            tools_skipped=sum(skip_reasons.values()),
            accepted_evidence=len(accepted_refs),
            rejected_evidence=len(rejected_keys),
            skip_reasons=dict(sorted(skip_reasons.items())),
            rejection_reasons=dict(sorted(rejection_reasons.items())),
        )

    def _api_retrieval_times(self, records: list[EvidenceRecord]) -> dict[str, str]:
        api_ids = _unique([record.api_run_id or "" for record in records])
        if not api_ids:
            return {}
        with session_scope(self.engine) as session:
            rows = session.exec(select(ApiRunRow).where(ApiRunRow.id.in_(api_ids))).all()
            return {row.id: row.retrieved_at.isoformat() for row in rows}

    @staticmethod
    def _evidence_system(record: EvidenceRecord) -> str:
        source_type = str(getattr(record.source_type, "value", record.source_type)).lower()
        organism = " ".join(
            str(item or "").lower()
            for item in (
                record.organism,
                record.species,
                record.value.get("organism"),
                record.value.get("species"),
            )
        )
        assertion = str(getattr(record.assertion_type, "value", record.assertion_type)).lower()
        text = " ".join((record.fact_type, record.display_text, organism)).lower()
        if source_type == "genetic_database" or any(
            term in organism for term in ("human", "homo sapiens")
        ):
            return "human"
        if source_type == "model_organism_database" or any(
            term in organism for term in ("mouse", "mus musculus", "rat", "drosophila", "zebrafish")
        ):
            return "animal/model-organism"
        if assertion in {"perturbation", "chemical_interaction", "chemical_tool"} or any(
            term in text for term in ("cell", "in vitro", "cellular", "culture")
        ):
            return "in-vitro/cellular"
        if source_type in {"structure_database", "interaction_database"} or any(
            term in text for term in ("predicted", "alphafold", "computational", "string")
        ):
            return "computational"
        return "not specified"

    @staticmethod
    def _claim_type(record: EvidenceRecord) -> str:
        assertion = str(getattr(record.assertion_type, "value", record.assertion_type)).lower()
        if assertion in {
            "variant_association",
            "disease_association",
            "ppi",
            "expression",
            "cell_type_expression",
        }:
            return "association"
        if assertion in {
            "protein_function",
            "pathway_membership",
            "protein_structure",
            "literature_summary",
            "orthology",
            "transcription_factor_association",
        }:
            return "mechanistic"
        if assertion in {
            "perturbation",
            "chemical_interaction",
            "chemical_tool",
            "knockout_phenotype",
        }:
            return "perturbational"
        return "not specified"

    @staticmethod
    def _stable_gap_id(assessment: EvidenceRequirementAssessment) -> str:
        return ":".join(
            (
                assessment.gene_symbol,
                assessment.evidence_need.value,
                "required" if assessment.required else "supporting",
                assessment.status.value,
            )
        ).lower()

    def _structured_gaps(
        self, assessments: list[EvidenceRequirementAssessment]
    ) -> list[EvidenceGap]:
        return [
            EvidenceGap(
                id=self._stable_gap_id(item),
                gene_symbol=item.gene_symbol,
                requirement_id=item.requirement_id,
                evidence_need=item.evidence_need,
                status=item.status,
                required=item.required,
                detail=item.detail,
            )
            for item in sorted(
                assessments,
                key=lambda item: (item.gene_symbol, item.requirement_id, item.evidence_need.value),
            )
            if item.status is not RequirementStatus.sufficient
        ]

    @staticmethod
    def _answer_sections(
        *,
        status: AnswerStatus,
        summary: str,
        categories: list[EvidenceCategoryBlock],
        comparison_decision: Any | None,
    ) -> list[AnswerSection]:
        sections = [
            AnswerSection(
                key="status",
                title="Evidence status",
                paragraphs=[
                    "All required evidence thresholds were met."
                    if status is AnswerStatus.answered
                    else "One or more required evidence needs remain below threshold or unsupported."
                ],
            ),
            AnswerSection(
                key="direct_answer",
                title="Direct answer",
                paragraphs=[summary],
            ),
        ]
        if comparison_decision is not None:
            sections.append(
                AnswerSection(
                    key="conditional_conclusion",
                    title="Comparison conclusion",
                    paragraphs=[comparison_decision.summary, *comparison_decision.limitations],
                )
            )
        findings = [
            f"{item.gene_symbol} · {item.evidence_need.value.replace('_', ' ')}: {item.summary}"
            for item in categories
            if item.unique_qualifying_count
        ]
        if findings:
            sections.append(
                AnswerSection(
                    key="key_findings",
                    title="Key findings",
                    paragraphs=findings[:8],
                )
            )
        dimension_lines = [
            f"{item.gene_symbol} · {item.category}: {item.unique_qualifying_count} unique item(s), "
            f"{item.distinct_source_count} source(s)."
            for item in categories
        ]
        if dimension_lines:
            sections.append(
                AnswerSection(
                    key="evidence_by_dimension",
                    title="Evidence by dimension",
                    paragraphs=dimension_lines,
                )
            )
        return sections

    def _citation_registry(
        self,
        records: list[EvidenceRecord],
        retrievals: list[RequirementRetrieval],
    ) -> list[CitationReference]:
        groups_by_record = {
            group.canonical_record.id: group
            for retrieval in retrievals
            for group in retrieval.canonical_groups
        }
        timestamps = self._api_retrieval_times(records)
        return [
            CitationReference(
                ordinal=index,
                evidence_record_id=record.id,
                public_evidence_ref=public_evidence_reference(record),
                source_id=record.source_id,
                source_name=record.source_name,
                title=record_title(record),
                public_identifier=public_identifier(record),
                source_url=record_source_url(record),
                evidence_need=(
                    groups_by_record[record.id].evidence_need
                    if record.id in groups_by_record
                    else None
                ),
                designation=(
                    groups_by_record[record.id].eligibility.designation
                    if record.id in groups_by_record
                    else "supporting"
                ),
                retrieved_at=timestamps.get(record.api_run_id or ""),
            )
            for index, record in enumerate(records, start=1)
        ]

    def _evidence_categories(
        self,
        records: list[EvidenceRecord],
        assessments: list[EvidenceRequirementAssessment],
    ) -> list[EvidenceCategoryBlock]:
        timestamps = self._api_retrieval_times(records)
        by_record = {record.id: record for record in records}
        ordinal_by_id = {record.id: index for index, record in enumerate(records, start=1)}
        blocks: list[EvidenceCategoryBlock] = []
        for assessment in sorted(
            assessments,
            key=lambda item: (item.gene_symbol, item.evidence_need.value, item.requirement_id),
        ):
            matched = [by_record[rid] for rid in assessment.evidence_record_ids if rid in by_record]
            systems = sorted({self._evidence_system(record) for record in matched}) or [
                "not specified"
            ]
            claim_types = sorted({self._claim_type(record) for record in matched}) or [
                "not specified"
            ]
            blocks.append(
                EvidenceCategoryBlock(
                    gene_symbol=assessment.gene_symbol,
                    category=assessment.evidence_need.value,
                    evidence_need=assessment.evidence_need,
                    evidence_system=", ".join(systems),
                    claim_type=", ".join(claim_types),
                    evidence_record_ids=[record.id for record in matched],
                    public_evidence_refs=[public_evidence_reference(record) for record in matched],
                    citation_ordinals=[
                        ordinal_by_id[record.id] for record in matched if record.id in ordinal_by_id
                    ],
                    source_names=sorted({record.source_name for record in matched}),
                    retrieval_timestamps=sorted(
                        {
                            timestamps[record.api_run_id]
                            for record in matched
                            if record.api_run_id in timestamps
                        }
                    ),
                    unique_qualifying_count=assessment.qualifying_count,
                    distinct_source_count=assessment.distinct_source_count,
                    direct_count=assessment.direct_count,
                    supporting_count=assessment.supporting_count,
                    status=assessment.status.value,
                    summary=assessment.detail,
                )
            )
        return blocks

    def _source_attempts(self, universes: dict[str, AgentEvidenceUniverse]) -> list[SourceAttempt]:
        run_ids = _unique(
            [run_id for universe in universes.values() for run_id in universe.dossier_run_ids]
        )
        if not run_ids:
            return []
        attempts: list[SourceAttempt] = []
        with session_scope(self.engine) as session:
            api_rows = session.exec(
                select(ApiRunRow).where(ApiRunRow.dossier_run_id.in_(run_ids))
            ).all()
            coverage_rows = session.exec(
                select(SourceCoverageResultRow).where(
                    SourceCoverageResultRow.dossier_run_id.in_(run_ids)
                )
            ).all()
            for row in api_rows:
                attempts.append(
                    SourceAttempt(
                        gene_symbol=row.gene_symbol,
                        dossier_run_id=row.dossier_run_id,
                        public_run_ref=public_run_reference(
                            row.dossier_run_id,
                            gene_symbol=row.gene_symbol,
                        ),
                        source_name=row.source_name,
                        status="success" if row.success else (row.error_type or "failed"),
                        retrieved_at=row.retrieved_at.isoformat(),
                        error_message="Source attempt failed." if row.error_message else None,
                    )
                )
            for row in coverage_rows:
                attempts.append(
                    SourceAttempt(
                        gene_symbol="",
                        dossier_run_id=row.dossier_run_id,
                        public_run_ref=public_run_reference(row.dossier_run_id),
                        source_name=row.source_name,
                        status=row.status,
                        error_message="Source attempt failed." if row.error_message else None,
                    )
                )
        return sorted(
            attempts,
            key=lambda item: (item.public_run_ref, item.source_name, item.retrieved_at or ""),
        )

    def _previous_source_failure(self, *, gene: str, sources: list[str]) -> str | None:
        wanted = {source.strip().lower() for source in sources if source.strip()}
        if not wanted:
            return None
        with session_scope(self.engine) as session:
            rows = [
                {
                    "source_name": row.source_name,
                    "success": row.success,
                    "status_code": row.status_code,
                    "error_type": row.error_type,
                }
                for row in session.exec(
                    select(ApiRunRow)
                    .where(ApiRunRow.gene_symbol == gene, ApiRunRow.source_name.in_(list(sources)))
                    .order_by(ApiRunRow.retrieved_at.desc())
                ).all()
            ]
        latest_by_source: dict[str, dict[str, Any]] = {}
        for row in rows:
            key = str(row["source_name"]).strip().lower()
            if key not in wanted or key in latest_by_source:
                continue
            latest_by_source[key] = row
        failed = [
            f"{row['source_name']}: {row['error_type'] or row['status_code'] or 'failed'}"
            for row in latest_by_source.values()
            if not row["success"]
        ]
        return "; ".join(sorted(failed)) if failed else None

    @staticmethod
    def _readiness_state(
        *,
        requirement: EvidenceRequirement,
        assessment: EvidenceRequirementAssessment,
        specs: list[Any],
        previous_failure: str | None,
    ) -> tuple[str, bool, bool, str]:
        if not requirement.required:
            reason = "Supporting requirement; disclose any gap without launching preparation acquisition."
            return "not_required", False, False, reason
        if assessment.status is RequirementStatus.sufficient:
            reason = "Configured sufficiency threshold is met by selected persisted evidence."
            return "sufficient_persisted_evidence", False, False, reason
        if not specs:
            reason = "No end-to-end validated acquisition capability is registered for this evidence need."
            return "unsupported_capability", False, False, reason
        if previous_failure:
            reason = f"Previous source failure is recorded and will not be hidden as simple missing evidence: {previous_failure}."
            return "previous_source_failure", True, False, reason
        if assessment.status is RequirementStatus.limited:
            reason = "Persisted evidence exists but is below the configured support threshold."
            return "insufficient_persisted_evidence_acquirable", True, True, reason
        reason = "No qualifying persisted evidence is present and an approved acquisition capability is registered."
        return "missing_acquirable", True, True, reason

    def _readiness_rows(
        self,
        *,
        retrievals: list[RequirementRetrieval],
        assessments: list[EvidenceRequirementAssessment],
        universes: dict[str, AgentEvidenceUniverse],
        plan: ScientificQuestionPlan,
    ) -> list[RequirementReadinessRow]:
        assessment_map = {(item.requirement_id, item.gene_symbol): item for item in assessments}
        rows: list[RequirementReadinessRow] = []
        for retrieval in retrievals:
            requirement = retrieval.requirement
            gene = retrieval.gene
            assessment = assessment_map[(requirement.id, gene)]
            specs = self._acquisition_specs(requirement, plan)
            sources = _unique(
                [
                    source
                    for spec in specs
                    for source in (
                        list(spec.sources) if spec.executor_kind == "source_workflow" else []
                    )
                ]
            )
            section_adapters = _unique(
                [
                    f"section_bundle:{','.join(spec.section_keys)}"
                    for spec in specs
                    if spec.executor_kind == "section_bundle"
                ]
            )
            adapters = [*sources, *section_adapters] or ["none"]
            previous_failure = self._previous_source_failure(gene=gene, sources=sources)
            state, possible, needed, reason = self._readiness_state(
                requirement=requirement,
                assessment=assessment,
                specs=specs,
                previous_failure=previous_failure,
            )
            qualifying_run_ids = _unique(
                [record.dossier_run_id for record in retrieval.qualifying_records]
            )
            universe = universes[gene]
            qualifying_tool_ids = [
                run_id for run_id in qualifying_run_ids if run_id in set(universe.tool_run_ids)
            ]
            source_freshness = (
                previous_failure
                if previous_failure
                else "eligible; no disqualifying source failure in selected deterministic policy"
                if possible
                else "not eligible for acquisition"
            )
            rows.append(
                RequirementReadinessRow(
                    gene=gene,
                    requirement_id=requirement.id,
                    evidence_need=requirement.evidence_need.value,
                    required=requirement.required,
                    minimum_support=requirement.minimum_support,
                    qualifying_persisted_record_count=len(retrieval.canonical_groups),
                    qualifying_dossier_run_refs=[
                        public_run_reference(run_id) for run_id in qualifying_run_ids
                    ],
                    qualifying_tool_run_refs=[
                        public_run_reference(run_id) for run_id in qualifying_tool_ids
                    ],
                    source_freshness_eligibility=source_freshness,
                    current_assessment=assessment.status.value,
                    registered_acquisition_capability=[spec.capability_id.value for spec in specs],
                    registered_source_adapter=adapters,
                    acquisition_possible=possible,
                    acquisition_needed=needed,
                    operational_state=state,
                    reason=reason,
                )
            )
        return sorted(rows, key=lambda item: (item.gene, item.requirement_id, item.evidence_need))

    def readiness_table(
        self,
        *,
        question: str,
        plan: ScientificQuestionPlan,
        request: ScientificAgentRequest,
    ) -> tuple[list[RequirementReadinessRow], dict[str, AgentEvidenceUniverse]]:
        universes = {
            gene: self._resolve_initial_universe(gene, request) for gene in plan.entities.genes
        }
        if not request.refresh_if_available:
            for requirement in sorted(
                plan.evidence_requirements, key=lambda item: (not item.required, item.id)
            ):
                for gene in requirement.genes:
                    universe = universes[gene]
                    reusable_runs = self._find_reusable_runs(
                        gene=gene,
                        requirement=requirement,
                        excluded_run_ids=set(universe.dossier_run_ids),
                        plan=plan,
                    )
                    for run_id, _capabilities in reusable_runs:
                        if run_id not in universe.reused_tool_run_ids:
                            universe.reused_tool_run_ids.append(run_id)
                            self._refresh_universe(universe)
        retrievals = self._retrieve_all(question=question, plan=plan, universes=universes)
        assessments = [
            self._assessment(
                item,
                acquisition_supported=bool(self._acquisition_specs(item.requirement, plan)),
            )
            for item in retrievals
        ]
        return self._readiness_rows(
            retrievals=retrievals,
            assessments=assessments,
            universes=universes,
            plan=plan,
        ), universes

    @staticmethod
    def _recommendations(
        gaps: list[EvidenceGap], records: list[EvidenceRecord]
    ) -> list[ExperimentRecommendation]:
        priority = {
            "human_genetic_association": 0,
            "repeat_instability_mechanism": 1,
            "therapeutic_perturbability": 2,
            "brain_expression": 3,
            "experimental_evidence": 4,
        }
        ordinal_by_id = {record.id: index for index, record in enumerate(records, start=1)}
        recommendations: list[ExperimentRecommendation] = []
        for gap in sorted(
            gaps, key=lambda item: (priority.get(item.evidence_need.value, 99), item.id)
        ):
            if len(recommendations) >= 3:
                break
            if gap.evidence_need.value == "human_genetic_association":
                description = f"Run or curate a human modifier-association analysis for {gap.gene_symbol} in HD cohorts."
                uncertainty = "Whether direct human genetic modifier evidence supports the target."
            elif gap.evidence_need.value == "repeat_instability_mechanism":
                description = f"Test {gap.gene_symbol} perturbation in a CAG-repeat instability assay with prespecified direction-of-effect readouts."
                uncertainty = "Whether the target changes somatic CAG-repeat expansion through a measurable mechanism."
            elif gap.evidence_need.value == "therapeutic_perturbability":
                description = f"Evaluate selective perturbation tools for {gap.gene_symbol} and document potency, selectivity, and assay context."
                uncertainty = (
                    "Whether the target is tractable with interpretable perturbation tools."
                )
            elif gap.evidence_need.value in {"brain_expression", "expression_context"}:
                description = f"Resolve {gap.gene_symbol} expression in HD-relevant brain regions and cell types."
                uncertainty = (
                    "Whether the target is present in the relevant tissue and cell context."
                )
            else:
                description = f"Collect targeted evidence for {gap.gene_symbol} {gap.evidence_need.value.replace('_', ' ')}."
                uncertainty = "Whether this evidence dimension is adequate for the decision."
            rationale_citation_ordinals = [
                ordinal_by_id[record.id]
                for record in records
                if record.gene_symbol.strip().upper() == gap.gene_symbol.strip().upper()
                and record_matches_need(record, gap.evidence_need)
                and record.id in ordinal_by_id
            ][:3]
            limitations = [
                "The recommendation addresses a deterministic evidence gap and does not predict the experiment's result."
            ]
            if not rationale_citation_ordinals:
                description = f"Gap-driven recommendation: {description}"
                limitations.append(
                    "This is a gap-driven recommendation; no compatible rationale EvidenceRecord was available."
                )
            recommendations.append(
                ExperimentRecommendation(
                    description=description,
                    gap_ids=[gap.id],
                    gap_labels=[f"{gap.gene_symbol} {gap.evidence_need.value.replace('_', ' ')}"],
                    decision_uncertainty=uncertainty,
                    rationale_citation_ordinals=rationale_citation_ordinals,
                    limitations=limitations,
                )
            )
        return recommendations

    def execute(self, request: ScientificAgentRequest) -> ScientificAgentResult:
        request_start = time.perf_counter()
        timings: dict[str, float] = {}

        def mark(stage: str, started: float) -> None:
            timings[stage] = round(time.perf_counter() - started, 6)

        self._request_indexed_universes = set()
        acquisition_allowed = (
            request.allow_tool_acquisition and request.research_mode is not ResearchMode.stored_only
        )
        if self.read_only and acquisition_allowed:
            raise ValueError(
                "Read-only scientific-agent mode requires allow_tool_acquisition=False."
            )
        if not self.read_only:
            init_db(self.engine)
        question = request.question.strip()
        context_gene = (request.context_gene or "").strip().upper() or None
        stage_start = time.perf_counter()
        plan_result = self.planner(
            question,
            context_gene=context_gene,
            settings=self.settings,
            known_genes=self._known_genes(),
        )
        mark("planning", stage_start)
        if plan_result.status is AnswerStatus.out_of_scope:
            timings["total_request"] = round(time.perf_counter() - request_start, 6)
            return ScientificAgentResult(
                status=AnswerStatus.out_of_scope,
                question=question,
                context_gene=context_gene,
                plan=plan_result.plan,
                summary=plan_result.message
                or "The request is outside supported biomedical research scope.",
                limitations=["No evidence acquisition or scientific synthesis was attempted."],
                evidence_gaps=[],
                agent_activity=[
                    "Planning classified the request as outside supported biomedical research scope."
                ],
                metadata={"timings": timings},
            )
        if plan_result.plan is None:
            status = plan_result.status or AnswerStatus.clarification_required
            timings["total_request"] = round(time.perf_counter() - request_start, 6)
            return ScientificAgentResult(
                status=status,
                question=question,
                context_gene=context_gene,
                summary=plan_result.message
                or "The scientific question could not be planned safely.",
                limitations=["No scientific answer was generated from model prior knowledge."],
                evidence_gaps=[plan_result.message or "A validated evidence plan is unavailable."],
                agent_activity=["Planning stopped before evidence retrieval."],
                metadata={"timings": timings},
            )

        stage_start = time.perf_counter()
        plan = self._augment_hd_plan(plan_result.plan)
        mark("gene_resolution", stage_start)
        if (
            len(plan.entities.genes) > MAX_GENES
            or len(plan.evidence_requirements) > MAX_REQUIREMENTS
        ):
            timings["total_request"] = round(time.perf_counter() - request_start, 6)
            return ScientificAgentResult(
                status=AnswerStatus.clarification_required,
                question=question,
                context_gene=context_gene,
                plan=plan,
                summary="The question exceeds the bounded gene or evidence-requirement limit.",
                limitations=["No acquisition was attempted."],
                evidence_gaps=[
                    "Narrow the question to at most 6 genes and 10 evidence requirements."
                ],
                agent_activity=["Validated plan exceeded execution bounds."],
                metadata={"timings": timings},
            )

        try:
            stage_start = time.perf_counter()
            universes = {
                gene: self._resolve_initial_universe(gene, request) for gene in plan.entities.genes
            }
            mark("persisted_evidence_inspection", stage_start)
        except ValueError as exc:
            timings["total_request"] = round(time.perf_counter() - request_start, 6)
            return ScientificAgentResult(
                status=AnswerStatus.clarification_required,
                question=question,
                context_gene=context_gene,
                plan=plan,
                summary=str(exc),
                limitations=["Evidence run ownership could not be validated."],
                evidence_gaps=[str(exc)],
                agent_activity=["Evidence-universe validation failed."],
                metadata={"timings": timings},
            )

        activity = [
            f"Planner method: {plan.planner_method.value}",
            f"Resolved genes: {', '.join(plan.entities.genes)}",
            f"Validated evidence requirements: {len(plan.evidence_requirements)}",
            "Initial requirement-level retrieval completed.",
        ]
        stage_start = time.perf_counter()
        initial_retrievals = self._retrieve_all(question=question, plan=plan, universes=universes)
        mark("chroma_retrieval_initial", stage_start)
        stage_start = time.perf_counter()
        initial_assessments = [
            self._assessment(
                item,
                acquisition_supported=bool(self._acquisition_specs(item.requirement, plan)),
            )
            for item in initial_retrievals
        ]
        mark("capability_gap_calculation", stage_start)

        tool_activity: list[ToolActivity] = []
        if not request.refresh_if_available:
            stage_start = time.perf_counter()
            initial_assessment_map = {
                (item.requirement_id, item.gene_symbol): item
                for item in initial_assessments
            }
            for requirement in sorted(
                plan.evidence_requirements, key=lambda item: (not item.required, item.id)
            ):
                for gene in requirement.genes:
                    assessment = initial_assessment_map[(requirement.id, gene)]
                    if assessment.status is RequirementStatus.sufficient:
                        continue
                    universe = universes[gene]
                    reusable_runs = self._find_reusable_runs(
                        gene=gene,
                        requirement=requirement,
                        excluded_run_ids=set(universe.dossier_run_ids),
                        plan=plan,
                    )
                    for run_id, reused_capabilities in reusable_runs:
                        if run_id in universe.reused_tool_run_ids:
                            continue
                        universe.reused_tool_run_ids.append(run_id)
                        self._refresh_universe(universe)
                        tool_activity.append(
                            ToolActivity(
                                gene_symbol=gene,
                                capability_ids=reused_capabilities,
                                executor_kind="reused_acquisition",
                                status="completed",
                                dossier_run_id=run_id,
                                public_run_ref=public_run_reference(run_id),
                                evidence_records_persisted=len(
                                    self._load_run_records(gene, run_id)
                                ),
                                execution_succeeded=True,
                                reused=True,
                            )
                        )
            mark("reuse_selection", stage_start)
            if tool_activity:
                stage_start = time.perf_counter()
                initial_retrievals = self._retrieve_all(
                    question=question, plan=plan, universes=universes
                )
                mark("chroma_retrieval_after_reuse", stage_start)
                stage_start = time.perf_counter()
                initial_assessments = [
                    self._assessment(
                        item,
                        acquisition_supported=bool(self._acquisition_specs(item.requirement, plan)),
                    )
                    for item in initial_retrievals
                ]
                mark("capability_gap_calculation_after_reuse", stage_start)
        readiness_rows = self._readiness_rows(
            retrievals=initial_retrievals,
            assessments=initial_assessments,
            universes=universes,
            plan=plan,
        )
        readiness_by_key = {(row.requirement_id, row.gene): row for row in readiness_rows}
        actionable: list[tuple[EvidenceRequirement, str]] = []
        for requirement in sorted(
            plan.evidence_requirements, key=lambda item: (not item.required, item.id)
        ):
            for gene in requirement.genes:
                row = readiness_by_key[(requirement.id, gene)]
                if (
                    row.acquisition_needed
                    or (
                        request.refresh_if_available
                        and (
                            requirement.required
                            or request.research_mode is ResearchMode.deep_research
                        )
                        and row.acquisition_possible
                    )
                    or (
                        request.research_mode is ResearchMode.deep_research
                        and bool(self._acquisition_specs(requirement, plan))
                    )
                ):
                    actionable.append((requirement, gene))

        scheduled: list[tuple[EvidenceRequirement, str, CapabilityId]] = []
        distinct_capabilities: list[CapabilityId] = []
        scheduled_keys: set[tuple[str, str, str, tuple[str, ...]]] = set()
        for requirement, gene in actionable:
            universe = universes[gene]
            if not request.refresh_if_available:
                reusable = self._find_reusable_run(
                    gene=gene,
                    requirement=requirement,
                    excluded_run_ids=set(universe.dossier_run_ids),
                    plan=plan,
                )
                if reusable:
                    run_id, reused_capabilities = reusable
                    universe.reused_tool_run_ids.append(run_id)
                    self._refresh_universe(universe)
                    tool_activity.append(
                        ToolActivity(
                            gene_symbol=gene,
                            capability_ids=reused_capabilities,
                            executor_kind="reused_acquisition",
                            status="completed",
                            dossier_run_id=run_id,
                            public_run_ref=public_run_reference(run_id),
                            evidence_records_persisted=len(self._load_run_records(gene, run_id)),
                            execution_succeeded=True,
                            reused=True,
                        )
                    )
                    continue
            if not acquisition_allowed:
                continue
            candidates = self._acquisition_specs(requirement, plan)
            if not candidates:
                continue
            capability = candidates[0].capability_id
            spec = CAPABILITY_REGISTRY[capability]
            source_key = tuple(spec.sources or spec.section_keys)
            scheduling_key = (gene, requirement.id, capability.value, source_key)
            if scheduling_key in scheduled_keys:
                continue
            scheduled_keys.add(scheduling_key)
            if capability not in distinct_capabilities:
                if len(distinct_capabilities) >= MAX_DISTINCT_ACQUISITIONS:
                    continue
                distinct_capabilities.append(capability)
            scheduled.append((requirement, gene, capability))

        grouped: dict[tuple[str, str], AcquisitionGroup] = {}
        for requirement, gene, capability in scheduled:
            spec = CAPABILITY_REGISTRY[capability]
            key = (gene, spec.executor_kind)
            group = grouped.setdefault(key, AcquisitionGroup(gene, spec.executor_kind))
            if capability not in group.capabilities:
                group.capabilities.append(capability)
            if requirement not in group.requirements:
                group.requirements.append(requirement)

        source_errors = validate_source_capabilities()
        if source_errors:
            activity.extend(f"Source capability disabled: {error}" for error in source_errors)
            grouped = {
                key: group
                for key, group in grouped.items()
                if group.executor_kind != "source_workflow"
            }

        execution_counts: dict[str, int] = defaultdict(int)
        ordered_groups = sorted(
            grouped.values(),
            key=lambda group: (group.executor_kind != "section_bundle", group.gene),
        )
        acquisition_started = time.perf_counter()
        for group in ordered_groups:
            if len([item for item in tool_activity if not item.reused]) >= MAX_TOTAL_EXECUTIONS:
                break
            if execution_counts[group.gene] >= MAX_EXECUTIONS_PER_GENE:
                continue
            universe = universes[group.gene]
            bootstrap = not universe.dossier_run_ids
            stage_start = time.perf_counter()
            result = self._execute_group(group, bootstrap_identity=bootstrap)
            mark(
                f"acquisition.{group.gene}.{group.executor_kind}.{','.join(cap.value for cap in group.capabilities)}",
                stage_start,
            )
            tool_activity.append(result)
            execution_counts[group.gene] += 1
            if result.dossier_run_id and result.execution_succeeded:
                universe.created_tool_run_ids.append(result.dossier_run_id)
                self._refresh_universe(universe)
        mark("acquisition_total", acquisition_started)

        if tool_activity:
            reused_count = len([item for item in tool_activity if item.reused])
            invoked_count = len(tool_activity) - reused_count
            activity.append(
                f"Acquisition stage: {invoked_count} execution(s), {reused_count} reused run(s)."
            )
        else:
            activity.append("Acquisition stage: no execution or historical overlay was needed.")

        stage_start = time.perf_counter()
        final_retrievals = self._retrieve_all(question=question, plan=plan, universes=universes)
        self._update_tool_scientific_counts(tool_activity, final_retrievals)
        mark("chroma_retrieval_final", stage_start)
        stage_start = time.perf_counter()
        assessments = [
            self._assessment(
                item,
                acquisition_supported=bool(self._acquisition_specs(item.requirement, plan)),
            )
            for item in final_retrievals
        ]
        required_assessments = [item for item in assessments if item.required]
        all_required_sufficient = bool(required_assessments) and all(
            item.status is RequirementStatus.sufficient for item in required_assessments
        )
        status = (
            AnswerStatus.answered if all_required_sufficient else AnswerStatus.insufficient_evidence
        )
        selected_records = self._select_final_records(final_retrievals)
        retrieval_method, embedding_backend = self._overall_retrieval_metadata(final_retrievals)

        gap_strings = [
            f"{item.gene_symbol}: {item.evidence_need.value} is {item.status.value} ({item.qualifying_count}/{item.minimum_support})."
            for item in assessments
            if item.status is not RequirementStatus.sufficient
        ]
        structured_gaps = self._structured_gaps(assessments)
        recommendations = self._recommendations(structured_gaps, selected_records)
        citation_registry = self._citation_registry(selected_records, final_retrievals)
        evidence_categories = self._evidence_categories(selected_records, assessments)
        dimensions: list[str] = []
        matrix = []
        comparison_decision = None
        is_comparison = plan.answer_mode.value == "comparison" or plan.intent.value == "comparison"
        if plan.analysis_lens == "hd_modifier_relevance" or is_comparison:
            matrix_requirements = (
                plan.evidence_requirements[: len(hd_modifier_requirements(plan.entities.genes))]
                if plan.analysis_lens == "hd_modifier_relevance"
                else plan.evidence_requirements
            )
            comparison_records = self._comparison_records(
                final_retrievals,
                matrix_requirements,
            )
            dimensions, matrix = build_hd_modifier_matrix(
                genes=plan.entities.genes,
                requirements=matrix_requirements,
                assessments=assessments,
                records=comparison_records,
                plan=plan,
                ordinal_by_id={
                    record.id: index for index, record in enumerate(selected_records, 1)
                },
            )
            comparison_decision = build_comparison_decision(
                plan=plan,
                matrix=matrix,
                assessments=assessments,
            )
            mark("comparison_matrix", stage_start)
        source_attempts = self._source_attempts(universes)
        retrieval_timestamps = sorted(
            {timestamp for block in evidence_categories for timestamp in block.retrieval_timestamps}
        )
        limitations = [
            "Evidence was selected only from the exact per-gene evidence universes reported in this response.",
            "The LLM was not used as a scientific source of truth.",
            "Directional conflict could not be assessed unless structured direction/effect metadata was present.",
        ]
        failures: list[ScientificFailure] = []
        mark("deterministic_assessment_gap_matrix", stage_start)
        if structured_gaps:
            limitations.append(
                "Limited, missing, and unsupported optional evidence remains disclosed in evidenceGaps."
            )
        if not selected_records:
            summary = (
                "Insufficient provenance-backed evidence is available in the selected universes."
            )
            generation_method = "abstain"
            grounding_validation_issues: list[dict[str, Any]] = []
            grounding_summary: dict[str, Any] = {
                "requestedSlotCount": 0,
                "acceptedSlotCount": 0,
                "fallbackSlotCount": 0,
                "diagnosticCounts": {},
            }
        else:
            stage_start = time.perf_counter()
            synthesis = try_grounded_synthesis(
                question=question,
                status=status,
                plan=plan,
                records=selected_records,
                assessments=assessments,
                gaps=structured_gaps,
                recommendations=recommendations,
                comparison_matrix=matrix,
                citation_registry=citation_registry,
                settings=self.settings,
            )
            summary = synthesis.summary
            generation_method = synthesis.generation_method
            recommendations = synthesis.recommendations
            mark("grounded_answer_generation_validation", stage_start)
            grounding_validation_issues = [
                issue.model_dump(mode="json") for issue in synthesis.validation_issues
            ]
            grounding_summary = {
                "requestedSlotCount": synthesis.requested_slot_count,
                "acceptedSlotCount": synthesis.accepted_slot_count,
                "fallbackSlotCount": synthesis.fallback_slot_count,
                "diagnosticCounts": synthesis.diagnostic_counts,
            }
            if synthesis.failure_type:
                failures.append(
                    ScientificFailure(
                        failure_type=synthesis.failure_type,
                        message=synthesis.failure_message or "Grounded answer provider failed.",
                    )
                )

        evidence_items, contextual_evidence = self._public_evidence_items(final_retrievals)
        request_records = {
            record.id: record
            for retrieval in final_retrievals
            for record in [
                *[item for group in retrieval.canonical_groups for item in group.backing_records],
                *[item for item, _eligibility in retrieval.contextual],
                *[item for item, _eligibility in retrieval.excluded],
            ]
        }
        private_identifiers = {
            str(value)
            for record in request_records.values()
            for value in (
                record.id,
                record.dossier_run_id,
                record.api_run_id,
                record.raw_artifact_id,
                record.raw_response_pointer,
            )
            if value
        }
        activity_summary = self._activity_summary(
            plan=plan,
            retrievals=final_retrievals,
            activities=tool_activity,
            readiness_rows=readiness_rows,
        )
        answer_sections = self._answer_sections(
            status=status,
            summary=summary,
            categories=evidence_categories,
            comparison_decision=comparison_decision,
        )

        activity.extend(
            [
                "Final requirement-level retrieval completed.",
                f"Required evidence status: {status.value}",
                f"Generation method: {generation_method}",
            ]
        )
        timings["total_request"] = round(time.perf_counter() - request_start, 6)
        return ScientificAgentResult(
            status=status,
            question=question,
            context_gene=context_gene,
            plan=plan,
            evidence_universes=universes,
            assessments=assessments,
            selected_records=selected_records,
            private_identifiers=private_identifiers,
            summary=summary,
            answer_sections=answer_sections,
            retrieval_method=retrieval_method,
            generation_method=generation_method,
            embedding_backend=embedding_backend,
            limitations=limitations,
            evidence_gaps=gap_strings,
            tool_activity=tool_activity,
            agent_activity=activity,
            comparison_dimensions=dimensions,
            comparison_matrix=matrix,
            comparison_decision=comparison_decision,
            evidence_categories=evidence_categories,
            evidence_items=evidence_items,
            contextual_evidence=contextual_evidence,
            structured_gaps=structured_gaps,
            recommendations=recommendations,
            citation_registry=citation_registry,
            source_attempts=source_attempts,
            retrieval_timestamps=retrieval_timestamps,
            failures=failures,
            activity_summary=activity_summary,
            cost_summary=CostSummary(
                estimated_model_cost_usd=None,
                external_tool_cost_usd=(
                    None if any(not item.reused for item in tool_activity) else 0.0
                ),
                actual_billed_cost_usd=None,
                cost_basis=[
                    "Model cost is unavailable unless provider-reported token usage is captured.",
                    "External public-data tools are reported separately from model usage.",
                    "Actual billed cost requires authoritative provider billing data.",
                ],
                provider_reported_usage={},
            ),
            metadata={
                "timings": timings,
                "readiness": [asdict(row) for row in readiness_rows],
                "groundingValidationIssues": grounding_validation_issues,
                "grounding": grounding_summary,
            },
        )


def summarize_agent_result_runs(result: ScientificAgentResult) -> dict[str, Any]:
    """Summarize acquisition/reuse provenance with public non-reversible handles."""
    created = _unique(
        [
            run_id
            for universe in result.evidence_universes.values()
            for run_id in universe.created_tool_run_ids
        ]
    )
    reused = _unique(
        [
            run_id
            for universe in result.evidence_universes.values()
            for run_id in universe.reused_tool_run_ids
        ]
    )
    tool_created = _unique(
        [
            item.dossier_run_id or ""
            for item in result.tool_activity
            if item.dossier_run_id and not item.reused
        ]
    )
    tool_reused = _unique(
        [
            item.dossier_run_id or ""
            for item in result.tool_activity
            if item.dossier_run_id and item.reused
        ]
    )
    return {
        "createdDossierRunRefs": [public_run_reference(run_id) for run_id in created],
        "reusedDossierRunRefs": [public_run_reference(run_id) for run_id in reused],
        "createdToolRunRefs": [
            public_run_reference(run_id) for run_id in (tool_created or created)
        ],
        "reusedToolRunRefs": [public_run_reference(run_id) for run_id in (tool_reused or reused)],
        "sourceAttempts": [item.model_dump(mode="json") for item in result.source_attempts],
        "evidenceRecordCounts": {
            item.public_run_ref or "unknown": item.evidence_records_persisted
            for item in result.tool_activity
            if item.public_run_ref
        },
        "toolActivity": [item.model_dump(mode="json") for item in result.tool_activity],
        "acquisitionStatus": "acquired"
        if tool_created or created
        else "reused"
        if tool_reused or reused
        else "none",
    }


__all__ = [
    "ACQUISITION_MANIFEST_KEY",
    "ACQUISITION_MANIFEST_SCHEMA",
    "ScientificAgentRequest",
    "ScientificAgentService",
    "RequirementReadinessRow",
    "summarize_agent_result_runs",
]
