"""Structured evidence planner with a conservative deterministic fallback."""

from __future__ import annotations

import logging
import re
from copy import deepcopy
from dataclasses import dataclass

from gene_dossier.config import Settings, get_settings
from gene_dossier.source_registry import list_source_names

from .capabilities import CAPABILITY_REGISTRY, NEED_CONTRIBUTORS, validated_capability_ids
from .models import (
    AnswerMode,
    AnswerStatus,
    CapabilityId,
    EvidenceNeed,
    EvidenceRequirement,
    PlannerMethod,
    ScientificEntities,
    ScientificIntent,
    ScientificQuestionPlan,
    ScientificQuestionPlanDraft,
    ScientificQueryPolicy,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PlanResult:
    plan: ScientificQuestionPlan | None
    status: AnswerStatus | None = None
    message: str | None = None


@dataclass(frozen=True)
class EvidenceNeedMetadata:
    label: str
    description: str


EVIDENCE_NEED_METADATA: dict[EvidenceNeed, EvidenceNeedMetadata] = {
    EvidenceNeed.identity_function: EvidenceNeedMetadata(
        label="Identity / Function",
        description="Gene identity, aliases, canonical identifiers, and core molecular function evidence.",
    ),
    EvidenceNeed.orthology_conservation: EvidenceNeedMetadata(
        label="Orthology / Conservation",
        description="Cross-species orthology, conservation, and homolog evidence.",
    ),
    EvidenceNeed.structure_domain: EvidenceNeedMetadata(
        label="Structure / Domains",
        description="Protein structure, domain, motif, and feature evidence.",
    ),
    EvidenceNeed.expression_context: EvidenceNeedMetadata(
        label="Expression Context",
        description="Tissue, cell-type, or biological-context expression evidence.",
    ),
    EvidenceNeed.brain_expression: EvidenceNeedMetadata(
        label="Brain Expression",
        description="Brain, neuronal, glial, or neuroanatomical expression evidence.",
    ),
    EvidenceNeed.experimental_evidence: EvidenceNeedMetadata(
        label="Experimental Evidence",
        description="Perturbation, model, or experimental evidence relevant to the question.",
    ),
    EvidenceNeed.transcriptional_regulation: EvidenceNeedMetadata(
        label="Transcriptional Regulation",
        description="Transcription factor, promoter, enhancer, and gene-regulatory evidence.",
    ),
    EvidenceNeed.protein_interaction: EvidenceNeedMetadata(
        label="Protein Interaction Evidence",
        description="Protein-protein interaction and interaction-partner evidence.",
    ),
    EvidenceNeed.pathway_membership: EvidenceNeedMetadata(
        label="Pathway Membership",
        description="Curated pathway membership and pathway-context evidence.",
    ),
    EvidenceNeed.hd_literature: EvidenceNeedMetadata(
        label="HD-Specific Literature",
        description="Literature evidence with explicit Huntington disease, huntingtin, or CAG-repeat context.",
    ),
    EvidenceNeed.disease_association: EvidenceNeedMetadata(
        label="Disease Association",
        description="Disease association evidence that is not direct human genetic modifier evidence.",
    ),
    EvidenceNeed.human_genetic_association: EvidenceNeedMetadata(
        label="Human Genetic Association",
        description="Direct human genetic association or modifier evidence.",
    ),
    EvidenceNeed.repeat_instability_mechanism: EvidenceNeedMetadata(
        label="Repeat Instability / Mechanistic Evidence",
        description="Evidence about CAG-repeat instability, somatic expansion, or related mechanisms.",
    ),
    EvidenceNeed.model_organism: EvidenceNeedMetadata(
        label="Model-Organism Evidence",
        description="Model-organism orthology, phenotype, or experimental model evidence.",
    ),
    EvidenceNeed.chemical_perturbation: EvidenceNeedMetadata(
        label="Chemical Perturbation Evidence",
        description="Chemical-gene interaction or perturbation evidence.",
    ),
    EvidenceNeed.therapeutic_perturbability: EvidenceNeedMetadata(
        label="Therapeutic / Perturbability Evidence",
        description="Chemical tool, perturbability, or therapeutic-modulation evidence.",
    ),
    EvidenceNeed.safety_tolerability: EvidenceNeedMetadata(
        label="Safety / Tolerability Evidence",
        description="Safety, tolerability, toxicity, or adverse-effect evidence.",
    ),
    EvidenceNeed.clinical_translational: EvidenceNeedMetadata(
        label="Clinical / Translational Evidence",
        description="Clinical development, translational biomarker, or human intervention evidence.",
    ),
}


def _registered_source_lookup() -> dict[str, str]:
    return {name.casefold(): name for name in list_source_names()}


def _explicit_registered_sources(question: str) -> list[str]:
    """Return only server-registered sources literally named by the scientist."""
    lower = question.casefold()
    matches = [
        canonical
        for normalized, canonical in _registered_source_lookup().items()
        if re.search(rf"(?<![a-z0-9]){re.escape(normalized)}(?![a-z0-9])", lower)
    ]
    return list(dict.fromkeys(matches))


def _validated_query_policy(
    policy: ScientificQueryPolicy,
    *,
    question: str,
) -> ScientificQueryPolicy:
    """Canonicalize planner policy without granting the model new source access."""
    registered = _registered_source_lookup()
    explicit = _explicit_registered_sources(question)
    sources: list[str] = []
    for requested in policy.source_restrictions:
        canonical = registered.get(requested.casefold())
        if canonical is None:
            raise ValueError(f"Unknown source restriction: {requested!r}.")
        if canonical not in explicit:
            raise ValueError(
                f"Source restriction {requested!r} was not explicitly named in the question."
            )
        sources.append(canonical)
    sources = list(dict.fromkeys([*sources, *explicit]))

    lower = question.casefold()
    ranking_requested = bool(
        re.search(
            r"\b(rank|ranking|winner|best|better|stronger|strongest|prefer|priority|prioritize)\b",
            lower,
        )
    )
    causal_required = bool(
        re.search(r"\b(caus(?:e|al|ality)|driv(?:e|er)|responsible for)\b", lower)
    )
    analyze_conflicts = bool(
        re.search(r"\b(conflict(?:ing)?|contradict(?:ory|ion)?|disagree(?:ment)?)\b", lower)
    )
    provenance_focus = bool(
        re.search(r"\b(source|provenance|where did|retriev(?:ed|al)|trace)\b", lower)
    )
    species_scope = policy.species_scope
    if re.search(r"\b(human|patient|clinical|people|participant)\b", lower):
        species_scope = "human"
    elif re.search(r"\b(mouse|mice|murine|model organism|drosophila|zebrafish)\b", lower):
        species_scope = "model_organism"

    criteria = list(policy.comparison_criteria)
    criterion_terms: tuple[tuple[tuple[str, ...], EvidenceNeed], ...] = (
        (("safety", "tolerability", "toxicity", "adverse"), EvidenceNeed.safety_tolerability),
        (("clinical", "translational", "biomarker"), EvidenceNeed.clinical_translational),
        (
            ("therapeutic", "tractability", "perturbability", "drug"),
            EvidenceNeed.therapeutic_perturbability,
        ),
        (("human genetic", "genetic modifier", "gwas"), EvidenceNeed.human_genetic_association),
        (
            ("repeat instability", "somatic cag", "somatic expansion"),
            EvidenceNeed.repeat_instability_mechanism,
        ),
        (("mechanism", "pathway"), EvidenceNeed.pathway_membership),
    )
    for terms, need in criterion_terms:
        if any(term in lower for term in terms):
            criteria.append(need)

    return ScientificQueryPolicy(
        source_restrictions=sources,
        species_scope=species_scope,
        provenance_focus=policy.provenance_focus or provenance_focus,
        analyze_conflicts=policy.analyze_conflicts or analyze_conflicts,
        # Only the scientist's wording may escalate causality. Model-set true values
        # otherwise reject all literature and empty the stored-evidence universe.
        causal_evidence_required=bool(causal_required),
        ranking_requested=policy.ranking_requested or ranking_requested,
        comparison_criteria=list(dict.fromkeys(criteria)),
    )


def _gene_is_explicit(question: str, gene: str) -> bool:
    return bool(
        re.search(rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])", question, re.IGNORECASE)
    )


_BIOMEDICAL_ACRONYM_STOPLIST = {
    "API",
    "CAG",
    "CNS",
    "CTD",
    "DNA",
    "GENE",
    "GENES",
    "GEO",
    "GO",
    "GWAS",
    "HD",
    "HUMAN",
    "KEGG",
    "LLM",
    "MOUSE",
    "MRNA",
    "NCBI",
    "PPI",
    "PROTEIN",
    "PROTEINS",
    "RAG",
    "RNA",
    "SNP",
    "SNPS",
    "TF",
    "WHAT",
    "WHEN",
    "WHERE",
    "WHICH",
    "WHY",
}


def _possible_unresolved_gene_tokens(
    question: str,
    *,
    resolved_genes: list[str],
    context_gene: str | None,
) -> list[str]:
    """Detect unresolved gene-like tokens only as a conservative ambiguity signal."""
    excluded = {gene.upper() for gene in resolved_genes}
    context = (context_gene or "").strip().upper()
    if context and _gene_is_explicit(question, context):
        excluded.add(context)
    uppercase_tokens = re.findall(
        r"(?<![A-Za-z0-9])([A-Z][A-Z0-9-]{1,11})(?![A-Za-z0-9])",
        question,
    )
    target_patterns = (
        r"\binteract(?:s|ed|ing)?\s+with\s+([A-Za-z][A-Za-z0-9-]{1,11})\b",
        r"\binteraction(?:s)?\s+(?:partners?\s+)?(?:of|for|with)\s+([A-Za-z][A-Za-z0-9-]{1,11})\b",
        r"\b(?:domains?|structure|expression|function)\s+(?:of|for|in)\s+([A-Za-z][A-Za-z0-9-]{1,11})\b",
    )
    target_tokens = [
        match.group(1).upper()
        for pattern in target_patterns
        for match in re.finditer(pattern, question, re.IGNORECASE)
    ]
    tokens = [token.upper() for token in uppercase_tokens] + target_tokens
    return list(
        dict.fromkeys(
            token
            for token in tokens
            if token not in excluded and token not in _BIOMEDICAL_ACRONYM_STOPLIST
        )
    )


def _validate_plan(
    draft: ScientificQuestionPlanDraft,
    *,
    question: str,
    context_gene: str | None,
) -> ScientificQuestionPlan:
    genes = list(draft.entities.genes)
    context = (context_gene or "").strip().upper() or None
    if not genes and context and draft.intent is not ScientificIntent.out_of_scope:
        genes = [context]
    if len(genes) > 6:
        raise ValueError("A scientific question may resolve at most 6 genes.")
    if len(genes) > 1 and any(not _gene_is_explicit(question, gene) for gene in genes):
        raise ValueError("Every multi-gene entity must appear explicitly in the question.")
    if len(genes) == 1 and not _gene_is_explicit(question, genes[0]) and genes[0] != context:
        raise ValueError(
            "The planned gene is neither explicit in the question nor the context gene."
        )
    if _possible_unresolved_gene_tokens(
        question,
        resolved_genes=genes,
        context_gene=context,
    ):
        raise ValueError(
            "The question appears to name a target that conflicts with the structured plan."
        )

    requirements: list[EvidenceRequirement] = []
    for requirement in draft.evidence_requirements:
        if any(gene not in genes for gene in requirement.genes):
            raise ValueError(f"Requirement {requirement.id!r} references an unresolved gene.")
        capabilities = validated_capability_ids(requirement)
        if not capabilities:
            # The planner identifies scientific evidence needs; the server owns
            # the validated mechanisms that may contribute to those needs.
            capabilities = list(NEED_CONTRIBUTORS[requirement.evidence_need])
        metadata = EVIDENCE_NEED_METADATA[requirement.evidence_need]
        requirements.append(
            requirement.model_copy(
                update={
                    "capability_ids": capabilities,
                    "label": metadata.label,
                    "description": metadata.description,
                }
            )
        )

    if draft.intent is not ScientificIntent.out_of_scope and not genes:
        raise ValueError("A supported biomedical plan requires a resolved gene.")
    if draft.intent is not ScientificIntent.out_of_scope and not requirements:
        raise ValueError("A supported biomedical plan requires at least one evidence requirement.")
    if draft.intent is not ScientificIntent.out_of_scope and not any(
        item.required for item in requirements
    ):
        raise ValueError(
            "A supported biomedical plan requires at least one required evidence requirement."
        )
    if draft.analysis_lens not in {"general", "hd_modifier_relevance"}:
        raise ValueError("Unknown analysis lens.")
    if draft.analysis_lens == "hd_modifier_relevance" and len(genes) < 2:
        raise ValueError("The HD modifier comparison lens requires multiple explicit genes.")

    primary = (draft.primary_gene or "").strip().upper() or None
    if primary and primary not in genes:
        raise ValueError("primary_gene must be one of the resolved genes.")
    diseases = list(draft.entities.diseases)
    if re.search(r"\b(huntington(?:'s)?(?: disease)?|hd)\b", question, re.IGNORECASE):
        diseases.append("Huntington disease")
    entities = draft.entities.model_copy(
        update={"genes": genes, "diseases": list(dict.fromkeys(diseases))}
    )
    query_policy = _validated_query_policy(draft.query_policy, question=question)
    return ScientificQuestionPlan(
        **draft.model_dump(
            exclude={
                "entities",
                "primary_gene",
                "evidence_requirements",
                "query_policy",
                "requires_multi_gene",
            }
        ),
        entities=entities,
        primary_gene=primary or (genes[0] if len(genes) == 1 else None),
        evidence_requirements=requirements,
        query_policy=query_policy,
        requires_multi_gene=len(genes) > 1,
        planner_method=PlannerMethod.llm_structured,
    )


def _planner_prompt(question: str, context_gene: str | None) -> str:
    capability_lines = "\n".join(
        f"- {capability.value}: {CAPABILITY_REGISTRY[capability].label}"
        for capability in CapabilityId
    )
    need_lines = "\n".join(
        f"- {need.value}: contributors={','.join(cap.value for cap in contributors)}"
        for need, contributors in NEED_CONTRIBUTORS.items()
    )
    return f"""You are an evidence-needs planner for a provenance-controlled biomedical research system.
Do not answer the scientific question. Return only the strict structured plan.
Identify explicit biomedical entities and what scientific evidence is needed.
Capability IDs are abstract evidence mechanisms. Never output function names, section keys, APIs, URLs,
endpoints, query parameters, or sources outside the supplied vocabulary.
Populate source_restrictions only when the user explicitly names a registered source. Never invent a source.
Represent explicit species, provenance, conflict, causality, ranking, and comparison-criterion requests in
query_policy. A request for the strongest or best evidence is a ranking request; it does not authorize a winner.
Do not invent a gene symbol. For a multi-gene plan every gene must occur literally in the question.
For multi-gene comparisons set primary_gene to null.
The context gene is only a hint when the question contains no explicit gene.
Every evidence requirement must explicitly set required true or false. Supporting requirements may be optional.
Use analysis_lens=hd_modifier_relevance only for an explicit HD modifier comparison; otherwise use general.
Keep disease_association separate from human_genetic_association. Generic disease, expression, pathway,
PPI, or literature evidence is not direct human genetic modifier evidence.
Safety/tolerability and clinical/translational evidence are distinct scientific needs and may be unsupported.
For HD modifier questions, request evidence needs such as hd_literature, human_genetic_association,
repeat_instability_mechanism, experimental_evidence, expression_context, pathway_membership, and
therapeutic_perturbability only when the question asks for them. Do not rank targets or declare a winner.
For therapeutic or next-experiment questions, plan evidence needs and gaps; do not propose experiments as evidence.

Supported capability IDs:
{capability_lines}

Supported scientific evidence needs and possible contributors:
{need_lines}

Context gene: {context_gene or "none"}
Registered sources that may be used only when literally named: {", ".join(list_source_names())}
Question: {question}
"""


def build_gemini_planner_schema() -> dict:
    """Return a Gemini-compatible wire schema for the canonical planner draft.

    Gemini rejects ``maxItems`` on the top-level array of nested requirement
    objects, while Pydantic still enforces ``ScientificQuestionPlanDraft`` after
    the wire response is parsed.
    """
    schema = deepcopy(ScientificQuestionPlanDraft.model_json_schema())
    evidence_requirements = schema.get("properties", {}).get("evidence_requirements")
    if isinstance(evidence_requirements, dict):
        evidence_requirements.pop("maxItems", None)
    return schema


def _try_structured_llm_plan(
    question: str,
    context_gene: str | None,
    settings: Settings,
) -> ScientificQuestionPlan | None:
    if not settings.has_llm():
        return None
    try:
        from gene_dossier.synthesis import build_chat_model_candidates
    except Exception as exc:  # noqa: BLE001
        logger.warning("scientific planner unavailable: %s", exc)
        return None

    prompt = _planner_prompt(question, context_gene)
    for candidate in build_chat_model_candidates(settings, purpose="planner"):
        try:
            if candidate.provider == "google_gemini":
                structured = candidate.model.with_structured_output(
                    schema=build_gemini_planner_schema(),
                    method="json_schema",
                )
            elif candidate.provider == "openai":
                structured = candidate.model.with_structured_output(
                    ScientificQuestionPlanDraft,
                    method="json_schema",
                    strict=True,
                )
            else:
                structured = candidate.model.with_structured_output(ScientificQuestionPlanDraft)
            raw = structured.invoke(prompt)
            draft = (
                raw
                if isinstance(raw, ScientificQuestionPlanDraft)
                else ScientificQuestionPlanDraft.model_validate(raw)
            )
            return _validate_plan(draft, question=question, context_gene=context_gene)
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "structured scientific planning failed via %s: %s", candidate.provider, exc
            )
    return None


def _single_requirement_plan(
    *,
    question: str,
    gene: str,
    intent: ScientificIntent,
    answer_mode: AnswerMode,
    need: EvidenceNeed,
    capability: CapabilityId,
    label: str,
    query_policy: ScientificQueryPolicy | None = None,
) -> ScientificQuestionPlan:
    requirement = EvidenceRequirement(
        id=need.value,
        label=label,
        description=f"Qualifying {label.lower()} evidence for {gene}.",
        genes=[gene],
        evidence_need=need,
        capability_ids=[capability],
        required=True,
        minimum_support=2,
        rationale="The question directly requests this evidence category.",
    )
    return ScientificQuestionPlan(
        intent=intent,
        entities=ScientificEntities(genes=[gene]),
        primary_gene=gene,
        objective=question.strip(),
        analysis_lens="general",
        answer_mode=answer_mode,
        evidence_requirements=[requirement],
        requires_multi_gene=False,
        query_policy=query_policy or ScientificQueryPolicy(),
        planner_method=PlannerMethod.deterministic_fallback,
    )


def deterministic_fallback_plan(
    question: str,
    context_gene: str | None,
    *,
    known_genes: list[str] | None = None,
) -> PlanResult:
    """Plan only unambiguous, established single-gene question shapes."""
    q = question.strip()
    lower = q.lower()
    context = (context_gene or "").strip().upper()
    explicit_known = [
        gene.strip().upper()
        for gene in (known_genes or [])
        if gene.strip() and _gene_is_explicit(q, gene.strip().upper())
    ]
    explicit_known = list(dict.fromkeys(explicit_known))
    if len(explicit_known) > 1:
        return PlanResult(
            None,
            AnswerStatus.clarification_required,
            "Structured planning is required to resolve multiple explicit genes safely.",
        )
    if not q:
        return PlanResult(
            None, AnswerStatus.clarification_required, "Please provide a biomedical question."
        )
    if any(
        term in lower
        for term in ("weather", "sports score", "recipe", "stock price", "write a poem")
    ):
        return PlanResult(
            None,
            AnswerStatus.out_of_scope,
            "The request is outside supported biomedical research scope.",
        )
    if any(term in lower for term in ("compare", "versus", " vs ", "genes")):
        return PlanResult(
            None,
            AnswerStatus.clarification_required,
            "Structured planning is required to resolve a multi-gene question safely.",
        )
    unresolved = _possible_unresolved_gene_tokens(
        q,
        resolved_genes=explicit_known,
        context_gene=context,
    )
    if unresolved:
        return PlanResult(
            None,
            AnswerStatus.clarification_required,
            "The question appears to name a gene that could not be safely resolved without structured planning.",
        )
    gene = explicit_known[0] if explicit_known else context
    if not gene:
        return PlanResult(
            None, AnswerStatus.clarification_required, "Select or name a context gene."
        )
    query_policy = _validated_query_policy(ScientificQueryPolicy(), question=q)

    options: list[tuple[tuple[str, ...], EvidenceNeed, CapabilityId, str, AnswerMode]] = [
        (
            ("interact", "interaction", "ppi", "protein partner"),
            EvidenceNeed.protein_interaction,
            CapabilityId.ppi,
            "Protein interaction",
            AnswerMode.fact,
        ),
        (
            ("domain", "structure", "alphafold"),
            EvidenceNeed.structure_domain,
            CapabilityId.structure_domain,
            "Structure and domain",
            AnswerMode.fact,
        ),
        (
            ("brain expression", "expressed in the brain", "brain region"),
            EvidenceNeed.brain_expression,
            CapabilityId.brain_expression,
            "Brain expression",
            AnswerMode.fact,
        ),
        (
            ("expression", "expressed", "tissue", "cell type"),
            EvidenceNeed.expression_context,
            CapabilityId.expression_context,
            "Expression",
            AnswerMode.fact,
        ),
        (
            ("transcription factor", " tf "),
            EvidenceNeed.transcriptional_regulation,
            CapabilityId.transcriptional_regulation,
            "Transcriptional regulation",
            AnswerMode.fact,
        ),
        (
            ("geo", "differential expression", "knockdown"),
            EvidenceNeed.experimental_evidence,
            CapabilityId.experimental_expression,
            "Experimental expression",
            AnswerMode.synthesis,
        ),
        (
            ("chemical perturbation", "chemical-gene", "ctd"),
            EvidenceNeed.chemical_perturbation,
            CapabilityId.chemical_perturbation,
            "Chemical perturbation",
            AnswerMode.synthesis,
        ),
        (
            ("pharmacolog", "chemical tool", "compound", "inhibitor", "drug"),
            EvidenceNeed.therapeutic_perturbability,
            CapabilityId.chemical_tools,
            "Therapeutic perturbability",
            AnswerMode.synthesis,
        ),
        (
            ("identity", "entrez", "uniprot", "alias", "function"),
            EvidenceNeed.identity_function,
            CapabilityId.identity_function,
            "Identity and function",
            AnswerMode.fact,
        ),
    ]
    for terms, need, capability, label, mode in options:
        if any(term in lower for term in terms):
            return PlanResult(
                _single_requirement_plan(
                    question=q,
                    gene=gene,
                    intent=ScientificIntent.single_gene_question,
                    answer_mode=mode,
                    need=need,
                    capability=capability,
                    label=label,
                    query_policy=query_policy,
                )
            )
    return PlanResult(
        None,
        AnswerStatus.clarification_required,
        "The question is biomedical, but structured planning is required to identify its evidence needs safely.",
    )


def plan_scientific_question(
    question: str,
    *,
    context_gene: str | None,
    settings: Settings | None = None,
    known_genes: list[str] | None = None,
) -> PlanResult:
    cfg = settings or get_settings()
    plan = _try_structured_llm_plan(question, context_gene, cfg)
    if plan is not None:
        if plan.intent is ScientificIntent.out_of_scope:
            return PlanResult(
                plan,
                AnswerStatus.out_of_scope,
                "The request is outside supported biomedical research scope.",
            )
        return PlanResult(plan)
    return deterministic_fallback_plan(question, context_gene, known_genes=known_genes)


__all__ = [
    "PlanResult",
    "build_gemini_planner_schema",
    "deterministic_fallback_plan",
    "plan_scientific_question",
]
