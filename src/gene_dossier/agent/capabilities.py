"""Server-owned scientific capability registry and evidence matching rules."""

from __future__ import annotations

from dataclasses import dataclass
from importlib.util import find_spec
from typing import Literal

from gene_dossier.models import EvidenceRecord
from gene_dossier.source_registry import get_source

from .models import CapabilityId, EvidenceNeed, EvidenceRequirement


ExecutorKind = Literal["section_bundle", "source_workflow", "retrieval_only"]


@dataclass(frozen=True)
class CapabilitySpec:
    capability_id: CapabilityId
    label: str
    executor_kind: ExecutorKind
    section_keys: tuple[str, ...] = ()
    sources: tuple[str, ...] = ()
    acquisition_enabled: bool = False
    can_refresh: bool = False
    limitations: str = ""


CAPABILITY_REGISTRY: dict[CapabilityId, CapabilitySpec] = {
    CapabilityId.identity_function: CapabilitySpec(
        CapabilityId.identity_function, "Identity and function", "section_bundle", ("1a",), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.orthology_conservation: CapabilitySpec(
        CapabilityId.orthology_conservation, "Orthology and conservation", "section_bundle", ("1b", "1e"), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.structure_domain: CapabilitySpec(
        CapabilityId.structure_domain, "Structure and domains", "section_bundle", ("1c", "1d"), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.expression_context: CapabilitySpec(
        CapabilityId.expression_context, "Expression context", "section_bundle", ("2a", "2b", "2c"), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.brain_expression: CapabilitySpec(
        CapabilityId.brain_expression, "Brain expression", "section_bundle", ("2a", "2b", "2c"), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.experimental_expression: CapabilitySpec(
        CapabilityId.experimental_expression, "Experimental expression", "section_bundle", ("3a",), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.transcriptional_regulation: CapabilitySpec(
        CapabilityId.transcriptional_regulation, "Transcriptional regulation", "section_bundle", ("4a",), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.ppi: CapabilitySpec(
        CapabilityId.ppi, "Protein interactions", "section_bundle", ("5a", "5b"), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.pathway: CapabilitySpec(
        CapabilityId.pathway, "Pathways", "source_workflow", sources=("Reactome",), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.hd_literature: CapabilitySpec(
        CapabilityId.hd_literature, "Huntington disease literature", "source_workflow", sources=("PubMed",), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.disease_association: CapabilitySpec(
        CapabilityId.disease_association, "Disease association", "retrieval_only", limitations="No controlled disease-association acquisition is enabled in this patch."
    ),
    CapabilityId.human_genetic_association: CapabilitySpec(
        CapabilityId.human_genetic_association, "Human genetic association", "retrieval_only", limitations="Generic disease, expression, pathway, and literature records do not establish a human genetic modifier association."
    ),
    CapabilityId.model_organism: CapabilitySpec(
        CapabilityId.model_organism, "Model-organism evidence", "source_workflow", sources=("MouseMine",), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.chemical_perturbation: CapabilitySpec(
        CapabilityId.chemical_perturbation, "Chemical perturbation", "section_bundle", ("6a",), acquisition_enabled=True, can_refresh=True
    ),
    CapabilityId.chemical_tools: CapabilitySpec(
        CapabilityId.chemical_tools, "Chemical tools", "section_bundle", ("7a",), acquisition_enabled=True, can_refresh=True
    ),
}


NEED_CONTRIBUTORS: dict[EvidenceNeed, tuple[CapabilityId, ...]] = {
    EvidenceNeed.identity_function: (CapabilityId.identity_function,),
    EvidenceNeed.orthology_conservation: (CapabilityId.orthology_conservation, CapabilityId.model_organism),
    EvidenceNeed.structure_domain: (CapabilityId.structure_domain,),
    EvidenceNeed.expression_context: (CapabilityId.expression_context, CapabilityId.brain_expression),
    EvidenceNeed.brain_expression: (CapabilityId.brain_expression, CapabilityId.expression_context),
    EvidenceNeed.experimental_evidence: (CapabilityId.experimental_expression, CapabilityId.model_organism, CapabilityId.chemical_perturbation),
    EvidenceNeed.transcriptional_regulation: (CapabilityId.transcriptional_regulation,),
    EvidenceNeed.protein_interaction: (CapabilityId.ppi,),
    EvidenceNeed.pathway_membership: (CapabilityId.pathway, CapabilityId.ppi),
    EvidenceNeed.hd_literature: (CapabilityId.hd_literature,),
    EvidenceNeed.disease_association: (CapabilityId.disease_association, CapabilityId.hd_literature),
    EvidenceNeed.human_genetic_association: (CapabilityId.human_genetic_association,),
    EvidenceNeed.repeat_instability_mechanism: (
        CapabilityId.hd_literature,
        CapabilityId.pathway,
        CapabilityId.model_organism,
        CapabilityId.experimental_expression,
        CapabilityId.identity_function,
    ),
    EvidenceNeed.model_organism: (CapabilityId.model_organism,),
    EvidenceNeed.chemical_perturbation: (CapabilityId.chemical_perturbation,),
    EvidenceNeed.therapeutic_perturbability: (CapabilityId.chemical_tools, CapabilityId.chemical_perturbation),
}


def validated_capability_ids(requirement: EvidenceRequirement) -> list[CapabilityId]:
    """Return planner-requested capabilities restricted to approved contributors."""
    allowed = NEED_CONTRIBUTORS[requirement.evidence_need]
    requested = requirement.capability_ids or list(allowed)
    return [capability for capability in dict.fromkeys(requested) if capability in allowed]


def acquisition_capabilities(requirement: EvidenceRequirement) -> list[CapabilitySpec]:
    return [
        CAPABILITY_REGISTRY[capability]
        for capability in validated_capability_ids(requirement)
        if CAPABILITY_REGISTRY[capability].acquisition_enabled
    ]


def validate_source_capabilities() -> list[str]:
    """Validate the registered client/normalizer contract for enabled source executors."""
    errors: list[str] = []
    for spec in CAPABILITY_REGISTRY.values():
        if not spec.acquisition_enabled or spec.executor_kind != "source_workflow":
            continue
        for source_name in spec.sources:
            source = get_source(source_name)
            if source is None:
                errors.append(f"{source_name}: source is not registered")
                continue
            if not source.client_module:
                errors.append(f"{source_name}: client module is missing")
            elif find_spec(f"gene_dossier.tools.{source.client_module}") is None:
                errors.append(f"{source_name}: client module cannot be imported")
            if not source.normalizer_module:
                errors.append(f"{source_name}: normalizer module is missing")
            elif find_spec(f"gene_dossier.normalize.{source.normalizer_module}") is None:
                errors.append(f"{source_name}: normalizer module cannot be imported")
            from gene_dossier.workflow import CLIENT_DISPATCH

            if source_name not in CLIENT_DISPATCH:
                errors.append(f"{source_name}: client is not wired into the deterministic workflow")
    return errors


def _record_text(record: EvidenceRecord) -> str:
    assertion = str(getattr(record.assertion_type, "value", record.assertion_type)).lower()
    source_type = str(getattr(record.source_type, "value", record.source_type)).lower()
    value = " ".join(f"{key} {val}" for key, val in (record.value or {}).items())
    return " ".join(
        [record.section, record.subsection or "", record.source_name, source_type, assertion, record.fact_type, value, record.display_text]
    ).lower()


def _assertion(record: EvidenceRecord) -> str:
    return str(getattr(record.assertion_type, "value", record.assertion_type)).lower()


def record_matches_capability(record: EvidenceRecord, capability: CapabilityId) -> bool:
    """Return whether a record was produced by the exact capability semantics.

    Requirement matching can intentionally accept evidence from several
    contributors. Acquisition manifests need the narrower inverse: proof that
    the capability being tagged actually supplied the record.
    """
    assertion = _assertion(record)
    source = record.source_name.strip().lower()
    text = _record_text(record)

    if capability is CapabilityId.identity_function:
        return assertion in {"gene_identity", "protein_function"} and source in {
            "ncbi gene",
            "uniprot",
        }
    if capability is CapabilityId.orthology_conservation:
        return source in {"ensembl", "ncbi datasets", "orthodb", "ucsc"} and (
            assertion == "orthology"
            or any(term in text for term in ("ortholog", "conservation", "homolog"))
        )
    if capability is CapabilityId.structure_domain:
        return assertion == "protein_structure" and source in {
            "alphafold",
            "cdd",
            "pdbe",
            "uniprot",
        }
    if capability is CapabilityId.expression_context:
        return assertion in {"expression", "cell_type_expression"}
    if capability is CapabilityId.brain_expression:
        return assertion in {"expression", "cell_type_expression"} and any(
            term in text
            for term in ("brain", "cortex", "stri", "neuron", "glia", "caudate", "putamen")
        )
    if capability is CapabilityId.experimental_expression:
        return assertion == "perturbation" and source in {"geo", "geo profiles"}
    if capability is CapabilityId.transcriptional_regulation:
        return assertion == "transcription_factor_association"
    if capability is CapabilityId.ppi:
        return assertion == "ppi" and source in {"string", "biogrid"}
    if capability is CapabilityId.pathway:
        return assertion == "pathway_membership" and source == "reactome"
    if capability is CapabilityId.hd_literature:
        return source == "pubmed" and record_matches_need(record, EvidenceNeed.hd_literature)
    if capability is CapabilityId.disease_association:
        return record_matches_need(record, EvidenceNeed.disease_association)
    if capability is CapabilityId.human_genetic_association:
        return record_matches_need(record, EvidenceNeed.human_genetic_association)
    if capability is CapabilityId.model_organism:
        return source == "mousemine" and assertion in {"knockout_phenotype", "orthology"}
    if capability is CapabilityId.chemical_perturbation:
        return source == "ctd" and assertion == "chemical_interaction"
    if capability is CapabilityId.chemical_tools:
        return assertion == "chemical_tool"
    return False


def record_matches_need(record: EvidenceRecord, need: EvidenceNeed) -> bool:
    """Metadata-first requirement matching; semantic ranking happens afterwards."""
    text = _record_text(record)
    assertion = _assertion(record)
    source = record.source_name.lower()

    if need is EvidenceNeed.identity_function:
        return assertion in {"gene_identity", "protein_function"}
    if need is EvidenceNeed.orthology_conservation:
        return assertion == "orthology" or any(term in text for term in ("ortholog", "conservation", "homolog"))
    if need is EvidenceNeed.structure_domain:
        return assertion == "protein_structure" or any(term in text for term in ("domain", "structure", "alphafold", "pdbe", "cdd"))
    if need is EvidenceNeed.expression_context:
        return assertion in {"expression", "cell_type_expression"}
    if need is EvidenceNeed.brain_expression:
        return assertion in {"expression", "cell_type_expression"} and any(
            term in text for term in ("brain", "cortex", "stri", "neuron", "glia", "caudate", "putamen")
        )
    if need is EvidenceNeed.experimental_evidence:
        return assertion in {"perturbation", "knockout_phenotype"} or source in {"geo", "geo profiles", "mousemine"}
    if need is EvidenceNeed.transcriptional_regulation:
        return assertion == "transcription_factor_association" or "transcription factor" in text
    if need is EvidenceNeed.protein_interaction:
        return assertion == "ppi" or source in {"string", "biogrid"}
    if need is EvidenceNeed.pathway_membership:
        return assertion == "pathway_membership" or source in {"reactome", "wikipathways"}
    if need is EvidenceNeed.hd_literature:
        return assertion == "literature_summary" and any(
            term in text for term in ("huntington", "huntingtin", "cag", "hd literature")
        )
    if need is EvidenceNeed.disease_association:
        return assertion == "disease_association"
    if need is EvidenceNeed.human_genetic_association:
        if assertion != "variant_association":
            return False
        return any(term in text for term in ("huntington", "cag", "modifier", "age at onset", "gwas", "somatic expansion"))
    if need is EvidenceNeed.repeat_instability_mechanism:
        return any(term in text for term in ("cag", "repeat instability", "repeat expansion", "somatic expansion", "mismatch repair"))
    if need is EvidenceNeed.model_organism:
        return assertion in {"knockout_phenotype", "orthology"} or source == "mousemine"
    if need is EvidenceNeed.chemical_perturbation:
        return assertion in {"chemical_interaction", "perturbation"} and (source == "ctd" or "chemical" in text)
    if need is EvidenceNeed.therapeutic_perturbability:
        return assertion in {"chemical_tool", "chemical_interaction"}
    return False


def qualifying_records(records: list[EvidenceRecord], requirement: EvidenceRequirement) -> list[EvidenceRecord]:
    genes = set(requirement.genes)
    return [
        record
        for record in records
        if record.gene_symbol.upper() in genes and record_matches_need(record, requirement.evidence_need)
    ]
