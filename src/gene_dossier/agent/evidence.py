"""Deterministic scientific evidence identity, eligibility, and public references."""

from __future__ import annotations

import base64
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Iterable

from gene_dossier.models import EvidenceRecord
from gene_dossier.source_ids import is_valid_source_id, parse_source_id, slugify

from .models import (
    EvidenceDesignation,
    EvidenceNeed,
    EvidenceRequirement,
    PublicEvidenceItem,
    ScientificQueryPolicy,
)


PUBLIC_EVIDENCE_REF_VERSION = "ev_v1"
PUBLIC_RUN_REF_VERSION = "run_v1"
_PUBLIC_EVIDENCE_RE = re.compile(r"^ev_v1_[a-z2-7]{32}$")
_PUBLIC_RUN_RE = re.compile(r"^run_v1_[a-z2-7]{32}$")
_SPACE_RE = re.compile(r"\s+")
_LOCAL_PATH_RE = re.compile(
    r"(?<!https:)(?<!http:)(?:/(?:Users|home|private|tmp|var)/[^\s,;]+)",
    re.IGNORECASE,
)
_HD_TERMS = ("huntington", "huntingtin", "cag", "somatic expansion")
_MECHANISM_TERMS = (
    "repeat instability",
    "repeat expansion",
    "somatic expansion",
    "mismatch repair",
    "dna repair",
    "repair mechanism",
)
_DIRECTION_KEYS = ("direction", "effect", "change", "modulation", "regulation", "action")
_EXTERNAL_ID_KEYS = (
    "pmid",
    "doi",
    "accession",
    "st_id",
    "pathway_id",
    "gds_accession",
    "gds_uid",
    "dataset_id",
    "variation_id",
    "canonical_spdi",
    "target_chembl_id",
    "assay_chembl_id",
    "activity_chembl_id",
    "molecule_chembl_id",
    "document_chembl_id",
    "activity_id",
    "assay_id",
    "aid",
    "chemical_id",
    "biogrid_interaction_id",
    "interaction_id",
    "partner_string_id",
    "mim_number",
    "phenotype_mim_number",
    "uid",
)
_TARGET_GENE_KEYS = (
    "target_gene",
    "target_gene_symbol",
    "gene_symbol_reported",
    "perturbed_gene",
    "query_gene",
    "official_symbol",
)


@dataclass(frozen=True)
class EvidenceEligibility:
    eligible: bool
    designation: EvidenceDesignation
    reason_code: str
    evidence_system: str
    claim_type: str
    directionality_known: bool
    direction: str | None


@dataclass(frozen=True)
class CanonicalEvidenceGroup:
    identity: str
    public_reference: str
    evidence_need: EvidenceNeed
    canonical_record: EvidenceRecord
    backing_records: tuple[EvidenceRecord, ...]
    eligibility: EvidenceEligibility

    @property
    def source_namespace(self) -> str:
        return _source_namespace(self.canonical_record)


@dataclass(frozen=True)
class CanonicalizationResult:
    qualifying: tuple[CanonicalEvidenceGroup, ...]
    contextual: tuple[tuple[EvidenceRecord, EvidenceEligibility], ...]
    excluded: tuple[tuple[EvidenceRecord, EvidenceEligibility], ...]


def _norm(value: Any) -> str:
    return _SPACE_RE.sub(" ", str(value or "").strip().casefold())


def _public_display_text(value: str) -> str:
    return _LOCAL_PATH_RE.sub("[local artifact]", str(value or "")).strip()


def _enum(value: Any) -> str:
    return str(getattr(value, "value", value) or "").strip().casefold()


def _source_namespace(record: EvidenceRecord) -> str:
    return slugify(record.source_name or "unknown", allow_underscore=True) or "unknown"


def _digest_reference(prefix: str, value: str) -> str:
    digest = hashlib.sha256(value.encode("utf-8")).digest()[:20]
    token = base64.b32encode(digest).decode("ascii").rstrip("=").lower()
    return f"{prefix}_{token}"


def public_run_reference(run_id: str, *, gene_symbol: str | None = None) -> str:
    """Return a stable, non-reversible public reference for a private run id."""
    del gene_symbol  # Retained for call-site compatibility; run identity is globally stable.
    material = "\x1f".join(("gene-dossier-public-run-v1", str(run_id).strip()))
    return _digest_reference(PUBLIC_RUN_REF_VERSION, material)


def is_public_evidence_reference(value: str) -> bool:
    return bool(_PUBLIC_EVIDENCE_RE.fullmatch(str(value or "")))


def is_public_run_reference(value: str) -> bool:
    return bool(_PUBLIC_RUN_RE.fullmatch(str(value or "")))


def _clean_external_value(key: str, value: Any) -> str | None:
    if value in (None, "", [], {}):
        return None
    if isinstance(value, (list, dict)):
        return None
    text = str(value).strip()
    if not text:
        return None
    if key == "doi":
        text = re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", text, flags=re.I)
    return text.casefold() if key == "doi" else text


def external_record_identifier(record: EvidenceRecord) -> tuple[str, str] | None:
    """Return a validated public record identifier, never a source registry id."""
    value = record.value or {}
    source = _source_namespace(record)
    source_order: tuple[str, ...]
    if source == "pubmed":
        source_order = ("pmid", "doi")
    elif source == "reactome":
        source_order = ("st_id", "pathway_id", "db_id", "doi")
    elif source == "geo":
        source_order = ("accession", "gds_accession", "gds_uid", "dataset_id")
    elif source == "clinvar":
        source_order = ("accession", "canonical_spdi", "variation_id", "uid")
    elif source == "chembl":
        source_order = (
            "activity_chembl_id",
            "activity_id",
            "assay_chembl_id",
            "molecule_chembl_id",
            "target_chembl_id",
            "document_chembl_id",
        )
    elif source == "string":
        source_order = ("interaction_id", "partner_string_id")
    else:
        source_order = _EXTERNAL_ID_KEYS
    for key in source_order:
        cleaned = _clean_external_value(key, value.get(key))
        if cleaned:
            return key, cleaned
    return None


def public_identifier(record: EvidenceRecord) -> str | None:
    found = external_record_identifier(record)
    if not found:
        return None
    key, value = found
    labels = {
        "pmid": "PMID",
        "doi": "DOI",
        "st_id": "Reactome",
        "gds_uid": "GEO",
        "gds_accession": "GEO",
        "canonical_spdi": "SPDI",
        "target_chembl_id": "ChEMBL",
        "mim_number": "OMIM",
    }
    return f"{labels.get(key, key.replace('_', ' ').title())}: {value}"


def _source_record_key(record: EvidenceRecord) -> tuple[str, str] | None:
    """Return source namespace plus a record-specific key when one is validated."""
    if not is_valid_source_id(record.source_id):
        return None
    parsed = parse_source_id(record.source_id)
    key = parsed["key"].strip()
    if not key or key == parsed["source"] or key == _source_namespace(record):
        return None
    return parsed["assertion"], key


def _stable_public_value(value: Any) -> Any:
    if isinstance(value, dict):
        private_keys = {
            "api_run_id",
            "raw_artifact_id",
            "dossier_run_id",
            "file_path",
            "path",
            "request_headers",
            "authorization",
            "retrieved_at",
            "created_at",
        }
        return {
            str(key): _stable_public_value(item)
            for key, item in sorted(value.items())
            if str(key).casefold() not in private_keys
        }
    if isinstance(value, list):
        return [_stable_public_value(item) for item in value]
    if isinstance(value, str):
        return _norm(value)
    return value


def canonical_evidence_identity(
    record: EvidenceRecord,
    *,
    evidence_need: EvidenceNeed | None = None,
) -> str:
    """Return one deterministic scientific identity independent of database ids."""
    source = _source_namespace(record)
    external = external_record_identifier(record)
    if external:
        kind, value = external
        return f"external:{source}:{kind}:{_norm(value)}"
    record_key = _source_record_key(record)
    if record_key:
        assertion, key = record_key
        return f"record:{source}:{assertion}:{_norm(key)}"
    payload = {
        "gene": _norm(record.official_symbol or record.gene_symbol),
        "source": source,
        "category": evidence_need.value if evidence_need else _primary_need(record).value,
        "assertion": _enum(record.assertion_type),
        "fact_type": _norm(record.fact_type),
        "organism": _norm(record.organism or record.species),
        "taxon_id": record.taxon_id,
        "content": _norm(record.display_text),
        "public_provenance": _stable_public_value(record.value or {}),
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return f"fingerprint:{hashlib.sha256(serialized.encode('utf-8')).hexdigest()}"


def public_evidence_reference(record: EvidenceRecord) -> str:
    identity = canonical_evidence_identity(record)
    return _digest_reference(
        PUBLIC_EVIDENCE_REF_VERSION,
        "\x1f".join(
            (
                "gene-dossier-public-evidence-v1",
                _norm(record.official_symbol or record.gene_symbol),
                identity,
            )
        ),
    )


def _record_text(record: EvidenceRecord, *, include_search: bool = False) -> str:
    value = record.value or {}
    ignored = set() if include_search else {"search_term", "query", "request_url"}
    values = " ".join(
        str(item)
        for key, item in value.items()
        if key not in ignored and item not in (None, "", [], {})
    )
    return _norm(
        " ".join(
            (
                record.section,
                record.subsection or "",
                record.fact_type,
                values,
                record.display_text,
            )
        )
    )


def _structured_categories(record: EvidenceRecord) -> set[str]:
    value = record.value or {}
    raw = value.get(
        "evidence_categories", value.get("evidence_category", value.get("evidence_need"))
    )
    if raw is None:
        return set()
    values = raw if isinstance(raw, list) else [raw]
    return {_enum(item) for item in values if _enum(item)}


def _structured_target_genes(record: EvidenceRecord) -> set[str]:
    value = record.value or {}
    genes: list[Any] = []
    for key in _TARGET_GENE_KEYS:
        item = value.get(key)
        genes.extend(item if isinstance(item, list) else [item])
    return {str(item).strip().upper() for item in genes if item not in (None, "")}


def _gene_aliases(record: EvidenceRecord, gene: str) -> set[str]:
    value = record.value or {}
    raw_aliases: list[Any] = []
    for key in ("aliases", "aliases_used", "gene_aliases"):
        item = value.get(key)
        raw_aliases.extend(item if isinstance(item, list) else [item])
    raw_aliases.extend((gene, record.official_symbol))
    return {str(item).strip() for item in raw_aliases if item not in (None, "")}


def _text_mentions_gene(record: EvidenceRecord, gene: str) -> bool:
    value = record.value or {}
    text = _norm(" ".join(str(value.get(key) or "") for key in ("title", "abstract", "summary")))
    if not text:
        text = _norm(record.display_text)
    return any(
        re.search(rf"(?<![a-z0-9]){re.escape(alias.casefold())}(?![a-z0-9])", text)
        for alias in _gene_aliases(record, gene)
        if alias
    )


def _is_human(record: EvidenceRecord) -> bool | None:
    if record.taxon_id is not None:
        return record.taxon_id == 9606
    organism = _norm(record.organism or record.species)
    if organism:
        if any(term in organism for term in ("homo sapiens", "human")):
            return True
        return False
    if _enum(record.source_type) == "genetic_database" and _source_namespace(record) in {
        "clinvar",
        "open-targets",
        "omim",
    }:
        return True
    return None


def _direction(record: EvidenceRecord) -> str | None:
    value = record.value or {}
    for key, item in value.items():
        if any(term in _norm(key) for term in _DIRECTION_KEYS) and item not in (None, "", [], {}):
            normalized = _norm(item)
            if any(
                term in normalized
                for term in ("increase", "upregulat", "enhance", "gain", "activate")
            ):
                return "increase"
            if any(
                term in normalized
                for term in ("decrease", "downregulat", "reduce", "suppress", "loss", "inhibit")
            ):
                return "decrease"
            return "reported"
    return None


def _primary_need(record: EvidenceRecord) -> EvidenceNeed:
    assertion = _enum(record.assertion_type)
    source = _source_namespace(record)
    if assertion in {"gene_identity", "protein_function"}:
        return EvidenceNeed.identity_function
    if assertion == "orthology":
        return EvidenceNeed.orthology_conservation
    if assertion == "protein_structure":
        return EvidenceNeed.structure_domain
    if assertion in {"expression", "cell_type_expression"}:
        return EvidenceNeed.expression_context
    if assertion == "transcription_factor_association":
        return EvidenceNeed.transcriptional_regulation
    if assertion == "ppi":
        return EvidenceNeed.protein_interaction
    if assertion == "pathway_membership":
        return EvidenceNeed.pathway_membership
    if assertion == "variant_association":
        return EvidenceNeed.human_genetic_association
    if assertion == "disease_association":
        return EvidenceNeed.disease_association
    if assertion == "knockout_phenotype" or source == "mousemine":
        return EvidenceNeed.model_organism
    if assertion == "chemical_tool":
        return EvidenceNeed.therapeutic_perturbability
    if assertion == "chemical_interaction":
        return EvidenceNeed.chemical_perturbation
    if assertion == "perturbation":
        return EvidenceNeed.experimental_evidence
    if assertion == "literature_summary":
        return EvidenceNeed.hd_literature
    return EvidenceNeed.identity_function


def _base_category_match(record: EvidenceRecord, need: EvidenceNeed) -> bool:
    assertion = _enum(record.assertion_type)
    source = _source_namespace(record)
    text = _record_text(record)
    if need is EvidenceNeed.identity_function:
        return assertion in {"gene_identity", "protein_function"}
    if need is EvidenceNeed.orthology_conservation:
        return assertion == "orthology"
    if need is EvidenceNeed.structure_domain:
        return assertion == "protein_structure"
    if need is EvidenceNeed.expression_context:
        return assertion in {"expression", "cell_type_expression"}
    if need is EvidenceNeed.brain_expression:
        return assertion in {"expression", "cell_type_expression"} and any(
            term in text
            for term in ("brain", "cortex", "stri", "neuron", "glia", "caudate", "putamen")
        )
    if need is EvidenceNeed.experimental_evidence:
        return assertion in {"perturbation", "knockout_phenotype"}
    if need is EvidenceNeed.transcriptional_regulation:
        return assertion == "transcription_factor_association"
    if need is EvidenceNeed.protein_interaction:
        return assertion == "ppi" and source in {"string", "biogrid"}
    if need is EvidenceNeed.pathway_membership:
        return assertion == "pathway_membership" and source in {"reactome", "wikipathways"}
    if need is EvidenceNeed.hd_literature:
        return assertion == "literature_summary" and any(term in text for term in _HD_TERMS)
    if need is EvidenceNeed.disease_association:
        return assertion == "disease_association"
    if need is EvidenceNeed.human_genetic_association:
        return assertion == "variant_association" and any(
            term in text for term in ("modifier", "gwas", "age at onset", "somatic expansion")
        )
    if need is EvidenceNeed.repeat_instability_mechanism:
        return any(term in text for term in _MECHANISM_TERMS)
    if need is EvidenceNeed.model_organism:
        return assertion in {"knockout_phenotype", "orthology"} and source == "mousemine"
    if need is EvidenceNeed.chemical_perturbation:
        return assertion == "chemical_interaction" and source == "ctd"
    if need is EvidenceNeed.therapeutic_perturbability:
        return assertion in {"chemical_tool", "chemical_interaction"}
    if need is EvidenceNeed.safety_tolerability:
        return bool(_structured_categories(record) & {need.value})
    if need is EvidenceNeed.clinical_translational:
        return bool(_structured_categories(record) & {need.value})
    return False


def evaluate_evidence(
    record: EvidenceRecord,
    requirement: EvidenceRequirement,
    *,
    gene: str,
    query_policy: ScientificQueryPolicy | None = None,
    disease_contexts: Iterable[str] = (),
) -> EvidenceEligibility:
    """Classify one persisted record without model judgment."""
    policy = query_policy or ScientificQueryPolicy()
    gene = gene.strip().upper()
    assertion = _enum(record.assertion_type)
    source = record.source_name.strip()
    source_type = _enum(record.source_type)
    direction = _direction(record)

    def decision(
        eligible: bool,
        designation: EvidenceDesignation,
        reason: str,
    ) -> EvidenceEligibility:
        return EvidenceEligibility(
            eligible=eligible,
            designation=designation,
            reason_code=reason,
            evidence_system=source_type or "unspecified",
            claim_type=assertion or _norm(record.fact_type) or "unspecified",
            directionality_known=direction is not None,
            direction=direction,
        )

    if record.gene_symbol.strip().upper() != gene:
        return decision(False, EvidenceDesignation.excluded, "off_gene_record")
    structured_targets = _structured_target_genes(record)
    if structured_targets and gene not in structured_targets:
        return decision(False, EvidenceDesignation.excluded, "structured_target_mismatch")
    if policy.source_restrictions and source.casefold() not in {
        item.casefold() for item in policy.source_restrictions
    }:
        return decision(False, EvidenceDesignation.excluded, "source_restriction_mismatch")
    explicit_categories = _structured_categories(record)
    if explicit_categories and requirement.evidence_need.value not in explicit_categories:
        return decision(False, EvidenceDesignation.excluded, "structured_category_mismatch")
    if not _base_category_match(record, requirement.evidence_need):
        return decision(False, EvidenceDesignation.excluded, "category_or_claim_mismatch")
    if (
        requirement.evidence_need
        in {
            EvidenceNeed.repeat_instability_mechanism,
            EvidenceNeed.experimental_evidence,
            EvidenceNeed.therapeutic_perturbability,
        }
        and _primary_need(record) is not requirement.evidence_need
        and requirement.evidence_need.value not in explicit_categories
    ):
        return decision(
            False, EvidenceDesignation.contextual, "cross_category_mapping_not_validated"
        )

    human = _is_human(record)
    if policy.species_scope == "human" and human is False:
        return decision(False, EvidenceDesignation.contextual, "nonhuman_context")
    if policy.species_scope == "model_organism" and human is True:
        return decision(False, EvidenceDesignation.contextual, "human_not_model_context")
    if requirement.evidence_need is EvidenceNeed.human_genetic_association and human is False:
        return decision(False, EvidenceDesignation.contextual, "nonhuman_genetic_evidence")

    source_namespace = _source_namespace(record)
    if source_namespace == "pubmed" and not _text_mentions_gene(record, gene):
        return decision(
            False, EvidenceDesignation.contextual, "target_not_validated_in_article_metadata"
        )
    if assertion in {"perturbation", "knockout_phenotype"} and not _text_mentions_gene(
        record, gene
    ):
        return decision(
            False, EvidenceDesignation.contextual, "target_not_validated_in_experiment_metadata"
        )
    if source_namespace == "ctd" and structured_targets and gene not in structured_targets:
        return decision(False, EvidenceDesignation.excluded, "chemical_target_mismatch")
    if source_namespace == "chembl":
        method = _norm((record.value or {}).get("target_selection_method"))
        if method and method != "matched":
            return decision(
                False, EvidenceDesignation.contextual, "chemical_target_not_safely_matched"
            )
    if requirement.evidence_need in {
        EvidenceNeed.chemical_perturbation,
        EvidenceNeed.therapeutic_perturbability,
    }:
        value = record.value or {}
        evidence_class = _norm(value.get("evidence_class"))
        target_validated = bool(structured_targets and gene in structured_targets)
        target_validated = (
            target_validated or _norm(value.get("target_selection_method")) == "matched"
        )
        target_validated = target_validated or evidence_class in {
            "direct_target_evidence",
            "explicit_chemical_tool",
            "perturbational_chemical",
        }
        if source_namespace == "ctd":
            target_validated = target_validated or gene in structured_targets
        if not target_validated:
            return decision(False, EvidenceDesignation.contextual, "chemical_target_not_validated")

    hd_context = any(
        "huntington" in _norm(disease)
        or bool(re.search(r"(?<![a-z0-9])hd(?![a-z0-9])", _norm(disease)))
        for disease in disease_contexts
    )
    if hd_context and requirement.evidence_need in {
        EvidenceNeed.experimental_evidence,
        EvidenceNeed.human_genetic_association,
        EvidenceNeed.repeat_instability_mechanism,
    }:
        if not any(term in _record_text(record) for term in _HD_TERMS):
            return decision(False, EvidenceDesignation.contextual, "disease_context_not_validated")

    if requirement.evidence_need is EvidenceNeed.human_genetic_association:
        if source_type != "genetic_database":
            return decision(False, EvidenceDesignation.contextual, "not_human_genetic_system")
        if not any(
            key in (record.value or {})
            for key in ("accession", "canonical_spdi", "variant_id", "gwas", "association_role")
        ):
            return decision(
                False, EvidenceDesignation.contextual, "human_genetic_metadata_incomplete"
            )

    if policy.causal_evidence_required:
        if assertion not in {"perturbation", "knockout_phenotype"}:
            return decision(
                False, EvidenceDesignation.contextual, "noncausal_evidence_for_causal_question"
            )
        if direction is None:
            return decision(
                False, EvidenceDesignation.contextual, "causal_direction_not_structured"
            )

    direct_assertions = {
        "gene_identity",
        "protein_function",
        "protein_structure",
        "ppi",
        "pathway_membership",
        "variant_association",
        "perturbation",
        "knockout_phenotype",
        "chemical_interaction",
        "chemical_tool",
    }
    designation = (
        EvidenceDesignation.direct
        if assertion in direct_assertions
        else EvidenceDesignation.supporting
    )
    return decision(True, designation, "qualifying")


def _grade_rank(record: EvidenceRecord) -> int:
    return {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}.get(
        str(getattr(record.evidence_grade, "value", record.evidence_grade)).upper(),
        99,
    )


def _representative_key(
    record: EvidenceRecord,
    eligibility: EvidenceEligibility,
) -> tuple[int, int, int, str, str, str, str]:
    created_at = getattr(record, "created_at", None)
    return (
        0 if eligibility.designation is EvidenceDesignation.direct else 1,
        _grade_rank(record),
        1 if record.manual_review_required else 0,
        record.source_name.casefold(),
        canonical_evidence_identity(record),
        _norm(record.source_id),
        created_at.isoformat() if isinstance(created_at, datetime) else "",
    )


def canonicalize_requirement_evidence(
    records: Iterable[EvidenceRecord],
    requirement: EvidenceRequirement,
    *,
    gene: str,
    query_policy: ScientificQueryPolicy | None = None,
    disease_contexts: Iterable[str] = (),
) -> CanonicalizationResult:
    """Apply eligibility and deduplicate one gene/dimension evidence set."""
    rows = list(records)
    grouped: dict[str, list[tuple[EvidenceRecord, EvidenceEligibility]]] = {}
    contextual: list[tuple[EvidenceRecord, EvidenceEligibility]] = []
    excluded: list[tuple[EvidenceRecord, EvidenceEligibility]] = []
    for record in rows:
        eligibility = evaluate_evidence(
            record,
            requirement,
            gene=gene,
            query_policy=query_policy,
            disease_contexts=disease_contexts,
        )
        if not eligibility.eligible:
            target = (
                contextual
                if eligibility.designation is EvidenceDesignation.contextual
                else excluded
            )
            target.append((record, eligibility))
            continue
        identity = canonical_evidence_identity(record, evidence_need=requirement.evidence_need)
        grouped.setdefault(identity, []).append((record, eligibility))

    qualifying: list[CanonicalEvidenceGroup] = []
    for identity, members in grouped.items():
        members.sort(key=lambda item: _representative_key(item[0], item[1]))
        representative, eligibility = members[0]
        qualifying.append(
            CanonicalEvidenceGroup(
                identity=identity,
                public_reference=public_evidence_reference(representative),
                evidence_need=requirement.evidence_need,
                canonical_record=representative,
                backing_records=tuple(item[0] for item in members),
                eligibility=eligibility,
            )
        )
    qualifying.sort(
        key=lambda item: (
            _representative_key(item.canonical_record, item.eligibility),
            item.identity,
        )
    )
    contextual.sort(key=lambda item: (item[1].reason_code, item[0].source_name, item[0].source_id))
    excluded.sort(key=lambda item: (item[1].reason_code, item[0].source_name, item[0].source_id))
    return CanonicalizationResult(tuple(qualifying), tuple(contextual), tuple(excluded))


def record_title(record: EvidenceRecord) -> str | None:
    value = record.value or {}
    for key in ("title", "display_name", "preferred_name", "variation_name", "chemical_name"):
        if value.get(key):
            return str(value[key]).strip()
    return None


def record_source_url(record: EvidenceRecord) -> str | None:
    value = record.value or {}
    for key in ("url", "source_url", "link", "detail_url", "browser_url"):
        candidate = value.get(key)
        if isinstance(candidate, str) and candidate.startswith(("http://", "https://")):
            return candidate
    found = external_record_identifier(record)
    if not found:
        return None
    key, value = found
    if key == "pmid":
        return f"https://pubmed.ncbi.nlm.nih.gov/{value}/"
    if key == "doi":
        return f"https://doi.org/{value}"
    if key == "st_id":
        return f"https://reactome.org/content/detail/{value}"
    return None


def persisted_timestamp(record: EvidenceRecord) -> str | None:
    created = record.created_at
    if isinstance(created, datetime):
        return created.isoformat()
    return str(created) if created else None


def public_item_from_group(group: CanonicalEvidenceGroup) -> PublicEvidenceItem:
    record = group.canonical_record
    return PublicEvidenceItem(
        public_evidence_ref=group.public_reference,
        gene_symbol=record.gene_symbol.upper(),
        source_name=record.source_name,
        public_identifier=public_identifier(record),
        title=record_title(record),
        source_url=record_source_url(record),
        evidence_need=group.evidence_need,
        designation=group.eligibility.designation,
        display_text=_public_display_text(record.display_text),
        retrieved_at=persisted_timestamp(record),
        backing_record_count=len(group.backing_records),
    )


def public_item_from_exclusion(
    record: EvidenceRecord,
    eligibility: EvidenceEligibility,
    *,
    evidence_need: EvidenceNeed,
) -> PublicEvidenceItem:
    return PublicEvidenceItem(
        public_evidence_ref=public_evidence_reference(record),
        gene_symbol=record.gene_symbol.upper(),
        source_name=record.source_name,
        public_identifier=public_identifier(record),
        title=record_title(record),
        source_url=record_source_url(record),
        evidence_need=evidence_need,
        designation=eligibility.designation,
        display_text=_public_display_text(record.display_text),
        retrieved_at=persisted_timestamp(record),
        exclusion_reason=eligibility.reason_code,
    )


def contains_private_value(payload: Any, private_values: Iterable[str]) -> str | None:
    """Recursively find the first exact known private value in a public payload."""
    values = {str(value) for value in private_values if value}
    if not values:
        return None
    if isinstance(payload, dict):
        for key, value in payload.items():
            found = contains_private_value(key, values) or contains_private_value(value, values)
            if found:
                return found
        return None
    if isinstance(payload, (list, tuple, set)):
        for item in payload:
            found = contains_private_value(item, values)
            if found:
                return found
        return None
    if isinstance(payload, str):
        return next((value for value in values if value in payload), None)
    return None


__all__ = [
    "CanonicalEvidenceGroup",
    "CanonicalizationResult",
    "EvidenceEligibility",
    "canonical_evidence_identity",
    "canonicalize_requirement_evidence",
    "contains_private_value",
    "evaluate_evidence",
    "external_record_identifier",
    "is_public_evidence_reference",
    "is_public_run_reference",
    "persisted_timestamp",
    "public_evidence_reference",
    "public_identifier",
    "public_item_from_exclusion",
    "public_item_from_group",
    "public_run_reference",
    "record_source_url",
    "record_title",
]
