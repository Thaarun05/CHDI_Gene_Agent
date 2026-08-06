"""Bundle-only Section 3a GEO Profiles (brain/neuron perturbation screening).

Owns GEO Profiles discovery, GDS metadata enrichment, bounded chart shortlist
acquisition, and section statuses. Accepted gene pointers are written only after
post-render evaluation in ``section_bundle`` — never from this node.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceCoverageResult,
    SourceStatus,
    SourceType,
)
from gene_dossier.section_1c import _append_evidence
from gene_dossier.section_3a_sources import (
    Section3aPaths,
    accept_gene_report,
    paths_for,
    sha256_bytes,
    write_json_atomic,
)
from gene_dossier.source_ids import make_source_id, slugify
from gene_dossier.tools import geo_profiles as gp
from gene_dossier.workflow import DossierState, WorkflowTransientContext


def _evidence(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    organism: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=make_source_id(
            gp.SOURCE_NAME, gene_symbol, AssertionType.perturbation, key
        ),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_PERTURBATIONS,
        subsection=SUBSECTION_3A,
        source_name=gp.SOURCE_NAME,
        source_type=SourceType.expression_database,
        assertion_type=AssertionType.perturbation,
        fact_type=fact_type,
        organism=organism,
        evidence_grade=EvidenceGrade.B,
        value=value,
        display_text=display_text,
    )

logger = logging.getLogger(__name__)

SECTION_PERTURBATIONS = "Perturbations in GEO that alter the gene"
SUBSECTION_3A = "GEO Profiles search focusing on brain and/or neurons"

STATUS_SUCCESS = "success"
STATUS_SOURCE_UNAVAILABLE = "source_unavailable"
STATUS_PARTIAL = "partial"
STATUS_UNAVAILABLE = "unavailable"
STATUS_NOT_ATTEMPTED = "not_attempted_optional"
STATUS_NO_RELEVANT = "no_relevant_profiles"
STATUS_FAILED = "failed"

GEO_INTRO_TEXT = (
    "The NCBI Gene Expression Omnibus (GEO) Profiles database indexes curated "
    "expression measurements for individual genes across GEO DataSets. Profiles "
    "summarize how a gene's reported expression varies across experimental "
    "conditions within a DataSet; DataSets follow MIAME-oriented curation and "
    "include sample annotations that support screening for brain or neuron "
    "contexts and experimental perturbations."
)

COMPARABILITY_NOTE = (
    "Expression values and scales are DataSet-specific and should not be "
    "compared directly between different GEO DataSets."
)

SCREENING_CAVEAT = (
    "Selection of GEO Profiles for presentation does not establish statistically "
    "significant differential expression, causality, or therapeutic relevance; "
    "charts are descriptive views of curated profile values within each DataSet."
)

SELECTION_POLICY = (
    "Profiles shown below were selected as potentially relevant based on neural "
    "context, perturbation/comparator design cues, metadata completeness, and "
    "(when figures are enabled) validated chart acquisition under a bounded "
    "diversity-aware shortlist."
)


@dataclass(frozen=True)
class Section3aConfig:
    output_root: str | Path | None = None
    force_refresh: bool = False
    max_discovery_profiles: int = gp.DEFAULT_MAX_DISCOVERY
    max_selected_profiles: int = gp.DEFAULT_MAX_SELECTED
    max_chart_candidates: int | None = None
    attempt_figures: bool = True

    def __post_init__(self) -> None:
        if int(self.max_discovery_profiles) < 1:
            raise ValueError("max_discovery_profiles must be >= 1")
        if int(self.max_selected_profiles) < 1:
            raise ValueError("max_selected_profiles must be >= 1")
        if self.max_chart_candidates is not None and int(self.max_chart_candidates) < 1:
            raise ValueError("max_chart_candidates must be >= 1 when set")


def build_intro_text(
    gene_symbol: str,
    *,
    exact_count: int | None,
    neural_count: int | None,
    subset_count: int | None,
) -> str:
    gene = (gene_symbol or "").strip() or "this gene"
    exact = "an unknown number of" if exact_count is None else str(exact_count)
    neural = "an unknown number of" if neural_count is None else str(neural_count)
    subset = "an unknown number of" if subset_count is None else str(subset_count)
    return (
        f"{GEO_INTRO_TEXT} For {gene}, GEO Profiles reports {exact} exact Gene Symbol "
        f"profile records, {neural} neural-context profile records, and {subset} "
        f"neural-context profile records flagged with value subset effect "
        f"(used only as a ranking signal, not as a formal differential-expression "
        f"call). {SELECTION_POLICY} {COMPARABILITY_NOTE}"
    )


def format_reporter_line(profile: dict[str, Any]) -> str:
    parts: list[str] = []
    gds = profile.get("gds_accession") or gp.format_gds_accession(profile.get("gds_uid"))
    if gds:
        parts.append(str(gds))
    gpl = str(profile.get("gpl") or "").strip()
    if gpl:
        if gpl.isdigit():
            gpl = f"GPL{gpl}"
        elif not gpl.upper().startswith("GPL"):
            gpl = f"GPL{gpl}"
        parts.append(gpl)
    idref = profile.get("idref")
    if idref:
        parts.append(str(idref))
    platform = profile.get("platform_technology") or profile.get("gdstype")
    if platform:
        parts.append(str(platform))
    return " · ".join(parts)


def presentation_profile(profile: dict[str, Any]) -> dict[str, Any]:
    gds_meta = dict(profile.get("gds_metadata") or {})
    uid = str(profile.get("profile_uid") or "")
    gds_acc = (
        profile.get("gds_accession")
        or gp.format_gds_accession(profile.get("gds_uid"))
        or ""
    )
    return {
        "profile_uid": uid,
        "profile_url": gp.PROFILE_PAGE_TMPL.format(uid=uid) if uid else None,
        "title": profile.get("title"),
        "taxon": profile.get("taxon"),
        "organism": profile.get("taxon") or gds_meta.get("organism"),
        "genename": profile.get("genename"),
        "gds_uid": gp.normalize_gds_uid(profile.get("gds_uid")),
        "gds_accession": gds_acc,
        "gpl": profile.get("gpl"),
        "idref": profile.get("idref"),
        "reporter_line": format_reporter_line({**profile, "gds_accession": gds_acc}),
        "sample_count": gds_meta.get("sample_count"),
        "pubmed_id": gds_meta.get("pubmed_id") or (gds_meta.get("pubmedids") or [None])[0],
        "gse": gds_meta.get("gse"),
        "subset_effect_flag": bool(profile.get("subset_effect_flag")),
        "final_score": profile.get("final_score"),
        "score_components": dict(profile.get("score_components") or {}),
        "diversity_keys": dict(profile.get("diversity_keys") or {}),
        "graph_status": profile.get("graph_status"),
        "graph_ok": bool(profile.get("graph_ok")),
        "selection_rank": profile.get("selection_rank"),
        "figure_relative_path": profile.get("figure_relative_path"),
        "local_artifact_path": profile.get("local_artifact_path"),
        "figure_sha256": profile.get("figure_sha256"),
        "image_width": profile.get("image_width"),
        "image_height": profile.get("image_height"),
        "acquisition_method": profile.get("acquisition_method"),
        "graph_requested_url": profile.get("graph_requested_url"),
        "graph_final_url": profile.get("graph_final_url"),
        "graph_url_origin": profile.get("graph_url_origin"),
        "validation_checks": dict(profile.get("validation_checks") or {}),
        "graph_error_type": profile.get("graph_error_type") or profile.get("error_type"),
    }


def _candidate_audit_row(profile: dict[str, Any]) -> dict[str, Any]:
    row = presentation_profile(profile)
    row.update(
        {
            "eligibility_status": profile.get("eligibility_status"),
            "rejection_reasons": list(profile.get("rejection_reasons") or []),
            "link_validation_status": profile.get("link_validation_status"),
            "in_chart_shortlist": bool(profile.get("in_chart_shortlist")),
            "selected": bool(profile.get("selected")),
            "base_score": profile.get("base_score"),
        }
    )
    return row


def _overall_status(scientific: str, visual: str) -> str:
    if scientific == STATUS_SOURCE_UNAVAILABLE:
        return STATUS_FAILED
    if scientific == STATUS_NO_RELEVANT:
        return STATUS_PARTIAL if visual in {STATUS_PARTIAL, STATUS_UNAVAILABLE} else scientific
    if visual in {STATUS_SUCCESS, STATUS_NOT_ATTEMPTED}:
        return STATUS_SUCCESS
    if visual == STATUS_PARTIAL:
        return STATUS_PARTIAL
    if visual == STATUS_UNAVAILABLE:
        return STATUS_PARTIAL
    return STATUS_PARTIAL


def evaluate_section_3a_scientific_complete(
    *,
    status: dict[str, Any],
    pdf_render_status: str,
) -> dict[str, Any]:
    rendering = dict(status.get("rendering_status") or status)
    summary = dict(status.get("summary") or {})
    scientific = str(rendering.get("scientific_status") or "")
    selected_count = int(
        summary.get("selected_profile_count")
        or len(summary.get("selected_profiles") or [])
        or 0
    )
    scientific_complete = (
        scientific == STATUS_SUCCESS
        and selected_count >= 1
        and pdf_render_status == STATUS_SUCCESS
    )
    return {
        "scientific_complete": scientific_complete,
        "scientific_status": scientific,
        "selected_profile_count": selected_count,
        "pdf_render_status": pdf_render_status,
    }


def evaluate_section_3a_visual_complete(
    *,
    status: dict[str, Any],
    embedded_figure_count: int,
    selected_count: int,
    pdf_render_status: str,
) -> dict[str, Any]:
    scientific = evaluate_section_3a_scientific_complete(
        status=status,
        pdf_render_status=pdf_render_status,
    )
    rendering = dict(status.get("rendering_status") or status)
    visual = str(rendering.get("visual_status") or "")
    visual_complete = (
        scientific.get("scientific_complete") is True
        and visual == STATUS_SUCCESS
        and selected_count >= 1
        and embedded_figure_count == selected_count
    )
    return {
        **scientific,
        "visual_complete": visual_complete,
        "visual_status": visual,
        "embedded_figure_count": int(embedded_figure_count),
        "selected_count": int(selected_count),
    }


def accept_scientific_complete_gene_report(
    paths: Section3aPaths,
    *,
    gene_symbol: str,
    attempt_dir: Path,
    acceptance: dict[str, Any],
    artifacts: dict[str, Any] | None = None,
) -> Path | None:
    """Pin scientific-complete acceptance; never replace a prior accepted pointer."""
    pointer = paths.accepted_gene_pointer(gene_symbol)
    if pointer.is_file():
        try:
            existing = json.loads(pointer.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            existing = {}
        prior = dict(existing.get("acceptance") or {})
        if prior.get("section_3a_scientific_complete") is True:
            return None
    return accept_gene_report(
        paths,
        gene_symbol=gene_symbol,
        attempt_dir=attempt_dir,
        acceptance=acceptance,
        artifacts=artifacts or {},
    )


def accept_visual_complete_gene_report(
    paths: Section3aPaths,
    *,
    gene_symbol: str,
    attempt_dir: Path,
    acceptance: dict[str, Any],
    artifacts: dict[str, Any] | None = None,
    promote_existing: bool = False,
) -> Path | None:
    """Pin visual-complete acceptance; never downgrade without promote."""
    pointer = paths.accepted_gene_pointer(gene_symbol)
    if pointer.is_file():
        try:
            existing = json.loads(pointer.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            existing = {}
        prior = dict(existing.get("acceptance") or {})
        if prior.get("section_3a_visual_complete") is True and not promote_existing:
            return None
    return accept_gene_report(
        paths,
        gene_symbol=gene_symbol,
        attempt_dir=attempt_dir,
        acceptance=acceptance,
        artifacts=artifacts or {},
    )


def node_generate_section_3a_derived_artifacts(
    state: DossierState,
    *,
    settings: Settings | None = None,
    persist_db: bool = True,
    transient: WorkflowTransientContext | None = None,
    config: Section3aConfig | None = None,
) -> DossierState:
    """Generate Section 3a attempt artifacts and return ``section_3a_status``.

    Does not update accepted gene pointers.
    """
    _ = transient
    cfg = settings or get_settings()
    section_cfg = config or Section3aConfig()
    run_type = state.get("run_type")
    selected_keys = list(state.get("selected_section_keys") or [])
    if run_type != "section_bundle" or "3a" not in selected_keys:
        return state

    gene = str(state.get("gene_symbol") or "").strip()
    if not gene:
        gene_ids = state.get("gene_ids") or {}
        gene = str(gene_ids.get("symbol") or gene_ids.get("gene_symbol") or "").strip()
    if not gene:
        errors = list(state.get("errors") or [])
        errors.append("Section 3a requires a resolved gene_symbol")
        return {**state, "errors": errors}

    collected = gp.collect_section_3a_profiles(
        gene,
        max_discovery_profiles=section_cfg.max_discovery_profiles,
        max_selected_profiles=section_cfg.max_selected_profiles,
        max_chart_candidates=section_cfg.max_chart_candidates,
        attempt_figures=section_cfg.attempt_figures,
        settings=cfg,
    )

    paths = paths_for(section_cfg.output_root or cfg.output_path)
    attempt_dir = paths.new_gene_attempt(
        gene,
        run_id=str(state.get("dossier_run_id") or ""),
    )
    figures_dir = attempt_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    selected_for_present: list[dict[str, Any]] = []
    figure_artifacts: dict[str, Any] = {}
    for profile in list(collected.get("selected_profiles") or []):
        present = presentation_profile(profile)
        image = profile.get("graph_image_bytes")
        if isinstance(image, (bytes, bytearray)) and present.get("graph_ok"):
            uid = str(present["profile_uid"])
            rel = f"figures/{uid}.png"
            out_path = attempt_dir / rel
            content = bytes(image)
            out_path.write_bytes(content)
            digest = sha256_bytes(content)
            present["figure_relative_path"] = rel
            present["local_artifact_path"] = str(out_path)
            present["figure_sha256"] = digest
            present["image_width"] = profile.get("image_width") or present.get("image_width")
            present["image_height"] = profile.get("image_height") or present.get("image_height")
            figure_artifacts[uid] = {
                "relative_path": rel,
                "local_artifact_path": str(out_path),
                "sha256": digest,
                "byte_size": len(content),
            }
            profile["figure_relative_path"] = rel
            profile["local_artifact_path"] = str(out_path)
            profile["figure_sha256"] = digest
        selected_for_present.append(present)

    candidates_audit = [
        _candidate_audit_row(c) for c in list(collected.get("candidates") or [])
    ]
    # Never persist image bytes in JSON audit.
    summary_payload = {
        "gene_symbol": gene,
        "exact_profile_count": collected.get("exact_profile_count"),
        "neural_profile_count": collected.get("neural_profile_count"),
        "subset_effect_profile_count": collected.get("subset_effect_profile_count"),
        "candidate_union_count": collected.get("candidate_union_count"),
        "candidate_retrieval_truncated": collected.get("candidate_retrieval_truncated"),
        "max_discovery_profiles": collected.get("max_discovery_profiles"),
        "max_chart_candidates": collected.get("max_chart_candidates"),
        "max_selected_profiles": collected.get("max_selected_profiles"),
        "attempt_figures": collected.get("attempt_figures"),
        "selected_profile_count": len(selected_for_present),
        "rejected_candidate_count": collected.get("rejected_candidate_count"),
        "scientific_status": collected.get("scientific_status"),
        "visual_status": collected.get("visual_status"),
        "exact_query": collected.get("exact_query"),
        "neural_query": collected.get("neural_query"),
        "subset_effect_query": collected.get("subset_effect_query"),
        "retrieved_at": collected.get("retrieved_at"),
        "selected_profiles": selected_for_present,
        "intro_text": build_intro_text(
            gene,
            exact_count=collected.get("exact_profile_count"),
            neural_count=collected.get("neural_profile_count"),
            subset_count=collected.get("subset_effect_profile_count"),
        ),
        "comparability_note": COMPARABILITY_NOTE,
        "screening_caveat": SCREENING_CAVEAT,
        "selection_policy": SELECTION_POLICY,
        "presentation_item_key": f"geo-profiles-{slugify(gene.lower())}",
    }
    write_json_atomic(attempt_dir / "summary.json", summary_payload)
    write_json_atomic(attempt_dir / "candidates.json", {"candidates": candidates_audit})
    write_json_atomic(
        attempt_dir / "selected_profiles.json",
        {"selected_profiles": selected_for_present},
    )

    scientific = str(collected.get("scientific_status") or STATUS_SOURCE_UNAVAILABLE)
    visual = str(collected.get("visual_status") or STATUS_UNAVAILABLE)
    overall = _overall_status(scientific, visual)
    # Outside-shortlist audit status must never rewrite visual_status.
    outside = sum(
        1
        for c in candidates_audit
        if c.get("graph_status") == gp.GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE
    )
    _ = outside

    evidence_records = list(state.get("evidence_records") or [])
    run_id = str(state.get("dossier_run_id") or "")
    summary_rec = _evidence(
        dossier_run_id=run_id,
        gene_symbol=gene,
        fact_type="section_3a_summary",
        key="section-3a-summary",
        value={
            "gene_symbol": gene,
            "scientific_status": scientific,
            "visual_status": visual,
            "overall_status": overall,
            "selected_profile_count": len(selected_for_present),
            "exact_profile_count": collected.get("exact_profile_count"),
            "neural_profile_count": collected.get("neural_profile_count"),
            "subset_effect_profile_count": collected.get("subset_effect_profile_count"),
            "attempt_dir": str(attempt_dir),
            "intro_text": summary_payload["intro_text"],
            "comparability_note": COMPARABILITY_NOTE,
            "screening_caveat": SCREENING_CAVEAT,
            "presentation_item_key": summary_payload["presentation_item_key"],
        },
        display_text=(
            f"{gene} GEO Profiles brain/neuron perturbation screening "
            f"({overall}; {len(selected_for_present)} selected profile records)."
        ),
    )
    _append_evidence(evidence_records, summary_rec, persist_db=persist_db)

    for present in selected_for_present:
        uid = str(present.get("profile_uid") or "")
        fig_path = present.get("figure_relative_path")
        abs_fig = str(attempt_dir / fig_path) if fig_path else None
        fig_rec = _evidence(
            dossier_run_id=run_id,
            gene_symbol=gene,
            fact_type="section_3a_profile",
            key=f"section-3a-profile-{uid}",
            value={
                **{k: v for k, v in present.items() if k != "score_components"},
                "local_artifact_path": abs_fig,
                "relative_path": abs_fig,
                "media_type": "image/png" if fig_path else None,
                "sha256": present.get("figure_sha256"),
                "width": present.get("image_width"),
                "height": present.get("image_height"),
            },
            display_text=str(present.get("title") or f"GEO Profile {uid}"),
            organism=str(present.get("taxon") or present.get("organism") or "") or None,
        )
        _append_evidence(evidence_records, fig_rec, persist_db=persist_db)

    coverage = list(state.get("coverage") or [])
    if scientific == STATUS_SUCCESS:
        cov_status = SourceStatus.success
    elif scientific == STATUS_NO_RELEVANT:
        cov_status = SourceStatus.partial
    else:
        cov_status = SourceStatus.failed
    coverage.append(
        SourceCoverageResult(
            dossier_run_id=run_id,
            source_name=gp.SOURCE_NAME,
            status=cov_status,
            evidence_record_count=1 + len(selected_for_present),
            error_message=None
            if scientific == STATUS_SUCCESS
            else f"section_3a scientific_status={scientific}",
            notes=f"visual_status={visual}; overall_status={overall}",
            report_sections_supported=[SECTION_PERTURBATIONS],
        )
    )

    section_status = {
        "section_key": "3a",
        "summary": summary_payload,
        "rendering_status": {
            "scientific_status": scientific,
            "visual_status": visual,
            "overall_status": overall,
        },
        "audit": {
            "gene_attempt_dir": str(attempt_dir),
            "artifacts": {
                "summary.json": "summary.json",
                "candidates.json": "candidates.json",
                "selected_profiles.json": "selected_profiles.json",
                **{
                    f"figures/{uid}.png": meta["relative_path"]
                    for uid, meta in figure_artifacts.items()
                },
            },
            "outside_shortlist_count": outside,
            "figure_artifacts": figure_artifacts,
            "search_status": collected.get("search_status"),
        },
    }
    write_json_atomic(attempt_dir / "section_3a_status.json", section_status)

    return {
        **state,
        "evidence_records": evidence_records,
        "coverage": coverage,
        "section_3a_status": section_status,
        "tool_results": [
            *(state.get("tool_results") or []),
            *list(collected.get("tool_results") or []),
        ],
    }


__all__ = [
    "COMPARABILITY_NOTE",
    "GEO_INTRO_TEXT",
    "SCREENING_CAVEAT",
    "SECTION_PERTURBATIONS",
    "SELECTION_POLICY",
    "STATUS_FAILED",
    "STATUS_NO_RELEVANT",
    "STATUS_NOT_ATTEMPTED",
    "STATUS_PARTIAL",
    "STATUS_SOURCE_UNAVAILABLE",
    "STATUS_SUCCESS",
    "STATUS_UNAVAILABLE",
    "SUBSECTION_3A",
    "Section3aConfig",
    "accept_scientific_complete_gene_report",
    "accept_visual_complete_gene_report",
    "build_intro_text",
    "evaluate_section_3a_scientific_complete",
    "evaluate_section_3a_visual_complete",
    "format_reporter_line",
    "node_generate_section_3a_derived_artifacts",
    "presentation_profile",
]
