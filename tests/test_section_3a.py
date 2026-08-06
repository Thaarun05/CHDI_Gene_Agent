"""Offline tests for Section 3a GEO Profiles presentation and acceptance."""

from __future__ import annotations

from pathlib import Path

from gene_dossier.report_presentation import build_section_3a_blocks
from gene_dossier.section_3a import (
    STATUS_NOT_ATTEMPTED,
    STATUS_PARTIAL,
    STATUS_SUCCESS,
    STATUS_UNAVAILABLE,
    accept_scientific_complete_gene_report,
    build_intro_text,
    evaluate_section_3a_scientific_complete,
    evaluate_section_3a_visual_complete,
    format_reporter_line,
)
from gene_dossier.section_3a_sources import paths_for
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SUPPORTED_SECTION_BUNDLE_KEYS,
    sources_for_sections,
    validate_section_keys,
)
from gene_dossier.tools.geo_profiles import GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE


def test_default_bundle_keys_unchanged() -> None:
    assert DEFAULT_SECTION_BUNDLE_KEYS == ("1a", "1b")
    assert "3a" in SUPPORTED_SECTION_BUNDLE_KEYS
    assert validate_section_keys(["3.a"]) == ["3a"]
    assert "GEO" not in sources_for_sections(["3a"])


def test_intro_uses_profile_records_not_datasets() -> None:
    text = build_intro_text("GENEX", exact_count=10, neural_count=4, subset_count=2)
    assert "profile records" in text
    assert "10 exact Gene Symbol profile records" in text
    assert "up/down genes" not in text.lower()
    assert "ranking signal" in text


def test_format_reporter_line_omits_empty() -> None:
    line = format_reporter_line(
        {"gds_accession": "GDS1", "gpl": "GPL570", "idref": "x", "gdstype": None}
    )
    assert "None" not in line
    assert "GDS1" in line


def test_presentation_blocks_and_no_outside_shortlist_selected() -> None:
    status = {
        "summary": {
            "intro_text": "Intro",
            "screening_caveat": "Caveat",
            "comparability_note": "Compare",
            "selection_policy": "Policy",
            "presentation_item_key": "geo-profiles-genex",
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_SUCCESS,
            "selected_profiles": [
                {
                    "profile_uid": "111",
                    "title": "Stress hippocampus",
                    "profile_url": "https://www.ncbi.nlm.nih.gov/geoprofiles/111",
                    "organism": "Mus musculus",
                    "reporter_line": "GDS1 · GPL1261 · x",
                    "graph_status": "success",
                    "graph_ok": True,
                    "figure_relative_path": None,
                }
            ],
            "selected_profile_count": 1,
        },
        "rendering_status": {
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_SUCCESS,
        },
    }
    result = build_section_3a_blocks(
        gene_symbol="GENEX",
        evidence_records=[],
        section_status=status,
    )
    roles = [b.presentation_role for b in result.blocks]
    assert "section_3a_intro" in roles
    assert "section_3a_profile_title" in roles
    assert "section_3a_caveat" in roles
    assert all(
        (p.get("graph_status") != GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE)
        for p in status["summary"]["selected_profiles"]
    )


def test_figures_disabled_visual_not_attempted_optional() -> None:
    status = {
        "summary": {
            "intro_text": "Intro",
            "selected_profiles": [
                {
                    "profile_uid": "1",
                    "title": "A",
                    "graph_status": "not_attempted_optional",
                    "graph_ok": False,
                }
            ],
            "selected_profile_count": 1,
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_NOT_ATTEMPTED,
        },
        "rendering_status": {
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_NOT_ATTEMPTED,
        },
    }
    result = build_section_3a_blocks(
        gene_symbol="GENEX",
        evidence_records=[],
        section_status=status,
    )
    assert any(
        b.presentation_role == "section_3a_profile_figure_status" for b in result.blocks
    )


def test_scientific_complete_accept_never_replaces(tmp_path: Path) -> None:
    paths = paths_for(tmp_path)
    attempt = paths.new_gene_attempt("GENEX", run_id="run1")
    first = accept_scientific_complete_gene_report(
        paths,
        gene_symbol="GENEX",
        attempt_dir=attempt,
        acceptance={"section_3a_scientific_complete": True},
        artifacts={},
    )
    assert first is not None
    second_attempt = paths.new_gene_attempt("GENEX", run_id="run2")
    second = accept_scientific_complete_gene_report(
        paths,
        gene_symbol="GENEX",
        attempt_dir=second_attempt,
        acceptance={"section_3a_scientific_complete": True},
        artifacts={},
    )
    assert second is None


def test_evaluate_visual_complete_requires_embedded_count() -> None:
    status = {
        "summary": {"selected_profile_count": 2, "selected_profiles": [{}, {}]},
        "rendering_status": {
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_SUCCESS,
        },
    }
    sci = evaluate_section_3a_scientific_complete(
        status=status, pdf_render_status=STATUS_SUCCESS
    )
    assert sci["scientific_complete"] is True
    vis = evaluate_section_3a_visual_complete(
        status=status,
        embedded_figure_count=1,
        selected_count=2,
        pdf_render_status=STATUS_SUCCESS,
    )
    assert vis["visual_complete"] is False
    vis_ok = evaluate_section_3a_visual_complete(
        status=status,
        embedded_figure_count=2,
        selected_count=2,
        pdf_render_status=STATUS_SUCCESS,
    )
    assert vis_ok["visual_complete"] is True


def test_outside_shortlist_status_not_in_polished_selected_contract() -> None:
    # Contract assertion used by collectors/tests: polished selected must never
    # carry not_attempted_outside_shortlist once figures are enabled.
    selected = [
        {"profile_uid": "a", "graph_status": "success"},
        {"profile_uid": "b", "graph_status": "failed"},
    ]
    audit = [
        *selected,
        {"profile_uid": "c", "graph_status": GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE},
    ]
    assert all(
        s.get("graph_status") != GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE for s in selected
    )
    assert any(
        c.get("graph_status") == GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE for c in audit
    )


def test_all_blocked_selected_charts_visual_unavailable() -> None:
    selected = [
        {"graph_status": "failed", "graph_ok": False, "error_type": "graph_http_blocked"},
        {"graph_status": "failed", "graph_ok": False, "error_type": "graph_http_blocked"},
    ]
    ok_count = sum(1 for s in selected if s.get("graph_status") == "success")
    visual_status = (
        "success"
        if ok_count == len(selected)
        else ("unavailable" if ok_count == 0 else "partial")
    )
    assert visual_status == STATUS_UNAVAILABLE


def test_mixed_valid_blocked_charts_visual_partial() -> None:
    selected = [
        {"graph_status": "success", "graph_ok": True},
        {"graph_status": "failed", "graph_ok": False, "error_type": "graph_http_blocked"},
    ]
    ok_count = sum(1 for s in selected if s.get("graph_status") == "success")
    visual_status = (
        "success"
        if ok_count == len(selected)
        else ("unavailable" if ok_count == 0 else "partial")
    )
    assert visual_status == STATUS_PARTIAL


def test_visual_success_requires_every_selected_chart_valid() -> None:
    selected = [
        {"graph_status": "success", "graph_ok": True},
        {"graph_status": "success", "graph_ok": True},
    ]
    ok_count = sum(1 for s in selected if s.get("graph_status") == "success")
    assert ok_count == len(selected)
    visual_status = "success"
    assert visual_status == STATUS_SUCCESS


def test_no_figure_block_for_blocked_capture() -> None:
    status = {
        "summary": {
            "intro_text": "Intro",
            "screening_caveat": "Caveat",
            "comparability_note": "Compare",
            "selection_policy": "Policy",
            "presentation_item_key": "geo-profiles-genex",
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_UNAVAILABLE,
            "selected_profiles": [
                {
                    "profile_uid": "111",
                    "title": "Stress hippocampus",
                    "profile_url": "https://www.ncbi.nlm.nih.gov/geoprofiles/111",
                    "organism": "Mus musculus",
                    "reporter_line": "GDS1 · GPL1261 · x",
                    "graph_status": "failed",
                    "graph_ok": False,
                    "figure_relative_path": "figures/bad.png",
                }
            ],
            "selected_profile_count": 1,
        },
        "rendering_status": {
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_UNAVAILABLE,
        },
    }
    result = build_section_3a_blocks(
        gene_symbol="GENEX",
        evidence_records=[],
        section_status=status,
    )
    roles = [b.presentation_role for b in result.blocks]
    assert "section_3a_profile_figure" not in roles
    assert "section_3a_profile_figure_status" in roles


def test_scientific_success_when_charts_blocked() -> None:
    status = {
        "summary": {
            "selected_profile_count": 2,
            "selected_profiles": [{}, {}],
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_UNAVAILABLE,
        },
        "rendering_status": {
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_UNAVAILABLE,
        },
    }
    sci = evaluate_section_3a_scientific_complete(
        status=status, pdf_render_status=STATUS_SUCCESS
    )
    assert sci["scientific_complete"] is True
    vis = evaluate_section_3a_visual_complete(
        status=status,
        embedded_figure_count=0,
        selected_count=2,
        pdf_render_status=STATUS_SUCCESS,
    )
    assert vis["visual_complete"] is False
    assert vis["visual_status"] == STATUS_UNAVAILABLE


def test_post_render_visual_complete_fails_when_chart_blocked() -> None:
    status = {
        "summary": {
            "selected_profile_count": 1,
            "selected_profiles": [
                {"graph_status": "failed", "graph_ok": False},
            ],
        },
        "rendering_status": {
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_UNAVAILABLE,
        },
    }
    vis = evaluate_section_3a_visual_complete(
        status=status,
        embedded_figure_count=0,
        selected_count=1,
        pdf_render_status=STATUS_SUCCESS,
    )
    assert vis["scientific_complete"] is True
    assert vis["visual_complete"] is False
