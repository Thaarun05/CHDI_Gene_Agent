"""Offline tests for Section 2b Barres / Allen / BrainRNASeq."""

from __future__ import annotations

from pathlib import Path

import pytest

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
    ToolResult,
)
from gene_dossier.report_presentation import build_barres_brain_expression_blocks
from gene_dossier.report_schema import ReportContentBlock, ReportSubsection
from gene_dossier.rancho_report import (
    render_section_2b_subsection_segments,
    split_section_2b_page_segments,
)
from gene_dossier.section_2b import (
    CATEGORY_NOTE,
    NOT_DETERMINED,
    PLOT_VERSION,
    SECTION_2B_PRODUCTION_DENYLIST,
    Section2bConfig,
    classify_expression,
    format_numeric_cell,
    group_celltype_expression,
    mouse_average_tpm_from_matrix,
    select_agilent_probe,
    select_replicate_columns,
    validate_probe_lookup_rows,
)
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SUPPORTED_SECTION_BUNDLE_KEYS,
    build_section_bundle_document,
    opaque_evidence_ref,
    render_section_bundle_html,
    sources_for_sections,
    validate_section_keys,
)
from gene_dossier.tools import brainrnaseq


def test_supported_keys_include_2b_default_unchanged():
    assert SUPPORTED_SECTION_BUNDLE_KEYS == (
        "1a",
        "1b",
        "1c",
        "1d",
        "1e",
        "2a",
        "2b",
        "2c",
        "3a",
    )
    assert DEFAULT_SECTION_BUNDLE_KEYS == ("1a", "1b")
    assert validate_section_keys(["2b", "2a", "1a"]) == ["1a", "2a", "2b"]
    assert validate_section_keys(["2.b"]) == ["2b"]
    assert sources_for_sections(["2b"]) == []
    assert "Allen Brain Atlas" not in sources_for_sections(["2b"])
    assert "BrainRNASeq" not in sources_for_sections(["1a", "2b"])
    assert "GTEx" not in sources_for_sections(["2a", "2b"])


def test_section_2b_config_validation():
    Section2bConfig()
    with pytest.raises(ValueError):
        Section2bConfig(donor_well_known_file_ids=(1, 2))
    with pytest.raises(ValueError):
        Section2bConfig(allen_probe_id=0)
    with pytest.raises(ValueError):
        Section2bConfig(category_policy="invented")


def test_probe_selection_explicit_single_and_unresolved():
    probes = [{"id": 1}, {"id": 2}, {"id": 3}]
    unresolved = select_agilent_probe(probes, allen_probe_id=None)
    assert unresolved["status"] == "probe_policy_unresolved"
    assert unresolved["selected_probe_id"] is None

    single = select_agilent_probe([{"id": 42}], allen_probe_id=None)
    assert single["status"] == "selected"
    assert single["selected_probe_id"] == 42

    explicit = select_agilent_probe(probes, allen_probe_id=2)
    assert explicit["status"] == "selected"
    assert explicit["selected_probe_id"] == 2

    missing = select_agilent_probe(probes, allen_probe_id=99)
    assert missing["status"] == "probe_policy_unresolved"


def test_validate_probe_lookup_no_truncation():
    payload = {
        "msg": [
            {
                "id": 10 + i,
                "name": f"p{i}",
                "gene": {"acronym": "GENEX", "entrez_id": 9999},
            }
            for i in range(5)
        ]
    }
    out = validate_probe_lookup_rows(
        payload, gene_symbol="GENEX", entrez_gene_id="9999"
    )
    assert out["validated_probe_count"] == 5
    assert out["discovered_probe_count"] == 5


def test_categories_unresolved_do_not_invent_labels():
    assert classify_expression(100.0) is None
    assert format_numeric_cell(None) == NOT_DETERMINED
    assert format_numeric_cell(12.345) == "12.345"
    assert format_numeric_cell(12.345678901) == "12.3456789"
    assert format_numeric_cell(10.0) == "10"


def test_group_celltype_sem_and_microglla_alias():
    row = {
        "gene_id": "GENEX",
        "astrocytes_fetal_1": "1",
        "astrocytes_fetal_2": "3",
        "astrocytes_mature_1": "2",
        "neurons_1": "4",
        "oligodendrocytes_1": "0.5",
        "endothelial_1": "0.2",
        "microglla_1": "1.5",
        "microglla_2": "2.5",
        "annotation": "ignore",
    }
    from gene_dossier.section_2b import HUMAN_CELLTYPE_GROUPS

    grouped = group_celltype_expression(row, groups=HUMAN_CELLTYPE_GROUPS)
    by_label = {g["label"]: g for g in grouped["groups"]}
    assert by_label["Astrocytes Fetal"]["mean_fpkm"] == 2.0
    assert by_label["Astrocytes Fetal"]["sem_fpkm"] is not None
    assert by_label["Astrocytes Mature"]["sem_fpkm"] is None  # n=1
    assert by_label["Microglia"]["mean_fpkm"] == 2.0
    assert grouped["microglia_prefix_matched"] == "microglla_"


def test_expression_from_raw_csv_no_http():
    raw = (
        "gene_id,id,neurons_1,astrocytes_mature_1,endothelial_1,microglia_1,"
        "oligodendrocytes_1,astrocytes_fetal_1\n"
        "GENEX,9999,12.5,3.1,0.4,1.2,0.8,2.0\n"
        "OTHER,1,0.1,0.2,0.0,0.0,0.0,0.0\n"
    )
    result = brainrnaseq.expression_from_raw_csv(raw, "GENEX", species="human")
    assert result.success is True
    assert result.endpoint_name == "expression_from_raw_csv"
    assert result.data["match_count"] == 1
    assert "raw_csv" in result.data


def test_mouse_average_tpm_full_matrix():
    raw = (
        "gene_id,astrocytes_1,neurons_1\n"
        "Genex,1,3\n"
        "Other,1,1\n"
    )
    out = mouse_average_tpm_from_matrix(raw, "Genex")
    assert out["success"] is True
    assert abs(out["mean_tpm"] - 625_000.0) < 1e-6
    assert out["sample_count"] == 2
    assert out["median_tpm"] == out["mean_tpm"]
    assert set(out["selected_replicate_columns"]) == {"astrocytes_1", "neurons_1"}


def test_select_replicate_columns_requires_numeric_suffix():
    from gene_dossier.section_2b import HUMAN_CELLTYPE_GROUPS

    prefixes = [p for _, prefs in HUMAN_CELLTYPE_GROUPS for p in prefs]
    selection = select_replicate_columns(
        [
            "gene_id",
            "astrocytes_fetal_1",
            "oligodendrocytes_average_count",
            "oligodendrocytes_standard_deviation",
            "oligodendrocytes_5",
            "neurons_x",
            "unnamed:0",
        ],
        allowed_prefixes=prefixes,
    )
    assert "astrocytes_fetal_1" in selection["selected"]
    assert "oligodendrocytes_5" in selection["selected"]
    assert "oligodendrocytes_average_count" not in selection["selected"]
    assert "oligodendrocytes_standard_deviation" not in selection["selected"]
    assert "neurons_x" not in selection["selected"]
    reasons = {e["column"]: e["reason"] for e in selection["excluded"]}
    assert reasons["oligodendrocytes_average_count"] == "summary_statistic_column"
    assert reasons["oligodendrocytes_standard_deviation"] == "summary_statistic_column"


def test_summary_columns_do_not_affect_group_stats():
    from gene_dossier.section_2b import HUMAN_CELLTYPE_GROUPS

    base = {
        "gene_id": "GENEX",
        "astrocytes_fetal_1": "2",
        "astrocytes_mature_1": "4",
        "neurons_1": "1",
        "oligodendrocytes_1": "1",
        "endothelial_1": "1",
        "microglia_1": "1",
    }
    polluted = {
        **base,
        "oligodendrocytes_average_count": "999",
        "oligodendrocytes_standard_deviation": "50",
    }
    clean = group_celltype_expression(base, groups=HUMAN_CELLTYPE_GROUPS)
    dirty = group_celltype_expression(polluted, groups=HUMAN_CELLTYPE_GROUPS)
    assert clean["groups"] == dirty["groups"]
    oligo = {g["label"]: g for g in dirty["groups"]}["Oligodendrocytes"]
    assert oligo["mean_fpkm"] == 1.0
    assert oligo["n"] == 1


def test_negative_group_values_excluded():
    from gene_dossier.section_2b import HUMAN_CELLTYPE_GROUPS

    row = {
        "gene_id": "GENEX",
        "astrocytes_fetal_1": "2",
        "astrocytes_fetal_2": "-1",
        "astrocytes_mature_1": "nan",
        "neurons_1": "1",
        "oligodendrocytes_1": "1",
        "endothelial_1": "1",
        "microglia_1": "1",
    }
    grouped = group_celltype_expression(row, groups=HUMAN_CELLTYPE_GROUPS)
    by_label = {g["label"]: g for g in grouped["groups"]}
    assert by_label["Astrocytes Fetal"]["n"] == 1
    assert by_label["Astrocytes Fetal"]["mean_fpkm"] == 2.0
    assert by_label["Astrocytes Mature"]["n"] == 0
    assert any(
        e.get("reason") == "negative_or_nonfinite_value"
        for e in grouped["excluded_columns"]
    )


def test_mouse_tpm_uses_replicate_columns_only():
    raw = (
        "gene_id,astrocytes_1,neurons_1,oligodendrocytes_average_count,"
        "oligodendrocytes_standard_deviation\n"
        "Genex,1,1,999,50\n"
        "Other,1,1,1,1\n"
    )
    out = mouse_average_tpm_from_matrix(raw, "Genex")
    assert out["success"] is True
    assert set(out["selected_replicate_columns"]) == {"astrocytes_1", "neurons_1"}
    # Without summary columns: each target TPM = 1/(1+1)*1e6 = 5e5; mean 5e5
    assert abs(out["mean_tpm"] - 500_000.0) < 1e-6


def test_mouse_tpm_rejects_duplicate_rows():
    raw = (
        "gene_id,astrocytes_1\n"
        "Genex,1\n"
        "Genex - Mus musculus,2\n"
        "Other,1\n"
    )
    out = mouse_average_tpm_from_matrix(raw, "Genex")
    assert out["success"] is False
    assert out["error_type"] == "multiple_gene_matches"


def test_mouse_tpm_malformed_not_coerced_to_zero():
    raw = (
        "gene_id,astrocytes_1,neurons_1\n"
        "Genex,1,1\n"
        "Other,bad,1\n"
        "Other2,-3,1\n"
    )
    out = mouse_average_tpm_from_matrix(raw, "Genex")
    assert out["success"] is True
    # astrocytes denominator skips bad/-3 → only Genex=1 → TPM=1e6
    # neurons denominator Genex+Other+Other2 = 3 → TPM=1/3*1e6
    assert abs(out["target_tpm_by_column"]["astrocytes_1"] - 1_000_000.0) < 1e-6
    assert abs(out["target_tpm_by_column"]["neurons_1"] - (1_000_000.0 / 3.0)) < 1e-6
    assert out["skipped_cells_count"] >= 2
    assert out["malformed_row_policy"] == "skip_malformed_gene_rows_from_denominator"


def test_mouse_tpm_zero_denominator_fails_sample():
    raw = (
        "gene_id,astrocytes_1,neurons_1\n"
        "Genex,0,1\n"
        "Other,0,1\n"
    )
    out = mouse_average_tpm_from_matrix(raw, "Genex")
    # astrocytes_1 all zeros → sample fails; neurons succeeds
    assert out["success"] is True
    assert out["sample_count"] == 1
    assert "astrocytes_1" not in out["denominators"]
    assert any(
        f["column"] == "astrocytes_1"
        and f["reason"] == "zero_negative_or_nonfinite_denominator"
        for f in out["failed_samples"]
    )


def test_plot_version_is_v2():
    assert PLOT_VERSION == "section_2b_celltype_fpkm_v2"


def test_genex_human_mouse_pipeline_without_validation_gene_constants():
    """GENEX/Genex fixtures prove match/group/TPM/figure without validation genes."""
    from gene_dossier.section_2b import (
        HUMAN_CELLTYPE_GROUPS,
        MOUSE_CELLTYPE_GROUPS,
        figure_status_from_panels,
        render_celltype_figure_png,
    )

    human_raw = (
        "gene_id,id,astrocytes_fetal_1,astrocytes_mature_1,neurons_1,"
        "oligodendrocytes_1,endothelial_1,microglia_1\n"
        "GENEX,9999,1.0,2.0,3.0,0.5,0.2,1.5\n"
        "OTHER,1,0,0,0,0,0,0\n"
    )
    mouse_raw = (
        "gene_id,id,astrocytes_1,neurons_1,opc_1,newly_formed_oligodendrocyte_1,"
        "myelinating_oligodendrocyte_1,endothelial_1,microglia_macrophage_1\n"
        "Genex,8888,5,10,2,1,3,0.5,4\n"
        "Other,2,1,1,1,1,1,1,1\n"
    )
    human = brainrnaseq.expression_from_raw_csv(human_raw, "GENEX", species="human")
    mouse = brainrnaseq.expression_from_raw_csv(mouse_raw, "Genex", species="mouse")
    assert human.data["match_count"] == 1
    assert mouse.data["match_count"] == 1
    h_groups = group_celltype_expression(
        human.data["matched_rows"][0], groups=HUMAN_CELLTYPE_GROUPS
    )["groups"]
    m_groups = group_celltype_expression(
        mouse.data["matched_rows"][0], groups=MOUSE_CELLTYPE_GROUPS
    )["groups"]
    tpm = mouse_average_tpm_from_matrix(mouse_raw, "Genex")
    assert tpm["success"] is True
    png, panels = render_celltype_figure_png(
        human_symbol="GENEX",
        mouse_symbol="Genex",
        human_groups=h_groups,
        mouse_groups=m_groups,
        dpi=72,
    )
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert panels == {"human": "success", "mouse": "success"}
    assert figure_status_from_panels(panels) == "success_two_panels"

    _, panels_h = render_celltype_figure_png(
        human_symbol="GENEX",
        mouse_symbol="Genex",
        human_groups=h_groups,
        mouse_groups=None,
        dpi=72,
    )
    assert figure_status_from_panels(panels_h) == "success_human_only"

    _, panels_m = render_celltype_figure_png(
        human_symbol="GENEX",
        mouse_symbol="Genex",
        human_groups=None,
        mouse_groups=m_groups,
        dpi=72,
    )
    assert figure_status_from_panels(panels_m) == "success_mouse_only"
    assert (
        figure_status_from_panels({"human": "unavailable", "mouse": "unavailable"})
        == "unavailable_no_panels"
    )


def test_presentation_block_order_and_page_break():
    summary = EvidenceRecord(
        source_id="src-sum",
        dossier_run_id="run",
        gene_symbol="GENEX",
        section="Tissue and cell expression",
        subsection="Barres Lab RNA-Seq brain specific expression data",
        source_name="Allen Brain Atlas",
        source_type=SourceType.expression_database,
        assertion_type=AssertionType.expression,
        fact_type="section_2b_summary_table",
        evidence_grade=EvidenceGrade.B,
        value={
            "average_human_brain_agilent_expression": None,
            "average_human_brain_rnaseq_expression_tpm": 10.0,
            "average_mouse_brain_rnaseq_expression_tpm": 20.0,
            "presentation_item_key": "barres-genex",
        },
        display_text="summary",
    )
    category = EvidenceRecord(
        source_id="src-cat",
        dossier_run_id="run",
        gene_symbol="GENEX",
        section="Tissue and cell expression",
        subsection="Barres Lab RNA-Seq brain specific expression data",
        source_name="BrainRNASeq",
        source_type=SourceType.expression_database,
        assertion_type=AssertionType.expression,
        fact_type="section_2b_category_status",
        evidence_grade=EvidenceGrade.B,
        value={"policy": "threshold_policy_unresolved", "note": CATEGORY_NOTE},
        display_text=CATEGORY_NOTE,
    )
    result = build_barres_brain_expression_blocks(
        gene_symbol="GENEX",
        evidence_records=[summary, category],
        section_status={
            "rendering_status": {"overall": "partial"},
            "summary": {
                "average_human_brain_rnaseq_expression_tpm": 10.0,
                "average_mouse_brain_rnaseq_expression_tpm": 20.0,
                "category_policy": "threshold_policy_unresolved",
                "presentation_item_key": "barres-genex",
                "brainrnaseq_url": "https://brainrnaseq.org/",
            },
        },
    )
    roles = [b.presentation_role for b in result.blocks]
    assert roles[0] == "section_2b_intro"
    assert roles[1] == "section_2b_summary_table"
    assert roles[2] == "section_2b_category_status"
    assert roles[3] == "section_2b_celltype_intro"
    assert result.blocks[3].presentation_page_break_before is True
    assert roles[4] == "section_2b_source_link"
    intro = result.blocks[0].text or ""
    assert "GENEX expression results are presented in the table below" in intro
    assert "downloaded, and average expression levels were calculated" in intro
    cell_intro = result.blocks[3].text or ""
    assert "relative GENEX RNA expression" in cell_intro
    assert "Barres Lab Brain RNA-Seq" in cell_intro
    link = result.blocks[4]
    assert link.links[0]["url"] == "https://brainrnaseq.org/"
    assert "Ben Barres lab (Link)" in (link.links[0].get("label") or "")
    table = result.blocks[1]
    assert len(table.table_headers) == 5
    assert table.table_rows[0][0] == NOT_DETERMINED
    assert table.table_rows[0][3] == NOT_DETERMINED
    assert "not determined" in (result.blocks[2].text or "").lower()
    assert all(
        "brainrnaseq.org/wp-content" not in (lnk.get("url") or "")
        for b in result.blocks
        for lnk in (b.links or [])
    )
    assert opaque_evidence_ref("2b", result.blocks[0], index=0) == "ev-2b-introduction"

    segments = split_section_2b_page_segments(list(result.blocks))
    assert len(segments) == 2
    html_segments = render_section_2b_subsection_segments(
        ReportSubsection(
            key="b",
            title="Barres Lab RNA-Seq brain specific expression data",
            toc_title="Barres",
            presentation_blocks=list(result.blocks),
        )
    )
    assert len(html_segments) == 2
    assert "b. Barres Lab" in html_segments[0]
    assert "b. Barres Lab" not in html_segments[1]
    assert "Additionally" not in html_segments[0]
    assert "Additionally" in html_segments[1]


def test_focused_2b_html_includes_major_heading():
    blocks = [
        ReportContentBlock(
            kind="narrative",
            text="intro",
            presentation_role="section_2b_intro",
            presentation_item_key="barres-x",
            evidence_ref="ev-2b-introduction",
        ),
        ReportContentBlock(
            kind="table",
            table_headers=["A", "B", "C", "D", "E"],
            table_rows=[["1", "2", "3", NOT_DETERMINED, NOT_DETERMINED]],
            presentation_role="section_2b_summary_table",
            presentation_item_key="barres-x",
            evidence_ref="ev-2b-summary-table",
        ),
        ReportContentBlock(
            kind="narrative",
            text="cell",
            presentation_role="section_2b_celltype_intro",
            presentation_item_key="barres-x",
            presentation_page_break_before=True,
            evidence_ref="ev-2b-celltype-introduction",
        ),
    ]
    document, _presentation, _audit = build_section_bundle_document(
        dossier_run_id="run-2b",
        gene_symbol="GENEX",
        section_keys=["2b"],
        evidence_records=[],
        section_status_by_key={
            "2b": {
                "rendering_status": {"overall": "partial"},
                "summary": {"presentation_item_key": "barres-genex"},
            }
        },
    )
    document.sections[0].subsections[0].presentation_blocks = blocks
    html = render_section_bundle_html(document, include_major_heading=True)
    assert "2. Expression pattern by cell and tissue" in html
    assert "subsection-2b" in html
    assert html.count('<h2 class="major-heading"') == 1


def test_assembled_2a_2b_does_not_repeat_major_heading_on_2b_pages():
    from gene_dossier.report_schema import ReportCover, ReportDocument, ReportMajorSection

    blocks_2a = [
        ReportContentBlock(
            kind="narrative",
            text="gtex",
            presentation_role="section_2a_gtex_intro",
            presentation_item_key="tissue-x",
            evidence_ref="ev-2a-gtex-introduction",
        ),
        ReportContentBlock(
            kind="link",
            text="",
            links=[{"label": "GTEx brain", "url": "https://www.gtexportal.org/home/gene/X"}],
            presentation_role="section_2a_gtex_brain_link",
            presentation_item_key="tissue-x",
            presentation_page_break_before=True,
            evidence_ref="ev-2a-gtex-brain-summary",
        ),
    ]
    blocks_2b = [
        ReportContentBlock(
            kind="narrative",
            text="intro",
            presentation_role="section_2b_intro",
            presentation_item_key="barres-x",
            evidence_ref="ev-2b-introduction",
        ),
        ReportContentBlock(
            kind="narrative",
            text="cell",
            presentation_role="section_2b_celltype_intro",
            presentation_item_key="barres-x",
            presentation_page_break_before=True,
            evidence_ref="ev-2b-celltype-introduction",
        ),
    ]
    document = ReportDocument(
        dossier_run_id="run",
        gene_symbol="GENEX",
        cover=ReportCover(gene_symbol="GENEX"),
        sections=[
            ReportMajorSection(
                number=2,
                key="2",
                title="Tissue and cell expression",
                toc_title="Tissue and cell expression",
                subsections=[
                    ReportSubsection(
                        key="a",
                        title="Tissue-specific information",
                        toc_title="Tissue-specific",
                        presentation_blocks=blocks_2a,
                    ),
                    ReportSubsection(
                        key="b",
                        title="Barres Lab RNA-Seq brain specific expression data",
                        toc_title="Barres",
                        presentation_blocks=blocks_2b,
                    ),
                ],
            )
        ],
    )
    html = render_section_bundle_html(document, include_major_heading=True)
    assert html.count('<h2 class="major-heading"') == 1
    assert "section-2b-page" in html or "subsection-2b" in html
    assert "section-2a-continuation" in html or "subsection-2a" in html


def test_no_gene_specific_branches_in_section_2b_source():
    production_paths = [
        Path("src/gene_dossier/section_2b.py"),
        Path("src/gene_dossier/tools/brainrnaseq.py"),
    ]
    gene_tokens = ("SREBF2", "CDH10", "Srebf2", "Cdh10")
    for path in production_paths:
        src = path.read_text(encoding="utf-8")
        for token in gene_tokens:
            assert token not in src, f"{token} found in {path}"
        assert "fetch_hba_expression" not in src
        assert "max_probes" not in src
    # Gold tokens live only inside SECTION_2B_PRODUCTION_DENYLIST in section_2b.
    section_src = Path("src/gene_dossier/section_2b.py").read_text(encoding="utf-8")
    assert "image-45fc289b" in SECTION_2B_PRODUCTION_DENYLIST
    assert "7.023964492" in SECTION_2B_PRODUCTION_DENYLIST
    assert "closest-to-gold" in SECTION_2B_PRODUCTION_DENYLIST
    # Ensure gold numeric literals are not used as assigned display values.
    assert "average_human_brain_agilent_expression\": 7.023964492" not in section_src
    assert "brainrnaseq.org/wp-content" not in section_src
    brs_src = Path("src/gene_dossier/tools/brainrnaseq.py").read_text(encoding="utf-8")
    for token in SECTION_2B_PRODUCTION_DENYLIST:
        assert token not in brs_src, f"gold token {token} found in brainrnaseq.py"


def test_separate_api_runs_on_browser_fallback(monkeypatch, tmp_path):
    """Direct 403 persists as failure; browser success is a separate ApiRun."""
    from gene_dossier.config import Settings
    from gene_dossier.section_2b import node_generate_section_2b_derived_artifacts
    from gene_dossier.tools import allen_brain, allen_human_rnaseq

    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")

    def _forbidden(species, gene_symbol="", settings=None):
        return ToolResult(
            source_name="BrainRNASeq",
            endpoint_name="download_csv",
            success=False,
            gene_symbol=gene_symbol or species,
            request_url=brainrnaseq.csv_url_for_species(species),
            request_params={"species": species},
            status_code=403,
            error_type="access_forbidden",
            error_message="403",
            data={"species": species, "retrieval_method": "httpx_direct"},
        )

    human_csv = (
        "gene_id,id,astrocytes_fetal_1,astrocytes_mature_1,neurons_1,"
        "oligodendrocytes_1,endothelial_1,microglia_1\n"
        "GENEX,9999,1,2,3,0.5,0.2,1.5\n"
    )
    mouse_csv = (
        "gene_id,id,astrocytes_1,neurons_1,opc_1,newly_formed_oligodendrocyte_1,"
        "myelinating_oligodendrocyte_1,endothelial_1,microglia_macrophage_1\n"
        "Genex,8888,5,10,2,1,3,0.5,4\n"
        "Other,2,1,1,1,1,1,1,1\n"
    )

    def _browser(species, gene_symbol="", settings=None):
        raw = human_csv if species == "human" else mouse_csv
        url = brainrnaseq.csv_url_for_species(species)
        return ToolResult(
            source_name="BrainRNASeq",
            endpoint_name="download_csv_browser",
            success=True,
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params={"species": species},
            status_code=None,
            data={
                "raw_csv": raw,
                "retrieval_method": "official_browser_download",
                "source_url": url,
                "final_url": url,
                "browser_channel": "chrome",
                "content_type": "text/csv",
                "byte_size": len(raw.encode("utf-8")),
                "sha256": "abc",
                "capture_via": "download_event",
                "http_status_observed": False,
            },
        )

    monkeypatch.setattr(brainrnaseq, "download_csv", _forbidden)
    monkeypatch.setattr(brainrnaseq, "download_csv_via_browser", _browser)
    monkeypatch.setattr(
        allen_brain,
        "probe_lookup",
        lambda gene, settings=None: ToolResult(
            source_name="Allen Brain Atlas",
            endpoint_name="probe_lookup",
            success=True,
            gene_symbol=gene,
            request_url="https://example.test/probe",
            request_params={},
            status_code=200,
            data={"msg": [{"id": 1}, {"id": 2}]},
        ),
    )
    monkeypatch.setattr(
        allen_brain,
        "microarray_expression",
        lambda probe_id, gene_symbol="", settings=None: ToolResult(
            source_name="Allen Brain Atlas",
            endpoint_name="microarray_expression",
            success=False,
            gene_symbol=gene_symbol,
            request_url="https://example.test/expr",
            request_params={},
            error_type="skipped",
            error_message="skipped in unit test",
        ),
    )
    monkeypatch.setattr(
        allen_human_rnaseq,
        "download_donor_zip",
        lambda donor_id, gene_symbol="", settings=None: ToolResult(
            source_name="Allen Human RNA-Seq",
            endpoint_name="download_donor_zip",
            success=False,
            gene_symbol=gene_symbol,
            request_url="https://example.test/zip",
            request_params={},
            error_type="skipped",
            error_message="skipped",
        ),
    )

    state = {
        "run_type": "section_bundle",
        "selected_section_keys": ["2b"],
        "dossier_run_id": "run-brs-fallback",
        "gene_symbol": "GENEX",
        "gene_ids": {"mouse_symbol": "Genex", "entrez_gene_id": "9999"},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }
    out = node_generate_section_2b_derived_artifacts(
        state, settings=settings, persist_db=False
    )
    status = out["section_2b_status"]
    audit = status["audit"]["brainrnaseq"]
    assert audit["human"]["download"]["success"] is False
    assert audit["human"]["download"]["error_type"] == "access_forbidden"
    assert audit["human"]["browser_fallback_attempted"] is True
    assert audit["human"]["selected_retrieval_method"] == "official_browser_download"
    assert audit["human"]["browser_fallback"]["success"] is True
    assert (
        audit["human"]["download"]["api_run_id"]
        != audit["human"]["browser_fallback"]["api_run_id"]
    )
    assert audit["mouse"]["selected_retrieval_method"] == "official_browser_download"
    assert status["rendering_status"]["celltype_figure"] == "success_two_panels"
    assert status["rendering_status"]["mouse_rnaseq_tpm"] == "success"
    assert status["rendering_status"]["agilent"] == "unavailable"
    assert status["rendering_status"]["overall"] == "partial"
    assert status["summary"]["human_brain_expression_category"] == NOT_DETERMINED
    assert status["summary"]["mouse_brain_expression_category"] == NOT_DETERMINED
    fig_audit = status["audit"]["celltype_figure"]
    assert fig_audit["plot_version"] == "section_2b_celltype_fpkm_v2"
    assert fig_audit["api_run_id"] is None
    assert fig_audit["derivation_type"] == "brainrnaseq_mean_sem_barplot"
    # No fake HTTP ApiRun for the derived matplotlib figure.
    assert not any(
        getattr(r, "endpoint_name", None) == "section_2b_celltype_figure"
        for r in out.get("api_runs") or []
    )
    fig_evidence = [
        e
        for e in out.get("evidence_records") or []
        if getattr(e, "fact_type", None) == "brainrnaseq_celltype_figure"
    ]
    assert fig_evidence
    assert fig_evidence[0].api_run_id is None
    assert fig_evidence[0].value.get("plot_version") == PLOT_VERSION
    assert "figure_api_run_id" not in fig_evidence[0].value
