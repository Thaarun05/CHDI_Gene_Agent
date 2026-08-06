"""Offline tests for Section 1e homologues / ortholog helpers."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.report_presentation import build_homologues_blocks, build_section_presentation
from gene_dossier.section_1e import (
    CAPTURE_ORIGIN_LIVE,
    MODEL_SPECIES_PRIORITY,
    REQUIRED_HEADERS,
    Section1eConfig,
    build_complete_narrative,
    build_incomplete_narrative,
    classify_membership,
    node_generate_section_1e_derived_artifacts,
    ortholog_ncbi_legacy_link,
    ortholog_ncbi_url,
    select_species_names,
)
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SUPPORTED_SECTION_BUNDLE_KEYS,
    assign_opaque_refs,
    render_section_bundle_html,
    sources_for_sections,
    validate_section_keys,
)
from gene_dossier.report_schema import (
    ReportContentBlock,
    ReportCover,
    ReportDocument,
    ReportMajorSection,
    ReportSubsection,
)


def test_default_bundle_keys_unchanged():
    assert DEFAULT_SECTION_BUNDLE_KEYS == (
        "1a",
        "1b",
        "1c",
        "1d",
        "1e",
        "2a",
        "2b",
        "2c",
        "3a",
        "4a",
    )
    assert "1e" in SUPPORTED_SECTION_BUNDLE_KEYS
    assert validate_section_keys(["1e", "1a"]) == ["1a", "1e"]
    assert sources_for_sections(["1e"]) == []


def test_section_1e_config_rejects_unknown_scope():
    with pytest.raises(ValueError, match="Unsupported"):
        Section1eConfig(ortholog_scope_tax_id=999999)
    cfg = Section1eConfig(ortholog_scope_tax_id=7776)
    assert cfg.ortholog_scope_label == "jawed vertebrates"
    assert cfg.ncbi_taxon_name == "Gnathostomata"
    assert Section1eConfig(ortholog_scope_tax_id=32523).ortholog_scope_label == "tetrapods"
    assert Section1eConfig(ortholog_scope_tax_id=32523).ncbi_taxon_name == "Tetrapoda"
    assert not hasattr(Section1eConfig(), "ortholog_capture_image")


def test_section_1e_scopes_include_ncbi_taxon_names():
    from gene_dossier.section_1e import SUPPORTED_SECTION_1E_SCOPES

    assert SUPPORTED_SECTION_1E_SCOPES[7776]["ncbi_taxon_name"] == "Gnathostomata"
    assert SUPPORTED_SECTION_1E_SCOPES[32523]["ncbi_taxon_name"] == "Tetrapoda"


def test_classify_membership_hierarchy():
    assert classify_membership(
        {"gene_groups": [{"id": "6721"}], "query_gene_ids": ["6721"]},
        query_gene_id="6721",
    ) == ("explicit", None)
    assert classify_membership(
        {"gene_groups": [{"id": "999"}], "query_gene_ids": ["6721"]},
        query_gene_id="6721",
    )[0] == "rejected"
    assert classify_membership(
        {"gene_groups": [], "query_gene_ids": ["6721"]},
        query_gene_id="6721",
    ) == ("endpoint_implicit", None)
    assert classify_membership(
        {"gene_groups": None, "query_gene_ids": ["1"]},
        query_gene_id="6721",
    )[0] == "unverified"


def test_select_species_names_priority_then_alpha():
    records = [
        {"tax_id": 9615, "common_name": "dog"},
        {"tax_id": 10090, "common_name": "house mouse"},
        {"tax_id": 7955, "common_name": "zebrafish"},
        {"tax_id": 9606, "common_name": "human"},
        {"tax_id": 9031, "common_name": "chicken"},
        {"tax_id": 12345, "common_name": "aardvark"},
        {"tax_id": 12346, "common_name": "zebra"},
    ]
    names = select_species_names(records, limit=10)
    assert "human" not in names
    assert names[0] == "house mouse"
    assert names[1] == "zebrafish"
    assert names[2] == "chicken"
    assert names.index("dog") < names.index("aardvark")
    assert MODEL_SPECIES_PRIORITY[0] == 10090


def test_narratives_avoid_homologene_and_hardcoded_counts():
    complete = build_complete_narrative(
        gene_symbol="SREBF2",
        scope_label="jawed vertebrates",
        ortholog_gene_count=42,
        species_names=["house mouse", "Norway rat"],
    )
    incomplete = build_incomplete_narrative(
        gene_symbol="CDH10", scope_label="tetrapods"
    )
    for text in (complete, incomplete):
        assert "HomoloGene" not in text
        assert "562" not in text
        assert "210" not in text
    assert "42 orthologous genes" in complete
    assert "house mouse and Norway rat" in complete
    assert "exact ortholog count is not reported" in incomplete


def test_ortholog_ncbi_url_datasets():
    url = ortholog_ncbi_url(entrez_gene_id="6721")
    assert url == "https://www.ncbi.nlm.nih.gov/datasets/gene/6721/#orthologs"
    with pytest.raises(ValueError, match="numeric"):
        ortholog_ncbi_url(entrez_gene_id="SREBF2")


def test_ortholog_ncbi_legacy_link_audit_only():
    url = ortholog_ncbi_legacy_link(
        entrez_gene_id="6721", scope_tax_id=7776, gene_symbol="SREBF2"
    )
    assert url == (
        "https://www.ncbi.nlm.nih.gov/gene/6721/ortholog/?scope=7776&term=SREBF2"
    )


def test_production_sources_forbid_gold_import_paths():
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "src" / "gene_dossier" / "section_1e.py",
        root / "scripts" / "run_section_bundle.py",
    ]
    forbidden = (
        "section_1e_gold_captures",
        "provided_ortholog_capture_image",
        "ortholog_capture_image",
        "capture_origin=provided",
    )
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for token in forbidden:
            assert token not in text, f"{token} found in {path}"


def test_cli_has_no_ortholog_capture_image_flag():
    root = Path(__file__).resolve().parents[1]
    text = (root / "scripts" / "run_section_bundle.py").read_text(encoding="utf-8")
    assert "--ortholog-capture-image" not in text
    assert "ortholog_capture_image" not in text


def test_headers_cover_required_is_subset():
    from gene_dossier.section_1e import _headers_cover_required

    visible = [
        "Checkbox",
        "Scientific name",
        "Symbol",
        "Length (aa)",
        "Architecture",
        "Action",
    ]
    assert _headers_cover_required(visible) is True
    assert _headers_cover_required(["Scientific name", "Symbol"]) is False


def test_parse_row_cells_by_header_index():
    from gene_dossier.section_1e import _parse_row_cells_by_header

    headers = ["Scientific name", "Symbol", "Length (aa)", "Architecture", "Action"]
    cells = ["Mus musculus", "Srebf2", "1138", "", "…"]
    parsed = _parse_row_cells_by_header(headers, cells)
    assert parsed["scientific name"] == "Mus musculus"
    assert parsed["symbol"] == "Srebf2"
    assert parsed["length (aa)"] == "1138"
    assert parsed["architecture"] == ""


def test_human_conditional_displayed_count():
    from gene_dossier.section_1e import (
        _count_claim_ready,
        _count_consistency_passed,
        _expected_displayed_gene_count,
        _human_reference_row_detected,
    )

    rows_with_human = [
        {"gene_id": "6721", "tax_id": "9606"},
        {"gene_id": "20788", "tax_id": "10090"},
    ]
    assert _human_reference_row_detected(
        rows_with_human, resolved_entrez_gene_id="6721"
    )
    assert _expected_displayed_gene_count(
        scoped_ortholog_gene_count=305, human_reference_row_detected=True
    ) == 306
    assert _expected_displayed_gene_count(
        scoped_ortholog_gene_count=305, human_reference_row_detected=False
    ) == 305
    assert (
        _count_consistency_passed(
            retrieval_complete=True,
            displayed_gene_count=306,
            expected_displayed_count=306,
        )
        is True
    )
    assert (
        _count_consistency_passed(
            retrieval_complete=True,
            displayed_gene_count=305,
            expected_displayed_count=306,
        )
        is False
    )
    assert (
        _count_consistency_passed(
            retrieval_complete=False,
            displayed_gene_count=999,
            expected_displayed_count=306,
        )
        is None
    )
    assert _count_claim_ready(pagination_complete=True, taxonomy_complete=True) is True
    assert _count_claim_ready(pagination_complete=True, taxonomy_complete=False) is False
    assert _count_claim_ready(pagination_complete=False, taxonomy_complete=True) is False


def test_official_capture_gate_and_provenance_from_bytes():
    from gene_dossier.section_1e import _official_capture_gate

    png = b"\x89PNG\r\n\x1a\n" + b"\x00" * 100
    digest = hashlib.sha256(png).hexdigest()
    meta = {
        "capture_origin": CAPTURE_ORIGIN_LIVE,
        "sha256": digest,
        "entrez_gene_id": "6721",
        "selected_scope_tax_id": 7776,
        "visible_headers": list(REQUIRED_HEADERS) + ["Action"],
        "displayed_gene_count": 10,
        "expected_displayed_count": 10,
    }
    assert _official_capture_gate(
        capture_api_success=True,
        capture_metadata=meta,
        captured_bytes=png,
        resolved_entrez_gene_id="6721",
        configured_scope_tax_id=7776,
        retrieval_complete=True,
    )
    bad = {**meta, "capture_origin": "provided"}
    assert not _official_capture_gate(
        capture_api_success=True,
        capture_metadata=bad,
        captured_bytes=png,
        resolved_entrez_gene_id="6721",
        configured_scope_tax_id=7776,
        retrieval_complete=True,
    )
    tampered = {**meta, "sha256": "0" * 64}
    assert not _official_capture_gate(
        capture_api_success=True,
        capture_metadata=tampered,
        captured_bytes=png,
        resolved_entrez_gene_id="6721",
        configured_scope_tax_id=7776,
        retrieval_complete=True,
    )


def test_hard_count_gate_rejects_when_retrieval_complete():
    """When pagination+taxonomy are complete, mismatched counts reject official."""
    from gene_dossier.section_1e import _official_capture_gate

    png = b"\x89PNG\r\n\x1a\n" + b"\x00" * 100
    digest = hashlib.sha256(png).hexdigest()
    meta = {
        "capture_origin": CAPTURE_ORIGIN_LIVE,
        "sha256": digest,
        "entrez_gene_id": "6721",
        "selected_scope_tax_id": 7776,
        "visible_headers": list(REQUIRED_HEADERS),
        "displayed_gene_count": 847,
        "expected_displayed_count": 727,
    }
    assert not _official_capture_gate(
        capture_api_success=True,
        capture_metadata=meta,
        captured_bytes=png,
        resolved_entrez_gene_id="6721",
        configured_scope_tax_id=7776,
        retrieval_complete=True,
    )
    # Matching counts pass.
    ok = {**meta, "displayed_gene_count": 727, "expected_displayed_count": 727}
    assert _official_capture_gate(
        capture_api_success=True,
        capture_metadata=ok,
        captured_bytes=png,
        resolved_entrez_gene_id="6721",
        configured_scope_tax_id=7776,
        retrieval_complete=True,
    )


def test_incomplete_retrieval_skips_hard_count_rejection():
    """Incomplete taxonomy/retrieval must not reject a live official capture."""
    from gene_dossier.section_1e import _count_consistency_passed, _official_capture_gate

    png = b"\x89PNG\r\n\x1a\n" + b"\x00" * 100
    digest = hashlib.sha256(png).hexdigest()
    meta = {
        "capture_origin": CAPTURE_ORIGIN_LIVE,
        "sha256": digest,
        "entrez_gene_id": "6721",
        "selected_scope_tax_id": 7776,
        "visible_headers": list(REQUIRED_HEADERS) + ["Feedback"],
        "displayed_gene_count": 847,
        # Undercounted scoped expectation from failed taxonomy.
        "expected_displayed_count": 727,
    }
    assert (
        _count_consistency_passed(
            retrieval_complete=False,
            displayed_gene_count=847,
            expected_displayed_count=727,
        )
        is None
    )
    assert _official_capture_gate(
        capture_api_success=True,
        capture_metadata=meta,
        captured_bytes=png,
        resolved_entrez_gene_id="6721",
        configured_scope_tax_id=7776,
        retrieval_complete=False,
    )

def test_datasets_gene_route_ok_ignores_fragment():
    from gene_dossier.section_1e import _datasets_gene_route_ok

    assert _datasets_gene_route_ok(
        "https://www.ncbi.nlm.nih.gov/datasets/gene/6721/",
        entrez_gene_id="6721",
    )
    assert _datasets_gene_route_ok(
        "https://www.ncbi.nlm.nih.gov/datasets/gene/6721/#orthologs",
        entrez_gene_id="6721",
    )
    assert not _datasets_gene_route_ok(
        "https://www.ncbi.nlm.nih.gov/gene/6721/ortholog/?scope=7776",
        entrez_gene_id="6721",
    )


def test_orthologs_tab_selected_uses_aria_not_fragment_only():
    from gene_dossier.section_1e import _orthologs_tab_selected

    class _FakeTab:
        def __init__(self, aria: str):
            self._aria = aria

        def count(self):
            return 1

        def get_attribute(self, name):
            return self._aria if name == "aria-selected" else None

    class _EmptyLoc:
        def count(self):
            return 0

        @property
        def first(self):
            return self

    class _FakePage:
        def __init__(self, aria: str, *, with_chrome: bool = False):
            self._aria = aria
            self._with_chrome = with_chrome

        def locator(self, sel):
            sel_l = sel.lower()
            if "orthologs" in sel_l and ("role='tab'" in sel.replace('"', "'") or "button" in sel_l):
                return type("L", (), {"first": _FakeTab(self._aria)})()
            if self._with_chrome:
                return type(
                    "L",
                    (),
                    {
                        "first": type("N", (), {"count": lambda self: 1})(),
                        "count": lambda self: 1,
                    },
                )()
            return _EmptyLoc()

    assert _orthologs_tab_selected(_FakePage("true")) is True
    assert _orthologs_tab_selected(_FakePage("false")) is False
    assert _orthologs_tab_selected(_FakePage("false", with_chrome=True)) is True


def _summary_record(**overrides):
    value = {
        "entrez_gene_id": "6721",
        "gene_symbol": "SREBF2",
        "scope_tax_id": 7776,
        "scope_label": "jawed vertebrates",
        "retrieval_complete": True,
        "count_status": "complete",
        "scoped_ortholog_gene_count": 3,
        "scoped_species_count": 3,
        "species_names": ["house mouse", "Norway rat", "zebrafish"],
        "ncbi_url": ortholog_ncbi_url(entrez_gene_id="6721"),
        "orthodb_url": "https://www.orthodb.org/?ncbi=6721",
        "table_status": "complete_api_fallback",
        "narrative": build_complete_narrative(
            gene_symbol="SREBF2",
            scope_label="jawed vertebrates",
            ortholog_gene_count=3,
            species_names=["house mouse", "Norway rat", "zebrafish"],
        ),
        "fallback_rows": [
            {
                "species": "house mouse",
                "gene": "Srebf2",
                "description": "sterol regulatory element binding transcription factor 2",
                "gene_id": "20788",
                "tax_id": "10090",
            }
        ],
        "presentation_item_key": "orthologs-srebf2",
    }
    value.update(overrides)
    return EvidenceRecord(
        source_id="src-1e-summary",
        dossier_run_id="run-1e",
        gene_symbol="SREBF2",
        official_symbol="SREBF2",
        section="Homologues",
        subsection="Homologues in model animals",
        source_name="NCBI Datasets",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="ortholog_collection_summary",
        species="human",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.B,
        value=value,
        display_text=str(value["narrative"]),
    )


def test_build_homologues_blocks_fallback_table():
    result = build_homologues_blocks(
        gene_symbol="SREBF2",
        evidence_records=[_summary_record()],
    )
    roles = [b.presentation_role for b in result.blocks]
    assert roles[0] == "section_1e_narrative"
    assert "section_1e_fallback_table" in roles
    assert "section_1e_attribution" in roles
    assert "HomoloGene" not in (result.blocks[0].text or "")
    polished, ref_map = assign_opaque_refs(section_key="1e", blocks=result.blocks)
    assert polished[0].evidence_ref == "ev-1e-ortholog-summary"
    assert "ev-1e-ortholog-fallback-table" in ref_map
    assert "ev-1e-ortholog-attribution" in ref_map


def test_official_capture_presentation_excludes_fallback(tmp_path: Path, monkeypatch):
    from gene_dossier.config import Settings, get_settings
    from gene_dossier.rancho_report import _rancho_css, _render_section_1e_blocks
    from gene_dossier.ucsc_figure import sha256_hex

    raw_root = tmp_path / "raw"
    raw_root.mkdir(parents=True)
    settings = Settings(raw_data_dir=raw_root, output_dir=tmp_path / "out")
    monkeypatch.setattr("gene_dossier.ucsc_figure.get_settings", lambda: settings)
    get_settings.cache_clear()

    try:
        from PIL import Image

        rel = "section_1e/ortholog-capture.png"
        img_path = raw_root / rel
        img_path.parent.mkdir(parents=True, exist_ok=True)
        Image.new("RGB", (700, 400), color=(240, 240, 240)).save(img_path)
        digest = sha256_hex(img_path.read_bytes())
    except Exception:
        pytest.skip("Pillow required for official capture presentation test")

    capture = EvidenceRecord(
        source_id="src-1e-capture",
        dossier_run_id="run-1e",
        gene_symbol="SREBF2",
        official_symbol="SREBF2",
        section="Homologues",
        subsection="Homologues in model animals",
        source_name="NCBI Gene",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="ortholog_table_capture",
        species="human",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.B,
        value={
            "status": "success",
            "relative_path": rel,
            "media_type": "image/png",
            "width": 700,
            "height": 400,
            "sha256": digest,
            "byte_size": img_path.stat().st_size,
            "presentation_item_key": "orthologs-srebf2",
            "visible_row_count": 10,
            "capture_origin": CAPTURE_ORIGIN_LIVE,
        },
        display_text="capture",
        raw_artifact_id="raw-cap",
    )
    summary = _summary_record(table_status="official_capture")
    result = build_homologues_blocks(
        gene_symbol="SREBF2",
        evidence_records=[summary, capture],
    )
    roles = [b.presentation_role for b in result.blocks]
    assert "section_1e_ortholog_capture" in roles
    assert "section_1e_fallback_table" not in roles
    fig = next(b for b in result.blocks if b.presentation_role == "section_1e_ortholog_capture")
    assert fig.kind == "figure"
    assert fig.figure_path == rel
    html = _render_section_1e_blocks(list(result.blocks))
    assert "section-1e-ortholog-capture" in html
    assert "section-1e-fallback-table" not in html
    assert "<th>Architecture</th>" not in html
    css = _rancho_css()
    assert "4.8in" in css
    assert "72%" in css
    assert "4.7in" not in css
    get_settings.cache_clear()


def test_dom_consistency_allows_human_and_name_warnings():
    from gene_dossier.section_1e import _dom_consistency

    result = _dom_consistency(
        row_meta=[
            {"gene_id": "6721", "tax_id": "9606", "scientific_name": "Homo sapiens"},
            {"gene_id": "20788", "tax_id": "10090", "scientific_name": "Mus musculus"},
            {"gene_id": "", "tax_id": "", "scientific_name": "mouse"},
        ],
        scoped_records=[
            {"gene_id": "20788", "tax_id": "10090", "scientific_name": "Mus musculus"},
        ],
        human_entrez_gene_id="6721",
    )
    assert result["status"] == "pass"
    assert result["gene_ok"] is True
    assert result["tax_ok"] is True


def test_classify_capture_failure_classes():
    from gene_dossier.section_1e import _CaptureFailure, _classify_capture_exception

    assert _classify_capture_exception(
        _CaptureFailure("blocked_or_consent_page", "403 Error - NCBI")
    )[0] == "blocked_or_consent_page"
    assert _classify_capture_exception(
        _CaptureFailure("table_selector_failure", "Scientific name/Symbol missing")
    )[0] == "table_selector_failure"
    assert _classify_capture_exception(
        ValueError("capture dimensions too small (10x10)")
    )[0] == "image_quality_rejection"


def test_section_1e_pdf_page_has_heading_without_ucsc(tmp_path: Path):
    pytest.importorskip("fitz")
    from gene_dossier.rancho_report import (
        SECTION_1C_PDF_PAGE_BREAK,
        _split_pdf_page_segments,
        render_rancho_pdf,
    )
    from gene_dossier.section_bundle import render_section_bundle_html

    tall = ReportContentBlock(
        kind="narrative",
        text="UCSC Genome Browser conservation view.",
        presentation_role=None,
    )
    narrative = ReportContentBlock(
        kind="narrative",
        text="NCBI Orthologs provides comparative information for human SREBF2.",
        presentation_role="section_1e_narrative",
        presentation_item_key="orthologs-srebf2",
    )
    document = ReportDocument(
        gene_symbol="SREBF2",
        dossier_run_id="run-1e-pdf",
        cover=ReportCover(gene_symbol="SREBF2", chromosome="22"),
        sections=[
            ReportMajorSection(
                number=1,
                key="1",
                title="Gene",
                toc_title="GENE",
                subsections=[
                    ReportSubsection(
                        key="b",
                        title="Gene conservation",
                        toc_title="GENE CONSERVATION",
                        presentation_blocks=[tall],
                        status="populated",
                    ),
                    ReportSubsection(
                        key="e",
                        title="Homologues in model animals",
                        toc_title="HOMOLOGUES IN MODEL ANIMALS",
                        presentation_blocks=[narrative],
                        status="populated",
                    ),
                ],
            )
        ],
    )
    html = render_section_bundle_html(document, include_page_chrome=False)
    assert "section-1e-page" in html
    assert SECTION_1C_PDF_PAGE_BREAK in html
    segments = _split_pdf_page_segments(html)
    assert len(segments) >= 2
    assert "Gene conservation" in segments[0]
    assert "Homologues in model animals" in segments[-1]
    assert "UCSC Genome Browser conservation view" not in segments[-1]

    pdf_path = tmp_path / "section_1.pdf"
    rendered = render_rancho_pdf(
        html, pdf_path, stamp_page_chrome=False, stamp_cover=False
    )
    assert rendered is not None and Path(rendered).is_file()

    import fitz

    with fitz.open(str(rendered)) as doc:
        assert doc.page_count >= 2
        texts = [doc.load_page(i).get_text() for i in range(doc.page_count)]
        e_pages = [i for i, text in enumerate(texts) if "Homologues in model animals" in text]
        b_pages = [i for i, text in enumerate(texts) if "Gene conservation" in text]
        assert e_pages and b_pages
        assert min(e_pages) > max(b_pages)
        assert "UCSC Genome Browser conservation view" not in texts[min(e_pages)]


def test_build_section_presentation_routes_1e():
    result = build_section_presentation(
        section_key="1e",
        gene_symbol="SREBF2",
        evidence_records=[_summary_record()],
    )
    assert result.blocks
    assert result.blocks[0].presentation_role == "section_1e_narrative"


def test_section_1e_starts_on_new_page_when_assembled(tmp_path: Path):
    from gene_dossier.rancho_report import SECTION_1C_PDF_PAGE_BREAK

    narrative = ReportContentBlock(
        kind="narrative",
        text="NCBI Orthologs provides comparative information for human SREBF2.",
        presentation_role="section_1e_narrative",
        presentation_item_key="orthologs-srebf2",
    )
    document = ReportDocument(
        gene_symbol="SREBF2",
        dossier_run_id="run-1e-page",
        cover=ReportCover(gene_symbol="SREBF2", chromosome="22"),
        sections=[
            ReportMajorSection(
                number=1,
                key="1",
                title="Gene",
                toc_title="GENE",
                subsections=[
                    ReportSubsection(
                        key="a",
                        title="Gene aliases",
                        toc_title="GENE ALIASES",
                        presentation_blocks=[
                            ReportContentBlock(
                                kind="narrative",
                                text="aliases",
                                presentation_role=None,
                            )
                        ],
                        status="populated",
                    ),
                    ReportSubsection(
                        key="e",
                        title="Homologues in model animals",
                        toc_title="HOMOLOGUES IN MODEL ANIMALS",
                        presentation_blocks=[narrative],
                        status="populated",
                    ),
                ],
            )
        ],
    )
    html = render_section_bundle_html(document, include_page_chrome=False)
    assert SECTION_1C_PDF_PAGE_BREAK in html
    a_marker = 'class="report-subsection subsection-a"'
    e_marker = 'class="report-subsection subsection-e"'
    assert html.index(a_marker) < html.index(SECTION_1C_PDF_PAGE_BREAK)
    assert html.index(SECTION_1C_PDF_PAGE_BREAK) < html.index(e_marker)
    assert "Homologues in model animals" in html


def _mock_ortholog_network(monkeypatch):
    from gene_dossier.models import ToolResult

    page = ToolResult(
        source_name="NCBI Datasets",
        endpoint_name="orthologs_by_gene_id",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/orthologs",
        request_params={},
        status_code=200,
        data={
            "reports": [
                {
                    "gene": {
                        "gene_id": "6721",
                        "symbol": "SREBF2",
                        "tax_id": 9606,
                        "taxname": "Homo sapiens",
                        "common_name": "human",
                        "description": "sterol regulatory element binding transcription factor 2",
                        "gene_groups": [{"id": "6721"}],
                    },
                    "query": ["6721"],
                },
                {
                    "gene": {
                        "gene_id": "20788",
                        "symbol": "Srebf2",
                        "tax_id": 10090,
                        "taxname": "Mus musculus",
                        "common_name": "house mouse",
                        "description": "sterol regulatory element binding transcription factor 2",
                        "gene_groups": [{"id": "6721"}],
                    },
                    "query": ["6721"],
                },
            ]
        },
    )

    def fake_iter(entrez, **kwargs):
        return [page], {
            "retrieval_complete": True,
            "pages_fetched": 1,
            "stop_reason": "exhausted",
        }

    def fake_tax(tax_ids, *, scope_tax_id, gene_symbol="", settings=None):
        tids = [str(tid) for tid in tax_ids]
        membership = {tid: True for tid in tids}
        return (
            membership,
            [],
            {
                "requested_tax_ids": tids,
                "requested_count": len(tids),
                "resolved_count": len(tids),
                "unknown_count": 0,
                "unresolved_count": 0,
                "failed_request_count": 0,
                "taxonomy_complete": True,
                "scope_tax_id": scope_tax_id,
            },
        )
    def fake_release(**kwargs):
        return ToolResult(
            source_name="OrthoDB",
            endpoint_name="orthodb_release_id",
            success=True,
            gene_symbol="SREBF2",
            request_url="https://example.test/release",
            request_params={},
            status_code=200,
            data="v12.2",
        )

    def fake_search(entrez, **kwargs):
        return ToolResult(
            source_name="OrthoDB",
            endpoint_name="genesearch",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://example.test/genesearch",
            request_params={},
            error_type="http_error",
            error_message="HTTP 500",
        )

    monkeypatch.setattr(
        "gene_dossier.section_1e.ncbi_datasets.iter_ortholog_pages", fake_iter
    )
    monkeypatch.setattr(
        "gene_dossier.section_1e.ncbi_taxonomy.resolve_taxonomy_memberships",
        fake_tax,
    )
    monkeypatch.setattr(
        "gene_dossier.section_1e.orthodb.fetch_release_id", fake_release
    )
    monkeypatch.setattr(
        "gene_dossier.section_1e.orthodb.fetch_gene_search", fake_search
    )


def test_node_skip_capture_partial_fallback(monkeypatch, tmp_path: Path):
    from gene_dossier.config import Settings

    settings = Settings(
        raw_data_dir=tmp_path / "raw",
        output_dir=tmp_path / "outputs",
    )
    Path(settings.raw_data_path).mkdir(parents=True, exist_ok=True)
    Path(settings.output_path).mkdir(parents=True, exist_ok=True)
    _mock_ortholog_network(monkeypatch)

    state = {
        "dossier_run_id": "run-1e-node",
        "gene_symbol": "SREBF2",
        "run_type": "section_bundle",
        "selected_section_keys": ["1e"],
        "gene_ids": {"entrez_gene_id": "6721"},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }
    out = node_generate_section_1e_derived_artifacts(
        state,
        settings=settings,
        persist_db=False,
        config=Section1eConfig(ortholog_scope_tax_id=7776),
        skip_table_capture=True,
    )
    status = out["section_1e_status"]
    summary = status["summary"]
    assert summary["scoped_ortholog_gene_count"] == 1
    assert summary["table_status"] == "complete_api_fallback"
    assert summary["ncbi_url"] == ortholog_ncbi_url(entrez_gene_id="6721")
    assert "HomoloGene" not in summary["narrative"]
    assert summary["orthodb_url"].endswith("6721")
    facts = {rec.fact_type for rec in out["evidence_records"]}
    assert "ortholog_collection_summary" in facts


def test_mock_capture_fixture_png_cannot_become_cli_official(monkeypatch, tmp_path: Path):
    """Fixture PNG via mock capture only — never CLI-supplied official_capture."""
    from gene_dossier.config import Settings
    from gene_dossier.models import ApiRun

    settings = Settings(
        raw_data_dir=tmp_path / "raw",
        output_dir=tmp_path / "outputs",
    )
    Path(settings.raw_data_path).mkdir(parents=True, exist_ok=True)
    Path(settings.output_path).mkdir(parents=True, exist_ok=True)
    _mock_ortholog_network(monkeypatch)

    fixture = tmp_path / "gold.png"
    try:
        from PIL import Image

        Image.new("RGB", (700, 400), color=(200, 200, 200)).save(fixture)
    except Exception:
        pytest.skip("Pillow required")

    def fake_capture(**kwargs):
        api = ApiRun(
            dossier_run_id=kwargs["dossier_run_id"],
            gene_symbol=kwargs["gene_symbol"],
            source_name="NCBI Gene",
            endpoint_name="capture_ncbi_ortholog_table",
            request_url=kwargs["page_url"],
            request_params={},
            success=True,
        )
        content = fixture.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        audit = {
            "status": "success",
            "capture_origin": "provided",
            "capture_metadata": {
                "capture_origin": "provided",
                "sha256": digest,
                "entrez_gene_id": kwargs["entrez_gene_id"],
                "selected_scope_tax_id": kwargs["scope_tax_id"],
                "visible_headers": list(REQUIRED_HEADERS),
                "displayed_gene_count": 2,
                "expected_displayed_count": 2,
            },
            "displayed_gene_count": 2,
            "human_reference_row_detected": True,
        }
        return api, None, None, audit, content

    monkeypatch.setattr(
        "gene_dossier.section_1e._capture_ncbi_ortholog_table", fake_capture
    )

    state = {
        "dossier_run_id": "run-1e-mock-gold",
        "gene_symbol": "SREBF2",
        "run_type": "section_bundle",
        "selected_section_keys": ["1e"],
        "gene_ids": {"entrez_gene_id": "6721"},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }
    out = node_generate_section_1e_derived_artifacts(
        state,
        settings=settings,
        persist_db=False,
        config=Section1eConfig(ortholog_scope_tax_id=7776),
        skip_table_capture=False,
    )
    summary = out["section_1e_status"]["summary"]
    assert summary["table_status"] != "official_capture"
    assert summary["table_status"] == "complete_api_fallback"


def test_incomplete_taxonomy_keeps_live_official_capture(monkeypatch, tmp_path: Path):
    """Incomplete taxonomy must not force API fallback when live capture succeeded."""
    from gene_dossier.config import Settings
    from gene_dossier.models import ApiRun, EvidenceRecord, AssertionType, EvidenceGrade, SourceType

    settings = Settings(
        raw_data_dir=tmp_path / "raw",
        output_dir=tmp_path / "outputs",
    )
    Path(settings.raw_data_path).mkdir(parents=True, exist_ok=True)
    Path(settings.output_path).mkdir(parents=True, exist_ok=True)
    _mock_ortholog_network(monkeypatch)

    # Override taxonomy to report incomplete retrieval (failed batches / unknowns).
    def incomplete_tax(tax_ids, *, scope_tax_id, gene_symbol="", settings=None):
        tids = [str(tid) for tid in tax_ids]
        # Only resolve human; leave mouse unknown → undercounted scoped set.
        membership = {tid: (True if tid == "9606" else None) for tid in tids}
        return (
            membership,
            [],
            {
                "requested_tax_ids": tids,
                "requested_count": len(tids),
                "resolved_count": sum(1 for v in membership.values() if v is not None),
                "unknown_count": sum(1 for v in membership.values() if v is None),
                "unresolved_count": sum(1 for v in membership.values() if v is None),
                "failed_request_count": 6,
                "taxonomy_complete": False,
                "scope_tax_id": scope_tax_id,
            },
        )

    monkeypatch.setattr(
        "gene_dossier.section_1e.ncbi_taxonomy.resolve_taxonomy_memberships",
        incomplete_tax,
    )

    fixture = tmp_path / "live.png"
    try:
        from PIL import Image

        Image.new("RGB", (700, 400), color=(210, 210, 210)).save(fixture)
    except Exception:
        pytest.skip("Pillow required")

    def fake_live_capture(**kwargs):
        content = fixture.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        api = ApiRun(
            dossier_run_id=kwargs["dossier_run_id"],
            gene_symbol=kwargs["gene_symbol"],
            source_name="NCBI Gene",
            endpoint_name="capture_ncbi_ortholog_table",
            request_url=kwargs["page_url"],
            request_params={},
            success=True,
        )
        # Simulate undercounted scoped expectation vs live displayed count.
        displayed = 847
        expected = 1  # would mismatch if hard-gated
        assert kwargs.get("count_claim_ready") is False
        meta = {
            "relative_path": "section_1e/live.png",
            "media_type": "image/png",
            "width": 700,
            "height": 400,
            "byte_size": len(content),
        }
        value = {
            "status": "success",
            "capture_origin": CAPTURE_ORIGIN_LIVE,
            "sha256": digest,
            "entrez_gene_id": kwargs["entrez_gene_id"],
            "selected_scope_tax_id": kwargs["scope_tax_id"],
            "visible_headers": list(REQUIRED_HEADERS) + ["Feedback"],
            "displayed_gene_count": displayed,
            "expected_displayed_count": expected,
            "relative_path": meta["relative_path"],
            "media_type": meta["media_type"],
            "width": meta["width"],
            "height": meta["height"],
            "byte_size": meta["byte_size"],
            "presentation_item_key": "orthologs-srebf2",
        }
        rec = EvidenceRecord(
            source_id="src-1e-live-cap",
            dossier_run_id=kwargs["dossier_run_id"],
            gene_symbol=kwargs["gene_symbol"],
            official_symbol=kwargs["gene_symbol"],
            section="Homologues",
            subsection="Homologues in model animals",
            source_name="NCBI Gene",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ortholog_table_capture",
            species="human",
            taxon_id=9606,
            evidence_grade=EvidenceGrade.B,
            value=value,
            display_text="live capture",
            raw_artifact_id="raw-live",
            api_run_id=api.id,
        )
        audit = {
            "status": "success",
            "capture_origin": CAPTURE_ORIGIN_LIVE,
            "capture_metadata": {
                "capture_origin": CAPTURE_ORIGIN_LIVE,
                "sha256": digest,
                "entrez_gene_id": kwargs["entrez_gene_id"],
                "selected_scope_tax_id": kwargs["scope_tax_id"],
                "visible_headers": list(REQUIRED_HEADERS) + ["Feedback"],
                "displayed_gene_count": displayed,
                "expected_displayed_count": expected,
            },
            "displayed_gene_count": displayed,
            "human_reference_row_detected": True,
            "count_claim_ready": False,
        }
        return api, meta, rec, audit, content

    monkeypatch.setattr(
        "gene_dossier.section_1e._capture_ncbi_ortholog_table", fake_live_capture
    )

    state = {
        "dossier_run_id": "run-1e-incomplete-tax",
        "gene_symbol": "SREBF2",
        "run_type": "section_bundle",
        "selected_section_keys": ["1e"],
        "gene_ids": {"entrez_gene_id": "6721"},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }
    out = node_generate_section_1e_derived_artifacts(
        state,
        settings=settings,
        persist_db=False,
        config=Section1eConfig(ortholog_scope_tax_id=7776),
        skip_table_capture=False,
    )
    summary = out["section_1e_status"]["summary"]
    assert summary["taxonomy_complete"] is False
    assert summary["count_claim_ready"] is False
    assert summary["table_status"] == "official_capture"
    assert summary["capture_origin"] == CAPTURE_ORIGIN_LIVE
    assert summary["scoped_ortholog_gene_count"] is None
    facts = {rec.fact_type for rec in out["evidence_records"]}
    assert "ortholog_table_capture" in facts


def test_taxonomy_retries_failed_batches(monkeypatch):
    """Failed taxonomy lookups are retried; audit tracks retries and completeness."""
    from gene_dossier.models import ToolResult
    from gene_dossier.tools import ncbi_taxonomy

    calls: list[list[str]] = []

    def flaky_fetch(tax_ids, **kwargs):
        batch = [str(t) for t in tax_ids]
        calls.append(batch)
        # Fail the first two attempts for any batch containing 10090.
        fail_count = sum(1 for b in calls if "10090" in b)
        if "10090" in batch and fail_count <= 2:
            return ToolResult(
                source_name="NCBI Datasets",
                endpoint_name="taxonomy_taxon",
                success=False,
                gene_symbol="SREBF2",
                request_url="https://example.test/taxonomy",
                request_params={"taxons": ",".join(batch), "retry_count": 1},
                status_code=503,
                error_type="http_error",
                error_message="HTTP 503",
            )
        nodes = []
        for tid in batch:
            nodes.append(
                {
                    "taxonomy": {
                        "tax_id": int(tid),
                        "lineage": [2759, 7776, int(tid)],
                    }
                }
            )
        return ToolResult(
            source_name="NCBI Datasets",
            endpoint_name="taxonomy_taxon",
            success=True,
            gene_symbol="SREBF2",
            request_url="https://example.test/taxonomy",
            request_params={"taxons": ",".join(batch)},
            status_code=200,
            data={"taxonomy_nodes": nodes},
        )

    monkeypatch.setattr(ncbi_taxonomy, "fetch_taxonomy_taxon", flaky_fetch)
    monkeypatch.setattr(ncbi_taxonomy.time, "sleep", lambda *_a, **_k: None)

    membership, results, audit = ncbi_taxonomy.resolve_taxonomy_memberships(
        ["9606", "10090"],
        scope_tax_id=7776,
        gene_symbol="SREBF2",
        batch_size=20,
        max_concurrency=1,
        max_attempts=1,
        retry_sleep_seconds=0.0,
        resolve_rounds=3,
    )
    assert membership["9606"] is True
    assert membership["10090"] is True
    assert audit["taxonomy_complete"] is True
    assert audit["resolved_count"] == 2
    assert audit["requested_count"] == 2
    assert audit["failed_request_count"] == 0
    assert audit["resolve_rounds_run"] >= 2
    assert len(calls) >= 2
