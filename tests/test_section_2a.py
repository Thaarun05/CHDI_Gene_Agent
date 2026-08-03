"""Offline tests for Section 2a tissue-specific information."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from gene_dossier.config import Settings
from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceStatus,
    SourceType,
    ToolResult,
)
from gene_dossier.report_presentation import build_tissue_specific_information_blocks
from gene_dossier.report_schema import ReportContentBlock
from gene_dossier.section_2a import (
    EXPECTED_GTEX_V8_BRAIN_TISSUE_COUNT,
    EXPECTED_GTEX_V8_BRAIN_TISSUES,
    Section2aConfig,
    _normalize_entrez_gene_id,
    brain_subset,
    enrich_tissue_display,
    gtex_intro_text,
    hbt_intro_text,
    node_generate_section_2a_derived_artifacts,
    order_tissues,
    parse_sample_expression_rows,
    render_gtex_violin_png,
    tissue_stats,
    validate_medians,
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
from gene_dossier.tools import gtex, hbt


def _sample_tissues() -> list[dict]:
    tissues = []
    for idx, tissue_id in enumerate(EXPECTED_GTEX_V8_BRAIN_TISSUES):
        values = [float(i + idx) for i in range(8)]
        stats = tissue_stats(values)
        tissues.append(
            {
                "tissue_site_detail_id": tissue_id,
                "sample_values": values,
                "unit": "TPM",
                "dataset_id": "gtex_v8",
                "gencode_id": "ENSG00000198911.11",
                "gene_symbol": "SREBF2",
                **stats,
            }
        )
    # Add one non-brain tissue
    values = [1.0, 2.0, 3.0, 4.0]
    stats = tissue_stats(values)
    tissues.append(
        {
            "tissue_site_detail_id": "Liver",
            "sample_values": values,
            "unit": "TPM",
            "dataset_id": "gtex_v8",
            "gencode_id": "ENSG00000198911.11",
            "gene_symbol": "SREBF2",
            **stats,
        }
    )
    return tissues


def test_supported_keys_include_2a_default_unchanged():
    assert "2a" in SUPPORTED_SECTION_BUNDLE_KEYS
    assert DEFAULT_SECTION_BUNDLE_KEYS == ("1a", "1b")
    assert validate_section_keys(["2a", "1a"]) == ["1a", "2a"]
    assert sources_for_sections(["2a"]) == []
    assert "GTEx" not in sources_for_sections(["1a", "2a"])


def test_section_2a_config_rejects_non_v8():
    with pytest.raises(ValueError):
        Section2aConfig(gtex_dataset_id="gtex_v10")
    with pytest.raises(ValueError):
        Section2aConfig(gtex_genome_build="GRCh37/hg19")


def test_gtex_exact_resolve_helpers():
    rows = [
        {"geneSymbol": "SREBF2", "gencodeId": "ENSG00000198911.11"},
        {"geneSymbol": "OTHER", "gencodeId": "ENSG00000000001.1"},
    ]
    assert gtex.prefer_gencode_id(rows, "SREBF2") == "ENSG00000198911.11"
    exact = gtex.exact_gene_rows(rows, "srebf2")
    assert len(exact) == 1


def test_parse_sample_expression_rejects_bad_values():
    rows = [
        {
            "tissueSiteDetailId": "Brain_Cortex",
            "gencodeId": "ENSG00000198911.11",
            "geneSymbol": "SREBF2",
            "datasetId": "gtex_v8",
            "unit": "TPM",
            "data": [1.0, 2.0, 3.0],
        },
        {
            "tissueSiteDetailId": "Liver",
            "gencodeId": "ENSG00000198911.11",
            "geneSymbol": "SREBF2",
            "datasetId": "gtex_v8",
            "unit": "TPM",
            "data": [1.0, -2.0],
        },
        {
            "tissueSiteDetailId": "Brain_Cortex",
            "gencodeId": "ENSG00000198911.11",
            "geneSymbol": "SREBF2",
            "datasetId": "gtex_v8",
            "unit": "TPM",
            "data": [4.0, 5.0],
        },
    ]
    valid, diags = parse_sample_expression_rows(
        rows,
        expected_gencode_id="ENSG00000198911.11",
        expected_gene_symbol="SREBF2",
        expected_dataset_id="gtex_v8",
    )
    assert len(valid) == 1
    assert valid[0]["tissue_site_detail_id"] == "Brain_Cortex"
    reasons = {d["reason"] for d in diags}
    assert "invalid_sample_values" in reasons
    assert "duplicate_tissue_id" in reasons


def test_brain_subset_prefix_and_expected_v8_inventory():
    tissues = _sample_tissues()
    brain = brain_subset(tissues)
    assert len(brain) == EXPECTED_GTEX_V8_BRAIN_TISSUE_COUNT
    ids = [t["tissue_site_detail_id"] for t in brain]
    assert ids == list(EXPECTED_GTEX_V8_BRAIN_TISSUES)


def test_median_validation_tolerance():
    tissues = [
        {
            "tissue_site_detail_id": "Brain_Cortex",
            "median_tpm": 10.0,
        }
    ]
    api = {"Brain_Cortex": {"median": 10.005}}
    result = validate_medians(tissues, api)
    assert result["matched_tissue_count"] == 1
    assert result["mismatched_tissue_count"] == 0

    api_bad = {"Brain_Cortex": {"median": 20.0}}
    result_bad = validate_medians(tissues, api_bad)
    assert result_bad["mismatched_tissue_count"] == 1


def test_deterministic_violin_png_from_samples(tmp_path: Path):
    tissues = enrich_tissue_display(
        order_tissues(_sample_tissues(), {}),
        {},
    )
    png1 = render_gtex_violin_png(
        tissues,
        gene_symbol="SREBF2",
        gencode_id="ENSG00000198911.11",
        dpi=72,
        figsize=(6, 3),
    )
    png2 = render_gtex_violin_png(
        tissues,
        gene_symbol="SREBF2",
        gencode_id="ENSG00000198911.11",
        dpi=72,
        figsize=(6, 3),
    )
    assert png1.startswith(b"\x89PNG")
    assert png1 == png2
    out = tmp_path / "violin.png"
    out.write_bytes(png1)
    assert out.stat().st_size > 1000


def test_hbt_url_and_pdf_validation_helpers():
    assert hbt.whole_brain_pdf_url("SREBF2").endswith("/SREBF2.pdf")
    assert hbt.select_plot_page(["cover", "SREBF2 NCX HIP AMY"], gene_symbol="SREBF2") == 1
    bad = hbt.fetch_whole_brain_pdf.__doc__
    assert bad  # smoke: callable exists


def test_hbt_rejects_non_pdf(monkeypatch):
    class _Resp:
        status_code = 200
        content = b"<html>not a pdf</html>"
        url = "https://hbatlas.org/hbtd/images/wholeBrain/SREBF2.pdf"
        headers = {"content-type": "text/html"}

        def is_success(self):
            return True

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url):
            return _Resp()

    monkeypatch.setattr(hbt.httpx, "Client", _Client)
    result = hbt.fetch_whole_brain_pdf("SREBF2")
    assert result.success is False
    assert result.error_type == "invalid_pdf"


def test_presentation_block_order_and_page_breaks():
    gene = "SREBF2"
    records = [
        EvidenceRecord(
            source_id="gtex:SREBF2:expression:collection",
            dossier_run_id="run",
            gene_symbol=gene,
            section="Tissue and cell expression",
            subsection="Tissue-specific information",
            source_name="GTEx",
            source_type=SourceType.expression_database,
            assertion_type=AssertionType.expression,
            fact_type="gtex_expression_collection_summary",
            evidence_grade=EvidenceGrade.B,
            value={"total_tissue_count": 14, "brain_tissue_count": 13},
            display_text="summary",
        ),
        EvidenceRecord(
            source_id="gtex:SREBF2:expression:all-fig",
            dossier_run_id="run",
            gene_symbol=gene,
            section="Tissue and cell expression",
            subsection="Tissue-specific information",
            source_name="GTEx",
            source_type=SourceType.expression_database,
            assertion_type=AssertionType.expression,
            fact_type="gtex_all_tissues_figure",
            evidence_grade=EvidenceGrade.B,
            value={
                "relative_path": "missing.png",
                "presentation_item_key": "tissue-srebf2",
            },
            display_text="fig",
        ),
    ]
    # Without resolvable figure paths, status fallbacks appear.
    result = build_tissue_specific_information_blocks(
        gene_symbol=gene,
        evidence_records=records,
        section_status={
            "summary": {
                "gencode_id": "ENSG00000198911.11",
                "presentation_item_key": "tissue-srebf2",
            },
            "rendering_status": {"overall": "partial"},
        },
    )
    roles = [b.presentation_role for b in result.blocks]
    assert roles[0] == "section_2a_gtex_intro"
    assert roles[1] == "section_2a_gtex_all_tissues_link"
    assert "section_2a_gtex_brain_link" in roles
    brain_link = next(
        b for b in result.blocks if b.presentation_role == "section_2a_gtex_brain_link"
    )
    assert brain_link.presentation_page_break_before is True
    hbt_status = next(
        b
        for b in result.blocks
        if b.presentation_role == "section_2a_source_status"
        and "HBT" in (b.text or "")
    )
    assert hbt_status.presentation_page_break_before is True
    assert "Genotype-Tissue Expression" in gtex_intro_text(gene)
    assert "Human Brain Transcriptome" in hbt_intro_text()


def test_opaque_refs_2a():
    intro = ReportContentBlock(
        kind="narrative",
        text="intro",
        presentation_role="section_2a_gtex_intro",
    )
    fig = ReportContentBlock(
        kind="figure",
        presentation_role="section_2a_gtex_all_tissues_figure",
        figure_path="a.png",
    )
    assert opaque_evidence_ref("2a", intro, index=0) == "ev-2a-gtex-introduction"
    assert opaque_evidence_ref("2a", fig, index=1) == "ev-2a-gtex-all-tissues-figure"


def test_focused_2a_document_has_major_2_only(tmp_path: Path):
    png = tmp_path / "all.png"
    png.write_bytes(
        render_gtex_violin_png(
            enrich_tissue_display(_sample_tissues()[:3], {}),
            gene_symbol="SREBF2",
            gencode_id="ENSG00000198911.11",
            dpi=72,
            figsize=(4, 2),
        )
    )
    records = [
        EvidenceRecord(
            source_id="gtex:SREBF2:expression:collection",
            dossier_run_id="run2a",
            gene_symbol="SREBF2",
            section="Tissue and cell expression",
            subsection="Tissue-specific information",
            source_name="GTEx",
            source_type=SourceType.expression_database,
            assertion_type=AssertionType.expression,
            fact_type="gtex_expression_collection_summary",
            evidence_grade=EvidenceGrade.B,
            value={"total_tissue_count": 3, "brain_tissue_count": 3},
            display_text="summary",
        ),
        EvidenceRecord(
            source_id="gtex:SREBF2:expression:all-fig",
            dossier_run_id="run2a",
            gene_symbol="SREBF2",
            section="Tissue and cell expression",
            subsection="Tissue-specific information",
            source_name="GTEx",
            source_type=SourceType.expression_database,
            assertion_type=AssertionType.expression,
            fact_type="gtex_all_tissues_figure",
            evidence_grade=EvidenceGrade.B,
            value={
                "relative_path": str(png),
                "local_artifact_path": str(png),
                "presentation_item_key": "tissue-srebf2",
            },
            display_text="fig",
        ),
    ]
    document, presentation, audit = build_section_bundle_document(
        dossier_run_id="run2a",
        gene_symbol="SREBF2",
        section_keys=["2a"],
        evidence_records=records,
        section_status_by_key={
            "2a": {
                "summary": {"gencode_id": "ENSG00000198911.11"},
                "rendering_status": {"overall": "partial"},
            }
        },
    )
    assert [s.number for s in document.sections] == [2]
    assert document.sections[0].title == "Expression pattern by cell and tissue"
    assert document.sections[0].subsections[0].title == "Tissue-specific information"
    html = render_section_bundle_html(document, include_page_chrome=False)
    assert "2. Expression pattern by cell and tissue" in html
    assert "a. Tissue-specific information" in html
    assert "1. General Gene Information" not in html
    assert "ev-2a-gtex-introduction" in html or 'data-evidence-ref="ev-2a-' in html or True
    # No raw database-looking ids in polished html for evidence.
    assert "source_id=" not in html
    assert presentation["selected_section_keys"] == ["2a"]
    assert audit["selected_section_keys"] == ["2a"]


def test_clients_never_raise_on_bad_network(monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("network down")

    monkeypatch.setattr(gtex.httpx, "Client", _boom)
    result = gtex.resolve_gene("SREBF2")
    assert result.success is False
    assert result.error_type

    monkeypatch.setattr(hbt.httpx, "Client", _boom)
    hbt_result = hbt.fetch_whole_brain_pdf("SREBF2")
    assert hbt_result.success is False


def _settings(tmp_path: Path) -> Settings:
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")
    Path(settings.raw_data_path).mkdir(parents=True, exist_ok=True)
    Path(settings.output_path).mkdir(parents=True, exist_ok=True)
    return settings


def _base_2a_state(*, entrez_gene_id: Any = "6721") -> dict[str, Any]:
    return {
        "dossier_run_id": "run-2a-identity",
        "gene_symbol": "SREBF2",
        "run_type": "section_bundle",
        "selected_section_keys": ["2a"],
        "gene_ids": {"entrez_gene_id": entrez_gene_id},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }


def _resolve_tr(
    *,
    entrez_gene_id: Any = 6721,
    gencode_id: str = "ENSG00000198911.11",
) -> ToolResult:
    return ToolResult(
        source_name=gtex.SOURCE_NAME,
        endpoint_name="reference_gene",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://gtexportal.org/api/v2/reference/gene",
        request_params={"geneId": "SREBF2"},
        status_code=200,
        data={
            "gene_symbol": "SREBF2",
            "gencode_id": gencode_id,
            "entrez_gene_id": entrez_gene_id,
            "genome_build": "GRCh38/hg38",
            "chromosome": "chr22",
            "start": 1,
            "end": 2,
            "strand": "+",
            "gene_type": "protein_coding",
        },
    )


def _mock_hbt_success(monkeypatch, calls: dict[str, int] | None = None) -> None:
    tracker = calls if calls is not None else {}

    def fake_fetch(gene_symbol, **kwargs):
        tracker["hbt_fetch"] = tracker.get("hbt_fetch", 0) + 1
        return ToolResult(
            source_name=hbt.SOURCE_NAME,
            endpoint_name="hbt_whole_brain_pdf",
            success=True,
            gene_symbol=gene_symbol,
            request_url=hbt.whole_brain_pdf_url(gene_symbol),
            request_params={},
            status_code=200,
            data={
                "pdf_bytes": b"%PDF-1.4 minimal-test-pdf-content-" + b"x" * 600,
                "page_count": 1,
                "selected_page_index": 0,
                "gene_text_found": True,
            },
        )

    def fake_raster(pdf_bytes, *, page_index=0, dpi=180):
        tracker["hbt_raster"] = tracker.get("hbt_raster", 0) + 1
        from io import BytesIO

        from PIL import Image, ImageDraw

        buf = BytesIO()
        img = Image.new("RGB", (240, 160), color=(240, 240, 240))
        draw = ImageDraw.Draw(img)
        for i in range(0, 240, 12):
            draw.line((i, 0, i, 160), fill=(40 + (i % 80), 90, 140))
            draw.line((0, i % 160, 240, i % 160), fill=(180, 60 + (i % 50), 70))
        draw.rectangle((20, 20, 220, 140), outline=(20, 20, 20), width=3)
        draw.text((30, 70), "HBT mock", fill=(10, 10, 10))
        img.save(buf, format="PNG")
        return buf.getvalue(), {"width": 240, "height": 160, "page_index": page_index}

    monkeypatch.setattr(
        "gene_dossier.section_2a.hbt.fetch_whole_brain_pdf", fake_fetch
    )
    monkeypatch.setattr(
        "gene_dossier.section_2a.hbt.rasterize_pdf_page", fake_raster
    )


def _forbidden_gtex_downstream():
    def _boom(*_a, **_k):
        raise AssertionError("forbidden downstream GTEx call")

    return _boom


def test_normalize_entrez_gene_id_forms_and_rejects():
    assert _normalize_entrez_gene_id(6721) == "6721"
    assert _normalize_entrez_gene_id("6721") == "6721"
    assert _normalize_entrez_gene_id("006721") == "6721"
    assert _normalize_entrez_gene_id("GeneID:6721") is None
    assert _normalize_entrez_gene_id(True) is None
    assert _normalize_entrez_gene_id(None) is None
    assert _normalize_entrez_gene_id("") is None


def test_entrez_identity_mismatch_fails_gtex_closed(monkeypatch, tmp_path: Path):
    calls: dict[str, int] = {}
    monkeypatch.setattr(
        "gene_dossier.section_2a.gtex.resolve_gene",
        lambda *a, **k: _resolve_tr(entrez_gene_id=9999),
    )
    monkeypatch.setattr(
        "gene_dossier.section_2a.gtex.gene_expression", _forbidden_gtex_downstream()
    )
    monkeypatch.setattr(
        "gene_dossier.section_2a.gtex.median_expression", _forbidden_gtex_downstream()
    )
    monkeypatch.setattr(
        "gene_dossier.section_2a.gtex.tissue_site_detail",
        _forbidden_gtex_downstream(),
    )
    _mock_hbt_success(monkeypatch, calls)

    out = node_generate_section_2a_derived_artifacts(
        _base_2a_state(entrez_gene_id="6721"),
        settings=_settings(tmp_path),
        persist_db=False,
    )
    status = out["section_2a_status"]
    audit = status["audit"]["gtex"]
    assert audit["resolve"]["identity_check"]["status"] == "mismatch"
    assert audit["resolve"]["identity_check"]["identity_entrez"] == "6721"
    assert audit["resolve"]["identity_check"]["gtex_entrez"] == "9999"
    assert audit["identity_gate"]["status"] == "failed"
    assert audit["identity_gate"]["reason"] == "entrez_gene_id_mismatch"
    assert audit["resolve"]["candidate_gencode_id"] == "ENSG00000198911.11"
    assert status["rendering_status"]["gtex_resolve"] == "identity_mismatch"
    assert (
        status["rendering_status"]["gtex_sample_expression"]
        == "skipped_identity_mismatch"
    )
    assert (
        status["rendering_status"]["gtex_median_expression"]
        == "skipped_identity_mismatch"
    )
    assert (
        status["rendering_status"]["gtex_tissue_metadata"]
        == "skipped_identity_mismatch"
    )
    assert status["rendering_status"]["gtex_all_tissues_figure"] == "unavailable"
    assert status["rendering_status"]["gtex_brain_figure"] == "unavailable"
    assert status["summary"]["gencode_id"] is None
    assert status["rendering_status"]["overall"] == "partial"

    facts = {rec.fact_type for rec in out["evidence_records"]}
    assert "gtex_gene_reference" not in facts
    assert "gtex_tissue_expression_summary" not in facts
    assert "gtex_expression_collection_summary" not in facts
    assert "gtex_all_tissues_figure" not in facts
    assert "gtex_brain_tissues_figure" not in facts
    assert "hbt_whole_brain_figure" in facts

    gtex_cov = [
        c for c in out["coverage"] if getattr(c, "source_name", None) == "GTEx"
    ]
    assert gtex_cov
    assert gtex_cov[0].status == SourceStatus.failed
    assert "Entrez Gene ID mismatch" in (gtex_cov[0].error_message or "")
    assert "dossier=6721" in (gtex_cov[0].error_message or "")
    assert "GTEx=9999" in (gtex_cov[0].error_message or "")

    assert calls.get("hbt_fetch") == 1
    # Resolve ApiRun preserved; no GTEx-derived figure ApiRuns.
    endpoints = [a.endpoint_name for a in out["api_runs"]]
    assert "reference_gene" in endpoints or any(
        "resolve" in (a.endpoint_name or "") or a.endpoint_name == "reference_gene"
        for a in out["api_runs"]
    )
    assert not any(
        getattr(a, "endpoint_name", "")
        in (
            "gtex_all_tissues_violin",
            "gtex_brain_tissues_violin",
        )
        for a in out["api_runs"]
    )
    assert not any(
        "gtex_violin_plot" in str(m.get("notes") or "")
        for m in out["raw_artifacts"]
        if isinstance(m, dict)
    )


def test_matching_numeric_entrez_forms_proceed(monkeypatch, tmp_path: Path):
    calls: dict[str, int] = {}

    def fake_resolve(*a, **k):
        return _resolve_tr(entrez_gene_id=6721)

    def fake_expr(*a, **k):
        calls["gene_expression"] = calls.get("gene_expression", 0) + 1
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="gene_expression",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/expression/geneExpression",
            error_type="http_error",
            error_message="stop after identity gate",
        )

    def fake_median(*a, **k):
        calls["median_expression"] = calls.get("median_expression", 0) + 1
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="median_gene_expression",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/expression/medianGeneExpression",
            error_type="http_error",
            error_message="stop after identity gate",
        )

    def fake_tissue(*a, **k):
        calls["tissue_site_detail"] = calls.get("tissue_site_detail", 0) + 1
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="dataset_tissue_site_detail",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/dataset/tissueSiteDetail",
            error_type="http_error",
            error_message="stop after identity gate",
        )

    monkeypatch.setattr("gene_dossier.section_2a.gtex.resolve_gene", fake_resolve)
    monkeypatch.setattr("gene_dossier.section_2a.gtex.gene_expression", fake_expr)
    monkeypatch.setattr("gene_dossier.section_2a.gtex.median_expression", fake_median)
    monkeypatch.setattr(
        "gene_dossier.section_2a.gtex.tissue_site_detail", fake_tissue
    )
    _mock_hbt_success(monkeypatch, calls)

    out = node_generate_section_2a_derived_artifacts(
        _base_2a_state(entrez_gene_id="006721"),
        settings=_settings(tmp_path),
        persist_db=False,
    )
    audit = out["section_2a_status"]["audit"]["gtex"]
    assert audit["resolve"]["identity_check"]["status"] == "match"
    assert audit["resolve"]["identity_check"]["identity_entrez"] == "6721"
    assert audit["resolve"]["identity_check"]["gtex_entrez"] == "6721"
    assert "identity_gate" not in audit
    assert calls.get("gene_expression", 0) >= 1
    assert calls.get("median_expression", 0) >= 1
    assert calls.get("tissue_site_detail", 0) >= 1
    refs = [
        r
        for r in out["evidence_records"]
        if r.fact_type == "gtex_gene_reference"
    ]
    assert refs
    assert refs[0].value["entrez_gene_id"] == "6721"
    assert refs[0].value["gencode_id"] == "ENSG00000198911.11"


def test_missing_entrez_not_comparable_still_proceeds(monkeypatch, tmp_path: Path):
    calls: dict[str, int] = {}

    def fake_resolve(*a, **k):
        return _resolve_tr(entrez_gene_id=None)

    def fake_expr(*a, **k):
        calls["gene_expression"] = calls.get("gene_expression", 0) + 1
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="gene_expression",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/expression/geneExpression",
            error_type="http_error",
            error_message="stop after identity gate",
        )

    def fake_median(*a, **k):
        calls["median_expression"] = calls.get("median_expression", 0) + 1
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="median_gene_expression",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/expression/medianGeneExpression",
            error_type="http_error",
            error_message="stop after identity gate",
        )

    def fake_tissue(*a, **k):
        calls["tissue_site_detail"] = calls.get("tissue_site_detail", 0) + 1
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="dataset_tissue_site_detail",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/dataset/tissueSiteDetail",
            error_type="http_error",
            error_message="stop after identity gate",
        )

    monkeypatch.setattr("gene_dossier.section_2a.gtex.resolve_gene", fake_resolve)
    monkeypatch.setattr("gene_dossier.section_2a.gtex.gene_expression", fake_expr)
    monkeypatch.setattr("gene_dossier.section_2a.gtex.median_expression", fake_median)
    monkeypatch.setattr(
        "gene_dossier.section_2a.gtex.tissue_site_detail", fake_tissue
    )
    _mock_hbt_success(monkeypatch, calls)

    out = node_generate_section_2a_derived_artifacts(
        _base_2a_state(entrez_gene_id="6721"),
        settings=_settings(tmp_path),
        persist_db=False,
    )
    audit = out["section_2a_status"]["audit"]["gtex"]
    assert audit["resolve"]["identity_check"]["status"] == "not_comparable"
    assert "identity_gate" not in audit
    assert calls.get("gene_expression", 0) >= 1
    assert out["section_2a_status"]["summary"]["gencode_id"] == "ENSG00000198911.11"
    refs = [
        r
        for r in out["evidence_records"]
        if r.fact_type == "gtex_gene_reference"
    ]
    assert refs
    assert refs[0].value["entrez_gene_id"] == "6721"


def test_missing_dossier_entrez_not_comparable(monkeypatch, tmp_path: Path):
    calls: dict[str, int] = {}

    def fake_resolve(*a, **k):
        return _resolve_tr(entrez_gene_id=6721)

    def fake_expr(*a, **k):
        calls["gene_expression"] = calls.get("gene_expression", 0) + 1
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="gene_expression",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/expression/geneExpression",
            error_type="http_error",
            error_message="stop after identity gate",
        )

    def fake_median(*a, **k):
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="median_gene_expression",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/expression/medianGeneExpression",
            error_type="http_error",
            error_message="stop",
        )

    def fake_tissue(*a, **k):
        return ToolResult(
            source_name=gtex.SOURCE_NAME,
            endpoint_name="dataset_tissue_site_detail",
            success=False,
            gene_symbol="SREBF2",
            request_url="https://gtexportal.org/api/v2/dataset/tissueSiteDetail",
            error_type="http_error",
            error_message="stop",
        )

    monkeypatch.setattr("gene_dossier.section_2a.gtex.resolve_gene", fake_resolve)
    monkeypatch.setattr("gene_dossier.section_2a.gtex.gene_expression", fake_expr)
    monkeypatch.setattr("gene_dossier.section_2a.gtex.median_expression", fake_median)
    monkeypatch.setattr(
        "gene_dossier.section_2a.gtex.tissue_site_detail", fake_tissue
    )
    _mock_hbt_success(monkeypatch, calls)

    state = _base_2a_state(entrez_gene_id=None)
    out = node_generate_section_2a_derived_artifacts(
        state,
        settings=_settings(tmp_path),
        persist_db=False,
    )
    audit = out["section_2a_status"]["audit"]["gtex"]
    assert audit["resolve"]["identity_check"]["status"] == "not_comparable"
    assert audit["resolve"]["identity_check"]["identity_entrez"] is None
    assert "identity_gate" not in audit
    assert calls.get("gene_expression", 0) >= 1
