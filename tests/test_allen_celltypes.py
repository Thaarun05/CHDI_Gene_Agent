"""Offline tests for the Allen Cell Types (10x) client.

No network. Fixtures are synthetic and gene-general (``GENEX`` / ``Genex``).
"""

from __future__ import annotations

import base64
import io
import json
import sys
import types
from pathlib import Path

import pytest

from gene_dossier.tools import allen_celltypes as ac

# 1x1 PNG so paired-scatter compositing (PIL) works in offline browser fakes.
_MIN_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg=="
)

SRC = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "gene_dossier"
    / "tools"
    / "allen_celltypes.py"
)

HUMAN_CSV = """feature,Astro_1,L2/3 IT_1,Pvalb_1,Vip_1
GENEX,0.0,4.5,1.5,0.5
GENEX-AS1,9.9,9.9,9.9,9.9
GENEXP1,8.8,8.8,8.8,8.8
OTHER,1.0,1.0,1.0,1.0
"""

MOUSE_CSV = """feature,1_CR,2_Lamp5,3_Pvalb
Genex,0.0,2.0,6.0
Other,1.0,1.0,1.0
"""


def leaf(alias: str, *, structure: str = "Cortical plate") -> dict:
    return {
        "leaf_attributes": [
            {
                "leaf": True,
                "cell_set_alias": alias,
                "original_label": alias,
                "cell_set_accession": f"CS_{alias}",
                "cell_set_designation": f"CTX-HPF {alias}",
                "cell_set_structure": structure,
            }
        ]
    }


def internal(alias: str, children: list[dict]) -> dict:
    return {
        "node_attributes": [{"cell_set_alias": alias, "cell_set_designation": alias}],
        "children": children,
    }


DEND = internal(
    "",
    [
        internal(
            "Gabaergic neurons",
            [internal("n17", [leaf("2_Lamp5")]), leaf("3_Pvalb")],
        ),
        internal("Glutamatergic neurons", [leaf("1_CR")]),
    ],
)


# --------------------------------------------------------------------------------------
# Source configuration and acquisition
# --------------------------------------------------------------------------------------
def test_human_m1_source_urls_are_exact_and_allowlisted():
    assert ac.human_m1_source_url("trimmed_means").endswith(
        "aibs_human_m1_10x/trimmed_means.csv"
    )
    assert ac.human_m1_source_url("taxonomy").endswith("aibs_human_m1_10x/dend.json")
    for url in ac.HUMAN_M1_SOURCES.values():
        assert ac.is_allowlisted_download_url(url)
        assert ac.DOWNLOAD_HOST in url


def test_only_the_two_required_human_sources_are_configured():
    # metadata.csv / tsne.csv / medians.csv and the per-cell matrix are out of scope.
    assert set(ac.HUMAN_M1_SOURCES) == {"trimmed_means", "taxonomy"}
    joined = " ".join(ac.HUMAN_M1_SOURCES.values())
    for excluded in ("metadata.csv", "tsne.csv", "medians.csv", "matrix.csv"):
        assert excluded not in joined


def test_off_bucket_url_rejected():
    assert not ac.is_allowlisted_download_url("https://example.com/trimmed_means.csv")
    assert not ac.is_allowlisted_download_url(
        "http://idk-etl-prod-download-bucket.s3.amazonaws.com/x.csv"
    )


def test_unknown_source_key_fails_closed_without_network():
    tr = ac.download_human_m1_source("medians")
    assert tr.success is False
    assert tr.error_type == "invalid_request"


def test_html_payload_rejected_for_csv_and_json():
    html = b"<!DOCTYPE html><html><body>oops</body></html>"
    for key in ("trimmed_means", "taxonomy"):
        result = ac.validate_human_m1_payload(key, html)
        assert result["ok"] is False
        assert result["error_type"] == "html_masquerading_as_data"


def test_malformed_taxonomy_json_rejected():
    result = ac.validate_human_m1_payload("taxonomy", b"{not json")
    assert result["ok"] is False
    assert result["error_type"] == "malformed_json"


def test_valid_payloads_accepted():
    assert ac.validate_human_m1_payload("trimmed_means", HUMAN_CSV.encode())["ok"] is True
    assert ac.validate_human_m1_payload("taxonomy", json.dumps(DEND).encode())["ok"] is True


def test_register_local_source_records_original_filename(tmp_path):
    src = tmp_path / "response.xls"
    src.write_text(MOUSE_CSV, encoding="utf-8")
    out = ac.register_local_source(src, canonical_name=ac.MOUSE_TRIMMED_MEANS_FILENAME)
    assert out["ok"] is True
    assert out["original_filename"] == "response.xls"
    assert out["canonical_name"] == "mouse_ctx_hpf_trimmed_means.csv"
    assert out["sha256"]
    assert out["retrieval_method"] == "local_registration"


def test_register_local_source_missing_file():
    out = ac.register_local_source("/nonexistent/response.xls", canonical_name="x.csv")
    assert out["ok"] is False
    assert out["error_type"] == "local_source_missing"


# --------------------------------------------------------------------------------------
# Gene row extraction
# --------------------------------------------------------------------------------------
def test_human_exact_gene_row_extraction():
    out = ac.extract_gene_row(io.StringIO(HUMAN_CSV), "GENEX")
    assert out.ok is True
    assert out.match_count == 1
    assert out.source_symbol == "GENEX"
    assert out.celltype_labels == ["Astro_1", "L2/3 IT_1", "Pvalb_1", "Vip_1"]
    assert out.values == [0.0, 4.5, 1.5, 0.5]


def test_mouse_exact_gene_row_extraction():
    out = ac.extract_gene_row(io.StringIO(MOUSE_CSV), "Genex")
    assert out.ok is True
    assert out.values == [0.0, 2.0, 6.0]


def test_substring_paralog_rows_never_match():
    out = ac.extract_gene_row(io.StringIO(HUMAN_CSV), "GENEX")
    # GENEX-AS1 (9.9) and GENEXP1 (8.8) must not contribute.
    assert 9.9 not in out.values
    assert 8.8 not in out.values


def test_case_insensitive_match_preserves_source_symbol():
    out = ac.extract_gene_row(io.StringIO(HUMAN_CSV), "genex")
    assert out.ok is True
    assert out.source_symbol == "GENEX"


def test_duplicate_exact_rows_fail_closed():
    text = "feature,A,B\nGENEX,1,2\nGENEX,3,4\n"
    out = ac.extract_gene_row(io.StringIO(text), "GENEX")
    assert out.ok is False
    assert out.error_type == "duplicate_gene_rows"
    assert out.match_count == 2


def test_missing_gene_fails_closed():
    out = ac.extract_gene_row(io.StringIO(HUMAN_CSV), "ABSENT")
    assert out.ok is False
    assert out.error_type == "gene_not_found"


def test_empty_matrix_fails_closed():
    out = ac.extract_gene_row(io.StringIO(""), "GENEX")
    assert out.ok is False
    assert out.error_type == "empty_matrix"


def test_matrix_celltype_labels():
    assert ac.matrix_celltype_labels(io.StringIO(MOUSE_CSV)) == [
        "1_CR",
        "2_Lamp5",
        "3_Pvalb",
    ]


# --------------------------------------------------------------------------------------
# Taxonomy
# --------------------------------------------------------------------------------------
def test_parse_dendrogram_collects_leaves_and_named_ancestors():
    tax = ac.parse_dendrogram(DEND)
    assert tax.leaf_count == 3
    assert set(tax.leaves) == {"1_CR", "2_Lamp5", "3_Pvalb"}
    # Anonymous structural node n17 is skipped; named ancestors are kept.
    assert tax.leaves["2_Lamp5"]["ancestors"] == ["Gabaergic neurons"]
    assert tax.leaves["1_CR"]["ancestors"] == ["Glutamatergic neurons"]


def test_parse_dendrogram_tolerates_garbage():
    assert ac.parse_dendrogram(None).leaf_count == 0
    assert ac.parse_dendrogram({"children": []}).leaf_count == 0


def test_reconcile_taxonomy_records_missing_cluster_without_zero_filling():
    result = ac.reconcile_taxonomy(
        taxonomy_leaves=["1_CR", "2_Lamp5", "3_Pvalb", "381_SMC-Peri"],
        matrix_labels=["1_CR", "2_Lamp5", "3_Pvalb"],
    )
    assert result["taxonomy_leaf_count"] == 4
    assert result["expression_cluster_count"] == 3
    assert result["analyzed_intersection_count"] == 3
    assert result["missing_expression_clusters"] == ["381_SMC-Peri"]


def test_reconcile_reports_matrix_only_columns():
    result = ac.reconcile_taxonomy(
        taxonomy_leaves=["a"], matrix_labels=["a", "extra"]
    )
    assert result["matrix_columns_absent_from_taxonomy"] == ["extra"]


# --------------------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------------------
def human_summary():
    out = ac.extract_gene_row(io.StringIO(HUMAN_CSV), "GENEX")
    return ac.summarize_celltype_expression(
        gene_symbol="GENEX",
        source_symbol=out.source_symbol,
        celltype_labels=out.celltype_labels,
        values=out.values,
    )


def test_nonzero_naming_replaces_detected():
    summary = human_summary()
    assert summary["valid_celltype_count"] == 4
    assert summary["nonzero_celltype_count"] == 3
    assert summary["nonzero_percentage"] == pytest.approx(75.0)
    for forbidden in ("detected_count", "detected_percentage", "detected_celltype_count"):
        assert forbidden not in summary


def test_summary_statistics():
    summary = human_summary()
    assert summary["maximum"] == pytest.approx(4.5)
    assert summary["minimum"] == pytest.approx(0.0)
    assert summary["median"] == pytest.approx(1.0)
    # max 4.5 / median 1.0
    assert summary["max_to_median_ratio"] == pytest.approx(4.5)


def test_top_ranking_is_descending():
    summary = human_summary()
    labels = [item["label"] for item in summary["top"]]
    assert labels == ["L2/3 IT_1", "Pvalb_1", "Vip_1", "Astro_1"]


def test_summary_attaches_taxonomy_ancestors():
    out = ac.extract_gene_row(io.StringIO(MOUSE_CSV), "Genex")
    tax = ac.parse_dendrogram(DEND)
    summary = ac.summarize_celltype_expression(
        gene_symbol="Genex",
        source_symbol=out.source_symbol,
        celltype_labels=out.celltype_labels,
        values=out.values,
        taxonomy=tax,
        count_noun="cluster",
    )
    assert summary["valid_cluster_count"] == 3
    assert summary["nonzero_cluster_count"] == 2
    top_labels = [item["label"] for item in summary["top"]]
    assert top_labels[0] == "3_Pvalb"
    assert summary["top"][0]["ancestors"] == ["Gabaergic neurons"]
    ancestors = {a["ancestor"] for a in summary["top_taxonomy_ancestors"]}
    assert "Gabaergic neurons" in ancestors


def test_summary_artifact_carries_provenance_and_terminology():
    artifact = ac.build_celltype_summary_artifact(
        dataset=ac.DATASET_HUMAN_M1,
        summary=human_summary(),
        reconciliation=ac.reconcile_taxonomy(taxonomy_leaves=[], matrix_labels=[]),
        match_count=1,
        source_checksums={"trimmed_means": "abc", "taxonomy": "def"},
        source_urls={"trimmed_means": ac.human_m1_source_url("trimmed_means")},
    )
    assert artifact["assay_terminology"] == "single-nucleus RNA sequencing"
    assert artifact["sampling_scope"] == "sampled human primary motor cortex cell types"
    assert artifact["target_row_match_count"] == 1
    assert artifact["calculation_version"] == ac.CALCULATION_VERSION
    json.dumps(artifact)


def test_mouse_dataset_terminology_is_not_snrnaseq():
    # DropViz-style single-cell wording must not be relabelled, and the mouse
    # Allen dataset must not silently inherit the human snRNA-seq term.
    assert ac.DATASET_ASSAY_TERMS[ac.DATASET_MOUSE_CTX_HPF] == "single-cell RNA sequencing"
    assert ac.DATASET_ASSAY_TERMS[ac.DATASET_HUMAN_M1] != ac.DATASET_ASSAY_TERMS[
        ac.DATASET_MOUSE_CTX_HPF
    ]


# --------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------
def test_figure_url_is_gene_general():
    url = ac.figure_url(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
    )
    assert ac.EXPLORER_DATASET_SLUGS[ac.DATASET_HUMAN_M1] in url
    assert ac.DATASET_HUMAN_M1 not in url
    assert "colorByFeatureValue=GENEX" in url
    assert ac.is_allowlisted_figure_url(url)

    heatmap = ac.figure_url(
        dataset=ac.DATASET_MOUSE_CTX_HPF,
        gene_symbol="Genex",
        visualization=ac.VISUALIZATION_HEATMAP,
    )
    assert ac.EXPLORER_DATASET_SLUGS[ac.DATASET_MOUSE_CTX_HPF] in heatmap
    assert "Heatmap" in heatmap
    assert "colorByFeatureValue=Genex" in heatmap


def test_figure_url_rejects_unknown_inputs():
    with pytest.raises(KeyError):
        ac.figure_url(dataset="nope", gene_symbol="GENEX", visualization="scatter")
    with pytest.raises(KeyError):
        ac.figure_url(
            dataset=ac.DATASET_HUMAN_M1, gene_symbol="GENEX", visualization="violin"
        )
    with pytest.raises(ValueError):
        ac.figure_url(
            dataset=ac.DATASET_HUMAN_M1, gene_symbol="  ", visualization="scatter"
        )


def test_landing_page_detected():
    slug = ac.EXPLORER_DATASET_SLUGS[ac.DATASET_HUMAN_M1]
    assert ac.is_landing_page("https://celltypes.brain-map.org/rnaseq", dataset=ac.DATASET_HUMAN_M1)
    assert ac.is_landing_page(
        f"https://celltypes.brain-map.org/rnaseq/{slug}",
        dataset=ac.DATASET_HUMAN_M1,
    )
    # Internal dataset key in the path is not a valid explorer deep link.
    assert ac.is_landing_page(
        f"https://celltypes.brain-map.org/rnaseq/{ac.DATASET_HUMAN_M1}?x=1",
        dataset=ac.DATASET_HUMAN_M1,
    )
    assert not ac.is_landing_page(
        ac.figure_url(
            dataset=ac.DATASET_HUMAN_M1,
            gene_symbol="GENEX",
            visualization=ac.VISUALIZATION_SCATTER,
        ),
        dataset=ac.DATASET_HUMAN_M1,
    )


def test_dataset_source_pages_are_dataset_specific():
    human = ac.dataset_source_page(ac.DATASET_HUMAN_M1)
    mouse = ac.dataset_source_page(ac.DATASET_MOUSE_CTX_HPF)
    assert human.endswith("/human-m1-10x")
    assert "mouse-whole-cortex-and-hippocampus-10x" in mouse
    assert human != ac.database_hub_url()
    assert mouse != ac.database_hub_url()


def test_valid_figure_capture_accepted():
    url = ac.figure_url(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
    )
    result = ac.validate_figure_capture(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        final_url=url,
        gene_visible=True,
        plot_width=800,
        plot_height=600,
    )
    assert result["valid"] is True
    assert result["rejection_reasons"] == []


def test_landing_page_capture_rejected():
    result = ac.validate_figure_capture(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        final_url="https://celltypes.brain-map.org/rnaseq",
        gene_visible=True,
        plot_width=800,
        plot_height=600,
    )
    assert result["valid"] is False
    assert "landing_page" in result["rejection_reasons"]


def test_capture_requires_requested_gene_visible():
    url = ac.figure_url(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
    )
    result = ac.validate_figure_capture(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        final_url=url,
        gene_visible=False,
        plot_width=800,
        plot_height=600,
    )
    assert result["valid"] is False
    assert "requested_gene_not_visible" in result["rejection_reasons"]


def test_capture_requires_nonzero_plot_dimensions():
    url = ac.figure_url(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
    )
    zero = ac.validate_figure_capture(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        final_url=url,
        gene_visible=True,
        plot_width=0,
        plot_height=0,
    )
    assert "plot_container_missing" in zero["rejection_reasons"]
    tiny = ac.validate_figure_capture(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        final_url=url,
        gene_visible=True,
        plot_width=10,
        plot_height=10,
    )
    assert "plot_container_too_small" in tiny["rejection_reasons"]


def test_capture_rejected_when_source_error_visible():
    url = ac.figure_url(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
    )
    result = ac.validate_figure_capture(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        final_url=url,
        gene_visible=True,
        plot_width=800,
        plot_height=600,
        source_error_visible=True,
    )
    assert "source_error_visible" in result["rejection_reasons"]


def test_off_host_final_url_rejected():
    result = ac.validate_figure_capture(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        final_url="https://evil.example.com/rnaseq/human-m1-10x?x=1",
        gene_visible=True,
        plot_width=800,
        plot_height=600,
    )
    assert "final_url_not_allowlisted" in result["rejection_reasons"]


# --------------------------------------------------------------------------------------
# Playwright figure capture (browser fully mocked)
# --------------------------------------------------------------------------------------
class _FakeNode:
    """One DOM element with a bounding box and a screenshot."""

    def __init__(
        self,
        *,
        width: int,
        height: int,
        png: bytes = _MIN_PNG,
        visible: bool = True,
        screenshot_error: Exception | None = None,
        x: int = 0,
        y: int = 0,
        text: str = "",
    ):
        self._box = {"x": x, "y": y, "width": width, "height": height}
        self._png = png
        self._visible = visible
        self._screenshot_error = screenshot_error
        self._text = text
        self.screenshot_calls = 0

    def is_visible(self) -> bool:
        return self._visible

    def bounding_box(self):
        return dict(self._box)

    def inner_text(self, timeout: int = 0) -> str:
        return self._text

    def screenshot(self, type="png"):  # noqa: A002
        self.screenshot_calls += 1
        if self._screenshot_error is not None:
            raise self._screenshot_error
        return self._png


class _FakeCaptureLocator:
    def __init__(self, nodes: list[_FakeNode], text: str | None = None):
        self._nodes = nodes
        self._text = text

    def count(self) -> int:
        return len(self._nodes)

    def nth(self, index: int) -> _FakeNode:
        return self._nodes[index]

    def inner_text(self, timeout: int = 0) -> str:
        return self._text or ""


class _FakeFrame:
    """One child frame with its own body text and elements."""

    def __init__(
        self,
        *,
        body_text: str = "",
        nodes_by_selector: dict[str, list[_FakeNode]] | None = None,
    ):
        self._body_text = body_text
        self._nodes = nodes_by_selector or {}

    def locator(self, selector):
        if selector == "body":
            return _FakeCaptureLocator([], text=self._body_text)
        return _FakeCaptureLocator(list(self._nodes.get(selector) or []))


class _FakeCapturePage:
    """Minimal Playwright page for :func:`capture_explorer_figure`."""

    def __init__(
        self,
        *,
        final_url: str,
        body_text: str,
        nodes_by_selector: dict[str, list[_FakeNode]] | None = None,
        goto_error: Exception | None = None,
        child_frames: list[_FakeFrame] | None = None,
    ):
        self.url = final_url
        self._body_text = body_text
        self._nodes = nodes_by_selector or {}
        self._goto_error = goto_error
        self.goto_calls: list[str] = []
        self.waited_selectors: list[str] = []
        self._main_frame = _FakeFrame(
            body_text=body_text, nodes_by_selector=nodes_by_selector
        )
        self.frames = [self._main_frame, *(child_frames or [])]

    def goto(self, url, wait_until=None, timeout=0):
        self.goto_calls.append(url)
        if self._goto_error is not None:
            raise self._goto_error

    def wait_for_selector(self, selector, state=None, timeout=0):
        self.waited_selectors.append(selector)
        if selector not in self._nodes:
            raise TimeoutError(f"no {selector}")

    def wait_for_timeout(self, ms):
        return None

    def locator(self, selector):
        if selector == "body":
            return _FakeCaptureLocator([], text=self._body_text)
        return _FakeCaptureLocator(list(self._nodes.get(selector) or []))


def _dataset_url(dataset: str, gene: str, visualization: str) -> str:
    return ac.figure_url(
        dataset=dataset, gene_symbol=gene, visualization=visualization
    )


def _plot_selector() -> str:
    return ac.PLOT_CONTAINER_SELECTORS[0]


def _paired_scatter_nodes(
    *,
    png: bytes = _MIN_PNG,
    screenshot_error: Exception | None = None,
) -> list[_FakeNode]:
    return [
        _FakeNode(
            width=450,
            height=700,
            png=png,
            x=0,
            y=0,
            screenshot_error=screenshot_error,
        ),
        _FakeNode(width=450, height=700, png=png, x=500, y=0),
    ]


def _paired_scatter_body(gene: str) -> str:
    return f"Cell type Gene Expression {gene} expression scatter"


def test_capture_accepts_validated_plot_and_returns_png():
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_HUMAN_M1, gene, ac.VISUALIZATION_SCATTER)
    left = _FakeNode(
        width=450, height=700, png=b"\x89PNG\r\n\x1a\nscatter", x=0, y=0
    )
    right = _FakeNode(
        width=450, height=700, png=b"\x89PNG\r\n\x1a\nscatter", x=500, y=0
    )
    page = _FakeCapturePage(
        final_url=url,
        body_text=f"Cell type Gene Expression {gene} scatter expression",
        nodes_by_selector={_plot_selector(): [left, right]},
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is True
    assert result["status"] == ac.CAPTURE_STATUS_SUCCESS
    assert result["png"] == b"\x89PNG\r\n\x1a\nscatter"
    assert result["plot_width"] == 450
    assert result["plot_height"] == 700
    assert result["plot_panel_count"] == 2
    assert result["validation"]["valid"] is True
    assert result["retrieval_method"] == ac.CAPTURE_RETRIEVAL_METHOD
    assert page.goto_calls == [url]
    assert left.screenshot_calls + right.screenshot_calls == 1


def test_capture_rejects_landing_page_without_screenshot():
    gene = "GENEX"
    node = _FakeNode(width=900, height=700)
    page = _FakeCapturePage(
        final_url="https://knowledge.brain-map.org/abcatlas",
        body_text=f"{gene} explorer",
        nodes_by_selector={_plot_selector(): [node]},
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is False
    assert result["status"] == ac.CAPTURE_STATUS_UNAVAILABLE
    assert result["png"] is None
    assert result["error_type"] == "figure_validation_failed"
    assert "landing_page" in result["validation"]["rejection_reasons"]
    assert node.screenshot_calls == 0


def test_capture_rejects_page_that_never_shows_requested_gene():
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_MOUSE_CTX_HPF, gene, ac.VISUALIZATION_HEATMAP)
    page = _FakeCapturePage(
        final_url=url,
        body_text="Gene Expression SOMEOTHERGENE heatmap taxonomy cluster dendrogram",
        nodes_by_selector={_plot_selector(): [_FakeNode(width=900, height=700)]},
    )
    # Force a URL that does not carry the requested gene so canvas-only labels
    # cannot be confirmed via colorByFeatureValue.
    page.url = (
        "https://celltypes.brain-map.org/rnaseq/mouse_ctx-hpf_10x"
        "?selectedVisualization=Heatmap&colorByFeature=Gene+Expression"
        "&colorByFeatureValue=SOMEOTHERGENE"
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_MOUSE_CTX_HPF,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_HEATMAP,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is False
    assert "requested_gene_not_visible" in result["validation"]["rejection_reasons"]


def test_capture_rejects_undersized_plot_container():
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_HUMAN_M1, gene, ac.VISUALIZATION_SCATTER)
    small = ac.MIN_PLOT_DIMENSION - 1
    page = _FakeCapturePage(
        final_url=url,
        body_text=gene,
        nodes_by_selector={_plot_selector(): [_FakeNode(width=small, height=small)]},
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is False
    assert "plot_container_missing" in result["validation"]["rejection_reasons"]


def test_capture_reports_navigation_failure_without_raising():
    page = _FakeCapturePage(
        final_url="",
        body_text="",
        goto_error=TimeoutError("navigation timed out"),
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is False
    assert result["error_type"] == "navigation_failed"
    assert result["status"] == ac.CAPTURE_STATUS_UNAVAILABLE


def test_capture_reports_screenshot_failure_without_raising():
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_HUMAN_M1, gene, ac.VISUALIZATION_SCATTER)
    nodes = _paired_scatter_nodes(
        screenshot_error=RuntimeError("element detached")
    )
    page = _FakeCapturePage(
        final_url=url,
        body_text=_paired_scatter_body(gene),
        nodes_by_selector={_plot_selector(): nodes},
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is False
    assert result["error_type"] == "screenshot_failed"
    assert result["png"] is None


def test_capture_rejects_unknown_dataset_before_touching_a_browser():
    calls: list[int] = []

    def _factory():
        calls.append(1)
        raise AssertionError("browser must not open for an invalid request")

    result = ac.capture_explorer_figure(
        dataset="not-a-dataset",
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
        page_factory=_factory,
    )

    assert result["ok"] is False
    assert result["error_type"] == "invalid_request"
    assert calls == []


def test_capture_picks_largest_visible_plot_element():
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_HUMAN_M1, gene, ac.VISUALIZATION_SCATTER)
    # Two visible panels satisfy paired-scatter composition; the larger wins.
    small = _FakeNode(width=420, height=400, png=b"\x89PNG\r\n\x1a\nsmall", x=0, y=0)
    large = _FakeNode(width=1100, height=800, png=b"\x89PNG\r\n\x1a\nlarge", x=500, y=0)
    hidden = _FakeNode(width=1600, height=1200, visible=False, x=0, y=0)
    page = _FakeCapturePage(
        final_url=url,
        body_text=_paired_scatter_body(gene),
        nodes_by_selector={_plot_selector(): [small, hidden, large]},
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is True
    assert result["png"] == b"\x89PNG\r\n\x1a\nlarge"
    assert result["plot_width"] == 1100
    assert result["plot_panel_count"] == 2


def test_capture_finds_the_plot_inside_the_viewer_iframe():
    """The explorer embeds its viewer in an iframe, so frames must be searched."""
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_HUMAN_M1, gene, ac.VISUALIZATION_SCATTER)
    nodes = _paired_scatter_nodes(png=b"\x89PNG\r\n\x1a\nframe")
    page = _FakeCapturePage(
        final_url=url,
        body_text="",
        nodes_by_selector={},
        child_frames=[
            _FakeFrame(
                body_text=_paired_scatter_body(gene),
                nodes_by_selector={_plot_selector(): nodes},
            )
        ],
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is True
    assert result["png"] == b"\x89PNG\r\n\x1a\nframe"
    assert result["plot_panel_count"] == 2


def test_capture_reports_a_source_error_raised_inside_the_iframe():
    """A 404 in the embedded viewer is a source error, not a silent empty plot."""
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_HUMAN_M1, gene, ac.VISUALIZATION_SCATTER)
    page = _FakeCapturePage(
        final_url=url,
        body_text="",
        nodes_by_selector={},
        child_frames=[
            _FakeFrame(
                body_text=f"404\n/{ac.DATASET_HUMAN_M1} does not exist.\n{gene}",
                nodes_by_selector={_plot_selector(): [_FakeNode(width=900, height=700)]},
            )
        ],
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        settle_ms=0,
        plot_timeout_ms=0,
        page_factory=lambda: page,
    )

    assert result["ok"] is False
    assert "source_error_visible" in result["validation"]["rejection_reasons"]


def test_capture_reports_missing_playwright_without_raising(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def _fake_import(name, *args, **kwargs):
        if name == "playwright.sync_api":
            raise ModuleNotFoundError("No module named 'playwright'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _fake_import)

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol="GENEX",
        visualization=ac.VISUALIZATION_SCATTER,
    )

    assert result["ok"] is False
    assert result["error_type"] == "playwright_unavailable"


def test_capture_retries_then_reports_failure(monkeypatch):
    gene = "GENEX"
    url = _dataset_url(ac.DATASET_HUMAN_M1, gene, ac.VISUALIZATION_SCATTER)
    pages = [
        _FakeCapturePage(
            final_url=url,
            body_text=gene,
            nodes_by_selector={},
        ),
        _FakeCapturePage(
            final_url=url,
            body_text=_paired_scatter_body(gene),
            nodes_by_selector={_plot_selector(): _paired_scatter_nodes()},
        ),
    ]

    class _FakeBrowserContext:
        def __init__(self, page):
            self._page = page

        def new_page(self):
            return self._page

        def close(self):
            return None

    class _FakeBrowser:
        def __init__(self, page):
            self._page = page

        def new_context(self, viewport=None, user_agent=None):
            return _FakeBrowserContext(self._page)

        def close(self):
            return None

    class _FakePW:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    fake_sync = types.ModuleType("playwright.sync_api")
    fake_sync.sync_playwright = lambda: _FakePW()
    fake_playwright = types.ModuleType("playwright")
    fake_playwright.sync_api = fake_sync
    monkeypatch.setitem(sys.modules, "playwright", fake_playwright)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", fake_sync)

    order = iter(pages)
    monkeypatch.setattr(
        ac,
        "_launch_playwright_chromium",
        lambda pw: (
            _FakeBrowser(next(order)),
            [{"channel": "chromium", "success": True}],
        ),
    )

    result = ac.capture_explorer_figure(
        dataset=ac.DATASET_HUMAN_M1,
        gene_symbol=gene,
        visualization=ac.VISUALIZATION_SCATTER,
        max_attempts=2,
        settle_ms=0,
        plot_timeout_ms=0,
    )

    # First attempt has no plot container; the retry supplies a valid one.
    assert result["ok"] is True
    assert [a["attempt"] for a in result["attempts"]] == [1, 2]
    assert result["browser_channel"] == "chromium"


# --------------------------------------------------------------------------------------
# Production hygiene
# --------------------------------------------------------------------------------------
def test_no_validation_gene_hardcoded_in_module():
    text = SRC.read_text(encoding="utf-8")
    for token in ("SREBF2", "Srebf2", "CDH10", "Cdh10"):
        assert token not in text, f"{token} must not appear in production source"


def test_module_does_not_fabricate_api_runs():
    assert "ApiRun(" not in SRC.read_text(encoding="utf-8")
