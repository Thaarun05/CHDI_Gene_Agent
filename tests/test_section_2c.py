"""Offline tests for Section 2c snRNA-Seq gene expression in cell type database.

Every fixture is synthetic and gene-general: the target symbols are ``GENEX`` /
``Genex`` so nothing here depends on a validation gene. No network, no browser.
"""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path
from typing import Any

import pytest

from gene_dossier.config import Settings, get_settings
from gene_dossier.report_presentation import build_section_2c_blocks
from gene_dossier.report_schema import PresentationRole, ReportSubsection
from gene_dossier.rancho_report import (
    render_section_2c_subsection_segments,
    split_section_2c_page_segments,
)
from gene_dossier.section_2c import (
    CALCULATION_VERSION,
    CHART_VERSION,
    GEO_SOURCE_KEY,
    SECTION_2C_INTRO_TEXT,
    STATUS_FIGURE_SUPPRESSED,
    STATUS_NOT_ATTEMPTED,
    STATUS_SOURCE_UNAVAILABLE,
    STATUS_SUCCESS,
    Section2cConfig,
    node_generate_section_2c_derived_artifacts,
)
from gene_dossier.section_2c_sources import accept_source, paths_for, sha256_bytes
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SUPPORTED_SECTION_BUNDLE_KEYS,
    build_section_bundle_document,
    render_section_bundle_html,
    sources_for_sections,
    validate_section_keys,
)
from gene_dossier.tools import allen_celltypes as ac
from gene_dossier.tools import dropviz_geo as dg

SRC = Path(__file__).resolve().parents[1] / "src" / "gene_dossier" / "section_2c.py"

HUMAN_SYMBOL = "GENEX"
MOUSE_SYMBOL = "Genex"

# ---------------------------------------------------------------------------
# Synthetic dataset-level sources
# ---------------------------------------------------------------------------
HUMAN_CSV = """feature,Astro_1,Inh Lamp5_1,Exc L2/3 IT_1,Vip_1,Oligo_1
GENEX,0.0,4.5,6.25,1.5,0.5
GENEX-AS1,9.9,9.9,9.9,9.9,9.9
OTHER,1.0,1.0,1.0,1.0,1.0
"""

MOUSE_CSV = """feature,1_CR,2_Lamp5,3_Pvalb,4_Oligo
Genex,0.0,2.0,6.0,3.0
Other,1.0,1.0,1.0,1.0
"""

# Column totals are all 100, so Genex is 5,000 / 20,000 / 10,000 / 0 per
# 100,000 and the descending order is pop_B > pop_C > pop_A > pop_D.
GEO_CSV = """gene,pop_A,pop_B,pop_C,pop_D
Gene1,10,0,5,20
Genex,5,20,10,0
Gene3,85,80,85,80
"""

# Fractional, non-integral values: value semantics cannot be resolved as counts.
GEO_CSV_FRACTIONAL = """gene,pop_A,pop_B,pop_C,pop_D
Gene1,0.15,0.02,0.51,0.29
Genex,0.35,0.78,0.44,0.02
Gene3,0.55,0.20,0.05,0.69
"""


def _leaf(alias: str, *, structure: str = "Cortical plate") -> dict[str, Any]:
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


def _internal(alias: str, children: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "node_attributes": [{"cell_set_alias": alias, "cell_set_designation": alias}],
        "children": children,
    }


HUMAN_DEND = _internal(
    "",
    [
        _internal("GABAergic", [_leaf("Inh Lamp5_1"), _leaf("Vip_1")]),
        _internal("Glutamatergic", [_leaf("Exc L2/3 IT_1")]),
        _internal("Non-neuronal", [_leaf("Astro_1"), _leaf("Oligo_1")]),
    ],
)

MOUSE_DEND = _internal(
    "",
    [
        _internal("GABAergic neurons", [_leaf("2_Lamp5"), _leaf("3_Pvalb")]),
        _internal("Glutamatergic neurons", [_leaf("1_CR")]),
        _internal("Non-neuronal", [_leaf("4_Oligo")]),
    ],
)


def _gzip(text: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as fh:
        fh.write(text.encode("utf-8"))
    return buf.getvalue()


def _accept(paths, *, source_key: str, content: bytes, official_url: str | None) -> None:
    attempt = paths.new_source_attempt(source_key)
    artifact = attempt / "payload.bin"
    artifact.write_bytes(content)
    accept_source(
        paths,
        source_key=source_key,
        attempt_dir=attempt,
        artifact_path=artifact,
        official_url=official_url or "",
        sha256=sha256_bytes(content),
        byte_size=len(content),
        validation={"ok": True, "fixture": True},
    )


def _seed_sources(root: Path, *, geo_csv: str = GEO_CSV, skip: set[str] | None = None):
    """Accept synthetic dataset-level sources under ``root``."""
    skipped = skip or set()
    paths = paths_for(root)
    paths.ensure()
    payloads = {
        ac.CACHE_KEY_HUMAN_TRIMMED_MEANS: (
            HUMAN_CSV.encode("utf-8"),
            ac.human_m1_source_url("trimmed_means"),
        ),
        ac.CACHE_KEY_HUMAN_TAXONOMY: (
            json.dumps(HUMAN_DEND).encode("utf-8"),
            ac.human_m1_source_url("taxonomy"),
        ),
        ac.CACHE_KEY_MOUSE_TRIMMED_MEANS: (MOUSE_CSV.encode("utf-8"), None),
        ac.CACHE_KEY_MOUSE_TAXONOMY: (json.dumps(MOUSE_DEND).encode("utf-8"), None),
        GEO_SOURCE_KEY: (_gzip(geo_csv), dg.supplementary_url()),
    }
    for source_key, (content, url) in payloads.items():
        if source_key in skipped:
            continue
        _accept(paths, source_key=source_key, content=content, official_url=url)
    return paths


def _run_section(
    tmp_path: Path,
    *,
    gene: str = HUMAN_SYMBOL,
    mouse_symbol: str | None = MOUSE_SYMBOL,
    geo_csv: str = GEO_CSV,
    skip: set[str] | None = None,
    config: Section2cConfig | None = None,
    run_id: str = "test-2c-run",
) -> dict[str, Any]:
    """Run the Section 2c node fully offline against synthetic accepted sources."""
    _seed_sources(tmp_path / "outputs", geo_csv=geo_csv, skip=skip)
    settings = Settings(
        raw_data_dir=tmp_path / "raw",
        output_dir=tmp_path / "outputs",
    )
    state: dict[str, Any] = {
        "run_type": "section_bundle",
        "selected_section_keys": ["2c"],
        "dossier_run_id": run_id,
        "gene_symbol": gene,
        "gene_ids": {"mouse_symbol": mouse_symbol} if mouse_symbol else {},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }
    return node_generate_section_2c_derived_artifacts(
        state,
        settings=settings,
        persist_db=False,
        config=config
        or Section2cConfig(
            output_root=tmp_path / "outputs",
            attempt_allen_figures=False,
        ),
    )


@pytest.fixture(scope="module")
def section_run(tmp_path_factory) -> dict[str, Any]:
    return _run_section(tmp_path_factory.mktemp("section2c"))


@pytest.fixture(scope="module")
def status(section_run) -> dict[str, Any]:
    return section_run["section_2c_status"]


def _evidence_payload(status: dict[str, Any]) -> dict[str, Any]:
    audit = status["audit"]
    path = Path(audit["artifact_root"]) / audit["artifacts"]["section_2c_evidence.json"]
    return json.loads(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Bundle wiring
# ---------------------------------------------------------------------------
def test_section_key_is_2c_and_opt_in():
    assert "2c" in SUPPORTED_SECTION_BUNDLE_KEYS
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
    assert "2c" in DEFAULT_SECTION_BUNDLE_KEYS
    assert validate_section_keys(["2.c"]) == ["2c"]
    assert validate_section_keys(["2c"]) == ["2c"]


def test_2c_is_ordered_after_2b():
    assert validate_section_keys(["2c", "2a", "2b"]) == ["2a", "2b", "2c"]
    keys = list(SUPPORTED_SECTION_BUNDLE_KEYS)
    assert keys.index("2c") == keys.index("2b") + 1


def test_section_owns_its_sources_so_generic_clients_stay_quiet():
    assert sources_for_sections(["2c"]) == []
    assert "Allen Brain Atlas" not in sources_for_sections(["1a", "2c"])
    assert "GEO" not in sources_for_sections(["1a", "2c"])


def test_subsection_title_matches_the_canonical_toc():
    from gene_dossier.report_schema import REPORT_SECTIONS

    major2 = next(sec for sec in REPORT_SECTIONS if sec.number == 2)
    sub = next(s for s in major2.subsections if s.key == "c")
    assert sub.title == "snRNA-Seq gene expression in cell type database"


def test_presentation_roles_registered():
    registered = set(PresentationRole.__args__)
    for role in (
        "section_2c_intro",
        "section_2c_human_narrative",
        "section_2c_human_table",
        "section_2c_human_scatter_figure",
        "section_2c_human_heatmap_figure",
        "section_2c_mouse_narrative",
        "section_2c_mouse_table",
        "section_2c_dropviz_narrative",
        "section_2c_dropviz_table",
        "section_2c_dropviz_rank_figure",
        "section_2c_therapeutic_narrative",
        "section_2c_source_link",
        "section_2c_source_status",
    ):
        assert role in registered


# ---------------------------------------------------------------------------
# Independent branch statuses
# ---------------------------------------------------------------------------
def test_all_three_branches_succeed_on_synthetic_sources(status):
    rendering = status["rendering_status"]
    assert rendering["allen_human_analysis_status"] == STATUS_SUCCESS
    assert rendering["allen_mouse_analysis_status"] == STATUS_SUCCESS
    assert rendering["dropviz_geo_matrix_status"] == STATUS_SUCCESS
    assert rendering["dropviz_rank_status"] == dg.RANK_STATUS_SUCCESS
    assert rendering["overall"] == STATUS_SUCCESS


def test_live_dropviz_client_is_never_invoked_in_production(status):
    assert status["rendering_status"]["dropviz_live_status"] == STATUS_NOT_ATTEMPTED
    assert "tools.dropviz" not in SRC.read_text(encoding="utf-8")
    assert "DropVizClient" not in SRC.read_text(encoding="utf-8")


def test_confidence_interval_stays_method_unresolved(status):
    assert (
        status["rendering_status"]["confidence_interval_status"]
        == dg.CONFIDENCE_INTERVAL_METHOD_UNRESOLVED
    )


def test_unavailable_allen_human_source_leaves_other_branches_intact(tmp_path):
    out = _run_section(
        tmp_path, skip={ac.CACHE_KEY_HUMAN_TRIMMED_MEANS}, run_id="human-missing"
    )
    rendering = out["section_2c_status"]["rendering_status"]
    assert rendering["allen_human_analysis_status"] == STATUS_SOURCE_UNAVAILABLE
    assert rendering["allen_mouse_analysis_status"] == STATUS_SUCCESS
    assert rendering["dropviz_rank_status"] == dg.RANK_STATUS_SUCCESS
    assert rendering["overall"] == "partial"

    summary = out["section_2c_status"]["summary"]
    # A failed branch never discards the successful ones.
    assert summary["mouse"]["top"]
    assert summary["dropviz"]["top_populations"]
    assert summary["human"]["narrative"]
    assert not summary["human"]["top"]


def test_unavailable_geo_source_keeps_allen_branches(tmp_path):
    out = _run_section(tmp_path, skip={GEO_SOURCE_KEY}, run_id="geo-missing")
    rendering = out["section_2c_status"]["rendering_status"]
    assert rendering["dropviz_geo_matrix_status"] == STATUS_SOURCE_UNAVAILABLE
    assert rendering["allen_human_analysis_status"] == STATUS_SUCCESS
    assert rendering["allen_mouse_analysis_status"] == STATUS_SUCCESS
    assert rendering["overall"] == "partial"
    assert out["section_2c_status"]["summary"]["human"]["top"]


# ---------------------------------------------------------------------------
# Quantitative content
# ---------------------------------------------------------------------------
def test_human_counts_use_canonical_metric_names(status):
    human = status["summary"]["human"]
    assert human["valid_celltype_count"] == 5
    assert human["nonzero_celltype_count"] == 4
    assert human["nonzero_percentage"] == pytest.approx(80.0)
    assert human["maximum"] == pytest.approx(6.25)


def test_mouse_counts_use_cluster_nouns(status):
    mouse = status["summary"]["mouse"]
    assert mouse["valid_cluster_count"] == 4
    assert mouse["nonzero_cluster_count"] == 3
    assert mouse["maximum"] == pytest.approx(6.0)
    assert mouse["top"][0]["label"] == "3_Pvalb"


def test_dropviz_ranking_is_descending_transcripts_per_100k(status):
    dropviz = status["summary"]["dropviz"]
    assert dropviz["value_semantics_status"] == dg.VALUE_SEMANTICS_RAW_COUNTS
    assert dropviz["ranking_metric"] == "transcripts_per_100k"
    assert dropviz["axis_label"] == "Transcripts per 100,000"
    labels = [row["population_label"] for row in dropviz["top_populations"]]
    assert labels == ["pop_B", "pop_C", "pop_A", "pop_D"]
    values = [row["ranking_value"] for row in dropviz["top_populations"]]
    assert values == sorted(values, reverse=True)
    assert values[0] == pytest.approx(20_000.0)
    assert dropviz["valid_population_count"] == 4
    assert dropviz["nonzero_population_count"] == 3
    assert dropviz["excluded_population_count"] == 0
    assert dropviz["target_row_match_count"] == 1


def test_dropviz_reports_both_share_metrics(status):
    dropviz = status["summary"]["dropviz"]
    assert dropviz["top_10_normalized_expression_share"] is not None
    assert dropviz["top_10_raw_target_count_share"] is not None
    assert 0.0 <= float(dropviz["top_10_normalized_expression_share"]) <= 100.0
    assert 0.0 <= float(dropviz["top_10_raw_target_count_share"]) <= 100.0


# ---------------------------------------------------------------------------
# Derived artifacts
# ---------------------------------------------------------------------------
def test_expected_derived_artifacts_are_persisted(status):
    artifacts = status["audit"]["artifacts"]
    for name in (
        "allen_human_celltype_summary.json",
        "allen_mouse_celltype_summary.json",
        "dropviz_geo_matrix_profile.json",
        "dropviz_population_expression_raw.csv",
        "dropviz_population_expression_ranked.csv",
        "dropviz_top_populations.json",
        "dropviz_population_audit.json",
        f"{CHART_VERSION}.png",
        "section_2c_evidence.json",
    ):
        assert artifacts.get(name), f"missing derived artifact {name}"


def test_derived_artifacts_declare_parents_and_fabricate_no_api_run(section_run):
    derived = [
        meta
        for meta in section_run["raw_artifacts"]
        if meta.get("artifact_class") == "derived"
    ]
    assert derived
    for meta in derived:
        # A local computation never fabricates an ApiRun.
        assert meta.get("api_run_id") in (None, "")
        assert "parent_raw_artifact_ids" in meta


def test_raw_csv_keeps_source_order_and_ranked_csv_is_descending(section_run):
    artifacts = section_run["section_2c_status"]["audit"]["artifacts"]
    root = Path(section_run["section_2c_status"]["audit"]["artifact_root"])
    raw = (root / artifacts["dropviz_population_expression_raw.csv"]).read_text()
    ranked = (root / artifacts["dropviz_population_expression_ranked.csv"]).read_text()

    raw_labels = [line.split(",")[0] for line in raw.strip().splitlines()[1:]]
    # The ranked export leads with an explicit rank column.
    ranked_lines = ranked.strip().splitlines()
    assert ranked_lines[0].split(",")[0] == "rank"
    ranked_labels = [line.split(",")[1] for line in ranked_lines[1:]]
    assert raw_labels == ["pop_A", "pop_B", "pop_C", "pop_D"]
    assert ranked_labels == ["pop_B", "pop_C", "pop_A", "pop_D"]
    assert [line.split(",")[0] for line in ranked_lines[1:]] == ["1", "2", "3", "4"]


def test_evidence_json_is_compact_and_omits_raw_matrices(status):
    payload = _evidence_payload(status)
    assert payload["section_key"] == "2c"
    assert payload["calculation_version"] == CALCULATION_VERSION
    assert payload["requested_gene_symbol"] == HUMAN_SYMBOL
    assert payload["resolved_mouse_symbol"] == MOUSE_SYMBOL
    assert set(payload["statuses"]) >= {
        "allen_human_analysis_status",
        "allen_human_scatter_status",
        "allen_human_heatmap_status",
        "allen_mouse_analysis_status",
        "allen_mouse_scatter_status",
        "allen_mouse_heatmap_status",
        "dropviz_background_status",
        "dropviz_geo_matrix_status",
        "dropviz_rank_status",
        "dropviz_live_status",
        "confidence_interval_status",
        "label_mapping_status",
        "value_semantics_status",
    }
    text = json.dumps(payload)
    for gene_row in ("Gene1", "Gene3", "OTHER"):
        assert gene_row not in text, "evidence JSON must not embed the raw matrix"


# ---------------------------------------------------------------------------
# Chart rules
# ---------------------------------------------------------------------------
def test_point_only_chart_is_emitted_when_ci_is_unresolved(section_run):
    status = section_run["section_2c_status"]
    assert status["rendering_status"]["dropviz_rank_figure_status"] == STATUS_SUCCESS
    assert status["summary"]["dropviz"]["figure_status"] == STATUS_SUCCESS
    evidence = _evidence_payload(status)
    figure = evidence["dropviz"]["rank_figure"]
    assert (
        figure["confidence_interval_status"]
        == dg.CONFIDENCE_INTERVAL_METHOD_UNRESOLVED
    )
    assert figure["confidence_intervals_drawn"] is False
    assert figure["sha256"]
    root = Path(status["audit"]["artifact_root"])
    png = root / status["audit"]["artifacts"][f"{CHART_VERSION}.png"]
    assert png.read_bytes().startswith(b"\x89PNG")


def test_no_chart_when_value_semantics_leave_ranking_unresolved(tmp_path):
    out = _run_section(tmp_path, geo_csv=GEO_CSV_FRACTIONAL, run_id="unresolved-semantics")
    status = out["section_2c_status"]
    rendering = status["rendering_status"]
    assert rendering["dropviz_rank_status"] == dg.RANK_STATUS_NORMALIZATION_UNRESOLVED
    assert rendering["dropviz_rank_figure_status"] == STATUS_FIGURE_SUPPRESSED
    assert f"{CHART_VERSION}.png" not in status["audit"]["artifacts"]
    assert rendering["overall"] == "partial"
    # The unresolved metric is reported rather than relabelled.
    assert status["summary"]["dropviz"]["ranking_metric"] is None
    fields = {issue["field"] for issue in status["summary"]["unresolved_issues"]}
    assert "value_semantics_status" in fields
    assert "dropviz_rank_status" in fields


def test_chart_axis_never_relabels_a_non_count_metric(tmp_path):
    out = _run_section(tmp_path, geo_csv=GEO_CSV_FRACTIONAL, run_id="axis-label")
    dropviz = out["section_2c_status"]["summary"]["dropviz"]
    assert dropviz["axis_label"] in (None, "")
    assert "Transcripts per 100,000" not in json.dumps(dropviz)


# ---------------------------------------------------------------------------
# Deterministic narrative
# ---------------------------------------------------------------------------
def test_intro_names_each_dataset_by_its_own_assay(status):
    intro = status["summary"]["intro_text"]
    assert intro == SECTION_2C_INTRO_TEXT
    assert "Drop-seq" in intro
    assert "single-nucleus" in intro


def test_human_branch_uses_single_nucleus_terminology(status):
    narrative = status["summary"]["human"]["narrative"]
    assert "single-nucleus RNA sequencing" in narrative
    assert "Drop-seq" not in narrative


def test_mouse_branch_uses_allen_dataset_terminology_not_snrnaseq(status):
    mouse = status["summary"]["mouse"]
    assert mouse["assay_terminology"] == ac.DATASET_ASSAY_TERMS[ac.DATASET_MOUSE_CTX_HPF]
    assert "single-nucleus" not in mouse["narrative"]


def test_dropviz_branch_is_described_as_dropseq_single_cell(status):
    narrative = status["summary"]["dropviz"]["narrative"]
    assert "single-cell RNA sequencing using Drop-seq" in narrative
    assert dg.GEO_ACCESSION in narrative


def test_nonzero_aggregate_expression_phrasing_replaces_detection(status):
    human = status["summary"]["human"]["narrative"]
    assert "nonzero aggregate expression across" in human
    assert "Within the Human M1 dataset, the gene showed" in human
    for narrative in (
        human,
        status["summary"]["mouse"]["narrative"],
        status["summary"]["dropviz"]["narrative"],
    ):
        assert "was detected in" not in narrative


def test_claims_never_generalize_to_the_whole_human_brain(status):
    blob = json.dumps(status["summary"])
    assert "throughout the human brain" not in blob
    assert "broadly expressed throughout" not in blob
    assert (
        "localization within sampled human primary motor cortex cell types"
        in status["summary"]["human"]["narrative"]
    )
    assert (
        "localization within sampled mouse cortex and hippocampus populations"
        in status["summary"]["mouse"]["narrative"]
    )
    assert (
        "broader adult mouse-brain population context"
        in status["summary"]["dropviz"]["narrative"]
    )


def test_therapeutic_narrative_is_cautious_and_reaches_no_verdict(status):
    therapeutic = status["summary"]["therapeutic_narrative"]
    assert "assuming the therapeutic reaches the relevant tissue" in therapeutic
    assert "may have greater potential for on-target pharmacology" in therapeutic
    assert "RNA abundance alone does not establish delivery" in therapeutic
    assert (
        "relative RNA-expression rankings and maximum-to-median ratios" in therapeutic
    )
    assert "Expression localization:" in therapeutic or "Expression localization for" in therapeutic
    assert "Expression breadth:" in therapeutic
    assert "Potential safety context:" in therapeutic
    # Dynamic localization names and nonzero counts remain.
    assert "Exc L2/3 IT_1" in therapeutic or "Inh Lamp5_1" in therapeutic
    assert "of" in therapeutic and "sampled" in therapeutic
    # Required caveat: RNA abundance is not protein, dependence, or tractability.
    assert "do not establish protein abundance" in therapeutic
    assert "not a target assessment" in therapeutic
    for banned in (
        "would reach",
        "would receive the greatest on-target exposure",
        "would be exposed to the same on-target activity",
        "concentration ratios",
    ):
        assert banned not in therapeutic
    for verdict in (
        "is a promising target",
        "is a viable target",
        "should be prioritized",
        "we recommend",
        "validated target",
    ):
        assert verdict not in therapeutic


def test_overlapping_taxonomy_ancestors_prose_is_dynamic_and_once(status, section_run):
    human_narrative = status["summary"]["human"]["narrative"]
    mouse_narrative = status["summary"]["mouse"]["narrative"]
    assert "overlapping taxonomy ancestors" in human_narrative
    assert "overlapping taxonomy ancestors" in mouse_narrative
    assert "highest-ranked cell types" in human_narrative
    assert "highest-ranked clusters" in mouse_narrative
    # Fixture dendrogram ancestors appear with dynamic counts.
    assert "GABAergic" in human_narrative
    assert "GABAergic neurons" in mouse_narrative or "Non-neuronal" in mouse_narrative
    result = _blocks(status, section_run["evidence_records"])
    blob = " ".join(b.text for b in result.blocks if b.text)
    human_count = sum(
        1
        for b in result.blocks
        if b.text and "highest-ranked cell types" in b.text
        and "overlapping taxonomy ancestors" in b.text
    )
    mouse_count = sum(
        1
        for b in result.blocks
        if b.text and "highest-ranked clusters" in b.text
        and "overlapping taxonomy ancestors" in b.text
    )
    assert human_count <= 1
    assert mouse_count <= 1
    # Ancestor counts are not summed into a single exclusive total.
    assert "mutually exclusive" not in blob.lower()
    assert "together account for" not in blob.lower()


def test_cge_most_is_cge_derived_in_prose_only():
    from gene_dossier.section_2c import (
        _overlapping_taxonomy_ancestors_sentence,
        allen_branch_narrative,
    )

    ancestors = [
        {"ancestor": "GABAergic", "top_member_count": 2},
        {"ancestor": "CGE (most)", "top_member_count": 2},
        {"ancestor": "Vip", "top_member_count": 1},
    ]
    analysis = {
        "dataset": ac.DATASET_HUMAN_M1,
        "dataset_label": "Human M1 10x",
        "assay_terminology": "single-nucleus RNA sequencing",
        "sampling_scope": "sampled human primary motor cortex cell types",
        "count_noun": "celltype",
        "summary": {
            "valid_celltype_count": 5,
            "nonzero_celltype_count": 4,
            "nonzero_percentage": 80.0,
            "median": 1.0,
            "max_to_median_ratio": 2.0,
            "top": [
                {"label": "Vip_1", "value": 2.0},
                {"label": "Lamp5_1", "value": 1.5},
            ],
            "top_taxonomy_ancestors": ancestors,
        },
        "taxonomy_reconciliation": {},
    }
    prose = allen_branch_narrative(analysis, gene_symbol=HUMAN_SYMBOL)
    assert "CGE-derived" in prose
    assert "CGE (most)" not in prose
    assert "overlapping taxonomy ancestors" in prose
    sentence = _overlapping_taxonomy_ancestors_sentence(
        ancestors,
        top_count=2,
        entity_noun="cell types",
    )
    assert sentence is not None
    assert "Among the 2 highest-ranked cell types" in sentence
    # Evidence values remain the source label.
    assert ancestors[1]["ancestor"] == "CGE (most)"


def test_figure_success_puts_taxonomy_only_in_scatter_interpretation():
    """Figure-led path uses short intro + scatter note; fallback uses full narrative."""
    from gene_dossier.section_2c import _branch_presentation

    analysis = {
        "ok": True,
        "dataset": ac.DATASET_HUMAN_M1,
        "dataset_label": "Human M1 10x",
        "assay_terminology": "single-nucleus RNA sequencing",
        "sampling_scope": "sampled human primary motor cortex cell types",
        "count_noun": "celltype",
        "summary": {
            "valid_celltype_count": 5,
            "nonzero_celltype_count": 4,
            "nonzero_percentage": 80.0,
            "median": 1.5,
            "max_to_median_ratio": 2.0,
            "top": [
                {"label": "Exc L2/3 IT_1", "value": 6.25},
                {"label": "Inh Lamp5_1", "value": 4.5},
            ],
            "top_taxonomy_ancestors": [
                {"ancestor": "GABAergic", "top_member_count": 1},
                {"ancestor": "Glutamatergic", "top_member_count": 1},
            ],
            "source_symbol": HUMAN_SYMBOL,
        },
        "taxonomy_reconciliation": {},
    }
    success = _branch_presentation(
        analysis=analysis,
        gene_symbol=HUMAN_SYMBOL,
        label_column="Human M1 transcriptomic cell type",
        count_noun="celltype",
        explorer_symbol=HUMAN_SYMBOL,
        figure_notes=[],
        figure_statuses={
            ac.VISUALIZATION_SCATTER: STATUS_SUCCESS,
            ac.VISUALIZATION_HEATMAP: STATUS_SUCCESS,
        },
        top_n=10,
    )
    assert "overlapping taxonomy ancestors" not in (success["narrative"] or "")
    assert "overlapping taxonomy ancestors" in (success["scatter_interpretation"] or "")
    fallback = _branch_presentation(
        analysis=analysis,
        gene_symbol=HUMAN_SYMBOL,
        label_column="Human M1 transcriptomic cell type",
        count_noun="celltype",
        explorer_symbol=HUMAN_SYMBOL,
        figure_notes=[],
        figure_statuses={
            ac.VISUALIZATION_SCATTER: STATUS_SOURCE_UNAVAILABLE,
            ac.VISUALIZATION_HEATMAP: STATUS_SOURCE_UNAVAILABLE,
        },
        top_n=10,
    )
    assert "overlapping taxonomy ancestors" in (fallback["narrative"] or "")
    assert fallback["scatter_interpretation"] is None


def test_comparability_and_dropviz_notes_appear_once(section_run, status):
    comparability = (
        "Absolute expression values should be interpreted within each dataset and "
        "should not be compared directly across Human M1, Mouse CTX-HPF, and DropViz "
        "because the assays, aggregation procedures, and expression scales differ."
    )
    population_note = (
        "Population identifiers encode brain-region and cluster numbers; descriptive "
        "cell-class mapping was unavailable."
    )
    tsne_note = (
        "Historical DropViz regional t-SNE views were not included because the saved "
        "application states are no longer reproducible; production results use the "
        "archived GSE116470 matrix."
    )
    assert status["summary"]["intro_text"].count(comparability) == 1
    assert status["summary"]["dropviz"]["population_identifier_note"] == population_note
    assert status["summary"]["dropviz"]["regional_tsne_limitation_note"] == tsne_note
    assert status["rendering_status"]["dropviz_regional_tsne_status"] == (
        "not_in_production_scope"
    )
    # In-narrative label-mapping sentence was replaced (limitation once only).
    assert "Population labels resolve to region and cluster identifiers only" not in (
        status["summary"]["dropviz"]["narrative"] or ""
    )
    result = _blocks(status, section_run["evidence_records"])
    texts = [b.text for b in result.blocks if b.text]
    blob = "\n".join(texts)
    assert blob.count(comparability) == 1
    assert blob.count(population_note) == 1
    assert blob.count(tsne_note) == 1
    assert "_state_id_" not in blob
    assert "dropviz.org" not in blob.lower()
    # Population note sits before the DropViz table.
    roles = [b.presentation_role for b in result.blocks]
    pop_idx = next(
        i for i, b in enumerate(result.blocks) if b.text == population_note
    )
    table_idx = roles.index("section_2c_dropviz_table")
    assert pop_idx < table_idx
    # t-SNE limitation sits after DropViz content and before therapeutic.
    tsne_idx = next(i for i, b in enumerate(result.blocks) if b.text == tsne_note)
    therapeutic_idx = roles.index("section_2c_therapeutic_narrative")
    assert table_idx < tsne_idx < therapeutic_idx


def test_narrative_reports_taxonomy_reconciliation_without_zero_filling(status):
    recon = status["summary"]["mouse"]["taxonomy_reconciliation"]
    assert recon["taxonomy_leaf_count"] == 4
    assert recon["expression_cluster_count"] == 4
    assert recon["missing_expression_clusters"] == []


def test_missing_taxonomy_cluster_is_reported_not_zero_filled(tmp_path):
    """A taxonomy leaf with no expression column must never become a zero."""
    dend = _internal(
        "",
        [
            _internal("GABAergic neurons", [_leaf("2_Lamp5"), _leaf("3_Pvalb")]),
            _internal("Glutamatergic neurons", [_leaf("1_CR")]),
            _internal("Non-neuronal", [_leaf("4_Oligo"), _leaf("5_SMC-Peri")]),
        ],
    )
    paths = _seed_sources(tmp_path / "outputs", skip={ac.CACHE_KEY_MOUSE_TAXONOMY})
    _accept(
        paths,
        source_key=ac.CACHE_KEY_MOUSE_TAXONOMY,
        content=json.dumps(dend).encode("utf-8"),
        official_url=None,
    )
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "outputs")
    out = node_generate_section_2c_derived_artifacts(
        {
            "run_type": "section_bundle",
            "selected_section_keys": ["2c"],
            "dossier_run_id": "missing-cluster",
            "gene_symbol": HUMAN_SYMBOL,
            "gene_ids": {"mouse_symbol": MOUSE_SYMBOL},
            "evidence_records": [],
            "api_runs": [],
            "raw_artifacts": [],
            "errors": [],
            "coverage": [],
        },
        settings=settings,
        persist_db=False,
        config=Section2cConfig(
            output_root=tmp_path / "outputs", attempt_allen_figures=False
        ),
    )
    mouse = out["section_2c_status"]["summary"]["mouse"]
    recon = mouse["taxonomy_reconciliation"]
    assert recon["missing_expression_clusters"] == ["5_SMC-Peri"]
    assert mouse["valid_cluster_count"] == 4
    assert "5_SMC-Peri" in mouse["narrative"]
    assert "recorded as missing rather than filled with a zero" in mouse["narrative"]


# ---------------------------------------------------------------------------
# Presentation + rendering
# ---------------------------------------------------------------------------
def _blocks(status: dict[str, Any], records=None):
    result = build_section_2c_blocks(
        gene_symbol=status["summary"]["gene_symbol"],
        evidence_records=list(records if records is not None else []),
        section_status=status,
    )
    return result


def test_presentation_renders_every_successful_branch(section_run, status):
    result = _blocks(status, section_run["evidence_records"])
    roles = [b.presentation_role for b in result.blocks]
    assert roles[0] == "section_2c_intro"
    assert "section_2c_human_narrative" in roles
    assert "section_2c_human_table" in roles
    assert "section_2c_mouse_narrative" in roles
    assert "section_2c_mouse_table" in roles
    assert "section_2c_dropviz_narrative" in roles
    assert "section_2c_dropviz_table" in roles
    assert roles[-1] == "section_2c_therapeutic_narrative"


def test_no_tsne_block_and_no_homepage_or_sample_image(section_run, status):
    result = _blocks(status, section_run["evidence_records"])
    blob = json.dumps(
        [
            {
                "role": b.presentation_role,
                "text": b.text,
                "figure_path": b.figure_path,
                "links": b.links,
            }
            for b in result.blocks
        ]
    )
    for banned in (
        "dropviz.org",
        "sample image",
        "homepage",
        "screenshot",
        "_state_id_",
    ):
        assert banned not in blob.lower()
    # Limitation prose may mention t-SNE, but no t-SNE figure/link is rendered.
    assert "section_2c_dropviz_tsne" not in blob
    assert not any(
        (b.presentation_role or "").endswith("tsne_figure") for b in result.blocks
    )

def test_quantitative_summary_renders_when_both_figures_unavailable(status):
    """A browser failure must never invalidate a valid structured analysis."""
    degraded = json.loads(json.dumps(status, default=str))
    degraded["rendering_status"]["allen_human_scatter_status"] = STATUS_SOURCE_UNAVAILABLE
    degraded["rendering_status"]["allen_human_heatmap_status"] = STATUS_SOURCE_UNAVAILABLE
    degraded["summary"]["human"]["figure_status_notes"] = [
        "The Allen Human M1 10x scatter figure is unavailable for this run.",
        "The Allen Human M1 10x heatmap figure is unavailable for this run.",
    ]

    result = _blocks(degraded)
    roles = [b.presentation_role for b in result.blocks]
    assert "section_2c_human_narrative" in roles
    assert "section_2c_human_table" in roles
    # No empty or broken figure blocks are emitted.
    assert "section_2c_human_scatter_figure" not in roles
    assert "section_2c_human_heatmap_figure" not in roles
    assert not [b for b in result.blocks if b.kind == "figure" and not b.figure_path]
    notes = [b.text for b in result.blocks if b.presentation_role == "section_2c_source_status"]
    assert any("unavailable" in note for note in notes)


def test_unavailable_branch_gets_narrative_and_diagnostic(tmp_path):
    out = _run_section(tmp_path, skip={GEO_SOURCE_KEY}, run_id="diag")
    status = out["section_2c_status"]
    result = _blocks(status, out["evidence_records"])
    roles = [b.presentation_role for b in result.blocks]
    assert "section_2c_dropviz_narrative" in roles
    assert "section_2c_dropviz_table" not in roles
    assert "section_2c_dropviz_rank_figure" not in roles
    assert any(
        d.severity == "warning" and "DropViz" in d.reason for d in result.diagnostics
    )


def test_page_segments_split_on_declared_breaks(section_run, status):
    result = _blocks(status, section_run["evidence_records"])
    segments = split_section_2c_page_segments(list(result.blocks))
    assert len(segments) == 4
    assert segments[0][0].presentation_role == "section_2c_intro"
    assert segments[1][0].presentation_role == "section_2c_mouse_narrative"
    assert segments[2][0].presentation_role == "section_2c_dropviz_narrative"
    # Without renderable figure bytes the therapeutic narrative keeps its page.
    assert segments[3][0].presentation_role == "section_2c_therapeutic_narrative"


def test_ranking_figure_opens_its_own_page(section_run, status, monkeypatch):
    """The PDF engine shrinks a chart to whatever space is left on a page, so
    the ranking figure must start a page and let the narrative follow it."""
    # Presentation resolves figure paths through the process-wide settings root,
    # so point it at this run's temporary raw store.
    artifact_root = Path(status["audit"]["artifact_root"])
    monkeypatch.setenv("RAW_DATA_DIR", str(artifact_root))
    get_settings.cache_clear()
    try:
        result = build_section_2c_blocks(
            gene_symbol=status["summary"]["gene_symbol"],
            evidence_records=list(section_run["evidence_records"]),
            section_status=status,
        )
        segments = split_section_2c_page_segments(list(result.blocks))
        assert [b.presentation_role for b in segments[-1]] == [
            "section_2c_dropviz_rank_figure",
            "section_2c_geo_attribution",
            "section_2c_source_status",
            "section_2c_therapeutic_narrative",
        ]
    finally:
        get_settings.cache_clear()


def test_subsection_heading_appears_once(section_run, status):
    result = _blocks(status, section_run["evidence_records"])
    sub = ReportSubsection(
        key="c",
        title="snRNA-Seq gene expression in cell type database",
        toc_title="SNRNA-SEQ GENE EXPRESSION IN CELL TYPE DATABASE",
        presentation_blocks=list(result.blocks),
        status="populated",
    )
    segments = render_section_2c_subsection_segments(sub)
    joined = "\n".join(segments)
    assert joined.count("c. snRNA-Seq gene expression in cell type database") == 1
    assert 'class="sub-heading"' in segments[0]


def test_focused_2c_html_has_one_major_heading_and_no_1c_renderer(section_run, status):
    document, presentation, audit = build_section_bundle_document(
        dossier_run_id="render-2c",
        gene_symbol=status["summary"]["gene_symbol"],
        section_keys=["2c"],
        evidence_records=list(section_run["evidence_records"]),
        api_runs=list(section_run["api_runs"]),
        raw_artifacts=list(section_run["raw_artifacts"]),
        section_status_by_key={"2c": status},
    )
    html = render_section_bundle_html(document)
    assert html.count('class="major-heading"') == 1
    assert "2. Expression pattern by cell and tissue" in html
    assert html.count("c. snRNA-Seq gene expression in cell type database") == 1
    # Section 1c's grouped renderer must not claim the 2c subsection.
    body = html.split("</style>", 1)[1]
    assert "section-1c" not in body
    assert "section_1c" not in body
    assert 'class="report-page section-bundle-body section-2-page"' in html
    assert audit["section_2c_status"] is status
    blocks = presentation["major_sections"][0]["subsections"][0]["blocks"]
    assert [b["presentation_role"] for b in blocks][0] == "section_2c_intro"


def test_assembled_major_two_orders_2c_after_2b_with_a_page_break(section_run, status):
    """2c must start on a clean page and never repeat the Major 2 heading."""
    document, _presentation, _audit = build_section_bundle_document(
        dossier_run_id="render-2c-assembled",
        gene_symbol=status["summary"]["gene_symbol"],
        section_keys=["2c"],
        evidence_records=list(section_run["evidence_records"]),
        api_runs=list(section_run["api_runs"]),
        raw_artifacts=list(section_run["raw_artifacts"]),
        section_status_by_key={"2c": status},
    )
    html = render_section_bundle_html(document)
    assert html.count('class="major-heading"') == 1
    # Count page sections only — the stylesheet also mentions the class name.
    assert (
        html.count('class="report-page section-bundle-body section-2c-continuation"')
        == 3
    )


# ---------------------------------------------------------------------------
# Gene generality / production hygiene
# ---------------------------------------------------------------------------
def test_no_validation_gene_or_golden_conclusion_in_production_source():
    text = SRC.read_text(encoding="utf-8")
    for token in ("SREBF2", "Srebf2", "CDH10", "Cdh10", "SN_4-3", "371_Oligo"):
        assert token not in text, f"{token} must not appear in production source"


def test_production_source_hardcodes_no_golden_numbers():
    text = SRC.read_text(encoding="utf-8")
    for token in ("19.4", "97.35", "94.83", "4.63", "0.68", "565", "387", "127"):
        assert token not in text, f"golden value {token} must not be hardcoded"


def test_an_arbitrary_synthetic_gene_renders_end_to_end(tmp_path):
    """Nothing in Section 2c depends on a particular gene symbol."""
    human_csv = "feature,Astro_1,Vip_1,Oligo_1\nZZTESTGENE,1.25,0.0,3.5\n"
    mouse_csv = "feature,1_CR,2_Lamp5\nZztestgene,0.5,4.0\n"
    geo_csv = "gene,pop_X,pop_Y\nZztestgene,3,9\nOther,97,91\n"
    dend_h = _internal("", [_internal("Non-neuronal", [_leaf("Astro_1"), _leaf("Oligo_1")]), _leaf("Vip_1")])
    dend_m = _internal("", [_internal("GABAergic neurons", [_leaf("2_Lamp5")]), _leaf("1_CR")])

    paths = paths_for(tmp_path / "outputs")
    paths.ensure()
    for source_key, content, url in (
        (
            ac.CACHE_KEY_HUMAN_TRIMMED_MEANS,
            human_csv.encode(),
            ac.human_m1_source_url("trimmed_means"),
        ),
        (
            ac.CACHE_KEY_HUMAN_TAXONOMY,
            json.dumps(dend_h).encode(),
            ac.human_m1_source_url("taxonomy"),
        ),
        (ac.CACHE_KEY_MOUSE_TRIMMED_MEANS, mouse_csv.encode(), None),
        (ac.CACHE_KEY_MOUSE_TAXONOMY, json.dumps(dend_m).encode(), None),
        (GEO_SOURCE_KEY, _gzip(geo_csv), dg.supplementary_url()),
    ):
        _accept(paths, source_key=source_key, content=content, official_url=url)

    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "outputs")
    out = node_generate_section_2c_derived_artifacts(
        {
            "run_type": "section_bundle",
            "selected_section_keys": ["2c"],
            "dossier_run_id": "arbitrary-gene",
            "gene_symbol": "ZZTESTGENE",
            "gene_ids": {"mouse_symbol": "Zztestgene"},
            "evidence_records": [],
            "api_runs": [],
            "raw_artifacts": [],
            "errors": [],
            "coverage": [],
        },
        settings=settings,
        persist_db=False,
        config=Section2cConfig(
            output_root=tmp_path / "outputs", attempt_allen_figures=False
        ),
    )
    status = out["section_2c_status"]
    assert status["rendering_status"]["overall"] == STATUS_SUCCESS
    assert status["summary"]["mouse_symbol"] == "Zztestgene"
    assert status["summary"]["dropviz"]["top_populations"][0]["population_label"] == "pop_Y"

    result = _blocks(status, out["evidence_records"])
    roles = [b.presentation_role for b in result.blocks]
    assert "section_2c_human_table" in roles
    assert "section_2c_dropviz_table" in roles
    assert "ZZTESTGENE" in status["summary"]["therapeutic_narrative"]


def test_gate_skips_non_bundle_runs(tmp_path):
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "outputs")
    base = {
        "run_type": "full_dossier",
        "selected_section_keys": ["2c"],
        "dossier_run_id": "gated",
        "gene_symbol": HUMAN_SYMBOL,
        "gene_ids": {},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }
    out = node_generate_section_2c_derived_artifacts(
        dict(base), settings=settings, persist_db=False, config=Section2cConfig()
    )
    assert "section_2c_status" not in out

    unselected = {**base, "run_type": "section_bundle", "selected_section_keys": ["2b"]}
    out = node_generate_section_2c_derived_artifacts(
        unselected, settings=settings, persist_db=False, config=Section2cConfig()
    )
    assert "section_2c_status" not in out


def test_figures_are_not_attempted_when_disabled(status):
    rendering = status["rendering_status"]
    for key in (
        "allen_human_scatter_status",
        "allen_human_heatmap_status",
        "allen_mouse_scatter_status",
        "allen_mouse_heatmap_status",
    ):
        assert rendering[key] == STATUS_NOT_ATTEMPTED
    # Skipping an optional capture is not the same as it being unavailable.
    assert rendering["overall"] == STATUS_SUCCESS
