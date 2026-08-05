"""Offline tests for the DropViz GEO (GSE116470) processed-matrix client.

No network. Every fixture is synthetic and gene-general: the target symbol is
``Genex`` so nothing here depends on a validation gene.
"""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

import pytest

from gene_dossier.tools import dropviz_geo as dg

SRC = Path(__file__).resolve().parents[1] / "src" / "gene_dossier" / "tools" / "dropviz_geo.py"

# Canonical synthetic matrix: column totals are all 100, so Genex is
# 5,000 / 20,000 / 0 per 100,000 and the order is pop_B > pop_A > pop_C.
FIXTURE_CSV = """gene,pop_A,pop_B,pop_C
Gene1,10,0,5
Genex,5,20,0
Gene3,85,80,95
"""


def gzip_bytes(text: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as fh:
        fh.write(text.encode("utf-8"))
    return buf.getvalue()


def scan_text(text: str, *, target: str | None = "Genex") -> dg.MatrixScan:
    return dg.scan_matrix(io.StringIO(text), target_gene=target)


def counts_semantics(scan: dg.MatrixScan) -> dict:
    return dg.classify_value_semantics(scan)


# --------------------------------------------------------------------------------------
# Download validation
# --------------------------------------------------------------------------------------
def test_valid_gzip_payload_accepted():
    result = dg.validate_gzip_payload(gzip_bytes(FIXTURE_CSV * 200))
    assert result["ok"] is True
    assert result["error_type"] is None


def test_html_masquerading_as_data_rejected():
    body = b"<!DOCTYPE html><html><body>Not found</body></html>" + b" " * 2000
    result = dg.validate_gzip_payload(body)
    assert result["ok"] is False
    assert result["error_type"] == "html_masquerading_as_data"


def test_truncated_gzip_rejected():
    full = gzip_bytes(FIXTURE_CSV * 500)
    truncated = full[: len(full) // 2]
    result = dg.validate_gzip_payload(truncated)
    assert result["ok"] is False
    assert result["error_type"] == "truncated_download"


def test_non_gzip_payload_rejected():
    result = dg.validate_gzip_payload(b"gene,pop_A\n" + b"x" * 4000)
    assert result["ok"] is False
    assert result["error_type"] == "invalid_gzip_magic"


def test_empty_download_rejected():
    assert dg.validate_gzip_payload(b"")["error_type"] == "empty_download"


def test_official_urls_are_ncbi_hosted():
    for url in dg.candidate_download_urls():
        assert dg.is_allowlisted_url(url)
    assert dg.supplementary_url().endswith(dg.METACELL_FILENAME)


def test_off_allowlist_url_rejected():
    assert not dg.is_allowlisted_url("https://example.com/GSE116470.csv.gz")
    assert not dg.is_allowlisted_url("ftp://ftp.ncbi.nlm.nih.gov/x.gz")


def test_download_refuses_off_allowlist_url_without_network():
    tr = dg.download_metacell_matrix(
        gene_symbol="Genex", urls=["https://evil.example.com/matrix.csv.gz"]
    )
    assert tr.success is False
    assert tr.error_type == "url_not_allowlisted"


def test_open_matrix_stream_roundtrip():
    lines = list(dg.open_matrix_stream(gzip_bytes(FIXTURE_CSV)))
    assert lines[0].startswith("gene,pop_A")
    assert len(lines) == 4


# --------------------------------------------------------------------------------------
# Matrix scan / orientation
# --------------------------------------------------------------------------------------
def test_scan_resolves_gene_column_and_labels():
    scan = scan_text(FIXTURE_CSV)
    assert scan.ok is True
    assert scan.gene_column_index == 0
    assert scan.gene_column_status == "resolved_by_header"
    assert scan.population_labels == ["pop_A", "pop_B", "pop_C"]
    assert scan.gene_row_count == 3


def test_streaming_column_totals():
    scan = scan_text(FIXTURE_CSV)
    assert scan.column_totals == [100.0, 100.0, 100.0]


def test_exact_gene_match_only():
    scan = scan_text(FIXTURE_CSV)
    assert scan.target_matches == 1
    assert scan.target_row == [5.0, 20.0, 0.0]


def test_no_substring_match():
    text = "gene,pop_A,pop_B\nGenex-AS1,7,7\nGenexP1,3,3\nOther,90,90\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    assert scan.target_matches == 0


def test_duplicate_exact_target_rows_detected():
    text = "gene,pop_A,pop_B\nGenex,5,5\nGenex,6,6\nOther,89,89\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    assert scan.target_matches == 2
    out = dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS
    )
    assert out["status"] == "target_gene_match_failed"


def test_ambiguous_orientation_fails_closed():
    text = "0.1,0.2,0.3\n1,2,3\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    assert scan.ok is False
    assert scan.error_type == "ambiguous_orientation"


def test_empty_matrix_fails_closed():
    scan = dg.scan_matrix(io.StringIO(""), target_gene="Genex")
    assert scan.ok is False
    assert scan.error_type == "empty_matrix"


def test_malformed_rows_counted_not_fatal():
    text = "gene,pop_A,pop_B\nGene1,1\nGenex,5,20\nGene3,94,80\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    assert scan.malformed_row_count == 1
    assert scan.target_matches == 1


def test_malformed_values_recorded():
    text = "gene,pop_A,pop_B\nGenex,abc,20\nOther,100,80\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    assert scan.malformed_value_count == 1


def test_duplicate_population_labels_detected():
    text = "gene,pop_A,pop_A,pop_B\nGenex,5,5,20\nOther,95,95,80\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    assert scan.duplicate_population_labels == ["pop_A"]


# --------------------------------------------------------------------------------------
# Value-semantics gate (the critical scientific checkpoint)
# --------------------------------------------------------------------------------------
def test_integral_nonnegative_classified_as_raw_counts():
    semantics = counts_semantics(scan_text(FIXTURE_CSV))
    assert semantics["value_semantics_status"] == dg.VALUE_SEMANTICS_RAW_COUNTS
    assert semantics["basis"] == "integral_nonnegative_values"


def test_small_integer_matrix_is_count_compatible_not_raw_counts():
    text = "gene,pop_A,pop_B\nGenex,1,2\nOther,3,4\n"
    semantics = counts_semantics(dg.scan_matrix(io.StringIO(text), target_gene="Genex"))
    assert semantics["value_semantics_status"] == dg.VALUE_SEMANTICS_COUNT_COMPATIBLE
    assert semantics["value_semantics_status"] in dg.COUNT_COMPATIBLE_SEMANTICS


def test_negative_values_classified_as_transformed():
    text = "gene,pop_A,pop_B\nGenex,-1.5,2.5\nOther,3.5,4.5\n"
    semantics = counts_semantics(dg.scan_matrix(io.StringIO(text), target_gene="Genex"))
    assert semantics["value_semantics_status"] == dg.VALUE_SEMANTICS_TRANSFORMED
    assert semantics["basis"] == "negative_values_present"


def test_nonnegative_decimals_do_not_become_raw_counts():
    # Explicitly guards the rule: nonnegativity alone must never imply counts.
    text = "gene,pop_A,pop_B\nGenex,0.31,0.72\nOther,1.24,2.51\n"
    semantics = counts_semantics(dg.scan_matrix(io.StringIO(text), target_gene="Genex"))
    assert semantics["value_semantics_status"] != dg.VALUE_SEMANTICS_RAW_COUNTS
    assert semantics["value_semantics_status"] == dg.VALUE_SEMANTICS_UNRESOLVED
    assert semantics["evidence"]["negative_value_count"] == 0


def test_near_constant_column_sums_classified_as_normalized():
    # Each column sums to 1.0 (a per-population normalization).
    text = "gene,pop_A,pop_B\nGenex,0.25,0.5\nOther,0.75,0.5\n"
    semantics = counts_semantics(dg.scan_matrix(io.StringIO(text), target_gene="Genex"))
    assert semantics["value_semantics_status"] == dg.VALUE_SEMANTICS_NORMALIZED
    assert semantics["basis"] == "near_constant_column_sums"


def test_documented_semantics_override_is_recorded():
    semantics = dg.classify_value_semantics(
        scan_text(FIXTURE_CSV),
        documented_semantics=dg.VALUE_SEMANTICS_NORMALIZED,
        documentation_reference="GEO GSE116470 series matrix description",
        documented_unit="CPM",
    )
    assert semantics["value_semantics_status"] == dg.VALUE_SEMANTICS_NORMALIZED
    assert semantics["basis"] == "source_documentation"
    assert semantics["documented_unit"] == "CPM"
    assert (
        semantics["evidence"]["source_documentation_reference"]
        == "GEO GSE116470 series matrix description"
    )


def test_profile_records_classification_evidence():
    scan = scan_text(FIXTURE_CSV)
    profile = dg.build_matrix_profile(
        scan,
        semantics=counts_semantics(scan),
        source_sha256="abc123",
        source_url=dg.supplementary_url(),
    )
    evidence = profile["value_semantics_evidence"]
    for key in (
        "integer_fraction",
        "minimum",
        "maximum",
        "decimal_value_count",
        "negative_value_count",
        "column_sum_distribution",
    ):
        assert key in evidence
    assert profile["source_filename"] == dg.METACELL_FILENAME
    assert profile["source_sha256"] == "abc123"
    assert profile["calculation_version"] == dg.CALCULATION_VERSION
    # Profile must be JSON-serializable for the artifact write.
    json.dumps(profile)


# --------------------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------------------
def build_counts_ranking():
    scan = scan_text(FIXTURE_CSV)
    return scan, dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS
    )


def test_transcripts_per_100k_calculation():
    _, out = build_counts_ranking()
    assert out["status"] == dg.RANK_STATUS_SUCCESS
    by_label = {r["population_label"]: r for r in out["records"]}
    assert by_label["pop_A"]["transcripts_per_100k"] == pytest.approx(5000.0)
    assert by_label["pop_B"]["transcripts_per_100k"] == pytest.approx(20000.0)
    assert by_label["pop_C"]["transcripts_per_100k"] == pytest.approx(0.0)


def test_raw_order_preserved_in_records():
    _, out = build_counts_ranking()
    assert [r["population_label"] for r in out["records"]] == ["pop_A", "pop_B", "pop_C"]


def test_descending_ranking_order():
    _, out = build_counts_ranking()
    ranked = dg.rank_records(out["records"])
    assert [r["population_label"] for r in ranked] == ["pop_B", "pop_A", "pop_C"]


def test_ranking_record_metric_and_unit_for_counts():
    _, out = build_counts_ranking()
    rec = out["records"][0]
    assert rec["ranking_metric"] == dg.RANKING_METRIC_PER_100K
    assert rec["expression_unit"] == dg.EXPRESSION_UNIT_PER_100K
    assert rec["population_total"] == 100.0
    assert out["presentation"]["axis_label"] == dg.AXIS_LABEL_PER_100K


def test_normalized_ranking_preserves_source_scale_and_unit():
    text = "gene,pop_A,pop_B\nGenex,0.25,0.5\nOther,0.75,0.5\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    out = dg.build_ranking_records(
        scan,
        value_semantics_status=dg.VALUE_SEMANTICS_NORMALIZED,
        documented_unit="CPM",
    )
    assert out["status"] == dg.RANK_STATUS_SUCCESS
    rec = {r["population_label"]: r for r in out["records"]}["pop_B"]
    assert rec["ranking_value"] == pytest.approx(0.5)
    assert rec["transcripts_per_100k"] is None
    assert rec["expression_unit"] == "CPM"
    assert out["presentation"]["axis_label"] == "CPM"


def test_normalized_without_documented_unit_is_unresolved():
    text = "gene,pop_A,pop_B\nGenex,0.25,0.5\nOther,0.75,0.5\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    out = dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_NORMALIZED
    )
    assert out["status"] == dg.RANK_STATUS_NORMALIZATION_UNRESOLVED
    assert out["records"] == []


@pytest.mark.parametrize(
    "status", [dg.VALUE_SEMANTICS_TRANSFORMED, dg.VALUE_SEMANTICS_UNRESOLVED]
)
def test_unresolved_semantics_produce_no_ranking(status):
    scan = scan_text(FIXTURE_CSV)
    out = dg.build_ranking_records(scan, value_semantics_status=status)
    assert out["status"] == dg.RANK_STATUS_NORMALIZATION_UNRESOLVED
    assert out["records"] == []
    assert out["presentation"]["chartable"] is False


def test_zero_total_population_excluded():
    text = "gene,pop_A,pop_Z\nGenex,5,0\nOther,95,0\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    out = dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS
    )
    reasons = {e["label"]: e["reason"] for e in out["excluded"]}
    assert reasons["pop_Z"] == "population_total_not_positive"
    assert [r["population_label"] for r in out["records"]] == ["pop_A"]


def test_duplicate_label_population_excluded_with_reason():
    text = "gene,pop_A,pop_A,pop_B\nGenex,5,5,20\nOther,95,95,80\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    out = dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS
    )
    reasons = {e["reason"] for e in out["excluded"]}
    assert "duplicate_population_label" in reasons
    assert [r["population_label"] for r in out["records"]] == ["pop_B"]


def test_malformed_target_value_excluded():
    text = "gene,pop_A,pop_B\nGenex,abc,20\nOther,100,80\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    out = dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS
    )
    reasons = {e["label"]: e["reason"] for e in out["excluded"]}
    assert reasons["pop_A"] == "malformed_target_value"


# --------------------------------------------------------------------------------------
# Summary metrics
# --------------------------------------------------------------------------------------
def test_nonzero_counts_and_percentage():
    _, out = build_counts_ranking()
    summary = dg.summarize_ranking(
        out["records"], value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS
    )
    assert summary["valid_population_count"] == 3
    assert summary["nonzero_population_count"] == 2
    assert summary["nonzero_percentage"] == pytest.approx(66.67, abs=0.01)
    # "detected" is reserved for a source-defined threshold.
    assert "detected_count" not in summary
    assert "detected_percentage" not in summary


def test_share_metrics_are_labelled_distinctly():
    _, out = build_counts_ranking()
    summary = dg.summarize_ranking(
        out["records"], value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS, top_n=1
    )
    # Top 1 of 3: normalized share = 20000 / 25000 = 80%; raw share = 20 / 25 = 80%.
    assert summary["top_10_normalized_expression_share"] == pytest.approx(80.0)
    assert summary["top_10_raw_target_count_share"] == pytest.approx(80.0)
    assert "normalized population-level expression signal" in (
        summary["top_10_normalized_expression_share_description"]
    )


def test_raw_share_omitted_for_normalized_semantics():
    text = "gene,pop_A,pop_B\nGenex,0.25,0.5\nOther,0.75,0.5\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    out = dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_NORMALIZED, documented_unit="CPM"
    )
    summary = dg.summarize_ranking(
        out["records"], value_semantics_status=dg.VALUE_SEMANTICS_NORMALIZED
    )
    assert summary["top_10_raw_target_count_share"] is None
    assert summary["top_10_normalized_expression_share"] is not None
    assert summary["expression_unit"] == "CPM"


def test_max_to_median_ratio():
    _, out = build_counts_ranking()
    summary = dg.summarize_ranking(
        out["records"], value_semantics_status=dg.VALUE_SEMANTICS_RAW_COUNTS
    )
    # values 20000, 5000, 0 -> median 5000
    assert summary["max_to_median_ratio"] == pytest.approx(4.0)


# --------------------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------------------
def test_population_label_parsing_extracts_encoded_fields_only():
    parsed = dg.parse_population_label("FC_1-2")
    assert parsed["region_abbreviation"] == "FC"
    assert parsed["cluster_number"] == 1
    assert parsed["subcluster_number"] == 2
    # Broad cell class is never inferred from biological knowledge.
    assert parsed["broad_cell_class"] is None
    assert parsed["neuronal_or_glial"] is None


def test_unstructured_label_is_unresolved():
    parsed = dg.parse_population_label("mystery")
    assert parsed["region_abbreviation"] is None
    assert dg.label_mapping_status([parsed]) == dg.LABEL_MAPPING_UNRESOLVED


def test_structured_labels_report_partial_mapping():
    parsed = [dg.parse_population_label(x) for x in ("FC_1-1", "HC_2-3")]
    assert dg.label_mapping_status(parsed) == dg.LABEL_MAPPING_PARTIAL


def test_empty_label_set_is_unresolved():
    assert dg.label_mapping_status([]) == dg.LABEL_MAPPING_UNRESOLVED


# --------------------------------------------------------------------------------------
# Chart
# --------------------------------------------------------------------------------------
def test_chart_renders_png_point_only_when_ci_unresolved():
    _, out = build_counts_ranking()
    ranked = dg.rank_records(out["records"])
    png = dg.render_top_populations_png(
        mouse_gene_symbol="Genex",
        top_populations=ranked,
        axis_label=out["presentation"]["axis_label"],
    )
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert len(png) > 1000


def test_chart_uses_documented_unit_axis_for_normalized_data():
    text = "gene,pop_A,pop_B\nGenex,0.25,0.5\nOther,0.75,0.5\n"
    scan = dg.scan_matrix(io.StringIO(text), target_gene="Genex")
    out = dg.build_ranking_records(
        scan, value_semantics_status=dg.VALUE_SEMANTICS_NORMALIZED, documented_unit="CPM"
    )
    assert out["presentation"]["axis_label"] == "CPM"
    png = dg.render_top_populations_png(
        mouse_gene_symbol="Genex",
        top_populations=dg.rank_records(out["records"]),
        axis_label=out["presentation"]["axis_label"],
    )
    assert png.startswith(b"\x89PNG\r\n\x1a\n")


# --------------------------------------------------------------------------------------
# Production hygiene
# --------------------------------------------------------------------------------------
def test_no_validation_gene_hardcoded_in_module():
    text = SRC.read_text(encoding="utf-8")
    for token in ("SREBF2", "Srebf2", "CDH10", "Cdh10"):
        assert token not in text, f"{token} must not appear in production source"


def test_module_declares_no_api_run_for_local_derivations():
    # Local derivation helpers must not fabricate HTTP provenance.
    text = SRC.read_text(encoding="utf-8")
    assert "ApiRun(" not in text
