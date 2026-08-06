"""Regression tests for the pre-Section-4a Harmonizome client public contract."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from gene_dossier.models import ToolResult
from gene_dossier.tools import harmonizome as hz


def test_fetch_tf_associations_outer_data_keys_legacy_payload() -> None:
    """Existing callers expect these top-level data keys and success semantics."""
    legacy_payload = {
        "symbol": "GENEX",
        "name": "Gene X",
        "synonyms": ["GX"],
        "associations": [
            {
                "gene": {"symbol": "TF1", "name": "TF one"},
                "dataset": {"name": "ENCODE Transcription Factor Targets"},
                "attribute": {"name": "TF1", "href": "/api/1.0/gene_set/TF1/ENCODE"},
                "thresholdValue": 1.0,
                "standardizedValue": 0.5,
            },
            {
                "gene": {"symbol": "OTHER"},
                "dataset": {"name": "Some Other Dataset"},
                "attribute": {"name": "OTHER"},
            },
        ],
    }
    fake = ToolResult(
        source_name=hz.SOURCE_NAME,
        endpoint_name="gene_associations",
        gene_symbol="GENEX",
        request_url="https://maayanlab.cloud/Harmonizome/api/1.0/gene/GENEX?showAssociations=true",
        request_params={"gene_symbol": "GENEX", "showAssociations": "true"},
        success=True,
        status_code=200,
        data=legacy_payload,
    )
    with patch.object(hz, "gene_associations", return_value=fake):
        result = hz.fetch_tf_associations("GENEX")
    assert result.success is True
    assert result.endpoint_name == "fetch_gene_associations"
    assert isinstance(result.data, dict)
    for key in (
        "gene_symbol",
        "name",
        "synonyms",
        "associations",
        "association_summaries",
        "tf_associations",
        "tf_summaries",
        "association_count",
        "tf_count",
        "tf_only",
        "raw",
    ):
        assert key in result.data
    assert result.data["tf_only"] is True
    assert result.data["tf_count"] == 1
    assert result.data["association_count"] == 1
    assert result.data["gene_symbol"] == "GENEX"


def test_fetch_gene_associations_failure_preserves_error_semantics() -> None:
    failed = ToolResult(
        source_name=hz.SOURCE_NAME,
        endpoint_name="gene_associations",
        gene_symbol="GENEX",
        request_url="https://example.invalid",
        request_params={},
        success=False,
        status_code=500,
        data={"raw_text": "boom"},
        error_type="http_error",
        error_message="HTTP 500",
    )
    with patch.object(hz, "gene_associations", return_value=failed):
        result = hz.fetch_gene_associations("GENEX")
    assert result.success is False
    assert result.error_type == "http_error"
    assert isinstance(result.data, dict)
    assert "raw" in result.data


def test_summarize_association_legacy_keys() -> None:
    row = {
        "gene": {"symbol": "A", "name": "Alpha"},
        "dataset": {"name": "ENCODE Transcription Factor Targets"},
        "attribute": {"name": "A", "href": "/x"},
        "thresholdValue": 1,
        "standardizedValue": 2,
    }
    summary = hz.summarize_association(row)
    assert summary["associated_gene_symbol"] == "A"
    assert summary["dataset_name"] == "ENCODE Transcription Factor Targets"
    assert summary["attribute_name"] == "A"
    assert summary["threshold_value"] == 1
