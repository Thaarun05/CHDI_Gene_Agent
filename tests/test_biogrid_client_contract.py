"""Regression tests for the pre-Section-5b BioGRID client public contract."""

from __future__ import annotations

from unittest.mock import patch

from gene_dossier.config import Settings
from gene_dossier.models import ToolResult
from gene_dossier.tools import biogrid as bg


def _mock_client(payload, *, status_code: int = 200, is_success: bool = True):
    import json as _json

    body = _json.dumps(payload).encode("utf-8") if not isinstance(payload, bytes) else payload

    class _Resp:
        def __init__(self):
            self.status_code = status_code
            self.is_success = is_success
            self.url = "https://webservice.thebiogrid.org/interactions/"
            self.content = body if isinstance(payload, (dict, list)) else body
            self.text = body.decode("utf-8", errors="replace")

        def json(self):
            if isinstance(payload, (dict, list)):
                return payload
            raise ValueError("not json")

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            _Client.last_params = dict(params or {})
            _Client.last_url = url
            return _Resp()

    _Client.last_params = {}
    _Client.last_url = ""
    return _Client


def test_build_interaction_params_legacy_defaults() -> None:
    params = bg.build_interaction_params("SREBF2", accesskey="SECRETKEY")
    assert params["searchNames"] == "true"
    assert params["geneList"] == "SREBF2"
    assert params["taxId"] == "9606"
    assert params["includeInteractors"] == "true"
    assert params["includeInteractorInteractions"] == "false"
    assert params["selfInteractionsExcluded"] == "true"
    assert params["interSpeciesExcluded"] == "true"
    assert params["format"] == "jsonExtended"
    assert params["max"] == "10000"
    assert params["accesskey"] == "SECRETKEY"


def test_safe_params_redacts_accesskey() -> None:
    safe = bg._safe_params({"geneList": "SREBF2", "accesskey": "SECRETKEY"})
    assert safe["accesskey"] == "***"
    assert "SECRETKEY" not in str(safe)


def test_interactions_requires_key() -> None:
    settings = Settings(biogrid_accesskey=None)
    result = bg.interactions("SREBF2", settings=settings)
    assert result.success is False
    assert result.error_type == "requires_key"
    assert result.source_name == bg.SOURCE_NAME
    assert result.endpoint_name == "interactions"


def test_interactions_success_outer_shape_and_legacy_filters() -> None:
    payload = {
        "1": {
            "BIOGRID_INTERACTION_ID": 1,
            "ENTREZ_GENE_A": "6721",
            "ENTREZ_GENE_B": "100",
            "BIOGRID_ID_A": 112599,
            "BIOGRID_ID_B": 1000,
            "OFFICIAL_SYMBOL_A": "SREBF2",
            "OFFICIAL_SYMBOL_B": "PARTNER",
            "EXPERIMENTAL_SYSTEM": "Two-hybrid",
            "EXPERIMENTAL_SYSTEM_TYPE": "physical",
            "ORGANISM_A": 9606,
            "ORGANISM_B": 9606,
        }
    }
    Client = _mock_client(payload)
    settings = Settings(biogrid_accesskey="SECRETKEY")
    with patch("gene_dossier.tools.biogrid.httpx.Client", Client):
        result = bg.interactions("SREBF2", settings=settings)
    assert result.success is True
    assert result.endpoint_name == "interactions"
    assert result.source_name == "BioGRID"
    assert result.data == payload
    assert Client.last_params["selfInteractionsExcluded"] == "true"
    assert Client.last_params["interSpeciesExcluded"] == "true"
    assert "SECRETKEY" not in (result.request_url or "")
    assert "accesskey=" in (result.request_url or "")
    assert result.request_params.get("accesskey") == "***"


def test_fetch_interactions_outer_data_keys() -> None:
    payload = {
        "9": {
            "BIOGRID_INTERACTION_ID": 9,
            "ENTREZ_GENE_A": "6721",
            "ENTREZ_GENE_B": "200",
            "BIOGRID_ID_A": 112599,
            "BIOGRID_ID_B": 2000,
            "OFFICIAL_SYMBOL_A": "SREBF2",
            "OFFICIAL_SYMBOL_B": "OTHER",
            "EXPERIMENTAL_SYSTEM": "Affinity Capture-MS",
            "EXPERIMENTAL_SYSTEM_TYPE": "physical",
            "PUBMED_ID": 1,
            "ORGANISM_A": 9606,
            "ORGANISM_B": 9606,
            "THROUGHPUT": "High Throughput",
            "SOURCE_DATABASE": "BIOGRID",
        }
    }
    Client = _mock_client(payload)
    settings = Settings(biogrid_accesskey="SECRETKEY")
    with patch("gene_dossier.tools.biogrid.httpx.Client", Client):
        result = bg.fetch_interactions("SREBF2", settings=settings)
    assert result.success is True
    assert result.endpoint_name == "fetch_interactions"
    assert isinstance(result.data, dict)
    assert set(result.data) >= {
        "gene_symbol",
        "tax_id",
        "interactions",
        "interaction_rows",
        "interaction_summaries",
        "interaction_count",
    }
    assert result.data["gene_symbol"] == "SREBF2"
    assert result.data["tax_id"] == 9606
    assert result.data["interaction_count"] == 1
    assert len(result.data["interaction_rows"]) == 1
    assert len(result.data["interaction_summaries"]) == 1
    summary = result.data["interaction_summaries"][0]
    assert summary["biogrid_interaction_id"] == 9
    assert summary["official_symbol_a"] == "SREBF2"


def test_interactions_as_list_dict_and_list() -> None:
    as_list = [{"BIOGRID_INTERACTION_ID": 1, "OFFICIAL_SYMBOL_A": "A"}]
    assert len(bg.interactions_as_list(as_list)) == 1
    as_dict = {"1": {"BIOGRID_INTERACTION_ID": 1, "OFFICIAL_SYMBOL_A": "A"}}
    rows = bg.interactions_as_list(as_dict)
    assert len(rows) == 1
    assert rows[0]["BIOGRID_INTERACTION_ID"] == 1


def test_fetch_interactions_propagates_failure() -> None:
    settings = Settings(biogrid_accesskey=None)
    result = bg.fetch_interactions("SREBF2", settings=settings)
    assert result.success is False
    assert result.endpoint_name == "fetch_interactions"
    assert "interactions" in (result.data or {})
