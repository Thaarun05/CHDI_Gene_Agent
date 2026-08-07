"""Regression tests for the pre-Section-5a STRING client public contract."""

from __future__ import annotations

from unittest.mock import patch

from gene_dossier.config import Settings
from gene_dossier.models import ToolResult
from gene_dossier.tools import string_db as sd


def test_prefer_string_id_exact_preferred_name() -> None:
    rows = [
        {"stringId": "9606.ENSP0000001", "preferredName": "OTHER"},
        {"stringId": "9606.ENSP0000002", "preferredName": "GENEX"},
    ]
    assert sd.prefer_string_id(rows, "GENEX") == "9606.ENSP0000002"
    assert sd.prefer_string_id(rows, "genex") == "9606.ENSP0000002"
    assert sd.prefer_string_id([], "GENEX") is None


def test_get_string_ids_success_outer_shape() -> None:
    import json as _json

    raw_rows = [
        {
            "queryItem": "GENEX",
            "preferredName": "GENEX",
            "stringId": "9606.ENSP0000999",
            "ncbiTaxonId": 9606,
            "taxonName": "Homo sapiens",
        }
    ]
    settings = Settings()
    body = _json.dumps(raw_rows).encode("utf-8")

    class _Resp:
        status_code = 200
        is_success = True
        url = "https://example/get_string_ids"
        history: list = []
        headers = {"content-type": "application/json"}
        content = body

        def json(self):
            return raw_rows

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            return _Resp()

    with patch("gene_dossier.tools.string_db.httpx.Client", _Client):
        result = sd.get_string_ids("GENEX", settings=settings)
    assert result.success is True
    assert result.source_name == sd.SOURCE_NAME
    assert result.endpoint_name == "get_string_ids"
    assert isinstance(result.data, dict)
    assert result.data["gene_symbol"] == "GENEX"
    assert result.data["species"] == 9606
    assert result.data["string_id"] == "9606.ENSP0000999"
    assert "raw" in result.data


def test_interaction_partners_passes_limit_and_required_score() -> None:
    captured: dict = {}

    class _Resp:
        status_code = 200
        is_success = True
        url = "https://example/interaction_partners"
        history: list = []
        headers = {"content-type": "application/json"}
        content = b"[]"

        def json(self):
            return []

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            captured["params"] = dict(params or {})
            return _Resp()

    with patch("gene_dossier.tools.string_db.httpx.Client", _Client):
        result = sd.interaction_partners(
            "9606.ENSP0000999",
            gene_symbol="GENEX",
            limit=25,
            required_score=700,
            settings=Settings(),
        )
    assert result.success is True
    assert result.endpoint_name == "interaction_partners"
    assert captured["params"]["limit"] == "25"
    assert captured["params"]["required_score"] == "700"
    assert captured["params"]["network_type"] == "functional"
    assert captured["params"]["identifiers"] == "9606.ENSP0000999"


def test_fetch_interaction_partners_outer_data_keys() -> None:
    resolve = ToolResult(
        source_name=sd.SOURCE_NAME,
        endpoint_name="get_string_ids",
        gene_symbol="GENEX",
        request_url="https://example/ids",
        request_params={},
        success=True,
        status_code=200,
        data={
            "gene_symbol": "GENEX",
            "species": 9606,
            "string_id": "9606.ENSP0000999",
            "raw": [],
        },
    )
    partners = ToolResult(
        source_name=sd.SOURCE_NAME,
        endpoint_name="interaction_partners",
        gene_symbol="GENEX",
        request_url="https://example/partners",
        request_params={
            "identifiers": "9606.ENSP0000999",
            "species": "9606",
            "limit": "100",
            "required_score": "400",
            "network_type": "functional",
        },
        success=True,
        status_code=200,
        data=[
            {
                "stringId_A": "9606.ENSP0000999",
                "stringId_B": "9606.ENSP0000001",
                "preferredName_A": "GENEX",
                "preferredName_B": "PART",
                "score": 0.9,
            }
        ],
    )
    with patch.object(sd, "get_string_ids", return_value=resolve), patch.object(
        sd, "interaction_partners", return_value=partners
    ) as mock_partners:
        result = sd.fetch_interaction_partners("GENEX", limit=50, required_score=400)
    assert result.success is True
    assert result.endpoint_name == "fetch_interaction_partners"
    assert isinstance(result.data, dict)
    for key in (
        "gene_symbol",
        "species",
        "string_id",
        "get_string_ids",
        "partners",
        "partner_count",
    ):
        assert key in result.data
    assert result.data["gene_symbol"] == "GENEX"
    assert result.data["string_id"] == "9606.ENSP0000999"
    assert result.data["partner_count"] == 1
    assert isinstance(result.data["partners"], list)
    mock_partners.assert_called_once()
    assert mock_partners.call_args.kwargs["limit"] == 50
    assert mock_partners.call_args.kwargs["required_score"] == 400


def test_fetch_interaction_partners_resolve_failure_preserves_error() -> None:
    failed = ToolResult(
        source_name=sd.SOURCE_NAME,
        endpoint_name="get_string_ids",
        gene_symbol="GENEX",
        request_url="https://example/ids",
        request_params={},
        success=False,
        status_code=404,
        data=[],
        error_type="no_results",
        error_message="No STRING ID for GENEX",
    )
    with patch.object(sd, "get_string_ids", return_value=failed):
        result = sd.fetch_interaction_partners("GENEX")
    assert result.success is False
    assert result.endpoint_name == "fetch_interaction_partners"
    assert result.error_type in {"no_results", "resolve_failed"}
    assert isinstance(result.data, dict)
    assert "get_string_ids" in result.data
