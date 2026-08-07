"""Offline tests for the version-pinned STRING client and network helpers."""

from __future__ import annotations

import io
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from PIL import Image

from gene_dossier.config import Settings
from gene_dossier.models import ToolResult
from gene_dossier.section_5a import (
    Section5aConfig,
    canonicalize_network,
    node_generate_section_5a_derived_artifacts,
    parse_network_edge,
)
from gene_dossier.tools import string_db as sd

FIXTURES = Path(__file__).parent / "fixtures" / "string_db"
QUERY_ID = "9606.ENSP00000354476"


def _load_json(name: str):
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


def _mock_client(body: bytes, *, content_type: str = "application/json", url: str = "https://example"):
    class _Resp:
        status_code = 200
        is_success = True
        history: list = []
        headers = {"content-type": content_type}
        content = body

        def __init__(self):
            self.url = url

        def json(self):
            return json.loads(body.decode("utf-8"))

    class _Client:
        def __init__(self, *a, **k):
            self.captured: dict = {}

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, request_url, params=None):
            self.captured["url"] = request_url
            self.captured["params"] = dict(params or {})
            _Client.last_captured = self.captured  # type: ignore[attr-defined]
            return _Resp()

    return _Client


def test_string_base_url_is_version_12_0() -> None:
    assert "version-12-0" in sd.STRING_BASE_URL
    assert sd.STRING_BASE_URL == "https://version-12-0.string-db.org"
    assert sd.STRING_VERSION == "12.0"


@pytest.mark.parametrize(
    "value,expected",
    [
        (9606, 9606),
        ("9606", 9606),
        (" 9606 ", 9606),
        (None, None),
        ("", None),
        ("human", None),
        (9606.5, None),
        (True, None),
        (0, None),
        ("9606.0", None),
    ],
)
def test_normalize_taxon_id(value, expected) -> None:
    assert sd.normalize_taxon_id(value) == expected


def test_resolve_string_identifier_not_found() -> None:
    body = json.dumps(_load_json("empty_get_string_ids.json")).encode("utf-8")
    Client = _mock_client(body)
    with patch("gene_dossier.tools.string_db.httpx.Client", Client):
        result = sd.resolve_string_identifier("MISSING", settings=Settings())
    assert result.success is False
    assert result.error_type == "not_found"
    assert (result.data or {}).get("identifier_status") == "not_found"


def test_resolve_string_identifier_ambiguous() -> None:
    body = json.dumps(_load_json("ambiguous_get_string_ids.json")).encode("utf-8")
    Client = _mock_client(body)
    with patch("gene_dossier.tools.string_db.httpx.Client", Client):
        result = sd.resolve_string_identifier("GENEX", settings=Settings())
    assert result.success is False
    assert result.error_type == "ambiguous"
    assert (result.data or {}).get("identifier_status") == "ambiguous"


def test_resolve_string_identifier_success_from_fixture() -> None:
    body = json.dumps(_load_json("get_string_ids_srebf2.json")).encode("utf-8")
    Client = _mock_client(body)
    with patch("gene_dossier.tools.string_db.httpx.Client", Client):
        result = sd.resolve_string_identifier("SREBF2", settings=Settings())
    assert result.success is True
    assert result.data["string_id"] == QUERY_ID
    assert result.data["identifier_status"] == "resolved"
    assert result.data["ncbi_taxon_id"] == 9606


def test_canonicalize_network_fixture_topology_and_string_taxon() -> None:
    rows = _load_json("network_srebf2.json")
    assert all(isinstance(r["ncbiTaxonId"], str) and r["ncbiTaxonId"] == "9606" for r in rows)
    canon = canonicalize_network(
        rows,
        query_string_id=QUERY_ID,
        species_taxon_id=9606,
        required_score=400,
    )
    stats = canon["stats"]
    assert stats["unique_node_count"] == 31
    assert stats["unique_edge_count"] == 238
    assert stats["direct_query_edge_count"] == 30
    assert stats["neighbor_neighbor_edge_count"] == 208
    assert stats["min_combined_score"] == 0.4
    assert stats["max_combined_score"] == 0.999
    # parse path stores integer taxon via normalize_taxon_id
    assert all(e["ncbi_taxon_id"] == 9606 for e in canon["canonical_edges"])


def test_parse_network_edge_self_edge_and_duplicate_undirected() -> None:
    self_row = {
        "stringId_A": QUERY_ID,
        "stringId_B": QUERY_ID,
        "preferredName_A": "SREBF2",
        "preferredName_B": "SREBF2",
        "ncbiTaxonId": "9606",
        "score": 0.9,
    }
    rec, warn = parse_network_edge(
        self_row,
        source_order=1,
        query_string_id=QUERY_ID,
        species_taxon_id=9606,
        required_score=400,
    )
    assert rec is None
    assert warn == "self-edge rejected"

    partner = "9606.ENSP00000000001"
    base = {
        "stringId_A": QUERY_ID,
        "stringId_B": partner,
        "preferredName_A": "SREBF2",
        "preferredName_B": "PART",
        "ncbiTaxonId": "9606",
        "score": 0.8,
        "nscore": 0.1,
        "fscore": 0.0,
        "pscore": 0.0,
        "ascore": 0.0,
        "escore": 0.5,
        "dscore": 0.0,
        "tscore": 0.2,
    }
    swapped = {
        **base,
        "stringId_A": partner,
        "stringId_B": QUERY_ID,
        "preferredName_A": "PART",
        "preferredName_B": "SREBF2",
    }
    canon = canonicalize_network(
        [base, swapped],
        query_string_id=QUERY_ID,
        species_taxon_id=9606,
        required_score=400,
    )
    assert canon["stats"]["unique_edge_count"] == 1
    assert canon["stats"]["duplicate_undirected_edge_count"] == 1
    assert any("duplicate undirected" in w for w in canon["warnings"])


def test_canonicalize_below_threshold_warning() -> None:
    row = {
        "stringId_A": QUERY_ID,
        "stringId_B": "9606.ENSP00000000002",
        "preferredName_A": "SREBF2",
        "preferredName_B": "LOW",
        "ncbiTaxonId": "9606",
        "score": 0.3,
    }
    canon = canonicalize_network(
        [row],
        query_string_id=QUERY_ID,
        species_taxon_id=9606,
        required_score=400,
    )
    assert canon["stats"]["unique_edge_count"] == 0
    assert any("below threshold" in w for w in canon["warnings"])


def test_required_score_700_propagates_to_network_image_and_link() -> None:
    captured: list[dict] = []

    class _Resp:
        status_code = 200
        is_success = True
        url = "https://version-12-0.string-db.org/api/json/network"
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
            captured.append({"url": url, "params": dict(params or {})})
            return _Resp()

    settings = Settings()
    with patch("gene_dossier.tools.string_db.httpx.Client", _Client):
        sd.fetch_network(
            QUERY_ID,
            gene_symbol="SREBF2",
            required_score=700,
            settings=settings,
        )
        sd.fetch_network_image(
            QUERY_ID,
            gene_symbol="SREBF2",
            required_score=700,
            settings=settings,
        )
        sd.fetch_network_link(
            QUERY_ID,
            gene_symbol="SREBF2",
            required_score=700,
            settings=settings,
        )
    assert len(captured) == 3
    for item in captured:
        assert item["params"]["required_score"] == "700"
        assert item["params"]["identifiers"] == QUERY_ID


def test_validate_network_image_bytes_rejects_html_tiny_blank() -> None:
    with pytest.raises(ValueError, match="HTML|html|too small"):
        sd.validate_network_image_bytes(
            b"<html>not an image</html>",
            content_type="text/html",
        )
    with pytest.raises(ValueError, match="too small"):
        sd.validate_network_image_bytes(b"PNG", content_type="image/png")

    blank = Image.new("RGB", (500, 400), (255, 255, 255))
    buf = io.BytesIO()
    blank.save(buf, format="PNG")
    with pytest.raises(ValueError, match="blank|low entropy"):
        sd.validate_network_image_bytes(
            buf.getvalue(),
            content_type="image/png",
            final_url="https://version-12-0.string-db.org/api/highres_image/network",
        )

    valid = (FIXTURES / "tiny_valid.png").read_bytes()
    info = sd.validate_network_image_bytes(
        valid,
        content_type="image/png",
        final_url="https://version-12-0.string-db.org/api/highres_image/network",
    )
    assert info["width"] == 500
    assert info["height"] == 400


def test_normal_collect_does_not_call_interaction_partners(tmp_path: Path) -> None:
    resolve = ToolResult(
        source_name=sd.SOURCE_NAME,
        endpoint_name="resolve_string_identifier",
        gene_symbol="SREBF2",
        request_url="https://example/ids",
        request_params={},
        success=True,
        status_code=200,
        data={
            "gene_symbol": "SREBF2",
            "string_id": QUERY_ID,
            "preferred_name": "SREBF2",
            "taxon_name": "Homo sapiens",
            "ncbi_taxon_id": 9606,
            "identifier_status": "resolved",
            "raw": _load_json("get_string_ids_srebf2.json"),
        },
    )
    network_rows = _load_json("network_srebf2.json")
    network = ToolResult(
        source_name=sd.SOURCE_NAME,
        endpoint_name="network",
        gene_symbol="SREBF2",
        request_url="https://example/network",
        request_params={"required_score": "400", "add_nodes": "30"},
        success=True,
        status_code=200,
        data=network_rows,
    )
    link = ToolResult(
        source_name=sd.SOURCE_NAME,
        endpoint_name="get_link",
        gene_symbol="SREBF2",
        request_url="https://example/get_link",
        request_params={},
        success=True,
        status_code=200,
        data=_load_json("get_link_srebf2.json"),
    )
    state = {
        "gene_symbol": "SREBF2",
        "dossier_run_id": "run-string-no-partners",
        "run_type": "section_bundle",
        "selected_section_keys": ["5a"],
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "tool_results": [],
        "coverage": [],
        "errors": [],
    }
    with (
        patch.object(sd, "resolve_string_identifier", return_value=resolve),
        patch.object(sd, "fetch_network", return_value=network),
        patch.object(sd, "fetch_network_link", return_value=link),
        patch.object(sd, "interaction_partners") as partners,
        patch(
            "gene_dossier.section_5a._persist_string_raw",
            return_value=(None, {"id": "raw-1", "sha256": "abc"}, "abc"),
        ),
    ):
        out = node_generate_section_5a_derived_artifacts(
            state,
            settings=Settings(output_path=tmp_path),
            persist_db=False,
            config=Section5aConfig(
                output_root=tmp_path,
                attempt_network_figure=False,
            ),
        )
    partners.assert_not_called()
    endpoints = [tr.endpoint_name for tr in out["tool_results"]]
    assert "interaction_partners" not in endpoints
    assert out["section_5a_status"]["rendering_status"]["scientific_status"] == "success"
    assert (
        out["section_5a_status"]["rendering_status"]["visual_status"]
        == "not_attempted_optional"
    )


def test_extract_link_url_from_fixture() -> None:
    url = sd.extract_link_url(_load_json("get_link_srebf2.json"))
    assert url is not None
    assert url.startswith("https://")
    assert "string-db.org" in url
