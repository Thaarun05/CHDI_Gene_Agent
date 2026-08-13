"""Tests for Reactome ContentService client retry / soft-fail behavior."""

from __future__ import annotations

from typing import Any

from gene_dossier.models import ToolResult
from gene_dossier.tools import reactome


def _mock_client(monkeypatch, responses: list):
    """Patch httpx.Client so successive GET calls return queued responses/exceptions."""
    queue = list(responses)

    class _Resp:
        def __init__(self, status_code: int, payload):
            self.status_code = status_code
            self._payload = payload
            self.text = str(payload)

        def json(self):
            if isinstance(self._payload, Exception):
                raise self._payload
            return self._payload

        @property
        def is_success(self):
            return 200 <= self.status_code < 300

    class _Client:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, url, headers=None):
            item = queue.pop(0)
            if isinstance(item, Exception):
                raise item
            status, payload = item
            return _Resp(status, payload)

    monkeypatch.setattr(reactome.httpx, "Client", _Client)
    monkeypatch.setattr(reactome.time, "sleep", lambda s: None)
    return queue


def test_fetch_pathways_retries_http_521_then_succeeds(monkeypatch):
    pathway = {
        "dbId": 1,
        "stId": "R-HSA-1655829",
        "displayName": "Cholesterol biosynthesis",
        "speciesName": "Homo sapiens",
    }
    _mock_client(
        monkeypatch,
        [
            (521, {"raw_text": "Web server is down"}),
            (200, [pathway]),
        ],
    )
    result = reactome.fetch_pathways(
        "Q12772",
        gene_symbol="SREBF2",
        max_attempts=3,
        retry_sleep_seconds=1.0,
    )
    assert result.success is True
    assert result.status_code == 200
    assert result.data["uniprot_accession"] == "Q12772"
    assert result.data["pathway_count"] == 1
    assert result.data["pathways"][0]["stId"] == "R-HSA-1655829"
    assert result.data["pathway_summaries"][0]["st_id"] == "R-HSA-1655829"
    attempts = result.data["request_attempts"]
    assert len(attempts) == 2
    assert attempts[0]["success"] is False
    assert attempts[0]["status_code"] == 521
    assert attempts[0]["error_message"] == "HTTP 521"
    assert attempts[1]["success"] is True


def test_fetch_pathways_all_attempts_fail_preserves_attempts(monkeypatch):
    _mock_client(
        monkeypatch,
        [
            (521, {"raw_text": "down-1"}),
            (521, {"raw_text": "down-2"}),
            (521, {"raw_text": "down-3"}),
        ],
    )
    result = reactome.fetch_pathways(
        "Q12772",
        gene_symbol="SREBF2",
        max_attempts=3,
        retry_sleep_seconds=0.5,
    )
    assert result.success is False
    assert result.error_type == "http_error"
    assert result.error_message == "HTTP 521"
    assert result.status_code == 521
    attempts = result.data["request_attempts"]
    assert len(attempts) == 3
    assert [a["attempt"] for a in attempts] == [1, 2, 3]
    assert all(a["status_code"] == 521 for a in attempts)
    assert all(a["success"] is False for a in attempts)


def test_fetch_pathways_does_not_retry_http_400(monkeypatch):
    queue = _mock_client(
        monkeypatch,
        [
            (400, {"messages": ["Bad request"]}),
            (200, [{"stId": "should-not-be-called"}]),
        ],
    )
    result = reactome.fetch_pathways(
        "Q12772",
        gene_symbol="SREBF2",
        max_attempts=3,
        retry_sleep_seconds=1.0,
    )
    assert result.success is False
    assert result.status_code == 400
    assert result.error_message == "HTTP 400"
    assert "request_attempts" not in (result.data or {})
    # Second queued success response must remain unused.
    assert len(queue) == 1


def test_fetch_pathways_successful_parse_shape(monkeypatch):
    pathways = [
        {
            "dbId": 42,
            "stId": "R-HSA-1655829",
            "stIdVersion": "R-HSA-1655829.1",
            "displayName": "Cholesterol biosynthesis",
            "speciesName": "Homo sapiens",
            "doi": "10.3180/REACT_1111.1",
            "hasDiagram": True,
            "schemaClass": "Pathway",
        }
    ]
    _mock_client(monkeypatch, [(200, pathways)])
    result = reactome.fetch_pathways("Q12772", gene_symbol="SREBF2", max_attempts=3)
    assert result.success is True
    assert result.data["pathway_count"] == 1
    assert result.data["pathways"] == pathways
    summary = result.data["pathway_summaries"][0]
    assert summary["st_id"] == "R-HSA-1655829"
    assert summary["display_name"] == "Cholesterol biosynthesis"
    assert summary["detail_url"].endswith("/R-HSA-1655829")
    assert "FLG=Q12772" in summary["browser_url"]
    assert "request_attempts" not in result.data


def test_pathways_by_uniprot_missing_accession_not_retried(monkeypatch):
    calls = {"n": 0}

    class _Client:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, url, headers=None):
            calls["n"] += 1
            raise AssertionError("HTTP should not be called for missing accession")

    monkeypatch.setattr(reactome.httpx, "Client", _Client)
    result = reactome.pathways_by_uniprot("", gene_symbol="SREBF2", max_attempts=3)
    assert result.success is False
    assert result.error_type == "invalid_request"
    assert calls["n"] == 0


def test_fetch_pathways_uses_validated_request_contract(monkeypatch):
    captured: dict[str, Any] = {}

    def fake_request_json(**kwargs: Any) -> ToolResult:
        captured.update(kwargs)
        return ToolResult(
            source_name="Reactome",
            endpoint_name=kwargs["endpoint_name"],
            success=True,
            gene_symbol=kwargs["gene_symbol"],
            request_url=(
                "https://reactome.org/ContentService/"
                "data/mapping/UniProt/Q9Y2M0/pathways"
            ),
            request_params=kwargs["request_params"],
            status_code=200,
            data=[
                {
                    "stId": "R-HSA-73894",
                    "displayName": "DNA Repair",
                    "speciesName": "Homo sapiens",
                }
            ],
        )

    monkeypatch.setattr(reactome, "_request_json", fake_request_json)
    result = reactome.fetch_pathways(
        "Q9Y2M0",
        gene_symbol="FAN1",
        max_attempts=1,
        retry_sleep_seconds=0,
    )

    assert captured["endpoint_name"] == "pathways_by_uniprot"
    assert captured["gene_symbol"] == "FAN1"
    assert captured["path"] == "data/mapping/UniProt/Q9Y2M0/pathways"
    assert captured["request_params"] == {
        "uniprot_accession": "Q9Y2M0",
        "max_attempts": 1,
        "retry_sleep_seconds": 0.0,
    }
    assert result.success
    assert result.source_name == "Reactome"
    assert result.endpoint_name == "fetch_pathways"
    assert result.request_url.endswith("/data/mapping/UniProt/Q9Y2M0/pathways")
    assert result.request_params == captured["request_params"]
    assert result.data["uniprot_accession"] == "Q9Y2M0"
