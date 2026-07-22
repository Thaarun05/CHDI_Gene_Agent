"""Tests for Allen Brain Atlas HBA probe / expression client."""

from __future__ import annotations

from gene_dossier.models import ToolResult
from gene_dossier.tools import allen_brain


def test_build_expression_criteria_quotes_probe_id():
    criteria = allen_brain.build_expression_criteria(1051154)
    assert criteria == "service::human_microarray_expression[probes$eq'1051154']"
    assert "probes$eq'1051154'" in criteria
    assert "probes$eq1051154]" not in criteria


def test_build_probe_lookup_criteria_unchanged():
    criteria = allen_brain.build_probe_lookup_criteria("SREBF2")
    assert criteria.startswith("model::Probe,rma::criteria,")
    assert "[probe_type$eq'DNA']" in criteria
    assert "products[abbreviation$eq'HumanMA']" in criteria
    assert "gene[acronym$eq'SREBF2']" in criteria
    assert "rma::options[only$eq'probes.id,probes.name,genes.acronym,genes.name,genes.entrez_id']" in criteria


def test_request_query_soft_fails_when_allen_success_false(monkeypatch):
    payload = {
        "success": False,
        "msg": (
            "Error in query: service::human_microarray_expression[probes$eq1051154], "
            "Informatics service request failed."
        ),
    }

    class _Resp:
        status_code = 200
        headers = {"content-type": "application/json"}

        def json(self):
            return payload

        @property
        def text(self):
            return str(payload)

        @property
        def is_success(self):
            return True

    class _Client:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, url, params=None):
            return _Resp()

    monkeypatch.setattr(allen_brain.httpx, "Client", _Client)
    result = allen_brain.microarray_expression(1051154, gene_symbol="SREBF2")
    assert result.success is False
    assert result.error_type == "api_error"
    assert isinstance(result.data, dict)
    assert result.data["success"] is False
    assert result.data["msg"] == payload["msg"]
    assert "Informatics service request failed" in (result.error_message or "")
    assert "probes$eq1051154" in (result.error_message or "")


def test_request_query_preview_when_msg_missing(monkeypatch):
    payload = {"success": False, "detail": "service unavailable"}

    class _Resp:
        status_code = 200
        headers = {"content-type": "application/json"}

        def json(self):
            return payload

        @property
        def text(self):
            return str(payload)

        @property
        def is_success(self):
            return True

    class _Client:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, url, params=None):
            return _Resp()

    monkeypatch.setattr(allen_brain.httpx, "Client", _Client)
    result = allen_brain.microarray_expression(1051154, gene_symbol="SREBF2")
    assert result.success is False
    assert result.data == payload
    assert "Allen API returned success=false" in (result.error_message or "")
    assert "preview:" in (result.error_message or "")
    assert "service unavailable" in (result.error_message or "")


def test_fetch_hba_expression_partial_success(monkeypatch):
    def fake_microarray(probe_id, *, gene_symbol="", settings=None):
        pid = int(probe_id)
        if pid == 1051154:
            return ToolResult(
                source_name="Allen Brain Atlas",
                endpoint_name="microarray_expression",
                success=False,
                gene_symbol=gene_symbol,
                request_url="https://example.test/fail",
                request_params={"probe_id": str(pid)},
                status_code=200,
                data={"success": False, "msg": f"failed probe {pid}"},
                error_type="api_error",
                error_message=f"failed probe {pid}",
            )
        return ToolResult(
            source_name="Allen Brain Atlas",
            endpoint_name="microarray_expression",
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/ok",
            request_params={"probe_id": str(pid)},
            status_code=200,
            data={
                "success": True,
                "msg": {"probes": [{"id": pid}], "samples": [], "expression": []},
            },
        )

    monkeypatch.setattr(allen_brain, "microarray_expression", fake_microarray)
    result = allen_brain.fetch_hba_expression(
        "SREBF2",
        probe_ids=[1051154, 1067243],
        max_probes=2,
        max_expression_attempts=1,
        retry_sleep_seconds=0,
    )
    assert result.success is True
    assert "1067243" in result.data["expressions"]
    assert "1051154" not in result.data["expressions"]
    assert len(result.data["failed_probes"]) == 1
    assert result.data["failed_probes"][0]["probe_id"] == 1051154
    assert result.data["failed_probes"][0]["data"]["success"] is False


def test_fetch_hba_expression_all_probes_fail(monkeypatch):
    def fake_microarray(probe_id, *, gene_symbol="", settings=None):
        pid = int(probe_id)
        return ToolResult(
            source_name="Allen Brain Atlas",
            endpoint_name="microarray_expression",
            success=False,
            gene_symbol=gene_symbol,
            request_url=f"https://example.test/{pid}",
            request_params={"probe_id": str(pid)},
            status_code=200,
            data={"success": False, "msg": f"failed probe {pid}"},
            error_type="api_error",
            error_message=f"failed probe {pid}",
        )

    monkeypatch.setattr(allen_brain, "microarray_expression", fake_microarray)
    result = allen_brain.fetch_hba_expression(
        "SREBF2",
        probe_ids=[1051154, 1067243],
        max_probes=2,
        max_expression_attempts=1,
        retry_sleep_seconds=0,
    )
    assert result.success is False
    assert result.error_type == "api_error"
    assert "failed probe 1051154" in (result.error_message or "")
    assert result.data["expressions"] == {}
    assert len(result.data["failed_probes"]) == 2
    assert all("data" in row for row in result.data["failed_probes"])


def test_fetch_hba_expression_retries_transient_then_succeeds(monkeypatch):
    calls: list[int] = []

    def fake_microarray(probe_id, *, gene_symbol="", settings=None):
        pid = int(probe_id)
        calls.append(pid)
        if len(calls) == 1:
            return ToolResult(
                source_name="Allen Brain Atlas",
                endpoint_name="microarray_expression",
                success=False,
                gene_symbol=gene_symbol,
                request_url="https://example.test/fail",
                request_params={"probe_id": str(pid)},
                status_code=200,
                data={
                    "success": False,
                    "msg": (
                        "Error in query: "
                        "service::human_microarray_expression[probes$eq'1051154'], "
                        "Informatics service request failed."
                    ),
                },
                error_type="api_error",
                error_message=(
                    "Error in query: "
                    "service::human_microarray_expression[probes$eq'1051154'], "
                    "Informatics service request failed."
                ),
            )
        return ToolResult(
            source_name="Allen Brain Atlas",
            endpoint_name="microarray_expression",
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/ok",
            request_params={"probe_id": str(pid)},
            status_code=200,
            data={
                "success": True,
                "msg": {"probes": [{"id": pid}], "samples": [], "expression": []},
            },
        )

    sleeps: list[float] = []
    monkeypatch.setattr(allen_brain, "microarray_expression", fake_microarray)
    monkeypatch.setattr(allen_brain.time, "sleep", lambda s: sleeps.append(s))
    result = allen_brain.fetch_hba_expression(
        "SREBF2",
        probe_ids=[1051154],
        max_probes=1,
        max_expression_attempts=3,
        retry_sleep_seconds=1.0,
    )
    assert result.success is True
    assert "1051154" in result.data["expressions"]
    assert result.data["expressions"]["1051154"]["success"] is True
    assert result.data["failed_probes"] == []
    attempts = result.data["expression_attempts"]["1051154"]
    assert len(attempts) == 2
    assert attempts[0]["success"] is False
    assert "Informatics service request failed" in (attempts[0]["error_message"] or "")
    assert attempts[1]["success"] is True
    assert calls == [1051154, 1051154]
    assert sleeps == [1.0]


def test_fetch_hba_expression_all_attempts_fail_preserves_attempts(monkeypatch):
    calls = 0

    def fake_microarray(probe_id, *, gene_symbol="", settings=None):
        nonlocal calls
        calls += 1
        pid = int(probe_id)
        return ToolResult(
            source_name="Allen Brain Atlas",
            endpoint_name="microarray_expression",
            success=False,
            gene_symbol=gene_symbol,
            request_url=f"https://example.test/{pid}/{calls}",
            request_params={"probe_id": str(pid)},
            status_code=200,
            data={
                "success": False,
                "msg": (
                    f"Error in query: probes$eq'{pid}', "
                    "Informatics service request failed."
                ),
            },
            error_type="api_error",
            error_message=(
                f"Error in query: probes$eq'{pid}', "
                "Informatics service request failed."
            ),
        )

    monkeypatch.setattr(allen_brain, "microarray_expression", fake_microarray)
    monkeypatch.setattr(allen_brain.time, "sleep", lambda s: None)
    result = allen_brain.fetch_hba_expression(
        "SREBF2",
        probe_ids=[1051154],
        max_probes=1,
        max_expression_attempts=3,
        retry_sleep_seconds=0.5,
    )
    assert result.success is False
    assert calls == 3
    assert result.data["expressions"] == {}
    assert len(result.data["failed_probes"]) == 1
    attempts = result.data["failed_probes"][0]["expression_attempts"]
    assert len(attempts) == 3
    assert [a["attempt"] for a in attempts] == [1, 2, 3]
    assert all(a["success"] is False for a in attempts)
    assert result.data["expression_attempts"]["1051154"] == attempts


def test_fetch_hba_expression_does_not_retry_invalid_request(monkeypatch):
    calls = 0

    def fake_microarray(probe_id, *, gene_symbol="", settings=None):
        nonlocal calls
        calls += 1
        return ToolResult(
            source_name="Allen Brain Atlas",
            endpoint_name="microarray_expression",
            success=False,
            gene_symbol=gene_symbol,
            request_url="https://example.test/invalid",
            request_params={},
            status_code=None,
            data=None,
            error_type="invalid_request",
            error_message="probe_id is required",
        )

    monkeypatch.setattr(allen_brain, "microarray_expression", fake_microarray)
    sleeps: list[float] = []
    monkeypatch.setattr(allen_brain.time, "sleep", lambda s: sleeps.append(s))
    result = allen_brain.fetch_hba_expression(
        "SREBF2",
        probe_ids=[1051154],
        max_probes=1,
        max_expression_attempts=3,
        retry_sleep_seconds=1.0,
    )
    assert result.success is False
    assert calls == 1
    assert sleeps == []
    assert result.data["failed_probes"][0]["error_type"] == "invalid_request"
    assert len(result.data["failed_probes"][0]["expression_attempts"]) == 1
