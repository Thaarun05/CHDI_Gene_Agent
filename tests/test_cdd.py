"""Tests for NCBI CDD Batch CD-Search status / cdsid parsing."""

from __future__ import annotations

from gene_dossier.models import ToolResult
from gene_dossier.tools import cdd

SAMPLE_RUNNING = """#Batch CD-search tool\tNIH/NLM/NCBI
#cdsid\tQM3-qcdsearch-18DA3DC3DD8814BD-15234ADDECADD29E
#datatype\thitsConcise Results
#status\t3\tmsg\tJob is still running
"""

SAMPLE_COMPLETED = """#Batch CD-search tool\tNIH/NLM/NCBI
#cdsid\tQM3-qcdsearch-259F25135530B294-255D2F6654B607F0
#datatype\thitsConcise Results
#status\t0
#Start time\t2026-07-22T20:28:23\tRun time\t0:00:00:01
"""

SAMPLE_HTML = """<!DOCTYPE HTML>
<html><body>
<script>var ctrlHandle = "QM3-qcdsearch-151F951D1B74B33D";</script>
<tr><td class="subttl">Search-ID:</td>
<td id="id_cdsid" class="subval">QM3-qcdsearch-151F951D1B74B33D</td></tr>
<input id="hid_DlRid" type="hidden" name="cdsid" value="QM3-qcdsearch-151F951D1B74B33D">
</body></html>
"""

SAMPLE_COMPLETED_WITH_WARNING = """#Batch CD-search tool\tNIH/NLM/NCBI
#cdsid\tQM3-qcdsearch-FC79C8B84F50055-18BB97B1104F6E23
#datatype\thitsConcise Results
#status\t0
#Start time\t2026-07-22T20:54:58\tRun time\t0:00:00:00
#status\tWarning: Too many queries.\tmsg\tError(s) occurred during search, the result may be incomplete
"""

SAMPLE_WARNING_ONLY = """#Batch CD-search tool\tNIH/NLM/NCBI
#cdsid\tQM3-qcdsearch-AAAAAAAAAAAAAAA-BBBBBBBBBBBBBBBB
#status\tWarning: Too many queries.\tmsg\tError(s) occurred during search, the result may be incomplete
"""


def test_parse_status_text_running_multifield_status():
    parsed = cdd.parse_status_text(SAMPLE_RUNNING)
    assert parsed["cdsid"] == "QM3-qcdsearch-18DA3DC3DD8814BD-15234ADDECADD29E"
    assert parsed["status"] == 3


def test_parse_status_text_completed_status_zero():
    parsed = cdd.parse_status_text(SAMPLE_COMPLETED)
    assert parsed["cdsid"] == "QM3-qcdsearch-259F25135530B294-255D2F6654B607F0"
    assert parsed["status"] == 0


def test_parse_status_text_keeps_numeric_status_when_warning_follows():
    parsed = cdd.parse_status_text(SAMPLE_COMPLETED_WITH_WARNING)
    assert parsed["cdsid"] == "QM3-qcdsearch-FC79C8B84F50055-18BB97B1104F6E23"
    assert parsed["status"] == 0
    assert "cdsid" in parsed and "status" in parsed and "fields" in parsed
    assert any("Too many queries" in line for line in parsed["status_lines"])
    assert "Too many queries" in parsed["fields"]["status"]


def test_parse_status_text_warning_only_leaves_status_none():
    parsed = cdd.parse_status_text(SAMPLE_WARNING_ONLY)
    assert parsed["status"] is None
    assert parsed["cdsid"] == "QM3-qcdsearch-AAAAAAAAAAAAAAA-BBBBBBBBBBBBBBBB"
    assert any("Too many queries" in line for line in parsed["status_lines"])
    assert "Too many queries" in parsed["fields"]["status"]


def test_parse_status_text_html_fallbacks():
    parsed = cdd.parse_status_text(SAMPLE_HTML)
    assert parsed["cdsid"] == "QM3-qcdsearch-151F951D1B74B33D"
    assert parsed["status"] is None


def test_parse_status_text_empty_body():
    parsed = cdd.parse_status_text("")
    assert parsed["cdsid"] is None
    assert parsed["status"] is None


def test_parse_status_text_does_not_invent_cdsid():
    parsed = cdd.parse_status_text("<html><body>no search id here</body></html>")
    assert parsed["cdsid"] is None


def test_submit_search_parse_failure_preserves_raw_and_preview(monkeypatch):
    raw = "<!DOCTYPE html>\n<html><body>Welcome to NCBI Batch CD-search</body></html>\n"

    def fake_get_text(*, endpoint_name, gene_symbol, params, settings):
        return ToolResult(
            source_name="CDD",
            endpoint_name=endpoint_name,
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/bwrpsb",
            request_params=params,
            status_code=200,
            data={"raw_text": raw, "content_type": "text/html"},
        )

    monkeypatch.setattr(cdd, "_get_text", fake_get_text)
    result = cdd.submit_search("Q12772", gene_symbol="SREBF2")
    assert result.success is False
    assert result.error_type == "parse_error"
    assert isinstance(result.data, dict)
    assert result.data["raw_text"] == raw
    assert "preview:" in (result.error_message or "")
    assert "Welcome to NCBI Batch CD-search" in (result.error_message or "")
    assert "tdata" in result.request_params
    assert result.request_params["tdata"] == "hits"


def test_submit_search_parses_text_response(monkeypatch):
    def fake_get_text(*, endpoint_name, gene_symbol, params, settings):
        assert params.get("tdata") == "hits"
        return ToolResult(
            source_name="CDD",
            endpoint_name=endpoint_name,
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/bwrpsb",
            request_params=params,
            status_code=200,
            data={"raw_text": SAMPLE_RUNNING, "content_type": "text/plain"},
        )

    monkeypatch.setattr(cdd, "_get_text", fake_get_text)
    result = cdd.submit_search("NP_004590.2", gene_symbol="SREBF2")
    assert result.success is True
    assert result.data["cdsid"] == "QM3-qcdsearch-18DA3DC3DD8814BD-15234ADDECADD29E"
    assert result.data["status"] == 3


def test_poll_status_includes_tdata_hits(monkeypatch):
    seen: dict = {}

    def fake_get_text(*, endpoint_name, gene_symbol, params, settings):
        seen.update(params)
        return ToolResult(
            source_name="CDD",
            endpoint_name=endpoint_name,
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/bwrpsb",
            request_params=params,
            status_code=200,
            data={"raw_text": SAMPLE_COMPLETED, "content_type": "text/plain"},
        )

    monkeypatch.setattr(cdd, "_get_text", fake_get_text)
    result = cdd.poll_status(
        "QM3-qcdsearch-259F25135530B294-255D2F6654B607F0", gene_symbol="SREBF2"
    )
    assert seen.get("tdata") == "hits"
    assert result.success is True
    assert result.data["status"] == 0
