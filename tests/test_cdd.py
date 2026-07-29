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

SAMPLE_INVALID_ID = """#Batch CD-search tool\tNIH/NLM/NCBI
#cdsid\tQM3-qcdsearch-BAD
#status\t1\tmsg\tInvalid search ID
"""

SAMPLE_ABUSE_INPUT = """#Batch CD-search tool\tNIH/NLM/NCBI
#cdsid\tQM3-qcdsearch-ABUSE
#status\t6\tmsg\tABUSEINPUT
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
    assert result.request_params["tdata"] == "hits"
    assert result.request_params["smode"] == "live"
    assert result.request_params["filter"] == "false"
    assert result.request_params["compbasedadj"] == "1"


def test_submit_search_parses_text_response(monkeypatch):
    def fake_get_text(*, endpoint_name, gene_symbol, params, settings):
        assert params.get("tdata") == "hits"
        assert params.get("smode") == "live"
        assert params.get("filter") == "false"
        assert params.get("compbasedadj") == "1"
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
    assert result.data["cdsid"] == "QM3-qcdsearch-18DA3DC3DD8814BD"
    assert result.data["master_cdsid"] == "QM3-qcdsearch-18DA3DC3DD8814BD"
    assert result.data["submit_request_cdsid"] == "QM3-qcdsearch-18DA3DC3DD8814BD-15234ADDECADD29E"
    assert result.data["status"] == 3


def test_submit_search_terminal_status_fails(monkeypatch):
    def fake_get_text(*, endpoint_name, gene_symbol, params, settings):
        return ToolResult(
            source_name="CDD",
            endpoint_name=endpoint_name,
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/bwrpsb",
            request_params=params,
            status_code=200,
            data={"raw_text": SAMPLE_INVALID_ID, "content_type": "text/plain"},
        )

    monkeypatch.setattr(cdd, "_get_text", fake_get_text)
    result = cdd.submit_search("Q12772", gene_symbol="SREBF2")
    assert result.success is False
    assert result.error_type == "terminal_status"
    assert "status 1" in (result.error_message or "")


def test_submit_search_status_6_terminal_status_fails(monkeypatch):
    def fake_get_text(*, endpoint_name, gene_symbol, params, settings):
        return ToolResult(
            source_name="CDD",
            endpoint_name=endpoint_name,
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/bwrpsb",
            request_params=params,
            status_code=200,
            data={"raw_text": SAMPLE_ABUSE_INPUT, "content_type": "text/plain"},
        )

    monkeypatch.setattr(cdd, "_get_text", fake_get_text)
    result = cdd.submit_search("Q12772", gene_symbol="SREBF2")
    assert result.success is False
    assert result.error_type == "terminal_status"
    assert "status 6" in (result.error_message or "")


def test_fetch_domains_terminal_poll_status_stops(monkeypatch):
    monkeypatch.setattr(
        cdd,
        "submit_search",
        lambda *args, **kwargs: ToolResult(
            source_name="CDD",
            endpoint_name="submit_search",
            success=True,
            gene_symbol="SREBF2",
            request_url="https://example.test/submit",
            request_params={},
            status_code=200,
            data={"cdsid": "QM3-qcdsearch-OK", "status": 3},
        ),
    )
    monkeypatch.setattr(
        cdd,
        "poll_status",
        lambda *args, **kwargs: ToolResult(
            source_name="CDD",
            endpoint_name="poll_status",
            success=True,
            gene_symbol="SREBF2",
            request_url="https://example.test/poll",
            request_params={},
            status_code=200,
            data={"cdsid": "QM3-qcdsearch-OK", "status": 5},
        ),
    )
    result = cdd.fetch_domains(
        "Q12772",
        gene_symbol="SREBF2",
        poll_interval_seconds=0,
        max_polls=3,
    )
    assert result.success is False
    assert result.error_type == "terminal_status"
    assert result.data["status"] == 5


def test_parse_features_text_rejects_hit_shaped_rows():
    rows = cdd.parse_features_text(
        """#Batch CD-search tool\tNIH/NLM/NCBI
Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort name\tIncomplete\tSuperfamily
Q#1 - NP_006718.2\tspecific\t206637\t164\t265\t1e-29\t112.7\tcd11304\tCadherin_repeat\t-\tcl46864
"""
    )
    assert rows == []


def test_parse_features_text_accepts_feature_columns():
    rows = cdd.parse_features_text(
        """Query\tFeature name\tFeature type\tQuery residues\tPSSM-ID\tAccession
Q#1 - NP_006718.2\tCa2+ binding site\tion binding site\t101,103\t206637\tcd11304
"""
    )
    assert rows[0]["Feature name"] == "Ca2+ binding site"


def test_parse_features_text_accepts_live_title_coordinates_columns():
    rows = cdd.parse_features_text(
        """Query\tType\tTitle\tcoordinates\tcomplete size\tmapped size\tsource domain
Q#1 - NP_006718.2[cadherin-10 isoform 1 preproprotein [Homo sapiens]]\tspecific\tCa2+ binding site\tE171,M172,N226,E228,D261,N263,D264\t7\t7\t206637
"""
    )
    assert rows[0]["Title"] == "Ca2+ binding site"
    summary = cdd.summarize_feature(rows[0])
    assert summary["feature_name"] == "Ca2+ binding site"
    assert summary["feature_type"] == "specific"
    assert summary["query_residues"] == "E171,M172,N226,E228,D261,N263,D264"
    assert summary["domain_accession"] == "206637"


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
