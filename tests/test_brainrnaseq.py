"""Tests for BrainRNASeq CSV download / gene filter client."""

from __future__ import annotations

from gene_dossier.tools import brainrnaseq


def _mock_client(monkeypatch, *, status_code: int, text: str, headers: dict | None = None):
    captured: dict = {}

    class _Headers(dict):
        def get(self, key, default=None):
            for k, v in self.items():
                if str(k).lower() == str(key).lower():
                    return v
            return default

    class _Resp:
        def __init__(self):
            self.status_code = status_code
            self.text = text
            self.headers = _Headers(headers or {})

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
            captured["url"] = url
            captured["headers"] = headers
            return _Resp()

    monkeypatch.setattr(brainrnaseq.httpx, "Client", _Client)
    return captured


SAMPLE_CSV = (
    "gene_id,id,neurons_1,astrocytes_mature_1\n"
    "SREBF2,6721,12.5,3.1\n"
    "OTHER,1,0.1,0.2\n"
)


def test_fetch_gene_expression_http_200_parses_matched_rows(monkeypatch):
    captured = _mock_client(
        monkeypatch,
        status_code=200,
        text=SAMPLE_CSV,
        headers={"content-type": "text/csv"},
    )
    result = brainrnaseq.fetch_gene_expression("SREBF2", species="human")
    assert result.success is True
    assert result.status_code == 200
    assert result.data["match_count"] == 1
    assert result.data["matched_rows"][0]["gene_id"] == "SREBF2"
    assert result.data["expression_summaries"][0]["gene_id"] == "SREBF2"
    assert "neurons_1" in result.data["expression_summaries"][0]["cell_type_values"]
    assert captured["headers"] == brainrnaseq.REQUEST_HEADERS


def test_download_csv_http_403_access_forbidden(monkeypatch):
    _mock_client(
        monkeypatch,
        status_code=403,
        text="<html>Forbidden by WAF</html>",
        headers={
            "content-type": "text/html",
            "server": "cloudflare",
            "cf-ray": "abc123",
        },
    )
    result = brainrnaseq.download_csv("human", gene_symbol="SREBF2")
    assert result.success is False
    assert result.status_code == 403
    assert result.error_type == "access_forbidden"
    assert "403" in (result.error_message or "")
    assert "forbidden" in (result.error_message or "").lower()
    assert result.request_url == brainrnaseq.HUMAN_CSV_URL


def test_download_csv_http_403_preserves_preview_and_details(monkeypatch):
    body = "<html>" + ("denied " * 200) + "</html>"
    _mock_client(
        monkeypatch,
        status_code=403,
        text=body,
        headers={
            "content-type": "text/html; charset=utf-8",
            "server": "cloudflare",
            "cf-ray": "ray-xyz",
            "set-cookie": "should-not-leak=1",
        },
    )
    result = brainrnaseq.download_csv("human", gene_symbol="SREBF2")
    assert result.success is False
    assert result.data["content_type"] == "text/html; charset=utf-8"
    assert result.data["response_headers"]["content-type"] == "text/html; charset=utf-8"
    assert result.data["response_headers"]["server"] == "cloudflare"
    assert result.data["response_headers"]["cf-ray"] == "ray-xyz"
    assert "set-cookie" not in result.data["response_headers"]
    assert "raw_csv" not in result.data
    preview = result.data["raw_text_preview"]
    assert preview
    assert len(preview) <= brainrnaseq.RAW_TEXT_PREVIEW_LIMIT
    assert result.request_url == brainrnaseq.HUMAN_CSV_URL

    forwarded = brainrnaseq.fetch_gene_expression("SREBF2", species="human")
    assert forwarded.success is False
    assert forwarded.error_type == "access_forbidden"
    assert forwarded.status_code == 403
    assert forwarded.data["raw_text_preview"] == preview


def test_download_csv_sends_polite_headers(monkeypatch):
    captured = _mock_client(
        monkeypatch,
        status_code=200,
        text=SAMPLE_CSV,
        headers={"content-type": "text/csv"},
    )
    result = brainrnaseq.download_csv("human", gene_symbol="SREBF2")
    assert result.success is True
    assert captured["headers"]["User-Agent"] == brainrnaseq.REQUEST_HEADERS["User-Agent"]
    assert captured["headers"]["Accept"] == brainrnaseq.REQUEST_HEADERS["Accept"]
    assert captured["headers"]["Referer"] == brainrnaseq.REQUEST_HEADERS["Referer"]
    assert result.request_params["request_headers"] == brainrnaseq.REQUEST_HEADERS
    assert "github.com" not in result.request_params["request_headers"]["User-Agent"].lower()


def test_exact_gene_matching_still_works():
    rows = [
        {"gene_id": "SREBF2", "id": "6721"},
        {"gene_id": "SREBF1", "id": "6720"},
        {"gene_id": "XSREBF2Y", "id": "9"},
    ]
    matched = brainrnaseq.filter_gene_rows(rows, "SREBF2")
    assert [r["gene_id"] for r in matched] == ["SREBF2"]
    assert brainrnaseq.row_matches_gene(rows[2], "SREBF2") is False
    assert (
        brainrnaseq.row_matches_gene(
            rows[2], "SREBF2", allow_substring_match=True
        )
        is True
    )
