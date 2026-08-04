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
    "gene_id,id,neurons_1,astrocytes_mature_1,endothelial_1,microglia_1,"
    "oligodendrocytes_1,astrocytes_fetal_1\n"
    "GENEX,9999,12.5,3.1,0.4,1.2,0.8,2.0\n"
    "OTHER,1,0.1,0.2,0.0,0.0,0.0,0.0\n"
)

SAMPLE_MOUSE_CSV = (
    "gene_id,id,astrocytes_1,neurons_1,opc_1,newly_formed_oligodendrocyte_1,"
    "myelinating_oligodendrocyte_1,endothelial_1,microglia_macrophage_1\n"
    "Genex,8888,5.0,10.0,2.0,1.0,3.0,0.5,4.0\n"
    "Other,2,1.0,1.0,1.0,1.0,1.0,1.0,1.0\n"
)


def test_fetch_gene_expression_http_200_parses_matched_rows(monkeypatch):
    captured = _mock_client(
        monkeypatch,
        status_code=200,
        text=SAMPLE_CSV,
        headers={"content-type": "text/csv"},
    )
    result = brainrnaseq.fetch_gene_expression("GENEX", species="human")
    assert result.success is True
    assert result.status_code == 200
    assert result.data["match_count"] == 1
    assert result.data["matched_rows"][0]["gene_id"] == "GENEX"
    assert result.data["expression_summaries"][0]["gene_id"] == "GENEX"
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
    result = brainrnaseq.download_csv("human", gene_symbol="GENEX")
    assert result.success is False
    assert result.status_code == 403
    assert result.error_type == "access_forbidden"
    assert "403" in (result.error_message or "")
    assert "forbidden" in (result.error_message or "").lower()
    assert result.request_url == brainrnaseq.HUMAN_CSV_URL
    assert brainrnaseq._should_fallback_to_browser(result) is True


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
    result = brainrnaseq.download_csv("human", gene_symbol="GENEX")
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

    forwarded = brainrnaseq.fetch_gene_expression("GENEX", species="human")
    assert forwarded.success is False
    assert forwarded.error_type == "access_forbidden"
    assert forwarded.status_code == 403
    assert forwarded.data["raw_text_preview"] == preview


def test_download_csv_rejects_html_soft_200(monkeypatch):
    _mock_client(
        monkeypatch,
        status_code=200,
        text="<!DOCTYPE html><html><body>challenge</body></html>",
        headers={"content-type": "text/html"},
    )
    result = brainrnaseq.download_csv("human", gene_symbol="GENEX")
    assert result.success is False
    assert result.error_type == "invalid_or_html_response"
    assert brainrnaseq._should_fallback_to_browser(result) is True


def test_validate_brainrnaseq_csv_bytes_human_and_mouse():
    ok_h = brainrnaseq.validate_brainrnaseq_csv_bytes(SAMPLE_CSV, species="human")
    assert ok_h["ok"] is True
    assert ok_h["row_count"] >= 1
    ok_m = brainrnaseq.validate_brainrnaseq_csv_bytes(SAMPLE_MOUSE_CSV, species="mouse")
    assert ok_m["ok"] is True
    bad = brainrnaseq.validate_brainrnaseq_csv_bytes(
        "<html>nope</html>", species="human"
    )
    assert bad["ok"] is False


def test_normalize_published_csv_url_www_host():
    www = (
        "https://www.brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-124.csv"
    )
    assert brainrnaseq.normalize_published_csv_url(www) == brainrnaseq.HUMAN_CSV_URL
    assert brainrnaseq.normalize_published_csv_url("https://evil.example/x.csv") is None


def test_download_csv_sends_polite_headers(monkeypatch):
    captured = _mock_client(
        monkeypatch,
        status_code=200,
        text=SAMPLE_CSV,
        headers={"content-type": "text/csv"},
    )
    result = brainrnaseq.download_csv("human", gene_symbol="GENEX")
    assert result.success is True
    assert captured["headers"]["User-Agent"] == brainrnaseq.REQUEST_HEADERS["User-Agent"]
    assert captured["headers"]["Accept"] == brainrnaseq.REQUEST_HEADERS["Accept"]
    assert captured["headers"]["Referer"] == brainrnaseq.REQUEST_HEADERS["Referer"]
    assert result.request_params["request_headers"] == brainrnaseq.REQUEST_HEADERS
    assert "github.com" not in result.request_params["request_headers"]["User-Agent"].lower()


def test_exact_gene_matching_still_works():
    rows = [
        {"gene_id": "GENEX", "id": "9999"},
        {"gene_id": "GENEY", "id": "9998"},
        {"gene_id": "XGENEXY", "id": "9"},
        {"gene_id": "GENEX - Homo sapiens", "id": "HGNC:1"},
        {"gene_id": "Genex - Mus musculus", "id": "2"},
        {"gene_id": "GENEX - antisense", "id": "3"},
        {"gene_id": "GENEX - pseudogene", "id": "4"},
    ]
    matched = brainrnaseq.filter_gene_rows(rows, "GENEX", species="human")
    assert [r["gene_id"] for r in matched] == ["GENEX", "GENEX - Homo sapiens"]
    assert brainrnaseq.row_matches_gene(rows[2], "GENEX") is False
    assert (
        brainrnaseq.row_matches_gene(
            rows[2], "GENEX", allow_substring_match=True
        )
        is True
    )
    assert brainrnaseq.row_matches_gene(rows[4], "Genex", species="mouse") is True
    assert brainrnaseq.row_matches_gene(rows[4], "GENEX", species="human") is False
    assert brainrnaseq.row_matches_gene(rows[3], "GENEX", species="human") is True
    assert brainrnaseq.row_matches_gene(rows[5], "GENEX") is False
    assert brainrnaseq.row_matches_gene(rows[6], "GENEX") is False
    assert brainrnaseq.normalize_gene_identifier("GENEX - antisense") is None
    assert brainrnaseq.normalize_gene_identifier("GENEX - Homo sapiens") == "GENEX"
    assert brainrnaseq.normalize_gene_identifier("Genex - Mus musculus") == "Genex"


def test_browser_download_chrome_success(monkeypatch):
    class _FakeResp:
        def __init__(self, url, body):
            self.url = url
            self._body = body
            self.headers = {"content-type": "text/csv"}
            self.status = 200

        def body(self):
            return self._body

    class _FakePage:
        def __init__(self):
            self.url = brainrnaseq.HUMAN_CSV_URL
            self._handlers = []

        def on(self, event, handler):
            self._handlers.append((event, handler))

        def expect_download(self, timeout=0):
            raise RuntimeError("no download event")

        def goto(self, url, wait_until=None, timeout=None):
            resp = _FakeResp(url, SAMPLE_CSV.encode("utf-8"))
            for event, handler in self._handlers:
                if event == "response":
                    handler(resp)
            return resp

    class _FakeContext:
        def new_page(self):
            return _FakePage()

        def close(self):
            return None

    class _FakeBrowser:
        def new_context(self, accept_downloads=True):
            return _FakeContext()

        def close(self):
            return None

    class _FakeChromium:
        def launch(self, headless=True, channel=None):
            if channel == "chrome":
                return _FakeBrowser()
            raise RuntimeError("should not fall through")

    class _FakePW:
        chromium = _FakeChromium()

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    monkeypatch.setattr(
        brainrnaseq,
        "_launch_playwright_chromium",
        lambda pw: (_FakeBrowser(), [{"channel": "chrome", "success": True}]),
    )

    import sys
    import types

    fake_sync = types.ModuleType("playwright.sync_api")
    fake_sync.sync_playwright = lambda: _FakePW()
    fake_playwright = types.ModuleType("playwright")
    fake_playwright.sync_api = fake_sync
    monkeypatch.setitem(sys.modules, "playwright", fake_playwright)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", fake_sync)

    result = brainrnaseq.download_csv_via_browser("human", gene_symbol="GENEX")
    assert result.success is True
    assert result.endpoint_name == "download_csv_browser"
    assert result.data["retrieval_method"] == "official_browser_download"
    assert result.data["browser_channel"] == "chrome"
    assert result.data["final_url"] == brainrnaseq.HUMAN_CSV_URL
    assert result.data["capture_via"] in {
        "response_body",
        "navigation_response_body",
    }
    assert result.status_code == 200
    assert result.data["http_status_observed"] is True
    assert "GENEX" in result.data["raw_csv"]


def test_browser_download_event_has_no_fabricated_status(monkeypatch):
    class _FakeDownload:
        def __init__(self, url, payload: bytes):
            self.url = url
            self._payload = payload

        def save_as(self, path: str) -> None:
            from pathlib import Path

            Path(path).write_bytes(self._payload)

    class _DownloadCtx:
        def __init__(self, download):
            self._download = download

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        @property
        def value(self):
            return self._download

    class _FakePage:
        def __init__(self):
            self.url = brainrnaseq.HUMAN_CSV_URL

        def on(self, event, handler):
            return None

        def expect_download(self, timeout=0):
            return _DownloadCtx(
                _FakeDownload(
                    brainrnaseq.HUMAN_CSV_URL, SAMPLE_CSV.encode("utf-8")
                )
            )

        def goto(self, url, wait_until=None, timeout=None):
            raise RuntimeError("Download is starting")

    class _FakeContext:
        def new_page(self):
            return _FakePage()

        def close(self):
            return None

    class _FakeBrowser:
        def new_context(self, accept_downloads=True):
            return _FakeContext()

        def close(self):
            return None

    monkeypatch.setattr(
        brainrnaseq,
        "_launch_playwright_chromium",
        lambda pw: (_FakeBrowser(), [{"channel": "chrome", "success": True}]),
    )

    import sys
    import types

    fake_sync = types.ModuleType("playwright.sync_api")
    fake_sync.sync_playwright = lambda: type(
        "PW",
        (),
        {
            "__enter__": lambda self: self,
            "__exit__": lambda self, *a: False,
            "chromium": object(),
        },
    )()
    # Minimal Error type for "Download is starting" handling.
    fake_sync.Error = RuntimeError
    fake_playwright = types.ModuleType("playwright")
    fake_playwright.sync_api = fake_sync
    monkeypatch.setitem(sys.modules, "playwright", fake_playwright)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", fake_sync)

    result = brainrnaseq.download_csv_via_browser("human", gene_symbol="GENEX")
    assert result.success is True
    assert result.data["capture_via"] == "download_event"
    assert result.status_code is None
    assert result.data["http_status_observed"] is False
    assert "GENEX" in result.data["raw_csv"]

def test_browser_download_chrome_fail_chromium_success(monkeypatch):
    class _FakeResp:
        def __init__(self, url, body):
            self.url = url
            self._body = body
            self.headers = {"content-type": "text/csv"}
            self.status = 200

        def body(self):
            return self._body

    class _FakePage:
        def __init__(self):
            self.url = brainrnaseq.MOUSE_CSV_URL
            self._handlers = []

        def on(self, event, handler):
            self._handlers.append((event, handler))

        def expect_download(self, timeout=0):
            raise RuntimeError("no download event")

        def goto(self, url, wait_until=None, timeout=None):
            resp = _FakeResp(url, SAMPLE_MOUSE_CSV.encode("utf-8"))
            for event, handler in self._handlers:
                if event == "response":
                    handler(resp)
            return resp

    class _FakeContext:
        def new_page(self):
            return _FakePage()

        def close(self):
            return None

    class _FakeBrowser:
        def new_context(self, accept_downloads=True):
            return _FakeContext()

        def close(self):
            return None

    monkeypatch.setattr(
        brainrnaseq,
        "_launch_playwright_chromium",
        lambda pw: (
            _FakeBrowser(),
            [
                {"channel": "chrome", "success": False, "error_type": "Error"},
                {"channel": "chromium", "success": True},
            ],
        ),
    )

    import sys
    import types

    class _FakePW:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    fake_sync = types.ModuleType("playwright.sync_api")
    fake_sync.sync_playwright = lambda: _FakePW()
    fake_playwright = types.ModuleType("playwright")
    fake_playwright.sync_api = fake_sync
    monkeypatch.setitem(sys.modules, "playwright", fake_playwright)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", fake_sync)

    result = brainrnaseq.download_csv_via_browser("mouse", gene_symbol="Genex")
    assert result.success is True
    assert result.data["browser_channel"] == "chromium"
    assert result.data["final_url"] == brainrnaseq.MOUSE_CSV_URL


def test_browser_download_both_channels_fail(monkeypatch):
    monkeypatch.setattr(
        brainrnaseq,
        "_launch_playwright_chromium",
        lambda pw: (_ for _ in ()).throw(RuntimeError("no browser")),
    )

    import sys
    import types

    class _FakePW:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    fake_sync = types.ModuleType("playwright.sync_api")
    fake_sync.sync_playwright = lambda: _FakePW()
    fake_playwright = types.ModuleType("playwright")
    fake_playwright.sync_api = fake_sync
    monkeypatch.setitem(sys.modules, "playwright", fake_playwright)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", fake_sync)

    # Patch to capture launch attempts then raise
    def _fail_launch(pw):
        attempts = [
            {"channel": "chrome", "success": False, "error_type": "Error"},
            {"channel": "chromium", "success": False, "error_type": "Error"},
        ]
        raise RuntimeError("launch failed")

    # Need launch to set attempts then raise — mimic real helper
    def _launch(pw):
        attempts = [
            {"channel": "chrome", "success": False},
            {"channel": "chromium", "success": False},
        ]
        exc = RuntimeError("both failed")
        # Attach attempts by raising after assigning via nonlocal pattern:
        # The real code catches and returns with attempts. Simulate raise path.
        raise RuntimeError("both failed")

    monkeypatch.setattr(brainrnaseq, "_launch_playwright_chromium", _launch)
    result = brainrnaseq.download_csv_via_browser("human", gene_symbol="GENEX")
    assert result.success is False
    assert result.error_type == "browser_launch_failed"
