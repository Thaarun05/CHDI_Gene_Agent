"""Offline tests for DropViz Playwright acquisition + rank extraction helpers."""

from __future__ import annotations

import json
import sys
import types
from pathlib import Path

from gene_dossier.tools import dropviz

FIXTURES = Path(__file__).parent / "fixtures" / "dropviz"
RANK_ZIP = FIXTURES / "rank_ok.zip"
DUP_ZIP = FIXTURES / "rank_dup_basename.zip"
RAW_CSV = FIXTURES / "clusters_top_raw.csv"


class _FakeDownload:
    def __init__(self, url: str, payload: bytes):
        self.url = url
        self._payload = payload

    def save_as(self, path: str) -> None:
        Path(path).write_bytes(self._payload)


class _DownloadCtx:
    def __init__(self, download: _FakeDownload):
        self._download = download

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    @property
    def value(self):
        return self._download


class _FakeLocator:
    def __init__(self, page: "_FakePage", selector: str):
        self._page = page
        self._selector = selector

    def click(self, timeout=0):
        self._page._clicked.append(self._selector)

    def count(self):
        return 1

    def fill(self, value, timeout=0):
        self._page._filled = value

    @property
    def first(self):
        return self


class _FakePage:
    def __init__(
        self,
        *,
        url: str = "https://dropviz.org/?_state_id_=abc",
        inventory: dict | None = None,
        shiny_ready: dict | None = None,
        download_payload: bytes | None = None,
        download_url: str = "https://dropviz.org/session/tmp/download/rank.zip",
        goto_final_url: str | None = None,
        body_text: str = "DropViz Genex Levels By Cluster",
    ):
        self.url = url
        self._goto_final_url = goto_final_url or url
        self._handlers: list = []
        self._clicked: list[str] = []
        self._filled = None
        self._inventory = inventory or {
            "title": "DropViz",
            "text": body_text,
            "html_preview": "",
            "downloads": [
                {
                    "id": "gene.expr.rank.cluster.dl",
                    "href": "/session/tmp/download/rank.zip",
                    "text": "",
                    "className": "shiny-download-link",
                },
                {
                    "id": "tsne.global.cluster.label.dl",
                    "href": "/session/tmp/download/tsne.zip",
                    "text": "",
                    "className": "shiny-download-link",
                },
            ],
            "images": [
                {
                    "id": "plot1",
                    "src": "data:image/png;base64,xxx",
                    "complete": True,
                    "naturalWidth": 100,
                    "naturalHeight": 80,
                    "alt": "",
                }
            ],
            "body_classes": [],
            "url": url,
        }
        self._shiny_ready = shiny_ready or {
            "ready": True,
            "reason": "ok",
            "shinyOk": True,
            "restore_error_present": False,
            "enabled_download_count": 2,
            "enabled_download_ids": [
                "gene.expr.rank.cluster.dl",
                "tsne.global.cluster.label.dl",
            ],
            "dormant_download_count": 0,
            "gene_specific_plot_count": 1,
            "homepage_asset_count": 0,
            "requested_gene_visible": True,
            "homepage_only": False,
        }
        self._download_payload = download_payload
        self._download_url = download_url
        self._body_text = body_text

    def on(self, event, handler):
        self._handlers.append((event, handler))

    def goto(self, url, wait_until=None, timeout=None):
        if self._goto_final_url:
            self.url = self._goto_final_url
        else:
            self.url = url
        return None

    def evaluate(self, script, arg=None):
        text = str(script)
        if "dropviz-script: download-meta" in text:
            return {
                "present": True,
                "disabled": False,
                "href": self._download_url,
            }
        if "dropviz-script: shiny-ready" in text:
            return dict(self._shiny_ready)
        if "dropviz-script: page-inventory" in text:
            inv = dict(self._inventory)
            inv["url"] = self.url
            return inv
        return {}

    def wait_for_timeout(self, ms):
        return None

    def wait_for_function(self, *args, **kwargs):
        return None

    def screenshot(self, path=None, full_page=False):
        if path:
            Path(path).write_bytes(b"\x89PNG\r\n\x1a\nfake")
        return None

    def locator(self, selector):
        return _FakeLocator(self, selector)

    def expect_download(self, timeout=0):
        payload = self._download_payload
        if payload is None:
            payload = RANK_ZIP.read_bytes()
        return _DownloadCtx(_FakeDownload(self._download_url, payload))

    @property
    def keyboard(self):
        class _K:
            def type(self, *a, **k):
                return None

            def press(self, *a, **k):
                return None

        return _K()


class _FakeContext:
    def __init__(self, page: _FakePage):
        self._page = page

    def new_page(self):
        return self._page

    def close(self):
        return None


class _FakeBrowser:
    def __init__(self, page: _FakePage):
        self._page = page

    def new_context(self, accept_downloads=True):
        return _FakeContext(self._page)

    def close(self):
        return None


def _install_playwright(monkeypatch, page: _FakePage, launch_attempts=None):
    attempts = launch_attempts or [{"channel": "chrome", "success": True}]
    browser = _FakeBrowser(page)

    class _FakePW:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        @property
        def chromium(self):
            return object()

    fake_sync = types.ModuleType("playwright.sync_api")
    fake_sync.sync_playwright = lambda: _FakePW()
    fake_sync.Error = RuntimeError
    fake_playwright = types.ModuleType("playwright")
    fake_playwright.sync_api = fake_sync
    monkeypatch.setitem(sys.modules, "playwright", fake_playwright)
    monkeypatch.setitem(sys.modules, "playwright.sync_api", fake_sync)
    monkeypatch.setattr(
        dropviz,
        "_launch_playwright_chromium",
        lambda pw: (browser, list(attempts)),
    )
    return browser, attempts


def test_allowlist_and_redirects():
    assert dropviz.is_allowed_dropviz_url("http://dropviz.org/?_state_id_=x")
    assert dropviz.is_allowed_dropviz_url("https://www.dropviz.org/")
    assert not dropviz.is_allowed_dropviz_url("https://evil.example/")
    assert not dropviz.is_allowed_dropviz_url("ftp://dropviz.org/")
    assert dropviz.redirect_is_allowed("http://dropviz.org/", "https://dropviz.org/")
    assert dropviz.redirect_is_allowed(
        "http://dropviz.org/", "https://www.dropviz.org/"
    )
    assert not dropviz.redirect_is_allowed(
        "http://dropviz.org/", "https://evil.example/"
    )
    assert (
        dropviz.normalize_dropviz_url("http://www.dropviz.org/?_state_id_=abc")
        == "http://dropviz.org/?_state_id_=abc"
    )
    assert (
        dropviz.normalize_dropviz_url("https://www.dropviz.org/?_state_id_=abc")
        == "https://dropviz.org/?_state_id_=abc"
    )


def test_view_classification_never_uses_state_id():
    assert (
        dropviz.classify_view(download_link_ids=["gene.expr.rank.cluster.dl"])
        == dropviz.VIEW_RANK
    )
    assert (
        dropviz.classify_view(download_link_ids=["tsne.global.cluster.label.dl"])
        == dropviz.VIEW_GLOBAL_TSNE
    )
    assert (
        dropviz.classify_view(download_link_ids=["tsne.local.label.dl"])
        == dropviz.VIEW_REGIONAL_TSNE
    )
    assert (
        dropviz.classify_view(
            download_link_ids=[
                "gene.expr.rank.cluster.dl",
                "tsne.global.cluster.label.dl",
            ]
        )
        == dropviz.VIEW_MIXED
    )
    assert (
        dropviz.classify_view(
            download_link_ids=[],
            page_text="http://dropviz.org/?_state_id_=05190dfa61f331d8",
        )
        == dropviz.VIEW_UNKNOWN
    )


def test_css_escape_and_download_selector():
    assert "\\" in dropviz.css_escape_id("gene.expr.rank.cluster.dl")
    sel = dropviz.download_selector("gene.expr.rank.cluster.dl")
    assert sel.startswith("#")
    assert "gene" in sel


def test_zip_magic_and_duplicate_basename_rejected():
    assert dropviz.is_zip_bytes(RANK_ZIP.read_bytes())
    assert not dropviz.is_zip_bytes(b"<!DOCTYPE html>")
    assert not dropviz.is_html_payload(RANK_ZIP.read_bytes())
    assert dropviz.is_html_payload(b"<!DOCTYPE html><html></html>")
    ok = dropviz.inventory_zip_basenames(RANK_ZIP)
    assert ok["ok"] is True
    assert "rank.Rdata" in ok["basenames"]
    bad = dropviz.inventory_zip_basenames(DUP_ZIP)
    assert bad["ok"] is False
    assert bad["error_type"] == "ambiguous_duplicate_basename"


def test_rank_raw_vs_ranked_and_ci_validation(tmp_path):
    derived = dropviz.derive_rank_outputs_from_raw_csv(RAW_CSV, tmp_path)
    assert derived["ok"] is True
    ranked = (tmp_path / "clusters_top_ranked.csv").read_text(encoding="utf-8")
    lines = [ln for ln in ranked.splitlines() if ln.strip()]
    assert "Cluster A" in lines[1]
    assert "sort_policy" in derived
    top = json.loads((tmp_path / "top_clusters.json").read_text(encoding="utf-8"))
    assert top["clusters"][0]["label"] == "Cluster A"
    assert top["clusters"][0]["target.sum.per.100k"] == 12.5
    missing = dropviz.derive_rank_outputs_from_raw_csv(tmp_path / "nope.csv", tmp_path)
    assert missing["ok"] is False
    assert missing["status"] == "missing_clusters_top"
    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text(
        "cx.disp,target.sum.per.100k,target.sum.L.per.100k,target.sum.R.per.100k\n"
        "X,5,10,20\n",
        encoding="utf-8",
    )
    bad = dropviz.derive_rank_outputs_from_raw_csv(bad_csv, tmp_path / "bad_out")
    assert bad["ok"] is False


def test_regional_fail_closed_ambiguous_expression():
    amb = dropviz.assess_regional_quantitative_fields(
        ["cell", "V1", "V2", "region", "expr_a", "expr_b"]
    )
    assert amb["ok"] is False
    assert amb["regional_quantitative_status"] == "unavailable"
    assert amb["regional_evidence_type"] == "figure_derived"
    assert amb["reason"] == "ambiguous_expression_field"
    good = dropviz.assess_regional_quantitative_fields(
        ["cell", "V1", "V2", "region", "cluster", "target.sum.per.100k"]
    )
    assert good["ok"] is True
    assert good["regional_quantitative_status"] == "available"


def test_detect_expired_state():
    assert dropviz.detect_state_failure("Error: state not found") == "state_not_found"
    assert dropviz.detect_state_failure("saved state expired") == "state_expired"
    assert dropviz.detect_state_failure("Normal DropViz content") is None


def test_production_denylist_in_dropviz_module():
    src = Path(dropviz.__file__).read_text(encoding="utf-8")
    forbidden = [
        "05190dfa61f331d8",
        "5c4cbc26b012914c",
        "719bad29fe0f17fc",
        "3814b82b16caaf25",
        "SREBF2",
        "CDH10",
        "Srebf2",
        "Cdh10",
    ]
    for token in forbidden:
        assert token not in src, f"production denylist hit: {token}"
    assert "golden" not in src.lower()


def test_inspect_chrome_success(tmp_path, monkeypatch):
    page = _FakePage()
    _install_playwright(monkeypatch, page)
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="http://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path / "state",
    )
    assert result.success is True
    assert result.source_name == "DropViz"
    data = result.data
    assert data["status"] == "success"
    assert data["payload"]["acquisition_status"] == "success"
    assert data["payload"]["view_type"] == dropviz.VIEW_MIXED
    assert data["audit"]["browser_channel"] == "chrome"
    assert (tmp_path / "state" / "manifest.json").is_file()
    assert (tmp_path / "state" / "network.json").is_file()
    assert "temporary_download_hrefs" in data["audit"]
    hrefs = data["audit"]["temporary_download_hrefs"]
    assert any("session" in h for h in hrefs)


def test_inspect_chrome_fail_chromium_success(tmp_path, monkeypatch):
    page = _FakePage()
    _install_playwright(
        monkeypatch,
        page,
        launch_attempts=[
            {"channel": "chrome", "success": False, "error_type": "Error"},
            {"channel": "chromium", "success": True},
        ],
    )
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path / "state",
    )
    assert result.success is True
    assert result.data["audit"]["browser_channel"] == "chromium"


def test_inspect_host_not_allowed(tmp_path):
    client = dropviz.DropVizClient()
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://evil.example/?_state_id_=x",
        output_dir=tmp_path,
    )
    assert result.success is False
    assert result.error_type == "host_not_allowed"


def test_inspect_redirect_not_allowed(tmp_path, monkeypatch):
    page = _FakePage(goto_final_url="https://evil.example/phish")
    _install_playwright(monkeypatch, page)
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path,
    )
    assert result.success is False
    assert result.error_type == "redirect_not_allowed"


def test_shiny_busy_disconnected_blank_image(tmp_path, monkeypatch):
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)

    busy = _FakePage(shiny_ready={"ready": False, "reason": "shiny_busy"})
    _install_playwright(monkeypatch, busy)
    r1 = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path / "busy",
    )
    assert r1.success is False
    assert r1.data["status"] == "shiny_busy"

    disc = _FakePage(shiny_ready={"ready": False, "reason": "shiny_disconnected"})
    _install_playwright(monkeypatch, disc)
    r2 = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path / "disc",
    )
    assert r2.success is False
    assert r2.data["status"] == "shiny_disconnected"

    blank = _FakePage(
        shiny_ready={"ready": False, "reason": "no_plot_or_download"},
        inventory={
            "title": "DropViz",
            "text": "loading",
            "downloads": [],
            "images": [
                {
                    "id": "x",
                    "src": "",
                    "complete": True,
                    "naturalWidth": 0,
                    "naturalHeight": 0,
                    "alt": "",
                }
            ],
            "body_classes": [],
            "url": "https://dropviz.org/",
        },
    )
    _install_playwright(monkeypatch, blank)
    r3 = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path / "blank",
    )
    assert r3.success is False
    assert r3.data["status"] == "blank_or_unready"


def test_expired_state_detection(tmp_path, monkeypatch):
    page = _FakePage(
        body_text="Error: state not found for bookmark",
        inventory={
            "title": "DropViz",
            "text": "Error: state not found for bookmark",
            "downloads": [],
            "images": [],
            "body_classes": [],
            "url": "https://dropviz.org/",
        },
        shiny_ready={"ready": True, "reason": "ok"},
    )
    _install_playwright(monkeypatch, page)
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=dead",
        output_dir=tmp_path,
    )
    assert result.success is False
    assert result.data["status"] == "state_not_found"


HOMEPAGE_TEXT = (
    "DropViz Exploring the Mouse Brain through Single Cell Expression Profiles "
    "Get Started Explore Cell Types Compare and Discover Genes Search By Gene"
)

RESTORE_ERROR_TEXT = (
    "Error in RestoreContext initialization: Session 0123456789abcdef not found"
)

# All 13 dormant download anchors exist in the DropViz homepage DOM.
HOMEPAGE_DORMANT_DOWNLOADS = [
    {
        "id": "gene.expr.rank.cluster.dl",
        "href": "",
        "text": "",
        "className": "shiny-download-link disabled",
    },
    {
        "id": "tsne.global.cluster.label.dl",
        "href": "",
        "text": "",
        "className": "shiny-download-link disabled",
    },
    {
        "id": "tsne.local.label.dl",
        "href": "",
        "text": "",
        "className": "shiny-download-link disabled",
    },
]

HOMEPAGE_SAMPLE_IMAGES = [
    {
        "id": None,
        "src": "http://dropviz.org/tsne-sample.jpg",
        "complete": True,
        "naturalWidth": 963,
        "naturalHeight": 857,
        "alt": "",
    },
    {
        "id": None,
        "src": "http://dropviz.org/rank-sample.jpg",
        "complete": True,
        "naturalWidth": 798,
        "naturalHeight": 796,
        "alt": "",
    },
    {
        "id": None,
        "src": "http://dropviz.org/StanleyCenter-web.png",
        "complete": True,
        "naturalWidth": 413,
        "naturalHeight": 160,
        "alt": "Stanley Center",
    },
]


def _homepage_evidence(*, restore_error: bool) -> dict:
    """Readiness evidence the real predicate returns for a homepage render."""
    return {
        "ready": False,
        "reason": "restore_error" if restore_error else "homepage_only",
        "shinyOk": True,
        "restore_error_present": restore_error,
        "enabled_download_count": 0,
        "enabled_download_ids": [],
        "dormant_download_count": len(HOMEPAGE_DORMANT_DOWNLOADS),
        "gene_specific_plot_count": 0,
        "homepage_asset_count": len(HOMEPAGE_SAMPLE_IMAGES),
        "requested_gene_visible": False,
        "homepage_only": True,
    }


def test_detect_restore_context_session_not_found():
    assert dropviz.detect_restore_error(RESTORE_ERROR_TEXT) is True
    assert dropviz.detect_state_failure(RESTORE_ERROR_TEXT) == "state_not_found"
    assert dropviz.detect_restore_error(HOMEPAGE_TEXT) is False
    assert dropviz.detect_state_failure(HOMEPAGE_TEXT) is None


def test_restore_error_page_reports_state_not_found(tmp_path, monkeypatch):
    page = _FakePage(
        body_text=RESTORE_ERROR_TEXT,
        inventory={
            "title": "DropViz",
            "text": RESTORE_ERROR_TEXT + " " + HOMEPAGE_TEXT,
            "downloads": HOMEPAGE_DORMANT_DOWNLOADS,
            "images": HOMEPAGE_SAMPLE_IMAGES,
            "body_classes": [],
            "url": "http://dropviz.org/?_state_id_=dead",
        },
        shiny_ready=_homepage_evidence(restore_error=True),
    )
    _install_playwright(monkeypatch, page)
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="http://dropviz.org/?_state_id_=dead",
        output_dir=tmp_path,
    )

    assert result.success is False
    assert result.data["status"] == "state_not_found"
    payload = result.data["payload"]
    assert payload["acquisition_status"] == "state_not_found"
    assert payload["view_type"] == dropviz.VIEW_UNKNOWN
    acceptance = payload["acceptance"]
    assert acceptance["restore_error_present"] is True
    assert acceptance["requested_gene_visible"] is False
    assert acceptance["usable_for_section_2c"] is False

    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["state_failure"] == "state_not_found"
    assert manifest["view_type"] == dropviz.VIEW_UNKNOWN


def test_homepage_without_restore_error_is_not_success(tmp_path, monkeypatch):
    """Dormant anchors plus sample art must never classify as mixed_view."""
    page = _FakePage(
        body_text=HOMEPAGE_TEXT,
        inventory={
            "title": "DropViz",
            "text": HOMEPAGE_TEXT,
            "downloads": HOMEPAGE_DORMANT_DOWNLOADS,
            "images": HOMEPAGE_SAMPLE_IMAGES,
            "body_classes": [],
            "url": "http://dropviz.org/",
        },
        shiny_ready=_homepage_evidence(restore_error=False),
    )
    _install_playwright(monkeypatch, page)
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="http://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path,
    )

    assert result.success is False
    assert result.data["status"] == "blank_or_unready"
    payload = result.data["payload"]
    assert payload["view_type"] == dropviz.VIEW_UNKNOWN
    assert payload["enabled_download_link_ids"] == []
    assert payload["acceptance"]["homepage_only"] is True
    # The dormant anchors are still inventoried for audit.
    assert "gene.expr.rank.cluster.dl" in payload["download_link_ids"]


def test_acceptance_requires_gene_evidence():
    ready_evidence = {
        "ready": True,
        "restore_error_present": False,
        "requested_gene_visible": True,
        "homepage_only": False,
        "gene_specific_plot_count": 3,
        "enabled_download_count": 2,
    }
    accepted = dropviz.acceptance_from_evidence(
        ready_evidence, gene_query_submitted=True
    )
    assert accepted["usable_for_section_2c"] is True
    assert accepted["gene_specific_plot_count"] == 3
    assert accepted["gene_query_submitted"] is True

    # A state failure overrides otherwise-ready evidence.
    blocked = dropviz.acceptance_from_evidence(
        ready_evidence, state_failure="state_not_found"
    )
    assert blocked["usable_for_section_2c"] is False


def test_gene_acceptance_mirrors_the_accepted_attempt():
    expired = {
        "state_url": "http://dropviz.org/?_state_id_=05190dfa61f331d8",
        "status": "state_not_found",
        "acceptance": {
            "requested_gene_visible": False,
            "homepage_only": True,
            "restore_error_present": True,
            "gene_query_submitted": False,
            "gene_specific_plot_count": 0,
            "enabled_download_count": 0,
            "usable_for_section_2c": False,
        },
    }
    dynamic = {
        "state_url": "http://dropviz.org/",
        "status": "partial_success",
        "acceptance": {
            "requested_gene_visible": True,
            "homepage_only": False,
            "restore_error_present": False,
            "gene_query_submitted": True,
            "gene_specific_plot_count": 4,
            "enabled_download_count": 6,
            "usable_for_section_2c": True,
        },
    }

    rolled = dropviz._aggregate_acceptance([expired, dynamic])

    # The dead saved state must not stamp its restore error onto the collection
    # that Section 2c will actually read.
    assert rolled["restore_error_present"] is False
    assert rolled["homepage_only"] is False
    assert rolled["requested_gene_visible"] is True
    assert rolled["usable_for_section_2c"] is True
    assert rolled["accepted_attempt"] == "http://dropviz.org/"
    assert [a["state_url"] for a in rolled["rejected_attempts"]] == [
        "http://dropviz.org/?_state_id_=05190dfa61f331d8"
    ]


def test_gene_acceptance_keeps_restore_error_when_nothing_is_usable():
    expired = {
        "state_url": "http://dropviz.org/?_state_id_=05190dfa61f331d8",
        "status": "state_not_found",
        "acceptance": {
            "requested_gene_visible": False,
            "homepage_only": True,
            "restore_error_present": True,
            "gene_query_submitted": False,
            "gene_specific_plot_count": 0,
            "enabled_download_count": 0,
            "usable_for_section_2c": False,
        },
    }

    rolled = dropviz._aggregate_acceptance([expired])

    assert rolled["usable_for_section_2c"] is False
    assert rolled["restore_error_present"] is True
    assert rolled["accepted_attempt"] is None


def test_download_success_html_rejected_zip_magic(tmp_path):
    client = dropviz.DropVizClient()
    page = _FakePage(download_payload=RANK_ZIP.read_bytes())
    ok = client.download_shiny_export(
        page=page,
        link_id="gene.expr.rank.cluster.dl",
        output_path=tmp_path / "rank.zip",
        gene_symbol="Genex",
    )
    assert ok.success is True
    assert ok.data["payload"]["acquisition_status"] == "success"
    assert ok.data["payload"]["stable_endpoint"] is False
    assert ok.data["audit"]["session_href_is_stable_endpoint"] is False
    assert any("session" in h for h in ok.data["audit"]["temporary_download_hrefs"])

    html_page = _FakePage(
        download_payload=b"<!DOCTYPE html><html><body>x</body></html>"
    )
    bad_html = client.download_shiny_export(
        page=html_page,
        link_id="gene.expr.rank.cluster.dl",
        output_path=tmp_path / "bad.html.zip",
        gene_symbol="Genex",
    )
    assert bad_html.success is False
    assert bad_html.error_type == "html_download_rejected"

    magic_page = _FakePage(download_payload=b"NOTZIP")
    bad_magic = client.download_shiny_export(
        page=magic_page,
        link_id="gene.expr.rank.cluster.dl",
        output_path=tmp_path / "notzip.bin",
        gene_symbol="Genex",
    )
    assert bad_magic.success is False
    assert bad_magic.error_type == "invalid_zip_magic"

    dup_page = _FakePage(download_payload=DUP_ZIP.read_bytes())
    bad_dup = client.download_shiny_export(
        page=dup_page,
        link_id="gene.expr.rank.cluster.dl",
        output_path=tmp_path / "dup.zip",
        gene_symbol="Genex",
    )
    assert bad_dup.success is False
    assert bad_dup.error_type == "ambiguous_duplicate_basename"


def test_rscript_unavailable_partial_success_preserves_acquisition(
    tmp_path, monkeypatch
):
    page = _FakePage()
    _install_playwright(monkeypatch, page)
    monkeypatch.setattr(dropviz, "rscript_available", lambda: False)
    monkeypatch.setattr(
        dropviz,
        "run_extract_dropviz_rank",
        lambda **kwargs: {
            "ok": False,
            "status": "rscript_unavailable",
            "error_type": "rscript_unavailable",
            "error_message": "Rscript not found on PATH",
            "api_run": None,
        },
    )
    monkeypatch.setattr(
        dropviz,
        "run_inspect_dropviz_rdata",
        lambda **kwargs: {
            "ok": False,
            "status": "rscript_unavailable",
            "error_type": "rscript_unavailable",
            "api_run": None,
        },
    )
    result = dropviz.collect_dropviz_gene(
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "out",
        saved_state_urls=["https://dropviz.org/?_state_id_=abc"],
        client=dropviz.DropVizClient(shiny_ready_timeout_ms=250),
    )
    assert result.success is True
    assert result.data["status"] == "partial_success"
    assert result.data["payload"]["acquisition_status"] == "success"
    assert result.data["payload"]["rank_extraction_status"] == "rscript_unavailable"
    assert (tmp_path / "out" / "state_1" / "manifest.json").is_file()
    assert (tmp_path / "out" / "summary.json").is_file()
    dumped = json.dumps(result.data)
    assert "ApiRun(" not in dumped


def test_no_fake_api_run_in_extraction_wrappers(tmp_path, monkeypatch):
    class _Completed:
        returncode = 0
        stdout = "ok"
        stderr = ""

    def _fake_run(cmd, **kwargs):
        out = Path(cmd[3])
        out.mkdir(parents=True, exist_ok=True)
        (out / "clusters_top_raw.csv").write_text(RAW_CSV.read_text(encoding="utf-8"))
        return _Completed()

    monkeypatch.setattr(dropviz.subprocess, "run", _fake_run)
    monkeypatch.setattr(dropviz.shutil, "which", lambda name: "/usr/bin/Rscript")
    result = dropviz.run_extract_dropviz_rank(
        zip_or_rdata_path=RANK_ZIP,
        output_dir=tmp_path / "extract",
    )
    assert result.get("api_run") is None
    assert result["ok"] is True
    assert (tmp_path / "extract" / "clusters_top_ranked.csv").is_file()


def test_collect_dynamic_fallback_only_for_state_failures(tmp_path, monkeypatch):
    page = _FakePage()
    _install_playwright(monkeypatch, page)
    calls = {"dyn": 0}
    real_drive = dropviz.DropVizClient.drive_genex_dynamic_ui

    def _track(self, **kwargs):
        calls["dyn"] += 1
        return real_drive(self, **kwargs)

    monkeypatch.setattr(dropviz.DropVizClient, "drive_genex_dynamic_ui", _track)
    monkeypatch.setattr(
        dropviz,
        "run_extract_dropviz_rank",
        lambda **kwargs: {
            "ok": False,
            "status": "extraction_failed",
            "error_type": "extraction_failed",
            "api_run": None,
        },
    )
    monkeypatch.setattr(
        dropviz,
        "run_inspect_dropviz_rdata",
        lambda **kwargs: {"ok": True, "status": "success", "api_run": None},
    )
    result = dropviz.collect_dropviz_gene(
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "nofallback",
        saved_state_urls=["https://dropviz.org/?_state_id_=abc"],
        client=dropviz.DropVizClient(shiny_ready_timeout_ms=250),
    )
    assert calls["dyn"] == 0
    assert result.data["audit"]["used_dynamic_fallback"] is False
    assert result.success is True
    assert result.data["status"] == "partial_success"


def test_collect_dynamic_fallback_on_expired_states(tmp_path, monkeypatch):
    page = _FakePage(
        body_text="Error: state not found",
        inventory={
            "title": "DropViz",
            "text": "Error: state not found",
            "downloads": [],
            "images": [],
            "body_classes": [],
            "url": "https://dropviz.org/",
        },
        shiny_ready={"ready": True, "reason": "ok"},
    )
    _install_playwright(monkeypatch, page)
    calls = {"dyn": 0}

    def _drive(self, **kwargs):
        calls["dyn"] += 1
        out = Path(kwargs["output_dir"])
        out.mkdir(parents=True, exist_ok=True)
        return dropviz._tool_result(
            endpoint_name="collect_dynamic_query",
            gene_symbol=kwargs["mouse_gene_symbol"],
            request_url=dropviz.BASE_URL,
            success=True,
            data=dropviz._envelope(
                status="success",
                payload={
                    "gene_symbol": kwargs["mouse_gene_symbol"],
                    "view_type": dropviz.VIEW_RANK,
                    "artifacts": [],
                    "downloads": [],
                    "extractions": [],
                    "acquisition_status": "success",
                    "rank_extraction_status": "not_attempted",
                    "regional_quantitative_status": "unavailable",
                    "regional_evidence_type": "figure_derived",
                },
                audit={},
            ),
        )

    monkeypatch.setattr(dropviz.DropVizClient, "collect_dynamic_query", _drive)
    result = dropviz.collect_dropviz_gene(
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "fallback",
        saved_state_urls=["https://dropviz.org/?_state_id_=dead"],
        client=dropviz.DropVizClient(shiny_ready_timeout_ms=250),
        skip_extraction=True,
    )
    assert calls["dyn"] == 1
    assert result.data["audit"]["used_dynamic_fallback"] is True
    assert result.data["audit"]["fallback_reason"] == "all_saved_states_unavailable"


GENEX_TABLE_CSV = (
    b'"Region","Class","Cluster","Genex Amount","Genex P-Val"\n'
    b'"Substantia Nigra","Neuron","Neuron_Th [#4]",2.83,1.8e-71\n'
    b'"Frontal Cortex","Neuron","Neuron_Layer5b [#4]",2.77,4.4e-63\n'
    b'"Cerebellum","Astrocyte","Astrocyte_Gja1 [#8]",1.05,2.5e-9\n'
)


def test_summarize_cluster_table_csv_ranks_by_amount(tmp_path):
    csv_path = tmp_path / "clusters_table.csv"
    csv_path.write_bytes(GENEX_TABLE_CSV)

    summary = dropviz.summarize_cluster_table_csv(csv_path, "Genex")

    assert summary["ok"] is True
    assert summary["row_count"] == 3
    assert summary["amount_column"] == "Genex Amount"
    assert summary["label_column"] == "Cluster"
    # The table export carries no confidence bounds; say so rather than imply it.
    assert summary["confidence_intervals_available"] is False
    amounts = [row["amount"] for row in summary["top_clusters"]]
    assert amounts == sorted(amounts, reverse=True)
    assert summary["top_clusters"][0]["cluster"] == "Neuron_Th [#4]"
    assert summary["top_clusters"][0]["region"] == "Substantia Nigra"


def test_summarize_cluster_table_csv_reports_missing_gene_columns(tmp_path):
    csv_path = tmp_path / "other.csv"
    csv_path.write_bytes(b'"Region","Class","Cluster","Other Amount"\n"A","B","C",1\n')

    summary = dropviz.summarize_cluster_table_csv(csv_path, "Genex")

    assert summary["ok"] is False
    assert summary["status"] == "gene_columns_missing"


class _QueryLocator:
    """Locator over the modelled DropViz Query DOM."""

    def __init__(self, page: "_FakeQueryPage", selector: str, index: int = 0):
        self._page = page
        self._selector = selector
        self._index = index

    def _matches(self) -> list[str]:
        return self._page.matches(self._selector)

    def count(self):
        return len(self._matches())

    def click(self, timeout=0, force=False):
        self._page.click_selector(self._selector)

    def type(self, text, delay=0):
        self._page.type_gene(text)

    def fill(self, value, timeout=0):
        self._page.type_gene(value)

    def inner_text(self):
        items = self._matches()
        return items[self._index] if self._index < len(items) else ""

    def nth(self, index):
        return _QueryLocator(self._page, self._selector, index)

    def screenshot(self, path=None):
        if path:
            Path(path).write_bytes(b"\x89PNG\r\n\x1a\nplot")

    def wait_for(self, state=None, timeout=0):
        return None

    @property
    def first(self):
        return _QueryLocator(self._page, self._selector, 0)


class _FakeQueryPage:
    """Models the DropViz Query panel: selectize, Update!, tabs, plots."""

    def __init__(
        self,
        *,
        gene: str = "Genex",
        enable_downloads: bool = True,
        rendered_top_region: str | None = "Substantia Nigra",
        table_export_fails: bool = False,
    ):
        self.url = "http://dropviz.org/"
        self._gene = gene
        self._enable_downloads = enable_downloads
        self._rendered_top_region = rendered_top_region
        self._table_export_fails = table_export_fails
        self.query_tab_open = False
        self.typed = ""
        self.selected: list[str] = []
        self.submitted = False
        self.active_main = "clusters"
        self.active_sub = "rank"
        self.clicked: list[str] = []
        self.last_download_id: str | None = None
        self.focus: str | None = None
        self.region_options = ["Substantia Nigra", "Frontal Cortex", "Cerebellum"]
        self.region_selected: list[str] = []
        self._handlers: list = []

    # -- DOM model ---------------------------------------------------------
    def matches(self, selector: str) -> list[str]:
        if selector == dropviz.QUERY_TAB_SELECTOR:
            return ["Query"]
        if selector == dropviz.GENE_SELECTIZE_INPUT_SELECTOR:
            return ["input"] if self.query_tab_open else []
        if selector == dropviz.REGION_SELECTIZE_INPUT_SELECTOR:
            return ["input"] if self.query_tab_open else []
        if selector.endswith(".selectize-dropdown-content .option"):
            if not self.typed:
                return []
            if f'select[id="{dropviz.REGION_SELECT_ID}"]' in selector:
                return list(self.region_options)
            if f'select[id="{dropviz.GENE_SELECT_ID}"]' in selector:
                return [self._gene, f"{self._gene}os"]
            # Unscoped selector: the gene dropdown wins in the real DOM too.
            return [self._gene, f"{self._gene}os"]
        if f'select[id="{dropviz.GENE_SELECT_ID}"]' in selector:
            return list(self.selected)
        if f'select[id="{dropviz.REGION_SELECT_ID}"]' in selector:
            return list(self.region_selected)
        if selector == dropviz.GENE_SELECTED_ITEM_FALLBACK_SELECTOR:
            return list(self.selected) + list(self.region_selected)
        if "data-value=" in selector or selector.startswith("#mainpanel"):
            return ["tab"]
        if selector.startswith("button[id=") or selector == dropviz.QUERY_UPDATE_BUTTON_SELECTOR:
            return ["Update!"]
        if "img" in selector:
            return ["img"]
        return ["node"]

    def click_selector(self, selector: str) -> None:
        self.clicked.append(selector)
        if selector.startswith('a[id="'):
            self.last_download_id = selector.split('a[id="')[1].rstrip('"]')
        if selector in {dropviz.QUERY_TAB_SELECTOR, dropviz.QUERY_TAB_FALLBACK_SELECTOR}:
            self.query_tab_open = True
        elif selector == dropviz.GENE_SELECTIZE_INPUT_SELECTOR:
            self.focus = "gene"
            self.typed = ""
        elif selector == dropviz.REGION_SELECTIZE_INPUT_SELECTOR:
            self.focus = "region"
            self.typed = ""
        elif selector.endswith(".selectize-dropdown-content .option"):
            if f'select[id="{dropviz.REGION_SELECT_ID}"]' in selector:
                self.region_selected = [self.typed]
            else:
                self.selected = [self._gene]
        elif selector == dropviz.QUERY_UPDATE_BUTTON_SELECTOR:
            self.submitted = True
        elif 'data-value="' in selector:
            value = selector.split('data-value="')[1].split('"')[0]
            if selector.startswith("#mainpanel"):
                self.active_main = value
            else:
                self.active_sub = value

    def type_gene(self, text: str) -> None:
        self.typed = text

    # -- Playwright surface ------------------------------------------------
    def on(self, event, handler):
        self._handlers.append((event, handler))

    def goto(self, url, wait_until=None, timeout=None):
        self.url = url

    def locator(self, selector):
        return _QueryLocator(self, selector)

    def wait_for_timeout(self, ms):
        return None

    def wait_for_function(self, *args, **kwargs):
        return None

    def screenshot(self, path=None, full_page=False):
        if path:
            Path(path).write_bytes(b"\x89PNG\r\n\x1a\npage")

    def expect_download(self, timeout=0):
        page = self

        class _LazyCtx:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            @property
            def value(self):
                link = page.last_download_id or ""
                if link.startswith("dt."):
                    if page._table_export_fails:
                        raise RuntimeError("socket hang up")
                    payload = GENEX_TABLE_CSV
                    name = "table.csv"
                else:
                    payload = RANK_ZIP.read_bytes()
                    name = "rank.zip"
                return _FakeDownload(
                    f"http://dropviz.org/session/tok/download/{name}", payload
                )

        return _LazyCtx()

    @property
    def keyboard(self):
        page = self

        class _K:
            def type(self, text, **kwargs):
                page.type_gene(text)

            def press(self, key, **kwargs):
                if key == "Enter" and page.typed:
                    page.selected = [page._gene]

        return _K()

    def evaluate(self, script, arg=None):
        text = str(script)
        if "dropviz-script: download-meta" in text:
            return {
                "present": True,
                "disabled": not (self.submitted and self._enable_downloads),
                "href": (
                    "/session/tok/download/rank.zip"
                    if self.submitted and self._enable_downloads
                    else ""
                ),
            }
        if "dropviz-script: cluster-table-top-row" in text:
            if not (self.submitted and self._rendered_top_region):
                return None
            return {
                "headers": ["Region", "Class", "Cluster"],
                "cells": [self._rendered_top_region, "Neuron", "Neuron_Th [#4]"],
                "region": self._rendered_top_region,
                "cluster": "Neuron_Th [#4]",
                "cell_class": "Neuron",
            }
        if "dropviz-script: plot-image" in text:
            if not self.submitted:
                return {"rendered": False, "reason": "not_painted"}
            return {
                "rendered": True,
                "reason": "ok",
                "src": "http://dropviz.org/session/tok/img.png",
                "width": 900,
                "height": 700,
            }
        if "dropviz-script: shiny-ready" in text:
            return {
                "ready": self.submitted,
                "reason": "ok" if self.submitted else "no_gene_specific_output",
                "restore_error_present": False,
                "enabled_download_count": 3 if self.submitted else 0,
                "enabled_download_ids": (
                    [v["download_id"] for v in dropviz.DYNAMIC_VIEW_PLAN]
                    if self.submitted
                    else []
                ),
                "gene_specific_plot_count": 3 if self.submitted else 0,
                "homepage_asset_count": 0,
                "requested_gene_visible": bool(self.selected),
                "homepage_only": not self.submitted,
            }
        return {}


def test_genex_dynamic_query_workflow(tmp_path):
    """The client must drive Query -> selectize -> Update! -> all three views."""
    page = _FakeQueryPage(gene="Genex")
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.collect_dynamic_query(
        page=page,
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "dynamic_genex",
        skip_extraction=True,
        plot_timeout_ms=1_000,
    )

    assert page.query_tab_open is True
    assert page.selected == ["Genex"]
    assert page.submitted is True

    payload = result.data["payload"]
    assert payload["dynamic_ui"] is True
    assert payload["acceptance"]["gene_query_submitted"] is True
    assert payload["acceptance"]["requested_gene_visible"] is True

    captured = {a["view"] for a in payload["artifacts"] if a.get("kind") == "plot_image"}
    assert captured == {"rank", "tsne_global", "tsne_local"}
    assert {d["view"] for d in payload["downloads"]} == {
        "rank",
        "tsne_global",
        "tsne_local",
        "clusters_table",
        "subclusters_table",
    }
    assert all(d["success"] for d in payload["downloads"])

    out = tmp_path / "dynamic_genex"
    assert (out / "manifest.json").is_file()
    assert (out / "network.json").is_file()
    assert (out / "query_diagnostics.json").is_file()
    assert (out / "images" / "rank_plot.png").is_file()
    assert (out / "images" / "global_tsne.png").is_file()
    assert (out / "images" / "local_tsne.png").is_file()
    assert (out / "downloads" / "rank.zip").is_file()


class _FakeAPIResponse:
    def __init__(self, payload: bytes, status: int = 200):
        self._payload = payload
        self.status = status

    @property
    def ok(self):
        return 200 <= self.status < 300

    def body(self):
        return self._payload


class _FakeAPIRequest:
    def __init__(self, payload: bytes, status: int = 200):
        self._payload = payload
        self._status = status
        self.calls: list[str] = []

    def get(self, url, timeout=0):
        self.calls.append(url)
        return _FakeAPIResponse(self._payload, self._status)


class _FakeContextWithRequest:
    def __init__(self, api: _FakeAPIRequest):
        self.request = api


class _PopupDownloadPage(_FakePage):
    """Models target=_blank links whose download never fires on this page."""

    def __init__(self, api: _FakeAPIRequest, **kwargs):
        super().__init__(**kwargs)
        self.context = _FakeContextWithRequest(api)

    def expect_download(self, timeout=0):
        raise TimeoutError("download event never fired (popup)")


def test_download_falls_back_to_same_session_context_request(tmp_path):
    """A popup-swallowed download must still be captured in-session."""
    api = _FakeAPIRequest(RANK_ZIP.read_bytes())
    page = _PopupDownloadPage(
        api,
        url="http://dropviz.org/",
        download_url="/session/tmp/download/rank.zip",
    )
    client = dropviz.DropVizClient()

    result = client.download_shiny_export(
        page=page,
        link_id="gene.expr.rank.cluster.dl",
        output_path=tmp_path / "rank.zip",
        gene_symbol="Genex",
    )

    assert result.success is True
    assert result.request_params["capture_via"] == "context_request"
    assert api.calls == ["http://dropviz.org/session/tmp/download/rank.zip"]
    assert (tmp_path / "rank.zip").read_bytes().startswith(dropviz.ZIP_MAGIC)
    # The temporary session URL stays in audit only.
    assert "session/tmp" not in json.dumps(result.data["payload"])


def test_download_fallback_rejects_off_allowlist_href(tmp_path):
    api = _FakeAPIRequest(RANK_ZIP.read_bytes())
    page = _PopupDownloadPage(api, url="http://dropviz.org/")
    page._download_url = "https://evil.example.com/rank.zip"

    def _evaluate(script, arg=None):
        text = str(script)
        if "dropviz-script: download-meta" in text:
            return {
                "present": True,
                "disabled": False,
                "href": "https://evil.example.com/rank.zip",
            }
        return {}

    page.evaluate = _evaluate  # type: ignore[method-assign]

    client = dropviz.DropVizClient()
    result = client.download_shiny_export(
        page=page,
        link_id="gene.expr.rank.cluster.dl",
        output_path=tmp_path / "rank.zip",
        gene_symbol="Genex",
    )

    assert result.success is False
    assert result.error_type == "download_failed"
    assert api.calls == []


def test_local_tsne_requires_region_selection(tmp_path):
    """Without a region the local t-SNE only renders a placeholder, so skip it."""
    page = _FakeQueryPage(gene="Genex")
    page.region_options = []  # no region matches the top cluster

    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.collect_dynamic_query(
        page=page,
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "no_region",
        skip_extraction=True,
        plot_timeout_ms=500,
    )

    payload = result.data["payload"]
    captured = {a["view"] for a in payload["artifacts"] if a.get("kind") == "plot_image"}
    assert "tsne_local" not in captured
    assert "rank" in captured

    diagnostics = json.loads(
        (tmp_path / "no_region" / "query_diagnostics.json").read_text()
    )
    assert diagnostics["region_phase"]["region_selected"] is False
    local = next(v for v in diagnostics["views"] if v["view"] == "tsne_local")
    assert local["download_status"] == "requires_region_selection"


def test_region_phase_uses_top_cluster_region(tmp_path):
    page = _FakeQueryPage(gene="Genex")
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    client.collect_dynamic_query(
        page=page,
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "region",
        skip_extraction=True,
        plot_timeout_ms=500,
    )

    diagnostics = json.loads(
        (tmp_path / "region" / "query_diagnostics.json").read_text()
    )
    # GENEX_TABLE_CSV ranks Substantia Nigra highest.
    assert diagnostics["region_phase"]["region"] == "Substantia Nigra"
    assert diagnostics["region_phase"]["region_source"] == "clusters_table_csv"
    assert diagnostics["region_phase"]["region_selected"] is True
    assert page.region_selected == ["Substantia Nigra"]


def test_region_falls_back_to_rendered_table_when_export_drops(tmp_path):
    """The live server drops CSV exports intermittently; the DOM still ranks."""
    page = _FakeQueryPage(
        gene="Genex", rendered_top_region="Frontal Cortex", table_export_fails=True
    )
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    client.collect_dynamic_query(
        page=page,
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "region",
        skip_extraction=True,
        plot_timeout_ms=500,
    )

    diagnostics = json.loads(
        (tmp_path / "region" / "query_diagnostics.json").read_text()
    )
    assert diagnostics["region_phase"]["region"] == "Frontal Cortex"
    assert diagnostics["region_phase"]["region_source"] == "rendered_cluster_table"
    assert page.region_selected == ["Frontal Cortex"]


def test_region_phase_reports_when_no_ranking_is_available(tmp_path):
    page = _FakeQueryPage(
        gene="Genex", rendered_top_region=None, table_export_fails=True
    )
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    client.collect_dynamic_query(
        page=page,
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "region",
        skip_extraction=True,
        plot_timeout_ms=500,
    )

    diagnostics = json.loads(
        (tmp_path / "region" / "query_diagnostics.json").read_text()
    )
    assert diagnostics["region_phase"] == {"attempted": False}
    local = [v for v in diagnostics["views"] if v["view"] == "tsne_local"]
    assert local[0]["skipped"] == "no_region_available"


def test_dynamic_query_without_rendered_plots_is_not_success(tmp_path):
    page = _FakeQueryPage(gene="Genex")

    # Model a query that never submits, so no gene-specific plot renders.
    def _no_submit(selector: str) -> None:
        _FakeQueryPage.click_selector(page, selector)
        page.submitted = False

    page.click_selector = _no_submit  # type: ignore[method-assign]

    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.collect_dynamic_query(
        page=page,
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "dyn_fail",
        skip_extraction=True,
        plot_timeout_ms=500,
    )

    assert result.success is False
    assert result.data["status"] == "plot_unavailable"
    assert result.data["payload"]["acceptance"]["usable_for_section_2c"] is False


def test_drive_genex_alias_delegates(tmp_path):
    page = _FakeQueryPage(gene="Genex")
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.drive_genex_dynamic_ui(
        page=page,
        mouse_gene_symbol="Genex",
        output_dir=tmp_path / "genex_alias",
        skip_extraction=True,
    )
    assert result.data["payload"]["dynamic_ui"] is True
    assert result.data["payload"]["gene_symbol"] == "Genex"


def test_session_href_not_stable_endpoint(tmp_path):
    client = dropviz.DropVizClient()
    page = _FakePage(
        download_payload=RANK_ZIP.read_bytes(),
        download_url="https://dropviz.org/session/abc123/download/rank.zip",
    )
    result = client.download_shiny_export(
        page=page,
        link_id="gene.expr.rank.cluster.dl",
        output_path=tmp_path / "rank.zip",
        gene_symbol="Genex",
    )
    assert result.success is True
    assert result.data["payload"]["stable_endpoint"] is False
    assert result.data["audit"]["session_href_is_stable_endpoint"] is False
    payload_json = json.dumps(result.data["payload"])
    assert "session/abc123" not in payload_json


def test_image_and_download_inventory_in_manifest(tmp_path, monkeypatch):
    page = _FakePage()
    _install_playwright(monkeypatch, page)
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path,
    )
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert "gene.expr.rank.cluster.dl" in manifest["download_link_ids"]
    assert manifest["image_inventory"][0]["naturalWidth"] > 0
    assert result.data["payload"]["screenshot_status"] == "success"


def test_tool_result_envelope_shape(tmp_path, monkeypatch):
    page = _FakePage()
    _install_playwright(monkeypatch, page)
    client = dropviz.DropVizClient(shiny_ready_timeout_ms=250)
    result = client.inspect_saved_state(
        gene_symbol="Genex",
        state_url="https://dropviz.org/?_state_id_=abc",
        output_dir=tmp_path,
    )
    assert set(result.data.keys()) == {"status", "payload", "audit"}
    payload = result.data["payload"]
    for key in (
        "gene_symbol",
        "state_url",
        "view_type",
        "artifacts",
        "downloads",
        "extractions",
        "acquisition_status",
        "rank_extraction_status",
        "regional_quantitative_status",
        "regional_evidence_type",
    ):
        assert key in payload
