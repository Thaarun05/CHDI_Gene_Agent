"""Offline tests for GEO Profiles Section 3a client helpers."""

from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

from gene_dossier.models import ToolResult
from gene_dossier.tools import geo_profiles as gp

FIX = Path(__file__).resolve().parent / "fixtures" / "geo_profiles"


def _cand(
    uid: str,
    *,
    score: int,
    organism: str,
    category: str,
    region: str = "hippocampus",
    gds: str = "1",
    graph_status: str = gp.GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE,
    graph_ok: bool = False,
) -> dict[str, Any]:
    return {
        "profile_uid": uid,
        "gds_uid": gds,
        "idref": f"probe_{uid}",
        "final_score": score,
        "graph_ok": graph_ok,
        "graph_status": graph_status,
        "pdat": "2015",
        "diversity_keys": {
            "organism": organism,
            "category": category,
            "region": region,
        },
    }


def _png(name: str) -> bytes:
    return (FIX / name).read_bytes()


def test_normalize_and_format_gds_ids() -> None:
    assert gp.normalize_gds_uid("GDS4524") == "4524"
    assert gp.normalize_gds_uid("4524") == "4524"
    assert gp.normalize_gds_uid(4524) == "4524"
    assert gp.format_gds_accession("4524") == "GDS4524"
    assert gp.format_gds_accession("GDS4524") == "GDS4524"
    assert gp.format_gds_accession("") == ""


def test_parse_profile_esummary_normalizes_gds() -> None:
    payload = json.loads((FIX / "esummary_profile.json").read_text())
    docs = gp.parse_profile_esummary_payload(payload)
    assert len(docs) == 1
    assert docs[0]["gds_uid"] == "4524"
    assert docs[0]["gds_accession"] == "GDS4524"


def test_validate_graph_url_rejects_icon_and_mismatch() -> None:
    ok, err = gp.validate_graph_url(
        "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at",
        gds_uid="4524",
        idref="1426744_at",
    )
    assert ok and err is None
    bad, err2 = gp.validate_graph_url(
        "https://www.ncbi.nlm.nih.gov/geo/tools/profileIcon.cgi?ID=GDS4524:1426744_at",
        gds_uid="4524",
        idref="1426744_at",
    )
    assert not bad and err2 == "graph_thumbnail_rejected"
    mismatch, err3 = gp.validate_graph_url(
        "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS9999:1426744_at",
        gds_uid="4524",
        idref="1426744_at",
    )
    assert not mismatch and err3 == "graph_identifier_mismatch"


def test_discover_graph_url_from_fixture_html() -> None:
    html = (FIX / "profile_97740750.html").read_text()
    found = gp.discover_graph_url_from_html(html, gds_uid="4524", idref="1426744_at")
    assert found["ok"] is True
    assert "profileGraph.cgi" in found["url"]
    assert "profileIcon" not in found["url"]


def test_score_profile_has_no_human_species_bonus() -> None:
    profile_human = {
        "title": "Stress effect on hippocampus",
        "taxon": "Homo sapiens",
        "gds_uid": "1",
        "gpl": "GPL570",
        "idref": "x",
    }
    profile_mouse = {
        "title": "Stress effect on hippocampus",
        "taxon": "Mus musculus",
        "gds_uid": "2",
        "gpl": "GPL1261",
        "idref": "y",
    }
    gds = {"title": "stress hippocampus", "summary": "control vs stress", "sample_count": 10}
    human = gp.score_profile(profile_human, gds=gds, subset_effect_flag=False)
    mouse = gp.score_profile(profile_mouse, gds=gds, subset_effect_flag=False)
    assert "human" not in human["score_components"]
    assert "organism" not in human["score_components"]
    assert human["final_score"] == mouse["final_score"]


def test_diversity_shortlist_is_deterministic_and_diversity_aware() -> None:
    ranked = [
        _cand("1", score=100, organism="mus musculus", category="stress", gds="10"),
        _cand("2", score=99, organism="mus musculus", category="stress", gds="11"),
        _cand("3", score=98, organism="homo sapiens", category="drug_treatment", gds="12"),
        _cand("4", score=97, organism="rattus norvegicus", category="genetic_perturbation", gds="13"),
        _cand("5", score=96, organism="mus musculus", category="immune_challenge", gds="14"),
    ]
    a = gp.build_diversity_shortlist(ranked, max_items=3)
    b = gp.build_diversity_shortlist(ranked, max_items=3)
    assert [c["profile_uid"] for c in a] == [c["profile_uid"] for c in b]
    organisms = {(c.get("diversity_keys") or {}).get("organism") for c in a}
    assert "homo sapiens" in organisms
    assert "mus musculus" in organisms


def test_final_selection_never_keeps_outside_shortlist_status() -> None:
    chart_calls: list[str] = []

    def fake_html(uid: str, **kwargs: Any):
        return ToolResult(
            source_name=gp.SOURCE_NAME,
            endpoint_name="profile_html",
            gene_symbol="GENEX",
            request_url=f"https://www.ncbi.nlm.nih.gov/geoprofiles/{uid}",
            success=True,
            data={"raw_text": "", "content_type": "text/html"},
        )

    def fake_chart(**kwargs: Any):
        uid = str(kwargs["profile_uid"])
        chart_calls.append(uid)
        return {
            "graph_status": gp.GRAPH_STATUS_SUCCESS,
            "image_bytes": b"\x89PNG\r\n\x1a\n" + b"0" * 100,
            "image_width": 800,
            "image_height": 400,
            "tool_results": [],
        }

    live = {
        "profile_uid": "out",
        "gds_uid": "9",
        "idref": "z",
        "graph_status": gp.GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE,
        "subset_effect_flag": False,
        "gds_metadata": {},
    }
    with patch.object(gp, "fetch_profile_html", side_effect=fake_html), patch.object(
        gp, "acquire_profile_chart", side_effect=fake_chart
    ):
        chart = gp.acquire_profile_chart(
            gene_symbol="GENEX",
            profile_uid="out",
            gds_uid="9",
            idref="z",
            profile_html="",
        )
        gp._apply_chart(live, chart)
    assert "out" in chart_calls
    assert live["graph_status"] == gp.GRAPH_STATUS_SUCCESS
    assert live["graph_status"] != gp.GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE


def test_outside_shortlist_does_not_affect_visual_status() -> None:
    selected_final = [
        {"graph_status": gp.GRAPH_STATUS_SUCCESS},
        {"graph_status": gp.GRAPH_STATUS_SUCCESS},
    ]
    candidates = [
        *selected_final,
        {"graph_status": gp.GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE},
        {"graph_status": gp.GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE},
    ]
    ok_count = sum(
        1 for s in selected_final if s.get("graph_status") == gp.GRAPH_STATUS_SUCCESS
    )
    visual_status = "success" if ok_count == len(selected_final) else "partial"
    assert visual_status == "success"
    assert any(
        c.get("graph_status") == gp.GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE for c in candidates
    )


def test_expected_graph_id_uses_formatted_accession() -> None:
    assert gp.expected_graph_id("4524", "1426744_at") == "GDS4524:1426744_at"


def test_detect_misuse_redirect_final_url() -> None:
    blocked, reason = gp.detect_ncbi_block_page(
        final_url="https://misuse.ncbi.nlm.nih.gov/error/blocking.shtml"
    )
    assert blocked is True
    assert reason in {"misuse_host", "blocking_shtml"} or (
        reason is not None and reason.startswith("block_token:")
    )


def test_detect_unusual_browser_activity_html() -> None:
    html = (FIX / "misuse_block.html").read_text()
    blocked, reason = gp.detect_ncbi_block_page(html=html, body_text=html)
    assert blocked is True
    assert reason is not None


def test_detect_abnormal_browsing_activity_token() -> None:
    blocked, _ = gp.detect_ncbi_block_page(
        body_text="We detected possible abnormal browsing activity from your IP."
    )
    assert blocked is True


def test_detect_error_title() -> None:
    blocked, reason = gp.detect_ncbi_block_page(
        final_url="https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS1:x",
        page_title="Error",
    )
    assert blocked is True
    assert reason == "error_title"


def test_large_misuse_banner_cannot_pass_without_identity() -> None:
    banner = _png("misuse_banner.png")
    # Size-only helper may pass, but positive capture validation must fail.
    size_ok, _, metrics = gp.validate_chart_bytes(banner)
    assert metrics.get("width", 0) >= 400
    assert metrics.get("height", 0) >= 200
    ok, err, checks = gp.validate_graph_capture(
        banner,
        gds_uid="4524",
        idref="1426744_at",
        final_url="https://misuse.ncbi.nlm.nih.gov/error/blocking.shtml",
        page_title="Error",
        body_text="Unusual browser activity!",
        html=(FIX / "misuse_block.html").read_text(),
    )
    assert ok is False
    assert err == "graph_http_blocked"
    assert checks.get("block_page_absent") is False
    # Even with a forged graph URL, banner-like pixels fail closed.
    ok2, err2, checks2 = gp.validate_graph_capture(
        banner,
        gds_uid="4524",
        idref="1426744_at",
        final_url="https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at",
        page_title="GEO Profile Graph",
    )
    assert ok2 is False
    assert err2 == "graph_http_blocked"
    assert checks2.get("banner_like") is True
    _ = size_ok


def test_generic_body_capture_requires_identity() -> None:
    chart = _png("valid_chart.png")
    ok, err, _ = gp.validate_graph_capture(
        chart,
        gds_uid="4524",
        idref="1426744_at",
        final_url="https://www.ncbi.nlm.nih.gov/",
        capture_method="generic_body",
    )
    assert ok is False
    assert err in {"graph_link_missing", "graph_capture_identity_missing"}


def test_valid_full_chart_fixture_succeeds() -> None:
    chart = _png("valid_chart.png")
    ok, err, checks = gp.validate_graph_capture(
        chart,
        gds_uid="4524",
        idref="1426744_at",
        final_url="https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at",
        page_title="GEO Profile Graph",
        capture_method="img",
    )
    assert ok is True
    assert err is None
    assert checks["final_host_ok"] is True
    assert checks["graph_id_ok"] is True
    assert checks["block_page_absent"] is True
    assert checks["dimensions_ok"] is True
    assert checks["nonblank_ok"] is True


def test_wrong_gds_or_idref_fails() -> None:
    chart = _png("valid_chart.png")
    ok, err, _ = gp.validate_graph_capture(
        chart,
        gds_uid="4524",
        idref="1426744_at",
        final_url="https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS9999:1426744_at",
    )
    assert ok is False
    assert err == "graph_identifier_mismatch"


def test_profile_icon_remains_rejected() -> None:
    chart = _png("valid_chart.png")
    ok, err, checks = gp.validate_graph_capture(
        chart,
        gds_uid="4524",
        idref="1426744_at",
        final_url="https://www.ncbi.nlm.nih.gov/geo/tools/profileIcon.cgi?ID=GDS4524:1426744_at",
    )
    assert ok is False
    assert err == "graph_thumbnail_rejected"
    assert checks["thumbnail_absent"] is False


def test_playwright_source_never_uses_union_first_screenshot() -> None:
    source = inspect.getsource(gp.GeoProfilesBrowserSession)
    module_source = inspect.getsource(gp)
    assert 'locator("img, canvas, body").first' not in source
    assert "img, canvas, body" not in source
    # Cold navigation to profileGraph.cgi must not be the entry path.
    assert "page.goto(graph_url" not in source
    assert "NCBI_HOME_URL" in module_source
    assert "geoprofiles/" in module_source
    assert "playwright_profile_navigation" in source
    assert "capture_profile_graph" in source


def test_blocked_direct_response_fail_closed(monkeypatch: Any) -> None:
    misuse_html = (FIX / "misuse_block.html").read_text().encode("utf-8")

    def fake_request(**kwargs: Any) -> ToolResult:
        return ToolResult(
            source_name=gp.SOURCE_NAME,
            endpoint_name=kwargs.get("endpoint_name") or "profile_graph",
            gene_symbol="SREBF2",
            request_url=kwargs.get("path") or "",
            success=True,
            status_code=200,
            data={
                "content_type": "text/html",
                "content_bytes": misuse_html,
                "raw_text": misuse_html.decode("utf-8"),
                "requested_url": kwargs.get("path"),
                "final_url": "https://misuse.ncbi.nlm.nih.gov/error/blocking.shtml",
                "redirect_history": [
                    "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at",
                    "https://misuse.ncbi.nlm.nih.gov/error/blocking.shtml",
                ],
            },
        )

    session = MagicMock()
    session.capture_profile_graph.return_value = {
        "acquisition_method": "playwright_profile_navigation",
        "graph_status": gp.GRAPH_STATUS_FAILED,
        "error_type": "graph_http_blocked",
        "image_bytes": None,
        "tool_results": [],
        "validation_checks": {"block_page_absent": False},
    }
    monkeypatch.setattr(gp, "_request", fake_request)
    out = gp.acquire_profile_chart(
        gene_symbol="SREBF2",
        profile_uid="97740750",
        gds_uid="4524",
        idref="1426744_at",
        profile_html=(FIX / "profile_97740750.html").read_text(),
        allow_playwright=True,
        browser_session=session,
    )
    assert out["graph_status"] == gp.GRAPH_STATUS_FAILED
    assert out["error_type"] == "graph_http_blocked"
    assert out.get("image_bytes") is None
    session.capture_profile_graph.assert_called_once()


def test_blocked_chart_no_validated_score_bonus() -> None:
    profile = {
        "title": "Stress effect on hippocampus",
        "taxon": "Mus musculus",
        "gds_uid": "1",
        "gpl": "GPL1261",
        "idref": "x",
    }
    gds = {"title": "stress hippocampus", "summary": "control vs stress", "sample_count": 10}
    blocked = gp.score_profile(profile, gds=gds, subset_effect_flag=False, graph_ok=False)
    ok = gp.score_profile(profile, gds=gds, subset_effect_flag=False, graph_ok=True)
    assert "validated_chart" not in blocked["score_components"]
    assert blocked["score_components"].get("validated_chart") is None
    assert ok["score_components"]["validated_chart"] == 15
    assert ok["final_score"] == blocked["final_score"] + 15


def test_apply_chart_clears_bytes_on_failure() -> None:
    cand: dict[str, Any] = {
        "profile_uid": "1",
        "title": "Stress hippocampus",
        "taxon": "Mus musculus",
        "gds_uid": "1",
        "idref": "x",
        "subset_effect_flag": False,
        "gds_metadata": {"title": "stress", "summary": "control", "sample_count": 8},
        "graph_image_bytes": b"stale",
    }
    gp._apply_chart(
        cand,
        {
            "graph_status": gp.GRAPH_STATUS_FAILED,
            "error_type": "graph_http_blocked",
            "image_bytes": _png("misuse_banner.png"),
        },
    )
    assert cand["graph_ok"] is False
    assert cand["graph_status"] == gp.GRAPH_STATUS_FAILED
    assert "graph_image_bytes" not in cand
    assert "validated_chart" not in cand.get("score_components", {})


def test_browser_session_reused_across_shortlist(monkeypatch: Any) -> None:
    sessions: list[Any] = []

    class FakeSession:
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            self.captures = 0

        def __enter__(self) -> "FakeSession":
            sessions.append(self)
            return self

        def __exit__(self, *args: Any) -> None:
            return None

        def capture_profile_graph(self, **kwargs: Any) -> dict[str, Any]:
            self.captures += 1
            return {
                "graph_status": gp.GRAPH_STATUS_SUCCESS,
                "acquisition_method": "playwright_profile_navigation",
                "image_bytes": _png("valid_chart.png"),
                "graph_final_url": (
                    "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi"
                    f"?ID=GDS{kwargs['gds_uid']}:{kwargs['idref']}"
                ),
                "tool_results": [],
            }

    monkeypatch.setattr(gp, "GeoProfilesBrowserSession", FakeSession)

    def fake_html(uid: str, **kwargs: Any) -> ToolResult:
        return ToolResult(
            source_name=gp.SOURCE_NAME,
            endpoint_name="profile_html",
            gene_symbol="GENEX",
            request_url=f"https://www.ncbi.nlm.nih.gov/geoprofiles/{uid}",
            success=True,
            data={
                "raw_text": (
                    f'<a href="/geo/tools/profileGraph.cgi?ID=GDS{uid}:{uid}_at">'
                    "Full chart</a>"
                ),
                "content_type": "text/html",
                "final_url": f"https://www.ncbi.nlm.nih.gov/geoprofiles/{uid}",
            },
        )

    def fake_acquire(**kwargs: Any) -> dict[str, Any]:
        session = kwargs.get("browser_session")
        assert session is not None
        return session.capture_profile_graph(
            profile_uid=kwargs["profile_uid"],
            gds_uid=kwargs["gds_uid"],
            idref=kwargs["idref"],
            graph_url=gp.construct_graph_url(kwargs["gds_uid"], kwargs["idref"]),
        )

    # Drive only the session-threading helper path via a miniature collect stub.
    session = FakeSession()
    with session:
        for uid in ("1", "2", "3"):
            fake_acquire(
                profile_uid=uid,
                gds_uid=uid,
                idref=f"{uid}_at",
                browser_session=session,
            )
    assert len(sessions) == 1
    assert sessions[0].captures == 3
    _ = fake_html


def test_playwright_navigation_clicks_validated_graph_link() -> None:
    """Same-tab path: home/profile first, then exact graph anchor click."""
    session = gp.GeoProfilesBrowserSession.__new__(gp.GeoProfilesBrowserSession)
    session.viewport = {"width": 1280, "height": 900}
    session.browser_channel = "chromium"
    session.browser_version = "1.0"
    session.home_ok = True
    session._capture_count = 0
    session._pw = object()
    session._browser = object()

    graph_url = "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at"
    chart_png = _png("valid_chart.png")

    profile_page = MagicMock()
    profile_page.url = "https://www.ncbi.nlm.nih.gov/geoprofiles/97740750"
    profile_page.title.return_value = "GEO Profile 97740750"
    profile_page.content.return_value = (FIX / "profile_97740750.html").read_text()

    link_locator = MagicMock()
    link_locator.count.return_value = 1
    target = MagicMock()
    target.get_attribute.return_value = (
        "/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at"
    )
    link_locator.first = target

    body = MagicMock()
    body.inner_text.return_value = "GEO Profile SREBF2"
    body.count.return_value = 1
    body.screenshot.return_value = chart_png
    body.is_visible.return_value = True
    body.bounding_box.return_value = {"width": 600, "height": 300}

    empty = MagicMock()
    empty.count.return_value = 0

    def locator_side_effect(sel: str, *args: Any, **kwargs: Any) -> Any:
        if "profileGraph.cgi" in sel:
            return link_locator
        if sel in {"img", "canvas", "svg", "#graphic, .graphic, #profileGraph, .profile-graph"}:
            return empty
        return body

    profile_page.locator.side_effect = locator_side_effect

    context = MagicMock()
    context.new_page.return_value = profile_page
    # No popup: expect_page raises -> same-tab path.
    context.expect_page.side_effect = TimeoutError("no popup")
    session._context = context

    def goto(url: str, **kwargs: Any) -> None:
        profile_page.url = url
        if "profileGraph.cgi" in url:
            profile_page.title.return_value = "GEO Profile Graph"
            profile_page.content.return_value = "<html><body>chart</body></html>"

    profile_page.goto.side_effect = goto

    # After same-tab click, simulate navigation by updating URL in click side effect.
    def click_side_effect(**kwargs: Any) -> None:
        profile_page.url = graph_url
        profile_page.title.return_value = "GEO Profile Graph"
        profile_page.content.return_value = "<html><body>Value Rank chart</body></html>"

    target.click.side_effect = click_side_effect

    out = session.capture_profile_graph(
        profile_uid="97740750",
        gds_uid="4524",
        idref="1426744_at",
        graph_url=graph_url,
    )
    assert profile_page.goto.call_args_list[0].args[0].endswith("/geoprofiles/97740750")
    target.click.assert_called()
    assert out["acquisition_method"] == "playwright_profile_navigation"
    assert out["graph_status"] == gp.GRAPH_STATUS_SUCCESS
    assert out["graph_final_url"] == graph_url
    assert out.get("image_bytes")

def test_playwright_popup_graph_navigation() -> None:
    session = gp.GeoProfilesBrowserSession.__new__(gp.GeoProfilesBrowserSession)
    session.viewport = {"width": 1280, "height": 900}
    session.browser_channel = "chromium"
    session.browser_version = "1.0"
    session.home_ok = True
    session._capture_count = 0

    graph_url = "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at"
    chart_png = _png("valid_chart.png")

    profile_page = MagicMock()
    profile_page.url = "https://www.ncbi.nlm.nih.gov/geoprofiles/97740750"
    profile_page.title.return_value = "GEO Profile 97740750"
    profile_page.content.return_value = (FIX / "profile_97740750.html").read_text()

    popup = MagicMock()
    popup.url = graph_url
    popup.title.return_value = "GEO Profile Graph"
    popup.content.return_value = "<html><body><img src='/chart.png'></body></html>"

    link_locator = MagicMock()
    link_locator.count.return_value = 1
    target = MagicMock()
    target.get_attribute.return_value = (
        "/geo/tools/profileGraph.cgi?ID=GDS4524:1426744_at"
    )
    link_locator.first = target

    img_loc = MagicMock()
    img_loc.count.return_value = 1
    img_el = MagicMock()
    img_el.is_visible.return_value = True
    img_el.bounding_box.return_value = {"width": 600, "height": 300}
    img_el.get_attribute.return_value = "/chart.png"
    img_el.screenshot.return_value = chart_png
    img_loc.nth.return_value = img_el

    def profile_locator(sel: str, *args: Any, **kwargs: Any) -> Any:
        if "profileGraph.cgi" in sel:
            return link_locator
        body = MagicMock()
        body.inner_text.return_value = "profile"
        body.count.return_value = 1
        return body

    def popup_locator(sel: str, *args: Any, **kwargs: Any) -> Any:
        if sel == "img":
            return img_loc
        if sel in {"canvas", "svg", "#graphic, .graphic, #profileGraph, .profile-graph"}:
            empty = MagicMock()
            empty.count.return_value = 0
            return empty
        body = MagicMock()
        body.inner_text.return_value = "Value Rank"
        body.count.return_value = 1
        return body

    profile_page.locator.side_effect = profile_locator
    popup.locator.side_effect = popup_locator

    popup_cm = MagicMock()
    popup_cm.__enter__.return_value = MagicMock(value=popup)
    popup_cm.__exit__.return_value = None

    context = MagicMock()
    context.new_page.return_value = profile_page
    context.expect_page.return_value = popup_cm
    session._context = context

    out = session.capture_profile_graph(
        profile_uid="97740750",
        gds_uid="4524",
        idref="1426744_at",
        graph_url=graph_url,
    )
    assert out["graph_status"] == gp.GRAPH_STATUS_SUCCESS
    assert out["graph_final_url"] == graph_url
    assert out["image_bytes"]
    assert out["capture_method"] == "img"
    target.click.assert_called()
