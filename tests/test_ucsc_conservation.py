"""Focused tests for UCSC Section 1b conservation parsing and presentation."""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from gene_dossier.config import Settings
from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.normalize.ucsc_conservation import build_conservation_evidence
from gene_dossier.report_presentation import (
    UCSC_STABLE_INTRO,
    build_conservation_blocks,
    build_section_presentation,
)
from gene_dossier.report_schema import ReportContentBlock, build_report_document
from gene_dossier.ucsc_coords import (
    api_to_display_start,
    display_to_api_start,
    interval_from_api,
    interval_from_display,
)
from gene_dossier.ucsc_figure import (
    UCSC_SECTION_1B_TRACK_PRESET_ID,
    build_safe_hgtracks_url,
    is_safe_ucsc_browser_url,
    redact_api_key,
    sanitize_params,
    sha256_hex,
    validate_image_bytes,
)
from gene_dossier.tools import ucsc as ucsc_client
from gene_dossier.ucsc_parse import (
    filter_exact_gene_rows,
    parse_known_gene_region,
    parse_search_response,
    parse_transcript_row,
    parse_ucsc_int_array,
    reconstruct_exons,
    select_canonical_transcript,
)

FIXTURES = Path(__file__).parent / "fixtures" / "ucsc"
SEARCH_JSON = FIXTURES / "srebf2_search_relevant.json"
TRACK_JSON = FIXTURES / "srebf2_known_gene_region.json"
FIGURE_PNG = FIXTURES / "srebf2_comprehensive_conservation.png"
EXPECTED_FIGURE_SHA = "3d165b72c20d11a0c921d16bf2cd17418a5169c2d0cec0537e297de5be0e3d6a"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_search_parser_nested_groups_and_counts():
    search = _load(SEARCH_JSON)
    inv = parse_search_response(search, gene_symbol="SREBF2", genome="hg38")
    assert inv.basic_release == 50
    assert inv.comprehensive_release == 50
    assert inv.basic_exact_gene_count == 23
    assert inv.comprehensive_exact_gene_count == 45
    assert inv.mane_exact_gene_count >= 1
    assert any(m.group_kind == "mane" and m.transcript_id == "ENST00000361204.9" for m in inv.matches)
    # description-only / hubs rejected
    assert all(m.group_kind != "reject" for m in inv.matches)
    assert not any("SCAP" in (m.displayed_gene_symbol or "") for m in inv.matches)
    assert not any(m.displayed_gene_symbol.upper().startswith("SREBF2-AS1") for m in inv.matches)


def test_older_gencode_does_not_override_current():
    search = _load(SEARCH_JSON)
    inv = parse_search_response(search, gene_symbol="SREBF2", genome="hg38")
    assert inv.basic_release == 50
    assert inv.comprehensive_release == 50


def test_coordinate_helpers_no_off_by_one():
    assert display_to_api_start(41833105) == 41833104
    assert api_to_display_start(41833104) == 41833105
    disp = interval_from_display("chr22", 41833105, 41907305)
    assert disp.api_start_0_based == 41833104
    assert disp.api_end_exclusive == 41907305
    assert disp.display_position == "chr22:41833105-41907305"
    api = interval_from_api("chr22", 41833104, 41907305)
    assert api.display_position == "chr22:41833105-41907305"


def test_known_gene_exact_filter_and_neighbors():
    track = _load(TRACK_JSON)
    region = parse_known_gene_region(track, gene_symbol="SREBF2", genome="hg38")
    assert region.regional_total_count == 48
    assert region.exact_gene_transcript_count == 45
    assert region.excluded_neighbor_count == 3
    assert set(region.excluded_neighbor_symbols) >= {"SREBF2-AS1", "MIR33A", "SHISA8"}
    assert region.current_gencode_release == 50
    assert region.release_source == "track_big_data_url"
    assert "gencodeV50" in (region.track_big_data_url or "")


def test_biggenepred_exon_reconstruction_and_trailing_commas():
    assert parse_ucsc_int_array("1,2,3,") == [1, 2, 3]
    exons, ok, diags = reconstruct_exons(
        chrom_start=100,
        chrom_end=200,
        block_count=2,
        block_sizes="10,20,",
        chrom_starts="0,50,",
    )
    assert ok
    assert exons[0].start == 100 and exons[0].end == 110
    assert exons[1].start == 150 and exons[1].end == 170
    assert not diags

    _, bad, diags = reconstruct_exons(
        chrom_start=100,
        chrom_end=120,
        block_count=2,
        block_sizes="10,20",
        chrom_starts="0,50",
    )
    assert not bad
    assert any(d.code == "exon_out_of_range" for d in diags)

    _, bad2, diags2 = reconstruct_exons(
        chrom_start=100,
        chrom_end=200,
        block_count=3,
        block_sizes="10,20",
        chrom_starts="0,50",
    )
    assert not bad2
    assert any(d.code == "malformed_exon_arrays" for d in diags2)


def test_canonical_selection_precedence_and_tie():
    track = _load(TRACK_JSON)
    region = parse_known_gene_region(track, gene_symbol="SREBF2", genome="hg38")
    selected, diags = select_canonical_transcript(region.exact_rows)
    assert selected is not None
    assert selected.name == "ENST00000361204.9"
    assert selected.is_mane_select
    assert selected.is_ensembl_canonical
    assert selected.is_gencode_primary
    assert selected.is_canonical_tier
    assert selected.rank == 1
    assert selected.block_count == 19
    assert not any(d.code == "ambiguous_canonical_transcript" for d in diags)

    # Tie: two identical MANE-like rows
    a = parse_transcript_row(selected.raw)
    b = parse_transcript_row({**selected.raw, "name": "ENST99999999999.1"})
    assert a and b
    tied, tdiags = select_canonical_transcript([a, b])
    assert tied is None
    assert any(d.code == "ambiguous_canonical_transcript" for d in tdiags)

    # Without MANE, Ensembl_canonical wins
    rows = []
    for row in region.exact_rows:
        tr = parse_transcript_row(
            {
                **row.raw,
                "tag": (row.tag_raw or "").replace("MANE_Select", "").strip(","),
            }
        )
        if tr:
            rows.append(tr)
    sel2, _ = select_canonical_transcript(rows)
    assert sel2 is not None
    assert sel2.is_ensembl_canonical or sel2.is_gencode_primary or sel2.is_canonical_tier


def test_no_production_srebf2_fallback_constants():
    import gene_dossier.tools.ucsc as ucsc

    assert not hasattr(ucsc, "DEFAULT_REGION_SREBF2")
    assert not hasattr(ucsc, "DEFAULT_CANONICAL_TRANSCRIPT_SREBF2")


def test_evidence_contract_and_builder_dynamic_text():
    search = _load(SEARCH_JSON)
    track = _load(TRACK_JSON)
    png = FIGURE_PNG.read_bytes()
    assert sha256_hex(png) == EXPECTED_FIGURE_SHA
    validated, err = validate_image_bytes(png)
    assert err is None and validated is not None

    figure_value = {
        "figure_artifact_id": "art1",
        "relative_path": "run/ucsc/figures/x.png",
        "local_artifact_path": "run/ucsc/figures/x.png",
        "media_type": validated.media_type,
        "width": validated.width,
        "height": validated.height,
        "byte_size": validated.byte_size,
        "sha256": validated.sha256,
        "genome": "hg38",
        "display_position": "chr22:41833105-41907305",
        "selected_transcript": "ENST00000361204.9",
        "retrieval_method": "attached_validated_ucsc_render",
        "origin_endpoint": "hgRenderTracks",
        "api_key_used": True,
        "api_key_persisted": False,
    }
    records, diags = build_conservation_evidence(
        dossier_run_id="run",
        gene_symbol="SREBF2",
        genome="hg38",
        search_payload=search,
        track_payload=track,
        figure_value=figure_value,
    )
    types = {r.fact_type for r in records}
    assert types == {
        "ucsc_gene_locus",
        "ucsc_transcript_inventory",
        "ucsc_canonical_transcript",
        "ucsc_conservation_figure",
    }
    tx = next(r for r in records if r.fact_type == "ucsc_canonical_transcript")
    assert tx.value["transcript_id"] == "ENST00000361204.9"
    assert tx.value["api_start_0_based"] == 41833104
    assert tx.value["display_start_1_based"] == 41833105
    assert tx.value["exon_count"] == 19
    assert "apiKey" not in str(tx.value.get("browser_url"))

    inv = next(r for r in records if r.fact_type == "ucsc_transcript_inventory")
    assert inv.value["exact_gene_transcript_count"] == 45
    assert inv.value["current_gencode_release"] == "V50"

    # Presentation
    name_rec = EvidenceRecord(
        source_id="ncbi-gene:srebf2:gene_identity:6721",
        dossier_run_id="run",
        gene_symbol="SREBF2",
        section="General Gene Information",
        source_name="NCBI Gene",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="entrez_gene_id",
        evidence_grade=EvidenceGrade.C,
        taxon_id=9606,
        organism="Homo sapiens",
        value={
            "gene_name": "sterol regulatory element binding transcription factor 2",
            "nomenclaturename": "sterol regulatory element binding transcription factor 2",
        },
        display_text="name",
    )
    result = build_conservation_blocks(
        gene_symbol="SREBF2",
        evidence_records=records + [name_rec],
    )
    assert len([b for b in result.blocks if b.kind == "narrative"]) == 1
    narrative = result.blocks[0].text or ""
    assert UCSC_STABLE_INTRO in narrative
    assert "45 SREBF2" in narrative
    assert "GENCODE V50" in narrative
    assert "V44" not in narrative
    assert "51 results" not in narrative
    assert "both MANE Select and Ensembl canonical" in narrative
    assert result.blocks[1].kind == "link"
    assert "ENST00000361204.9" in (result.blocks[1].text or "")
    assert "chr22:41833105-41907305" in (result.blocks[1].text or "")
    url = result.blocks[1].links[0]["url"]
    assert is_safe_ucsc_browser_url(url)
    assert "apiKey" not in url

    # Figure omitted when path missing → diagnostic, no fabricated image
    assert any(b.kind == "figure" for b in result.blocks) or any(
        d.field == "figure" for d in result.diagnostics
    )


def test_safe_url_and_secret_redaction():
    url = build_safe_hgtracks_url(
        genome="hg38",
        display_position="chr22:41833105-41907305",
        transcript_id="ENST00000361204.9",
    )
    assert url and is_safe_ucsc_browser_url(url)
    assert not is_safe_ucsc_browser_url(url + "&apiKey=SECRET")
    assert not is_safe_ucsc_browser_url("https://evil.example/cgi-bin/hgTracks?db=hg38")
    assert "REDACTED" in redact_api_key("https://x?apiKey=super-secret&db=hg38")
    assert "apiKey" not in sanitize_params({"db": "hg38", "apiKey": "secret"})


def test_figure_validation_rejects_html_captcha_and_bootstrap():
    ok, err = validate_image_bytes(FIGURE_PNG.read_bytes())
    assert ok and err is None
    bad, err2 = validate_image_bytes(b"<html>cf-turnstile protecting itself from bots</html>")
    assert bad is None and err2 is not None
    bad2, err3 = validate_image_bytes(b"<html><script src='hgTracks.js'></script></html>")
    assert bad2 is None and err3 is not None
    bad3, err4 = validate_image_bytes(b"")
    assert bad3 is None and err4 is not None


class _FakeResponse:
    def __init__(self, *, status_code: int, content: bytes, text: str | None = None):
        self.status_code = status_code
        self.content = content
        self._text = text

    @property
    def is_success(self) -> bool:
        return 200 <= self.status_code < 300

    @property
    def text(self) -> str:
        return self._text if self._text is not None else self.content.decode("utf-8", errors="replace")


class _FakeClient:
    def __init__(self, responses: list[_FakeResponse]):
        self._responses = list(responses)
        self.calls: list[tuple[str, dict | None]] = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def close(self):
        return None

    def get(self, url, params=None):
        self.calls.append((str(url), dict(params) if params else None))
        return self._responses.pop(0)


def _synthetic_live_png(*, width: int = 1400, height: int = 400) -> bytes:
    """Build a minimal valid PNG with nontrivial IHDR dimensions (no Pillow)."""
    import struct
    import zlib

    def _chunk(tag: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + tag
            + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
        )

    # Single black row, filtered with filter type 0, repeated for height.
    raw = b"".join(b"\x00" + (b"\x00" * width * 3) for _ in range(height))
    return (
        b"\x89PNG\r\n\x1a\n"
        + _chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + _chunk(b"IDAT", zlib.compress(raw, 9))
        + _chunk(b"IEND", b"")
    )


def test_live_fetch_accepts_direct_image_without_serializing_bytes(monkeypatch, caplog):
    import logging

    png = _synthetic_live_png()
    fake = _FakeClient([_FakeResponse(status_code=200, content=png)])
    monkeypatch.setattr(ucsc_client.httpx, "Client", lambda timeout=None: fake)
    settings = Settings(ucsc_browser_api_key="fake-ucsc-secret-key-xyz")
    from gene_dossier.ucsc_figure import install_ucsc_api_key_log_redaction

    install_ucsc_api_key_log_redaction()
    logging.getLogger().setLevel(logging.DEBUG)
    with caplog.at_level(logging.DEBUG):
        logging.getLogger("httpx").info(
            "HTTP Request: GET https://genome.ucsc.edu/cgi-bin/hgRenderTracks?apiKey=%s",
            "fake-ucsc-secret-key-xyz",
        )
        logging.getLogger("httpx._client").info(
            "HTTP Request: GET https://genome.ucsc.edu/cgi-bin/hgRenderTracks?apiKey=fake-ucsc-secret-key-xyz"
        )
        logging.getLogger("httpcore").debug(
            "send_request_headers.complete return_value=apiKey=fake-ucsc-secret-key-xyz"
        )
        logging.getLogger("httpcore.connection").debug(
            "connect_tcp.started url=https://genome.ucsc.edu/"
            "cgi-bin/hgRenderTracks?apiKey=fake-ucsc-secret-key-xyz"
        )
        logging.getLogger("httpcore.http11").debug(
            "send_request_headers apiKey=%s",
            "fake-ucsc-secret-key-xyz",
        )
        meta, payload, attempts = ucsc_client.fetch_conservation_figure(
            gene_symbol="SREBF2",
            genome="hg38",
            display_position="chr22:41833105-41907305",
            selected_transcript="ENST00000361204.9",
            settings=settings,
        )
    assert meta is not None and meta["status"] == "ok"
    assert payload is not None
    assert payload.width == 1400
    assert "content=" not in repr(payload)
    assert meta["track_preset_id"] == UCSC_SECTION_1B_TRACK_PRESET_ID
    assert "apiKey" not in str(meta)
    assert attempts[-1]["endpoint_name"] == "hgRenderTracks"
    assert fake.calls[0][1]["apiKey"] == "fake-ucsc-secret-key-xyz"
    joined = "\n".join(r.getMessage() for r in caplog.records)
    assert "fake-ucsc-secret-key-xyz" not in joined
    assert "apiKey=REDACTED" in joined


def test_child_logger_api_key_redaction_at_debug(caplog):
    """Descendant httpcore/httpx loggers must redact even when filters are not inherited."""
    import logging

    from gene_dossier.ucsc_figure import install_ucsc_api_key_log_redaction

    install_ucsc_api_key_log_redaction()
    logging.getLogger().setLevel(logging.DEBUG)
    secret = "fake-secret"
    with caplog.at_level(logging.DEBUG):
        logging.getLogger("httpcore.connection").debug(
            "connect_tcp.started url=https://genome.ucsc.edu/"
            f"cgi-bin/hgRenderTracks?apiKey={secret}"
        )
        logging.getLogger("httpcore.http11").debug(
            "send_request_headers apiKey=%s",
            secret,
        )
        logging.getLogger("httpx._client").info(
            f"HTTP Request: GET https://example/?apiKey={secret}"
        )
        logging.getLogger("httpcore.proxy").debug(
            "proxy.connect apiKey=%s",
            secret,
        )
    joined = "\n".join(r.getMessage() for r in caplog.records)
    assert secret not in joined
    assert "apiKey=REDACTED" in joined
    assert all(secret not in r.getMessage() for r in caplog.records)


def test_live_fetch_accepts_html_wrapper_and_tracks_second_request(monkeypatch, caplog):
    import logging

    png = _synthetic_live_png()
    wrapper = b"<html><img src='logo.png'><img src='/trash/figure123.png?db=hg38&amp;x=1'></html>"
    fake = _FakeClient(
        [
            _FakeResponse(status_code=200, content=wrapper, text=wrapper.decode()),
            _FakeResponse(status_code=200, content=png),
        ]
    )
    monkeypatch.setattr(ucsc_client.httpx, "Client", lambda timeout=None: fake)
    settings = Settings(ucsc_browser_api_key="fake-ucsc-secret-key-xyz")
    from gene_dossier.ucsc_figure import install_ucsc_api_key_log_redaction

    install_ucsc_api_key_log_redaction()
    logging.getLogger().setLevel(logging.DEBUG)
    with caplog.at_level(logging.DEBUG):
        meta, payload, attempts = ucsc_client.fetch_conservation_figure(
            gene_symbol="SREBF2",
            genome="hg38",
            display_position="chr22:41833105-41907305",
            selected_transcript="ENST00000361204.9",
            settings=settings,
        )
    assert meta is not None and meta["status"] == "ok"
    assert meta["wrapper_used"] is True
    assert payload is not None
    assert payload.wrapper_request_index == 0
    assert attempts[0]["endpoint_name"] == "hgRenderTracks"
    assert attempts[1]["endpoint_name"] == "hgRenderTracks_image_asset"
    assert attempts[1]["request_url"] == "https://genome.ucsc.edu/trash/figure123.png"
    assert attempts[1]["request_params"] == {"db": "hg38", "x": "1"}
    assert attempts[1]["parent_request_index"] == 0
    joined = "\n".join(r.getMessage() for r in caplog.records)
    assert "fake-ucsc-secret-key-xyz" not in joined


def test_live_fetch_rejects_non_2xx_before_wrapper_parse(monkeypatch):
    fake = _FakeClient([_FakeResponse(status_code=403, content=b"<html><img src='/trash/x.png'></html>")])
    monkeypatch.setattr(ucsc_client.httpx, "Client", lambda timeout=None: fake)
    settings = Settings(ucsc_browser_api_key="secret")
    meta, payload, attempts = ucsc_client.fetch_conservation_figure(
        gene_symbol="SREBF2",
        genome="hg38",
        display_position="chr22:41833105-41907305",
        selected_transcript="ENST00000361204.9",
        settings=settings,
    )
    assert meta is not None and meta["status"] == "failed"
    assert meta["error_type"] == "http_error"
    assert payload is None
    assert len(fake.calls) == 1


def test_live_fetch_retries_image_asset_on_429(monkeypatch):
    png = _synthetic_live_png()
    wrapper = b"<html><img src='/trash/figure123.png'></html>"
    fake = _FakeClient(
        [
            _FakeResponse(status_code=200, content=wrapper, text=wrapper.decode()),
            _FakeResponse(status_code=429, content=b"retry"),
            _FakeResponse(status_code=200, content=png),
        ]
    )
    monkeypatch.setattr(ucsc_client.httpx, "Client", lambda timeout=None: fake)
    monkeypatch.setattr(ucsc_client, "sleep", lambda *_: None)
    settings = Settings(ucsc_browser_api_key="secret")
    meta, payload, attempts = ucsc_client.fetch_conservation_figure(
        gene_symbol="SREBF2",
        genome="hg38",
        display_position="chr22:41833105-41907305",
        selected_transcript="ENST00000361204.9",
        settings=settings,
    )
    assert meta is not None and meta["status"] == "ok"
    assert payload is not None
    assert sum(1 for a in attempts if a["endpoint_name"] == "hgRenderTracks_image_asset") == 2


def test_missing_key_is_nonfatal(monkeypatch):
    settings = Settings(ucsc_browser_api_key=None)
    meta, payload, attempts = ucsc_client.fetch_conservation_figure(
        gene_symbol="SREBF2",
        genome="hg38",
        display_position="chr22:41833105-41907305",
        selected_transcript="ENST00000361204.9",
        settings=settings,
    )
    assert meta is not None and meta["status"] == "failed"
    assert meta["error_type"] == "missing_key"
    assert payload is None
    assert attempts == []


def test_extract_prefers_trash_over_logo():
    from gene_dossier.ucsc_figure import extract_ucsc_image_url_from_html

    html = '<img src="/images/logo.png"><img src="/trash/hgt_abc.png?db=hg38">'
    assert extract_ucsc_image_url_from_html(html) == (
        "https://genome.ucsc.edu/trash/hgt_abc.png?db=hg38"
    )


def test_section_1b_block_order_and_no_audit_dump(tmp_path: Path):
    search = _load(SEARCH_JSON)
    track = _load(TRACK_JSON)
    # Place managed figure under a fake raw root via monkeypatch path in value absolute file
    fig_path = tmp_path / "fig.png"
    fig_path.write_bytes(FIGURE_PNG.read_bytes())
    records, _ = build_conservation_evidence(
        dossier_run_id="run",
        gene_symbol="SREBF2",
        genome="hg38",
        search_payload=search,
        track_payload=track,
        figure_value={
            "relative_path": str(fig_path),  # will fail portable check if absolute home-like
            "sha256": EXPECTED_FIGURE_SHA,
            "media_type": "image/png",
            "width": 80,
            "height": 40,
            "retrieval_method": "attached_validated_ucsc_render",
            "api_key_persisted": False,
        },
    )
    # Use absolute path under tmp — resolve_artifact_path may reject; force via patching value after
    for rec in records:
        if rec.fact_type == "ucsc_conservation_figure":
            rec.value["relative_path"] = str(fig_path)
            rec.value["local_artifact_path"] = str(fig_path)

    # Monkeypatch resolve to accept tmp file
    import gene_dossier.report_presentation as rp

    original = rp._resolve_figure_path

    def _fake_resolve(value):
        path = Path(value.get("relative_path"))
        if path.is_file() and sha256_hex(path.read_bytes()) == value.get("sha256"):
            return str(path), []
        return original(value)

    rp._resolve_figure_path = _fake_resolve  # type: ignore[assignment]
    try:
        result = build_section_presentation(
            section_key="1b",
            gene_symbol="SREBF2",
            evidence_records=records,
        )
    finally:
        rp._resolve_figure_path = original  # type: ignore[assignment]

    kinds = [b.kind for b in result.blocks]
    assert kinds[:3] == ["narrative", "link", "figure"]
    assert sum(1 for k in kinds if k == "narrative") == 1
    assert UCSC_STABLE_INTRO in (result.blocks[0].text or "")
    assert "transcript models" in (result.blocks[0].text or "")
    html_bits = " ".join((b.text or "") + (b.figure_caption or "") for b in result.blocks)
    assert "Key findings" not in html_bits
    assert "Limitations" not in html_bits
    assert "Supporting evidence" not in html_bits


def test_generic_dispatch_preserves_1a_and_adds_1b():
    search = _load(SEARCH_JSON)
    track = _load(TRACK_JSON)
    ucsc_records, _ = build_conservation_evidence(
        dossier_run_id="run",
        gene_symbol="SREBF2",
        genome="hg38",
        search_payload=search,
        track_payload=track,
        figure_value=None,
    )
    human = EvidenceRecord(
        source_id="ncbi-gene:srebf2:gene_identity:6721",
        dossier_run_id="run",
        gene_symbol="SREBF2",
        official_symbol="SREBF2",
        section="General Gene Information",
        source_name="NCBI Gene",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="entrez_gene_id",
        evidence_grade=EvidenceGrade.C,
        taxon_id=9606,
        organism="Homo sapiens",
        value={
            "entrez_gene_id": "6721",
            "gene_symbol": "SREBF2",
            "gene_name": "sterol regulatory element binding transcription factor 2",
            "nomenclaturename": "sterol regulatory element binding transcription factor 2",
            "aliases": ["BHLHD2"],
            "taxon_id": 9606,
        },
        display_text="SREBF2 Entrez Gene ID is 6721.",
    )
    doc = build_report_document(
        dossier_run_id="run",
        gene_symbol="SREBF2",
        evidence_records=ucsc_records + [human],
        report_sections=None,
    )
    major = next(s for s in doc.sections if s.key == "1")
    sub_a = next(s for s in major.subsections if s.key == "a")
    sub_b = next(s for s in major.subsections if s.key == "b")
    assert sub_a.presentation_blocks
    assert sub_a.presentation_blocks[0].presentation_role == "gene_aliases_table"
    assert sub_b.presentation_blocks
    assert sub_b.presentation_blocks[0].kind == "narrative"
    # unknown section unchanged: empty presentation stays empty
    sub_c = next(s for s in major.subsections if s.key == "c")
    assert sub_c.presentation_blocks == []


def test_inputs_not_mutated():
    search = _load(SEARCH_JSON)
    track = _load(TRACK_JSON)
    search_copy = json.loads(json.dumps(search))
    track_copy = json.loads(json.dumps(track))
    build_conservation_evidence(
        dossier_run_id="run",
        gene_symbol="SREBF2",
        genome="hg38",
        search_payload=search,
        track_payload=track,
        figure_value=None,
    )
    assert search == search_copy
    assert track == track_copy
