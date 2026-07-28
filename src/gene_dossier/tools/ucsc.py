"""UCSC Genome Browser client for locus, track, and live figure retrieval."""

from __future__ import annotations

from dataclasses import dataclass, field
from time import sleep
from typing import Any
from urllib.parse import urlencode, urlparse

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult
from gene_dossier.ucsc_coords import display_to_api_end, display_to_api_start
from gene_dossier.ucsc_figure import (
    UCSC_SECTION_1B_TRACK_PRESET,
    UCSC_SECTION_1B_TRACK_PRESET_ID,
    UCSC_SECTION_1B_TRACK_PRESET_VERSION,
    build_safe_hgtracks_url,
    extract_ucsc_image_url_from_html,
    install_ucsc_api_key_log_redaction,
    redact_api_key,
    sanitize_params,
    split_url_for_provenance,
    validate_live_render_image_bytes,
)
from gene_dossier.ucsc_parse import (
    parse_known_gene_region,
    parse_position,
    parse_search_response,
    select_canonical_transcript,
)

SOURCE_NAME = "UCSC"
UCSC_API_BASE = "https://api.genome.ucsc.edu"
UCSC_BROWSER_BASE = "https://genome.ucsc.edu/cgi-bin/hgTracks"
UCSC_RENDER_URL = "https://genome.ucsc.edu/cgi-bin/hgRenderTracks"

DEFAULT_GENOME = "hg38"
DEFAULT_TRACK = "knownGene"

# Install once at import so -v / DEBUG never leaks apiKey values.
install_ucsc_api_key_log_redaction()


@dataclass(frozen=True)
class UCSCLiveFigurePayload:
    """Workflow-local live image payload. Must never enter serializable state."""

    content: bytes = field(repr=False)
    media_type: str
    width: int
    height: int
    sha256: str
    byte_size: int
    request_chain: tuple[dict[str, Any], ...] = field(repr=False)
    image_request_index: int
    track_preset_id: str
    track_preset_version: int
    track_params: dict[str, str]
    genome: str
    display_position: str
    selected_transcript: str | None
    wrapper_request_index: int | None = None


@dataclass(frozen=True)
class UCSCWorkflowExecution:
    """Serializable ToolResult plus workflow-local live figure payload."""

    tool_result: ToolResult
    live_figure: UCSCLiveFigurePayload | None = None


def _tool_result(
    *,
    endpoint_name: str,
    gene_symbol: str,
    request_url: str,
    request_params: dict[str, Any],
    success: bool,
    status_code: int | None = None,
    data: Any | None = None,
    error_type: str | None = None,
    error_message: str | None = None,
) -> ToolResult:
    """Build a uniform :class:`ToolResult` for this source."""
    return ToolResult(
        source_name=SOURCE_NAME,
        endpoint_name=endpoint_name,
        success=success,
        gene_symbol=gene_symbol,
        request_url=redact_api_key(request_url),
        request_params=sanitize_params(request_params),
        status_code=status_code,
        data=data,
        error_type=error_type,
        error_message=redact_api_key(error_message) if error_message else None,
    )


def _base_url(url: str) -> str:
    parsed = urlparse(url)
    if not parsed.scheme or not parsed.netloc:
        return url
    return f"{parsed.scheme}://{parsed.netloc}{parsed.path}"


def browser_url(
    chrom: str,
    display_start: int,
    display_end: int,
    *,
    genome: str = DEFAULT_GENOME,
    transcript_id: str | None = None,
) -> str:
    """Build a credential-free UCSC hgTracks URL for a display interval."""
    url = build_safe_hgtracks_url(
        genome=genome,
        display_position=f"{chrom}:{display_start}-{display_end}",
        transcript_id=transcript_id,
    )
    return url or (
        f"{UCSC_BROWSER_BASE}?"
        f"{urlencode({'db': genome, 'position': f'{chrom}:{display_start}-{display_end}'})}"
    )


def summarize_known_gene(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key knownGene track fields (not evidence)."""
    return {
        "name": row.get("name"),
        "chrom": row.get("chrom"),
        "chrom_start": row.get("chromStart", row.get("txStart")),
        "chrom_end": row.get("chromEnd", row.get("txEnd")),
        "tx_start": row.get("txStart"),
        "tx_end": row.get("txEnd"),
        "strand": row.get("strand"),
        "exon_starts": row.get("exonStarts"),
        "exon_ends": row.get("exonEnds"),
        "gene_name": row.get("geneName"),
        "gene_name2": row.get("geneName2"),
        "tag": row.get("tag"),
        "tier": row.get("tier"),
        "rank": row.get("rank"),
    }


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a UCSC API path and return :class:`ToolResult`."""
    url = f"{UCSC_API_BASE}/{path.lstrip('/')}"
    query = {k: str(v) for k, v in params.items() if v is not None}
    request_url = f"{url}?{urlencode(query)}" if query else url
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=query,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            status_code=response.status_code,
            data=payload,
            error_type="http_error",
            error_message=f"HTTP {response.status_code}",
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type="http_error",
            error_message=redact_api_key(str(exc)),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type=type(exc).__name__,
            error_message=redact_api_key(str(exc)),
        )


def _render_request(
    *,
    gene_symbol: str,
    settings: Settings,
    params: dict[str, Any],
    client: httpx.Client | None = None,
) -> tuple[httpx.Response | None, list[dict[str, Any]], str | None]:
    """Call hgRenderTracks with bounded retry and sanitized provenance."""
    secret = settings.ucsc_browser_api_key
    if secret is None or not secret.get_secret_value().strip():
        return None, [], "UCSC_BROWSER_API_KEY is not configured"

    outgoing_params = dict(params)
    outgoing_params["apiKey"] = secret.get_secret_value()
    sanitized_params = sanitize_params(params)
    attempts: list[dict[str, Any]] = []
    last_error: str | None = None

    owns_client = client is None
    session = client or httpx.Client(timeout=settings.http_timeout_seconds)
    try:
        for attempt_no in range(3):
            try:
                response = session.get(UCSC_RENDER_URL, params=outgoing_params)
            except httpx.HTTPError as exc:
                last_error = redact_api_key(str(exc))
                attempts.append(
                    {
                        "endpoint_name": "hgRenderTracks",
                        "request_url": UCSC_RENDER_URL,
                        "request_params": dict(sanitized_params),
                        "status_code": None,
                        "success": False,
                        "error_type": type(exc).__name__,
                        "error_message": last_error,
                        "parent_request_index": None,
                    }
                )
                break

            attempts.append(
                {
                    "endpoint_name": "hgRenderTracks",
                    "request_url": UCSC_RENDER_URL,
                    "request_params": dict(sanitized_params),
                    "status_code": response.status_code,
                    "success": response.is_success,
                    "error_type": None if response.is_success else "http_error",
                    "error_message": None if response.is_success else f"HTTP {response.status_code}",
                    "parent_request_index": None,
                }
            )
            if response.status_code == 429 or 500 <= response.status_code < 600:
                last_error = f"HTTP {response.status_code}"
                if attempt_no < 2:
                    sleep(0.5 * (2**attempt_no))
                    continue
            return response, attempts, last_error
        return None, attempts, last_error or "hgRenderTracks request failed"
    finally:
        if owns_client:
            session.close()


def _fetch_image_asset(
    *,
    client: httpx.Client,
    image_url: str,
    wrapper_request_index: int,
    attempts: list[dict[str, Any]],
) -> httpx.Response | None:
    """GET a UCSC image asset with the same transient retry policy as hgRenderTracks."""
    base_url, asset_params = split_url_for_provenance(image_url)
    last_resp: httpx.Response | None = None
    for attempt_no in range(3):
        try:
            response = client.get(image_url)
        except httpx.HTTPError as exc:
            attempts.append(
                {
                    "endpoint_name": "hgRenderTracks_image_asset",
                    "request_url": base_url,
                    "request_params": dict(asset_params),
                    "status_code": None,
                    "success": False,
                    "error_type": type(exc).__name__,
                    "error_message": redact_api_key(str(exc)),
                    "parent_request_index": wrapper_request_index,
                }
            )
            return None
        attempts.append(
            {
                "endpoint_name": "hgRenderTracks_image_asset",
                "request_url": base_url,
                "request_params": dict(asset_params),
                "status_code": response.status_code,
                "success": response.is_success,
                "error_type": None if response.is_success else "http_error",
                "error_message": None if response.is_success else f"HTTP {response.status_code}",
                "parent_request_index": wrapper_request_index,
            }
        )
        last_resp = response
        if response.status_code == 429 or 500 <= response.status_code < 600:
            if attempt_no < 2:
                sleep(0.5 * (2**attempt_no))
                continue
        return response
    return last_resp


def _figure_failure(
    *,
    error_type: str,
    error_message: str,
    genome: str,
    display_position: str,
    selected_transcript: str | None,
    params: dict[str, Any],
) -> dict[str, Any]:
    return {
        "status": "failed",
        "error_type": error_type,
        "error_message": redact_api_key(error_message),
        "genome": genome,
        "display_position": display_position,
        "selected_transcript": selected_transcript,
        "track_preset_id": UCSC_SECTION_1B_TRACK_PRESET_ID,
        "track_preset_version": UCSC_SECTION_1B_TRACK_PRESET_VERSION,
        "track_params": sanitize_params(params),
        "pixel_width": int(UCSC_SECTION_1B_TRACK_PRESET["pix"]),
    }


def _figure_success_meta(
    *,
    validated: Any,
    genome: str,
    display_position: str,
    selected_transcript: str | None,
    params: dict[str, Any],
    wrapper_used: bool,
) -> dict[str, Any]:
    return {
        "status": "ok",
        "media_type": validated.media_type,
        "width": validated.width,
        "height": validated.height,
        "sha256": validated.sha256,
        "byte_size": validated.byte_size,
        "genome": genome,
        "display_position": display_position,
        "selected_transcript": selected_transcript,
        "retrieval_method": "programmatic_browser_render",
        "origin_endpoint": "hgRenderTracks",
        "api_key_used": True,
        "api_key_persisted": False,
        "track_preset_id": UCSC_SECTION_1B_TRACK_PRESET_ID,
        "track_preset_version": UCSC_SECTION_1B_TRACK_PRESET_VERSION,
        "track_params": sanitize_params(params),
        "pixel_width": int(UCSC_SECTION_1B_TRACK_PRESET["pix"]),
        "wrapper_used": wrapper_used,
    }


def fetch_conservation_figure(
    *,
    gene_symbol: str,
    genome: str,
    display_position: str,
    selected_transcript: str | None,
    settings: Settings,
) -> tuple[dict[str, Any] | None, UCSCLiveFigurePayload | None, list[dict[str, Any]]]:
    """Fetch a live UCSC conservation figure and return metadata + transient bytes."""
    params = {
        "db": genome,
        "position": display_position,
        **UCSC_SECTION_1B_TRACK_PRESET,
    }
    requested_pix = int(UCSC_SECTION_1B_TRACK_PRESET["pix"])
    with httpx.Client(timeout=settings.http_timeout_seconds) as client:
        response, attempts, render_error = _render_request(
            gene_symbol=gene_symbol,
            settings=settings,
            params=params,
            client=client,
        )
        if response is None:
            return (
                _figure_failure(
                    error_type="missing_key"
                    if render_error and "not configured" in render_error
                    else "http_error",
                    error_message=render_error or "hgRenderTracks request failed",
                    genome=genome,
                    display_position=display_position,
                    selected_transcript=selected_transcript,
                    params=params,
                ),
                None,
                attempts,
            )

        if not response.is_success:
            return (
                _figure_failure(
                    error_type="http_error",
                    error_message=f"HTTP {response.status_code}",
                    genome=genome,
                    display_position=display_position,
                    selected_transcript=selected_transcript,
                    params=params,
                ),
                None,
                attempts,
            )

        validated, image_error = validate_live_render_image_bytes(
            response.content, requested_pix=requested_pix
        )
        if validated is not None:
            payload = UCSCLiveFigurePayload(
                content=validated.content,
                media_type=validated.media_type,
                width=validated.width,
                height=validated.height,
                sha256=validated.sha256,
                byte_size=validated.byte_size,
                request_chain=tuple(attempts),
                image_request_index=len(attempts) - 1,
                track_preset_id=UCSC_SECTION_1B_TRACK_PRESET_ID,
                track_preset_version=UCSC_SECTION_1B_TRACK_PRESET_VERSION,
                track_params=sanitize_params(params),
                genome=genome,
                display_position=display_position,
                selected_transcript=selected_transcript,
                wrapper_request_index=None,
            )
            return (
                _figure_success_meta(
                    validated=validated,
                    genome=genome,
                    display_position=display_position,
                    selected_transcript=selected_transcript,
                    params=params,
                    wrapper_used=False,
                ),
                payload,
                attempts,
            )

        if image_error is None or image_error.code != "html_figure_wrapper":
            return (
                _figure_failure(
                    error_type=image_error.code if image_error else "invalid_figure",
                    error_message=image_error.message if image_error else "Invalid figure response",
                    genome=genome,
                    display_position=display_position,
                    selected_transcript=selected_transcript,
                    params=params,
                ),
                None,
                attempts,
            )

        wrapper_request_index = len(attempts) - 1
        image_url = extract_ucsc_image_url_from_html(response.text)
        if not image_url:
            return (
                _figure_failure(
                    error_type="html_wrapper_missing_image",
                    error_message="HTML wrapper did not contain an approved UCSC image URL",
                    genome=genome,
                    display_position=display_position,
                    selected_transcript=selected_transcript,
                    params=params,
                ),
                None,
                attempts,
            )

        image_resp = _fetch_image_asset(
            client=client,
            image_url=image_url,
            wrapper_request_index=wrapper_request_index,
            attempts=attempts,
        )
        if image_resp is None or not image_resp.is_success:
            return (
                _figure_failure(
                    error_type="http_error",
                    error_message=(
                        f"Image asset HTTP {image_resp.status_code}"
                        if image_resp is not None
                        else "Image asset request failed"
                    ),
                    genome=genome,
                    display_position=display_position,
                    selected_transcript=selected_transcript,
                    params=params,
                ),
                None,
                attempts,
            )

        validated, image_error = validate_live_render_image_bytes(
            image_resp.content, requested_pix=requested_pix
        )
        if validated is None:
            return (
                _figure_failure(
                    error_type=image_error.code if image_error else "invalid_figure",
                    error_message=image_error.message if image_error else "Invalid image asset bytes",
                    genome=genome,
                    display_position=display_position,
                    selected_transcript=selected_transcript,
                    params=params,
                ),
                None,
                attempts,
            )

        payload = UCSCLiveFigurePayload(
            content=validated.content,
            media_type=validated.media_type,
            width=validated.width,
            height=validated.height,
            sha256=validated.sha256,
            byte_size=validated.byte_size,
            request_chain=tuple(attempts),
            image_request_index=len(attempts) - 1,
            track_preset_id=UCSC_SECTION_1B_TRACK_PRESET_ID,
            track_preset_version=UCSC_SECTION_1B_TRACK_PRESET_VERSION,
            track_params=sanitize_params(params),
            genome=genome,
            display_position=display_position,
            selected_transcript=selected_transcript,
            wrapper_request_index=wrapper_request_index,
        )
        return (
            _figure_success_meta(
                validated=validated,
                genome=genome,
                display_position=display_position,
                selected_transcript=selected_transcript,
                params=params,
                wrapper_used=True,
            ),
            payload,
            attempts,
        )


def search(
    gene_symbol: str,
    *,
    genome: str = DEFAULT_GENOME,
    settings: Settings | None = None,
) -> ToolResult:
    """Search UCSC for ``gene_symbol`` on ``genome``."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    params = {"genome": genome, "search": symbol}
    return _request_json(
        endpoint_name="search",
        gene_symbol=symbol,
        path="search",
        params=params,
        settings=cfg,
    )


def get_track_data(
    chrom: str,
    api_start_0_based: int,
    api_end_exclusive: int,
    *,
    gene_symbol: str = "",
    genome: str = DEFAULT_GENOME,
    track: str = DEFAULT_TRACK,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch UCSC track data for an API (0-based half-open) interval."""
    cfg = settings or get_settings()
    if api_start_0_based is None or api_end_exclusive is None:
        return _tool_result(
            endpoint_name="get_track_data",
            gene_symbol=gene_symbol,
            request_url=f"{UCSC_API_BASE}/getData/track",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="UCSC track requests require both start and end",
        )
    params = {
        "genome": genome,
        "track": track,
        "chrom": chrom,
        "start": int(api_start_0_based),
        "end": int(api_end_exclusive),
    }
    return _request_json(
        endpoint_name="get_track_data",
        gene_symbol=gene_symbol or chrom,
        path="getData/track",
        params=params,
        settings=cfg,
    )


def fetch_gene_region_execution(
    gene_symbol: str,
    *,
    genome: str = DEFAULT_GENOME,
    track: str = DEFAULT_TRACK,
    settings: Settings | None = None,
) -> UCSCWorkflowExecution:
    """Search → track → transcript context → live figure, with workflow-local bytes.

    The returned `tool_result.data` is JSON-serializable. Live figure bytes, when
    present, are carried only in `UCSCWorkflowExecution.live_figure`.
    """
    cfg = settings or get_settings()
    searched = search(gene_symbol, genome=genome, settings=cfg)
    if not searched.success:
        return UCSCWorkflowExecution(
            tool_result=_tool_result(
                endpoint_name="fetch_gene_region",
                gene_symbol=gene_symbol,
                request_url=searched.request_url,
                request_params=searched.request_params,
                success=False,
                status_code=searched.status_code,
                data={"search": searched.data},
                error_type=searched.error_type or "search_failed",
                error_message=searched.error_message or "UCSC search failed",
            )
        )

    payload = searched.data if isinstance(searched.data, dict) else {}
    inventory = parse_search_response(payload, gene_symbol=gene_symbol, genome=genome)
    display = inventory.selected_display_interval
    if display is None:
        return UCSCWorkflowExecution(
            tool_result=_tool_result(
                endpoint_name="fetch_gene_region",
                gene_symbol=gene_symbol,
                request_url=searched.request_url,
                request_params=searched.request_params,
                success=True,
                status_code=searched.status_code,
                data={
                    "gene_symbol": gene_symbol,
                    "genome": genome,
                    "selection_method": "ambiguous",
                    "search": searched.data,
                    "chrom": None,
                    "start": None,
                    "end": None,
                    "api_start_0_based": None,
                    "api_end_exclusive": None,
                    "display_start_1_based": None,
                    "display_end_1_based": None,
                    "track_data": None,
                    "browser_url": None,
                },
            )
        )

    api_start = display_to_api_start(display.display_start_1_based)
    api_end = display_to_api_end(display.display_end_1_based)
    track_res = get_track_data(
        display.chrom,
        api_start,
        api_end,
        gene_symbol=gene_symbol,
        genome=genome,
        track=track,
        settings=cfg,
    )
    browser = browser_url(
        display.chrom,
        display.display_start_1_based,
        display.display_end_1_based,
        genome=genome,
    )
    if not track_res.success:
        return UCSCWorkflowExecution(
            tool_result=_tool_result(
                endpoint_name="fetch_gene_region",
                gene_symbol=gene_symbol,
                request_url=track_res.request_url,
                request_params=track_res.request_params,
                success=False,
                status_code=track_res.status_code,
                data={
                    "gene_symbol": gene_symbol,
                    "genome": genome,
                    "selection_method": "matched",
                    "search": searched.data,
                    "chrom": display.chrom,
                    "start": api_start,
                    "end": api_end,
                    "api_start_0_based": api_start,
                    "api_end_exclusive": api_end,
                    "display_start_1_based": display.display_start_1_based,
                    "display_end_1_based": display.display_end_1_based,
                    "position": display.display_position,
                    "track_data": track_res.data,
                    "browser_url": browser,
                },
                error_type=track_res.error_type or "track_failed",
                error_message=track_res.error_message or "UCSC track data failed",
            )
        )

    track_payload = track_res.data if isinstance(track_res.data, dict) else {}
    known = track_payload.get(track) or track_payload.get("knownGene") or []
    if not isinstance(known, list):
        known = []
    transcript_summaries = [
        summarize_known_gene(row) for row in known if isinstance(row, dict)
    ]
    region = parse_known_gene_region(
        track_payload,
        gene_symbol=gene_symbol,
        genome=genome,
        search_inventory=inventory,
    )
    selected_transcript_row, _sel_diags = select_canonical_transcript(region.exact_rows)
    selected_transcript_id = selected_transcript_row.name if selected_transcript_row else None
    figure_meta, live_figure, figure_attempts = fetch_conservation_figure(
        gene_symbol=gene_symbol,
        genome=genome,
        display_position=display.display_position,
        selected_transcript=selected_transcript_id,
        settings=cfg,
    )
    data = {
        "gene_symbol": gene_symbol,
        "genome": genome,
        "selection_method": "matched",
        "search": searched.data,
        "chrom": display.chrom,
        "start": api_start,
        "end": api_end,
        "api_start_0_based": api_start,
        "api_end_exclusive": api_end,
        "display_start_1_based": display.display_start_1_based,
        "display_end_1_based": display.display_end_1_based,
        "position": display.display_position,
        "display_position": display.display_position,
        "track_data": track_res.data,
        "transcript_summaries": transcript_summaries,
        "transcript_count": len(transcript_summaries),
        "selected_transcript": (
            {
                "name": selected_transcript_row.name,
                "gene_name": selected_transcript_row.gene_name,
                "chrom": selected_transcript_row.chrom,
                "chrom_start": selected_transcript_row.chrom_start,
                "chrom_end": selected_transcript_row.chrom_end,
                "tag": selected_transcript_row.tag_raw,
                "tier": selected_transcript_row.tier_raw,
                "rank": selected_transcript_row.rank,
            }
            if selected_transcript_row is not None
            else None
        ),
        "browser_url": browser,
        "figure": figure_meta if figure_meta and figure_meta.get("status") == "ok" else None,
        "figure_failure": figure_meta if figure_meta and figure_meta.get("status") != "ok" else None,
        "figure_attempts": figure_attempts,
    }
    return UCSCWorkflowExecution(
        tool_result=_tool_result(
        endpoint_name="fetch_gene_region",
        gene_symbol=gene_symbol,
        request_url=track_res.request_url,
        request_params={
            "genome": genome,
            "track": track,
            "chrom": display.chrom,
            "start": api_start,
            "end": api_end,
            "selection_method": "matched",
        },
        success=True,
        status_code=track_res.status_code,
        data=data,
        ),
        live_figure=live_figure,
    )


def fetch_gene_region(
    gene_symbol: str,
    *,
    genome: str = DEFAULT_GENOME,
    track: str = DEFAULT_TRACK,
    settings: Settings | None = None,
) -> ToolResult:
    """Compatibility wrapper returning only the serializable ToolResult."""
    return fetch_gene_region_execution(
        gene_symbol,
        genome=genome,
        track=track,
        settings=settings,
    ).tool_result


__all__ = [
    "SOURCE_NAME",
    "UCSC_API_BASE",
    "UCSC_BROWSER_BASE",
    "UCSC_RENDER_URL",
    "DEFAULT_GENOME",
    "DEFAULT_TRACK",
    "UCSCLiveFigurePayload",
    "UCSCWorkflowExecution",
    "parse_position",
    "browser_url",
    "summarize_known_gene",
    "fetch_conservation_figure",
    "fetch_gene_region_execution",
    "search",
    "get_track_data",
    "fetch_gene_region",
]
