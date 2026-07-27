"""UCSC Genome Browser API client.

Searches for a gene region, fetches knownGene track data, and builds browser
URLs. Conservation evidence normalization lives in
``normalize/ucsc_conservation.py``.

Key endpoints (validated)::

    GET https://api.genome.ucsc.edu/search?genome=hg38&search={symbol}
    GET https://api.genome.ucsc.edu/getData/track
        ?genome=hg38&track=knownGene&chrom={chrom}&start={api_start}&end={api_end}

Browser view (not JSON)::

    https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={display}

Track requests use 0-based half-open API coordinates derived from the
1-based display position via ``ucsc_coords`` helpers. Never reuse the display
start as the track API start.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult
from gene_dossier.ucsc_coords import display_to_api_end, display_to_api_start
from gene_dossier.ucsc_figure import (
    build_safe_hgtracks_url,
    redact_api_key,
    sanitize_params,
)
from gene_dossier.ucsc_parse import parse_position, parse_search_response

SOURCE_NAME = "UCSC"
UCSC_API_BASE = "https://api.genome.ucsc.edu"
UCSC_BROWSER_BASE = "https://genome.ucsc.edu/cgi-bin/hgTracks"
UCSC_RENDER_URL = "https://genome.ucsc.edu/cgi-bin/hgRenderTracks"

DEFAULT_GENOME = "hg38"
DEFAULT_TRACK = "knownGene"


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


def fetch_gene_region(
    gene_symbol: str,
    *,
    genome: str = DEFAULT_GENOME,
    track: str = DEFAULT_TRACK,
    settings: Settings | None = None,
) -> ToolResult:
    """Search → resolve display locus → fetch knownGene with API coordinates.

    On success, ``data`` includes search payload, track payload, and display /
    API coordinate fields. If no exact-gene positional match exists, returns
    success with ``selection_method="ambiguous"`` and no track fetch.

    Never raises.
    """
    cfg = settings or get_settings()
    searched = search(gene_symbol, genome=genome, settings=cfg)
    if not searched.success:
        return _tool_result(
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

    payload = searched.data if isinstance(searched.data, dict) else {}
    inventory = parse_search_response(payload, gene_symbol=gene_symbol, genome=genome)
    display = inventory.selected_display_interval
    if display is None:
        return _tool_result(
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
        return _tool_result(
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

    track_payload = track_res.data if isinstance(track_res.data, dict) else {}
    known = track_payload.get(track) or track_payload.get("knownGene") or []
    if not isinstance(known, list):
        known = []
    transcript_summaries = [
        summarize_known_gene(row) for row in known if isinstance(row, dict)
    ]

    return _tool_result(
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
            "transcript_summaries": transcript_summaries,
            "transcript_count": len(transcript_summaries),
            "browser_url": browser,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "UCSC_API_BASE",
    "UCSC_BROWSER_BASE",
    "UCSC_RENDER_URL",
    "DEFAULT_GENOME",
    "DEFAULT_TRACK",
    "parse_position",
    "browser_url",
    "summarize_known_gene",
    "search",
    "get_track_data",
    "fetch_gene_region",
]
