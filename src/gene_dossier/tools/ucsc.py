"""UCSC Genome Browser API client.

Searches for a gene region, fetches knownGene track data, and builds browser
URLs. Does **not** normalize into evidence records — that belongs in
``normalize/gene_identity.py``.

Key endpoints (validated)::

    GET https://api.genome.ucsc.edu/search?genome=hg38&search={symbol}
    GET https://api.genome.ucsc.edu/getData/track
        ?genome=hg38&track=knownGene&chrom={chrom}&start={start}&end={end}

Browser view (not JSON)::

    https://genome.ucsc.edu/cgi-bin/hgTracks
        ?db=hg38&position={chrom}:{start}-{end}&knownGene=pack&cons100way=full

NOTE: Track requests require both ``start`` and ``end``.

For SREBF2 / hg38, validated region ``chr22:41833105-41907305``
(canonical transcript ``ENST00000361204.9``).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import re
from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "UCSC"
UCSC_API_BASE = "https://api.genome.ucsc.edu"
UCSC_BROWSER_BASE = "https://genome.ucsc.edu/cgi-bin/hgTracks"

DEFAULT_GENOME = "hg38"
DEFAULT_TRACK = "knownGene"

# Validated SREBF2 / hg38 anchors.
DEFAULT_REGION_SREBF2 = "chr22:41833105-41907305"
DEFAULT_CANONICAL_TRANSCRIPT_SREBF2 = "ENST00000361204.9"

_POSITION_RE = re.compile(
    r"^(?P<chrom>chr[\w.]+):(?P<start>\d+)-(?P<end>\d+)$",
    re.IGNORECASE,
)


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
        request_url=request_url,
        request_params=request_params,
        status_code=status_code,
        data=data,
        error_type=error_type,
        error_message=error_message,
    )


def parse_position(position: str) -> dict[str, Any] | None:
    """Parse ``chrom:start-end`` into ``{chrom, start, end}``."""
    text = (position or "").strip().replace(",", "")
    match = _POSITION_RE.match(text)
    if not match:
        return None
    start = int(match.group("start"))
    end = int(match.group("end"))
    if end < start:
        return None
    return {
        "chrom": match.group("chrom"),
        "start": start,
        "end": end,
        "position": f"{match.group('chrom')}:{start}-{end}",
    }


def browser_url(
    chrom: str,
    start: int,
    end: int,
    *,
    genome: str = DEFAULT_GENOME,
    known_gene: str = "pack",
    cons100way: str = "full",
) -> str:
    """Build a UCSC hgTracks browser URL for the region (not a JSON API)."""
    params = {
        "db": genome,
        "position": f"{chrom}:{start}-{end}",
        "knownGene": known_gene,
        "cons100way": cons100way,
    }
    return f"{UCSC_BROWSER_BASE}?{urlencode(params)}"


def summarize_position_match(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key UCSC search hit fields (not evidence)."""
    position = row.get("position")
    parsed = parse_position(str(position)) if position else None
    return {
        "db": row.get("db"),
        "name": row.get("name"),
        "description": row.get("description"),
        "position": position,
        "chrom": parsed["chrom"] if parsed else None,
        "start": parsed["start"] if parsed else None,
        "end": parsed["end"] if parsed else None,
    }


def summarize_known_gene(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key knownGene track fields (not evidence)."""
    return {
        "name": row.get("name"),
        "chrom": row.get("chrom"),
        "tx_start": row.get("txStart"),
        "tx_end": row.get("txEnd"),
        "strand": row.get("strand"),
        "exon_starts": row.get("exonStarts"),
        "exon_ends": row.get("exonEnds"),
        "gene_name": row.get("geneName"),
        "gene_name2": row.get("geneName2"),
    }


def prefer_position_match(
    matches: list[Any],
    gene_symbol: str,
    *,
    genome: str = DEFAULT_GENOME,
) -> dict[str, Any] | None:
    """Prefer a search hit whose name matches ``gene_symbol`` on ``genome``.

    Does not blindly trust the first hit. Returns ``None`` when no safe match.
    """
    target = gene_symbol.strip().upper()
    if not target:
        return None
    scored: list[tuple[int, dict[str, Any]]] = []
    for row in matches:
        if not isinstance(row, dict):
            continue
        position = row.get("position")
        if not position or not parse_position(str(position)):
            continue
        name = str(row.get("name") or "").strip().upper()
        db = str(row.get("db") or "").strip()
        if name != target and target not in name.split():
            # Also allow description mentioning the symbol as weaker signal only
            # when name is empty — still require exact name preference first.
            continue
        rank = 0
        if db and db != genome:
            rank += 2
        if name != target:
            rank += 1
        scored.append((rank, row))
    if not scored:
        return None
    scored.sort(key=lambda item: item[0])
    best_rank = scored[0][0]
    best = [row for rank, row in scored if rank == best_rank]
    if len(best) != 1:
        return None
    return best[0]


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
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
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
    start: int,
    end: int,
    *,
    gene_symbol: str = "",
    genome: str = DEFAULT_GENOME,
    track: str = DEFAULT_TRACK,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch UCSC track data for a region (``start`` and ``end`` both required)."""
    cfg = settings or get_settings()
    if start is None or end is None:
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
        "start": int(start),
        "end": int(end),
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
    """Search → select coordinates → fetch knownGene track + browser URL.

    On success, ``data`` includes selected chrom/start/end, track summaries, and
    a browser URL. If no safe position match exists, returns success with
    ``selection_method="ambiguous"`` and no track fetch (does not guess).

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
    matches = payload.get("positionMatches") or []
    if not isinstance(matches, list):
        matches = []
    match_summaries = [
        summarize_position_match(m) for m in matches if isinstance(m, dict)
    ]
    selected = prefer_position_match(matches, gene_symbol, genome=genome)
    if selected is None:
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
                "position_match_summaries": match_summaries,
                "chrom": None,
                "start": None,
                "end": None,
                "track_data": None,
                "browser_url": None,
            },
        )

    parsed = parse_position(str(selected.get("position")))
    if not parsed:
        return _tool_result(
            endpoint_name="fetch_gene_region",
            gene_symbol=gene_symbol,
            request_url=searched.request_url,
            request_params=searched.request_params,
            success=False,
            status_code=searched.status_code,
            data={
                "gene_symbol": gene_symbol,
                "search": searched.data,
                "selected_match": selected,
            },
            error_type="parse_error",
            error_message="Could not parse UCSC position from selected match",
        )

    chrom = parsed["chrom"]
    start = parsed["start"]
    end = parsed["end"]
    track_res = get_track_data(
        chrom,
        start,
        end,
        gene_symbol=gene_symbol,
        genome=genome,
        track=track,
        settings=cfg,
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
                "selected_match": selected,
                "chrom": chrom,
                "start": start,
                "end": end,
                "track_data": track_res.data,
                "browser_url": browser_url(chrom, start, end, genome=genome),
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
            "chrom": chrom,
            "start": start,
            "end": end,
            "selection_method": "matched",
        },
        success=True,
        status_code=track_res.status_code,
        data={
            "gene_symbol": gene_symbol,
            "genome": genome,
            "selection_method": "matched",
            "search": searched.data,
            "selected_match": selected,
            "position_match_summaries": match_summaries,
            "chrom": chrom,
            "start": start,
            "end": end,
            "position": parsed["position"],
            "track_data": track_res.data,
            "transcript_summaries": transcript_summaries,
            "transcript_count": len(transcript_summaries),
            "browser_url": browser_url(chrom, start, end, genome=genome),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "UCSC_API_BASE",
    "UCSC_BROWSER_BASE",
    "DEFAULT_GENOME",
    "DEFAULT_TRACK",
    "DEFAULT_REGION_SREBF2",
    "DEFAULT_CANONICAL_TRANSCRIPT_SREBF2",
    "parse_position",
    "browser_url",
    "summarize_position_match",
    "summarize_known_gene",
    "prefer_position_match",
    "search",
    "get_track_data",
    "fetch_gene_region",
]
