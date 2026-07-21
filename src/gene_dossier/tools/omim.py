"""OMIM client (entry search + entry detail).

Fetches OMIM gene/entry metadata. Priority C scaffold: requires ``OMIM_API_KEY``.
Does **not** normalize into evidence records — that belongs in
``normalize/variants.py``.

Key endpoints (validated)::

    GET https://api.omim.org/api/entry/search
        ?search={symbol}&include=geneMap&format=json&apiKey={{omim_api_key}}
    GET https://api.omim.org/api/entry
        ?mimNumber={mim}&include=geneMap,clinicalSynopsis,text&format=json
        &apiKey={{omim_api_key}}

For SREBF2, validated MIM number is ``600481``. NOTE: no strong OMIM
disease/phenotype relationship was confirmed from the validated response — do
not invent phenotype links.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "OMIM"
OMIM_API_BASE = "https://api.omim.org/api"

DEFAULT_SEARCH_INCLUDE = "geneMap"
DEFAULT_ENTRY_INCLUDE = "geneMap,clinicalSynopsis,text"
DEFAULT_MIM_SREBF2 = "600481"


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


def _safe_params(params: dict[str, Any]) -> dict[str, Any]:
    """Redact apiKey for provenance logging."""
    out = {k: v for k, v in params.items() if k.lower() != "apikey"}
    if any(k.lower() == "apikey" for k in params):
        out["apiKey"] = "***"
    return out


def summarize_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Extract key OMIM entry fields (not evidence)."""
    titles = entry.get("titles") or {}
    if not isinstance(titles, dict):
        titles = {}
    gene_map = entry.get("geneMap") or {}
    if not isinstance(gene_map, dict):
        gene_map = {}
    phenotype_maps = gene_map.get("phenotypeMapList") or []
    if not isinstance(phenotype_maps, list):
        phenotype_maps = []
    return {
        "mim_number": entry.get("mimNumber"),
        "preferred_title": titles.get("preferredTitle"),
        "alternative_titles": titles.get("alternativeTitles"),
        "chromosome": gene_map.get("chromosome"),
        "cyto_location": gene_map.get("cytoLocation"),
        "computed_cyto_location": gene_map.get("computedCytoLocation"),
        "gene_symbols": gene_map.get("geneSymbols"),
        "gene_name": gene_map.get("geneName"),
        "gene_ids": gene_map.get("geneIDs"),
        "ensembl_ids": gene_map.get("ensemblIDs"),
        "phenotype_map_list": phenotype_maps,
        "phenotype_map_count": len(phenotype_maps),
    }


def extract_entry_list(payload: Any) -> list[dict[str, Any]]:
    """Return entry dicts from an OMIM JSON payload."""
    if not isinstance(payload, dict):
        return []
    omim = payload.get("omim") or payload
    if not isinstance(omim, dict):
        return []
    entry_list = omim.get("entryList") or []
    if not isinstance(entry_list, list):
        return []
    out: list[dict[str, Any]] = []
    for item in entry_list:
        if isinstance(item, dict):
            entry = item.get("entry")
            if isinstance(entry, dict):
                out.append(entry)
            elif "mimNumber" in item:
                out.append(item)
    return out


def extract_mim_numbers(payload: Any) -> list[str]:
    """Unique MIM numbers from an OMIM search/entry payload."""
    out: list[str] = []
    seen: set[str] = set()
    for entry in extract_entry_list(payload):
        mim = entry.get("mimNumber")
        if mim is None:
            continue
        key = str(mim)
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def prefer_mim_number(
    entries: list[dict[str, Any]],
    gene_symbol: str,
) -> str | None:
    """Prefer a MIM whose geneMap symbols match ``gene_symbol``; else ``None``."""
    target = gene_symbol.strip().upper()
    if not target:
        return None
    matches: list[str] = []
    for entry in entries:
        gene_map = entry.get("geneMap") or {}
        if not isinstance(gene_map, dict):
            continue
        symbols_raw = gene_map.get("geneSymbols") or ""
        symbols = {s.strip().upper() for s in str(symbols_raw).replace(",", " ").split() if s.strip()}
        if target in symbols:
            mim = entry.get("mimNumber")
            if mim is not None:
                matches.append(str(mim))
    if len(matches) == 1:
        return matches[0]
    return None


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET an OMIM API path; return :class:`ToolResult` (never raises)."""
    if not settings.has_key("OMIM_API_KEY"):
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=f"{OMIM_API_BASE}/{path.lstrip('/')}",
            request_params={k: v for k, v in params.items() if k.lower() != "apikey"},
            success=False,
            error_type="requires_key",
            error_message="OMIM_API_KEY missing",
        )

    query = dict(params)
    query["apiKey"] = str(settings.omim_api_key)
    query.setdefault("format", "json")
    safe = _safe_params(query)
    url = f"{OMIM_API_BASE}/{path.lstrip('/')}"
    request_url = f"{url}?{urlencode(safe)}"
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
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
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
            request_params=safe,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def search_entries(
    gene_symbol: str,
    *,
    include: str = DEFAULT_SEARCH_INCLUDE,
    settings: Settings | None = None,
) -> ToolResult:
    """Search OMIM entries for ``gene_symbol``."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol:
        return _tool_result(
            endpoint_name="search_entries",
            gene_symbol=gene_symbol,
            request_url=f"{OMIM_API_BASE}/entry/search",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="gene_symbol is required",
        )
    params = {
        "search": symbol,
        "include": include,
        "format": "json",
    }
    return _request_json(
        endpoint_name="search_entries",
        gene_symbol=symbol,
        path="entry/search",
        params=params,
        settings=cfg,
    )


def get_entry(
    mim_number: str | int,
    *,
    gene_symbol: str = "",
    include: str = DEFAULT_ENTRY_INCLUDE,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch an OMIM entry by MIM number."""
    cfg = settings or get_settings()
    mim = str(mim_number).strip()
    if not mim:
        return _tool_result(
            endpoint_name="get_entry",
            gene_symbol=gene_symbol,
            request_url=f"{OMIM_API_BASE}/entry",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="mimNumber is required",
        )
    params = {
        "mimNumber": mim,
        "include": include,
        "format": "json",
    }
    return _request_json(
        endpoint_name="get_entry",
        gene_symbol=gene_symbol or mim,
        path="entry",
        params=params,
        settings=cfg,
    )


def fetch_gene_entry(
    gene_symbol: str,
    *,
    mim_number: str | int | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Search OMIM (if needed), select a MIM, then fetch entry detail.

    If ``mim_number`` is provided, skips search. Otherwise prefers a geneMap
    symbol match; if ambiguous/unmatched, returns search results without
    guessing a MIM.

    Includes caveat that absence of phenotype maps must not be overinterpreted.
    """
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    search_payload: Any = None
    selected_mim: str | None = str(mim_number).strip() if mim_number is not None else None
    selection_method = "provided" if selected_mim else None
    last_url = f"{OMIM_API_BASE}/entry/search"
    last_params: dict[str, Any] = {}
    last_status: int | None = None

    if not selected_mim:
        searched = search_entries(symbol, settings=cfg)
        last_url = searched.request_url
        last_params = searched.request_params
        last_status = searched.status_code
        if not searched.success:
            return _tool_result(
                endpoint_name="fetch_gene_entry",
                gene_symbol=symbol,
                request_url=searched.request_url,
                request_params=searched.request_params,
                success=False,
                status_code=searched.status_code,
                data={"search": searched.data},
                error_type=searched.error_type or "search_failed",
                error_message=searched.error_message or "OMIM search failed",
            )
        search_payload = searched.data
        entries = extract_entry_list(searched.data)
        selected_mim = prefer_mim_number(entries, symbol)
        if selected_mim is None:
            return _tool_result(
                endpoint_name="fetch_gene_entry",
                gene_symbol=symbol,
                request_url=searched.request_url,
                request_params=searched.request_params,
                success=True,
                status_code=searched.status_code,
                data={
                    "gene_symbol": symbol,
                    "selection_method": "ambiguous",
                    "selected_mim": None,
                    "search": searched.data,
                    "search_summaries": [summarize_entry(e) for e in entries],
                    "entry": None,
                    "entry_summary": None,
                    "caveat": (
                        "No unique geneMap symbol match; not guessing a MIM. "
                        "Do not invent OMIM disease/phenotype relationships."
                    ),
                },
            )
        selection_method = "matched"

    entry_res = get_entry(selected_mim, gene_symbol=symbol, settings=cfg)
    last_url = entry_res.request_url
    last_params = entry_res.request_params
    last_status = entry_res.status_code
    if not entry_res.success:
        return _tool_result(
            endpoint_name="fetch_gene_entry",
            gene_symbol=symbol,
            request_url=entry_res.request_url,
            request_params=entry_res.request_params,
            success=False,
            status_code=entry_res.status_code,
            data={
                "gene_symbol": symbol,
                "selected_mim": selected_mim,
                "selection_method": selection_method,
                "search": search_payload,
                "entry": entry_res.data,
            },
            error_type=entry_res.error_type or "entry_failed",
            error_message=entry_res.error_message or "OMIM entry fetch failed",
        )

    entries = extract_entry_list(entry_res.data)
    entry_summary = summarize_entry(entries[0]) if entries else None
    return _tool_result(
        endpoint_name="fetch_gene_entry",
        gene_symbol=symbol,
        request_url=last_url,
        request_params={
            "gene_symbol": symbol,
            "selected_mim": selected_mim,
            "selection_method": selection_method,
            **(last_params or {}),
        },
        success=True,
        status_code=last_status,
        data={
            "gene_symbol": symbol,
            "selected_mim": selected_mim,
            "selection_method": selection_method,
            "search": search_payload,
            "entry": entry_res.data,
            "entry_summary": entry_summary,
            "caveat": (
                "OMIM gene entries may lack disease/phenotype maps; do not invent "
                "phenotype relationships from absence of data."
            ),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "OMIM_API_BASE",
    "DEFAULT_SEARCH_INCLUDE",
    "DEFAULT_ENTRY_INCLUDE",
    "DEFAULT_MIM_SREBF2",
    "summarize_entry",
    "extract_entry_list",
    "extract_mim_numbers",
    "prefer_mim_number",
    "search_entries",
    "get_entry",
    "fetch_gene_entry",
]
