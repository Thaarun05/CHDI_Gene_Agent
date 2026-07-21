"""Commercial antibody discovery client (SerpAPI Google search).

Priority C scaffold: there is no single reliable official REST API for commercial
antibodies. SerpAPI Google search is used for discovery only. Vendor product
pages are the final provenance; snippets can be stale. Does **not** normalize
into evidence records and does **not** invent catalog numbers or product claims.

Key endpoint (validated)::

    GET https://serpapi.com/search.json
        ?engine=google
        &q=SREBF2 OR SREBP2 antibody catalog Abcam Novus ...
        &api_key={{serpapi_api_key}}

Uses ``SERPAPI_API_KEY`` when present (redacted in provenance). Without the key,
calls fail with ``requires_key`` so coverage can mark the source appropriately.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import re
from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Antibodies"
SERPAPI_SEARCH_URL = "https://serpapi.com/search.json"

DEFAULT_VENDORS = (
    "Abcam",
    "Novus",
    "R&D",
    "R&D Systems",
    "Santa Cruz",
    "OriGene",
    "LSBio",
    "Biorbyt",
    "Sino Biological",
    "Proteintech",
)

# Gene-specific query aliases (do not apply SREBF2 terms globally).
GENE_SPECIFIC_ANTIBODY_TERMS: dict[str, list[str]] = {
    "SREBF2": ["SREBF2", "SREBP2", '"SREBP2 Antibody"'],
}

_CATALOG_RE = re.compile(
    r"(?:catalog\s*(?:number|no\.?|#)|cat\.?\s*(?:no\.?|#)|sku)\s*[:#-]?\s*"
    r"([A-Za-z0-9][A-Za-z0-9._/-]{2,})",
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


def _safe_params(params: dict[str, Any]) -> dict[str, Any]:
    """Redact api_key for provenance logging."""
    out = {k: v for k, v in params.items() if k != "api_key"}
    if "api_key" in params:
        out["api_key"] = "***"
    return out


def build_antibody_query(
    gene_symbol: str,
    *,
    aliases: list[str] | None = None,
    vendors: list[str] | tuple[str, ...] | None = None,
) -> str:
    """Build a Google discovery query for commercial antibodies."""
    symbol = gene_symbol.strip()
    terms: list[str] = []
    seen: set[str] = set()

    def _add(token: str) -> None:
        cleaned = token.strip()
        if not cleaned:
            return
        key = cleaned.lower()
        if key in seen:
            return
        seen.add(key)
        terms.append(cleaned)

    if symbol:
        _add(symbol)
    if aliases:
        for alias in aliases:
            _add(str(alias))
    for mapped in GENE_SPECIFIC_ANTIBODY_TERMS.get(symbol.upper(), []):
        _add(mapped)

    _add("antibody")
    _add("catalog")
    vendor_list = list(vendors) if vendors is not None else list(DEFAULT_VENDORS)
    for vendor in vendor_list:
        _add(vendor)

    return " ".join(terms)


def guess_vendor_name(title: str, source: str, link: str) -> str | None:
    """Best-effort vendor guess from title/source/link (discovery only)."""
    blob = f"{title} {source} {link}".lower()
    for vendor in DEFAULT_VENDORS:
        if vendor.lower() in blob:
            return vendor
    return source.strip() or None


def guess_catalog_number(title: str, snippet: str) -> str | None:
    """Extract an explicit catalog-like token when present; else ``None``.

    Does **not** invent catalog numbers.
    """
    for text in (title, snippet):
        match = _CATALOG_RE.search(text or "")
        if match:
            return match.group(1)
    return None


def summarize_antibody_hit(row: dict[str, Any]) -> dict[str, Any]:
    """Extract discovery fields from one Google organic result (not evidence)."""
    title = str(row.get("title") or "")
    source = str(row.get("source") or "")
    link = str(row.get("link") or "")
    snippet = str(row.get("snippet") or "")
    return {
        "title": title or None,
        "source": source or None,
        "link": link or None,
        "snippet": snippet or None,
        "antibody_product_name": title or None,
        "vendor_name": guess_vendor_name(title, source, link),
        "catalog_number": guess_catalog_number(title, snippet),
        "antibody_description": snippet or None,
        "product_url": link or None,
        "caveat": (
            "SerpAPI/Google snippets are discovery metadata and may be stale; "
            "vendor product pages are final provenance."
        ),
    }


def search_antibodies(
    gene_symbol: str,
    *,
    query: str | None = None,
    aliases: list[str] | None = None,
    vendors: list[str] | tuple[str, ...] | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Run a SerpAPI Google search for commercial antibody listings."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol and not query:
        return _tool_result(
            endpoint_name="search_antibodies",
            gene_symbol=gene_symbol,
            request_url=SERPAPI_SEARCH_URL,
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="gene_symbol or query is required",
        )

    if not cfg.has_key("SERPAPI_API_KEY"):
        return _tool_result(
            endpoint_name="search_antibodies",
            gene_symbol=symbol or gene_symbol,
            request_url=SERPAPI_SEARCH_URL,
            request_params={"engine": "google"},
            success=False,
            error_type="requires_key",
            error_message="SERPAPI_API_KEY missing",
        )

    q = (
        query
        if query is not None
        else build_antibody_query(symbol, aliases=aliases, vendors=vendors)
    )
    params = {
        "engine": "google",
        "q": q,
        "api_key": str(cfg.serpapi_api_key),
    }
    safe = _safe_params(params)
    request_url = f"{SERPAPI_SEARCH_URL}?{urlencode(safe)}"
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(SERPAPI_SEARCH_URL, params=params)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            if isinstance(payload, dict) and payload.get("error"):
                return _tool_result(
                    endpoint_name="search_antibodies",
                    gene_symbol=symbol or gene_symbol,
                    request_url=request_url,
                    request_params=safe,
                    success=False,
                    status_code=response.status_code,
                    data=payload,
                    error_type="api_error",
                    error_message=str(payload.get("error")),
                )
            return _tool_result(
                endpoint_name="search_antibodies",
                gene_symbol=symbol or gene_symbol,
                request_url=request_url,
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="search_antibodies",
            gene_symbol=symbol or gene_symbol,
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
            endpoint_name="search_antibodies",
            gene_symbol=symbol or gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="search_antibodies",
            gene_symbol=symbol or gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="search_antibodies",
            gene_symbol=symbol or gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_antibodies(
    gene_symbol: str,
    *,
    aliases: list[str] | None = None,
    vendors: list[str] | tuple[str, ...] | None = None,
    query: str | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Search for commercial antibodies and attach discovery summaries.

    On success, ``data`` includes raw SerpAPI payload and antibody summaries.
    Catalog numbers are only filled when explicitly present in text.
    """
    cfg = settings or get_settings()
    q = (
        query
        if query is not None
        else build_antibody_query(gene_symbol, aliases=aliases, vendors=vendors)
    )
    result = search_antibodies(
        gene_symbol,
        query=q,
        aliases=aliases,
        vendors=vendors,
        settings=cfg,
    )
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_antibodies",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"query": q, "raw": result.data},
            error_type=result.error_type or "search_failed",
            error_message=result.error_message or "Antibody search failed",
        )

    organic = []
    if isinstance(result.data, dict):
        organic = result.data.get("organic_results") or []
    if not isinstance(organic, list):
        organic = []
    summaries = [
        summarize_antibody_hit(row) for row in organic if isinstance(row, dict)
    ]
    return _tool_result(
        endpoint_name="fetch_antibodies",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "query": q,
            "raw": result.data,
            "antibody_summaries": summaries,
            "hit_count": len(summaries),
            "caveat": (
                "Discovery only. Do not invent product claims; confirm on vendor "
                "product pages before citing antibodies in the dossier."
            ),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "SERPAPI_SEARCH_URL",
    "DEFAULT_VENDORS",
    "GENE_SPECIFIC_ANTIBODY_TERMS",
    "build_antibody_query",
    "guess_vendor_name",
    "guess_catalog_number",
    "summarize_antibody_hit",
    "search_antibodies",
    "fetch_antibodies",
]
