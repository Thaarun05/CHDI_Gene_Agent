"""Patents client via SerpAPI Google Patents.

Searches Google Patents for gene-related patents. Priority C scaffold: requires
``SERPAPI_API_KEY``. Does **not** normalize into evidence records — relevance
grading (``direct`` / ``pathway`` / ``marker-list`` / ``weak``) belongs in a
later normalizer; the client only provides provisional hints.

Key endpoint (validated)::

    GET https://serpapi.com/search.json
        ?engine=google_patents
        &q=(SREBF2 OR SREBP2 OR "SREBP-2")
        &api_key={{serpapi_api_key}}

Build link::

    https://patents.google.com/{patent_id}

NOTE: Some hits only mention the gene in large gene lists — do not overinterpret
those as gene-specific patents.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import quote, urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Patents"
SERPAPI_SEARCH_URL = "https://serpapi.com/search.json"
GOOGLE_PATENTS_BASE = "https://patents.google.com"

# Gene-specific patent query tokens (do not apply SREBF2 aliases globally).
GENE_SPECIFIC_PATENT_TERMS: dict[str, list[str]] = {
    "SREBF2": [
        "SREBF2",
        "SREBP2",
        '"SREBP-2"',
        '"sterol regulatory element-binding protein 2"',
    ],
}


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


def patent_url(patent_id: str) -> str:
    """Build a Google Patents URL for a patent ID."""
    # Keep path separators in IDs like ``patent/US1234567A1``.
    return f"{GOOGLE_PATENTS_BASE}/{quote(str(patent_id).strip(), safe='/')}"


def build_patent_query(
    gene_symbol: str,
    aliases: list[str] | None = None,
) -> str:
    """Build a Google Patents ``q`` string for ``gene_symbol``.

    Includes the symbol, optional workflow aliases, and gene-specific mapped
    terms only when present in :data:`GENE_SPECIFIC_PATENT_TERMS`.
    """
    symbol = gene_symbol.strip()
    tokens: list[str] = []
    seen: set[str] = set()

    def _add(token: str) -> None:
        cleaned = token.strip()
        if not cleaned:
            return
        key = cleaned.lower()
        if key in seen:
            return
        seen.add(key)
        tokens.append(cleaned)

    if symbol:
        _add(symbol)
        _add(f'"{symbol}"')
    if aliases:
        for alias in aliases:
            _add(str(alias))
            if " " in str(alias).strip() or "-" in str(alias):
                _add(f'"{str(alias).strip()}"')
    for mapped in GENE_SPECIFIC_PATENT_TERMS.get(symbol.upper(), []):
        _add(mapped)

    if not tokens:
        return ""
    if len(tokens) == 1:
        return tokens[0]
    return "(" + " OR ".join(tokens) + ")"


def summarize_patent(
    row: dict[str, Any],
    *,
    gene_symbol: str = "",
) -> dict[str, Any]:
    """Extract key SerpAPI Google Patents fields (not evidence)."""
    patent_id = row.get("patent_id")
    return {
        "title": row.get("title"),
        "publication_number": row.get("publication_number"),
        "patent_id": patent_id,
        "assignee": row.get("assignee"),
        "inventor": row.get("inventor"),
        "priority_date": row.get("priority_date"),
        "filing_date": row.get("filing_date"),
        "publication_date": row.get("publication_date"),
        "grant_date": row.get("grant_date"),
        "snippet": row.get("snippet"),
        "patent_link": row.get("patent_link"),
        "pdf": row.get("pdf"),
        "google_patents_url": patent_url(str(patent_id)) if patent_id else None,
        # Provisional triage only — final relevance belongs in the normalizer.
        "relevance_hint": suggest_relevance_hint(row, gene_symbol=gene_symbol),
    }


def suggest_relevance_hint(
    row: dict[str, Any],
    *,
    gene_symbol: str = "",
) -> str:
    """Provisional relevance hint for triage (not normalized evidence).

    Levels: ``direct``, ``pathway``, ``marker-list``, ``weak``.
    """
    symbol = gene_symbol.strip().upper()
    title = str(row.get("title") or "")
    snippet = str(row.get("snippet") or "")
    text = f"{title} {snippet}"
    upper = text.upper()

    if symbol and symbol in title.upper():
        return "direct"
    # Heuristic: dense comma-separated gene-like tokens → marker list mention.
    comma_parts = [p.strip() for p in snippet.replace(";", ",").split(",") if p.strip()]
    if len(comma_parts) >= 8 and symbol and symbol in upper:
        return "marker-list"
    if symbol and symbol in upper:
        return "pathway"
    return "weak"


def search_patents(
    gene_symbol: str,
    *,
    query: str | None = None,
    aliases: list[str] | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Run a SerpAPI Google Patents search for ``gene_symbol``."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol and not query:
        return _tool_result(
            endpoint_name="search_patents",
            gene_symbol=gene_symbol,
            request_url=SERPAPI_SEARCH_URL,
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="gene_symbol or query is required",
        )

    if not cfg.has_key("SERPAPI_API_KEY"):
        return _tool_result(
            endpoint_name="search_patents",
            gene_symbol=symbol or gene_symbol,
            request_url=SERPAPI_SEARCH_URL,
            request_params={"engine": "google_patents"},
            success=False,
            error_type="requires_key",
            error_message="SERPAPI_API_KEY missing",
        )

    q = query if query is not None else build_patent_query(symbol, aliases=aliases)
    if not q.strip():
        return _tool_result(
            endpoint_name="search_patents",
            gene_symbol=symbol,
            request_url=SERPAPI_SEARCH_URL,
            request_params={"engine": "google_patents"},
            success=False,
            error_type="invalid_request",
            error_message="Empty patent search query",
        )

    params = {
        "engine": "google_patents",
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
            # SerpAPI may return HTTP 200 with an error object.
            if isinstance(payload, dict) and payload.get("error"):
                return _tool_result(
                    endpoint_name="search_patents",
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
                endpoint_name="search_patents",
                gene_symbol=symbol or gene_symbol,
                request_url=request_url,
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="search_patents",
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
            endpoint_name="search_patents",
            gene_symbol=symbol or gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="search_patents",
            gene_symbol=symbol or gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="search_patents",
            gene_symbol=symbol or gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_patents(
    gene_symbol: str,
    *,
    aliases: list[str] | None = None,
    query: str | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Search patents and attach light summaries + provisional relevance hints.

    On success, ``data`` includes::

        {
          "gene_symbol": ...,
          "query": ...,
          "raw": <SerpAPI payload>,
          "patent_summaries": [...],
          "patent_count": N,
          "caveat": "...",
        }
    """
    cfg = settings or get_settings()
    q = query if query is not None else build_patent_query(gene_symbol, aliases=aliases)
    result = search_patents(
        gene_symbol, query=q, aliases=aliases, settings=cfg
    )
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_patents",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"query": q, "raw": result.data},
            error_type=result.error_type or "search_failed",
            error_message=result.error_message or "Patent search failed",
        )

    organic = []
    if isinstance(result.data, dict):
        organic = result.data.get("organic_results") or []
    if not isinstance(organic, list):
        organic = []
    summaries = [
        summarize_patent(row, gene_symbol=gene_symbol)
        for row in organic
        if isinstance(row, dict)
    ]
    return _tool_result(
        endpoint_name="fetch_patents",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "query": q,
            "raw": result.data,
            "patent_summaries": summaries,
            "patent_count": len(summaries),
            "caveat": (
                "Do not overinterpret patents that only mention the gene in "
                "large marker lists; relevance_hint is provisional."
            ),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "SERPAPI_SEARCH_URL",
    "GOOGLE_PATENTS_BASE",
    "GENE_SPECIFIC_PATENT_TERMS",
    "patent_url",
    "build_patent_query",
    "summarize_patent",
    "suggest_relevance_hint",
    "search_patents",
    "fetch_patents",
]
