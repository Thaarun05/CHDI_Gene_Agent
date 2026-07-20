"""PubMed client (E-utilities ESearch + ESummary + EFetch).

Searches literature for a gene in the context of Huntington's disease. Does **not**
normalize into evidence records — that belongs in ``normalize/literature.py``.

Default search term (per platform spec)::

    {gene_symbol}[Title/Abstract] AND "Huntington Disease"[MeSH Terms]

Rules:
- Do not assume every hit strongly supports a gene–HD link.
- Use ``NCBI_API_KEY`` when present.
- Never raise: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "PubMed"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

DEFAULT_RETMAX = 50
HD_MESH = '"Huntington Disease"[MeSH Terms]'


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


def _with_api_key(params: dict[str, Any], settings: Settings) -> dict[str, Any]:
    """Return a copy of ``params`` with ``api_key`` added when configured."""
    out = dict(params)
    if settings.has_key("NCBI_API_KEY"):
        out["api_key"] = settings.ncbi_api_key
    return out


def _safe_params(params: dict[str, Any]) -> dict[str, Any]:
    """Redact API key for provenance logging."""
    out = {k: v for k, v in params.items() if k != "api_key"}
    if "api_key" in params:
        out["api_key"] = "***"
    return out


def build_hd_search_term(gene_symbol: str) -> str:
    """Build the default gene × Huntington Disease PubMed search term."""
    return f'{gene_symbol}[Title/Abstract] AND {HD_MESH}'


def _request(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
    expect_json: bool = True,
) -> ToolResult:
    """GET an E-utilities endpoint; return :class:`ToolResult` (never raises)."""
    query = _with_api_key(params, settings)
    url = f"{EUTILS_BASE}/{path}"
    safe = _safe_params(query)
    request_url = f"{url}?{urlencode(query)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query)
        if expect_json:
            try:
                payload: Any = response.json()
            except ValueError:
                payload = {"raw_text": response.text[:4000]}
        else:
            payload = {"raw_text": response.text, "content_type": response.headers.get("content-type")}

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


def esearch(
    gene_symbol: str,
    *,
    term: str | None = None,
    retmax: int = DEFAULT_RETMAX,
    settings: Settings | None = None,
) -> ToolResult:
    """Search PubMed. Default term is gene Title/Abstract + HD MeSH."""
    cfg = settings or get_settings()
    search_term = term if term is not None else build_hd_search_term(gene_symbol)
    params = {
        "db": "pubmed",
        "term": search_term,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": "relevance",
    }
    return _request(
        endpoint_name="esearch",
        gene_symbol=gene_symbol,
        path="esearch.fcgi",
        params=params,
        settings=cfg,
        expect_json=True,
    )


def esummary(
    pubmed_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PubMed ESummary metadata for one or more PMIDs."""
    cfg = settings or get_settings()
    if isinstance(pubmed_ids, (list, tuple)):
        id_str = ",".join(str(p) for p in pubmed_ids)
    else:
        id_str = str(pubmed_ids)
    params = {
        "db": "pubmed",
        "id": id_str,
        "retmode": "json",
    }
    return _request(
        endpoint_name="esummary",
        gene_symbol=gene_symbol or id_str,
        path="esummary.fcgi",
        params=params,
        settings=cfg,
        expect_json=True,
    )


def efetch(
    pubmed_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    rettype: str = "abstract",
    retmode: str = "xml",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PubMed records (default: abstract XML) for one or more PMIDs."""
    cfg = settings or get_settings()
    if isinstance(pubmed_ids, (list, tuple)):
        id_str = ",".join(str(p) for p in pubmed_ids)
    else:
        id_str = str(pubmed_ids)
    params = {
        "db": "pubmed",
        "id": id_str,
        "rettype": rettype,
        "retmode": retmode,
    }
    return _request(
        endpoint_name="efetch",
        gene_symbol=gene_symbol or id_str,
        path="efetch.fcgi",
        params=params,
        settings=cfg,
        expect_json=False,
    )


def extract_id_list(esearch_result: ToolResult) -> list[str]:
    """Return PMIDs from a successful ESearch :class:`ToolResult`."""
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return []
    result = esearch_result.data.get("esearchresult") or {}
    ids = result.get("idlist") or []
    return [str(i) for i in ids]


def search_hd_literature(
    gene_symbol: str,
    *,
    retmax: int = DEFAULT_RETMAX,
    fetch_abstracts: bool = True,
    settings: Settings | None = None,
) -> ToolResult:
    """Run HD-focused ESearch, then ESummary (+ optional EFetch abstracts).

    On success, ``data`` is::

        {
          "gene_symbol": ...,
          "search_term": ...,
          "pmids": [...],
          "count": int | None,
          "esearch": <raw>,
          "esummary": <raw> | null,
          "efetch": <raw xml wrapper> | null,
          "caveat": "Do not assume every hit strongly supports the gene–HD link.",
        }

    Never raises.
    """
    cfg = settings or get_settings()
    term = build_hd_search_term(gene_symbol)
    search = esearch(gene_symbol, term=term, retmax=retmax, settings=cfg)
    if not search.success:
        return _tool_result(
            endpoint_name="search_hd_literature",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=False,
            status_code=search.status_code,
            data={"esearch": search.data, "search_term": term},
            error_type=search.error_type or "esearch_failed",
            error_message=search.error_message or "PubMed ESearch failed",
        )

    pmids = extract_id_list(search)
    count_raw = None
    if isinstance(search.data, dict):
        count_raw = (search.data.get("esearchresult") or {}).get("count")

    if not pmids:
        return _tool_result(
            endpoint_name="search_hd_literature",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=True,
            status_code=search.status_code,
            data={
                "gene_symbol": gene_symbol,
                "search_term": term,
                "pmids": [],
                "count": int(count_raw) if count_raw is not None else 0,
                "esearch": search.data,
                "esummary": None,
                "efetch": None,
                "caveat": "Do not assume every hit strongly supports the gene–HD link.",
            },
        )

    ids_for_fetch = pmids[:retmax]
    summary = esummary(ids_for_fetch, gene_symbol=gene_symbol, settings=cfg)
    fetch: ToolResult | None = None
    if fetch_abstracts:
        fetch = efetch(ids_for_fetch, gene_symbol=gene_symbol, settings=cfg)

    # Success if we got summaries; efetch failure is recorded but does not void the search.
    ok = summary.success
    error_type = None if ok else (summary.error_type or "esummary_failed")
    error_message = None if ok else (summary.error_message or "PubMed ESummary failed")

    return _tool_result(
        endpoint_name="search_hd_literature",
        gene_symbol=gene_symbol,
        request_url=summary.request_url if summary.success else search.request_url,
        request_params={
            "search_term": term,
            "pmids": ids_for_fetch,
            "fetch_abstracts": fetch_abstracts,
        },
        success=ok,
        status_code=summary.status_code or search.status_code,
        data={
            "gene_symbol": gene_symbol,
            "search_term": term,
            "pmids": ids_for_fetch,
            "count": int(count_raw) if count_raw is not None else len(pmids),
            "esearch": search.data,
            "esummary": summary.data,
            "efetch": fetch.data if fetch is not None else None,
            "efetch_success": fetch.success if fetch is not None else None,
            "caveat": "Do not assume every hit strongly supports the gene–HD link.",
        },
        error_type=error_type,
        error_message=error_message,
    )


__all__ = [
    "SOURCE_NAME",
    "EUTILS_BASE",
    "DEFAULT_RETMAX",
    "HD_MESH",
    "build_hd_search_term",
    "esearch",
    "esummary",
    "efetch",
    "extract_id_list",
    "search_hd_literature",
]
