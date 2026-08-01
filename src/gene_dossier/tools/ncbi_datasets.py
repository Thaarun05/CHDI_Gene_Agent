"""NCBI Datasets API client (gene / ortholog reports).

Fetches ortholog reports for an NCBI Gene ID (and optional gene package metadata).
Does **not** normalize into evidence records — that belongs in
``normalize/gene_identity.py``.

Key endpoint (validated)::

    GET https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}/orthologs
        ?returned_content=COMPLETE&page_size=1000

NOTE: Use ``page_size=1000`` so ortholog pages are not truncated.

For SREBF2, the expected Entrez Gene ID is ``6721``.

Optional ``NCBI_API_KEY`` is sent as an ``api-key`` header when present (redacted
in provenance fields).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import quote, urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "NCBI Datasets"
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"

DEFAULT_PAGE_SIZE = 1000
DEFAULT_RETURNED_CONTENT = "COMPLETE"
DEFAULT_GENE_ID_SREBF2 = "6721"


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


def _headers(settings: Settings) -> dict[str, str]:
    """Build request headers; include api-key when configured."""
    headers = {"Accept": "application/json"}
    if settings.has_key("NCBI_API_KEY"):
        headers["api-key"] = str(settings.ncbi_api_key)
    return headers


def _safe_params(params: dict[str, Any], *, api_key_used: bool) -> dict[str, Any]:
    """Copy params for provenance; redact api-key presence."""
    out = dict(params)
    if api_key_used:
        out["api-key"] = "***"
    return out


def _report_query_ids(report: dict[str, Any], payload: dict[str, Any] | None = None) -> list[str]:
    """Return query Gene IDs attached to a report or page payload."""
    raw = report.get("query")
    if raw is None and isinstance(payload, dict):
        raw = payload.get("query")
    if raw is None:
        return []
    if isinstance(raw, list):
        return [str(item).strip() for item in raw if str(item).strip()]
    text = str(raw).strip()
    return [text] if text else []


def summarize_ortholog_report(
    report: dict[str, Any],
    *,
    payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Extract key fields from one ortholog ``reports[*]`` entry (not evidence)."""
    gene = report.get("gene") or {}
    if not isinstance(gene, dict):
        gene = {}
    organism = gene.get("organism") or {}
    if not isinstance(organism, dict):
        organism = {}
    annotations = gene.get("annotations")
    genomic_locations: list[Any] = []
    if isinstance(annotations, list):
        for item in annotations:
            if not isinstance(item, dict):
                continue
            locs = item.get("genomic_locations")
            if isinstance(locs, list):
                genomic_locations.extend(locs)
    return {
        "gene_id": gene.get("gene_id"),
        "symbol": gene.get("symbol"),
        "description": gene.get("description"),
        "tax_id": gene.get("tax_id"),
        "taxname": gene.get("taxname"),
        "common_name": gene.get("common_name"),
        "scientific_name": (
            gene.get("taxname")
            or organism.get("scientific_name")
            or gene.get("scientific_name")
        ),
        "type": gene.get("type"),
        "chromosomes": gene.get("chromosomes"),
        "swiss_prot_accessions": gene.get("swiss_prot_accessions"),
        "ensembl_gene_ids": gene.get("ensembl_gene_ids"),
        "protein_count": gene.get("protein_count"),
        "gene_groups": gene.get("gene_groups"),
        "annotations": annotations,
        "genomic_locations": genomic_locations or None,
        "genomic_ranges": gene.get("genomic_ranges"),
        "query_gene_ids": _report_query_ids(report, payload),
    }


def extract_next_page_token(payload: Any) -> str | None:
    """Return the next-page token from a Datasets ortholog payload, if any."""
    if not isinstance(payload, dict):
        return None
    for key in ("next_page_token", "nextPageToken", "page_token", "pageToken"):
        raw = payload.get(key)
        if raw is None:
            continue
        text = str(raw).strip()
        if text:
            return text
    return None


def extract_reports(payload: Any) -> list[dict[str, Any]]:
    """Return ``reports`` list from a Datasets gene/ortholog payload."""
    if not isinstance(payload, dict):
        return []
    reports = payload.get("reports")
    if isinstance(reports, list):
        return [r for r in reports if isinstance(r, dict)]
    return []


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a Datasets v2 path and return :class:`ToolResult`."""
    url = f"{DATASETS_BASE}/{path.lstrip('/')}"
    query = {k: str(v) for k, v in params.items() if v is not None}
    headers = _headers(settings)
    api_key_used = "api-key" in headers
    safe = _safe_params(query, api_key_used=api_key_used)
    request_url = f"{url}?{urlencode(query)}" if query else url
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query, headers=headers)
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


def gene_by_id(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    returned_content: str = DEFAULT_RETURNED_CONTENT,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch NCBI Datasets gene package/report for an Entrez Gene ID."""
    cfg = settings or get_settings()
    gid = str(gene_id).strip()
    if not gid:
        return _tool_result(
            endpoint_name="gene_by_id",
            gene_symbol=gene_symbol,
            request_url=f"{DATASETS_BASE}/gene/id/",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="gene_id is required",
        )
    path = f"gene/id/{quote(gid, safe='')}"
    params = {"returned_content": returned_content}
    return _request_json(
        endpoint_name="gene_by_id",
        gene_symbol=gene_symbol or gid,
        path=path,
        params=params,
        settings=cfg,
    )


DEFAULT_MAX_ORTHOLOG_PAGES = 100


def orthologs_by_gene_id(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    returned_content: str = DEFAULT_RETURNED_CONTENT,
    page_size: int = DEFAULT_PAGE_SIZE,
    page_token: str | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch one NCBI Datasets ortholog page for an Entrez Gene ID."""
    cfg = settings or get_settings()
    gid = str(gene_id).strip()
    if not gid:
        return _tool_result(
            endpoint_name="orthologs_by_gene_id",
            gene_symbol=gene_symbol,
            request_url=f"{DATASETS_BASE}/gene/id/",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="gene_id is required",
        )
    path = f"gene/id/{quote(gid, safe='')}/orthologs"
    params: dict[str, Any] = {
        "returned_content": returned_content,
        "page_size": page_size,
    }
    if page_token:
        params["page_token"] = page_token
    return _request_json(
        endpoint_name="orthologs_by_gene_id",
        gene_symbol=gene_symbol or gid,
        path=path,
        params=params,
        settings=cfg,
    )


def iter_ortholog_pages(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    returned_content: str = DEFAULT_RETURNED_CONTENT,
    page_size: int = DEFAULT_PAGE_SIZE,
    max_pages: int = DEFAULT_MAX_ORTHOLOG_PAGES,
    settings: Settings | None = None,
) -> tuple[list[ToolResult], dict[str, Any]]:
    """Fetch ortholog pages until exhaustion or a defensive stop.

    Returns ``(page_tool_results, pagination_audit)``. Each successful or failed
    HTTP attempt is a separate :class:`ToolResult` for one ApiRun/raw artifact.
    """
    cfg = settings or get_settings()
    gid = str(gene_id).strip()
    pages: list[ToolResult] = []
    seen_tokens: set[str] = set()
    token: str | None = None
    retrieval_complete = False
    stop_reason: str | None = None
    empty_token_streak = 0

    for page_index in range(1, max_pages + 1):
        if token is not None:
            if token in seen_tokens:
                stop_reason = "repeated_page_token"
                break
            seen_tokens.add(token)
        result = orthologs_by_gene_id(
            gid,
            gene_symbol=gene_symbol,
            returned_content=returned_content,
            page_size=page_size,
            page_token=token,
            settings=cfg,
        )
        pages.append(result)
        if not result.success:
            stop_reason = result.error_type or "page_request_failed"
            break
        payload = result.data if isinstance(result.data, dict) else {}
        reports = extract_reports(payload)
        next_token = extract_next_page_token(payload)
        if next_token and not reports:
            empty_token_streak += 1
            if empty_token_streak >= 2:
                stop_reason = "empty_pages_with_token"
                break
        else:
            empty_token_streak = 0
        # Drift check: page-level or report-level query should match requested gene.
        page_query = _report_query_ids({}, payload)
        if not page_query and reports:
            page_query = _report_query_ids(reports[0], payload)
        if page_query and gid not in page_query:
            stop_reason = "query_gene_id_drift"
            break
        if not next_token:
            retrieval_complete = True
            stop_reason = "exhausted"
            break
        token = next_token
    else:
        stop_reason = "max_pages_exceeded"

    audit = {
        "requested_gene_id": gid,
        "page_count": len(pages),
        "successful_pages": sum(1 for page in pages if page.success),
        "retrieval_complete": retrieval_complete,
        "stop_reason": stop_reason,
        "seen_page_tokens": sorted(seen_tokens),
        "max_pages": max_pages,
        "page_size": page_size,
    }
    return pages, audit


def fetch_orthologs(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    page_size: int = DEFAULT_PAGE_SIZE,
    include_gene_package: bool = False,
    max_pages: int = DEFAULT_MAX_ORTHOLOG_PAGES,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch orthologs (validated) with light summaries across pages.

    On success, ``data`` includes::

        {
          "gene_id": ...,
          "gene_symbol": ...,
          "orthologs": <first_or_combined_raw>,
          "ortholog_pages": [<raw page payloads>...],
          "ortholog_summaries": [...],
          "ortholog_count": N,
          "retrieval_complete": bool,
          "pagination": {...},
          "gene_package": <optional raw>,
        }

    Never raises. Prefer :func:`iter_ortholog_pages` when each page must become
    its own ApiRun/raw artifact (Section 1e).
    """
    cfg = settings or get_settings()
    gid = str(gene_id).strip()
    pages, pagination = iter_ortholog_pages(
        gid,
        gene_symbol=gene_symbol,
        page_size=page_size,
        max_pages=max_pages,
        settings=cfg,
    )
    if not pages:
        return _tool_result(
            endpoint_name="fetch_orthologs",
            gene_symbol=gene_symbol or gid,
            request_url=f"{DATASETS_BASE}/gene/id/{quote(gid, safe='')}/orthologs",
            request_params={"gene_id": gid, "page_size": page_size},
            success=False,
            data={"gene_id": gid, "pagination": pagination},
            error_type="orthologs_failed",
            error_message="NCBI Datasets orthologs returned no pages",
        )

    first = pages[0]
    if not first.success and len(pages) == 1:
        return _tool_result(
            endpoint_name="fetch_orthologs",
            gene_symbol=gene_symbol or gid,
            request_url=first.request_url,
            request_params=first.request_params,
            success=False,
            status_code=first.status_code,
            data={"gene_id": gid, "orthologs": first.data, "pagination": pagination},
            error_type=first.error_type or "orthologs_failed",
            error_message=first.error_message or "NCBI Datasets orthologs failed",
        )

    summaries: list[dict[str, Any]] = []
    raw_pages: list[Any] = []
    for page in pages:
        if not page.success or not isinstance(page.data, dict):
            continue
        raw_pages.append(page.data)
        for report in extract_reports(page.data):
            summaries.append(summarize_ortholog_report(report, payload=page.data))

    gene_payload: Any = None
    last_url = pages[-1].request_url
    last_params = pages[-1].request_params
    last_status = pages[-1].status_code
    if include_gene_package:
        gene = gene_by_id(gid, gene_symbol=gene_symbol, settings=cfg)
        last_url = gene.request_url
        last_params = gene.request_params
        last_status = gene.status_code
        if not gene.success:
            return _tool_result(
                endpoint_name="fetch_orthologs",
                gene_symbol=gene_symbol or gid,
                request_url=gene.request_url,
                request_params=gene.request_params,
                success=False,
                status_code=gene.status_code,
                data={
                    "gene_id": gid,
                    "orthologs": raw_pages[0] if raw_pages else first.data,
                    "ortholog_pages": raw_pages,
                    "ortholog_summaries": summaries,
                    "ortholog_count": len(summaries),
                    "retrieval_complete": bool(pagination.get("retrieval_complete")),
                    "pagination": pagination,
                    "gene_package": gene.data,
                },
                error_type=gene.error_type or "gene_package_failed",
                error_message=gene.error_message
                or "NCBI Datasets gene package failed",
            )
        gene_payload = gene.data

    success = bool(raw_pages)
    return _tool_result(
        endpoint_name="fetch_orthologs",
        gene_symbol=gene_symbol or gid,
        request_url=last_url,
        request_params={
            "gene_id": gid,
            "page_size": page_size,
            "include_gene_package": include_gene_package,
            **(last_params or {}),
        },
        success=success,
        status_code=last_status,
        data={
            "gene_id": gid,
            "gene_symbol": gene_symbol or None,
            "orthologs": raw_pages[0] if raw_pages else first.data,
            "ortholog_pages": raw_pages,
            "ortholog_summaries": summaries,
            "ortholog_count": len(summaries),
            "retrieval_complete": bool(pagination.get("retrieval_complete")),
            "pagination": pagination,
            "gene_package": gene_payload,
        },
        error_type=None if success else (first.error_type or "orthologs_failed"),
        error_message=None
        if success
        else (first.error_message or "NCBI Datasets orthologs failed"),
    )


__all__ = [
    "SOURCE_NAME",
    "DATASETS_BASE",
    "DEFAULT_PAGE_SIZE",
    "DEFAULT_RETURNED_CONTENT",
    "DEFAULT_GENE_ID_SREBF2",
    "DEFAULT_MAX_ORTHOLOG_PAGES",
    "summarize_ortholog_report",
    "extract_reports",
    "extract_next_page_token",
    "gene_by_id",
    "orthologs_by_gene_id",
    "iter_ortholog_pages",
    "fetch_orthologs",
]
