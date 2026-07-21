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


def summarize_ortholog_report(report: dict[str, Any]) -> dict[str, Any]:
    """Extract key fields from one ortholog ``reports[*]`` entry (not evidence)."""
    gene = report.get("gene") or {}
    if not isinstance(gene, dict):
        gene = {}
    organism = gene.get("organism") or {}
    if not isinstance(organism, dict):
        organism = {}
    return {
        "gene_id": gene.get("gene_id"),
        "symbol": gene.get("symbol"),
        "description": gene.get("description"),
        "tax_id": gene.get("tax_id"),
        "common_name": gene.get("common_name"),
        "scientific_name": organism.get("scientific_name")
        or gene.get("scientific_name"),
        "chromosomes": gene.get("chromosomes"),
        "genomic_ranges": gene.get("genomic_ranges"),
        "annotations": gene.get("annotations"),
    }


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


def orthologs_by_gene_id(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    returned_content: str = DEFAULT_RETURNED_CONTENT,
    page_size: int = DEFAULT_PAGE_SIZE,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch NCBI Datasets ortholog reports for an Entrez Gene ID."""
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
    params = {
        "returned_content": returned_content,
        "page_size": page_size,
    }
    return _request_json(
        endpoint_name="orthologs_by_gene_id",
        gene_symbol=gene_symbol or gid,
        path=path,
        params=params,
        settings=cfg,
    )


def fetch_orthologs(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    page_size: int = DEFAULT_PAGE_SIZE,
    include_gene_package: bool = False,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch orthologs (validated) with light summaries.

    On success, ``data`` includes::

        {
          "gene_id": ...,
          "gene_symbol": ...,
          "orthologs": <raw>,
          "ortholog_summaries": [...],
          "ortholog_count": N,
          "gene_package": <optional raw>,
        }

    Never raises.
    """
    cfg = settings or get_settings()
    gid = str(gene_id).strip()
    orth = orthologs_by_gene_id(
        gid, gene_symbol=gene_symbol, page_size=page_size, settings=cfg
    )
    if not orth.success:
        return _tool_result(
            endpoint_name="fetch_orthologs",
            gene_symbol=gene_symbol or gid,
            request_url=orth.request_url,
            request_params=orth.request_params,
            success=False,
            status_code=orth.status_code,
            data={"gene_id": gid, "orthologs": orth.data},
            error_type=orth.error_type or "orthologs_failed",
            error_message=orth.error_message or "NCBI Datasets orthologs failed",
        )

    reports = extract_reports(orth.data)
    summaries = [summarize_ortholog_report(r) for r in reports]

    gene_payload: Any = None
    last_url = orth.request_url
    last_params = orth.request_params
    last_status = orth.status_code
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
                    "orthologs": orth.data,
                    "ortholog_summaries": summaries,
                    "ortholog_count": len(summaries),
                    "gene_package": gene.data,
                },
                error_type=gene.error_type or "gene_package_failed",
                error_message=gene.error_message
                or "NCBI Datasets gene package failed",
            )
        gene_payload = gene.data

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
        success=True,
        status_code=last_status,
        data={
            "gene_id": gid,
            "gene_symbol": gene_symbol or None,
            "orthologs": orth.data,
            "ortholog_summaries": summaries,
            "ortholog_count": len(summaries),
            "gene_package": gene_payload,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "DATASETS_BASE",
    "DEFAULT_PAGE_SIZE",
    "DEFAULT_RETURNED_CONTENT",
    "DEFAULT_GENE_ID_SREBF2",
    "summarize_ortholog_report",
    "extract_reports",
    "gene_by_id",
    "orthologs_by_gene_id",
    "fetch_orthologs",
]
