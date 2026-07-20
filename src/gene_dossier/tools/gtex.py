"""GTEx Portal API v2 client.

Resolves a gene to a GTEx GENCODE ID, then fetches median tissue expression and
single-tissue eQTLs. Does **not** normalize into evidence records.

Key endpoints (validated)::

    GET /reference/gene?geneId={symbol}&genomeBuild=GRCh38/hg38
    GET /expression/medianGeneExpression?gencodeId={id}&datasetId=gtex_v8
    GET /association/singleTissueEqtl?gencodeId={id}&datasetId=gtex_v8&itemsPerPage=...

For SREBF2, the expected GTEx GENCODE ID is ``ENSG00000198911.11``.

GTEx is human-only. Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "GTEx"
GTEX_BASE = "https://gtexportal.org/api/v2"

DEFAULT_DATASET = "gtex_v8"
DEFAULT_GENOME_BUILD = "GRCh38/hg38"
DEFAULT_EQTL_PAGE_SIZE = 1000  # keep MVP payloads bounded; caller can raise


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


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a GTEx API path and return JSON as :class:`ToolResult`."""
    url = f"{GTEX_BASE}{path}"
    # httpx encodes list values as repeated keys (needed for tissueSiteDetailId filters).
    request_url = f"{url}?{urlencode(params, doseq=True)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=params)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:2000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=params,
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
            request_params=params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def _data_list(payload: Any) -> list[Any]:
    """Return GTEx ``data`` array from a response payload."""
    if isinstance(payload, dict):
        data = payload.get("data")
        if isinstance(data, list):
            return data
    return []


def prefer_gencode_id(gene_rows: list[Any], gene_symbol: str) -> str | None:
    """Pick a GENCODE ID from gene lookup rows, preferring exact geneSymbol match."""
    target = gene_symbol.strip().upper()
    exact: list[str] = []
    fallback: list[str] = []
    for row in gene_rows:
        if not isinstance(row, dict):
            continue
        gid = row.get("gencodeId")
        if not gid:
            continue
        symbol = str(row.get("geneSymbol") or row.get("geneSymbolUpper") or "").upper()
        if symbol == target:
            exact.append(str(gid))
        else:
            fallback.append(str(gid))
    if exact:
        return exact[0]
    return fallback[0] if fallback else None


def resolve_gene(
    gene_symbol: str,
    *,
    genome_build: str = DEFAULT_GENOME_BUILD,
    settings: Settings | None = None,
) -> ToolResult:
    """Resolve ``gene_symbol`` to a GTEx GENCODE ID via ``/reference/gene``."""
    cfg = settings or get_settings()
    params = {"geneId": gene_symbol, "genomeBuild": genome_build}
    result = _request_json(
        endpoint_name="resolve_gene",
        gene_symbol=gene_symbol,
        path="/reference/gene",
        params=params,
        settings=cfg,
    )
    if not result.success:
        return result

    rows = _data_list(result.data)
    gencode_id = prefer_gencode_id(rows, gene_symbol)
    if not gencode_id:
        return _tool_result(
            endpoint_name="resolve_gene",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data=result.data,
            error_type="no_results",
            error_message=f"No GTEx GENCODE ID for {gene_symbol}",
        )
    return _tool_result(
        endpoint_name="resolve_gene",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "gencode_id": gencode_id,
            "genome_build": genome_build,
            "raw": result.data,
        },
    )


def median_expression(
    gencode_id: str,
    *,
    gene_symbol: str = "",
    dataset_id: str = DEFAULT_DATASET,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch median gene expression across GTEx tissues."""
    # TODO: add sample-level /expression/geneExpression fallback for box/violin plots
    # (validated in CHDI API map §2a.3). Median endpoint remains the MVP default.
    cfg = settings or get_settings()
    params = {"gencodeId": gencode_id, "datasetId": dataset_id}
    return _request_json(
        endpoint_name="median_expression",
        gene_symbol=gene_symbol or gencode_id,
        path="/expression/medianGeneExpression",
        params=params,
        settings=cfg,
    )


def single_tissue_eqtl(
    gencode_id: str,
    *,
    gene_symbol: str = "",
    dataset_id: str = DEFAULT_DATASET,
    items_per_page: int = DEFAULT_EQTL_PAGE_SIZE,
    tissue_site_detail_ids: list[str] | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch single-tissue eQTL associations for a GENCODE ID."""
    # TODO: implement eQTL pagination (page / next pageToken). Validated map uses
    # itemsPerPage=100000 for full retrieval; MVP caps page size to bound payloads.
    cfg = settings or get_settings()
    params: dict[str, Any] = {
        "gencodeId": gencode_id,
        "datasetId": dataset_id,
        "itemsPerPage": str(items_per_page),
    }
    if tissue_site_detail_ids:
        params["tissueSiteDetailId"] = list(tissue_site_detail_ids)
    return _request_json(
        endpoint_name="single_tissue_eqtl",
        gene_symbol=gene_symbol or gencode_id,
        path="/association/singleTissueEqtl",
        params=params,
        settings=cfg,
    )


def fetch_expression_and_eqtl(
    gene_symbol: str,
    *,
    gencode_id: str | None = None,
    dataset_id: str = DEFAULT_DATASET,
    items_per_page: int = DEFAULT_EQTL_PAGE_SIZE,
    settings: Settings | None = None,
) -> ToolResult:
    """Resolve gene (if needed), then fetch median expression + eQTLs.

    On success, ``data`` includes ``gencode_id``, ``median_expression``, and ``eqtl``.
    """
    cfg = settings or get_settings()
    resolved_id = gencode_id
    resolve_payload: Any = None
    if not resolved_id:
        resolved = resolve_gene(gene_symbol, settings=cfg)
        if not resolved.success:
            return _tool_result(
                endpoint_name="fetch_expression_and_eqtl",
                gene_symbol=gene_symbol,
                request_url=resolved.request_url,
                request_params=resolved.request_params,
                success=False,
                status_code=resolved.status_code,
                data={"resolve_gene": resolved.data},
                error_type=resolved.error_type or "resolve_failed",
                error_message=resolved.error_message or "GTEx gene resolve failed",
            )
        resolve_payload = resolved.data
        resolved_id = (resolved.data or {}).get("gencode_id")
        if not resolved_id:
            return _tool_result(
                endpoint_name="fetch_expression_and_eqtl",
                gene_symbol=gene_symbol,
                request_url=resolved.request_url,
                request_params=resolved.request_params,
                success=False,
                status_code=resolved.status_code,
                data={"resolve_gene": resolved.data},
                error_type="no_results",
                error_message=f"No GTEx GENCODE ID for {gene_symbol}",
            )

    median = median_expression(
        resolved_id, gene_symbol=gene_symbol, dataset_id=dataset_id, settings=cfg
    )
    eqtl = single_tissue_eqtl(
        resolved_id,
        gene_symbol=gene_symbol,
        dataset_id=dataset_id,
        items_per_page=items_per_page,
        settings=cfg,
    )

    # Partial success is still useful: mark success if either call worked.
    ok = median.success or eqtl.success
    return _tool_result(
        endpoint_name="fetch_expression_and_eqtl",
        gene_symbol=gene_symbol,
        request_url=median.request_url if median.success else eqtl.request_url,
        request_params={
            "gencode_id": resolved_id,
            "dataset_id": dataset_id,
            "items_per_page": items_per_page,
        },
        success=ok,
        status_code=median.status_code or eqtl.status_code,
        data={
            "gene_symbol": gene_symbol,
            "gencode_id": resolved_id,
            "dataset_id": dataset_id,
            "resolve_gene": resolve_payload,
            "median_expression": median.data,
            "median_expression_success": median.success,
            "eqtl": eqtl.data,
            "eqtl_success": eqtl.success,
            "note": "GTEx is human-only.",
        },
        error_type=None if ok else (median.error_type or eqtl.error_type),
        error_message=None
        if ok
        else (median.error_message or eqtl.error_message or "GTEx expression/eQTL failed"),
    )


__all__ = [
    "SOURCE_NAME",
    "GTEX_BASE",
    "DEFAULT_DATASET",
    "DEFAULT_GENOME_BUILD",
    "DEFAULT_EQTL_PAGE_SIZE",
    "resolve_gene",
    "median_expression",
    "single_tissue_eqtl",
    "fetch_expression_and_eqtl",
    "prefer_gencode_id",
]
