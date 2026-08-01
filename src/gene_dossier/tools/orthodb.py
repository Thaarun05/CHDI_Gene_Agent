"""OrthoDB v12 client (Section 1e supporting source).

Endpoints::

    GET https://data.orthodb.org/v12/orthodb_release_id
    GET https://data.orthodb.org/v12/genesearch?ncbi={entrez_gene_id}

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "OrthoDB"
ORTHODB_BASE = "https://data.orthodb.org/v12"


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


def public_gene_url(entrez_gene_id: str | int) -> str:
    return f"https://www.orthodb.org/?ncbi={str(entrez_gene_id).strip()}"


def _request(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any] | None,
    settings: Settings,
) -> ToolResult:
    url = f"{ORTHODB_BASE}/{path.lstrip('/')}"
    query = {k: str(v) for k, v in (params or {}).items() if v is not None}
    request_url = f"{url}?{urlencode(query)}" if query else url
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query or None)
        text = response.text
        try:
            payload: Any = response.json()
        except ValueError:
            payload = text.strip().strip('"')
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
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_release_id(
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch OrthoDB data/API release identifier."""
    cfg = settings or get_settings()
    return _request(
        endpoint_name="orthodb_release_id",
        gene_symbol=gene_symbol,
        path="orthodb_release_id",
        params=None,
        settings=cfg,
    )


def fetch_gene_search(
    entrez_gene_id: str | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Search OrthoDB by NCBI Gene ID (``ncbi`` query parameter)."""
    cfg = settings or get_settings()
    gid = str(entrez_gene_id).strip()
    if not gid:
        return _tool_result(
            endpoint_name="genesearch",
            gene_symbol=gene_symbol,
            request_url=f"{ORTHODB_BASE}/genesearch",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="entrez_gene_id is required",
        )
    return _request(
        endpoint_name="genesearch",
        gene_symbol=gene_symbol or gid,
        path="genesearch",
        params={"ncbi": gid},
        settings=cfg,
    )


def _organism_tax_id(payload: dict[str, Any]) -> int | None:
    organism = payload.get("organism")
    if not isinstance(organism, dict):
        return None
    org_id = str(organism.get("id") or "")
    if "_" in org_id:
        head = org_id.split("_", 1)[0]
        try:
            return int(head)
        except ValueError:
            pass
    xref = str(payload.get("organism_xref") or "")
    if "/taxonomy/" in xref:
        tail = xref.rsplit("/taxonomy/", 1)[-1]
        try:
            return int(tail.strip().split("/")[0])
        except ValueError:
            return None
    return None


def _find_entrez_xrefs(payload: Any, *, target: str) -> list[str]:
    """Collect string values that look like explicit Entrez/NCBI Gene IDs."""
    found: list[str] = []

    def walk(obj: Any, path: str = "") -> None:
        if isinstance(obj, dict):
            for key, value in obj.items():
                lower = str(key).lower()
                if lower in {
                    "ncbi_gene_id",
                    "entrez_gene_id",
                    "entrezid",
                    "ncbigeneid",
                    "geneid",
                }:
                    text = str(value).strip()
                    if text:
                        found.append(text)
                walk(value, f"{path}.{key}")
        elif isinstance(obj, list):
            for index, value in enumerate(obj[:50]):
                walk(value, f"{path}[{index}]")

    walk(payload)
    return [item for item in found if item == target]


def validate_gene_search_for_human(
    payload: Any,
    *,
    entrez_gene_id: str,
) -> tuple[bool, dict[str, Any]]:
    """Fail-closed validation for OrthoDB supporting evidence (not section abort)."""
    diagnostics: dict[str, Any] = {
        "requested_entrez_gene_id": str(entrez_gene_id).strip(),
    }
    if not isinstance(payload, dict):
        diagnostics["reason"] = "payload_not_object"
        return False, diagnostics
    status = str(payload.get("status") or "").lower()
    diagnostics["status"] = status
    if status != "ok":
        diagnostics["reason"] = "status_not_ok"
        return False, diagnostics
    tax_id = _organism_tax_id(payload)
    diagnostics["organism_tax_id"] = tax_id
    organism = payload.get("organism") if isinstance(payload.get("organism"), dict) else {}
    diagnostics["organism_name"] = organism.get("name")
    if tax_id != 9606 and str(organism.get("name") or "").lower() != "homo sapiens":
        diagnostics["reason"] = "organism_not_human"
        return False, diagnostics
    try:
        matched = int(str(payload.get("nb_genes_matched_the_query") or "0"))
    except ValueError:
        matched = None
    diagnostics["nb_genes_matched_the_query"] = matched
    if matched is None:
        diagnostics["reason"] = "match_count_unparseable"
        return False, diagnostics
    if matched != 1:
        diagnostics["reason"] = "ambiguous_or_empty_match"
        return False, diagnostics
    explicit = _find_entrez_xrefs(payload, target=str(entrez_gene_id).strip())
    diagnostics["explicit_entrez_xrefs"] = explicit
    # When OrthoDB omits an explicit Entrez xref field, accept a unique human hit
    # for a request that was issued with ncbi={entrez} (recorded by the caller).
    diagnostics["validation_mode"] = (
        "explicit_xref" if explicit else "unique_human_ncbi_query"
    )
    diagnostics["reason"] = "ok"
    return True, diagnostics


__all__ = [
    "SOURCE_NAME",
    "ORTHODB_BASE",
    "public_gene_url",
    "fetch_release_id",
    "fetch_gene_search",
    "validate_gene_search_for_human",
]
