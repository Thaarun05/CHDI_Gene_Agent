"""NCATS Inxight Drugs client (substance search + relationships).

Does **not** normalize into evidence records — that belongs in Section 7a
orchestration / ``normalize/``.

Key endpoints (validated)::

    GET https://drugs.ncats.io/api/v1/substances/search
        ?facet=Primary+Target/{label}&top=...&skip=...
    GET https://drugs.ncats.io/api/v1/substances({uuid})/relationships

Facet matches alone are not confirmed target relationships.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "NCATS Inxight"
BASE = "https://drugs.ncats.io/api/v1"

DEFAULT_TOP = 10
PRIMARY_TARGET_FACET_PREFIX = "Primary Target/"


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


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any] | list[tuple[str, Any]] | None,
    settings: Settings,
) -> ToolResult:
    url = f"{BASE}/{path.lstrip('/')}"
    if params is None:
        query_pairs: list[tuple[str, str]] = []
        request_params: dict[str, Any] = {}
    elif isinstance(params, list):
        query_pairs = [(str(k), str(v)) for k, v in params if v is not None]
        request_params = {}
        for k, v in query_pairs:
            if k in request_params:
                existing = request_params[k]
                if isinstance(existing, list):
                    existing.append(v)
                else:
                    request_params[k] = [existing, v]
            else:
                request_params[k] = v
    else:
        query_pairs = [(str(k), str(v)) for k, v in params.items() if v is not None]
        request_params = {k: v for k, v in query_pairs}

    request_url = f"{url}?{urlencode(query_pairs)}" if query_pairs else url
    headers = {"Accept": "application/json"}
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query_pairs, headers=headers)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=request_params,
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
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def primary_target_facet(label: str) -> str:
    """Build an Inxight ``Primary Target/{label}`` facet value."""
    cleaned = str(label).strip()
    if cleaned.lower().startswith("primary target/"):
        return cleaned
    return f"{PRIMARY_TARGET_FACET_PREFIX}{cleaned}"


def search_substances(
    *,
    gene_symbol: str = "",
    query: str | None = None,
    primary_target: str | None = None,
    facets: list[str] | None = None,
    top: int = DEFAULT_TOP,
    skip: int = 0,
    extra_params: dict[str, Any] | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Search Inxight substances (supports Primary Target facet filters).

    Prefer ``primary_target`` (protein name / symbol) which becomes
    ``facet=Primary Target/{value}``. Additional ``facets`` and free-text
    ``query`` (``q``) may be supplied.
    """
    cfg = settings or get_settings()
    pairs: list[tuple[str, Any]] = []
    if query:
        pairs.append(("q", str(query).strip()))
    facet_values: list[str] = []
    if primary_target:
        facet_values.append(primary_target_facet(primary_target))
    if facets:
        for facet in facets:
            cleaned = str(facet).strip()
            if cleaned:
                facet_values.append(cleaned)
    for facet in facet_values:
        pairs.append(("facet", facet))
    pairs.append(("top", int(top)))
    pairs.append(("skip", int(skip)))
    if extra_params:
        for key, value in extra_params.items():
            if value is None:
                continue
            pairs.append((str(key), value))

    return _request_json(
        endpoint_name="search_substances",
        gene_symbol=gene_symbol or str(primary_target or query or ""),
        path="substances/search",
        params=pairs,
        settings=cfg,
    )


def substance_relationships(
    uuid: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch GINAS relationships for a substance UUID."""
    cfg = settings or get_settings()
    sid = str(uuid).strip()
    if not sid:
        return _tool_result(
            endpoint_name="substance_relationships",
            gene_symbol=gene_symbol,
            request_url=f"{BASE}/substances()/relationships",
            request_params={"uuid": ""},
            success=False,
            error_type="invalid_request",
            error_message="NCATS Inxight relationships require a substance uuid",
        )
    # GINAS-style path: /substances({uuid})/relationships
    path = f"substances({sid})/relationships"
    return _request_json(
        endpoint_name="substance_relationships",
        gene_symbol=gene_symbol or sid,
        path=path,
        params=None,
        settings=cfg,
    )


__all__ = [
    "SOURCE_NAME",
    "BASE",
    "DEFAULT_TOP",
    "PRIMARY_TARGET_FACET_PREFIX",
    "primary_target_facet",
    "search_substances",
    "substance_relationships",
]
