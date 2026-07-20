"""STRING database client.

Resolves a gene/protein to a STRING ID, then fetches interaction partners.
Does **not** normalize into evidence records — that belongs in ``normalize/ppi.py``.

Key endpoints (validated)::

    GET https://string-db.org/api/json/get_string_ids
        ?identifiers={symbol}&species=9606&echo_query=1&caller_identity=...
    GET https://string-db.org/api/json/interaction_partners
        ?identifiers={string_id}&species=9606&limit=100&required_score=400
        &network_type=functional&caller_identity=...

Uses ``settings.caller_identity`` (default ``gene_dossier_platform``).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "STRING"
STRING_BASE = "https://string-db.org/api/json"

SPECIES_HUMAN = 9606
DEFAULT_LIMIT = 100
DEFAULT_REQUIRED_SCORE = 400
DEFAULT_NETWORK_TYPE = "functional"


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
    """GET a STRING JSON API path and return :class:`ToolResult`."""
    url = f"{STRING_BASE}/{path.lstrip('/')}"
    request_url = f"{url}?{urlencode(params)}"
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


def prefer_string_id(rows: list[Any], gene_symbol: str) -> str | None:
    """Pick a STRING ID from ``get_string_ids`` rows, preferring exact preferredName."""
    target = gene_symbol.strip().upper()
    exact: list[str] = []
    fallback: list[str] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        sid = row.get("stringId") or row.get("string_id")
        if not sid:
            continue
        name = str(row.get("preferredName") or "").upper()
        if name == target:
            exact.append(str(sid))
        else:
            fallback.append(str(sid))
    if exact:
        return exact[0]
    return fallback[0] if fallback else None


def get_string_ids(
    gene_symbol: str,
    *,
    species: int = SPECIES_HUMAN,
    settings: Settings | None = None,
) -> ToolResult:
    """Map ``gene_symbol`` to STRING identifier(s)."""
    cfg = settings or get_settings()
    params = {
        "identifiers": gene_symbol,
        "species": str(species),
        "echo_query": "1",
        "caller_identity": cfg.caller_identity,
    }
    result = _request_json(
        endpoint_name="get_string_ids",
        gene_symbol=gene_symbol,
        path="get_string_ids",
        params=params,
        settings=cfg,
    )
    if not result.success:
        return result

    rows = result.data if isinstance(result.data, list) else []
    string_id = prefer_string_id(rows, gene_symbol)
    if not string_id:
        return _tool_result(
            endpoint_name="get_string_ids",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data=result.data,
            error_type="no_results",
            error_message=f"No STRING ID for {gene_symbol}",
        )
    return _tool_result(
        endpoint_name="get_string_ids",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "species": species,
            "string_id": string_id,
            "raw": result.data,
        },
    )


def interaction_partners(
    string_id: str,
    *,
    gene_symbol: str = "",
    species: int = SPECIES_HUMAN,
    limit: int = DEFAULT_LIMIT,
    required_score: int = DEFAULT_REQUIRED_SCORE,
    network_type: str = DEFAULT_NETWORK_TYPE,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch STRING interaction partners for a resolved STRING ID."""
    cfg = settings or get_settings()
    params = {
        "identifiers": string_id,
        "species": str(species),
        "limit": str(limit),
        "required_score": str(required_score),
        "network_type": network_type,
        "caller_identity": cfg.caller_identity,
    }
    return _request_json(
        endpoint_name="interaction_partners",
        gene_symbol=gene_symbol or string_id,
        path="interaction_partners",
        params=params,
        settings=cfg,
    )


def fetch_interaction_partners(
    gene_symbol: str,
    *,
    species: int = SPECIES_HUMAN,
    limit: int = DEFAULT_LIMIT,
    required_score: int = DEFAULT_REQUIRED_SCORE,
    network_type: str = DEFAULT_NETWORK_TYPE,
    string_id: str | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Resolve STRING ID (if needed) and fetch interaction partners.

    On success, ``data`` includes ``string_id`` and ``partners`` (raw partners list).
    """
    cfg = settings or get_settings()
    resolved_id = string_id
    resolve_payload: Any = None

    if not resolved_id:
        resolved = get_string_ids(gene_symbol, species=species, settings=cfg)
        if not resolved.success:
            return _tool_result(
                endpoint_name="fetch_interaction_partners",
                gene_symbol=gene_symbol,
                request_url=resolved.request_url,
                request_params=resolved.request_params,
                success=False,
                status_code=resolved.status_code,
                data={"get_string_ids": resolved.data},
                error_type=resolved.error_type or "resolve_failed",
                error_message=resolved.error_message or "STRING ID resolve failed",
            )
        resolve_payload = resolved.data
        resolved_id = (resolved.data or {}).get("string_id")
        if not resolved_id:
            return _tool_result(
                endpoint_name="fetch_interaction_partners",
                gene_symbol=gene_symbol,
                request_url=resolved.request_url,
                request_params=resolved.request_params,
                success=False,
                status_code=resolved.status_code,
                data={"get_string_ids": resolved.data},
                error_type="no_results",
                error_message=f"No STRING ID for {gene_symbol}",
            )

    partners = interaction_partners(
        resolved_id,
        gene_symbol=gene_symbol,
        species=species,
        limit=limit,
        required_score=required_score,
        network_type=network_type,
        settings=cfg,
    )
    if not partners.success:
        return _tool_result(
            endpoint_name="fetch_interaction_partners",
            gene_symbol=gene_symbol,
            request_url=partners.request_url,
            request_params=partners.request_params,
            success=False,
            status_code=partners.status_code,
            data={
                "string_id": resolved_id,
                "get_string_ids": resolve_payload,
                "partners": partners.data,
            },
            error_type=partners.error_type or "partners_failed",
            error_message=partners.error_message or "STRING interaction_partners failed",
        )

    return _tool_result(
        endpoint_name="fetch_interaction_partners",
        gene_symbol=gene_symbol,
        request_url=partners.request_url,
        request_params=partners.request_params,
        success=True,
        status_code=partners.status_code,
        data={
            "gene_symbol": gene_symbol,
            "species": species,
            "string_id": resolved_id,
            "get_string_ids": resolve_payload,
            "partners": partners.data,
            "partner_count": len(partners.data) if isinstance(partners.data, list) else None,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "STRING_BASE",
    "SPECIES_HUMAN",
    "DEFAULT_LIMIT",
    "DEFAULT_REQUIRED_SCORE",
    "DEFAULT_NETWORK_TYPE",
    "get_string_ids",
    "interaction_partners",
    "fetch_interaction_partners",
    "prefer_string_id",
]
