"""STRING database client (version-pinned v12).

Resolves a gene/protein to a STRING ID and fetches networks / partners.
Does **not** normalize into evidence records — that belongs in
``normalize/ppi.py`` / ``section_5a``.

Production base::

    https://version-12-0.string-db.org

Legacy public helpers (``get_string_ids``, ``interaction_partners``,
``fetch_interaction_partners``, ``prefer_string_id``) preserve their
behavioral output contracts. Section 5a helpers use
``caller_identity=gene_dossier``.

Note: omitting ``limit`` on ``interaction_partners`` is **not** unlimited
(STRING may still apply a server-side default, e.g. 10 rows).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from typing import Any
from urllib.parse import urlencode, urlparse

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "STRING"
STRING_BASE_URL = "https://version-12-0.string-db.org"
STRING_API_JSON = f"{STRING_BASE_URL}/api/json"
STRING_API_IMAGE = f"{STRING_BASE_URL}/api/highres_image"
STRING_VERSION = "12.0"
# Backward-compatible alias used by older imports/tests.
STRING_BASE = STRING_API_JSON

SPECIES_HUMAN = 9606
DEFAULT_LIMIT = 100
DEFAULT_REQUIRED_SCORE = 400
DEFAULT_NETWORK_TYPE = "functional"
DEFAULT_ADD_NODES = 30
DEFAULT_ADD_COLOR_NODES = 10
DEFAULT_ADD_WHITE_NODES = 20
DEFAULT_NETWORK_FLAVOR = "evidence"
SECTION_5A_CALLER_IDENTITY = "gene_dossier"

_STRING_HOST_RE = re.compile(r"(^|\.)string-db\.org$", re.IGNORECASE)


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


def normalize_taxon_id(value: Any) -> int | None:
    """Strict integer-like NCBI taxon parser.

    Accepts ``9606``, ``"9606"``, ``" 9606 "``. Rejects ``None``, ``""``,
    ``"human"``, floats like ``9606.5``, and arbitrary strings.
    """
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value > 0 else None
    if isinstance(value, float):
        if value.is_integer() and value > 0:
            return int(value)
        return None
    text = str(value).strip()
    if not text or not re.fullmatch(r"[0-9]+", text):
        return None
    try:
        parsed = int(text)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def normalize_string_request_identity(
    *,
    method: str,
    url_path: str,
    identifiers: str,
    species: int | str,
    query_params: dict[str, Any] | None = None,
) -> str:
    """Stable identity: method|versioned_path|identifiers|species|sorted params."""
    params = {
        str(k): str(v)
        for k, v in dict(query_params or {}).items()
        if str(k) not in {"identifiers", "species"}
    }
    sorted_params = tuple(sorted(params.items()))
    path = str(url_path or "").rstrip("/")
    return "|".join(
        [
            str(method or "GET").upper(),
            path,
            str(identifiers or "").strip(),
            str(species),
            repr(sorted_params),
        ]
    )


def _attach_meta(payload: Any, meta: dict[str, Any]) -> Any:
    if isinstance(payload, dict):
        cleaned = {k: v for k, v in payload.items() if k != "_string_meta"}
        return {**cleaned, "_string_meta": meta}
    return {"value": payload, "_string_meta": meta}


def string_meta(data: Any) -> dict[str, Any]:
    """Extract ``_string_meta`` from a ToolResult payload when present."""
    if isinstance(data, dict):
        meta = data.get("_string_meta")
        if isinstance(meta, dict):
            return dict(meta)
    return {}


def lookup_string_meta(tool_result: ToolResult | None, transient: Any | None = None) -> dict[str, Any]:
    """Resolve provenance meta from payload wrap or transient request cache."""
    if tool_result is None:
        return {}
    meta = string_meta(tool_result.data)
    if meta:
        return meta
    if transient is None:
        return {}
    params = dict(tool_result.request_params or {})
    # Reconstruct identity from the request URL path when possible.
    url = str(tool_result.request_url or "")
    path = url.split("?", 1)[0]
    identity = normalize_string_request_identity(
        method="GET",
        url_path=path,
        identifiers=str(params.get("identifiers") or tool_result.gene_symbol),
        species=str(params.get("species") or ""),
        query_params=params,
    )
    cached = transient.get_cached_request(identity)
    if isinstance(cached, dict) and isinstance(cached.get("meta"), dict):
        return dict(cached["meta"])
    return {}


def unwrap_string_payload(data: Any) -> Any:
    """Return the raw JSON payload without ``_string_meta`` wrapping."""
    if isinstance(data, dict) and "_string_meta" in data and "value" in data and len(data) <= 2:
        return data.get("value")
    if isinstance(data, dict) and "_string_meta" in data:
        return {k: v for k, v in data.items() if k != "_string_meta"}
    return data


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
    transient: Any | None = None,
    preserve_list: bool = False,
    base_url: str = STRING_API_JSON,
) -> ToolResult:
    """GET a STRING JSON API path and return :class:`ToolResult`.

    SHA-256 is computed on response bytes **before** parse. Exact bytes are
    retained on the transient request cache when available.
    """
    url = f"{base_url.rstrip('/')}/{path.lstrip('/')}"
    request_params = dict(params)
    requested_url = f"{url}?{urlencode(request_params)}"
    identity = normalize_string_request_identity(
        method="GET",
        url_path=url,
        identifiers=str(request_params.get("identifiers") or gene_symbol),
        species=str(request_params.get("species") or ""),
        query_params=request_params,
    )
    if transient is not None:
        cached = transient.get_cached_request(identity)
        if isinstance(cached, dict) and cached.get("tool_result") is not None:
            return cached["tool_result"]

    try:
        retrieved_at = datetime.now(timezone.utc).isoformat()
        with httpx.Client(
            timeout=settings.http_timeout_seconds, follow_redirects=True
        ) as client:
            response = client.get(url, params=request_params)
        content = bytes(response.content)
        raw_sha = hashlib.sha256(content).hexdigest()
        final_url = str(response.url)
        redirect_history = [str(r.url) for r in response.history]
        content_type = response.headers.get("content-type")
        decoding_method = "utf-8"
        replacement_count = 0
        decode_warning = None
        payload: Any
        try:
            text = content.decode("utf-8")
            payload = json.loads(text)
        except UnicodeDecodeError:
            text = content.decode("utf-8", errors="replace")
            replacement_count = text.count("\ufffd")
            decoding_method = "utf-8_replace"
            if replacement_count:
                decode_warning = (
                    f"UTF-8 decoding replaced {replacement_count} invalid byte "
                    "sequence(s) with U+FFFD; exact response bytes are preserved."
                )
            try:
                payload = json.loads(text)
            except Exception:  # noqa: BLE001
                payload = {"raw_text": text[:4000]}
                decoding_method = "utf-8_replace_unparsed"
        except ValueError:
            try:
                payload = {"raw_text": content.decode("utf-8", errors="replace")[:4000]}
            except Exception:  # noqa: BLE001
                payload = {"raw_text": str(content[:4000])}
            decoding_method = "utf-8_non_json"

        meta = {
            "requested_url": requested_url,
            "final_url": final_url,
            "redirect_history": redirect_history,
            "status_code": response.status_code,
            "content_type": content_type,
            "retrieved_at": retrieved_at,
            "response_body_sha256": raw_sha,
            "response_byte_length": len(content),
            "decoding_method": decoding_method,
            "utf8_replacement_char_count": replacement_count,
            "decode_warning": decode_warning,
            "request_identity": identity,
            "method": "GET",
            "string_version": STRING_VERSION,
            "string_base_url": STRING_BASE_URL,
        }

        success = bool(response.is_success)
        if preserve_list and isinstance(payload, list):
            data_out: Any = payload
        else:
            data_out = _attach_meta(payload, meta)

        result = _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=success,
            status_code=response.status_code,
            data=data_out,
            error_type=None if success else "http_error",
            error_message=None if success else f"HTTP {response.status_code}",
        )
        if transient is not None:
            transient.put_cached_request(
                identity,
                {
                    "tool_result": result,
                    "response_bytes": content,
                    "meta": meta,
                },
            )
        return result
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def _request_bytes(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
    transient: Any | None = None,
    base_url: str = STRING_API_IMAGE,
) -> ToolResult:
    """GET binary content (e.g. high-res PNG) with SHA before any interpretation."""
    url = f"{base_url.rstrip('/')}/{path.lstrip('/')}"
    request_params = dict(params)
    requested_url = f"{url}?{urlencode(request_params)}"
    identity = normalize_string_request_identity(
        method="GET",
        url_path=url,
        identifiers=str(request_params.get("identifiers") or gene_symbol),
        species=str(request_params.get("species") or ""),
        query_params=request_params,
    )
    if transient is not None:
        cached = transient.get_cached_request(identity)
        if isinstance(cached, dict) and cached.get("tool_result") is not None:
            return cached["tool_result"]

    try:
        retrieved_at = datetime.now(timezone.utc).isoformat()
        with httpx.Client(
            timeout=settings.http_timeout_seconds, follow_redirects=True
        ) as client:
            response = client.get(url, params=request_params)
        content = bytes(response.content)
        raw_sha = hashlib.sha256(content).hexdigest()
        final_url = str(response.url)
        content_type = response.headers.get("content-type")
        meta = {
            "requested_url": requested_url,
            "final_url": final_url,
            "redirect_history": [str(r.url) for r in response.history],
            "status_code": response.status_code,
            "content_type": content_type,
            "retrieved_at": retrieved_at,
            "response_body_sha256": raw_sha,
            "response_byte_length": len(content),
            "request_identity": identity,
            "method": "GET",
            "string_version": STRING_VERSION,
            "string_base_url": STRING_BASE_URL,
        }
        success = bool(response.is_success)
        data = {
            "content": content if success else None,
            "byte_size": len(content),
            "sha256": raw_sha,
            "content_type": content_type,
            "final_url": final_url,
            "_string_meta": meta,
        }
        result = _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=success,
            status_code=response.status_code,
            data=data,
            error_type=None if success else "http_error",
            error_message=None if success else f"HTTP {response.status_code}",
        )
        if transient is not None:
            transient.put_cached_request(
                identity,
                {
                    "tool_result": result,
                    "response_bytes": content,
                    "meta": meta,
                },
            )
        return result
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
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
    transient: Any | None = None,
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
        transient=transient,
        preserve_list=True,
    )
    if not result.success:
        return result

    raw = unwrap_string_payload(result.data)
    rows = raw if isinstance(raw, list) else []
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
            "raw": rows,
            "_string_meta": string_meta(result.data) or None,
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
    transient: Any | None = None,
) -> ToolResult:
    """Fetch STRING interaction partners for a resolved STRING ID.

    Omitting ``limit`` is **not** unlimited — STRING may still apply a
    server-side default (observed as ~10 rows without an explicit limit).
    """
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
        transient=transient,
        preserve_list=True,
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
    transient: Any | None = None,
) -> ToolResult:
    """Resolve STRING ID (if needed) and fetch interaction partners.

    On success, ``data`` includes ``gene_symbol``, ``species``, ``string_id``,
    ``get_string_ids``, ``partners``, and ``partner_count``.
    """
    cfg = settings or get_settings()
    resolved_id = string_id
    resolve_payload: Any = None

    if not resolved_id:
        resolved = get_string_ids(
            gene_symbol, species=species, settings=cfg, transient=transient
        )
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
        transient=transient,
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

    partner_rows = unwrap_string_payload(partners.data)
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
            "partners": partner_rows,
            "partner_count": len(partner_rows) if isinstance(partner_rows, list) else None,
        },
    )


def resolve_string_identifier(
    gene_symbol: str,
    *,
    species: int = SPECIES_HUMAN,
    settings: Settings | None = None,
    transient: Any | None = None,
    caller_identity: str = SECTION_5A_CALLER_IDENTITY,
) -> ToolResult:
    """Fail-closed single human STRING ID resolution for Section 5a."""
    cfg = settings or get_settings()
    params = {
        "identifiers": gene_symbol,
        "species": str(species),
        "echo_query": "1",
        "caller_identity": caller_identity,
    }
    result = _request_json(
        endpoint_name="resolve_string_identifier",
        gene_symbol=gene_symbol,
        path="get_string_ids",
        params=params,
        settings=cfg,
        transient=transient,
        preserve_list=True,
    )
    if not result.success:
        return _tool_result(
            endpoint_name="resolve_string_identifier",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data=result.data,
            error_type=result.error_type or "source_unavailable",
            error_message=result.error_message or "STRING get_string_ids failed",
        )

    raw = unwrap_string_payload(result.data)
    rows = [r for r in (raw if isinstance(raw, list) else []) if isinstance(r, dict)]
    target = gene_symbol.strip().upper()
    species_prefix = f"{species}."
    matches: list[dict[str, Any]] = []
    for row in rows:
        preferred = str(row.get("preferredName") or "").strip().upper()
        taxon = normalize_taxon_id(row.get("ncbiTaxonId"))
        sid = str(row.get("stringId") or row.get("string_id") or "").strip()
        if preferred != target:
            continue
        if taxon != int(species):
            continue
        if not sid.startswith(species_prefix):
            continue
        matches.append(row)

    if not matches:
        return _tool_result(
            endpoint_name="resolve_string_identifier",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"raw": rows, "identifier_status": "not_found"},
            error_type="not_found",
            error_message=f"No STRING ID for {gene_symbol}",
        )

    unique_ids = sorted(
        {
            str(m.get("stringId") or m.get("string_id"))
            for m in matches
            if m.get("stringId") or m.get("string_id")
        }
    )
    if len(unique_ids) != 1:
        return _tool_result(
            endpoint_name="resolve_string_identifier",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={
                "raw": rows,
                "matches": matches,
                "identifier_status": "ambiguous",
            },
            error_type="ambiguous",
            error_message=f"Ambiguous STRING ID mapping for {gene_symbol}",
        )

    row = matches[0]
    sid = unique_ids[0]
    return _tool_result(
        endpoint_name="resolve_string_identifier",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "query_item": row.get("queryItem"),
            "string_id": sid,
            "preferred_name": row.get("preferredName"),
            "ncbi_taxon_id": normalize_taxon_id(row.get("ncbiTaxonId")),
            "taxon_name": row.get("taxonName"),
            "annotation": row.get("annotation"),
            "identifier_status": "resolved",
            "raw": rows,
            "_string_meta": string_meta(result.data) or None,
        },
    )


def fetch_network(
    string_id: str,
    *,
    gene_symbol: str = "",
    species: int = SPECIES_HUMAN,
    add_nodes: int = DEFAULT_ADD_NODES,
    required_score: int = DEFAULT_REQUIRED_SCORE,
    network_type: str = DEFAULT_NETWORK_TYPE,
    settings: Settings | None = None,
    transient: Any | None = None,
    caller_identity: str = SECTION_5A_CALLER_IDENTITY,
) -> ToolResult:
    """Fetch bounded functional network JSON for a resolved STRING ID."""
    cfg = settings or get_settings()
    params = {
        "identifiers": string_id,
        "species": str(species),
        "add_nodes": str(int(add_nodes)),
        "required_score": str(int(required_score)),
        "network_type": network_type,
        "caller_identity": caller_identity,
    }
    return _request_json(
        endpoint_name="network",
        gene_symbol=gene_symbol or string_id,
        path="network",
        params=params,
        settings=cfg,
        transient=transient,
        preserve_list=True,
    )


def fetch_network_image(
    string_id: str,
    *,
    gene_symbol: str = "",
    species: int = SPECIES_HUMAN,
    add_color_nodes: int = DEFAULT_ADD_COLOR_NODES,
    add_white_nodes: int = DEFAULT_ADD_WHITE_NODES,
    required_score: int = DEFAULT_REQUIRED_SCORE,
    network_type: str = DEFAULT_NETWORK_TYPE,
    network_flavor: str = DEFAULT_NETWORK_FLAVOR,
    hide_disconnected_nodes: int = 1,
    settings: Settings | None = None,
    transient: Any | None = None,
    caller_identity: str = SECTION_5A_CALLER_IDENTITY,
) -> ToolResult:
    """Fetch official high-resolution evidence-style network PNG bytes."""
    cfg = settings or get_settings()
    params = {
        "identifiers": string_id,
        "species": str(species),
        "add_color_nodes": str(int(add_color_nodes)),
        "add_white_nodes": str(int(add_white_nodes)),
        "required_score": str(int(required_score)),
        "network_flavor": network_flavor,
        "network_type": network_type,
        "hide_disconnected_nodes": str(int(hide_disconnected_nodes)),
        "caller_identity": caller_identity,
    }
    return _request_bytes(
        endpoint_name="highres_image_network",
        gene_symbol=gene_symbol or string_id,
        path="network",
        params=params,
        settings=cfg,
        transient=transient,
        base_url=STRING_API_IMAGE,
    )


def fetch_network_link(
    string_id: str,
    *,
    gene_symbol: str = "",
    species: int = SPECIES_HUMAN,
    add_color_nodes: int = DEFAULT_ADD_COLOR_NODES,
    add_white_nodes: int = DEFAULT_ADD_WHITE_NODES,
    required_score: int = DEFAULT_REQUIRED_SCORE,
    network_type: str = DEFAULT_NETWORK_TYPE,
    network_flavor: str = DEFAULT_NETWORK_FLAVOR,
    settings: Settings | None = None,
    transient: Any | None = None,
    caller_identity: str = SECTION_5A_CALLER_IDENTITY,
) -> ToolResult:
    """Fetch stable STRING web-network URL via ``/api/json/get_link``."""
    cfg = settings or get_settings()
    params = {
        "identifiers": string_id,
        "species": str(species),
        "add_color_nodes": str(int(add_color_nodes)),
        "add_white_nodes": str(int(add_white_nodes)),
        "required_score": str(int(required_score)),
        "network_flavor": network_flavor,
        "network_type": network_type,
        "caller_identity": caller_identity,
    }
    return _request_json(
        endpoint_name="get_link",
        gene_symbol=gene_symbol or string_id,
        path="get_link",
        params=params,
        settings=cfg,
        transient=transient,
        preserve_list=True,
    )


def extract_network_rows(payload: Any) -> list[dict[str, Any]]:
    """Normalize a ``/network`` ToolResult payload into a list of edge dicts."""
    raw = unwrap_string_payload(payload)
    if isinstance(raw, dict) and isinstance(raw.get("edges"), list):
        raw = raw["edges"]
    if not isinstance(raw, list):
        return []
    return [row for row in raw if isinstance(row, dict)]


def extract_link_url(payload: Any) -> str | None:
    """Validate and extract the single HTTPS STRING URL from ``get_link``."""
    raw = unwrap_string_payload(payload)
    if isinstance(raw, dict) and isinstance(raw.get("value"), list):
        raw = raw["value"]
    if isinstance(raw, str):
        candidates = [raw]
    elif isinstance(raw, list):
        candidates = [str(x) for x in raw if x]
    else:
        return None
    if len(candidates) != 1:
        return None
    url = candidates[0].strip()
    if not url.lower().startswith("https://"):
        return None
    host = (urlparse(url).hostname or "").lower()
    if not host or not _STRING_HOST_RE.search(host):
        return None
    return url


def validate_network_image_bytes(
    content: bytes,
    *,
    final_url: str | None = None,
    content_type: str | None = None,
    min_width: int = 400,
    min_height: int = 300,
) -> dict[str, Any]:
    """Validate official STRING network PNG bytes (fail closed)."""
    if not content or len(content) < 64:
        raise ValueError("image payload too small")
    ct = (content_type or "").lower()
    if ct and "html" in ct:
        raise ValueError("image response looks like HTML")
    if ct and not any(token in ct for token in ("image/png", "image/", "octet-stream")):
        raise ValueError(f"unsupported image content type: {content_type}")
    if final_url:
        host = (urlparse(final_url).hostname or "").lower()
        if host != "version-12-0.string-db.org":
            raise ValueError(f"image final host is not version-pinned STRING v12: {host}")
    from PIL import Image, ImageStat
    import io

    with Image.open(io.BytesIO(content)) as img:
        width, height = img.size
        if width < min_width or height < min_height:
            raise ValueError(f"image too small: {width}x{height}")
        converted = img.convert("RGBA")
        stat = ImageStat.Stat(converted)
        extrema = converted.getextrema()
        alpha_range = extrema[3] if len(extrema) > 3 else (255, 255)
        variance = sum(stat.var[:3])
        if alpha_range == (0, 0):
            raise ValueError("image is fully transparent")
        if variance < 2.0:
            raise ValueError("image appears blank or low entropy")
        media = "image/png" if (img.format or "").upper() == "PNG" else f"image/{(img.format or 'unknown').lower()}"
    return {
        "media_type": media,
        "width": width,
        "height": height,
        "byte_size": len(content),
        "pixel_variance": variance,
        "sha256": hashlib.sha256(content).hexdigest(),
    }


__all__ = [
    "SOURCE_NAME",
    "STRING_BASE",
    "STRING_BASE_URL",
    "STRING_API_JSON",
    "STRING_API_IMAGE",
    "STRING_VERSION",
    "SPECIES_HUMAN",
    "DEFAULT_LIMIT",
    "DEFAULT_REQUIRED_SCORE",
    "DEFAULT_NETWORK_TYPE",
    "DEFAULT_ADD_NODES",
    "DEFAULT_ADD_COLOR_NODES",
    "DEFAULT_ADD_WHITE_NODES",
    "DEFAULT_NETWORK_FLAVOR",
    "SECTION_5A_CALLER_IDENTITY",
    "normalize_taxon_id",
    "normalize_string_request_identity",
    "string_meta",
    "lookup_string_meta",
    "unwrap_string_payload",
    "prefer_string_id",
    "get_string_ids",
    "interaction_partners",
    "fetch_interaction_partners",
    "resolve_string_identifier",
    "fetch_network",
    "fetch_network_image",
    "fetch_network_link",
    "extract_network_rows",
    "extract_link_url",
    "validate_network_image_bytes",
]
