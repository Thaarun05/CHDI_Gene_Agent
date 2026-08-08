"""BioGRID interactions client.

Fetches curated protein/genetic interactions for a gene with required filters.
Does **not** normalize into evidence records — that belongs in ``normalize/ppi.py``.

Legacy endpoint (validated defaults)::

    GET https://webservice.thebiogrid.org/interactions/
        ?searchNames=true
        &geneList={symbol}
        &taxId=9606
        &includeInteractors=true
        &includeInteractorInteractions=false
        &selfInteractionsExcluded=true
        &interSpeciesExcluded=true
        &format=jsonExtended
        &max=10000
        &accesskey={{biogrid_accesskey}}

Section 5b uses explicit helpers with selfInteractionsExcluded=false and
interSpeciesExcluded=false, plus /version (plain text) and start pagination.

Requires ``BIOGRID_ACCESSKEY``. Access key is redacted in provenance fields.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "BioGRID"
BIOGRID_BASE_URL = "https://webservice.thebiogrid.org"
BIOGRID_INTERACTIONS_URL = f"{BIOGRID_BASE_URL}/interactions/"
BIOGRID_VERSION_URL = f"{BIOGRID_BASE_URL}/version"

TAX_ID_HUMAN = 9606
DEFAULT_MAX = 10000
DEFAULT_FORMAT = "jsonExtended"


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
    """Redact accesskey for provenance logging."""
    out = {k: v for k, v in params.items() if k != "accesskey"}
    if "accesskey" in params:
        out["accesskey"] = "***"
    return out


def build_interaction_params(
    gene_symbol: str,
    *,
    tax_id: int = TAX_ID_HUMAN,
    max_results: int = DEFAULT_MAX,
    search_names: bool = True,
    include_interactors: bool = True,
    include_interactor_interactions: bool = False,
    self_interactions_excluded: bool = True,
    inter_species_excluded: bool = True,
    fmt: str = DEFAULT_FORMAT,
    accesskey: str,
    start: int | None = None,
) -> dict[str, str]:
    """Build filtered BioGRID interactions query params (accesskey included)."""
    params = {
        "searchNames": "true" if search_names else "false",
        "geneList": gene_symbol.strip(),
        "taxId": str(tax_id),
        "includeInteractors": "true" if include_interactors else "false",
        "includeInteractorInteractions": (
            "true" if include_interactor_interactions else "false"
        ),
        "selfInteractionsExcluded": "true" if self_interactions_excluded else "false",
        "interSpeciesExcluded": "true" if inter_species_excluded else "false",
        "format": fmt,
        "max": str(max_results),
        "accesskey": accesskey,
    }
    if start is not None:
        params["start"] = str(int(start))
    return params


def summarize_interaction(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key BioGRID interaction fields (not evidence)."""
    return {
        "biogrid_interaction_id": row.get("BIOGRID_INTERACTION_ID"),
        "entrez_gene_a": row.get("ENTREZ_GENE_A"),
        "entrez_gene_b": row.get("ENTREZ_GENE_B"),
        "official_symbol_a": row.get("OFFICIAL_SYMBOL_A"),
        "official_symbol_b": row.get("OFFICIAL_SYMBOL_B"),
        "synonyms_a": row.get("SYNONYMS_A"),
        "synonyms_b": row.get("SYNONYMS_B"),
        "experimental_system": row.get("EXPERIMENTAL_SYSTEM"),
        "experimental_system_type": row.get("EXPERIMENTAL_SYSTEM_TYPE"),
        "author": row.get("AUTHOR"),
        "pubmed_id": row.get("PUBMED_ID"),
        "organism_a": row.get("ORGANISM_A"),
        "organism_b": row.get("ORGANISM_B"),
        "throughput": row.get("THROUGHPUT"),
        "source_database": row.get("SOURCE_DATABASE"),
    }


def interactions_as_list(payload: Any) -> list[dict[str, Any]]:
    """Normalize BioGRID ``jsonExtended`` payload to a list of interaction dicts.

    BioGRID often returns a dict keyed by interaction ID; sometimes a list.
    """
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, dict)]
    if isinstance(payload, dict):
        rows: list[dict[str, Any]] = []
        for key, value in payload.items():
            if key == "_biogrid_meta":
                continue
            if not isinstance(value, dict):
                continue
            if (
                "BIOGRID_INTERACTION_ID" in value
                or "OFFICIAL_SYMBOL_A" in value
                or "ENTREZ_GENE_A" in value
            ):
                row = dict(value)
                row.setdefault("BIOGRID_INTERACTION_ID", key)
                rows.append(row)
        if rows:
            return rows
        return [
            v
            for k, v in payload.items()
            if k != "_biogrid_meta" and isinstance(v, dict)
        ]
    return []


def interactions(
    gene_symbol: str,
    *,
    tax_id: int = TAX_ID_HUMAN,
    max_results: int = DEFAULT_MAX,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch filtered BioGRID interactions for ``gene_symbol``.

    Always requires ``geneList`` and the validated filter set. Never calls the
    unfiltered endpoint. Legacy defaults exclude self and cross-species.
    """
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol:
        return _tool_result(
            endpoint_name="interactions",
            gene_symbol=gene_symbol,
            request_url=BIOGRID_INTERACTIONS_URL,
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="BioGRID requires a non-empty geneList filter",
        )

    if not cfg.has_key("BIOGRID_ACCESSKEY"):
        return _tool_result(
            endpoint_name="interactions",
            gene_symbol=symbol,
            request_url=BIOGRID_INTERACTIONS_URL,
            request_params={"geneList": symbol, "taxId": str(tax_id)},
            success=False,
            error_type="requires_key",
            error_message="BIOGRID_ACCESSKEY missing",
        )

    params = build_interaction_params(
        symbol,
        tax_id=tax_id,
        max_results=max_results,
        accesskey=str(cfg.biogrid_accesskey),
    )
    safe = _safe_params(params)
    request_url = f"{BIOGRID_INTERACTIONS_URL}?{urlencode(safe)}"
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(BIOGRID_INTERACTIONS_URL, params=params)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name="interactions",
                gene_symbol=symbol,
                request_url=request_url,
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="interactions",
            gene_symbol=symbol,
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
            endpoint_name="interactions",
            gene_symbol=symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="interactions",
            gene_symbol=symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="interactions",
            gene_symbol=symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_interactions(
    gene_symbol: str,
    *,
    tax_id: int = TAX_ID_HUMAN,
    max_results: int = DEFAULT_MAX,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch BioGRID interactions and attach light summaries.

    On success, ``data`` includes::

        {
          "gene_symbol": ...,
          "tax_id": ...,
          "interactions": <raw payload>,
          "interaction_rows": [...],
          "interaction_summaries": [...],
          "interaction_count": N,
        }
    """
    cfg = settings or get_settings()
    result = interactions(
        gene_symbol, tax_id=tax_id, max_results=max_results, settings=cfg
    )
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_interactions",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"interactions": result.data},
            error_type=result.error_type or "interactions_failed",
            error_message=result.error_message or "BioGRID interactions failed",
        )

    rows = interactions_as_list(result.data)
    summaries = [summarize_interaction(row) for row in rows]
    return _tool_result(
        endpoint_name="fetch_interactions",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "tax_id": tax_id,
            "interactions": result.data,
            "interaction_rows": rows,
            "interaction_summaries": summaries,
            "interaction_count": len(summaries),
        },
    )


# --------------------------------------------------------------------------------------
# Section 5b helpers (do not alter legacy interactions / fetch_interactions schema)
# --------------------------------------------------------------------------------------


def normalize_biogrid_request_identity(
    *,
    method: str,
    url_path: str,
    gene_symbol: str,
    query_params: dict[str, Any] | None = None,
) -> str:
    """Stable identity: method|path|gene|sorted non-secret params."""
    params = {
        str(k): str(v)
        for k, v in dict(query_params or {}).items()
        if str(k).lower() not in {"accesskey", "genelist"}
    }
    sorted_params = tuple(sorted(params.items()))
    path = str(url_path or "").rstrip("/")
    return "|".join(
        [
            str(method or "GET").upper(),
            path,
            str(gene_symbol or "").strip().upper(),
            repr(sorted_params),
        ]
    )


def _attach_meta(payload: Any, meta: dict[str, Any]) -> Any:
    if isinstance(payload, dict):
        cleaned = {k: v for k, v in payload.items() if k != "_biogrid_meta"}
        return {**cleaned, "_biogrid_meta": meta}
    return {"value": payload, "_biogrid_meta": meta}


def biogrid_meta(data: Any) -> dict[str, Any]:
    if isinstance(data, dict):
        meta = data.get("_biogrid_meta")
        if isinstance(meta, dict):
            return dict(meta)
    return {}


def lookup_biogrid_meta(
    tool_result: ToolResult | None, transient: Any | None = None
) -> dict[str, Any]:
    if tool_result is None:
        return {}
    meta = biogrid_meta(tool_result.data)
    if meta:
        return meta
    if transient is None:
        return {}
    params = dict(tool_result.request_params or {})
    url = str(tool_result.request_url or "")
    path = url.split("?", 1)[0]
    identity = normalize_biogrid_request_identity(
        method="GET",
        url_path=path,
        gene_symbol=str(tool_result.gene_symbol or params.get("geneList") or ""),
        query_params=params,
    )
    cached = transient.get_cached_request(identity)
    if isinstance(cached, dict) and isinstance(cached.get("meta"), dict):
        return dict(cached["meta"])
    return {}


def unwrap_biogrid_payload(data: Any) -> Any:
    if isinstance(data, dict) and "_biogrid_meta" in data and "value" in data and len(data) <= 2:
        return data.get("value")
    if isinstance(data, dict) and "_biogrid_meta" in data:
        return {k: v for k, v in data.items() if k != "_biogrid_meta"}
    return data


def fetch_version(
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
    transient: Any | None = None,
) -> ToolResult:
    """GET /version as plain text (exact bytes → SHA → trimmed version string)."""
    cfg = settings or get_settings()
    symbol = (gene_symbol or "").strip()
    if not cfg.has_key("BIOGRID_ACCESSKEY"):
        return _tool_result(
            endpoint_name="version",
            gene_symbol=symbol,
            request_url=BIOGRID_VERSION_URL,
            request_params={},
            success=False,
            error_type="requires_key",
            error_message="BIOGRID_ACCESSKEY missing",
        )

    live_params = {"accesskey": str(cfg.biogrid_accesskey)}
    safe = _safe_params(live_params)
    requested_url = f"{BIOGRID_VERSION_URL}?{urlencode(safe)}"
    identity = normalize_biogrid_request_identity(
        method="GET",
        url_path=BIOGRID_VERSION_URL,
        gene_symbol=symbol or "_version",
        query_params=safe,
    )
    if transient is not None:
        cached = transient.get_cached_request(identity)
        if isinstance(cached, dict) and cached.get("tool_result") is not None:
            return cached["tool_result"]

    try:
        retrieved_at = datetime.now(timezone.utc).isoformat()
        with httpx.Client(timeout=cfg.http_timeout_seconds, follow_redirects=True) as client:
            response = client.get(BIOGRID_VERSION_URL, params=live_params)
        content = bytes(response.content)
        raw_sha = hashlib.sha256(content).hexdigest()
        text = content.decode("utf-8", errors="replace")
        version = text.strip()
        meta = {
            "requested_url": requested_url,
            "final_url": str(response.url),
            "redirect_history": [str(r.url) for r in response.history],
            "status_code": response.status_code,
            "content_type": response.headers.get("content-type"),
            "retrieved_at": retrieved_at,
            "response_body_sha256": raw_sha,
            "response_byte_length": len(content),
            "decoding_method": "utf-8_text_trim",
            "endpoint_name": "version",
            "request_identity": identity,
        }
        data = _attach_meta({"version": version, "raw_text": text}, meta)
        if response.is_success and version:
            result = _tool_result(
                endpoint_name="version",
                gene_symbol=symbol,
                request_url=requested_url,
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=data,
            )
        else:
            result = _tool_result(
                endpoint_name="version",
                gene_symbol=symbol,
                request_url=requested_url,
                request_params=safe,
                success=False,
                status_code=response.status_code,
                data=data,
                error_type="http_error" if not response.is_success else "invalid_version",
                error_message=(
                    f"HTTP {response.status_code}"
                    if not response.is_success
                    else "empty BioGRID version"
                ),
            )
        if transient is not None:
            transient.put_cached_request(
                identity,
                {"tool_result": result, "response_bytes": content, "meta": meta},
            )
        return result
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name="version",
            gene_symbol=symbol,
            request_url=requested_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_interactions_page(
    gene_symbol: str,
    *,
    start: int = 0,
    tax_id: int = TAX_ID_HUMAN,
    max_results: int = DEFAULT_MAX,
    self_interactions_excluded: bool = False,
    inter_species_excluded: bool = False,
    include_interactors: bool = True,
    include_interactor_interactions: bool = False,
    settings: Settings | None = None,
    transient: Any | None = None,
) -> ToolResult:
    """Section-5b page fetch with explicit filters and SHA-before-parse cache."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol:
        return _tool_result(
            endpoint_name="interactions_page",
            gene_symbol=gene_symbol,
            request_url=BIOGRID_INTERACTIONS_URL,
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="BioGRID requires a non-empty geneList filter",
        )
    if not cfg.has_key("BIOGRID_ACCESSKEY"):
        return _tool_result(
            endpoint_name="interactions_page",
            gene_symbol=symbol,
            request_url=BIOGRID_INTERACTIONS_URL,
            request_params={
                "geneList": symbol,
                "taxId": str(tax_id),
                "start": str(start),
            },
            success=False,
            error_type="requires_key",
            error_message="BIOGRID_ACCESSKEY missing",
        )

    live_params = build_interaction_params(
        symbol,
        tax_id=tax_id,
        max_results=max_results,
        include_interactors=include_interactors,
        include_interactor_interactions=include_interactor_interactions,
        self_interactions_excluded=self_interactions_excluded,
        inter_species_excluded=inter_species_excluded,
        accesskey=str(cfg.biogrid_accesskey),
        start=start,
    )
    safe = _safe_params(live_params)
    requested_url = f"{BIOGRID_INTERACTIONS_URL}?{urlencode(safe)}"
    identity = normalize_biogrid_request_identity(
        method="GET",
        url_path=BIOGRID_INTERACTIONS_URL.rstrip("/"),
        gene_symbol=symbol,
        query_params=safe,
    )
    if transient is not None:
        cached = transient.get_cached_request(identity)
        if isinstance(cached, dict) and cached.get("tool_result") is not None:
            return cached["tool_result"]

    try:
        retrieved_at = datetime.now(timezone.utc).isoformat()
        with httpx.Client(timeout=cfg.http_timeout_seconds, follow_redirects=True) as client:
            response = client.get(BIOGRID_INTERACTIONS_URL, params=live_params)
        content = bytes(response.content)
        raw_sha = hashlib.sha256(content).hexdigest()
        try:
            payload: Any = json.loads(content.decode("utf-8"))
            decoding_method = "utf-8_json"
        except Exception:  # noqa: BLE001
            try:
                payload = json.loads(content.decode("utf-8", errors="replace"))
                decoding_method = "utf-8_replace_json"
            except Exception:  # noqa: BLE001
                payload = {"raw_text": content.decode("utf-8", errors="replace")[:4000]}
                decoding_method = "utf-8_non_json"

        rows = interactions_as_list(payload)
        meta = {
            "requested_url": requested_url,
            "final_url": str(response.url),
            "redirect_history": [str(r.url) for r in response.history],
            "status_code": response.status_code,
            "content_type": response.headers.get("content-type"),
            "retrieved_at": retrieved_at,
            "response_body_sha256": raw_sha,
            "response_byte_length": len(content),
            "decoding_method": decoding_method,
            "endpoint_name": "interactions_page",
            "request_identity": identity,
            "page_start": int(start),
            "page_row_count": len(rows),
            "max": int(max_results),
        }
        data = _attach_meta(
            {
                "gene_symbol": symbol,
                "tax_id": tax_id,
                "page_start": int(start),
                "interactions": payload,
                "interaction_rows": rows,
                "interaction_count": len(rows),
            },
            meta,
        )
        if response.is_success:
            result = _tool_result(
                endpoint_name="interactions_page",
                gene_symbol=symbol,
                request_url=requested_url,
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=data,
            )
        else:
            result = _tool_result(
                endpoint_name="interactions_page",
                gene_symbol=symbol,
                request_url=requested_url,
                request_params=safe,
                success=False,
                status_code=response.status_code,
                data=data,
                error_type="http_error",
                error_message=f"HTTP {response.status_code}",
            )
        if transient is not None:
            transient.put_cached_request(
                identity,
                {"tool_result": result, "response_bytes": content, "meta": meta},
            )
        return result
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name="interactions_page",
            gene_symbol=symbol,
            request_url=requested_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_all_interactions_section_5b(
    gene_symbol: str,
    *,
    tax_id: int = TAX_ID_HUMAN,
    max_results: int = DEFAULT_MAX,
    settings: Settings | None = None,
    transient: Any | None = None,
) -> tuple[list[ToolResult], list[dict[str, Any]], list[str]]:
    """Paginate Section-5b interactions. Never silently truncates.

    Returns ``(page_tool_results, annotated_rows, errors)``.
    Each annotated row includes ``_page_start`` and ``_within_page_index``.
    """
    pages: list[ToolResult] = []
    annotated: list[dict[str, Any]] = []
    errors: list[str] = []
    start = 0
    while True:
        page = fetch_interactions_page(
            gene_symbol,
            start=start,
            tax_id=tax_id,
            max_results=max_results,
            self_interactions_excluded=False,
            inter_species_excluded=False,
            settings=settings,
            transient=transient,
        )
        pages.append(page)
        if not page.success:
            errors.append(page.error_message or "interactions_page failed")
            break
        payload = unwrap_biogrid_payload(page.data) or {}
        rows = list(payload.get("interaction_rows") or [])
        for idx, row in enumerate(rows):
            item = dict(row)
            item["_page_start"] = start
            item["_within_page_index"] = idx
            annotated.append(item)
        if len(rows) < max_results:
            break
        start += max_results
    return pages, annotated, errors


__all__ = [
    "SOURCE_NAME",
    "BIOGRID_BASE_URL",
    "BIOGRID_INTERACTIONS_URL",
    "BIOGRID_VERSION_URL",
    "TAX_ID_HUMAN",
    "DEFAULT_MAX",
    "DEFAULT_FORMAT",
    "build_interaction_params",
    "summarize_interaction",
    "interactions_as_list",
    "interactions",
    "fetch_interactions",
    "normalize_biogrid_request_identity",
    "biogrid_meta",
    "lookup_biogrid_meta",
    "unwrap_biogrid_payload",
    "fetch_version",
    "fetch_interactions_page",
    "fetch_all_interactions_section_5b",
]
