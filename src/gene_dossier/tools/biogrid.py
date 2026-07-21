"""BioGRID interactions client.

Fetches curated protein/genetic interactions for a gene with required filters.
Does **not** normalize into evidence records — that belongs in ``normalize/ppi.py``.

Key endpoint (validated)::

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

NOTE: Do **not** call BioGRID without filters — the unfiltered endpoint returns
millions of interactions.

Requires ``BIOGRID_ACCESSKEY``. Access key is redacted in provenance fields.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "BioGRID"
BIOGRID_INTERACTIONS_URL = "https://webservice.thebiogrid.org/interactions/"

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
) -> dict[str, str]:
    """Build filtered BioGRID interactions query params (accesskey included)."""
    return {
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
        # Fallback: all dict values (still filtered to dicts only).
        return [v for v in payload.values() if isinstance(v, dict)]
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
    unfiltered endpoint.
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
            # Live call uses real accesskey; provenance uses redacted ``safe``.
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


__all__ = [
    "SOURCE_NAME",
    "BIOGRID_INTERACTIONS_URL",
    "TAX_ID_HUMAN",
    "DEFAULT_MAX",
    "DEFAULT_FORMAT",
    "build_interaction_params",
    "summarize_interaction",
    "interactions_as_list",
    "interactions",
    "fetch_interactions",
]
