"""Open Targets Platform GraphQL client.

Fetches disease associations and tractability / chemical-probe annotations for an
Ensembl gene ID. Does **not** normalize into evidence records — that belongs in
``normalize/variants.py`` (associations) and related normalizers.

Key endpoint (validated)::

    POST https://api.platform.opentargets.org/api/v4/graphql

Queries:
- associatedDiseases (page size 1000)
- tractability + chemicalProbes

For SREBF2, the expected Ensembl gene ID is ``ENSG00000198911``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Open Targets"
GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Validated against CHDI_Data_APIs_Gene_Report_SREBF2.md §9.6 / §7.5.
DISEASE_ASSOCIATIONS_QUERY = """
query GeneDiseaseEvidence($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    associatedDiseases(page: { index: 0, size: 1000 }) {
      count
      rows {
        score
        disease {
          id
          name
        }
        datatypeScores {
          id
          score
        }
        datasourceScores {
          id
          score
        }
      }
    }
  }
}
""".strip()

TRACTABILITY_QUERY = """
query GeneTractability($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    tractability {
      modality
      label
      value
    }
    chemicalProbes {
      id
      drugId
      drugFromSourceId
      targetFromSourceId
      mechanismOfAction
      isHighQuality
      origin
      probeMinerScore
      probesDrugsScore
      control
      urls {
        niceName
        url
      }
    }
  }
}
""".strip()

DEFAULT_ENSEMBL_ID_SREBF2 = "ENSG00000198911"


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


def _post_graphql(
    *,
    endpoint_name: str,
    gene_symbol: str,
    query: str,
    variables: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """POST a GraphQL query to Open Targets; return :class:`ToolResult`."""
    body = {"query": query, "variables": variables}
    # Provenance: store operation identity + variables, not the full query text twice.
    request_params = {
        "endpoint": GRAPHQL_URL,
        "operation": endpoint_name,
        "variables": variables,
    }
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
    }
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.post(GRAPHQL_URL, json=body, headers=headers)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if not response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=GRAPHQL_URL,
                request_params=request_params,
                success=False,
                status_code=response.status_code,
                data=payload,
                error_type="http_error",
                error_message=f"HTTP {response.status_code}",
            )

        # GraphQL can return 200 with top-level errors.
        if isinstance(payload, dict) and payload.get("errors"):
            errors = payload.get("errors") or []
            first = errors[0] if errors else {}
            message = (
                first.get("message")
                if isinstance(first, dict)
                else "GraphQL error"
            )
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=GRAPHQL_URL,
                request_params=request_params,
                success=False,
                status_code=response.status_code,
                data=payload,
                error_type="graphql_error",
                error_message=str(message),
            )

        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=GRAPHQL_URL,
            request_params=request_params,
            success=True,
            status_code=response.status_code,
            data=payload,
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=GRAPHQL_URL,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=GRAPHQL_URL,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=GRAPHQL_URL,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def summarize_disease_row(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key fields from one associatedDiseases row (not evidence)."""
    disease = row.get("disease") or {}
    if not isinstance(disease, dict):
        disease = {}
    return {
        "disease_id": disease.get("id"),
        "disease_name": disease.get("name"),
        "score": row.get("score"),
        "datatype_scores": list(row.get("datatypeScores") or []),
        "datasource_scores": list(row.get("datasourceScores") or []),
    }


def summarize_tractability_item(item: dict[str, Any]) -> dict[str, Any]:
    """Extract key fields from one tractability entry (not evidence)."""
    return {
        "modality": item.get("modality"),
        "label": item.get("label"),
        "value": item.get("value"),
    }


def disease_associations(
    ensembl_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch Open Targets disease associations for an Ensembl gene ID.

    On success, ``data`` is the raw GraphQL JSON response.
    """
    cfg = settings or get_settings()
    ensembl_id = ensembl_id.strip()
    return _post_graphql(
        endpoint_name="disease_associations",
        gene_symbol=gene_symbol or ensembl_id,
        query=DISEASE_ASSOCIATIONS_QUERY,
        variables={"ensemblId": ensembl_id},
        settings=cfg,
    )


def tractability(
    ensembl_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch Open Targets tractability and chemical probes for an Ensembl gene ID.

    On success, ``data`` is the raw GraphQL JSON response.
    """
    cfg = settings or get_settings()
    ensembl_id = ensembl_id.strip()
    return _post_graphql(
        endpoint_name="tractability",
        gene_symbol=gene_symbol or ensembl_id,
        query=TRACTABILITY_QUERY,
        variables={"ensemblId": ensembl_id},
        settings=cfg,
    )


def fetch_disease_associations(
    ensembl_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Disease associations with light row summaries for later normalization.

    On success, ``data`` includes::

        {
          "ensembl_id": ...,
          "gene_symbol": ...,
          "approved_symbol": ...,
          "approved_name": ...,
          "count": int | None,
          "association_summaries": [...],
          "raw": <GraphQL payload>,
        }
    """
    cfg = settings or get_settings()
    ensembl_id = ensembl_id.strip()
    result = disease_associations(
        ensembl_id, gene_symbol=gene_symbol, settings=cfg
    )
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_disease_associations",
            gene_symbol=gene_symbol or ensembl_id,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"ensembl_id": ensembl_id, "raw": result.data},
            error_type=result.error_type or "disease_associations_failed",
            error_message=result.error_message
            or "Open Targets disease associations failed",
        )

    target = ((result.data or {}).get("data") or {}).get("target") or {}
    if not isinstance(target, dict):
        target = {}
    associated = target.get("associatedDiseases") or {}
    if not isinstance(associated, dict):
        associated = {}
    rows = associated.get("rows") or []
    summaries = [
        summarize_disease_row(row) for row in rows if isinstance(row, dict)
    ]
    return _tool_result(
        endpoint_name="fetch_disease_associations",
        gene_symbol=gene_symbol or ensembl_id,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "ensembl_id": ensembl_id,
            "gene_symbol": gene_symbol or target.get("approvedSymbol"),
            "approved_symbol": target.get("approvedSymbol"),
            "approved_name": target.get("approvedName"),
            "target_id": target.get("id"),
            "count": associated.get("count"),
            "association_summaries": summaries,
            "raw": result.data,
        },
    )


def fetch_tractability(
    ensembl_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Tractability + chemical probes with light summaries.

    On success, ``data`` includes::

        {
          "ensembl_id": ...,
          "gene_symbol": ...,
          "approved_symbol": ...,
          "approved_name": ...,
          "tractability_summaries": [...],
          "chemical_probes": [...],
          "raw": <GraphQL payload>,
        }
    """
    cfg = settings or get_settings()
    ensembl_id = ensembl_id.strip()
    result = tractability(ensembl_id, gene_symbol=gene_symbol, settings=cfg)
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_tractability",
            gene_symbol=gene_symbol or ensembl_id,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"ensembl_id": ensembl_id, "raw": result.data},
            error_type=result.error_type or "tractability_failed",
            error_message=result.error_message
            or "Open Targets tractability failed",
        )

    target = ((result.data or {}).get("data") or {}).get("target") or {}
    if not isinstance(target, dict):
        target = {}
    tract_rows = target.get("tractability") or []
    tract_summaries = [
        summarize_tractability_item(item)
        for item in tract_rows
        if isinstance(item, dict)
    ]
    probes = target.get("chemicalProbes") or []
    return _tool_result(
        endpoint_name="fetch_tractability",
        gene_symbol=gene_symbol or ensembl_id,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "ensembl_id": ensembl_id,
            "gene_symbol": gene_symbol or target.get("approvedSymbol"),
            "approved_symbol": target.get("approvedSymbol"),
            "approved_name": target.get("approvedName"),
            "target_id": target.get("id"),
            "tractability_summaries": tract_summaries,
            "chemical_probes": probes if isinstance(probes, list) else [],
            "raw": result.data,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "GRAPHQL_URL",
    "DISEASE_ASSOCIATIONS_QUERY",
    "TRACTABILITY_QUERY",
    "DEFAULT_ENSEMBL_ID_SREBF2",
    "summarize_disease_row",
    "summarize_tractability_item",
    "disease_associations",
    "tractability",
    "fetch_disease_associations",
    "fetch_tractability",
]
