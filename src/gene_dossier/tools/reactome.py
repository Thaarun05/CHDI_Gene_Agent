"""Reactome ContentService client.

Fetches pathway membership for a UniProt accession and optional pathway detail.
Does **not** normalize into evidence records — that belongs in ``normalize/pathways.py``.

Key endpoints (validated)::

    GET https://reactome.org/ContentService/data/mapping/UniProt/{accession}/pathways
    GET https://reactome.org/ContentService/data/query/{stId}

Build links::

    https://reactome.org/content/detail/{stId}
    https://reactome.org/PathwayBrowser/#/{stId}&FLG={accession}

For SREBF2, the expected UniProt accession is ``Q12772``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Reactome"
CONTENT_SERVICE_BASE = "https://reactome.org/ContentService"
DETAIL_BASE = "https://reactome.org/content/detail"
BROWSER_BASE = "https://reactome.org/PathwayBrowser/#"


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


def pathway_detail_url(st_id: str) -> str:
    """Return the Reactome content-detail URL for a stable pathway ID."""
    return f"{DETAIL_BASE}/{quote(st_id, safe='')}"


def pathway_browser_url(st_id: str, uniprot_accession: str) -> str:
    """Return the Pathway Browser URL flagged for a UniProt accession."""
    return f"{BROWSER_BASE}/{quote(st_id, safe='')}&FLG={quote(uniprot_accession, safe='')}"


def summarize_pathway(
    pathway: dict[str, Any],
    *,
    uniprot_accession: str | None = None,
) -> dict[str, Any]:
    """Extract key fields from a Reactome pathway mapping row (not evidence)."""
    st_id = pathway.get("stId") or pathway.get("st_id")
    summary: dict[str, Any] = {
        "db_id": pathway.get("dbId") or pathway.get("db_id"),
        "st_id": st_id,
        "st_id_version": pathway.get("stIdVersion") or pathway.get("st_id_version"),
        "display_name": pathway.get("displayName") or pathway.get("display_name"),
        "species_name": pathway.get("speciesName") or pathway.get("species_name"),
        "doi": pathway.get("doi"),
        "has_diagram": pathway.get("hasDiagram") or pathway.get("has_diagram"),
        "release_date": pathway.get("releaseDate") or pathway.get("release_date"),
        "last_updated_date": pathway.get("lastUpdatedDate")
        or pathway.get("last_updated_date"),
        "schema_class": pathway.get("schemaClass") or pathway.get("schema_class"),
    }
    if isinstance(st_id, str) and st_id.strip():
        summary["detail_url"] = pathway_detail_url(st_id)
        if uniprot_accession:
            summary["browser_url"] = pathway_browser_url(st_id, uniprot_accession)
    return summary


def summarize_pathway_detail(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract key fields from a Reactome ``/data/query/{stId}`` payload."""
    summations = payload.get("summation") or []
    texts: list[str] = []
    if isinstance(summations, list):
        for item in summations:
            if isinstance(item, dict) and item.get("text"):
                texts.append(str(item["text"]))

    pmids: list[int | str] = []
    lit_refs = payload.get("literatureReference") or []
    if isinstance(lit_refs, list):
        for ref in lit_refs:
            if not isinstance(ref, dict):
                continue
            pmid = ref.get("pubMedIdentifier")
            if pmid is not None:
                pmids.append(pmid)

    st_id = payload.get("stId") or payload.get("st_id")
    return {
        "st_id": st_id,
        "display_name": payload.get("displayName") or payload.get("display_name"),
        "summation_texts": texts,
        "pubmed_ids": pmids,
        "detail_url": pathway_detail_url(str(st_id)) if st_id else None,
    }


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    request_params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a Reactome ContentService path and return :class:`ToolResult`."""
    url = f"{CONTENT_SERVICE_BASE}/{path.lstrip('/')}"
    headers = {"Accept": "application/json"}
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, headers=headers)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:2000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=url,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=url,
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
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def pathways_by_uniprot(
    uniprot_accession: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch Reactome pathways mapped to a UniProt accession.

    On success, ``data`` is the raw pathway list from ContentService.
    """
    cfg = settings or get_settings()
    accession = uniprot_accession.strip()
    path = f"data/mapping/UniProt/{quote(accession, safe='')}/pathways"
    params = {"uniprot_accession": accession}
    return _request_json(
        endpoint_name="pathways_by_uniprot",
        gene_symbol=gene_symbol or accession,
        path=path,
        request_params=params,
        settings=cfg,
    )


def pathway_detail(
    st_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch Reactome pathway detail for a stable ID (e.g. ``R-HSA-1655829``)."""
    cfg = settings or get_settings()
    stable_id = st_id.strip()
    path = f"data/query/{quote(stable_id, safe='')}"
    params = {"st_id": stable_id}
    return _request_json(
        endpoint_name="pathway_detail",
        gene_symbol=gene_symbol or stable_id,
        path=path,
        request_params=params,
        settings=cfg,
    )


def fetch_pathways(
    uniprot_accession: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch pathways for a UniProt accession and attach light summaries + links.

    On success, ``data`` includes::

        {
          "uniprot_accession": ...,
          "gene_symbol": ...,
          "pathways": <raw list>,
          "pathway_summaries": [...],
          "pathway_count": N,
        }

    Never raises.
    """
    cfg = settings or get_settings()
    accession = uniprot_accession.strip()
    result = pathways_by_uniprot(
        accession, gene_symbol=gene_symbol, settings=cfg
    )
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_pathways",
            gene_symbol=gene_symbol or accession,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"uniprot_accession": accession, "pathways": result.data},
            error_type=result.error_type or "pathways_failed",
            error_message=result.error_message or "Reactome pathways fetch failed",
        )

    raw = result.data if isinstance(result.data, list) else []
    summaries = [
        summarize_pathway(row, uniprot_accession=accession)
        for row in raw
        if isinstance(row, dict)
    ]
    return _tool_result(
        endpoint_name="fetch_pathways",
        gene_symbol=gene_symbol or accession,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "uniprot_accession": accession,
            "gene_symbol": gene_symbol or None,
            "pathways": result.data,
            "pathway_summaries": summaries,
            "pathway_count": len(summaries),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "CONTENT_SERVICE_BASE",
    "DETAIL_BASE",
    "BROWSER_BASE",
    "pathway_detail_url",
    "pathway_browser_url",
    "summarize_pathway",
    "summarize_pathway_detail",
    "pathways_by_uniprot",
    "pathway_detail",
    "fetch_pathways",
]
