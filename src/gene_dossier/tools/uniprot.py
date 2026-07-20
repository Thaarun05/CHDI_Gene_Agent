"""UniProt REST client (uniprotkb/search).

Resolves a reviewed protein accession for a gene symbol. Does **not** normalize into
evidence records — that belongs in ``normalize/protein.py`` / ``gene_identity.py``.

Default human query (per platform spec)::

    (gene_exact:{symbol}) AND (organism_id:9606) AND (reviewed:true)

For SREBF2, the expected UniProt accession is ``Q12772``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "UniProt"
UNIPROT_BASE = "https://rest.uniprot.org"

ORGANISM_HUMAN = 9606
ORGANISM_MOUSE = 10090
ORGANISM_RAT = 10116

DEFAULT_FIELDS = (
    "accession,id,gene_names,protein_name,organism_name,organism_id,xref_ensembl,"
    "cc_function,cc_subcellular_location,cc_disease,ft_domain,ft_repeat"
)


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


def build_search_query(
    gene_symbol: str,
    *,
    organism_id: int = ORGANISM_HUMAN,
    reviewed: bool = True,
) -> str:
    """Build a UniProtKB search query string."""
    parts = [f"(gene_exact:{gene_symbol})", f"(organism_id:{organism_id})"]
    if reviewed:
        parts.append("(reviewed:true)")
    return " AND ".join(parts)


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a UniProt REST path and return JSON as :class:`ToolResult`."""
    url = f"{UNIPROT_BASE}{path}"
    request_url = f"{url}?{urlencode(params)}"
    headers = {"Accept": "application/json"}
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds, headers=headers) as client:
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


def extract_results(search_payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Return the ``results`` list from a UniProt search JSON body."""
    results = search_payload.get("results")
    return results if isinstance(results, list) else []


def summarize_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Extract key fields from one UniProtKB search hit."""
    genes = entry.get("genes") or []
    gene_names: list[str] = []
    for g in genes:
        if not isinstance(g, dict):
            continue
        primary = (g.get("geneName") or {}).get("value")
        if primary:
            gene_names.append(str(primary))
        for syn in g.get("synonyms") or []:
            if isinstance(syn, dict) and syn.get("value"):
                gene_names.append(str(syn["value"]))

    protein = entry.get("proteinDescription") or {}
    recommended = ((protein.get("recommendedName") or {}).get("fullName") or {}).get("value")

    ensembl_ids: list[str] = []
    for xref in entry.get("uniProtKBCrossReferences") or []:
        if isinstance(xref, dict) and xref.get("database") == "Ensembl" and xref.get("id"):
            ensembl_ids.append(str(xref["id"]))

    organism = entry.get("organism") or {}
    return {
        "accession": entry.get("primaryAccession") or entry.get("accession"),
        "uni_protkb_id": entry.get("uniProtkbId") or entry.get("id"),
        "gene_names": gene_names,
        "protein_name": recommended,
        "organism_name": organism.get("scientificName") or organism.get("commonName"),
        "organism_id": (organism.get("taxonId") or organism.get("taxonID")),
        "ensembl_xrefs": ensembl_ids,
    }


def search_reviewed(
    gene_symbol: str,
    *,
    organism_id: int = ORGANISM_HUMAN,
    fields: str = DEFAULT_FIELDS,
    size: int = 25,
    settings: Settings | None = None,
) -> ToolResult:
    """Search reviewed UniProtKB entries for an exact gene symbol + organism.

    On success, ``data`` includes::

        {
          "gene_symbol": ...,
          "organism_id": ...,
          "query": ...,
          "selected_accession": "Q12772" | null,
          "entries": [summarized hits...],
          "raw": <full search json>,
        }
    """
    cfg = settings or get_settings()
    query = build_search_query(gene_symbol, organism_id=organism_id, reviewed=True)
    params = {
        "query": query,
        "format": "json",
        "fields": fields,
        "size": str(size),
    }
    result = _request_json(
        endpoint_name="search_reviewed",
        gene_symbol=gene_symbol,
        path="/uniprotkb/search",
        params=params,
        settings=cfg,
    )
    if not result.success:
        return result

    raw = result.data if isinstance(result.data, dict) else {}
    entries = [summarize_entry(e) for e in extract_results(raw) if isinstance(e, dict)]
    selected = entries[0]["accession"] if entries else None

    return _tool_result(
        endpoint_name="search_reviewed",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "organism_id": organism_id,
            "query": query,
            "selected_accession": selected,
            "entries": entries,
            "raw": raw,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "UNIPROT_BASE",
    "ORGANISM_HUMAN",
    "ORGANISM_MOUSE",
    "ORGANISM_RAT",
    "DEFAULT_FIELDS",
    "build_search_query",
    "search_reviewed",
    "extract_results",
    "summarize_entry",
]
