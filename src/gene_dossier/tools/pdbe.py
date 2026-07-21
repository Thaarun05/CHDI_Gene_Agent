"""PDBe (Protein Data Bank in Europe) client.

Fetches experimentally solved structures mapped to a UniProt accession, plus
optional UniProt/SIFTS mapping and entry summaries. Does **not** normalize into
evidence records — that belongs in ``normalize/protein.py``.

Key endpoints (validated)::

    GET https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot}
    GET https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}
    GET https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}

For SREBF2 / ``Q12772``, a validated PDB example is ``1ukl``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "PDBe"
PDBE_API_BASE = "https://www.ebi.ac.uk/pdbe/api"

DEFAULT_ACCESSION_SREBF2 = "Q12772"
DEFAULT_PDB_SREBF2 = "1ukl"


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
    request_params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a PDBe API path and return :class:`ToolResult`."""
    url = f"{PDBE_API_BASE}/{path.lstrip('/')}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        # PDBe may return 404 / empty-data messages for unmapped IDs.
        if response.is_success:
            empty_msg = ""
            if isinstance(payload, str):
                empty_msg = payload
            elif isinstance(payload, dict):
                empty_msg = str(
                    payload.get("message") or payload.get("error") or ""
                )
            if "does not contain any data" in empty_msg.lower():
                return _tool_result(
                    endpoint_name=endpoint_name,
                    gene_symbol=gene_symbol,
                    request_url=url,
                    request_params=request_params,
                    success=False,
                    status_code=response.status_code,
                    data=payload,
                    error_type="no_data",
                    error_message=empty_msg or "PDBe endpoint returned no data",
                )
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


def summarize_best_structure(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key best_structures fields (not evidence)."""
    return {
        "pdb_id": row.get("pdb_id"),
        "chain_id": row.get("chain_id"),
        "unp_start": row.get("unp_start"),
        "unp_end": row.get("unp_end"),
        "pdb_start": row.get("pdb_start"),
        "pdb_end": row.get("pdb_end"),
        "coverage": row.get("coverage"),
        "resolution": row.get("resolution"),
        "experimental_method": row.get("experimental_method"),
    }


def summarize_entry_summary(entry: dict[str, Any]) -> dict[str, Any]:
    """Extract key PDB entry summary fields (not evidence)."""
    return {
        "title": entry.get("title"),
        "experimental_method": entry.get("experimental_method"),
        "experimental_method_class": entry.get("experimental_method_class"),
        "entry_authors": entry.get("entry_authors"),
        "deposition_date": entry.get("deposition_date"),
        "release_date": entry.get("release_date"),
        "revision_date": entry.get("revision_date"),
    }


def extract_best_structures(
    payload: Any, uniprot_accession: str
) -> list[dict[str, Any]]:
    """Return best-structure rows for ``uniprot_accession`` from a PDBe payload."""
    if not isinstance(payload, dict):
        return []
    rows = payload.get(uniprot_accession) or payload.get(uniprot_accession.upper())
    if isinstance(rows, list):
        return [r for r in rows if isinstance(r, dict)]
    return []


def extract_pdb_ids(best_rows: list[dict[str, Any]]) -> list[str]:
    """Unique PDB IDs from best-structure rows (order preserved)."""
    out: list[str] = []
    seen: set[str] = set()
    for row in best_rows:
        pdb_id = row.get("pdb_id")
        if not pdb_id:
            continue
        key = str(pdb_id).lower()
        if key not in seen:
            seen.add(key)
            out.append(str(pdb_id).lower())
    return out


def best_structures(
    uniprot_accession: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PDBe best structures mapped to a UniProt accession."""
    cfg = settings or get_settings()
    accession = uniprot_accession.strip()
    if not accession:
        return _tool_result(
            endpoint_name="best_structures",
            gene_symbol=gene_symbol,
            request_url=f"{PDBE_API_BASE}/mappings/best_structures/",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="UniProt accession is required",
        )
    path = f"mappings/best_structures/{quote(accession, safe='')}"
    return _request_json(
        endpoint_name="best_structures",
        gene_symbol=gene_symbol or accession,
        path=path,
        request_params={"uniprot_accession": accession},
        settings=cfg,
    )


def uniprot_mapping(
    pdb_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PDBe UniProt/SIFTS mappings for a PDB entry."""
    cfg = settings or get_settings()
    pdb = pdb_id.strip().lower()
    if not pdb:
        return _tool_result(
            endpoint_name="uniprot_mapping",
            gene_symbol=gene_symbol,
            request_url=f"{PDBE_API_BASE}/mappings/uniprot/",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="PDB ID is required",
        )
    path = f"mappings/uniprot/{quote(pdb, safe='')}"
    return _request_json(
        endpoint_name="uniprot_mapping",
        gene_symbol=gene_symbol or pdb,
        path=path,
        request_params={"pdb_id": pdb},
        settings=cfg,
    )


def entry_summary(
    pdb_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PDBe entry summary for a PDB ID."""
    cfg = settings or get_settings()
    pdb = pdb_id.strip().lower()
    if not pdb:
        return _tool_result(
            endpoint_name="entry_summary",
            gene_symbol=gene_symbol,
            request_url=f"{PDBE_API_BASE}/pdb/entry/summary/",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="PDB ID is required",
        )
    path = f"pdb/entry/summary/{quote(pdb, safe='')}"
    return _request_json(
        endpoint_name="entry_summary",
        gene_symbol=gene_symbol or pdb,
        path=path,
        request_params={"pdb_id": pdb},
        settings=cfg,
    )


def fetch_structures(
    uniprot_accession: str,
    *,
    gene_symbol: str = "",
    max_entries: int = 5,
    include_mappings: bool = True,
    include_summaries: bool = True,
    settings: Settings | None = None,
) -> ToolResult:
    """Best structures for a UniProt accession, plus optional mapping/summary.

    On success, ``data`` includes best-structure summaries and, for up to
    ``max_entries`` unique PDB IDs, optional UniProt mapping and entry summary
    payloads. Never raises.
    """
    cfg = settings or get_settings()
    accession = uniprot_accession.strip()
    best = best_structures(accession, gene_symbol=gene_symbol, settings=cfg)
    if not best.success:
        return _tool_result(
            endpoint_name="fetch_structures",
            gene_symbol=gene_symbol or accession,
            request_url=best.request_url,
            request_params=best.request_params,
            success=False,
            status_code=best.status_code,
            data={"best_structures": best.data},
            error_type=best.error_type or "best_structures_failed",
            error_message=best.error_message or "PDBe best_structures failed",
        )

    rows = extract_best_structures(best.data, accession)
    structure_summaries = [summarize_best_structure(r) for r in rows]
    pdb_ids = extract_pdb_ids(rows)[: max(0, max_entries)]

    mappings: dict[str, Any] = {}
    summaries: dict[str, Any] = {}
    entry_summaries: list[dict[str, Any]] = []
    last_url = best.request_url
    last_params: dict[str, Any] = {
        "uniprot_accession": accession,
        "pdb_ids": pdb_ids,
        "include_mappings": include_mappings,
        "include_summaries": include_summaries,
    }
    last_status = best.status_code

    for pdb in pdb_ids:
        if include_mappings:
            mapped = uniprot_mapping(pdb, gene_symbol=gene_symbol, settings=cfg)
            last_url = mapped.request_url
            last_status = mapped.status_code
            if not mapped.success:
                return _tool_result(
                    endpoint_name="fetch_structures",
                    gene_symbol=gene_symbol or accession,
                    request_url=mapped.request_url,
                    request_params={**last_params, "failed_pdb_id": pdb},
                    success=False,
                    status_code=mapped.status_code,
                    data={
                        "uniprot_accession": accession,
                        "best_structures": best.data,
                        "structure_summaries": structure_summaries,
                        "pdb_ids": pdb_ids,
                        "uniprot_mappings": mappings,
                        "entry_summaries_raw": summaries,
                        "failed_pdb_id": pdb,
                        "failed_step": "uniprot_mapping",
                    },
                    error_type=mapped.error_type or "uniprot_mapping_failed",
                    error_message=mapped.error_message
                    or f"PDBe UniProt mapping failed for {pdb}",
                )
            mappings[pdb] = mapped.data

        if include_summaries:
            summary = entry_summary(pdb, gene_symbol=gene_symbol, settings=cfg)
            last_url = summary.request_url
            last_status = summary.status_code
            if not summary.success:
                return _tool_result(
                    endpoint_name="fetch_structures",
                    gene_symbol=gene_symbol or accession,
                    request_url=summary.request_url,
                    request_params={**last_params, "failed_pdb_id": pdb},
                    success=False,
                    status_code=summary.status_code,
                    data={
                        "uniprot_accession": accession,
                        "best_structures": best.data,
                        "structure_summaries": structure_summaries,
                        "pdb_ids": pdb_ids,
                        "uniprot_mappings": mappings,
                        "entry_summaries_raw": summaries,
                        "failed_pdb_id": pdb,
                        "failed_step": "entry_summary",
                    },
                    error_type=summary.error_type or "entry_summary_failed",
                    error_message=summary.error_message
                    or f"PDBe entry summary failed for {pdb}",
                )
            summaries[pdb] = summary.data
            # Payload shape: {pdb_id: [entry, ...]}
            entries = []
            if isinstance(summary.data, dict):
                entries = summary.data.get(pdb) or summary.data.get(pdb.upper()) or []
            if isinstance(entries, list):
                for entry in entries:
                    if isinstance(entry, dict):
                        s = summarize_entry_summary(entry)
                        s["pdb_id"] = pdb
                        entry_summaries.append(s)

    return _tool_result(
        endpoint_name="fetch_structures",
        gene_symbol=gene_symbol or accession,
        request_url=last_url,
        request_params=last_params,
        success=True,
        status_code=last_status,
        data={
            "uniprot_accession": accession,
            "gene_symbol": gene_symbol or None,
            "best_structures": best.data,
            "structure_summaries": structure_summaries,
            "structure_count": len(structure_summaries),
            "pdb_ids": pdb_ids,
            "uniprot_mappings": mappings,
            "entry_summaries_raw": summaries,
            "entry_summaries": entry_summaries,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "PDBE_API_BASE",
    "DEFAULT_ACCESSION_SREBF2",
    "DEFAULT_PDB_SREBF2",
    "summarize_best_structure",
    "summarize_entry_summary",
    "extract_best_structures",
    "extract_pdb_ids",
    "best_structures",
    "uniprot_mapping",
    "entry_summary",
    "fetch_structures",
]
