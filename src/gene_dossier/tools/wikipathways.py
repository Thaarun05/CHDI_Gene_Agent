"""WikiPathways client (bulk JSON + local filtering).

WikiPathways deprecated the old live query web services. The validated workflow
is: download bulk JSON, then filter locally for the gene. Does **not** normalize
into evidence records — that belongs in ``normalize/pathways.py``.

Key endpoints (validated)::

    GET https://www.wikipathways.org/json/findPathwaysByText.json
    GET https://www.wikipathways.org/json/findPathwaysByXref.json

Xref ``pathwayInfo`` fields include id, url, name, species, revision,
description, ncbigene, uniprot, ensembl.

For SREBF2, matched identifiers include ``6721``, ``Q12772``,
``ENSG00000198911``, and ``SREBF2``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any, Iterable

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "WikiPathways"
WIKIPATHWAYS_JSON_BASE = "https://www.wikipathways.org/json"
TEXT_BULK_URL = f"{WIKIPATHWAYS_JSON_BASE}/findPathwaysByText.json"
XREF_BULK_URL = f"{WIKIPATHWAYS_JSON_BASE}/findPathwaysByXref.json"

SPECIES_HUMAN = "Homo sapiens"


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
    url: str,
    request_params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a WikiPathways JSON bulk file and return :class:`ToolResult`."""
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

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


def extract_pathway_info(payload: Any) -> list[dict[str, Any]]:
    """Return ``pathwayInfo`` rows from a bulk JSON payload."""
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, dict)]
    if not isinstance(payload, dict):
        return []
    info = payload.get("pathwayInfo")
    if isinstance(info, list):
        return [row for row in info if isinstance(row, dict)]
    # Some dumps may nest under other keys; keep empty rather than guess wrongly.
    return []


def summarize_pathway(info: dict[str, Any]) -> dict[str, Any]:
    """Extract key WikiPathways pathwayInfo fields (not evidence)."""
    return {
        "id": info.get("id"),
        "url": info.get("url"),
        "name": info.get("name"),
        "species": info.get("species"),
        "revision": info.get("revision"),
        "description": info.get("description"),
        "ncbigene": list(info.get("ncbigene") or [])
        if isinstance(info.get("ncbigene"), list)
        else info.get("ncbigene"),
        "uniprot": list(info.get("uniprot") or [])
        if isinstance(info.get("uniprot"), list)
        else info.get("uniprot"),
        "ensembl": list(info.get("ensembl") or [])
        if isinstance(info.get("ensembl"), list)
        else info.get("ensembl"),
    }


def _as_str_set(values: Iterable[Any] | None) -> set[str]:
    """Normalize identifier iterables to a lowercase string set."""
    out: set[str] = set()
    if values is None:
        return out
    for value in values:
        if value is None:
            continue
        text = str(value).strip()
        if text:
            out.add(text.lower())
    return out


def _field_values(info: dict[str, Any], key: str) -> set[str]:
    """Collect string values from a pathwayInfo field (list or scalar)."""
    raw = info.get(key)
    if raw is None:
        return set()
    if isinstance(raw, list):
        return _as_str_set(raw)
    return _as_str_set([raw])


def pathway_matches(
    info: dict[str, Any],
    *,
    gene_symbol: str | None = None,
    entrez_ids: list[str | int] | None = None,
    uniprot_ids: list[str] | None = None,
    ensembl_ids: list[str] | None = None,
    species: str | None = SPECIES_HUMAN,
    allow_text_symbol_match: bool = False,
) -> bool:
    """True if ``info`` matches any provided identifier (local filter).

    Prefer exact identifier matches on ``ncbigene``, ``uniprot``, ``ensembl``,
    ``genes``, and ``symbols``. Do **not** match ``gene_symbol`` by loose
    substring in name/description unless ``allow_text_symbol_match=True``.

    Matching is OR across identifier types. When ``species`` is set, the pathway
    species must match (case-insensitive).
    """
    if species:
        pathway_species = str(info.get("species") or "").strip().lower()
        if pathway_species and pathway_species != species.strip().lower():
            return False

    wanted_symbols = _as_str_set([gene_symbol] if gene_symbol else [])
    wanted_entrez = _as_str_set(entrez_ids)
    wanted_uniprot = _as_str_set(uniprot_ids)
    wanted_ensembl = _as_str_set(ensembl_ids)

    if not (wanted_symbols or wanted_entrez or wanted_uniprot or wanted_ensembl):
        return False

    if wanted_entrez and (_field_values(info, "ncbigene") & wanted_entrez):
        return True
    if wanted_uniprot and (_field_values(info, "uniprot") & wanted_uniprot):
        return True
    if wanted_ensembl and (_field_values(info, "ensembl") & wanted_ensembl):
        return True

    if wanted_symbols:
        structured_symbols = _field_values(info, "genes") | _field_values(
            info, "symbols"
        )
        if wanted_symbols & structured_symbols:
            return True
        if allow_text_symbol_match:
            name = str(info.get("name") or "").lower()
            desc = str(info.get("description") or "").lower()
            for sym in wanted_symbols:
                if sym and (sym in name or sym in desc):
                    return True

    return False


def filter_pathways(
    pathway_info: list[dict[str, Any]],
    *,
    gene_symbol: str | None = None,
    entrez_ids: list[str | int] | None = None,
    uniprot_ids: list[str] | None = None,
    ensembl_ids: list[str] | None = None,
    species: str | None = SPECIES_HUMAN,
    allow_text_symbol_match: bool = False,
) -> list[dict[str, Any]]:
    """Filter pathwayInfo rows locally for the requested gene identifiers."""
    return [
        row
        for row in pathway_info
        if pathway_matches(
            row,
            gene_symbol=gene_symbol,
            entrez_ids=entrez_ids,
            uniprot_ids=uniprot_ids,
            ensembl_ids=ensembl_ids,
            species=species,
            allow_text_symbol_match=allow_text_symbol_match,
        )
    ]


def download_text_bulk(
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Download WikiPathways ``findPathwaysByText.json`` bulk file."""
    cfg = settings or get_settings()
    return _request_json(
        endpoint_name="download_text_bulk",
        gene_symbol=gene_symbol or "WikiPathways",
        url=TEXT_BULK_URL,
        request_params={"bulk": "findPathwaysByText.json"},
        settings=cfg,
    )


def download_xref_bulk(
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Download WikiPathways ``findPathwaysByXref.json`` bulk file."""
    cfg = settings or get_settings()
    return _request_json(
        endpoint_name="download_xref_bulk",
        gene_symbol=gene_symbol or "WikiPathways",
        url=XREF_BULK_URL,
        request_params={"bulk": "findPathwaysByXref.json"},
        settings=cfg,
    )


def fetch_pathways(
    gene_symbol: str,
    *,
    entrez_ids: list[str | int] | None = None,
    uniprot_ids: list[str] | None = None,
    ensembl_ids: list[str] | None = None,
    species: str | None = SPECIES_HUMAN,
    include_text_bulk: bool = False,
    allow_text_symbol_match: bool = False,
    settings: Settings | None = None,
) -> ToolResult:
    """Download xref bulk JSON and filter locally for the gene.

    Xref bulk is the primary structured source. Text bulk is optional
    (``include_text_bulk=True``) and filtered the same way when present.
    Name/description symbol substring matching is off by default
    (``allow_text_symbol_match=False``).

    On success, ``data`` includes matched pathway summaries plus raw bulk
    payloads for artifact storage. Never raises.
    """
    cfg = settings or get_settings()
    xref = download_xref_bulk(gene_symbol=gene_symbol, settings=cfg)
    if not xref.success:
        return _tool_result(
            endpoint_name="fetch_pathways",
            gene_symbol=gene_symbol,
            request_url=xref.request_url,
            request_params=xref.request_params,
            success=False,
            status_code=xref.status_code,
            data={"xref_bulk": xref.data},
            error_type=xref.error_type or "xref_bulk_failed",
            error_message=xref.error_message
            or "WikiPathways xref bulk download failed",
        )

    xref_rows = extract_pathway_info(xref.data)
    matched = filter_pathways(
        xref_rows,
        gene_symbol=gene_symbol,
        entrez_ids=entrez_ids,
        uniprot_ids=uniprot_ids,
        ensembl_ids=ensembl_ids,
        species=species,
        allow_text_symbol_match=allow_text_symbol_match,
    )
    # Deduplicate by pathway id when possible.
    by_id: dict[str, dict[str, Any]] = {}
    no_id: list[dict[str, Any]] = []
    for row in matched:
        pid = row.get("id")
        if pid:
            by_id[str(pid)] = row
        else:
            no_id.append(row)
    matched_unique = list(by_id.values()) + no_id
    summaries = [summarize_pathway(row) for row in matched_unique]

    text_payload: Any = None
    text_summaries: list[dict[str, Any]] = []
    last_url = xref.request_url
    last_params: dict[str, Any] = {
        "gene_symbol": gene_symbol,
        "entrez_ids": [str(i) for i in (entrez_ids or [])],
        "uniprot_ids": list(uniprot_ids or []),
        "ensembl_ids": list(ensembl_ids or []),
        "species": species,
        "include_text_bulk": include_text_bulk,
        "allow_text_symbol_match": allow_text_symbol_match,
    }
    last_status = xref.status_code

    if include_text_bulk:
        text = download_text_bulk(gene_symbol=gene_symbol, settings=cfg)
        last_url = text.request_url
        last_status = text.status_code
        if not text.success:
            return _tool_result(
                endpoint_name="fetch_pathways",
                gene_symbol=gene_symbol,
                request_url=text.request_url,
                request_params=last_params,
                success=False,
                status_code=text.status_code,
                data={
                    "xref_bulk": xref.data,
                    "pathway_summaries": summaries,
                    "pathway_count": len(summaries),
                    "text_bulk": text.data,
                },
                error_type=text.error_type or "text_bulk_failed",
                error_message=text.error_message
                or "WikiPathways text bulk download failed",
            )
        text_payload = text.data
        text_matched = filter_pathways(
            extract_pathway_info(text.data),
            gene_symbol=gene_symbol,
            entrez_ids=entrez_ids,
            uniprot_ids=uniprot_ids,
            ensembl_ids=ensembl_ids,
            species=species,
            allow_text_symbol_match=allow_text_symbol_match,
        )
        text_summaries = [summarize_pathway(row) for row in text_matched]

    return _tool_result(
        endpoint_name="fetch_pathways",
        gene_symbol=gene_symbol,
        request_url=last_url,
        request_params=last_params,
        success=True,
        status_code=last_status,
        data={
            "gene_symbol": gene_symbol,
            "entrez_ids": [str(i) for i in (entrez_ids or [])],
            "uniprot_ids": list(uniprot_ids or []),
            "ensembl_ids": list(ensembl_ids or []),
            "species": species,
            "allow_text_symbol_match": allow_text_symbol_match,
            "xref_bulk": xref.data,
            "matched_pathways": matched_unique,
            "pathway_summaries": summaries,
            "pathway_count": len(summaries),
            "text_bulk": text_payload,
            "text_pathway_summaries": text_summaries,
            "text_pathway_count": len(text_summaries),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "WIKIPATHWAYS_JSON_BASE",
    "TEXT_BULK_URL",
    "XREF_BULK_URL",
    "SPECIES_HUMAN",
    "extract_pathway_info",
    "summarize_pathway",
    "pathway_matches",
    "filter_pathways",
    "download_text_bulk",
    "download_xref_bulk",
    "fetch_pathways",
]
