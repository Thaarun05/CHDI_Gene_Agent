"""Ensembl REST client (lookup/symbol).

Resolves a gene symbol to Ensembl gene ID, location, biotype, and canonical
transcript when available. Does **not** normalize into evidence records.

Endpoint::

    GET https://rest.ensembl.org/lookup/symbol/{species}/{symbol}
        ?content-type=application/json
        [&expand=1]   # include Transcript list when requesting canonical transcript

For SREBF2 / homo_sapiens, the expected Ensembl gene ID is ``ENSG00000198911``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Ensembl"
ENSEMBL_BASE = "https://rest.ensembl.org"

SPECIES_HUMAN = "homo_sapiens"
SPECIES_MOUSE = "mus_musculus"
SPECIES_RAT = "rattus_norvegicus"


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
    """GET an Ensembl REST path and return JSON as :class:`ToolResult`."""
    url = f"{ENSEMBL_BASE}{path}"
    # Ensembl also accepts content-type as a query param (validated Postman style).
    query = {"content-type": "application/json", **params}
    request_url = f"{url}?{urlencode(query)}"
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds, headers=headers) as client:
            response = client.get(url, params=query)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:2000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=query,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
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
            request_params=query,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def extract_canonical_transcript(lookup_payload: dict[str, Any]) -> str | None:
    """Return the canonical transcript ID from an expanded lookup payload, if present."""
    # Prefer explicit field when Ensembl provides it.
    for key in ("canonical_transcript", "canonicalTranscript"):
        value = lookup_payload.get(key)
        if isinstance(value, str) and value.strip():
            # Sometimes annotated as "ENST...stable_id.version"
            return value.split(".")[0] if value else None

    transcripts = lookup_payload.get("Transcript") or lookup_payload.get("transcripts") or []
    if not isinstance(transcripts, list):
        return None
    for tr in transcripts:
        if not isinstance(tr, dict):
            continue
        is_canonical = tr.get("is_canonical") in (1, True, "1")
        if is_canonical:
            tid = tr.get("id") or tr.get("stable_id")
            if tid:
                return str(tid)
    return None


def summarize_lookup(payload: dict[str, Any]) -> dict[str, Any]:
    """Extract the key identity fields from a lookup/symbol JSON body."""
    return {
        "ensembl_gene_id": payload.get("id"),
        "display_name": payload.get("display_name") or payload.get("displayName"),
        "biotype": payload.get("biotype"),
        "seq_region_name": payload.get("seq_region_name"),
        "start": payload.get("start"),
        "end": payload.get("end"),
        "strand": payload.get("strand"),
        "assembly_name": payload.get("assembly_name"),
        "species": payload.get("species"),
        "canonical_transcript": extract_canonical_transcript(payload),
        "description": payload.get("description"),
    }


def lookup_symbol(
    gene_symbol: str,
    *,
    species: str = SPECIES_HUMAN,
    expand: bool = True,
    settings: Settings | None = None,
) -> ToolResult:
    """Lookup a gene symbol via Ensembl REST ``/lookup/symbol/{species}/{symbol}``.

    When ``expand=True`` (default), requests transcript expansion so a canonical
    transcript can be recovered when available.

    On success, ``data`` includes both the raw payload and a ``summary`` dict with
    ``ensembl_gene_id``, location, biotype, and ``canonical_transcript``.
    """
    cfg = settings or get_settings()
    path = f"/lookup/symbol/{species}/{gene_symbol}"
    params: dict[str, Any] = {}
    if expand:
        params["expand"] = "1"
    result = _request_json(
        endpoint_name="lookup_symbol",
        gene_symbol=gene_symbol,
        path=path,
        params=params,
        settings=cfg,
    )
    if not result.success:
        return result

    raw = result.data if isinstance(result.data, dict) else {}
    summary = summarize_lookup(raw)
    return _tool_result(
        endpoint_name="lookup_symbol",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params={**result.request_params, "species": species},
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "species": species,
            "summary": summary,
            "raw": raw,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "ENSEMBL_BASE",
    "SPECIES_HUMAN",
    "SPECIES_MOUSE",
    "SPECIES_RAT",
    "lookup_symbol",
    "summarize_lookup",
    "extract_canonical_transcript",
]
