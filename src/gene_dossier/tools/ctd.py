"""CTD (Comparative Toxicogenomics Database) client.

Fetches chemical–gene interaction rows via the batch query tool. CTD returns
**TSV**, not JSON — preserve the raw text for ``raw_store`` artifacts before any
normalization in ``normalize/chemicals.py``.

Key endpoint (validated)::

    GET https://ctdbase.org/tools/batchQuery.go
        ?inputType=gene&inputTerms={symbol}&report=cgixns&actionTypes=ANY&format=tsv

Expected TSV columns::

    Input, ChemicalName, ChemicalID, CasRN, GeneSymbol, GeneID,
    Organism, OrganismID, Interaction, InteractionActions, PubMedIDs

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import csv
import io
from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "CTD"
BATCH_QUERY_URL = "https://ctdbase.org/tools/batchQuery.go"

DEFAULT_INPUT_TYPE = "gene"
DEFAULT_REPORT = "cgixns"
DEFAULT_ACTION_TYPES = "ANY"
DEFAULT_FORMAT = "tsv"

# Validated header names from the API map (order may vary; match by name).
EXPECTED_COLUMNS = (
    "Input",
    "ChemicalName",
    "ChemicalID",
    "CasRN",
    "GeneSymbol",
    "GeneID",
    "Organism",
    "OrganismID",
    "Interaction",
    "InteractionActions",
    "PubMedIDs",
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


def build_batch_query_params(
    gene_symbol: str,
    *,
    input_type: str = DEFAULT_INPUT_TYPE,
    report: str = DEFAULT_REPORT,
    action_types: str = DEFAULT_ACTION_TYPES,
    fmt: str = DEFAULT_FORMAT,
) -> dict[str, str]:
    """Build query parameters for the CTD batch query tool."""
    return {
        "inputType": input_type,
        "inputTerms": gene_symbol.strip(),
        "report": report,
        "actionTypes": action_types,
        "format": fmt,
    }


def parse_tsv(raw_tsv: str) -> list[dict[str, str]]:
    """Parse CTD batch-query TSV into row dicts (still not evidence records).

    Skips blank lines. Returns an empty list when there is no header/body.
    """
    text = (raw_tsv or "").lstrip("\ufeff").strip()
    if not text:
        return []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    rows: list[dict[str, str]] = []
    for row in reader:
        if not isinstance(row, dict):
            continue
        # Drop None keys from uneven rows; coerce values to str.
        cleaned = {
            str(k): ("" if v is None else str(v))
            for k, v in row.items()
            if k is not None
        }
        if any(v.strip() for v in cleaned.values()):
            rows.append(cleaned)
    return rows


def summarize_interaction_row(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key CTD interaction fields for later normalization."""
    return {
        "input": row.get("Input"),
        "chemical_name": row.get("ChemicalName"),
        "chemical_id": row.get("ChemicalID"),
        "cas_rn": row.get("CasRN"),
        "gene_symbol": row.get("GeneSymbol"),
        "gene_id": row.get("GeneID"),
        "organism": row.get("Organism"),
        "organism_id": row.get("OrganismID"),
        "interaction": row.get("Interaction"),
        "interaction_actions": row.get("InteractionActions"),
        "pubmed_ids": row.get("PubMedIDs"),
    }


def batch_query(
    gene_symbol: str,
    *,
    input_type: str = DEFAULT_INPUT_TYPE,
    report: str = DEFAULT_REPORT,
    action_types: str = DEFAULT_ACTION_TYPES,
    fmt: str = DEFAULT_FORMAT,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch raw CTD batch-query TSV for ``gene_symbol``.

    On success, ``data`` is::

        {
          "raw_tsv": <str>,
          "content_type": <str | None>,
        }

    Never raises.
    """
    cfg = settings or get_settings()
    params = build_batch_query_params(
        gene_symbol,
        input_type=input_type,
        report=report,
        action_types=action_types,
        fmt=fmt,
    )
    request_url = f"{BATCH_QUERY_URL}?{urlencode(params)}"
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(BATCH_QUERY_URL, params=params)
        text = response.text
        content_type = response.headers.get("content-type")
        payload = {"raw_tsv": text, "content_type": content_type}

        if response.is_success:
            return _tool_result(
                endpoint_name="batch_query",
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="batch_query",
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
            endpoint_name="batch_query",
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="batch_query",
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="batch_query",
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_chemical_gene_interactions(
    gene_symbol: str,
    *,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch CTD chemical–gene interactions and attach parsed row views.

    On success, ``data`` includes::

        {
          "gene_symbol": ...,
          "raw_tsv": <str>,          # preserve for raw_store
          "content_type": ...,
          "rows": [<raw TSV dicts>],
          "interaction_summaries": [...],
          "row_count": N,
        }

    Parsing helpers are convenience only — not evidence normalization.
    Never raises.
    """
    cfg = settings or get_settings()
    result = batch_query(gene_symbol, settings=cfg)
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_chemical_gene_interactions",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data=result.data,
            error_type=result.error_type or "batch_query_failed",
            error_message=result.error_message or "CTD batch query failed",
        )

    raw_tsv = ""
    content_type = None
    if isinstance(result.data, dict):
        raw_tsv = str(result.data.get("raw_tsv") or "")
        content_type = result.data.get("content_type")

    rows = parse_tsv(raw_tsv)
    summaries = [summarize_interaction_row(row) for row in rows]
    return _tool_result(
        endpoint_name="fetch_chemical_gene_interactions",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": gene_symbol,
            "raw_tsv": raw_tsv,
            "content_type": content_type,
            "rows": rows,
            "interaction_summaries": summaries,
            "row_count": len(rows),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "BATCH_QUERY_URL",
    "DEFAULT_INPUT_TYPE",
    "DEFAULT_REPORT",
    "DEFAULT_ACTION_TYPES",
    "DEFAULT_FORMAT",
    "EXPECTED_COLUMNS",
    "build_batch_query_params",
    "parse_tsv",
    "summarize_interaction_row",
    "batch_query",
    "fetch_chemical_gene_interactions",
]
