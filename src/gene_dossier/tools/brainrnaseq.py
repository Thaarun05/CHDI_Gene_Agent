"""BrainRNASeq client (CSV bulk download + local gene filter).

Downloads published human/mouse cell-type expression CSVs and filters rows for a
gene. Priority C scaffold: CSV download, not a JSON API. Does **not** normalize
into evidence records — that belongs in ``normalize/expression.py``.

Key endpoints (validated)::

    GET https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-124.csv
        (human)
    GET https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-120.csv
        (mouse)

NOTE: Parse rows where ``gene_id`` / ``id`` matches the requested symbol
(case-insensitive). Prefer exact matches by default to avoid short-symbol
false positives.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import csv
import io
from typing import Any, Literal

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "BrainRNASeq"
HUMAN_CSV_URL = (
    "https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-124.csv"
)
MOUSE_CSV_URL = (
    "https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-120.csv"
)

SpeciesChoice = Literal["human", "mouse"]

# Cell-type column prefixes noted in the validated API map.
HUMAN_CELLTYPE_PREFIXES = (
    "astrocytes_fetal_",
    "astrocytes_mature_",
    "endothelial_",
    "microglia_",
    "neurons_",
    "oligodendrocytes_",
)
MOUSE_CELLTYPE_PREFIXES = (
    "astrocytes_",
    "endothelial_",
    "microglia_macrophage_",
    "myelinating_oligodendrocyte_",
    "neurons_",
    "newly_formed_oligodendrocyte_",
    "opc_",
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


def csv_url_for_species(species: SpeciesChoice) -> str:
    """Return the validated BrainRNASeq CSV URL for ``species``."""
    if species == "mouse":
        return MOUSE_CSV_URL
    return HUMAN_CSV_URL


def parse_csv(raw_csv: str) -> list[dict[str, str]]:
    """Parse BrainRNASeq CSV text into row dicts."""
    text = (raw_csv or "").lstrip("\ufeff").strip()
    if not text:
        return []
    reader = csv.DictReader(io.StringIO(text))
    rows: list[dict[str, str]] = []
    for row in reader:
        if not isinstance(row, dict):
            continue
        cleaned = {
            str(k).strip(): ("" if v is None else str(v).strip())
            for k, v in row.items()
            if k is not None
        }
        if any(v for v in cleaned.values()):
            rows.append(cleaned)
    return rows


def row_matches_gene(
    row: dict[str, Any],
    gene_symbol: str,
    *,
    allow_substring_match: bool = False,
) -> bool:
    """True if ``gene_id`` / ``id`` matches ``gene_symbol``.

    Exact case-insensitive match by default. Optional substring match is off by
    default to avoid short-symbol false positives.
    """
    target = gene_symbol.strip().lower()
    if not target:
        return False
    for key in ("gene_id", "id", "Gene", "gene", "symbol"):
        value = str(row.get(key) or "").strip().lower()
        if not value:
            continue
        if value == target:
            return True
        if allow_substring_match and target in value:
            return True
    return False


def filter_gene_rows(
    rows: list[dict[str, Any]],
    gene_symbol: str,
    *,
    allow_substring_match: bool = False,
) -> list[dict[str, Any]]:
    """Filter CSV rows for ``gene_symbol``."""
    return [
        row
        for row in rows
        if row_matches_gene(
            row, gene_symbol, allow_substring_match=allow_substring_match
        )
    ]


def summarize_expression_row(
    row: dict[str, Any],
    *,
    species: SpeciesChoice = "human",
) -> dict[str, Any]:
    """Extract gene identifiers plus cell-type expression columns (not evidence)."""
    prefixes = (
        MOUSE_CELLTYPE_PREFIXES if species == "mouse" else HUMAN_CELLTYPE_PREFIXES
    )
    cell_types: dict[str, Any] = {}
    for key, value in row.items():
        key_s = str(key)
        if any(key_s.startswith(prefix) for prefix in prefixes):
            cell_types[key_s] = value
    return {
        "gene_id": row.get("gene_id") or row.get("Gene") or row.get("gene"),
        "id": row.get("id"),
        "species": species,
        "cell_type_values": cell_types,
        "cell_type_count": len(cell_types),
    }


def download_csv(
    species: SpeciesChoice = "human",
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Download a BrainRNASeq CSV bulk file (raw text preserved)."""
    cfg = settings or get_settings()
    url = csv_url_for_species(species)
    request_params = {"species": species, "url": url}
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(url)
        text = response.text
        payload = {
            "raw_csv": text,
            "content_type": response.headers.get("content-type"),
            "species": species,
        }
        if response.is_success:
            return _tool_result(
                endpoint_name="download_csv",
                gene_symbol=gene_symbol or species,
                request_url=url,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
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
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_gene_expression(
    gene_symbol: str,
    *,
    species: SpeciesChoice = "human",
    allow_substring_match: bool = False,
    settings: Settings | None = None,
) -> ToolResult:
    """Download CSV and return rows matching ``gene_symbol``.

    On success, ``data`` includes raw CSV, matched rows, and light summaries.
    Never raises.
    """
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol:
        return _tool_result(
            endpoint_name="fetch_gene_expression",
            gene_symbol=gene_symbol,
            request_url=csv_url_for_species(species),
            request_params={"species": species},
            success=False,
            error_type="invalid_request",
            error_message="gene_symbol is required",
        )

    downloaded = download_csv(species, gene_symbol=symbol, settings=cfg)
    if not downloaded.success:
        return _tool_result(
            endpoint_name="fetch_gene_expression",
            gene_symbol=symbol,
            request_url=downloaded.request_url,
            request_params=downloaded.request_params,
            success=False,
            status_code=downloaded.status_code,
            data=downloaded.data,
            error_type=downloaded.error_type or "download_failed",
            error_message=downloaded.error_message or "BrainRNASeq CSV download failed",
        )

    raw_csv = ""
    if isinstance(downloaded.data, dict):
        raw_csv = str(downloaded.data.get("raw_csv") or "")
    rows = parse_csv(raw_csv)
    matched = filter_gene_rows(
        rows, symbol, allow_substring_match=allow_substring_match
    )
    summaries = [
        summarize_expression_row(row, species=species) for row in matched
    ]
    return _tool_result(
        endpoint_name="fetch_gene_expression",
        gene_symbol=symbol,
        request_url=downloaded.request_url,
        request_params={
            "species": species,
            "gene_symbol": symbol,
            "allow_substring_match": allow_substring_match,
        },
        success=True,
        status_code=downloaded.status_code,
        data={
            "gene_symbol": symbol,
            "species": species,
            "raw_csv": raw_csv,
            "content_type": (
                downloaded.data.get("content_type")
                if isinstance(downloaded.data, dict)
                else None
            ),
            "matched_rows": matched,
            "expression_summaries": summaries,
            "match_count": len(matched),
            "row_count_total": len(rows),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "HUMAN_CSV_URL",
    "MOUSE_CSV_URL",
    "HUMAN_CELLTYPE_PREFIXES",
    "MOUSE_CELLTYPE_PREFIXES",
    "csv_url_for_species",
    "parse_csv",
    "row_matches_gene",
    "filter_gene_rows",
    "summarize_expression_row",
    "download_csv",
    "fetch_gene_expression",
]
