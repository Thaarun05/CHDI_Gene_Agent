"""CTD (Comparative Toxicogenomics Database) client.

Legacy gene-scoped batch query (``batchQuery.go``) is preserved for workflow /
registry callers. Section 6a uses the official bulk chemical–gene interaction
download instead (``CTD_chem_gene_ixns.tsv.gz``) because the website returns
ALTCHA HTML to automated batchQuery requests.

Never raises from public fetch helpers: failures return :class:`ToolResult`.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
from typing import Any, Iterator
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "CTD"
BATCH_QUERY_URL = "https://ctdbase.org/tools/batchQuery.go"
CHEM_GENE_IXNS_BULK_URL = "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz"
CHEM_GENE_IXNS_SOURCE_KEY = "ctd_chem_gene_ixns"

DEFAULT_INPUT_TYPE = "gene"
DEFAULT_REPORT = "cgixns"
DEFAULT_ACTION_TYPES = "ANY"
DEFAULT_FORMAT = "tsv"

# Legacy batch-query columns (closed set for that endpoint).
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

# Section 6a bulk required columns (open schema: extras allowed).
BULK_REQUIRED_COLUMNS = (
    "ChemicalName",
    "ChemicalID",
    "CasRN",
    "GeneSymbol",
    "GeneID",
    "GeneForms",
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
    """Fetch raw CTD batch-query TSV for ``gene_symbol``. Never raises."""
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
    except Exception as exc:  # noqa: BLE001
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
    """Fetch CTD chemical–gene interactions and attach parsed row views."""
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


# ---------------------------------------------------------------------------
# Section 6a bulk chemical–gene interactions (official TSV.gz)
# ---------------------------------------------------------------------------


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def parse_bulk_comment_metadata(text_prefix: str) -> dict[str, Any]:
    """Extract CTD bulk comment metadata (e.g. Report created timestamp)."""
    report_created = None
    comment_lines: list[str] = []
    for line in (text_prefix or "").splitlines():
        stripped = line.strip()
        if not stripped.startswith("#"):
            break
        comment_lines.append(stripped)
        lower = stripped.lstrip("#").strip()
        if lower.lower().startswith("report created:"):
            report_created = lower.split(":", 1)[1].strip()
    return {"ctd_report_created": report_created, "comment_lines": comment_lines}


def validate_bulk_header(fieldnames: list[str] | None) -> dict[str, Any]:
    """Source-level header validation: required columns must exist; extras OK."""
    names = [str(n) for n in (fieldnames or []) if n is not None]
    missing = [c for c in BULK_REQUIRED_COLUMNS if c not in names]
    extras = [c for c in names if c not in BULK_REQUIRED_COLUMNS]
    return {
        "ok": not missing,
        "fieldnames": names,
        "missing_required_columns": missing,
        "extra_columns": extras,
    }


def iter_bulk_tsv_rows(
    gzip_bytes: bytes,
) -> tuple[dict[str, Any], Iterator[dict[str, str]]]:
    """Validate gzip + header, then yield data rows.

    Source-level validation only. Does not inspect row content.
    Raises ``ValueError`` when gzip/header validation fails.

    Official CTD dumps declare columns in a comment after ``# Fields:``::

        # Fields:
        # ChemicalName\\tChemicalID\\t...

    Fixture / alternate dumps may use a normal first non-comment TSV header
    line instead. Both forms are accepted.

    Rows are streamed from the gzip; the full decompressed TSV is not retained.
    """
    if not gzip_bytes:
        raise ValueError("empty_gzip")
    try:
        gz = gzip.GzipFile(fileobj=io.BytesIO(gzip_bytes))
        text_stream = io.TextIOWrapper(gz, encoding="utf-8", errors="replace", newline="")
    except OSError as exc:
        raise ValueError(f"gzip_decompress_failed:{exc}") from exc

    comment_lines: list[str] = []
    comment_fieldnames: list[str] | None = None
    expect_fields = False
    first_body_line: str | None = None
    header_probe_decompressed_bytes = 0

    try:
        while True:
            line = text_stream.readline()
            if line == "":
                break
            header_probe_decompressed_bytes += len(
                line.encode("utf-8", errors="replace")
            )
            if not comment_lines and line.startswith("\ufeff"):
                line = line.lstrip("\ufeff")
            stripped = line.strip()
            if stripped.startswith("#") or (comment_fieldnames is None and first_body_line is None and not stripped):
                if stripped.startswith("#"):
                    comment_lines.append(stripped)
                    inner = stripped.lstrip("#").strip()
                    if inner.lower().startswith("fields:"):
                        expect_fields = True
                        continue
                    if expect_fields and "\t" in inner:
                        comment_fieldnames = [c.strip() for c in inner.split("\t")]
                        expect_fields = False
                        continue
                    expect_fields = False
                continue
            # First non-comment, non-empty body line
            first_body_line = line.rstrip("\n\r")
            break
    except OSError as exc:
        text_stream.close()
        raise ValueError(f"gzip_decompress_failed:{exc}") from exc

    meta = parse_bulk_comment_metadata("\n".join(comment_lines) + "\n")

    if comment_fieldnames:
        fieldnames = comment_fieldnames
        header_origin = "comment_fields"
        pending_data_line = first_body_line
    else:
        if not first_body_line:
            text_stream.close()
            raise ValueError("empty_tsv_body")
        fieldnames = [c.strip() for c in first_body_line.split("\t")]
        header_origin = "tsv_header_line"
        pending_data_line = None

    header = validate_bulk_header(fieldnames)
    if not header["ok"]:
        text_stream.close()
        raise ValueError(
            "missing_required_columns:" + ",".join(header["missing_required_columns"])
        )
    meta = {
        **meta,
        **header,
        # Bytes read while locating the header / first body line only — not the
        # full decompressed TSV size (stream continues for data rows).
        "header_probe_decompressed_bytes": header_probe_decompressed_bytes,
        "header_origin": header_origin,
    }

    def _rows() -> Iterator[dict[str, str]]:
        try:
            if pending_data_line is not None:
                yield from _parse_data_line(pending_data_line, fieldnames)
            while True:
                line = text_stream.readline()
                if line == "":
                    break
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                yield from _parse_data_line(line.rstrip("\n\r"), fieldnames)
        finally:
            text_stream.close()

    return meta, _rows()


def _parse_data_line(line: str, fieldnames: list[str]) -> Iterator[dict[str, str]]:
    values = line.split("\t")
    cleaned: dict[str, str] = {}
    for idx, name in enumerate(fieldnames):
        cleaned[name] = values[idx].strip() if idx < len(values) else ""
    if len(values) > len(fieldnames):
        for extra_i, val in enumerate(values[len(fieldnames) :]):
            cleaned[f"_extra_{extra_i}"] = val.strip()
    if any(cleaned.values()):
        yield cleaned


def download_chem_gene_ixns_bulk(
    *,
    settings: Settings | None = None,
    url: str = CHEM_GENE_IXNS_BULK_URL,
) -> ToolResult:
    """Download the official CTD chemical–gene interactions gzip. Never raises."""
    cfg = settings or get_settings()
    params: dict[str, Any] = {"url": url}
    try:
        timeout = max(float(cfg.http_timeout_seconds), 300.0)
        with httpx.Client(timeout=timeout, follow_redirects=True) as client:
            response = client.get(url)
        content = response.content or b""
        content_type = response.headers.get("content-type")
        payload = {
            "content": content,
            "content_type": content_type,
            "byte_size": len(content),
            "sha256": sha256_bytes(content) if content else None,
            "final_url": str(response.url),
            "source_key": CHEM_GENE_IXNS_SOURCE_KEY,
        }
        if response.is_success and content:
            try:
                with gzip.GzipFile(fileobj=io.BytesIO(content)) as gz:
                    head = gz.read(256 * 1024).decode("utf-8", errors="replace")
                payload["probe_metadata"] = parse_bulk_comment_metadata(head)
            except Exception as probe_exc:  # noqa: BLE001
                return _tool_result(
                    endpoint_name="download_chem_gene_ixns_bulk",
                    gene_symbol="",
                    request_url=url,
                    request_params=params,
                    success=False,
                    status_code=response.status_code,
                    data=payload,
                    error_type="gzip_probe_failed",
                    error_message=str(probe_exc),
                )
            return _tool_result(
                endpoint_name="download_chem_gene_ixns_bulk",
                gene_symbol="",
                request_url=url,
                request_params=params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="download_chem_gene_ixns_bulk",
            gene_symbol="",
            request_url=url,
            request_params=params,
            success=False,
            status_code=response.status_code,
            data=payload,
            error_type="http_error" if not response.is_success else "empty_body",
            error_message=(
                f"HTTP {response.status_code}"
                if not response.is_success
                else "empty CTD bulk body"
            ),
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name="download_chem_gene_ixns_bulk",
            gene_symbol="",
            request_url=url,
            request_params=params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="download_chem_gene_ixns_bulk",
            gene_symbol="",
            request_url=url,
            request_params=params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name="download_chem_gene_ixns_bulk",
            gene_symbol="",
            request_url=url,
            request_params=params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


__all__ = [
    "SOURCE_NAME",
    "BATCH_QUERY_URL",
    "CHEM_GENE_IXNS_BULK_URL",
    "CHEM_GENE_IXNS_SOURCE_KEY",
    "DEFAULT_INPUT_TYPE",
    "DEFAULT_REPORT",
    "DEFAULT_ACTION_TYPES",
    "DEFAULT_FORMAT",
    "EXPECTED_COLUMNS",
    "BULK_REQUIRED_COLUMNS",
    "build_batch_query_params",
    "parse_tsv",
    "summarize_interaction_row",
    "batch_query",
    "fetch_chemical_gene_interactions",
    "sha256_bytes",
    "parse_bulk_comment_metadata",
    "validate_bulk_header",
    "iter_bulk_tsv_rows",
    "download_chem_gene_ixns_bulk",
]
