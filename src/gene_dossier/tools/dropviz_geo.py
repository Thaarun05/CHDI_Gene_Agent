"""DropViz processed-matrix client (NCBI GEO GSE116470, offline analysis).

Section 2c uses the *published* DropViz metacell matrix from GEO rather than the
live DropViz Shiny server, which is unreliable and whose plot exports are broken
server-side. This module acquires the supplementary matrix, profiles it, and —
only when the value semantics are established — derives a population-level
expression ranking for one target gene.

Key endpoint::

    GET https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116470/suppl/
        GSE116470_metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv.gz

Scientific gate: a filename containing "metacells" does not establish that the
values are counts. :func:`classify_value_semantics` inspects the actual numbers
and returns one of ``raw_counts``, ``count_compatible``, ``normalized_expression``,
``transformed_expression`` or ``unresolved``. Transcripts-per-100,000 is computed
only for count-compatible data; a normalized source scale is preserved and
labelled with its own unit; anything unresolved yields no quantitative chart.

Never raises: transport and parse failures return
:class:`~gene_dossier.models.ToolResult` or structured error dicts.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import math
import statistics
from dataclasses import dataclass, field
from typing import Any, Iterable, Iterator, Sequence
from urllib.parse import urlparse

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "GEO"

GEO_ACCESSION = "GSE116470"
METACELL_FILENAME = (
    "GSE116470_metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv.gz"
)

# Official NCBI hosts only. The FTP-over-HTTPS path is primary; the GEO download
# CGI is an equally official fallback when the suppl path layout changes.
ALLOWED_HOSTS = frozenset({"ftp.ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov"})
_SUPPL_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116470/suppl"
_DOWNLOAD_CGI = "https://www.ncbi.nlm.nih.gov/geo/download/"

REQUEST_HEADERS = {
    "User-Agent": "GeneDossier/0.1.0 (research; provenance-first gene dossier client)",
    "Accept": "application/gzip,application/octet-stream,*/*;q=0.8",
}
SAFE_RESPONSE_HEADER_KEYS = ("content-type", "content-length", "server", "last-modified")

GZIP_MAGIC = b"\x1f\x8b"

# Bumped whenever the derivation changes in a way that alters published numbers.
CALCULATION_VERSION = "dropviz_geo_v1"

# Expected shape of the published matrix. Used as validation *signal* only:
# a mismatch is recorded, never forced onto the data.
EXPECTED_GENE_ROWS = 32_307
EXPECTED_POPULATION_COLUMNS = 565
SHAPE_TOLERANCE_FRACTION = 0.05

GENE_COLUMN_CANDIDATES = ("", "gene", "genes", "symbol", "gene_symbol", "feature", "id")

# -- value-semantics vocabulary -------------------------------------------------------
VALUE_SEMANTICS_RAW_COUNTS = "raw_counts"
VALUE_SEMANTICS_COUNT_COMPATIBLE = "count_compatible"
VALUE_SEMANTICS_NORMALIZED = "normalized_expression"
VALUE_SEMANTICS_TRANSFORMED = "transformed_expression"
VALUE_SEMANTICS_UNRESOLVED = "unresolved"

COUNT_COMPATIBLE_SEMANTICS = frozenset(
    {VALUE_SEMANTICS_RAW_COUNTS, VALUE_SEMANTICS_COUNT_COMPATIBLE}
)

# Classification thresholds, recorded verbatim in the profile so the decision is
# auditable rather than implicit.
INTEGRAL_FRACTION_THRESHOLD = 0.999
RAW_COUNT_MIN_MAXIMUM = 10.0
NORMALIZED_COLUMN_SUM_REL_TOL = 0.01

RANKING_METRIC_PER_100K = "transcripts_per_100k"
EXPRESSION_UNIT_PER_100K = "transcripts_per_100k"
AXIS_LABEL_PER_100K = "Transcripts per 100,000"

RANK_STATUS_SUCCESS = "success"
RANK_STATUS_NORMALIZATION_UNRESOLVED = "normalization_unresolved"

CONFIDENCE_INTERVAL_METHOD_UNRESOLVED = "method_unresolved"

LABEL_MAPPING_RESOLVED = "resolved"
LABEL_MAPPING_PARTIAL = "partial"
LABEL_MAPPING_UNRESOLVED = "unresolved"

TOP_N_DEFAULT = 10

CHART_TITLE = "DropViz populations with highest expression"
CHART_SOURCE_NOTE = "Derived from GSE116470 processed metacell counts"


# --------------------------------------------------------------------------------------
# URL resolution
# --------------------------------------------------------------------------------------
def supplementary_url(filename: str = METACELL_FILENAME) -> str:
    """Return the official NCBI suppl-path URL for one GEO supplementary file."""
    return f"{_SUPPL_BASE}/{filename}"


def download_cgi_url(
    filename: str = METACELL_FILENAME, *, accession: str = GEO_ACCESSION
) -> str:
    """Return the official GEO download-CGI URL for one supplementary file."""
    return f"{_DOWNLOAD_CGI}?acc={accession}&format=file&file={filename}"


def candidate_download_urls(
    filename: str = METACELL_FILENAME, *, accession: str = GEO_ACCESSION
) -> tuple[str, ...]:
    """Official download URLs to try in order (both NCBI-hosted)."""
    return (supplementary_url(filename), download_cgi_url(filename, accession=accession))


def is_allowlisted_url(url: str | None) -> bool:
    """True when ``url`` is https and hosted on an approved NCBI host."""
    if not url:
        return False
    parsed = urlparse(url.strip())
    if parsed.scheme not in {"http", "https"}:
        return False
    return (parsed.hostname or "").lower() in ALLOWED_HOSTS


# --------------------------------------------------------------------------------------
# Acquisition
# --------------------------------------------------------------------------------------
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


def _safe_response_headers(headers: Any) -> dict[str, str]:
    out: dict[str, str] = {}
    if headers is None:
        return out
    for key in SAFE_RESPONSE_HEADER_KEYS:
        value = headers.get(key)
        if value is not None and str(value).strip():
            out[key] = str(value)
    return out


def _looks_like_html(content: bytes) -> bool:
    head = content[:512].lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html") or b"<html" in head[:200]


def validate_gzip_payload(content: bytes) -> dict[str, Any]:
    """Validate downloaded bytes as a complete gzip member.

    Returns ``{"ok": bool, "error_type": str | None, "error_message": str | None,
    "decompressed_head": bytes}``. A truncated archive is detected by attempting
    a bounded decompression of the head, which raises ``EOFError`` when the
    stream ends prematurely.
    """
    if not content:
        return {"ok": False, "error_type": "empty_download", "error_message": "no bytes"}
    if _looks_like_html(content):
        return {
            "ok": False,
            "error_type": "html_masquerading_as_data",
            "error_message": "response body is HTML, not a gzip archive",
        }
    if not content.startswith(GZIP_MAGIC):
        return {
            "ok": False,
            "error_type": "invalid_gzip_magic",
            "error_message": "payload does not start with the gzip magic bytes",
        }
    # Decompress the whole stream in bounded chunks: cheap enough for a ~40 MB
    # archive and the only reliable truncation check. A byte-size floor is not
    # usable here because a well-compressed archive can be legitimately small.
    try:
        head = b""
        with gzip.GzipFile(fileobj=io.BytesIO(content)) as fh:
            while True:
                chunk = fh.read(1 << 20)
                if not chunk:
                    break
                if len(head) < 4096:
                    head += chunk[:4096]
    except EOFError:
        return {
            "ok": False,
            "error_type": "truncated_download",
            "error_message": "gzip stream ended before the final block",
        }
    except (OSError, gzip.BadGzipFile) as exc:
        return {
            "ok": False,
            "error_type": "malformed_gzip",
            "error_message": str(exc),
        }
    return {"ok": True, "error_type": None, "error_message": None, "decompressed_head": head}


def download_metacell_matrix(
    *,
    gene_symbol: str = "",
    filename: str = METACELL_FILENAME,
    accession: str = GEO_ACCESSION,
    settings: Settings | None = None,
    urls: Sequence[str] | None = None,
) -> ToolResult:
    """Download the GEO processed metacell matrix (never raises).

    Each attempted URL is one real network request; the caller records one
    ApiRun per attempt. On success ``data`` carries ``content`` (bytes) for the
    section node to persist as the immutable raw artifact.
    """
    cfg = settings or get_settings()
    attempts: list[dict[str, Any]] = []
    candidates = tuple(urls) if urls else candidate_download_urls(filename, accession=accession)

    last: ToolResult | None = None
    for url in candidates:
        if not is_allowlisted_url(url):
            last = _tool_result(
                endpoint_name="geo_supplementary_download",
                gene_symbol=gene_symbol,
                request_url=url,
                request_params={"accession": accession, "filename": filename},
                success=False,
                error_type="url_not_allowlisted",
                error_message=f"host is not an approved NCBI host: {url!r}",
            )
            attempts.append({"url": url, "error_type": "url_not_allowlisted"})
            continue

        request_params = {
            "accession": accession,
            "filename": filename,
            "url": url,
            "request_headers": dict(REQUEST_HEADERS),
            "follow_redirects": True,
        }
        try:
            with httpx.Client(
                timeout=cfg.http_timeout_seconds, follow_redirects=True
            ) as client:
                response = client.get(url, headers=REQUEST_HEADERS)
        except httpx.HTTPError as exc:
            attempts.append({"url": url, "error_type": "network_error", "detail": str(exc)})
            last = _tool_result(
                endpoint_name="geo_supplementary_download",
                gene_symbol=gene_symbol,
                request_url=url,
                request_params=request_params,
                success=False,
                error_type="network_error",
                error_message=str(exc),
            )
            continue

        content = response.content or b""
        meta = {
            "accession": accession,
            "supplementary_filename": filename,
            "resolved_url": str(response.url) if response.url else url,
            "byte_size": len(content),
            "sha256": hashlib.sha256(content).hexdigest() if content else None,
            "content_type": response.headers.get("content-type"),
            "response_headers": _safe_response_headers(response.headers),
            "compression": "gzip",
        }
        if not response.is_success:
            attempts.append({"url": url, "status_code": response.status_code})
            last = _tool_result(
                endpoint_name="geo_supplementary_download",
                gene_symbol=gene_symbol,
                request_url=str(response.url) if response.url else url,
                request_params=request_params,
                success=False,
                status_code=response.status_code,
                data={**meta, "attempts": list(attempts)},
                error_type="http_error",
                error_message=f"HTTP {response.status_code}",
            )
            continue

        check = validate_gzip_payload(content)
        if not check["ok"]:
            attempts.append({"url": url, "error_type": check["error_type"]})
            last = _tool_result(
                endpoint_name="geo_supplementary_download",
                gene_symbol=gene_symbol,
                request_url=str(response.url) if response.url else url,
                request_params=request_params,
                success=False,
                status_code=response.status_code,
                data={
                    **meta,
                    "attempts": list(attempts),
                    "raw_text_preview": content[:400].decode("utf-8", "replace"),
                },
                error_type=str(check["error_type"]),
                error_message=str(check["error_message"]),
            )
            continue

        attempts.append({"url": url, "status_code": response.status_code, "ok": True})
        return _tool_result(
            endpoint_name="geo_supplementary_download",
            gene_symbol=gene_symbol,
            request_url=str(response.url) if response.url else url,
            request_params=request_params,
            success=True,
            status_code=response.status_code,
            data={**meta, "content": content, "attempts": list(attempts)},
        )

    return last or _tool_result(
        endpoint_name="geo_supplementary_download",
        gene_symbol=gene_symbol,
        request_url=supplementary_url(filename),
        request_params={"accession": accession, "filename": filename},
        success=False,
        error_type="no_candidate_urls",
        error_message="no official download URL was attempted",
    )


# --------------------------------------------------------------------------------------
# Matrix scan (single streaming pass)
# --------------------------------------------------------------------------------------
@dataclass
class MatrixScan:
    """Everything one streaming pass over the matrix yields."""

    ok: bool
    error_type: str | None = None
    error_message: str | None = None
    delimiter: str = ","
    gene_column_index: int = 0
    gene_column_name: str = ""
    gene_column_status: str = "unresolved"
    population_labels: list[str] = field(default_factory=list)
    duplicate_population_labels: list[str] = field(default_factory=list)
    empty_population_labels: int = 0
    gene_row_count: int = 0
    duplicate_gene_symbols: int = 0
    malformed_row_count: int = 0
    malformed_value_count: int = 0
    missing_value_count: int = 0
    column_totals: list[float] = field(default_factory=list)
    column_valid_counts: list[int] = field(default_factory=list)
    target_matches: int = 0
    target_row: list[float | None] = field(default_factory=list)
    target_source_symbol: str | None = None
    value_count: int = 0
    integral_value_count: int = 0
    decimal_value_count: int = 0
    negative_value_count: int = 0
    nonfinite_value_count: int = 0
    minimum: float | None = None
    maximum: float | None = None
    truncated: bool = False


def normalize_gene_symbol(symbol: str) -> str:
    """Documented symbol normalization: strip surrounding whitespace/quotes only.

    Casing is preserved for the comparison key via ``casefold`` at match time;
    no substring or alias expansion is performed.
    """
    return str(symbol or "").strip().strip('"').strip("'")


def _resolve_gene_column(header: Sequence[str]) -> tuple[int, str, str]:
    """Return ``(index, name, status)`` for the gene identifier column."""
    if not header:
        return 0, "", "unresolved"
    first = normalize_gene_symbol(header[0]).casefold()
    if first in GENE_COLUMN_CANDIDATES:
        return 0, header[0], "resolved_by_header"
    for idx, name in enumerate(header):
        if normalize_gene_symbol(name).casefold() in GENE_COLUMN_CANDIDATES and name != "":
            return idx, name, "resolved_by_header"
    return 0, header[0], "assumed_first_column"


def _parse_value(raw: str) -> tuple[float | None, str | None]:
    """Parse one matrix cell. Returns ``(value, problem)``."""
    text = (raw or "").strip()
    if text == "" or text.upper() in {"NA", "NAN", "NULL", "."}:
        return None, "missing"
    try:
        value = float(text)
    except ValueError:
        return None, "malformed"
    if not math.isfinite(value):
        return None, "nonfinite"
    return value, None


def scan_matrix(
    stream: Iterable[str],
    *,
    target_gene: str | None = None,
    max_rows: int | None = None,
) -> MatrixScan:
    """Single streaming pass: profile the matrix and collect column totals.

    Also captures the exact target-gene row when ``target_gene`` is given.
    Genes are expected as rows and populations as columns; an ambiguous
    orientation fails closed.
    """
    reader = csv.reader(stream)
    try:
        header = next(reader)
    except StopIteration:
        return MatrixScan(ok=False, error_type="empty_matrix", error_message="no header row")

    if len(header) < 2:
        return MatrixScan(
            ok=False,
            error_type="malformed_matrix",
            error_message=f"header has {len(header)} field(s); expected a gene column plus populations",
        )

    gene_idx, gene_name, gene_status = _resolve_gene_column(header)
    labels_raw = [h for i, h in enumerate(header) if i != gene_idx]
    labels = [normalize_gene_symbol(x) for x in labels_raw]

    # A fully numeric header means the file has no population labels, which makes
    # the orientation ambiguous rather than merely unlabelled.
    numeric_header = 0
    for label in labels:
        try:
            float(label)
        except ValueError:
            continue
        numeric_header += 1
    if labels and numeric_header == len(labels):
        return MatrixScan(
            ok=False,
            error_type="ambiguous_orientation",
            error_message="header row is entirely numeric; population labels not identifiable",
        )

    seen_labels: dict[str, int] = {}
    duplicates: list[str] = []
    empty_labels = 0
    for label in labels:
        if not label:
            empty_labels += 1
            continue
        seen_labels[label] = seen_labels.get(label, 0) + 1
        if seen_labels[label] == 2:
            duplicates.append(label)

    ncols = len(labels)
    scan = MatrixScan(
        ok=True,
        delimiter=",",
        gene_column_index=gene_idx,
        gene_column_name=gene_name,
        gene_column_status=gene_status,
        population_labels=labels,
        duplicate_population_labels=duplicates,
        empty_population_labels=empty_labels,
        column_totals=[0.0] * ncols,
        column_valid_counts=[0] * ncols,
    )

    target_key = normalize_gene_symbol(target_gene).casefold() if target_gene else None
    seen_genes: set[str] = set()
    expected_fields = len(header)

    for row in reader:
        if not row or all(not str(cell).strip() for cell in row):
            continue
        if len(row) != expected_fields:
            scan.malformed_row_count += 1
            continue

        symbol = normalize_gene_symbol(row[gene_idx])
        scan.gene_row_count += 1
        key = symbol.casefold()
        if key in seen_genes:
            scan.duplicate_gene_symbols += 1
        else:
            seen_genes.add(key)

        is_target = target_key is not None and key == target_key
        if is_target:
            scan.target_matches += 1
            scan.target_source_symbol = symbol
            target_values: list[float | None] = [None] * ncols

        col = 0
        for i, cell in enumerate(row):
            if i == gene_idx:
                continue
            value, problem = _parse_value(cell)
            if problem == "missing":
                scan.missing_value_count += 1
            elif problem == "malformed":
                scan.malformed_value_count += 1
            elif problem == "nonfinite":
                scan.nonfinite_value_count += 1

            if value is not None:
                scan.value_count += 1
                if float(value).is_integer():
                    scan.integral_value_count += 1
                else:
                    scan.decimal_value_count += 1
                if value < 0:
                    scan.negative_value_count += 1
                scan.minimum = value if scan.minimum is None else min(scan.minimum, value)
                scan.maximum = value if scan.maximum is None else max(scan.maximum, value)
                scan.column_totals[col] += value
                scan.column_valid_counts[col] += 1
            if is_target:
                target_values[col] = value
            col += 1

        if is_target:
            scan.target_row = target_values

        if max_rows is not None and scan.gene_row_count >= max_rows:
            scan.truncated = True
            break

    return scan


def open_matrix_stream(content: bytes) -> Iterator[str]:
    """Yield decoded text lines from gzip-compressed matrix bytes."""
    with gzip.GzipFile(fileobj=io.BytesIO(content)) as fh:
        wrapper = io.TextIOWrapper(fh, encoding="utf-8", newline="")
        for line in wrapper:
            yield line


# --------------------------------------------------------------------------------------
# Value semantics
# --------------------------------------------------------------------------------------
def _column_sum_distribution(totals: Sequence[float]) -> dict[str, float | None]:
    finite = [t for t in totals if math.isfinite(t)]
    if not finite:
        return {"minimum": None, "median": None, "maximum": None, "relative_spread": None}
    lo, hi = min(finite), max(finite)
    med = statistics.median(finite)
    spread = None if med == 0 else (hi - lo) / abs(med)
    return {"minimum": lo, "median": med, "maximum": hi, "relative_spread": spread}


def classify_value_semantics(
    scan: MatrixScan,
    *,
    documented_semantics: str | None = None,
    documentation_reference: str | None = None,
    documented_unit: str | None = None,
) -> dict[str, Any]:
    """Classify what the matrix values actually are, with recorded evidence.

    An explicit ``documented_semantics`` (backed by a citation) wins; otherwise
    the decision is made from the observed numbers. Nonnegativity alone never
    implies counts — integrality is the discriminator.
    """
    total = scan.value_count
    integral_fraction = (scan.integral_value_count / total) if total else 0.0
    decimal_fraction = (scan.decimal_value_count / total) if total else 0.0
    col_dist = _column_sum_distribution(scan.column_totals)

    evidence = {
        "value_count": total,
        "integer_fraction": integral_fraction,
        "decimal_value_fraction": decimal_fraction,
        "decimal_value_count": scan.decimal_value_count,
        "negative_value_count": scan.negative_value_count,
        "nonfinite_value_count": scan.nonfinite_value_count,
        "minimum": scan.minimum,
        "maximum": scan.maximum,
        "column_sum_distribution": col_dist,
        "thresholds": {
            "integral_fraction_threshold": INTEGRAL_FRACTION_THRESHOLD,
            "raw_count_min_maximum": RAW_COUNT_MIN_MAXIMUM,
            "normalized_column_sum_rel_tol": NORMALIZED_COLUMN_SUM_REL_TOL,
        },
        "source_documentation_reference": documentation_reference,
    }

    if documented_semantics:
        status = str(documented_semantics)
        return {
            "value_semantics_status": status,
            "basis": "source_documentation",
            "documented_unit": documented_unit,
            "evidence": evidence,
        }

    if total == 0:
        return {
            "value_semantics_status": VALUE_SEMANTICS_UNRESOLVED,
            "basis": "no_parsable_values",
            "documented_unit": None,
            "evidence": evidence,
        }

    if scan.negative_value_count > 0:
        return {
            "value_semantics_status": VALUE_SEMANTICS_TRANSFORMED,
            "basis": "negative_values_present",
            "documented_unit": None,
            "evidence": evidence,
        }

    if integral_fraction >= INTEGRAL_FRACTION_THRESHOLD:
        maximum = scan.maximum or 0.0
        status = (
            VALUE_SEMANTICS_RAW_COUNTS
            if maximum >= RAW_COUNT_MIN_MAXIMUM
            else VALUE_SEMANTICS_COUNT_COMPATIBLE
        )
        return {
            "value_semantics_status": status,
            "basis": "integral_nonnegative_values",
            "documented_unit": None,
            "evidence": evidence,
        }

    spread = col_dist.get("relative_spread")
    if spread is not None and spread <= NORMALIZED_COLUMN_SUM_REL_TOL:
        # Near-constant column sums indicate a per-population normalization whose
        # unit must come from documentation, not from us.
        return {
            "value_semantics_status": VALUE_SEMANTICS_NORMALIZED,
            "basis": "near_constant_column_sums",
            "documented_unit": documented_unit,
            "evidence": evidence,
        }

    return {
        "value_semantics_status": VALUE_SEMANTICS_UNRESOLVED,
        "basis": "non_integral_values_without_documented_scale",
        "documented_unit": None,
        "evidence": evidence,
    }


def build_matrix_profile(
    scan: MatrixScan,
    *,
    semantics: dict[str, Any],
    source_filename: str = METACELL_FILENAME,
    source_sha256: str | None = None,
    source_url: str | None = None,
) -> dict[str, Any]:
    """Assemble ``dropviz_geo_matrix_profile.json`` content."""
    population_count = len(scan.population_labels)
    shape_ok_rows = _within_tolerance(scan.gene_row_count, EXPECTED_GENE_ROWS)
    shape_ok_cols = _within_tolerance(population_count, EXPECTED_POPULATION_COLUMNS)
    return {
        "calculation_version": CALCULATION_VERSION,
        "accession": GEO_ACCESSION,
        "source_filename": source_filename,
        "source_sha256": source_sha256,
        "source_url": source_url,
        "delimiter": scan.delimiter,
        "header_present": True,
        "gene_identifier_column": {
            "index": scan.gene_column_index,
            "name": scan.gene_column_name,
            "status": scan.gene_column_status,
        },
        "orientation": "genes_rows_populations_columns",
        "gene_row_count": scan.gene_row_count,
        "population_column_count": population_count,
        "expected_gene_rows": EXPECTED_GENE_ROWS,
        "expected_population_columns": EXPECTED_POPULATION_COLUMNS,
        "shape_matches_expected": bool(shape_ok_rows and shape_ok_cols),
        "duplicate_population_labels": scan.duplicate_population_labels,
        "empty_population_labels": scan.empty_population_labels,
        "duplicate_gene_symbols": scan.duplicate_gene_symbols,
        "malformed_row_count": scan.malformed_row_count,
        "malformed_value_count": scan.malformed_value_count,
        "missing_value_count": scan.missing_value_count,
        "value_semantics_status": semantics["value_semantics_status"],
        "value_semantics_basis": semantics["basis"],
        "documented_unit": semantics.get("documented_unit"),
        "value_semantics_evidence": semantics["evidence"],
    }


def _within_tolerance(actual: int, expected: int) -> bool:
    if expected <= 0:
        return False
    return abs(actual - expected) / expected <= SHAPE_TOLERANCE_FRACTION


# --------------------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------------------
def ranking_presentation(
    value_semantics_status: str, *, documented_unit: str | None = None
) -> dict[str, Any]:
    """Metric, unit and axis label implied by the established value semantics."""
    if value_semantics_status in COUNT_COMPATIBLE_SEMANTICS:
        return {
            "ranking_metric": RANKING_METRIC_PER_100K,
            "expression_unit": EXPRESSION_UNIT_PER_100K,
            "axis_label": AXIS_LABEL_PER_100K,
            "chartable": True,
        }
    if value_semantics_status == VALUE_SEMANTICS_NORMALIZED:
        unit = documented_unit or ""
        return {
            "ranking_metric": unit or "documented_source_metric",
            "expression_unit": unit,
            "axis_label": unit,
            "chartable": bool(unit),
        }
    return {
        "ranking_metric": None,
        "expression_unit": None,
        "axis_label": None,
        "chartable": False,
    }


def build_ranking_records(
    scan: MatrixScan,
    *,
    value_semantics_status: str,
    documented_unit: str | None = None,
) -> dict[str, Any]:
    """Build per-population ranking records plus exclusions.

    Returns ``{"status", "records", "excluded", "presentation"}``. ``records``
    preserve source column order; sorting is the caller's derivation step.
    """
    presentation = ranking_presentation(
        value_semantics_status, documented_unit=documented_unit
    )
    counts_mode = value_semantics_status in COUNT_COMPATIBLE_SEMANTICS
    if not counts_mode and not presentation["chartable"]:
        # Either the semantics are transformed/unresolved, or a normalized scale
        # arrived without a documented unit to label it with.
        reason = (
            "normalized values have no documented unit"
            if value_semantics_status == VALUE_SEMANTICS_NORMALIZED
            else f"value semantics are {value_semantics_status}"
        )
        return {
            "status": RANK_STATUS_NORMALIZATION_UNRESOLVED,
            "records": [],
            "excluded": [],
            "presentation": presentation,
            "reason": reason,
        }

    if scan.target_matches != 1:
        return {
            "status": "target_gene_match_failed",
            "records": [],
            "excluded": [],
            "presentation": presentation,
            "reason": f"expected exactly one exact gene row, found {scan.target_matches}",
        }

    duplicates = set(scan.duplicate_population_labels)

    records: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for idx, label in enumerate(scan.population_labels):
        target_value = scan.target_row[idx] if idx < len(scan.target_row) else None
        total = scan.column_totals[idx] if idx < len(scan.column_totals) else None

        if not label:
            excluded.append({"index": idx, "label": label, "reason": "empty_population_label"})
            continue
        if label in duplicates:
            excluded.append({"index": idx, "label": label, "reason": "duplicate_population_label"})
            continue
        if target_value is None:
            excluded.append({"index": idx, "label": label, "reason": "malformed_target_value"})
            continue

        if counts_mode:
            if total is None or not math.isfinite(total):
                excluded.append({"index": idx, "label": label, "reason": "population_total_missing"})
                continue
            if total <= 0:
                excluded.append(
                    {"index": idx, "label": label, "reason": "population_total_not_positive"}
                )
                continue
            per_100k = target_value / total * 100_000.0
            records.append(
                {
                    "population_label": label,
                    "target_source_value": target_value,
                    "ranking_value": per_100k,
                    "ranking_metric": presentation["ranking_metric"],
                    "expression_unit": presentation["expression_unit"],
                    "population_total": total,
                    "transcripts_per_100k": per_100k,
                }
            )
        else:
            # Documented normalized scale: rank on the source value as published
            # and leave transcripts_per_100k unpopulated.
            records.append(
                {
                    "population_label": label,
                    "target_source_value": target_value,
                    "ranking_value": target_value,
                    "ranking_metric": presentation["ranking_metric"],
                    "expression_unit": presentation["expression_unit"],
                    "population_total": None,
                    "transcripts_per_100k": None,
                }
            )

    if not records:
        return {
            "status": "no_valid_populations",
            "records": [],
            "excluded": excluded,
            "presentation": presentation,
            "reason": "every population column was excluded",
        }

    return {
        "status": RANK_STATUS_SUCCESS,
        "records": records,
        "excluded": excluded,
        "presentation": presentation,
        "reason": None,
    }


def rank_records(records: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return records sorted by descending ranking value (stable on label)."""
    return sorted(
        records,
        key=lambda r: (-float(r["ranking_value"]), str(r["population_label"])),
    )


def summarize_ranking(
    records: Sequence[dict[str, Any]],
    *,
    value_semantics_status: str,
    top_n: int = TOP_N_DEFAULT,
) -> dict[str, Any]:
    """Objective breadth / concentration metrics over valid populations."""
    ranked = rank_records(records)
    values = [float(r["ranking_value"]) for r in ranked]
    valid_count = len(ranked)
    nonzero_count = sum(1 for v in values if v > 0)
    nonzero_percentage = (nonzero_count / valid_count * 100.0) if valid_count else 0.0
    top = ranked[:top_n]

    total_ranking = sum(values)
    top_ranking = sum(float(r["ranking_value"]) for r in top)
    normalized_share = (top_ranking / total_ranking * 100.0) if total_ranking > 0 else None

    raw_share: float | None = None
    total_raw = sum(float(r["target_source_value"]) for r in ranked)
    if value_semantics_status in COUNT_COMPATIBLE_SEMANTICS and total_raw > 0:
        top_raw = sum(float(r["target_source_value"]) for r in top)
        raw_share = top_raw / total_raw * 100.0

    maximum = max(values) if values else None
    median = statistics.median(values) if values else None
    mean = statistics.fmean(values) if values else None
    max_to_median = (maximum / median) if (maximum is not None and median) else None

    return {
        "calculation_version": CALCULATION_VERSION,
        "value_semantics_status": value_semantics_status,
        "ranking_metric": ranked[0]["ranking_metric"] if ranked else None,
        "expression_unit": ranked[0]["expression_unit"] if ranked else None,
        "valid_population_count": valid_count,
        "nonzero_population_count": nonzero_count,
        "nonzero_percentage": round(nonzero_percentage, 2),
        "maximum": maximum,
        "median": median,
        "mean": mean,
        "max_to_median_ratio": max_to_median,
        "top_populations": top,
        "bottom_populations": ranked[-top_n:] if valid_count > top_n else [],
        "top_10_normalized_expression_share": (
            round(normalized_share, 2) if normalized_share is not None else None
        ),
        "top_10_normalized_expression_share_description": (
            "Share of the normalized population-level expression signal represented "
            "by the ten highest-ranked populations."
        ),
        "top_10_raw_target_count_share": (round(raw_share, 2) if raw_share is not None else None),
        "top_10_raw_target_count_share_description": (
            "Share of the raw target-gene counts summed across valid populations."
            if raw_share is not None
            else None
        ),
        "total_target_source_value": total_raw,
    }


# --------------------------------------------------------------------------------------
# Population labels
# --------------------------------------------------------------------------------------
def parse_population_label(label: str) -> dict[str, Any]:
    """Extract only what the GEO label explicitly encodes.

    DropViz metacell columns look like ``FC_1-1`` / ``FC_1_1`` (region, cluster,
    subcluster). Broad cell class is NOT inferred from biological knowledge.
    """
    text = str(label or "").strip()
    out: dict[str, Any] = {
        "display_label": text,
        "region_abbreviation": None,
        "cluster_number": None,
        "subcluster_number": None,
        "broad_cell_class": None,
        "neuronal_or_glial": None,
    }
    if not text:
        return out

    region, sep, remainder = text.partition("_")
    if sep and region:
        out["region_abbreviation"] = region
        parts = remainder.replace("-", "_").split("_")
        if parts and parts[0].isdigit():
            out["cluster_number"] = int(parts[0])
        if len(parts) > 1 and parts[1].isdigit():
            out["subcluster_number"] = int(parts[1])
    return out


def label_mapping_status(parsed: Sequence[dict[str, Any]]) -> str:
    """Roll parsed labels up into resolved / partial / unresolved."""
    if not parsed:
        return LABEL_MAPPING_UNRESOLVED
    with_region = sum(1 for p in parsed if p.get("region_abbreviation"))
    with_cluster = sum(1 for p in parsed if p.get("cluster_number") is not None)
    if with_region == len(parsed) and with_cluster == len(parsed):
        # Region and cluster are resolved; broad class is never inferred, so a
        # fully resolved label set still stops short of cell-class mapping.
        return LABEL_MAPPING_PARTIAL
    if with_region or with_cluster:
        return LABEL_MAPPING_PARTIAL
    return LABEL_MAPPING_UNRESOLVED


# --------------------------------------------------------------------------------------
# Chart
# --------------------------------------------------------------------------------------
def render_top_populations_png(
    *,
    mouse_gene_symbol: str,
    top_populations: Sequence[dict[str, Any]],
    axis_label: str,
    confidence_interval_status: str = CONFIDENCE_INTERVAL_METHOD_UNRESOLVED,
    dpi: int = 180,
) -> bytes:
    """Deterministic point chart of the highest-expressing populations.

    Confidence whiskers are drawn only when a validated interval method supplied
    ``lower``/``upper`` on every plotted record.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = list(top_populations)
    labels = [str(r.get("population_label", "")) for r in rows]
    values = [float(r.get("ranking_value", 0.0)) for r in rows]

    draw_ci = confidence_interval_status != CONFIDENCE_INTERVAL_METHOD_UNRESOLVED and all(
        r.get("lower") is not None and r.get("upper") is not None for r in rows
    )

    height = max(2.6, 0.42 * len(rows) + 1.7)
    fig, ax = plt.subplots(figsize=(7.6, height), dpi=dpi)
    positions = list(range(len(rows) - 1, -1, -1))

    if draw_ci:
        for pos, row in zip(positions, rows):
            ax.plot(
                [float(row["lower"]), float(row["upper"])],
                [pos, pos],
                color="#333333",
                linewidth=1.0,
                solid_capstyle="butt",
            )
    ax.scatter(values, positions, color="#111111", s=18, zorder=3)

    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(axis_label, fontsize=9)
    ax.set_title(f"{CHART_TITLE}\n{mouse_gene_symbol}", fontsize=11, linespacing=1.35)
    ax.set_xlim(left=0)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.text(0.01, 0.01, CHART_SOURCE_NOTE, fontsize=7, color="#555555")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.18, facecolor="white")
    plt.close(fig)
    return buf.getvalue()


__all__ = [
    "ALLOWED_HOSTS",
    "CALCULATION_VERSION",
    "CONFIDENCE_INTERVAL_METHOD_UNRESOLVED",
    "COUNT_COMPATIBLE_SEMANTICS",
    "GEO_ACCESSION",
    "LABEL_MAPPING_PARTIAL",
    "LABEL_MAPPING_RESOLVED",
    "LABEL_MAPPING_UNRESOLVED",
    "METACELL_FILENAME",
    "MatrixScan",
    "RANK_STATUS_NORMALIZATION_UNRESOLVED",
    "RANK_STATUS_SUCCESS",
    "SOURCE_NAME",
    "VALUE_SEMANTICS_COUNT_COMPATIBLE",
    "VALUE_SEMANTICS_NORMALIZED",
    "VALUE_SEMANTICS_RAW_COUNTS",
    "VALUE_SEMANTICS_TRANSFORMED",
    "VALUE_SEMANTICS_UNRESOLVED",
    "build_matrix_profile",
    "build_ranking_records",
    "candidate_download_urls",
    "classify_value_semantics",
    "download_cgi_url",
    "download_metacell_matrix",
    "is_allowlisted_url",
    "label_mapping_status",
    "normalize_gene_symbol",
    "open_matrix_stream",
    "parse_population_label",
    "rank_records",
    "ranking_presentation",
    "render_top_populations_png",
    "scan_matrix",
    "summarize_ranking",
    "supplementary_url",
    "validate_gzip_payload",
]
