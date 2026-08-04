"""Allen Human Brain Atlas RNA-seq donor ZIP client.

Downloads well-known donor packages and parses TPM matrices for a gene.
Does **not** normalize into evidence records — that belongs in Section 2b.

Key endpoints::

    GET https://human.brain-map.org/api/v2/well_known_file_download/{id}

Validated donor well-known file IDs: ``278447594``, ``278448166``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`
or structured error dicts from pure helpers.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
import statistics
import zipfile
from typing import Any

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Allen Brain Atlas"
HUMAN_BRAIN_MAP_BASE = "https://human.brain-map.org/api/v2"
WELL_KNOWN_FILE_DOWNLOAD = f"{HUMAN_BRAIN_MAP_BASE}/well_known_file_download"

# Locked Section 2b donor packages (two distinct brains).
DEFAULT_DONOR_WELL_KNOWN_FILE_IDS: tuple[int, ...] = (278447594, 278448166)

REQUIRED_ZIP_MEMBERS: tuple[str, ...] = (
    "Contents.txt",
    "Genes.csv",
    "SampleAnnot.csv",
    "RNAseqTPM.csv",
    "RNAseqCounts.csv",
    "Ontology.csv",
)
# Persist as derived artifacts parented to the raw ZIP (not as fake HTTP).
PERSISTED_DERIVED_MEMBERS: tuple[str, ...] = (
    "Contents.txt",
    "Genes.csv",
    "SampleAnnot.csv",
    "RNAseqTPM.csv",
)
# Presence-validated only; do not duplicate full contents unless used.
VALIDATED_ONLY_MEMBERS: tuple[str, ...] = ("RNAseqCounts.csv", "Ontology.csv")

ZIP_MAGIC = b"PK"
REQUEST_HEADERS = {
    "User-Agent": "GeneDossier/0.1.0 (research; provenance-first gene dossier client)",
    "Accept": "application/zip,application/octet-stream,*/*;q=0.8",
}

# Scale-detection tolerances (persisted in audit exactly as used).
SCALE_FRACTION_OF_MILLION_TARGET = 1.0
SCALE_CONVENTIONAL_TPM_TARGET = 1_000_000.0
SCALE_REL_TOL = 0.05
SCALE_ABS_TOL_FRACTION = 0.05
SCALE_ABS_TOL_CONVENTIONAL = 50_000.0

RECOGNIZED_REPLICATE_SAMPLE_VALUES = frozenset({"yes", "no"})
BIOLOGICAL_REPLICATE_VALUE = "no"

SAFE_RESPONSE_HEADER_KEYS = ("content-type", "content-length", "server")


def download_url_for_donor(well_known_file_id: int | str) -> str:
    """Return the Allen well-known-file download URL for one donor package."""
    return f"{WELL_KNOWN_FILE_DOWNLOAD}/{int(well_known_file_id)}"


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
    return (
        head.startswith(b"<!doctype html")
        or head.startswith(b"<html")
        or b"<html" in head[:200]
    )


def download_donor_zip(
    well_known_file_id: int | str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Download one Allen human RNA-seq donor ZIP (never raises).

    On success ``data`` includes ``content`` (bytes), ``byte_size``, ``sha256``,
    and ``well_known_file_id``. The section node persists ``content`` as the
    immutable raw HTTP artifact.
    """
    cfg = settings or get_settings()
    try:
        donor_id = int(well_known_file_id)
    except (TypeError, ValueError):
        return _tool_result(
            endpoint_name="well_known_file_download",
            gene_symbol=gene_symbol,
            request_url=WELL_KNOWN_FILE_DOWNLOAD,
            request_params={"well_known_file_id": well_known_file_id},
            success=False,
            error_type="invalid_request",
            error_message=f"invalid well_known_file_id: {well_known_file_id!r}",
        )

    url = download_url_for_donor(donor_id)
    request_params = {
        "well_known_file_id": donor_id,
        "url": url,
        "request_headers": dict(REQUEST_HEADERS),
        "follow_redirects": True,
    }
    try:
        with httpx.Client(
            timeout=cfg.http_timeout_seconds,
            follow_redirects=True,
        ) as client:
            response = client.get(url, headers=REQUEST_HEADERS)
        content = response.content or b""
        content_type = response.headers.get("content-type")
        response_headers = _safe_response_headers(response.headers)
        meta = {
            "well_known_file_id": donor_id,
            "byte_size": len(content),
            "sha256": hashlib.sha256(content).hexdigest() if content else None,
            "content_type": content_type,
            "response_headers": response_headers,
        }
        if not response.is_success:
            return _tool_result(
                endpoint_name="well_known_file_download",
                gene_symbol=gene_symbol or str(donor_id),
                request_url=str(response.url) if response.url else url,
                request_params=request_params,
                success=False,
                status_code=response.status_code,
                data={**meta, "raw_text_preview": content[:400].decode("utf-8", "replace")},
                error_type="http_error",
                error_message=f"HTTP {response.status_code}",
            )
        if not content.startswith(ZIP_MAGIC):
            err = "html_masquerading_as_zip" if _looks_like_html(content) else "invalid_zip_magic"
            return _tool_result(
                endpoint_name="well_known_file_download",
                gene_symbol=gene_symbol or str(donor_id),
                request_url=str(response.url) if response.url else url,
                request_params=request_params,
                success=False,
                status_code=response.status_code,
                data={
                    **meta,
                    "raw_text_preview": content[:400].decode("utf-8", "replace"),
                },
                error_type=err,
                error_message=(
                    "Donor download looked like HTML, not a ZIP"
                    if err == "html_masquerading_as_zip"
                    else "Donor download missing ZIP magic bytes"
                ),
            )
        if _looks_like_html(content):
            return _tool_result(
                endpoint_name="well_known_file_download",
                gene_symbol=gene_symbol or str(donor_id),
                request_url=str(response.url) if response.url else url,
                request_params=request_params,
                success=False,
                status_code=response.status_code,
                data={
                    **meta,
                    "raw_text_preview": content[:400].decode("utf-8", "replace"),
                },
                error_type="html_masquerading_as_zip",
                error_message="Donor download looked like HTML, not a ZIP",
            )
        return _tool_result(
            endpoint_name="well_known_file_download",
            gene_symbol=gene_symbol or str(donor_id),
            request_url=str(response.url) if response.url else url,
            request_params=request_params,
            success=True,
            status_code=response.status_code,
            data={
                **meta,
                "content": content,
            },
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name="well_known_file_download",
            gene_symbol=gene_symbol or str(donor_id),
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="well_known_file_download",
            gene_symbol=gene_symbol or str(donor_id),
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="well_known_file_download",
            gene_symbol=gene_symbol or str(donor_id),
            request_url=url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def unpack_zip_members(zip_bytes: bytes) -> dict[str, Any]:
    """Unpack required members by basename; reject missing or duplicate names.

    Returns ``{"success": True, "members": {basename: bytes}, ...}`` or a
    structured failure dict (never raises).
    """
    if not zip_bytes or not zip_bytes.startswith(ZIP_MAGIC):
        return {
            "success": False,
            "error_type": "invalid_zip_magic",
            "error_message": "ZIP magic bytes missing",
            "members": {},
        }
    if _looks_like_html(zip_bytes):
        return {
            "success": False,
            "error_type": "html_masquerading_as_zip",
            "error_message": "Payload looks like HTML, not a ZIP",
            "members": {},
        }
    try:
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
            basenames: dict[str, list[str]] = {}
            for info in zf.infolist():
                if info.is_dir():
                    continue
                name = info.filename.replace("\\", "/")
                base = name.rsplit("/", 1)[-1]
                if not base:
                    continue
                basenames.setdefault(base, []).append(name)

            members: dict[str, bytes] = {}
            diagnostics: list[dict[str, Any]] = []
            for required in REQUIRED_ZIP_MEMBERS:
                paths = basenames.get(required) or []
                if not paths:
                    return {
                        "success": False,
                        "error_type": "missing_zip_member",
                        "error_message": f"Missing required ZIP member {required!r}",
                        "members": {},
                        "basenames_seen": sorted(basenames),
                    }
                if len(paths) > 1:
                    return {
                        "success": False,
                        "error_type": "duplicate_zip_member",
                        "error_message": (
                            f"Ambiguous duplicate ZIP member {required!r}: {paths}"
                        ),
                        "members": {},
                        "duplicate_paths": paths,
                    }
                members[required] = zf.read(paths[0])
                diagnostics.append(
                    {
                        "basename": required,
                        "zip_path": paths[0],
                        "byte_size": len(members[required]),
                    }
                )
            return {
                "success": True,
                "members": members,
                "member_diagnostics": diagnostics,
                "basenames_seen": sorted(basenames),
            }
    except zipfile.BadZipFile as exc:
        return {
            "success": False,
            "error_type": "bad_zip",
            "error_message": str(exc),
            "members": {},
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "success": False,
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "members": {},
        }


def _parse_csv_text(raw: str | bytes) -> list[dict[str, str]]:
    if isinstance(raw, bytes):
        text = raw.decode("utf-8-sig", errors="replace")
    else:
        text = (raw or "").lstrip("\ufeff")
    text = text.strip()
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
        if any(cleaned.values()):
            rows.append(cleaned)
    return rows


def _normalize_entrez(value: Any) -> str | None:
    if value is None or isinstance(value, bool):
        return None
    text = str(value).strip()
    if not text or text.upper() == "NA" or not text.isdigit():
        return None
    return str(int(text))


def scale_detection_tolerances() -> dict[str, float]:
    """Exact tolerances used by :func:`detect_tpm_scale` (for audit persistence)."""
    return {
        "fraction_of_million_target": SCALE_FRACTION_OF_MILLION_TARGET,
        "conventional_tpm_target": SCALE_CONVENTIONAL_TPM_TARGET,
        "relative_tolerance": SCALE_REL_TOL,
        "absolute_tolerance_fraction": SCALE_ABS_TOL_FRACTION,
        "absolute_tolerance_conventional": SCALE_ABS_TOL_CONVENTIONAL,
    }


def detect_tpm_scale(column_sums: list[float]) -> dict[str, Any]:
    """Classify TPM matrix scale from per-sample column sums.

    - median ≈ 1 → ``fraction_of_million`` (multiply values by 1e6)
    - median ≈ 1e6 → ``conventional_tpm`` (identity)
    - otherwise fail
    """
    tolerances = scale_detection_tolerances()
    finite = [float(s) for s in column_sums if math.isfinite(float(s))]
    if not finite:
        return {
            "success": False,
            "scale_mode": None,
            "scale_multiplier": None,
            "median_column_sum": None,
            "column_sums": list(column_sums),
            "tolerances": tolerances,
            "error_type": "empty_column_sums",
            "error_message": "No finite sample-column sums for scale detection",
        }
    median_sum = float(statistics.median(finite))
    if math.isclose(
        median_sum,
        SCALE_FRACTION_OF_MILLION_TARGET,
        rel_tol=SCALE_REL_TOL,
        abs_tol=SCALE_ABS_TOL_FRACTION,
    ):
        return {
            "success": True,
            "scale_mode": "fraction_of_million",
            "scale_multiplier": 1_000_000.0,
            "median_column_sum": median_sum,
            "column_sums": [float(s) for s in column_sums],
            "tolerances": tolerances,
        }
    if math.isclose(
        median_sum,
        SCALE_CONVENTIONAL_TPM_TARGET,
        rel_tol=SCALE_REL_TOL,
        abs_tol=SCALE_ABS_TOL_CONVENTIONAL,
    ):
        return {
            "success": True,
            "scale_mode": "conventional_tpm",
            "scale_multiplier": 1.0,
            "median_column_sum": median_sum,
            "column_sums": [float(s) for s in column_sums],
            "tolerances": tolerances,
        }
    return {
        "success": False,
        "scale_mode": None,
        "scale_multiplier": None,
        "median_column_sum": median_sum,
        "column_sums": [float(s) for s in column_sums],
        "tolerances": tolerances,
        "error_type": "unrecognized_tpm_scale",
        "error_message": (
            f"Median sample-column sum {median_sum} is neither ~1 "
            f"(fraction_of_million) nor ~1e6 (conventional_tpm)"
        ),
    }


def _brain_identity(sample_rows: list[dict[str, str]]) -> dict[str, Any]:
    identities: list[str] = []
    for row in sample_rows:
        brain = str(row.get("brain") or "").strip()
        if brain:
            identities.append(brain)
    distinct = sorted(set(identities))
    if len(distinct) != 1:
        return {
            "success": False,
            "brain_identity": None,
            "distinct_brains": distinct,
            "error_type": "brain_identity_invalid",
            "error_message": (
                "Expected exactly one nonempty brain identity per package; "
                f"got {distinct!r}"
            ),
        }
    return {
        "success": True,
        "brain_identity": distinct[0],
        "distinct_brains": distinct,
    }


def filter_biological_samples(
    sample_rows: list[dict[str, str]],
    *,
    exclude_technical_replicates: bool = True,
) -> dict[str, Any]:
    """Keep biological samples; reject unrecognized replicate_sample values."""
    retained: list[dict[str, str]] = []
    technical: list[dict[str, str]] = []
    unknown: list[dict[str, str]] = []
    for row in sample_rows:
        raw = str(row.get("replicate_sample") or "").strip()
        key = raw.lower()
        if key not in RECOGNIZED_REPLICATE_SAMPLE_VALUES:
            unknown.append(row)
            continue
        if exclude_technical_replicates and key != BIOLOGICAL_REPLICATE_VALUE:
            technical.append(row)
            continue
        retained.append(row)

    if unknown:
        return {
            "success": False,
            "retained_samples": [],
            "technical_replicate_count": len(technical),
            "unknown_replicate_count": len(unknown),
            "unknown_replicate_values": sorted(
                {
                    str(r.get("replicate_sample") or "").strip()
                    for r in unknown
                }
            ),
            "error_type": "unrecognized_replicate_sample",
            "error_message": (
                "Unrecognized replicate_sample values must not be treated as "
                f"biological: {[str(r.get('replicate_sample')) for r in unknown[:5]]}"
            ),
        }
    if not retained:
        return {
            "success": False,
            "retained_samples": [],
            "technical_replicate_count": len(technical),
            "unknown_replicate_count": 0,
            "error_type": "no_biological_samples",
            "error_message": "No biological samples retained after replicate filtering",
        }

    well_ids = [str(r.get("well_id") or "").strip() for r in retained]
    if any(not w for w in well_ids):
        return {
            "success": False,
            "retained_samples": [],
            "technical_replicate_count": len(technical),
            "unknown_replicate_count": 0,
            "error_type": "missing_well_id",
            "error_message": "Retained sample missing well_id",
        }
    if len(well_ids) != len(set(well_ids)):
        return {
            "success": False,
            "retained_samples": [],
            "technical_replicate_count": len(technical),
            "unknown_replicate_count": 0,
            "error_type": "duplicate_well_id",
            "error_message": "Duplicate well_id among retained biological samples",
        }
    return {
        "success": True,
        "retained_samples": retained,
        "technical_replicate_count": len(technical),
        "unknown_replicate_count": 0,
        "retained_well_ids": well_ids,
        "retained_sample_names": [
            str(r.get("RNAseq_sample_name") or r.get("sample_name") or "").strip()
            for r in retained
        ],
    }


def _tpm_matrix_has_header(first_cell: str) -> bool:
    """Allen donor TPM matrices have no header; first cell is a gene symbol."""
    text = (first_cell or "").strip().strip('"').lower()
    return text in {"gene_symbol", "gene", "symbol", "gene_id"}


def compute_sample_column_sums(
    tpm_bytes: bytes,
    *,
    sample_column_count: int,
) -> dict[str, Any]:
    """Sum each sample column across all gene rows (scale-detection input).

    Allen Human Brain Atlas RNAseqTPM.csv has **no header row**: column 0 is the
    gene symbol and columns 1..N align positionally with SampleAnnot row order.
    """
    text = tpm_bytes.decode("utf-8-sig", errors="replace")
    reader = csv.reader(io.StringIO(text))
    rows = [row for row in reader if row]
    if not rows:
        return {
            "success": False,
            "error_type": "empty_tpm_matrix",
            "error_message": "RNAseqTPM.csv is empty",
            "column_sums": [],
        }
    start_idx = 0
    if _tpm_matrix_has_header(rows[0][0] if rows[0] else ""):
        start_idx = 1
        n_cols = len(rows[0]) - 1
    else:
        n_cols = len(rows[0]) - 1
    if n_cols != sample_column_count:
        return {
            "success": False,
            "error_type": "tpm_sample_count_mismatch",
            "error_message": (
                f"RNAseqTPM sample columns ({n_cols}) != SampleAnnot rows "
                f"({sample_column_count})"
            ),
            "column_sums": [],
            "tpm_sample_column_count": n_cols,
            "sample_annot_row_count": sample_column_count,
        }
    sums = [0.0] * n_cols
    gene_rows = 0
    for row in rows[start_idx:]:
        if not row:
            continue
        if len(row) != n_cols + 1:
            return {
                "success": False,
                "error_type": "tpm_row_width_mismatch",
                "error_message": (
                    f"TPM row width {len(row)} != expected {n_cols + 1}"
                ),
                "column_sums": [],
            }
        gene_rows += 1
        for i in range(n_cols):
            try:
                value = float(row[i + 1])
            except (TypeError, ValueError):
                return {
                    "success": False,
                    "error_type": "non_numeric_tpm",
                    "error_message": f"Non-numeric TPM at gene row {gene_rows} col {i}",
                    "column_sums": [],
                }
            if not math.isfinite(value) or value < 0:
                return {
                    "success": False,
                    "error_type": "invalid_tpm_value",
                    "error_message": (
                        f"Non-finite or negative TPM at gene row {gene_rows} col {i}"
                    ),
                    "column_sums": [],
                }
            sums[i] += value
    return {
        "success": True,
        "column_sums": sums,
        "gene_row_count": gene_rows,
        "has_header_row": start_idx == 1,
        "column_alignment": "sample_annot_row_order",
    }


def find_gene_tpm_row(
    tpm_bytes: bytes,
    *,
    gene_symbol: str,
    entrez_gene_id: str | None,
    genes_rows: list[dict[str, str]],
) -> dict[str, Any]:
    """Locate the single TPM matrix row for the gene (symbol + Entrez gate)."""
    expected_entrez = _normalize_entrez(entrez_gene_id)
    symbol_target = (gene_symbol or "").strip().lower()
    if not symbol_target:
        return {
            "success": False,
            "error_type": "invalid_request",
            "error_message": "gene_symbol is required",
        }

    gene_matches = [
        row
        for row in genes_rows
        if str(row.get("gene_symbol") or "").strip().lower() == symbol_target
    ]
    if expected_entrez is not None:
        entrez_matches = [
            row
            for row in gene_matches
            if _normalize_entrez(row.get("entrez_id")) == expected_entrez
        ]
        if not entrez_matches:
            entrez_matches = [
                row
                for row in genes_rows
                if _normalize_entrez(row.get("entrez_id")) == expected_entrez
            ]
        gene_matches = entrez_matches

    if len(gene_matches) != 1:
        return {
            "success": False,
            "error_type": "gene_identity_mismatch",
            "error_message": (
                f"Expected exactly one Genes.csv match for {gene_symbol!r} "
                f"(entrez={expected_entrez}); got {len(gene_matches)}"
            ),
            "match_count": len(gene_matches),
        }
    gene_meta = gene_matches[0]
    matrix_symbol = str(gene_meta.get("gene_symbol") or "").strip()

    text = tpm_bytes.decode("utf-8-sig", errors="replace")
    reader = csv.reader(io.StringIO(text))
    rows = [row for row in reader if row]
    if not rows:
        return {
            "success": False,
            "error_type": "empty_tpm_matrix",
            "error_message": "RNAseqTPM.csv is empty",
        }
    start_idx = 1 if _tpm_matrix_has_header(rows[0][0] if rows[0] else "") else 0

    matched_rows: list[list[str]] = []
    for row in rows[start_idx:]:
        if not row:
            continue
        if str(row[0]).strip().strip('"').lower() == matrix_symbol.lower():
            matched_rows.append(row)
    if len(matched_rows) != 1:
        return {
            "success": False,
            "error_type": "tpm_gene_row_count",
            "error_message": (
                f"Expected exactly one TPM row for {matrix_symbol!r}; "
                f"got {len(matched_rows)}"
            ),
            "match_count": len(matched_rows),
        }
    row = matched_rows[0]
    values: list[float] = []
    for raw in row[1:]:
        try:
            number = float(raw)
        except (TypeError, ValueError):
            return {
                "success": False,
                "error_type": "non_numeric_tpm",
                "error_message": "Selected TPM row contains non-numeric values",
            }
        if not math.isfinite(number) or number < 0:
            return {
                "success": False,
                "error_type": "invalid_tpm_value",
                "error_message": "Selected TPM values must be finite and nonnegative",
            }
        values.append(number)
    return {
        "success": True,
        "gene_symbol": matrix_symbol,
        "entrez_id": _normalize_entrez(gene_meta.get("entrez_id")),
        "gene_metadata": gene_meta,
        "raw_tpm_values": values,
        "column_alignment": "sample_annot_row_order",
    }


def parse_donor_package(
    zip_bytes: bytes,
    *,
    gene_symbol: str,
    entrez_gene_id: str | None,
    well_known_file_id: int | str | None = None,
    exclude_technical_replicates: bool = True,
) -> dict[str, Any]:
    """Unpack, validate, scale-detect, and extract gene TPM for one donor ZIP."""
    unpacked = unpack_zip_members(zip_bytes)
    if not unpacked.get("success"):
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            "error_type": unpacked.get("error_type"),
            "error_message": unpacked.get("error_message"),
            "unpack": unpacked,
        }
    members: dict[str, bytes] = unpacked["members"]
    sample_rows = _parse_csv_text(members["SampleAnnot.csv"])
    genes_rows = _parse_csv_text(members["Genes.csv"])
    if not sample_rows:
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            "error_type": "empty_sample_annot",
            "error_message": "SampleAnnot.csv has no rows",
        }

    brain = _brain_identity(sample_rows)
    if not brain.get("success"):
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            **brain,
        }

    filtered = filter_biological_samples(
        sample_rows,
        exclude_technical_replicates=exclude_technical_replicates,
    )
    if not filtered.get("success"):
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            "brain_identity": brain.get("brain_identity"),
            **filtered,
        }

    sums_result = compute_sample_column_sums(
        members["RNAseqTPM.csv"],
        sample_column_count=len(sample_rows),
    )
    if not sums_result.get("success"):
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            "brain_identity": brain.get("brain_identity"),
            **sums_result,
        }

    scale = detect_tpm_scale(list(sums_result["column_sums"]))
    if not scale.get("success"):
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            "brain_identity": brain.get("brain_identity"),
            "scale_detection": scale,
            "sample_column_sum_diagnostics": {
                "column_sums": sums_result["column_sums"],
                "gene_row_count": sums_result.get("gene_row_count"),
                "sample_column_names": sums_result.get("sample_column_names"),
            },
            "error_type": scale.get("error_type"),
            "error_message": scale.get("error_message"),
        }

    gene_row = find_gene_tpm_row(
        members["RNAseqTPM.csv"],
        gene_symbol=gene_symbol,
        entrez_gene_id=entrez_gene_id,
        genes_rows=genes_rows,
    )
    if not gene_row.get("success"):
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            "brain_identity": brain.get("brain_identity"),
            "scale_detection": scale,
            "sample_column_sum_diagnostics": {
                "column_sums": sums_result["column_sums"],
                "gene_row_count": sums_result.get("gene_row_count"),
                "sample_column_names": sums_result.get("sample_column_names"),
                "tolerances": scale.get("tolerances"),
            },
            **gene_row,
        }

    # Map retained samples to TPM columns by SampleAnnot row order (no TPM header).
    retained = list(filtered["retained_samples"])
    retained_indices: list[int] = []
    for sample in retained:
        try:
            idx = sample_rows.index(sample)
        except ValueError:
            return {
                "success": False,
                "well_known_file_id": well_known_file_id,
                "brain_identity": brain.get("brain_identity"),
                "error_type": "sample_index_missing",
                "error_message": "Retained sample not found in SampleAnnot rows",
                "scale_detection": scale,
            }
        retained_indices.append(idx)

    selected_raw: list[float] = []
    selected_meta: list[dict[str, Any]] = []
    multiplier = float(scale["scale_multiplier"])
    tpm_values = list(gene_row["raw_tpm_values"])
    if len(tpm_values) != len(sample_rows):
        return {
            "success": False,
            "well_known_file_id": well_known_file_id,
            "brain_identity": brain.get("brain_identity"),
            "error_type": "tpm_sample_count_mismatch",
            "error_message": (
                f"Gene TPM columns ({len(tpm_values)}) != SampleAnnot rows "
                f"({len(sample_rows)})"
            ),
            "scale_detection": scale,
        }
    for sample, idx in zip(retained, retained_indices):
        raw_val = float(tpm_values[idx])
        scaled = raw_val * multiplier
        if not math.isfinite(scaled) or scaled < 0:
            return {
                "success": False,
                "well_known_file_id": well_known_file_id,
                "brain_identity": brain.get("brain_identity"),
                "error_type": "invalid_selected_tpm",
                "error_message": "Selected TPM after scaling must be finite and nonnegative",
                "scale_detection": scale,
            }
        selected_raw.append(scaled)
        selected_meta.append(
            {
                "RNAseq_sample_name": str(
                    sample.get("RNAseq_sample_name") or ""
                ).strip(),
                "well_id": str(sample.get("well_id") or "").strip(),
                "replicate_sample": str(sample.get("replicate_sample") or "").strip(),
                "sample_annot_index": idx,
                "tpm": scaled,
                "raw_matrix_value": raw_val,
            }
        )

    mean_tpm = float(statistics.mean(selected_raw)) if selected_raw else None
    derived_members = {
        name: members[name] for name in PERSISTED_DERIVED_MEMBERS if name in members
    }
    return {
        "success": True,
        "well_known_file_id": (
            int(well_known_file_id) if well_known_file_id is not None else None
        ),
        "brain_identity": brain["brain_identity"],
        "sample_annot_row_count": len(sample_rows),
        "retained_sample_count": len(retained),
        "technical_replicate_count": filtered.get("technical_replicate_count"),
        "retained_well_ids": filtered.get("retained_well_ids"),
        "scale_detection": scale,
        "sample_column_sum_diagnostics": {
            "column_sums": sums_result["column_sums"],
            "gene_row_count": sums_result.get("gene_row_count"),
            "has_header_row": sums_result.get("has_header_row"),
            "column_alignment": sums_result.get("column_alignment"),
            "tolerances": scale.get("tolerances"),
            "median_column_sum": scale.get("median_column_sum"),
            "scale_mode": scale.get("scale_mode"),
            "scale_multiplier": scale.get("scale_multiplier"),
        },
        "gene_symbol": gene_row.get("gene_symbol"),
        "entrez_id": gene_row.get("entrez_id"),
        "selected_sample_tpms": selected_meta,
        "mean_tpm": mean_tpm,
        "unit": "TPM",
        "derived_members": derived_members,
        "validated_only_members_present": {
            name: name in members for name in VALIDATED_ONLY_MEMBERS
        },
        "unpack": {
            "member_diagnostics": unpacked.get("member_diagnostics"),
            "basenames_seen": unpacked.get("basenames_seen"),
        },
    }


def pooled_mean_across_donors(
    donor_results: list[dict[str, Any]],
) -> dict[str, Any]:
    """Pooled mean over all retained biological samples (not equal-weight donors)."""
    if len(donor_results) != 2:
        return {
            "success": False,
            "error_type": "donor_count",
            "error_message": "Polished human RNA-seq TPM requires exactly two donor packages",
        }
    if not all(r.get("success") for r in donor_results):
        return {
            "success": False,
            "error_type": "incomplete_donors",
            "error_message": "One or more donor packages failed validation",
            "donor_success": [bool(r.get("success")) for r in donor_results],
        }
    brains = [str(r.get("brain_identity") or "") for r in donor_results]
    if any(not b for b in brains) or len(set(brains)) != 2:
        return {
            "success": False,
            "error_type": "brain_identity_not_distinct",
            "error_message": f"Donor brain identities must be distinct; got {brains!r}",
            "brain_identities": brains,
        }
    all_tpms: list[float] = []
    for result in donor_results:
        for sample in result.get("selected_sample_tpms") or []:
            value = float(sample["tpm"])
            if not math.isfinite(value) or value < 0:
                return {
                    "success": False,
                    "error_type": "invalid_selected_tpm",
                    "error_message": "Pooled TPM values must be finite and nonnegative",
                }
            all_tpms.append(value)
    if not all_tpms:
        return {
            "success": False,
            "error_type": "no_biological_samples",
            "error_message": "No retained biological TPM values to pool",
        }
    return {
        "success": True,
        "mean_tpm": float(statistics.mean(all_tpms)),
        "retained_sample_count": len(all_tpms),
        "per_donor_retained_counts": [
            int(r.get("retained_sample_count") or 0) for r in donor_results
        ],
        "brain_identities": brains,
        "unit": "TPM",
    }


__all__ = [
    "SOURCE_NAME",
    "HUMAN_BRAIN_MAP_BASE",
    "WELL_KNOWN_FILE_DOWNLOAD",
    "DEFAULT_DONOR_WELL_KNOWN_FILE_IDS",
    "REQUIRED_ZIP_MEMBERS",
    "PERSISTED_DERIVED_MEMBERS",
    "VALIDATED_ONLY_MEMBERS",
    "REQUEST_HEADERS",
    "SCALE_FRACTION_OF_MILLION_TARGET",
    "SCALE_CONVENTIONAL_TPM_TARGET",
    "SCALE_REL_TOL",
    "SCALE_ABS_TOL_FRACTION",
    "SCALE_ABS_TOL_CONVENTIONAL",
    "RECOGNIZED_REPLICATE_SAMPLE_VALUES",
    "download_url_for_donor",
    "download_donor_zip",
    "unpack_zip_members",
    "scale_detection_tolerances",
    "detect_tpm_scale",
    "filter_biological_samples",
    "compute_sample_column_sums",
    "find_gene_tpm_row",
    "parse_donor_package",
    "pooled_mean_across_donors",
]
