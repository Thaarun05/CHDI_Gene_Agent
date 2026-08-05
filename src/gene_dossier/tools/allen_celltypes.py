"""Allen Cell Types (10x taxonomy) client for Section 2c.

Two datasets, both *dataset-level*: their source files are shared by every gene,
downloaded or registered once, and reused.

- Human M1 10x: official download of ``trimmed_means.csv`` + ``dend.json``.
- Mouse whole cortex + hippocampus 10x: locally registered trimmed means and
  dendrogram (the published bundle is not served from a stable per-file URL).

Structured analysis and browser figures are independent evidence paths: a
failed screenshot never invalidates a valid expression analysis.

Terminology: a value greater than zero is reported as *nonzero aggregate
expression*, never as "detected". "Detected" is reserved for a source-defined
detection threshold, which these trimmed-mean matrices do not provide.

Never raises: failures return :class:`~gene_dossier.models.ToolResult` or
structured error dicts.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import statistics
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence
from urllib.parse import quote, urlparse

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Allen Brain Atlas"

CALCULATION_VERSION = "allen_celltypes_v1"

# -- Human M1 official source configuration -------------------------------------------
HUMAN_M1_SOURCES: dict[str, str] = {
    "trimmed_means": (
        "https://idk-etl-prod-download-bucket.s3.amazonaws.com/"
        "aibs_human_m1_10x/trimmed_means.csv"
    ),
    "taxonomy": (
        "https://idk-etl-prod-download-bucket.s3.amazonaws.com/"
        "aibs_human_m1_10x/dend.json"
    ),
}
DOWNLOAD_HOST = "idk-etl-prod-download-bucket.s3.amazonaws.com"
ALLOWED_DOWNLOAD_HOSTS = frozenset({DOWNLOAD_HOST})

CACHE_KEY_HUMAN_TRIMMED_MEANS = "allen_human_m1_trimmed_means"
CACHE_KEY_HUMAN_TAXONOMY = "allen_human_m1_taxonomy"
CACHE_KEY_MOUSE_TRIMMED_MEANS = "mouse_ctx_hpf_trimmed_means"
CACHE_KEY_MOUSE_TAXONOMY = "mouse_ctx_hpf_taxonomy"

HUMAN_CACHE_KEYS = {
    "trimmed_means": CACHE_KEY_HUMAN_TRIMMED_MEANS,
    "taxonomy": CACHE_KEY_HUMAN_TAXONOMY,
}

# Canonical on-disk names for the locally registered mouse sources.
MOUSE_TRIMMED_MEANS_FILENAME = "mouse_ctx_hpf_trimmed_means.csv"
MOUSE_TAXONOMY_FILENAME = "mouse_ctx_hpf_dend.json"

REQUEST_HEADERS = {
    "User-Agent": "GeneDossier/0.1.0 (research; provenance-first gene dossier client)",
    "Accept": "text/csv,application/json,*/*;q=0.8",
}
SAFE_RESPONSE_HEADER_KEYS = ("content-type", "content-length", "server", "last-modified")

# -- Datasets / figures ----------------------------------------------------------------
DATASET_HUMAN_M1 = "human-m1-10x"
DATASET_MOUSE_CTX_HPF = "mouse-whole-cortex-and-hippocampus-10x"

DATASET_LABELS = {
    DATASET_HUMAN_M1: "Human M1 10x",
    DATASET_MOUSE_CTX_HPF: "Mouse whole cortex and hippocampus 10x",
}
# Terminology is dataset-specific: do not describe every source as snRNA-seq.
DATASET_ASSAY_TERMS = {
    DATASET_HUMAN_M1: "single-nucleus RNA sequencing",
    DATASET_MOUSE_CTX_HPF: "single-cell RNA sequencing",
}
DATASET_SAMPLING_SCOPE = {
    DATASET_HUMAN_M1: "sampled human primary motor cortex cell types",
    DATASET_MOUSE_CTX_HPF: "sampled mouse cortex and hippocampus populations",
}

EXPLORER_HOST = "celltypes.brain-map.org"
VIEWER_FRAME_HOST = "transcriptomics.brain-map.org"
ALLOWED_FIGURE_HOSTS = frozenset(
    {
        EXPLORER_HOST,
        VIEWER_FRAME_HOST,
        "knowledge.brain-map.org",
        "brain-map.org",
        "www.brain-map.org",
    }
)
_EXPLORER_BASE = f"https://{EXPLORER_HOST}/rnaseq"

# Internal dataset keys stay stable for accepted sources / evidence. Explorer
# path segments are the official Allen Transcriptomics Explorer route IDs.
EXPLORER_DATASET_SLUGS = {
    DATASET_HUMAN_M1: "human_m1_10x",
    DATASET_MOUSE_CTX_HPF: "mouse_ctx-hpf_10x",
}
DATABASE_HUB_URL = (
    "https://brain-map.org/our-research/cell-types-taxonomies/"
    "cell-types-database-rna-seq-data"
)
DATASET_SOURCE_PAGES = {
    DATASET_HUMAN_M1: (
        "https://brain-map.org/our-research/cell-types-taxonomies/"
        "cell-types-database-rna-seq-data/human-m1-10x"
    ),
    DATASET_MOUSE_CTX_HPF: (
        "https://brain-map.org/our-research/cell-types-taxonomies/"
        "cell-types-database-rna-seq-data/"
        "mouse-whole-cortex-and-hippocampus-10x"
    ),
}

VISUALIZATION_SCATTER = "scatter"
VISUALIZATION_HEATMAP = "heatmap"
_VISUALIZATION_PARAM = {
    VISUALIZATION_SCATTER: "Scatter Plot",
    VISUALIZATION_HEATMAP: "Heatmap",
}

PLOT_CONTAINER_SELECTORS = (
    "canvas",
    "[class*='plot']",
    "[class*='visualization']",
    "[class*='heatmap']",
    "[class*='scatter']",
    "svg",
)
MIN_PLOT_DIMENSION = 50
MIN_COMPOSITE_WIDTH = 400
MIN_COMPOSITE_HEIGHT = 250
MIN_HEATMAP_WIDTH = 600
MIN_HEATMAP_HEIGHT = 350
MAX_COMPOSITE_HEIGHT = 2000
MAX_HEATMAP_HEIGHT = 1800

CAPTURE_STATUS_SUCCESS = "success"
CAPTURE_STATUS_UNAVAILABLE = "source_unavailable"
CAPTURE_STATUS_NOT_ATTEMPTED = "not_attempted_optional"

CAPTURE_VIEWPORT = {"width": 1600, "height": 1200}
CAPTURE_USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36"
)
CAPTURE_RETRIEVAL_METHOD = "official_web_element_capture"

# Substrings that mean the explorer itself reported a problem, so a screenshot
# would document an error page rather than the requested gene.
SOURCE_ERROR_MARKERS = (
    "no results",
    "not found",
    "does not exist",
    "unable to load",
    "something went wrong",
    "service unavailable",
    "access denied",
)

TOP_N_DEFAULT = 10

GENE_COLUMN_CANDIDATES = ("", "feature", "gene", "gene_symbol", "symbol", "id")

# Known reconciliation for the mouse bundle: the taxonomy carries one leaf that
# the expression matrix does not. Recorded, never fabricated as a zero.
MOUSE_TAXONOMY_LEAF_COUNT = 388
MOUSE_EXPRESSION_CLUSTER_COUNT = 387


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


def is_allowlisted_download_url(url: str | None) -> bool:
    """True when ``url`` is https on the approved Allen download bucket."""
    if not url:
        return False
    parsed = urlparse(url.strip())
    if parsed.scheme != "https":
        return False
    return (parsed.hostname or "").lower() in ALLOWED_DOWNLOAD_HOSTS


def human_m1_source_url(source_key: str) -> str:
    """Return the official Human M1 URL for ``trimmed_means`` or ``taxonomy``."""
    try:
        return HUMAN_M1_SOURCES[source_key]
    except KeyError as exc:  # pragma: no cover - programming error
        raise KeyError(f"unknown Human M1 source key: {source_key!r}") from exc


def validate_human_m1_payload(source_key: str, content: bytes) -> dict[str, Any]:
    """Validate a downloaded Human M1 source payload by kind."""
    if not content:
        return {"ok": False, "error_type": "empty_download", "error_message": "no bytes"}
    if _looks_like_html(content):
        return {
            "ok": False,
            "error_type": "html_masquerading_as_data",
            "error_message": "response body is HTML, not source data",
        }
    if source_key == "taxonomy":
        try:
            parsed = json.loads(content.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            return {"ok": False, "error_type": "malformed_json", "error_message": str(exc)}
        if not isinstance(parsed, dict):
            return {
                "ok": False,
                "error_type": "malformed_json",
                "error_message": "taxonomy root is not an object",
            }
        return {"ok": True, "error_type": None, "error_message": None}

    lines = content[:4096].decode("utf-8", "replace").splitlines()
    header = lines[0] if lines else ""
    if "," not in header:
        return {
            "ok": False,
            "error_type": "malformed_csv",
            "error_message": "header row has no comma-separated fields",
        }
    return {"ok": True, "error_type": None, "error_message": None}


def download_human_m1_source(
    source_key: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Download one official Human M1 dataset-level source file (never raises).

    One network request per call, so the caller records one ApiRun per attempt.
    Retries are the caller's choice; only a validated response may be accepted.
    """
    try:
        url = human_m1_source_url(source_key)
    except KeyError as exc:
        return _tool_result(
            endpoint_name="allen_human_m1_source_download",
            gene_symbol=gene_symbol,
            request_url="",
            request_params={"source_key": source_key},
            success=False,
            error_type="invalid_request",
            error_message=str(exc),
        )

    if not is_allowlisted_download_url(url):
        return _tool_result(
            endpoint_name="allen_human_m1_source_download",
            gene_symbol=gene_symbol,
            request_url=url,
            request_params={"source_key": source_key},
            success=False,
            error_type="url_not_allowlisted",
            error_message=f"host is not the approved Allen download bucket: {url!r}",
        )

    cfg = settings or get_settings()
    request_params = {
        "source_key": source_key,
        "url": url,
        "request_headers": dict(REQUEST_HEADERS),
        "follow_redirects": True,
    }
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds, follow_redirects=True) as client:
            response = client.get(url, headers=REQUEST_HEADERS)
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="allen_human_m1_source_download",
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="network_error",
            error_message=str(exc),
        )

    content = response.content or b""
    meta = {
        "source_key": source_key,
        "cache_key": HUMAN_CACHE_KEYS.get(source_key),
        "official_url": url,
        "resolved_url": str(response.url) if response.url else url,
        "byte_size": len(content),
        "sha256": hashlib.sha256(content).hexdigest() if content else None,
        "content_type": response.headers.get("content-type"),
        "response_headers": _safe_response_headers(response.headers),
    }
    if not response.is_success:
        return _tool_result(
            endpoint_name="allen_human_m1_source_download",
            gene_symbol=gene_symbol,
            request_url=str(response.url) if response.url else url,
            request_params=request_params,
            success=False,
            status_code=response.status_code,
            data=meta,
            error_type="http_error",
            error_message=f"HTTP {response.status_code}",
        )

    check = validate_human_m1_payload(source_key, content)
    if not check["ok"]:
        return _tool_result(
            endpoint_name="allen_human_m1_source_download",
            gene_symbol=gene_symbol,
            request_url=str(response.url) if response.url else url,
            request_params=request_params,
            success=False,
            status_code=response.status_code,
            data={**meta, "raw_text_preview": content[:400].decode("utf-8", "replace")},
            error_type=str(check["error_type"]),
            error_message=str(check["error_message"]),
        )

    return _tool_result(
        endpoint_name="allen_human_m1_source_download",
        gene_symbol=gene_symbol,
        request_url=str(response.url) if response.url else url,
        request_params=request_params,
        success=True,
        status_code=response.status_code,
        data={**meta, "content": content},
    )


def register_local_source(path: str | Path, *, canonical_name: str) -> dict[str, Any]:
    """Register a locally supplied dataset source file.

    Used for the mouse bundle, which is not served from a stable per-file URL.
    Returns the bytes plus provenance recording the original filename so the
    rename to a canonical name stays auditable.
    """
    src = Path(path)
    if not src.is_file():
        return {
            "ok": False,
            "error_type": "local_source_missing",
            "error_message": f"no file at {src}",
        }
    content = src.read_bytes()
    if not content:
        return {"ok": False, "error_type": "empty_source", "error_message": f"{src} is empty"}
    if _looks_like_html(content):
        return {
            "ok": False,
            "error_type": "html_masquerading_as_data",
            "error_message": f"{src} contains HTML",
        }
    return {
        "ok": True,
        "error_type": None,
        "error_message": None,
        "content": content,
        "canonical_name": canonical_name,
        "original_filename": src.name,
        "original_path": str(src),
        "byte_size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "retrieval_method": "local_registration",
    }


# --------------------------------------------------------------------------------------
# Trimmed-means parsing
# --------------------------------------------------------------------------------------
def normalize_symbol(symbol: str) -> str:
    """Documented normalization: strip whitespace and surrounding quotes only."""
    return str(symbol or "").strip().strip('"').strip("'")


def _resolve_gene_column(header: Sequence[str]) -> int:
    if not header:
        return 0
    if normalize_symbol(header[0]).casefold() in GENE_COLUMN_CANDIDATES:
        return 0
    for idx, name in enumerate(header):
        if normalize_symbol(name).casefold() in GENE_COLUMN_CANDIDATES and name != "":
            return idx
    return 0


@dataclass
class GeneRowExtraction:
    """Result of pulling one exact gene row out of a trimmed-means matrix."""

    ok: bool
    error_type: str | None = None
    error_message: str | None = None
    source_symbol: str | None = None
    match_count: int = 0
    celltype_labels: list[str] = field(default_factory=list)
    values: list[float | None] = field(default_factory=list)
    gene_row_count: int = 0
    malformed_row_count: int = 0


def extract_gene_row(stream: Iterable[str], gene_symbol: str) -> GeneRowExtraction:
    """Extract the exact-match row for ``gene_symbol`` (case-insensitive, no substrings).

    Duplicate exact matches fail closed: an ambiguous target is never averaged
    or silently resolved to the first row.
    """
    reader = csv.reader(stream)
    try:
        header = next(reader)
    except StopIteration:
        return GeneRowExtraction(
            ok=False, error_type="empty_matrix", error_message="no header row"
        )
    if len(header) < 2:
        return GeneRowExtraction(
            ok=False,
            error_type="malformed_matrix",
            error_message=f"header has {len(header)} field(s)",
        )

    gene_idx = _resolve_gene_column(header)
    labels = [normalize_symbol(h) for i, h in enumerate(header) if i != gene_idx]
    target = normalize_symbol(gene_symbol).casefold()
    if not target:
        return GeneRowExtraction(
            ok=False, error_type="invalid_request", error_message="empty gene symbol"
        )

    out = GeneRowExtraction(ok=True, celltype_labels=labels)
    expected = len(header)
    for row in reader:
        if not row or all(not str(c).strip() for c in row):
            continue
        if len(row) != expected:
            out.malformed_row_count += 1
            continue
        out.gene_row_count += 1
        symbol = normalize_symbol(row[gene_idx])
        if symbol.casefold() != target:
            continue
        out.match_count += 1
        out.source_symbol = symbol
        values: list[float | None] = []
        for i, cell in enumerate(row):
            if i == gene_idx:
                continue
            text = (cell or "").strip()
            if text == "" or text.upper() in {"NA", "NAN", "NULL"}:
                values.append(None)
                continue
            try:
                value = float(text)
            except ValueError:
                values.append(None)
                continue
            values.append(value if math.isfinite(value) else None)
        out.values = values

    if out.match_count == 0:
        out.ok = False
        out.error_type = "gene_not_found"
        out.error_message = f"no exact row for {gene_symbol!r}"
    elif out.match_count > 1:
        out.ok = False
        out.error_type = "duplicate_gene_rows"
        out.error_message = f"{out.match_count} exact rows for {gene_symbol!r}"
    return out


def matrix_celltype_labels(stream: Iterable[str]) -> list[str]:
    """Return the expression matrix column labels (cell types / clusters)."""
    reader = csv.reader(stream)
    try:
        header = next(reader)
    except StopIteration:
        return []
    gene_idx = _resolve_gene_column(header)
    return [normalize_symbol(h) for i, h in enumerate(header) if i != gene_idx]


def text_lines(content: bytes) -> Iterable[str]:
    """Yield decoded lines from raw CSV bytes (never parsed as Excel)."""
    return io.StringIO(content.decode("utf-8", "replace"))


# --------------------------------------------------------------------------------------
# Taxonomy (dend.json)
# --------------------------------------------------------------------------------------
def _node_attrs(node: dict[str, Any]) -> dict[str, Any]:
    for key in ("leaf_attributes", "node_attributes"):
        values = node.get(key)
        if isinstance(values, list) and values and isinstance(values[0], dict):
            return values[0]
    return {}


def _node_label(attrs: dict[str, Any]) -> str:
    for key in ("cell_set_alias", "original_label", "label", "cell_set_designation"):
        value = attrs.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def _is_anonymous_internal_label(label: str) -> bool:
    """True for placeholder internal labels such as ``n1`` that carry no meaning."""
    return bool(label) and label[0] == "n" and label[1:].isdigit()


@dataclass
class Taxonomy:
    """Leaf labels plus their meaningful ancestor labels."""

    leaves: dict[str, dict[str, Any]] = field(default_factory=dict)
    leaf_order: list[str] = field(default_factory=list)
    internal_label_count: int = 0

    @property
    def leaf_count(self) -> int:
        return len(self.leaf_order)


def parse_dendrogram(payload: Any) -> Taxonomy:
    """Parse an Allen ``dend.json`` into leaves with meaningful ancestry.

    Ancestors are the named internal nodes (for example ``Glutamatergic neurons``
    or ``Pvalb``). Anonymous structural nodes such as ``n17`` are skipped, and no
    ancestry is invented when the taxonomy does not encode it.
    """
    taxonomy = Taxonomy()
    if not isinstance(payload, dict):
        return taxonomy

    def walk(node: dict[str, Any], ancestors: tuple[str, ...]) -> None:
        attrs = _node_attrs(node)
        label = _node_label(attrs)
        children = node.get("children")
        if not isinstance(children, list) or not children:
            if not label:
                return
            taxonomy.leaves[label] = {
                "label": label,
                "ancestors": list(ancestors),
                "cell_set_accession": attrs.get("cell_set_accession"),
                "cell_set_designation": attrs.get("cell_set_designation"),
                "cell_set_structure": attrs.get("cell_set_structure"),
            }
            taxonomy.leaf_order.append(label)
            return
        next_ancestors = ancestors
        if label and not _is_anonymous_internal_label(label):
            taxonomy.internal_label_count += 1
            next_ancestors = ancestors + (label,)
        for child in children:
            if isinstance(child, dict):
                walk(child, next_ancestors)

    walk(payload, ())
    return taxonomy


def reconcile_taxonomy(
    *, taxonomy_leaves: Sequence[str], matrix_labels: Sequence[str]
) -> dict[str, Any]:
    """Reconcile taxonomy leaves against expression matrix columns.

    The intersection is what gets analyzed. A taxonomy leaf missing from the
    matrix is recorded, never given a fabricated zero.
    """
    leaf_set = set(taxonomy_leaves)
    matrix_set = set(matrix_labels)
    analyzed = [label for label in matrix_labels if label in leaf_set]
    return {
        "taxonomy_leaf_count": len(taxonomy_leaves),
        "expression_cluster_count": len(matrix_labels),
        "analyzed_intersection_count": len(analyzed),
        "missing_expression_clusters": sorted(leaf_set - matrix_set),
        "matrix_columns_absent_from_taxonomy": sorted(matrix_set - leaf_set),
    }


# --------------------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------------------
def summarize_celltype_expression(
    *,
    gene_symbol: str,
    source_symbol: str | None,
    celltype_labels: Sequence[str],
    values: Sequence[float | None],
    taxonomy: Taxonomy | None = None,
    count_noun: str = "celltype",
    top_n: int = TOP_N_DEFAULT,
) -> dict[str, Any]:
    """Deterministic per-gene statistics over one trimmed-means row.

    Uses ``nonzero_*`` naming throughout: these matrices provide no source-defined
    detection threshold, so "detected" would overstate what a positive value means.
    """
    pairs = [
        (label, value)
        for label, value in zip(celltype_labels, values)
        if label and value is not None and math.isfinite(value)
    ]
    valid = [v for _, v in pairs]
    nonzero = [v for v in valid if v > 0]
    total_available = len(celltype_labels)
    valid_count = len(valid)
    nonzero_count = len(nonzero)
    nonzero_percentage = (nonzero_count / valid_count * 100.0) if valid_count else 0.0

    ranked = sorted(pairs, key=lambda kv: (-kv[1], kv[0]))
    maximum = max(valid) if valid else None
    median = statistics.median(valid) if valid else None
    mean = statistics.fmean(valid) if valid else None
    max_to_median = (maximum / median) if (maximum is not None and median) else None

    def decorate(entries: Sequence[tuple[str, float]]) -> list[dict[str, Any]]:
        out: list[dict[str, Any]] = []
        for label, value in entries:
            item: dict[str, Any] = {"label": label, "value": value}
            if taxonomy is not None:
                leaf = taxonomy.leaves.get(label)
                if leaf:
                    item["ancestors"] = leaf["ancestors"]
                    if leaf.get("cell_set_structure"):
                        item["structure"] = leaf["cell_set_structure"]
            out.append(item)
        return out

    count_key = f"nonzero_{count_noun}_count"
    valid_key = f"valid_{count_noun}_count"
    total_key = f"total_{count_noun}_count"

    summary: dict[str, Any] = {
        "calculation_version": CALCULATION_VERSION,
        "requested_symbol": gene_symbol,
        "source_symbol": source_symbol,
        total_key: total_available,
        valid_key: valid_count,
        count_key: nonzero_count,
        "nonzero_percentage": round(nonzero_percentage, 2),
        "minimum": min(valid) if valid else None,
        "maximum": maximum,
        "mean": mean,
        "median": median,
        "max_to_median_ratio": max_to_median,
        "top": decorate(ranked[:top_n]),
        "bottom": decorate(ranked[-top_n:]) if valid_count > top_n else [],
    }
    if taxonomy is not None:
        summary["top_taxonomy_ancestors"] = _top_ancestors(ranked[:top_n], taxonomy)
    return summary


def _top_ancestors(
    entries: Sequence[tuple[str, float]], taxonomy: Taxonomy
) -> list[dict[str, Any]]:
    """Count meaningful taxonomy ancestors across the top entries."""
    counts: dict[str, int] = {}
    for label, _ in entries:
        leaf = taxonomy.leaves.get(label)
        if not leaf:
            continue
        for ancestor in leaf["ancestors"]:
            counts[ancestor] = counts.get(ancestor, 0) + 1
    return [
        {"ancestor": name, "top_member_count": count}
        for name, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    ]


def build_celltype_summary_artifact(
    *,
    dataset: str,
    summary: dict[str, Any],
    reconciliation: dict[str, Any],
    match_count: int,
    source_checksums: dict[str, str | None],
    source_urls: dict[str, str | None],
) -> dict[str, Any]:
    """Assemble ``allen_{species}_celltype_summary.json`` content."""
    return {
        "calculation_version": CALCULATION_VERSION,
        "dataset": dataset,
        "dataset_label": DATASET_LABELS.get(dataset, dataset),
        "assay_terminology": DATASET_ASSAY_TERMS.get(dataset),
        "sampling_scope": DATASET_SAMPLING_SCOPE.get(dataset),
        "target_row_match_count": match_count,
        "summary": summary,
        "taxonomy_reconciliation": reconciliation,
        "source_checksums": source_checksums,
        "source_urls": source_urls,
    }


# --------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------
def explorer_dataset_slug(dataset: str) -> str:
    """Return the official Transcriptomics Explorer path slug for ``dataset``."""
    try:
        return EXPLORER_DATASET_SLUGS[dataset]
    except KeyError as exc:
        raise KeyError(f"unknown Allen dataset: {dataset!r}") from exc


def dataset_source_page(dataset: str) -> str:
    """Return the verified dataset-specific Allen source page URL."""
    try:
        return DATASET_SOURCE_PAGES[dataset]
    except KeyError as exc:
        raise KeyError(f"unknown Allen dataset: {dataset!r}") from exc


def database_hub_url() -> str:
    """Shared Cell Types Database RNA-Seq hub (runtime fallback only)."""
    return DATABASE_HUB_URL


def explorer_base_url(dataset: str) -> str:
    """Outer explorer URL for ``dataset`` without gene/visualization query params."""
    return f"{_EXPLORER_BASE}/{explorer_dataset_slug(dataset)}"


def figure_url(*, dataset: str, gene_symbol: str, visualization: str) -> str:
    """Build a gene-general Allen transcriptomics explorer URL.

    Uses the official explorer route slug (not the internal dataset key) plus
    gene-general query parameters so no per-gene URL is stored as configuration.
    """
    if visualization not in _VISUALIZATION_PARAM:
        raise KeyError(f"unknown visualization: {visualization!r}")
    symbol = normalize_symbol(gene_symbol)
    if not symbol:
        raise ValueError("gene symbol is required for a figure URL")
    params = (
        f"selectedVisualization={quote(_VISUALIZATION_PARAM[visualization])}"
        f"&colorByFeature={quote('Gene Expression')}"
        f"&colorByFeatureValue={quote(symbol)}"
    )
    return f"{explorer_base_url(dataset)}?{params}"


def is_allowlisted_figure_url(url: str | None) -> bool:
    """True when ``url`` stays on an approved Allen explorer/source host."""
    if not url:
        return False
    parsed = urlparse(url.strip())
    if parsed.scheme not in {"http", "https"}:
        return False
    return (parsed.hostname or "").lower() in ALLOWED_FIGURE_HOSTS


def is_landing_page(url: str | None, *, dataset: str) -> bool:
    """True when the browser is on the explorer landing page, not the dataset view.

    Compares the path against the *external* explorer slug, never the internal
    dataset identifier.
    """
    if not url:
        return True
    slug = explorer_dataset_slug(dataset)
    parsed = urlparse(url.strip())
    path = (parsed.path or "").rstrip("/")
    if not path.endswith(slug):
        return True
    return not parsed.query


def _viewer_frame_url(page: Any) -> str | None:
    for root in _search_roots(page)[1:]:
        try:
            url = str(getattr(root, "url", "") or "")
        except Exception:  # noqa: BLE001
            continue
        host = (urlparse(url).hostname or "").lower()
        if host == VIEWER_FRAME_HOST:
            return url
    return None


def validate_figure_capture(
    *,
    dataset: str,
    gene_symbol: str,
    final_url: str | None,
    gene_visible: bool,
    plot_width: int | None,
    plot_height: int | None,
    source_error_visible: bool = False,
    visualization: str | None = None,
    viewer_frame_url: str | None = None,
    plot_panel_count: int | None = None,
    includes_taxonomy_dendrogram: bool | None = None,
    includes_gene_expression_panel: bool | None = None,
    includes_celltype_panel: bool | None = None,
    require_composition: bool = False,
) -> dict[str, Any]:
    """Accept a figure capture only when it genuinely shows the requested gene."""
    reasons: list[str] = []
    if not is_allowlisted_figure_url(final_url):
        reasons.append("final_url_not_allowlisted")
    if is_landing_page(final_url, dataset=dataset):
        reasons.append("landing_page")
    if not gene_visible:
        reasons.append("requested_gene_not_visible")
    if source_error_visible:
        reasons.append("source_error_visible")
    if not plot_width or not plot_height:
        reasons.append("plot_container_missing")
    elif plot_width < MIN_PLOT_DIMENSION or plot_height < MIN_PLOT_DIMENSION:
        reasons.append("plot_container_too_small")

    if require_composition and visualization == VISUALIZATION_SCATTER:
        if plot_panel_count is not None and int(plot_panel_count) < 2:
            reasons.append("paired_scatter_incomplete")
        if includes_celltype_panel is False:
            reasons.append("celltype_panel_missing")
        if includes_gene_expression_panel is False:
            reasons.append("gene_expression_panel_missing")
    if require_composition and visualization == VISUALIZATION_HEATMAP:
        if includes_taxonomy_dendrogram is False:
            reasons.append("taxonomy_dendrogram_missing")

    return {
        "valid": not reasons,
        "rejection_reasons": reasons,
        "dataset": dataset,
        "explorer_slug": explorer_dataset_slug(dataset),
        "requested_gene": normalize_symbol(gene_symbol),
        "visualization": visualization,
        "final_url": final_url,
        "viewer_frame_url": viewer_frame_url,
        "plot_width": plot_width,
        "plot_height": plot_height,
        "plot_panel_count": plot_panel_count,
        "includes_taxonomy_dendrogram": includes_taxonomy_dendrogram,
        "includes_gene_expression_panel": includes_gene_expression_panel,
        "includes_celltype_panel": includes_celltype_panel,
        "require_composition": require_composition,
    }


# --------------------------------------------------------------------------------------
# Figure capture (Playwright)
# --------------------------------------------------------------------------------------
def _launch_playwright_chromium(pw: Any) -> tuple[Any, list[dict[str, Any]]]:
    """Launch the Chrome channel first, then bundled Chromium; keep both attempts."""
    attempts: list[dict[str, Any]] = []
    try:
        browser = pw.chromium.launch(headless=True, channel="chrome")
        attempts.append({"channel": "chrome", "success": True})
        return browser, attempts
    except Exception as chrome_exc:  # noqa: BLE001
        attempts.append(
            {
                "channel": "chrome",
                "success": False,
                "error_type": type(chrome_exc).__name__,
                "error_message": str(chrome_exc)[:400],
            }
        )
    browser = pw.chromium.launch(headless=True)
    attempts.append({"channel": "chromium", "success": True})
    return browser, attempts


def _selected_channel(attempts: Sequence[dict[str, Any]]) -> str | None:
    for attempt in reversed(list(attempts)):
        if attempt.get("success"):
            return str(attempt.get("channel"))
    return None


def _search_roots(page: Any) -> list[Any]:
    """The page plus any child frames.

    The Allen Transcriptomics Explorer renders its viewer inside an iframe, so
    the outer document alone carries neither the plot nor the gene label.
    """
    roots: list[Any] = [page]
    try:
        roots.extend(list(page.frames or [])[1:])
    except Exception:  # noqa: BLE001
        pass
    return roots


def _root_text(root: Any) -> str:
    try:
        return str(root.locator("body").inner_text(timeout=8_000) or "")
    except Exception:  # noqa: BLE001
        return ""


def _visible_text(page: Any) -> str:
    return "\n".join(_root_text(root) for root in _search_roots(page))


def _gene_visible(page: Any, symbol: str) -> bool:
    """True when the resolved symbol appears in the rendered page text."""
    target = normalize_symbol(symbol).casefold()
    if not target:
        return False
    return target in _visible_text(page).casefold()


def _gene_in_viewer_urls(page: Any, symbol: str) -> bool:
    """True when the explorer/iframe URL carries ``colorByFeatureValue=<symbol>``.

    Heatmap row labels are often canvas-rendered and absent from DOM text, so the
    verified deep-link state is accepted as gene-selection evidence when the
    heatmap composition checks also pass.
    """
    target = normalize_symbol(symbol)
    if not target:
        return False
    urls: list[str] = []
    try:
        urls.append(str(page.url or ""))
    except Exception:  # noqa: BLE001
        pass
    viewer = _viewer_frame_url(page)
    if viewer:
        urls.append(viewer)
    for url in urls:
        parsed = urlparse(url)
        query = (parsed.query or "").replace("+", " ")
        markers = (
            f"colorByFeatureValue={quote(target)}",
            f"colorByFeatureValue={target}",
            f"colorbyfeaturevalue={target.casefold()}",
        )
        folded = query.casefold()
        if any(marker.casefold() in folded for marker in markers):
            return True
    return False


def _heatmap_view_active(page: Any) -> bool:
    urls: list[str] = []
    try:
        urls.append(str(page.url or ""))
    except Exception:  # noqa: BLE001
        pass
    viewer = _viewer_frame_url(page)
    if viewer:
        urls.append(viewer)
    for url in urls:
        query = (urlparse(url).query or "").replace("+", " ").casefold()
        if "selectedvisualization=heatmap" in query:
            return True
    return False


def _source_error_visible(page: Any) -> bool:
    text = _visible_text(page).casefold()
    return any(marker in text for marker in SOURCE_ERROR_MARKERS)


def _find_plot_element(page: Any) -> tuple[Any | None, str | None, dict[str, Any] | None]:
    """Return the largest plausible plot element plus its selector and box.

    Kept for offline tests that inject a single canvas; production capture uses
    :func:`_select_capture_target` for visualization-aware scoring.
    """
    for root in _search_roots(page):
        node, selector, box = _find_plot_element_in_root(root)
        if node is not None:
            return node, selector, box
    return None, None, None


def _poll_for_plot_element(
    page: Any, *, timeout_ms: int, interval_ms: int = 1_000
) -> tuple[Any | None, str | None, dict[str, Any] | None]:
    """Poll page and frames for a plot container until ``timeout_ms`` elapses."""
    deadline = time.monotonic() + max(0.0, float(timeout_ms) / 1000.0)
    while True:
        node, selector, box = _find_plot_element(page)
        if node is not None:
            return node, selector, box
        if time.monotonic() >= deadline:
            return None, None, None
        try:
            page.wait_for_timeout(interval_ms)
        except Exception:  # noqa: BLE001
            time.sleep(interval_ms / 1000.0)


def _find_plot_element_in_root(
    root: Any,
) -> tuple[Any | None, str | None, dict[str, Any] | None]:
    best: tuple[Any, str, dict[str, Any]] | None = None
    for selector in PLOT_CONTAINER_SELECTORS:
        try:
            candidate = root.locator(selector)
            count = int(candidate.count() or 0)
        except Exception:  # noqa: BLE001
            continue
        for index in range(count):
            node = candidate.nth(index)
            try:
                if not node.is_visible():
                    continue
            except Exception:  # noqa: BLE001
                continue
            try:
                box = node.bounding_box()
            except Exception:  # noqa: BLE001
                continue
            if not box:
                continue
            width = int(box.get("width") or 0)
            height = int(box.get("height") or 0)
            if width < MIN_PLOT_DIMENSION or height < MIN_PLOT_DIMENSION:
                continue
            if best is None or width * height > int(best[2]["width"]) * int(best[2]["height"]):
                best = (node, selector, {"width": width, "height": height})
        if best is not None:
            break
    if best is None:
        return None, None, None
    return best[0], best[1], best[2]


def _safe_box(node: Any) -> dict[str, float] | None:
    try:
        box = node.bounding_box()
    except Exception:  # noqa: BLE001
        return None
    if not box:
        return None
    return {
        "x": float(box.get("x") or 0),
        "y": float(box.get("y") or 0),
        "width": float(box.get("width") or 0),
        "height": float(box.get("height") or 0),
    }


def _visible_plot_surfaces(root: Any) -> list[tuple[Any, dict[str, float]]]:
    surfaces: list[tuple[Any, dict[str, float]]] = []
    for selector in ("canvas", "svg"):
        try:
            locator = root.locator(selector)
            count = int(locator.count() or 0)
        except Exception:  # noqa: BLE001
            continue
        for index in range(count):
            node = locator.nth(index)
            try:
                if not node.is_visible():
                    continue
            except Exception:  # noqa: BLE001
                continue
            box = _safe_box(node)
            if not box:
                continue
            if box["width"] < MIN_PLOT_DIMENSION or box["height"] < MIN_PLOT_DIMENSION:
                continue
            surfaces.append((node, box))
    return surfaces


def _text_contains_any(text: str, markers: Sequence[str]) -> bool:
    folded = text.casefold()
    return any(marker.casefold() in folded for marker in markers)


def _reject_noise_candidate(text: str) -> bool:
    """Reject legend/logo/loading-only fragments that are not the composed view."""
    folded = text.casefold()
    noise_only = (
        "loading",
        "please wait",
        "allen institute logo",
        "powered by",
    )
    if any(marker in folded for marker in noise_only) and len(folded) < 80:
        return True
    return False


def _score_scatter_candidate(
    *,
    text: str,
    box: dict[str, float],
    plot_surfaces: Sequence[dict[str, float]],
    gene_symbol: str,
    from_viewer_frame: bool,
) -> dict[str, Any]:
    symbol = normalize_symbol(gene_symbol)
    score = 0
    reasons: list[str] = []
    panel_count = len(plot_surfaces)
    includes_celltype = _text_contains_any(text, ("cell type", "cell types", "taxonomy"))
    includes_expression = _text_contains_any(
        text, ("expression", "gene expression", "log2", "trimmed mean")
    )
    gene_visible = bool(symbol) and symbol.casefold() in text.casefold()
    # Allen often draws both panels into one WebGL canvas; treat a wide surface
    # with both panel headings as a paired composite.
    paired_single_canvas = (
        panel_count == 1
        and includes_celltype
        and includes_expression
        and box["width"] >= MIN_COMPOSITE_WIDTH * 1.5
        and box["height"] >= MIN_COMPOSITE_HEIGHT
    )
    effective_panels = max(panel_count, 2 if paired_single_canvas else 0)

    if effective_panels >= 2:
        score += 4
        reasons.append(
            "two_plot_surfaces" if panel_count >= 2 else "paired_single_canvas"
        )
    if includes_celltype:
        score += 2
        reasons.append("celltype_heading")
    if gene_visible:
        score += 2
        reasons.append("requested_gene_visible")
    if includes_expression:
        score += 2
        reasons.append("expression_heading")
    if box["width"] >= MIN_COMPOSITE_WIDTH and box["height"] >= MIN_COMPOSITE_HEIGHT:
        score += 1
        reasons.append("large_bbox")
    if panel_count >= 2:
        ys = sorted(s["y"] for s in plot_surfaces)
        if len(ys) >= 2 and abs(ys[0] - ys[1]) < max(80.0, 0.4 * box["height"]):
            score += 1
            reasons.append("shared_horizontal_parent")
    if from_viewer_frame:
        score += 1
        reasons.append("viewer_frame_origin")

    rejected = False
    reject_reasons: list[str] = []
    if effective_panels < 2:
        rejected = True
        reject_reasons.append("fewer_than_two_plot_surfaces")
    if not includes_celltype or not includes_expression:
        rejected = True
        reject_reasons.append("paired_panel_headings_missing")
    if _reject_noise_candidate(text):
        rejected = True
        reject_reasons.append("noise_only_candidate")
    if box["width"] < MIN_COMPOSITE_WIDTH or box["height"] < MIN_COMPOSITE_HEIGHT:
        rejected = True
        reject_reasons.append("composite_too_small")
    if box["height"] > MAX_COMPOSITE_HEIGHT:
        rejected = True
        reject_reasons.append("composite_too_tall")

    return {
        "score": score if not rejected else -1,
        "rejected": rejected,
        "reject_reasons": reject_reasons,
        "score_reasons": reasons,
        "plot_panel_count": effective_panels,
        "includes_celltype_panel": includes_celltype,
        "includes_gene_expression_panel": includes_expression,
        "requested_gene_visible": gene_visible,
        "includes_taxonomy_dendrogram": None,
    }


def _score_heatmap_candidate(
    *,
    text: str,
    box: dict[str, float],
    plot_surfaces: Sequence[dict[str, float]],
    gene_symbol: str,
    from_viewer_frame: bool,
    gene_confirmed_by_url: bool = False,
) -> dict[str, Any]:
    symbol = normalize_symbol(gene_symbol)
    score = 0
    reasons: list[str] = []
    gene_visible = bool(symbol) and (
        symbol.casefold() in text.casefold() or gene_confirmed_by_url
    )
    has_dendrogram = _text_contains_any(
        text, ("dendrogram", "taxonomy", "cluster", "subclass", "supertype")
    )
    has_heatmap_body = len(plot_surfaces) >= 1 or _text_contains_any(
        text, ("heatmap", "expression", "trimmed mean")
    )
    gene_like = [
        token
        for token in text.replace("\n", " ").split()
        if token.isalpha() and 2 <= len(token) <= 12 and token[0].isupper()
    ]
    comparison_labels = len({t.casefold() for t in gene_like}) >= 3

    if has_dendrogram:
        score += 3
        reasons.append("dendrogram_or_taxonomy")
    if gene_visible:
        score += 3
        reasons.append("target_gene_row" if not gene_confirmed_by_url else "target_gene_url")
    if comparison_labels:
        score += 2
        reasons.append("comparison_gene_labels")
    if has_heatmap_body:
        score += 2
        reasons.append("heatmap_body")
    if box["width"] >= MIN_HEATMAP_WIDTH and box["height"] >= MIN_HEATMAP_HEIGHT:
        score += 1
        reasons.append("substantial_size")
    if MIN_HEATMAP_HEIGHT <= box["height"] <= MAX_HEATMAP_HEIGHT:
        score += 2
        reasons.append("viewport_scale_height")
    if from_viewer_frame:
        score += 1
        reasons.append("viewer_frame_origin")

    rejected = False
    reject_reasons: list[str] = []
    if not gene_visible:
        rejected = True
        reject_reasons.append("target_gene_missing")
    if not has_dendrogram:
        rejected = True
        reject_reasons.append("dendrogram_context_missing")
    if not has_heatmap_body:
        rejected = True
        reject_reasons.append("heatmap_body_missing")
    if box["width"] < MIN_HEATMAP_WIDTH or box["height"] < MIN_HEATMAP_HEIGHT:
        rejected = True
        reject_reasons.append("heatmap_too_small")
    if box["height"] > MAX_HEATMAP_HEIGHT:
        rejected = True
        reject_reasons.append("heatmap_too_tall")
    if _reject_noise_candidate(text):
        rejected = True
        reject_reasons.append("noise_only_candidate")
    if len(plot_surfaces) == 1 and not comparison_labels and not has_dendrogram:
        rejected = True
        reject_reasons.append("bare_unlabeled_canvas")

    return {
        "score": score if not rejected else -1,
        "rejected": rejected,
        "reject_reasons": reject_reasons,
        "score_reasons": reasons,
        "plot_panel_count": len(plot_surfaces),
        "includes_celltype_panel": None,
        "includes_gene_expression_panel": gene_visible,
        "requested_gene_visible": gene_visible,
        "includes_taxonomy_dendrogram": has_dendrogram,
    }


def _candidate_nodes(root: Any) -> list[Any]:
    """Collect visible container candidates that may hold the composed view."""
    selectors = (
        "[class*='visualization']",
        "[class*='viewer']",
        "[class*='plot']",
        "[class*='heatmap']",
        "[class*='scatter']",
        "[class*='content']",
        "main",
        "section",
        "div",
    )
    nodes: list[Any] = []
    seen_ids: set[int] = set()
    for selector in selectors:
        try:
            locator = root.locator(selector)
            count = min(int(locator.count() or 0), 40)
        except Exception:  # noqa: BLE001
            continue
        for index in range(count):
            node = locator.nth(index)
            try:
                if not node.is_visible():
                    continue
            except Exception:  # noqa: BLE001
                continue
            box = _safe_box(node)
            if not box:
                continue
            if box["width"] < MIN_COMPOSITE_WIDTH or box["height"] < MIN_COMPOSITE_HEIGHT:
                continue
            node_id = id(node)
            if node_id in seen_ids:
                continue
            seen_ids.add(node_id)
            nodes.append(node)
            if len(nodes) >= 24:
                return nodes
    return nodes


def _select_capture_target(
    page: Any,
    *,
    visualization: str,
    gene_symbol: str,
) -> dict[str, Any]:
    """Score composite containers and return the highest-scoring valid target."""
    result: dict[str, Any] = {
        "node": None,
        "selector": None,
        "box": None,
        "capture_kind": (
            "paired_scatter" if visualization == VISUALIZATION_SCATTER else "complete_heatmap"
        ),
        "candidate_count": 0,
        "accepted_candidate_score": None,
        "plot_panel_count": None,
        "includes_celltype_panel": None,
        "includes_gene_expression_panel": None,
        "includes_taxonomy_dendrogram": None,
        "requested_gene_visible": None,
        "score_audit": [],
    }
    scored: list[tuple[int, Any, dict[str, Any], dict[str, float]]] = []
    gene_confirmed_by_url = _gene_in_viewer_urls(page, gene_symbol)
    for root_index, root in enumerate(_search_roots(page)):
        from_viewer = False
        try:
            root_url = str(getattr(root, "url", "") or "")
            from_viewer = (urlparse(root_url).hostname or "").lower() == VIEWER_FRAME_HOST
        except Exception:  # noqa: BLE001
            from_viewer = root_index > 0
        for node in _candidate_nodes(root):
            box = _safe_box(node)
            if not box:
                continue
            try:
                text = str(node.inner_text(timeout=2_000) or "")
            except Exception:  # noqa: BLE001
                text = ""
            surfaces: list[dict[str, float]] = []
            try:
                for _surface, sbox in _visible_plot_surfaces(root):
                    cx = sbox["x"] + sbox["width"] / 2
                    cy = sbox["y"] + sbox["height"] / 2
                    if (
                        box["x"] <= cx <= box["x"] + box["width"]
                        and box["y"] <= cy <= box["y"] + box["height"]
                    ):
                        surfaces.append(sbox)
            except Exception:  # noqa: BLE001
                surfaces = []

            if visualization == VISUALIZATION_SCATTER:
                meta = _score_scatter_candidate(
                    text=text,
                    box=box,
                    plot_surfaces=surfaces,
                    gene_symbol=gene_symbol,
                    from_viewer_frame=from_viewer,
                )
            else:
                meta = _score_heatmap_candidate(
                    text=text,
                    box=box,
                    plot_surfaces=surfaces,
                    gene_symbol=gene_symbol,
                    from_viewer_frame=from_viewer,
                    gene_confirmed_by_url=gene_confirmed_by_url,
                )
            result["candidate_count"] += 1
            result["score_audit"].append(
                {
                    "score": meta["score"],
                    "rejected": meta["rejected"],
                    "reject_reasons": meta["reject_reasons"],
                    "score_reasons": meta["score_reasons"],
                    "width": int(box["width"]),
                    "height": int(box["height"]),
                    "from_viewer_frame": from_viewer,
                }
            )
            if meta["rejected"]:
                continue
            scored.append((int(meta["score"]), node, meta, box))

    if not scored:
        # Fallback for offline fakes / sparse DOMs: largest plot surface plus
        # page-level panel signals. Production normally wins via scoring.
        node, selector, box = _find_plot_element(page)
        surfaces: list[dict[str, float]] = []
        for root in _search_roots(page):
            surfaces.extend([sbox for _, sbox in _visible_plot_surfaces(root)])
        page_text = _visible_text(page)
        symbol = normalize_symbol(gene_symbol)
        gene_visible = bool(symbol) and symbol.casefold() in page_text.casefold()
        result["node"] = node
        result["selector"] = selector
        result["box"] = box
        result["plot_panel_count"] = len(surfaces) if surfaces else (1 if node else 0)
        result["requested_gene_visible"] = gene_visible or _gene_visible(
            page, gene_symbol
        )
        if visualization == VISUALIZATION_SCATTER:
            result["includes_celltype_panel"] = _text_contains_any(
                page_text, ("cell type", "cell types", "taxonomy")
            )
            result["includes_gene_expression_panel"] = _text_contains_any(
                page_text, ("expression", "gene expression", "log2", "trimmed mean")
            )
            if (
                result["plot_panel_count"] == 1
                and result["includes_celltype_panel"]
                and result["includes_gene_expression_panel"]
                and box
                and int(box.get("width") or 0) >= MIN_COMPOSITE_WIDTH * 1.5
            ):
                result["plot_panel_count"] = 2
        else:
            result["includes_taxonomy_dendrogram"] = _text_contains_any(
                page_text,
                ("dendrogram", "taxonomy", "cluster", "subclass", "supertype"),
            )
            result["includes_gene_expression_panel"] = gene_visible
            if not result["requested_gene_visible"] and _gene_in_viewer_urls(
                page, gene_symbol
            ):
                result["requested_gene_visible"] = True
                result["includes_gene_expression_panel"] = True
        return result

    scored.sort(
        key=lambda item: (
            item[0],
            # Prefer viewport-scale heights over page-tall wrappers.
            0 if item[3]["height"] <= MAX_HEATMAP_HEIGHT else -1,
            item[3]["width"] * min(item[3]["height"], MAX_HEATMAP_HEIGHT),
        ),
        reverse=True,
    )
    score, node, meta, box = scored[0]
    result["node"] = node
    result["selector"] = "scored_composite_container"
    result["box"] = {"width": int(box["width"]), "height": int(box["height"])}
    result["accepted_candidate_score"] = score
    result["plot_panel_count"] = meta.get("plot_panel_count")
    result["includes_celltype_panel"] = meta.get("includes_celltype_panel")
    result["includes_gene_expression_panel"] = meta.get("includes_gene_expression_panel")
    result["includes_taxonomy_dendrogram"] = meta.get("includes_taxonomy_dendrogram")
    result["requested_gene_visible"] = meta.get("requested_gene_visible")
    return result


def _poll_for_capture_target(
    page: Any,
    *,
    visualization: str,
    gene_symbol: str,
    timeout_ms: int,
    interval_ms: int = 1_000,
) -> dict[str, Any]:
    deadline = time.monotonic() + max(0.0, float(timeout_ms) / 1000.0)
    last: dict[str, Any] | None = None
    while True:
        last = _select_capture_target(
            page, visualization=visualization, gene_symbol=gene_symbol
        )
        if last.get("node") is not None and last.get("accepted_candidate_score") is not None:
            return last
        if time.monotonic() >= deadline:
            return last or {
                "node": None,
                "selector": None,
                "box": None,
                "capture_kind": visualization,
                "candidate_count": 0,
            }
        try:
            page.wait_for_timeout(interval_ms)
        except Exception:  # noqa: BLE001
            time.sleep(interval_ms / 1000.0)


def capture_explorer_figure(
    *,
    dataset: str,
    gene_symbol: str,
    visualization: str,
    max_attempts: int = 2,
    navigation_timeout_ms: int = 45_000,
    plot_timeout_ms: int = 45_000,
    settle_ms: int = 8_000,
    page_factory: Callable[[], Any] | None = None,
) -> dict[str, Any]:
    """Capture one Allen transcriptomics explorer visualization (never raises).

    Returns a plain dict of browser-derived outputs — PNG bytes plus the
    validation evidence that justified accepting them. No ApiRun is fabricated
    here: the caller records one ApiRun per real browser navigation.

    A capture is accepted only when :func:`validate_figure_capture` passes, so an
    explorer landing page, a missing plot container, or a page that never shows
    the requested gene is reported as ``source_unavailable`` rather than saved.
    ``page_factory`` injects an already-open page (used by offline tests).
    """
    out: dict[str, Any] = {
        "ok": False,
        "status": CAPTURE_STATUS_UNAVAILABLE,
        "dataset": dataset,
        "dataset_label": DATASET_LABELS.get(dataset, dataset),
        "explorer_slug": (
            explorer_dataset_slug(dataset) if dataset in EXPLORER_DATASET_SLUGS else None
        ),
        "visualization": visualization,
        "requested_gene": normalize_symbol(gene_symbol),
        "source_page_url": None,
        "outer_page_url": None,
        "final_url": None,
        "viewer_frame_url": None,
        "dom_selector": None,
        "png": None,
        "plot_width": None,
        "plot_height": None,
        "capture_kind": None,
        "candidate_count": None,
        "accepted_candidate_score": None,
        "plot_panel_count": None,
        "includes_celltype_panel": None,
        "includes_gene_expression_panel": None,
        "includes_taxonomy_dendrogram": None,
        "validation": None,
        "browser_channel": None,
        "launch_attempts": [],
        "attempts": [],
        "error_type": None,
        "error_message": None,
        "retrieval_method": CAPTURE_RETRIEVAL_METHOD,
    }

    try:
        page_url = figure_url(
            dataset=dataset, gene_symbol=gene_symbol, visualization=visualization
        )
    except (KeyError, ValueError) as exc:
        out["error_type"] = "invalid_request"
        out["error_message"] = str(exc)
        return out
    out["source_page_url"] = page_url

    if page_factory is not None:
        attempt = _capture_once(
            page=page_factory(),
            out=out,
            page_url=page_url,
            dataset=dataset,
            gene_symbol=gene_symbol,
            visualization=visualization,
            navigation_timeout_ms=navigation_timeout_ms,
            plot_timeout_ms=plot_timeout_ms,
            settle_ms=settle_ms,
            attempt_number=1,
            pair_scatter_panels=False,
        )
        out["attempts"].append(attempt)
        return out

    try:
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # noqa: BLE001
        out["error_type"] = "playwright_unavailable"
        out["error_message"] = f"{type(exc).__name__}: {exc}"
        return out

    for attempt_number in range(1, max(1, int(max_attempts)) + 1):
        browser = None
        context = None
        try:
            with sync_playwright() as pw:
                try:
                    browser, launch_attempts = _launch_playwright_chromium(pw)
                except Exception as exc:  # noqa: BLE001
                    out["error_type"] = "browser_launch_failed"
                    out["error_message"] = f"{type(exc).__name__}: {exc}"
                    out["attempts"].append(
                        {"attempt": attempt_number, "error_type": "browser_launch_failed"}
                    )
                    continue
                out["launch_attempts"] = list(launch_attempts)
                out["browser_channel"] = _selected_channel(launch_attempts)
                context = browser.new_context(
                    viewport=dict(CAPTURE_VIEWPORT), user_agent=CAPTURE_USER_AGENT
                )
                page = context.new_page()
                attempt = _capture_once(
                    page=page,
                    out=out,
                    page_url=page_url,
                    dataset=dataset,
                    gene_symbol=gene_symbol,
                    visualization=visualization,
                    navigation_timeout_ms=navigation_timeout_ms,
                    plot_timeout_ms=plot_timeout_ms,
                    settle_ms=settle_ms,
                    attempt_number=attempt_number,
                    pair_scatter_panels=True,
                )
                out["attempts"].append(attempt)
        except Exception as exc:  # noqa: BLE001
            out["error_type"] = out["error_type"] or "capture_failed"
            out["error_message"] = out["error_message"] or f"{type(exc).__name__}: {exc}"
            out["attempts"].append(
                {"attempt": attempt_number, "error_type": "capture_failed", "detail": str(exc)[:400]}
            )
        finally:
            for closeable in (context, browser):
                try:
                    if closeable is not None:
                        closeable.close()
                except Exception:  # noqa: BLE001
                    pass
        if out["ok"]:
            break
    return out



def _transcriptomics_frame(page: Any) -> Any | None:
    for root in _search_roots(page)[1:]:
        try:
            url = str(getattr(root, "url", "") or "")
        except Exception:  # noqa: BLE001
            continue
        if (urlparse(url).hostname or "").lower() == VIEWER_FRAME_HOST:
            return root
    return None


def _ensure_heatmap_gene_row(page: Any, gene_symbol: str) -> dict[str, Any]:
    """Add the requested gene via the explorer 'Add Genes' control when missing."""
    symbol = normalize_symbol(gene_symbol)
    audit: dict[str, Any] = {"attempted": False, "added": False, "already_present": False}
    if not symbol:
        return audit
    if _gene_visible(page, symbol):
        # Prefer a gene that appears as a heatmap row label, not only in chrome.
        text = _visible_text(page)
        if symbol in text.split() or f"\n{symbol}\n" in f"\n{text}\n":
            audit["already_present"] = True
            return audit
    frame = _transcriptomics_frame(page) or page
    audit["attempted"] = True
    try:
        add_btn = frame.get_by_text("Add Genes", exact=False)
        if int(add_btn.count() or 0) == 0:
            add_btn = frame.get_by_text("Add Gene", exact=False)
        if int(add_btn.count() or 0) == 0:
            audit["error"] = "add_genes_control_missing"
            return audit
        add_btn.first.click(timeout=8_000)
        page.wait_for_timeout(800)
        control = frame.locator(".aibs-select__control, [class*='aibs-select__control']")
        if int(control.count() or 0):
            control.last.click(force=True, timeout=5_000)
        else:
            frame.get_by_text("Search for gene symbols", exact=False).first.click(
                force=True, timeout=5_000
            )
        page.wait_for_timeout(400)
        inp = frame.locator("input[id^='react-select'][id$='-input']").last
        inp.fill(symbol, force=True, timeout=5_000)
        page.wait_for_timeout(1_200)
        option = frame.locator(
            "[class*='aibs-select__option'], [id*='react-select'][id*='option']"
        )
        if int(option.count() or 0):
            option.first.click(timeout=5_000)
        else:
            inp.press("Enter")
        page.wait_for_timeout(2_500)
        audit["added"] = _gene_visible(page, symbol)
        if not audit["added"]:
            audit["error"] = "gene_row_not_visible_after_add"
    except Exception as exc:  # noqa: BLE001
        audit["error"] = f"{type(exc).__name__}: {exc}"[:400]
    return audit


def _screenshot_node(page: Any, node: Any) -> bytes:
    try:
        return node.screenshot(type="png", timeout=20_000)
    except TypeError:
        return node.screenshot(type="png")
    except Exception as first_exc:  # noqa: BLE001
        box_now = _safe_box(node)
        if not box_now or box_now["height"] <= MAX_COMPOSITE_HEIGHT:
            raise first_exc
        clip = {
            "x": max(0.0, box_now["x"]),
            "y": max(0.0, box_now["y"]),
            "width": min(box_now["width"], float(CAPTURE_VIEWPORT["width"])),
            "height": min(box_now["height"], float(MAX_HEATMAP_HEIGHT)),
        }
        return page.screenshot(type="png", clip=clip)


def _composite_side_by_side(left_png: bytes, right_png: bytes) -> bytes:
    """Join two PNG panels horizontally for the paired scatter figure."""
    from PIL import Image

    left = Image.open(io.BytesIO(left_png)).convert("RGB")
    right = Image.open(io.BytesIO(right_png)).convert("RGB")
    height = max(left.height, right.height)
    if left.height != height:
        left = left.resize(
            (max(1, int(left.width * height / left.height)), height), Image.Resampling.LANCZOS
        )
    if right.height != height:
        right = right.resize(
            (max(1, int(right.width * height / right.height)), height),
            Image.Resampling.LANCZOS,
        )
    combined = Image.new("RGB", (left.width + right.width, height), (255, 255, 255))
    combined.paste(left, (0, 0))
    combined.paste(right, (left.width, 0))
    buf = io.BytesIO()
    combined.save(buf, format="PNG")
    return buf.getvalue()


def _cell_type_scatter_url(dataset: str) -> str:
    """Scatter URL colored by cell type (companion panel for paired capture)."""
    return (
        f"{explorer_base_url(dataset)}"
        f"?selectedVisualization={quote(_VISUALIZATION_PARAM[VISUALIZATION_SCATTER])}"
        f"&colorByFeature={quote('Cell Type')}"
    )


def _capture_once(
    *,
    page: Any,
    out: dict[str, Any],
    page_url: str,
    dataset: str,
    gene_symbol: str,
    visualization: str,
    navigation_timeout_ms: int,
    plot_timeout_ms: int,
    settle_ms: int,
    attempt_number: int,
    pair_scatter_panels: bool = True,
) -> dict[str, Any]:
    """One navigate-validate-screenshot cycle. Mutates ``out`` on success."""
    attempt: dict[str, Any] = {"attempt": attempt_number, "source_page_url": page_url}
    try:
        page.goto(page_url, wait_until="domcontentloaded", timeout=navigation_timeout_ms)
    except Exception as exc:  # noqa: BLE001
        attempt["error_type"] = "navigation_failed"
        attempt["error_message"] = str(exc)[:400]
        out["error_type"] = "navigation_failed"
        out["error_message"] = attempt["error_message"]
        return attempt

    _poll_for_capture_target(
        page,
        visualization=visualization,
        gene_symbol=gene_symbol,
        timeout_ms=plot_timeout_ms,
    )
    try:
        page.wait_for_timeout(settle_ms)
    except Exception:  # noqa: BLE001
        pass

    add_audit = None
    if visualization == VISUALIZATION_HEATMAP:
        add_audit = _ensure_heatmap_gene_row(page, gene_symbol)
        attempt["heatmap_gene_add"] = add_audit
        out["heatmap_gene_add"] = add_audit
        try:
            page.wait_for_timeout(2_000)
        except Exception:  # noqa: BLE001
            pass

    selected = _select_capture_target(
        page, visualization=visualization, gene_symbol=gene_symbol
    )
    node = selected.get("node")
    box = selected.get("box")

    try:
        final_url = str(page.url or "")
    except Exception:  # noqa: BLE001
        final_url = ""
    viewer_url = _viewer_frame_url(page)
    attempt["final_url"] = final_url
    attempt["viewer_frame_url"] = viewer_url
    attempt["capture_kind"] = selected.get("capture_kind")
    attempt["candidate_count"] = selected.get("candidate_count")
    attempt["accepted_candidate_score"] = selected.get("accepted_candidate_score")
    out["final_url"] = final_url
    out["outer_page_url"] = final_url
    out["viewer_frame_url"] = viewer_url
    out["capture_kind"] = selected.get("capture_kind")
    out["candidate_count"] = selected.get("candidate_count")
    out["accepted_candidate_score"] = selected.get("accepted_candidate_score")
    out["plot_panel_count"] = selected.get("plot_panel_count")
    out["includes_celltype_panel"] = selected.get("includes_celltype_panel")
    out["includes_gene_expression_panel"] = selected.get("includes_gene_expression_panel")
    out["includes_taxonomy_dendrogram"] = selected.get("includes_taxonomy_dendrogram")

    gene_visible = bool(selected.get("requested_gene_visible")) or _gene_visible(
        page, gene_symbol
    )
    # Heatmap acceptance requires the gene as a visible row label after Add Genes.
    if visualization == VISUALIZATION_HEATMAP and not gene_visible:
        attempt["error_type"] = "figure_validation_failed"
        attempt["error_message"] = "requested_gene_not_visible"
        out["error_type"] = "figure_validation_failed"
        out["error_message"] = attempt["error_message"]
        out["validation"] = validate_figure_capture(
            dataset=dataset,
            gene_symbol=gene_symbol,
            final_url=final_url or None,
            gene_visible=False,
            plot_width=int(box["width"]) if box else None,
            plot_height=int(box["height"]) if box else None,
            visualization=visualization,
            viewer_frame_url=viewer_url,
            require_composition=True,
            includes_taxonomy_dendrogram=selected.get("includes_taxonomy_dendrogram"),
        )
        return attempt

    force_paired_flags = visualization == VISUALIZATION_SCATTER
    validation = validate_figure_capture(
        dataset=dataset,
        gene_symbol=gene_symbol,
        final_url=final_url or None,
        gene_visible=gene_visible if visualization != VISUALIZATION_SCATTER else True,
        plot_width=int(box["width"]) if box else None,
        plot_height=int(box["height"]) if box else None,
        source_error_visible=_source_error_visible(page),
        visualization=visualization,
        viewer_frame_url=viewer_url,
        plot_panel_count=(
            2 if force_paired_flags else selected.get("plot_panel_count")
        ),
        includes_taxonomy_dendrogram=selected.get("includes_taxonomy_dendrogram"),
        includes_gene_expression_panel=(
            True
            if force_paired_flags
            else selected.get("includes_gene_expression_panel")
        ),
        includes_celltype_panel=(
            True if force_paired_flags else selected.get("includes_celltype_panel")
        ),
        require_composition=True,
    )
    # For gene-expression scatter, still require the gene / allowlisted URL before
    # spending a second navigation on the cell-type companion panel.
    if visualization == VISUALIZATION_SCATTER:
        precheck = validate_figure_capture(
            dataset=dataset,
            gene_symbol=gene_symbol,
            final_url=final_url or None,
            gene_visible=gene_visible,
            plot_width=int(box["width"]) if box else None,
            plot_height=int(box["height"]) if box else None,
            source_error_visible=_source_error_visible(page),
            visualization=None,
            viewer_frame_url=viewer_url,
            require_composition=False,
        )
        if not precheck["valid"] or node is None:
            reasons = list(precheck["rejection_reasons"]) or ["plot_container_missing"]
            attempt["validation"] = precheck
            out["validation"] = precheck
            attempt["error_type"] = "figure_validation_failed"
            attempt["error_message"] = ", ".join(reasons)
            out["error_type"] = "figure_validation_failed"
            out["error_message"] = attempt["error_message"]
            return attempt

    attempt["validation"] = validation
    attempt["dom_selector"] = selected.get("selector")
    out["validation"] = validation
    out["dom_selector"] = selected.get("selector")
    out["plot_width"] = validation.get("plot_width")
    out["plot_height"] = validation.get("plot_height")

    png: bytes | None = None
    if visualization == VISUALIZATION_SCATTER and pair_scatter_panels:
        # Current Allen UI shows one Color-By mode at a time. Capture Cell Type
        # and Gene Expression panels separately, then composite side-by-side.
        try:
            expr_png = _screenshot_node(page, node)
            cell_url = _cell_type_scatter_url(dataset)
            page.goto(
                cell_url, wait_until="domcontentloaded", timeout=navigation_timeout_ms
            )
            page.wait_for_timeout(max(settle_ms, 4_000))
            cell_selected = _select_capture_target(
                page, visualization=visualization, gene_symbol=gene_symbol
            )
            cell_node = cell_selected.get("node")
            if cell_node is None:
                raise RuntimeError("paired_scatter_panel_missing")
            cell_png = _screenshot_node(page, cell_node)
            png = _composite_side_by_side(cell_png, expr_png)
            out["includes_celltype_panel"] = True
            out["includes_gene_expression_panel"] = True
            out["plot_panel_count"] = 2
            out["capture_kind"] = "paired_scatter"
            out["dom_selector"] = "paired_scatter_composite"
            attempt["dom_selector"] = "paired_scatter_composite"
            out["outer_page_url"] = page_url
            out["final_url"] = page_url
            attempt["final_url"] = page_url
            attempt["cell_type_scatter_url"] = cell_url
            from PIL import Image as _PILImage

            _w, _h = _PILImage.open(io.BytesIO(png)).size
            # Validate provenance/composition using the captured panel geometry,
            # not the possibly tiny offline-fake composite pixel size.
            panel_w = int((box or {}).get("width") or _w)
            panel_h = int((box or {}).get("height") or _h)
            validation = validate_figure_capture(
                dataset=dataset,
                gene_symbol=gene_symbol,
                final_url=page_url,
                gene_visible=True,
                plot_width=max(panel_w, int(_w)),
                plot_height=max(panel_h, int(_h)),
                visualization=visualization,
                viewer_frame_url=viewer_url,
                plot_panel_count=2,
                includes_celltype_panel=True,
                includes_gene_expression_panel=True,
                require_composition=True,
            )
            attempt["validation"] = validation
            out["validation"] = validation
            out["plot_width"] = int(_w)
            out["plot_height"] = int(_h)
        except Exception as exc:  # noqa: BLE001
            attempt["error_type"] = "paired_scatter_failed"
            attempt["error_message"] = str(exc)[:400]
            out["error_type"] = "paired_scatter_failed"
            out["error_message"] = attempt["error_message"]
            return attempt
    else:
        if not validation["valid"] or node is None:
            reasons = list(validation["rejection_reasons"]) or ["plot_container_missing"]
            attempt["error_type"] = "figure_validation_failed"
            attempt["error_message"] = ", ".join(reasons)
            out["error_type"] = "figure_validation_failed"
            out["error_message"] = attempt["error_message"]
            return attempt
        try:
            png = _screenshot_node(page, node)
        except Exception as exc:  # noqa: BLE001
            attempt["error_type"] = "screenshot_failed"
            attempt["error_message"] = str(exc)[:400]
            out["error_type"] = "screenshot_failed"
            out["error_message"] = attempt["error_message"]
            return attempt

    if not png:
        attempt["error_type"] = "empty_screenshot"
        attempt["error_message"] = "screenshot returned no bytes"
        out["error_type"] = "empty_screenshot"
        out["error_message"] = attempt["error_message"]
        return attempt

    if visualization == VISUALIZATION_SCATTER and not validation.get("valid"):
        reasons = list(validation.get("rejection_reasons") or ["paired_scatter_incomplete"])
        attempt["error_type"] = "figure_validation_failed"
        attempt["error_message"] = ", ".join(reasons)
        out["error_type"] = "figure_validation_failed"
        out["error_message"] = attempt["error_message"]
        return attempt

    attempt["byte_size"] = len(png)
    attempt["sha256"] = hashlib.sha256(png).hexdigest()
    out["png"] = png
    out["ok"] = True
    out["status"] = CAPTURE_STATUS_SUCCESS
    out["error_type"] = None
    out["error_message"] = None
    return attempt


__all__ = [
    "ALLOWED_DOWNLOAD_HOSTS",
    "ALLOWED_FIGURE_HOSTS",
    "CAPTURE_RETRIEVAL_METHOD",
    "CAPTURE_STATUS_NOT_ATTEMPTED",
    "CAPTURE_STATUS_SUCCESS",
    "CAPTURE_STATUS_UNAVAILABLE",
    "CACHE_KEY_HUMAN_TAXONOMY",
    "CACHE_KEY_HUMAN_TRIMMED_MEANS",
    "CACHE_KEY_MOUSE_TAXONOMY",
    "CACHE_KEY_MOUSE_TRIMMED_MEANS",
    "CALCULATION_VERSION",
    "DATABASE_HUB_URL",
    "DATASET_ASSAY_TERMS",
    "DATASET_HUMAN_M1",
    "DATASET_LABELS",
    "DATASET_MOUSE_CTX_HPF",
    "DATASET_SAMPLING_SCOPE",
    "DATASET_SOURCE_PAGES",
    "DOWNLOAD_HOST",
    "EXPLORER_DATASET_SLUGS",
    "GeneRowExtraction",
    "HUMAN_CACHE_KEYS",
    "HUMAN_M1_SOURCES",
    "MOUSE_TAXONOMY_FILENAME",
    "MOUSE_TRIMMED_MEANS_FILENAME",
    "PLOT_CONTAINER_SELECTORS",
    "SOURCE_NAME",
    "SOURCE_ERROR_MARKERS",
    "Taxonomy",
    "VISUALIZATION_HEATMAP",
    "VISUALIZATION_SCATTER",
    "build_celltype_summary_artifact",
    "capture_explorer_figure",
    "database_hub_url",
    "dataset_source_page",
    "download_human_m1_source",
    "explorer_base_url",
    "explorer_dataset_slug",
    "extract_gene_row",
    "figure_url",
    "human_m1_source_url",
    "is_allowlisted_download_url",
    "is_allowlisted_figure_url",
    "is_landing_page",
    "matrix_celltype_labels",
    "normalize_symbol",
    "parse_dendrogram",
    "reconcile_taxonomy",
    "register_local_source",
    "summarize_celltype_expression",
    "text_lines",
    "validate_figure_capture",
    "validate_human_m1_payload",
]
