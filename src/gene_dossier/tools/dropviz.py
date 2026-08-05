"""DropViz client (Playwright acquisition + local R extraction).

DropViz is a session-bound Shiny portal for adult mouse brain single-cell
expression (http://dropviz.org/). There is no stable REST CSV/ZIP URL —
Playwright against saved-state or live UI is the primary acquisition path.

Scope (Section 2c collection only):
- inspect saved-state / dynamic Genex UI
- classify views from rendered content / download links
- download Shiny ZIP exports with provenance manifests
- extract rank tables via Rscript (missing R → partial_success)

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
Temporary Shiny session hrefs belong only in ``data.audit``, never payload.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import re
import shutil
import subprocess
import zipfile
from pathlib import Path
from typing import Any, Sequence
from urllib.parse import urlparse, urljoin

from gene_dossier.models import ToolResult

logger = logging.getLogger(__name__)

SOURCE_NAME = "DropViz"
BASE_URL = "http://dropviz.org/"
ALLOWED_HOSTS = frozenset({"dropviz.org", "www.dropviz.org"})

# Discovery aids only — never treat as guaranteed present.
CANDIDATE_TSNE_LINK_IDS = (
    "tsne.global.cluster.label.dl",
    "tsne.global.subcluster.label.dl",
    "tsne.local.label.dl",
)
CANDIDATE_RANK_LINK_IDS = (
    "gene.expr.rank.cluster.dl",
    "gene.expr.rank.subcluster.dl",
)

# Selectors verified against broadinstitute/dropviz app.R:
#   selectizeInput("user.genes", multiple=TRUE), actionButton("go", "Update!"),
#   navbarPage(id="top-nav"), tabsetPanel(id="mainpanel"/"clusterpanel"/
#   "subclusterpanel"), and downloadLink ids rendered by plotDownload().
QUERY_TAB_SELECTOR = '#top-nav a[data-value="Query"]'
QUERY_TAB_FALLBACK_SELECTOR = 'a[id="select.go.3"], a[id="select.analysis.tab"]'
GENE_SELECTIZE_INPUT_SELECTOR = 'input[id="user.genes-selectized"]'
# selectizeInput("tissue", "Limit By Region"). The local/regional t-SNE output
# renders a "choose a region" placeholder until this filter is committed.
REGION_SELECTIZE_INPUT_SELECTOR = 'input[id="tissue-selectized"]'
REGION_SELECT_ID = "tissue"
GENE_SELECT_ID = "user.genes"
GENE_SELECTIZE_DROPDOWN_OPTION = ".selectize-dropdown-content .option"


def selectize_dropdown_option_selector(select_id: str) -> str:
    """Options of one selectize control; the page hosts several dropdowns."""
    return (
        f'.form-group:has(select[id="{select_id}"]) '
        ".selectize-dropdown-content .option"
    )


def selectize_item_selector(select_id: str) -> str:
    """Committed chips of one selectize control."""
    return f'.form-group:has(select[id="{select_id}"]) .selectize-input .item'
# Scoped to the user.genes control: the Query page hosts several other
# selectize widgets whose chips would otherwise look like a gene selection.
GENE_SELECTED_ITEM_SELECTOR = (
    '.form-group:has(select[id="user.genes"]) .selectize-input .item'
)
GENE_SELECTED_ITEM_FALLBACK_SELECTOR = ".selectize-input .item"
QUERY_UPDATE_BUTTON_SELECTOR = 'button[id="go"], #go'

# Ordered capture plan: each entry names the tab path, the Shiny image output
# that must render, and the downloadLink that only enables once it has.
DYNAMIC_VIEW_PLAN: tuple[dict[str, str], ...] = (
    {
        "key": "rank",
        "main_tab": "clusters",
        "sub_tab_container": "clusterpanel",
        "sub_tab": "rank",
        "image_selector": '[id="gene.expr.rank.cluster.output"] img',
        "download_id": "gene.expr.rank.cluster.dl",
        "image_name": "rank_plot.png",
        "download_name": "rank.zip",
        "extract_dir": "rank",
    },
    {
        "key": "tsne_global",
        "main_tab": "clusters",
        "sub_tab_container": "clusterpanel",
        "sub_tab": "tsne",
        "image_selector": '[id="tsne.global.cluster.label"] img',
        "download_id": "tsne.global.cluster.label.dl",
        "image_name": "global_tsne.png",
        "download_name": "tsne_global.zip",
        "extract_dir": "tsne_global",
    },
    {
        "key": "tsne_local",
        "main_tab": "subclusters",
        "sub_tab_container": "subclusterpanel",
        "sub_tab": "tsne",
        "image_selector": '[id="tsne.local.label"] img',
        "download_id": "tsne.local.label.dl",
        "image_name": "local_tsne.png",
        "download_name": "tsne_local.zip",
        "extract_dir": "tsne_local",
        "requires_region": True,
    },
)

# tableDownload() exports stream reliably from the live server, unlike the
# plotDownload() ZIPs, so they are captured as a separate structured channel.
DYNAMIC_TABLE_EXPORTS: tuple[dict[str, str], ...] = (
    {
        "key": "clusters_table",
        "main_tab": "clusters",
        "sub_tab_container": "clusterpanel",
        "sub_tab": "Table",
        "download_id": "dt.clusters.dl",
        "download_name": "clusters_table.csv",
    },
    {
        "key": "subclusters_table",
        "main_tab": "subclusters",
        "sub_tab_container": "subclusterpanel",
        "sub_tab": "Table",
        "download_id": "dt.subclusters.dl",
        "download_name": "subclusters_table.csv",
    },
)

VIEW_RANK = "rank_expression_view"
VIEW_GLOBAL_TSNE = "global_tsne_view"
VIEW_REGIONAL_TSNE = "regional_tsne_view"
VIEW_MIXED = "mixed_view"
VIEW_UNKNOWN = "unknown_view"

ZIP_MAGIC = b"PK\x03\x04"
RANK_SORT_POLICY = (
    "descending by target.sum.per.100k; raw source row order preserved in "
    "clusters_top_raw.csv; CI validated lower_bound <= estimate <= upper_bound"
)

_STATE_FAILURE_STATUSES = frozenset(
    {
        "state_expired",
        "state_not_found",
        "saved_state_unavailable",
    }
)

# Production denylist markers — tests assert gene-specific state IDs / screenshots
# never appear in this module.
_PRODUCTION_DENYLIST_MARKERS = (
    # Intentionally empty of gene-specific state IDs / saved-state URLs.
)


def _tool_result(
    *,
    endpoint_name: str,
    gene_symbol: str,
    request_url: str,
    request_params: dict[str, Any] | None = None,
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
        request_params=request_params or {},
        status_code=status_code,
        data=data,
        error_type=error_type,
        error_message=error_message,
    )


def _envelope(
    *,
    status: str,
    payload: dict[str, Any],
    audit: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Nest status/payload/audit under ToolResult.data."""
    return {
        "status": status,
        "payload": payload,
        "audit": audit or {},
    }


def sha256_file(path: Path | str) -> str:
    """Return SHA-256 hex digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(content: bytes) -> str:
    """Return SHA-256 hex digest of bytes."""
    return hashlib.sha256(content).hexdigest()


def css_escape_id(element_id: str) -> str:
    """Escape ``.`` (and similar) for CSS selectors on Shiny download IDs."""
    return re.sub(r"([.!#$%&'*+,/:;<=>?@\[\]^`{|}~])", r"\\\1", element_id)


def download_selector(link_id: str) -> str:
    """Return a CSS selector for a Shiny downloadLink id."""
    return f"#{css_escape_id(link_id)}"


def is_allowed_dropviz_url(url: str) -> bool:
    """True when URL host is allowlisted (http or https)."""
    text = (url or "").strip()
    if not text:
        return False
    try:
        parsed = urlparse(text)
    except Exception:  # noqa: BLE001
        return False
    if parsed.scheme not in {"http", "https"}:
        return False
    host = (parsed.hostname or "").lower()
    return host in ALLOWED_HOSTS


def normalize_dropviz_url(url: str) -> str | None:
    """Return a canonical DropViz URL, or None if not allowlisted.

    Preserves http/https scheme (DropViz often serves bookmarks over http).
    Normalizes host to ``dropviz.org`` (strips www).
    """
    text = (url or "").strip()
    if not text:
        return None
    try:
        parsed = urlparse(text)
    except Exception:  # noqa: BLE001
        return None
    if parsed.scheme not in {"http", "https"}:
        return None
    host = (parsed.hostname or "").lower()
    if host not in ALLOWED_HOSTS:
        return None
    path = parsed.path or "/"
    query = f"?{parsed.query}" if parsed.query else ""
    fragment = f"#{parsed.fragment}" if parsed.fragment else ""
    return f"{parsed.scheme}://dropviz.org{path}{query}{fragment}"


def redirect_is_allowed(from_url: str, to_url: str) -> bool:
    """Allow HTTP→HTTPS (or host www) redirects only within ALLOWED_HOSTS."""
    if not is_allowed_dropviz_url(from_url):
        return False
    if not is_allowed_dropviz_url(to_url):
        return False
    return True


def is_html_payload(content: bytes | str) -> bool:
    """True when body looks like HTML rather than a binary export."""
    if isinstance(content, bytes):
        head = content[:200].lstrip()
        try:
            text = head.decode("utf-8", errors="ignore")
        except Exception:  # noqa: BLE001
            return False
    else:
        text = str(content or "")
    head = text.lstrip("\ufeff").lstrip()[:200].lower()
    return head.startswith("<!doctype") or head.startswith("<html")


def is_zip_bytes(content: bytes) -> bool:
    """True when bytes start with ZIP local-file magic."""
    return isinstance(content, (bytes, bytearray)) and content[:4] == ZIP_MAGIC


def inventory_zip_basenames(zip_path: Path | str) -> dict[str, Any]:
    """Inventory ZIP member basenames; reject ambiguous duplicate basenames.

    Returns ``{"ok": True, "members": [...], "basenames": [...]}`` or
    ``{"ok": False, "error_type": "ambiguous_duplicate_basename", ...}``.
    """
    path = Path(zip_path)
    try:
        with zipfile.ZipFile(path, "r") as zf:
            members = zf.namelist()
    except zipfile.BadZipFile as exc:
        return {
            "ok": False,
            "error_type": "invalid_zip",
            "error_message": str(exc)[:400],
            "members": [],
            "basenames": [],
        }
    basenames: list[str] = []
    seen: dict[str, str] = {}
    duplicates: list[str] = []
    for name in members:
        base = Path(name).name
        if not base or name.endswith("/"):
            continue
        basenames.append(base)
        if base in seen and seen[base] != name:
            duplicates.append(base)
        seen[base] = name
    if duplicates:
        return {
            "ok": False,
            "error_type": "ambiguous_duplicate_basename",
            "error_message": (
                "ZIP contains duplicate basenames that would collide on unpack: "
                + ", ".join(sorted(set(duplicates)))
            ),
            "members": members,
            "basenames": basenames,
            "duplicate_basenames": sorted(set(duplicates)),
        }
    return {
        "ok": True,
        "members": members,
        "basenames": basenames,
        "basename_to_member": seen,
    }


def classify_view(
    *,
    download_link_ids: Sequence[str] | None = None,
    page_text: str | None = None,
) -> str:
    """Classify DropViz view from rendered download IDs / content — never state_id."""
    ids = {str(x) for x in (download_link_ids or []) if x}
    text = (page_text or "").lower()

    has_rank = any(
        i.startswith("gene.expr.rank.") or i.startswith("gene.expr.heatmap.")
        for i in ids
    ) or ("levels by cluster" in text and "transcripts per 100" in text)
    # Also detect rank via known IDs even if only discovery candidates.
    has_rank = has_rank or bool(ids & set(CANDIDATE_RANK_LINK_IDS))

    has_global_tsne = any(
        i.startswith("tsne.global.") for i in ids
    ) or ("global region space" in text and "t-sne" in text)
    has_global_tsne = has_global_tsne or bool(
        ids & {"tsne.global.cluster.label.dl", "tsne.global.subcluster.label.dl"}
    )

    has_regional_tsne = (
        "tsne.local.label.dl" in ids
        or "local cluster space" in text
        or "regional" in text and "t-sne" in text and "local" in text
    )

    has_tsne = has_global_tsne or has_regional_tsne or any(
        i.startswith("tsne.") for i in ids
    )

    if has_rank and has_tsne:
        return VIEW_MIXED
    if has_rank:
        return VIEW_RANK
    if has_regional_tsne and not has_global_tsne:
        return VIEW_REGIONAL_TSNE
    if has_global_tsne:
        return VIEW_GLOBAL_TSNE
    if has_tsne:
        return VIEW_GLOBAL_TSNE
    return VIEW_UNKNOWN


def find_rank_download_id(download_link_ids: Sequence[str]) -> str | None:
    """Prefer cluster rank download, then subcluster."""
    ids = list(download_link_ids or [])
    for candidate in CANDIDATE_RANK_LINK_IDS:
        if candidate in ids:
            return candidate
    for link_id in ids:
        if "rank" in link_id and link_id.endswith(".dl"):
            return link_id
    return None


def find_tsne_download_ids(download_link_ids: Sequence[str]) -> list[str]:
    """Return present candidate t-SNE download IDs in preferred order."""
    ids = set(download_link_ids or [])
    return [c for c in CANDIDATE_TSNE_LINK_IDS if c in ids]


def _finite_nonnegative(value: Any) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and number >= 0.0


def validate_rank_row(row: dict[str, Any]) -> dict[str, Any]:
    """Validate one clusters.top-style row; returns ok / reason."""
    label = (
        row.get("cx.disp")
        or row.get("cluster.disp")
        or row.get("subcluster.disp")
        or row.get("label")
        or ""
    )
    label = str(label).strip()
    if not label:
        return {"ok": False, "reason": "empty_label"}

    estimate = row.get("target.sum.per.100k")
    lower = row.get("target.sum.L.per.100k", row.get("lower_bound"))
    upper = row.get("target.sum.R.per.100k", row.get("upper_bound"))
    if not _finite_nonnegative(estimate):
        return {"ok": False, "reason": "non_finite_or_negative_estimate"}
    if lower is not None and not _finite_nonnegative(lower):
        return {"ok": False, "reason": "non_finite_or_negative_lower_bound"}
    if upper is not None and not _finite_nonnegative(upper):
        return {"ok": False, "reason": "non_finite_or_negative_upper_bound"}
    if lower is not None and upper is not None:
        try:
            lo_f, est_f, hi_f = float(lower), float(estimate), float(upper)
        except (TypeError, ValueError):
            return {"ok": False, "reason": "non_numeric_ci"}
        if not (lo_f <= est_f <= hi_f):
            return {"ok": False, "reason": "ci_bounds_violated"}
    return {"ok": True, "label": label}


def derive_rank_outputs_from_raw_csv(
    raw_csv_path: Path | str,
    output_dir: Path | str,
) -> dict[str, Any]:
    """Write ranked CSV + top_clusters.json from source-order raw CSV.

    Preserves ``clusters_top_raw.csv`` as-is (caller may already have written it).
    Sorting policy is recorded in the returned manifest.
    """
    raw_path = Path(raw_csv_path)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not raw_path.is_file():
        return {
            "ok": False,
            "status": "missing_clusters_top",
            "error_type": "missing_clusters_top",
            "error_message": f"Raw CSV not found: {raw_path}",
        }

    with raw_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if "target.sum.per.100k" not in fieldnames:
        return {
            "ok": False,
            "status": "missing_clusters_top",
            "error_type": "missing_expression_column",
            "error_message": "target.sum.per.100k column missing from raw CSV",
            "fieldnames": fieldnames,
        }

    validated: list[dict[str, Any]] = []
    invalid: list[dict[str, Any]] = []
    for idx, row in enumerate(rows):
        check = validate_rank_row(row)
        if check.get("ok"):
            validated.append(row)
        else:
            invalid.append({"row_index": idx, "reason": check.get("reason"), "row": row})

    if not validated:
        return {
            "ok": False,
            "status": "rank_validation_failed",
            "error_type": "rank_validation_failed",
            "error_message": "No valid rank rows after validation",
            "invalid_rows": invalid,
        }

    def _sort_key(row: dict[str, Any]) -> float:
        try:
            return -float(row.get("target.sum.per.100k") or 0.0)
        except (TypeError, ValueError):
            return 0.0

    ranked = sorted(validated, key=_sort_key)
    ranked_path = out_dir / "clusters_top_ranked.csv"
    with ranked_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ranked)

    top_clusters = []
    for row in ranked:
        label = (
            row.get("cx.disp")
            or row.get("cluster.disp")
            or row.get("subcluster.disp")
            or row.get("label")
            or ""
        )
        top_clusters.append(
            {
                "label": str(label),
                "target.sum.per.100k": _as_float(row.get("target.sum.per.100k")),
                "target.sum.L.per.100k": _as_float(
                    row.get("target.sum.L.per.100k", row.get("lower_bound"))
                ),
                "target.sum.R.per.100k": _as_float(
                    row.get("target.sum.R.per.100k", row.get("upper_bound"))
                ),
                "gene": row.get("gene"),
                "region": row.get("region.disp") or row.get("region.abbrev"),
            }
        )
    top_path = out_dir / "top_clusters.json"
    top_path.write_text(
        json.dumps({"clusters": top_clusters, "sort_policy": RANK_SORT_POLICY}, indent=2)
        + "\n",
        encoding="utf-8",
    )

    manifest = {
        "ok": True,
        "status": "success",
        "sort_policy": RANK_SORT_POLICY,
        "raw_csv": str(raw_path),
        "ranked_csv": str(ranked_path),
        "top_clusters_json": str(top_path),
        "raw_row_count": len(rows),
        "valid_row_count": len(validated),
        "invalid_row_count": len(invalid),
        "invalid_rows": invalid,
        "parent_sha256": sha256_file(raw_path),
        "ranked_sha256": sha256_file(ranked_path),
        "top_clusters_sha256": sha256_file(top_path),
    }
    manifest_path = out_dir / "rank_derived_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    manifest["manifest_path"] = str(manifest_path)
    return manifest


def _as_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def assess_regional_quantitative_fields(field_names: Sequence[str]) -> dict[str, Any]:
    """Fail-closed regional stats field assessment.

    Requires unambiguous cell/observation id, two coords, one region, one
    expression field (+ optional cluster). Ambiguous expression → unavailable.
    """
    names = [str(n) for n in field_names]
    lower_map = {n.lower(): n for n in names}

    def _pick(candidates: Sequence[str]) -> str | None:
        for c in candidates:
            if c in lower_map:
                return lower_map[c]
        return None

    cell = _pick(("cell", "cell_id", "barcode", "observation", "obs_id", "cell.id"))
    x = _pick(("v1", "x", "tsne_1", "tsne1", "umap_1", "umap1"))
    y = _pick(("v2", "y", "tsne_2", "tsne2", "umap_2", "umap2"))
    region = _pick(("region", "region.disp", "region.abbrev", "exp.label", "tissue"))
    cluster = _pick(
        ("cluster", "subcluster", "cx.disp", "cluster.disp", "subcluster.disp", "cx")
    )

    expression_candidates = [
        n
        for n in names
        if any(
            key in n.lower()
            for key in (
                "expr",
                "expression",
                "target.sum",
                "per.100k",
                "transcript",
                "umi",
                "count",
                "heat",
            )
        )
    ]
    # Exclude coordinate-like / id fields mistaken for expression.
    expression_candidates = [
        n
        for n in expression_candidates
        if n.lower()
        not in {
            "cell",
            "cell_id",
            "barcode",
            "v1",
            "v2",
            "x",
            "y",
            "region",
            "cluster",
            "subcluster",
        }
    ]

    if len(expression_candidates) != 1:
        return {
            "ok": False,
            "regional_quantitative_status": "unavailable",
            "regional_evidence_type": "figure_derived",
            "reason": "ambiguous_expression_field",
            "expression_candidates": expression_candidates,
            "fields": {
                "cell": cell,
                "x": x,
                "y": y,
                "region": region,
                "cluster": cluster,
            },
        }

    if not (cell and x and y and region):
        return {
            "ok": False,
            "regional_quantitative_status": "unavailable",
            "regional_evidence_type": "figure_derived",
            "reason": "missing_required_fields",
            "expression_candidates": expression_candidates,
            "fields": {
                "cell": cell,
                "x": x,
                "y": y,
                "region": region,
                "cluster": cluster,
                "expression": expression_candidates[0],
            },
        }

    return {
        "ok": True,
        "regional_quantitative_status": "available",
        "regional_evidence_type": "table_derived",
        "reason": None,
        "fields": {
            "cell": cell,
            "x": x,
            "y": y,
            "region": region,
            "cluster": cluster,
            "expression": expression_candidates[0],
        },
    }


def summarize_cluster_table_csv(
    csv_path: Path | str,
    gene_symbol: str,
    *,
    top_n: int = 10,
) -> dict[str, Any]:
    """Rank a DropViz cluster/subcluster table export by the gene's amount.

    The table export carries ``<Gene> Amount`` and ``<Gene> P-Val`` columns. It
    is *not* the rank plot's transcripts-per-100,000 with confidence bounds, so
    the summary labels its own units rather than borrowing the rank schema.
    """
    path = Path(csv_path)
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False,
            "status": "table_parse_failed",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:300],
            "source_csv": str(path),
        }

    if not rows:
        return {
            "ok": False,
            "status": "table_empty",
            "error_message": "Cluster table export contained no rows",
            "source_csv": str(path),
        }

    fields = list(rows[0].keys())
    amount_col = next(
        (f for f in fields if f.lower() == f"{gene_symbol.lower()} amount"), None
    )
    pval_col = next(
        (f for f in fields if f.lower() == f"{gene_symbol.lower()} p-val"), None
    )
    label_col = next(
        (f for f in fields if f.lower() in {"cluster", "subcluster"}), None
    )
    if not amount_col or not label_col:
        return {
            "ok": False,
            "status": "gene_columns_missing",
            "error_message": (
                f"No '{gene_symbol} Amount' / cluster column in {fields}"
            ),
            "source_csv": str(path),
            "columns": fields,
        }

    ranked: list[dict[str, Any]] = []
    for row in rows:
        amount = _as_float(row.get(amount_col))
        if amount is None:
            continue
        ranked.append(
            {
                "region": row.get("Region"),
                "class": row.get("Class"),
                "cluster": row.get(label_col),
                "amount": amount,
                "p_value": _as_float(row.get(pval_col)) if pval_col else None,
            }
        )
    ranked.sort(key=lambda item: item["amount"], reverse=True)

    return {
        "ok": True,
        "status": "success",
        "gene_symbol": gene_symbol,
        "source_csv": str(path),
        "sha256": sha256_file(path),
        "row_count": len(rows),
        "ranked_row_count": len(ranked),
        "amount_column": amount_col,
        "p_value_column": pval_col,
        "label_column": label_col,
        "units": "DropViz table export 'Amount' (not transcripts per 100,000)",
        "confidence_intervals_available": False,
        "top_clusters": ranked[:top_n],
    }


def _repo_scripts_dir() -> Path:
    """Return ``scripts/`` next to the installed package root when possible."""
    # src/gene_dossier/tools/dropviz.py → repo root is parents[3]
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "scripts"
        if (candidate / "extract_dropviz_rank.R").is_file():
            return candidate
    return here.parents[3] / "scripts"


def rscript_available() -> bool:
    """True when ``Rscript`` is on PATH."""
    return shutil.which("Rscript") is not None


def run_extract_dropviz_rank(
    *,
    zip_or_rdata_path: Path | str,
    output_dir: Path | str,
    rscript_path: str | None = None,
) -> dict[str, Any]:
    """Invoke ``extract_dropviz_rank.R``; never raises."""
    script = _repo_scripts_dir() / "extract_dropviz_rank.R"
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exe = rscript_path or shutil.which("Rscript")
    if not exe:
        return {
            "ok": False,
            "status": "rscript_unavailable",
            "error_type": "rscript_unavailable",
            "error_message": "Rscript not found on PATH",
            "output_dir": str(out_dir),
        }
    if not script.is_file():
        return {
            "ok": False,
            "status": "extraction_failed",
            "error_type": "extraction_script_missing",
            "error_message": f"Missing R script: {script}",
        }
    try:
        completed = subprocess.run(
            [exe, str(script), str(zip_or_rdata_path), str(out_dir)],
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
    except FileNotFoundError:
        return {
            "ok": False,
            "status": "rscript_unavailable",
            "error_type": "rscript_unavailable",
            "error_message": "Rscript executable disappeared during invoke",
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "ok": False,
            "status": "extraction_failed",
            "error_type": "extraction_timeout",
            "error_message": str(exc)[:400],
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False,
            "status": "extraction_failed",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:400],
        }

    stdout = (completed.stdout or "")[-2000:]
    stderr = (completed.stderr or "")[-2000:]
    if completed.returncode != 0:
        # Detect honest missing clusters.top
        combined = (stdout + "\n" + stderr).lower()
        if "clusters.top" in combined and (
            "not found" in combined or "missing" in combined
        ):
            status = "missing_clusters_top"
        else:
            status = "extraction_failed"
        return {
            "ok": False,
            "status": status,
            "error_type": status,
            "error_message": stderr or stdout or f"Rscript exit {completed.returncode}",
            "returncode": completed.returncode,
            "stdout": stdout,
            "stderr": stderr,
        }

    raw_csv = out_dir / "clusters_top_raw.csv"
    if not raw_csv.is_file():
        return {
            "ok": False,
            "status": "missing_clusters_top",
            "error_type": "missing_clusters_top",
            "error_message": "Rscript succeeded but clusters_top_raw.csv missing",
            "stdout": stdout,
            "stderr": stderr,
        }

    # Ensure ranked outputs exist (R may have written them; Python re-derives).
    derived = derive_rank_outputs_from_raw_csv(raw_csv, out_dir)
    return {
        "ok": bool(derived.get("ok")),
        "status": derived.get("status", "success" if derived.get("ok") else "extraction_failed"),
        "returncode": completed.returncode,
        "stdout": stdout,
        "stderr": stderr,
        "derived": derived,
        "raw_csv": str(raw_csv),
        "output_dir": str(out_dir),
        "api_run": None,  # never invent fake ApiRuns for derived R extraction
    }


def run_inspect_dropviz_rdata(
    *,
    zip_or_rdata_path: Path | str,
    output_dir: Path | str,
    rscript_path: str | None = None,
) -> dict[str, Any]:
    """Invoke ``inspect_dropviz_rdata.R``; never raises."""
    script = _repo_scripts_dir() / "inspect_dropviz_rdata.R"
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exe = rscript_path or shutil.which("Rscript")
    if not exe:
        return {
            "ok": False,
            "status": "rscript_unavailable",
            "error_type": "rscript_unavailable",
            "error_message": "Rscript not found on PATH",
            "output_dir": str(out_dir),
            "api_run": None,
        }
    if not script.is_file():
        return {
            "ok": False,
            "status": "extraction_failed",
            "error_type": "extraction_script_missing",
            "error_message": f"Missing R script: {script}",
            "api_run": None,
        }
    try:
        completed = subprocess.run(
            [exe, str(script), str(zip_or_rdata_path), str(out_dir)],
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
    except FileNotFoundError:
        return {
            "ok": False,
            "status": "rscript_unavailable",
            "error_type": "rscript_unavailable",
            "error_message": "Rscript executable disappeared during invoke",
            "api_run": None,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False,
            "status": "extraction_failed",
            "error_type": type(exc).__name__,
            "error_message": str(exc)[:400],
            "api_run": None,
        }

    inventory_path = out_dir / "rdata_inventory.json"
    inventory: dict[str, Any] | None = None
    if inventory_path.is_file():
        try:
            inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            inventory = None

    regional = None
    if inventory and isinstance(inventory.get("objects"), list):
        # Attempt regional field assessment on first data.frame-like object.
        for obj in inventory["objects"]:
            cols = obj.get("columns") or obj.get("names") or []
            if cols:
                regional = assess_regional_quantitative_fields(cols)
                break

    ok = completed.returncode == 0
    return {
        "ok": ok,
        "status": "success" if ok else "extraction_failed",
        "returncode": completed.returncode,
        "stdout": (completed.stdout or "")[-2000:],
        "stderr": (completed.stderr or "")[-2000:],
        "inventory_path": str(inventory_path) if inventory_path.is_file() else None,
        "inventory": inventory,
        "regional": regional,
        "api_run": None,
        "output_dir": str(out_dir),
    }


def _launch_playwright_chromium(pw: Any) -> tuple[Any, list[dict[str, Any]]]:
    """Launch Chrome channel first, then bundled Chromium; keep both attempts."""
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
    try:
        browser = pw.chromium.launch(headless=True)
        attempts.append({"channel": "chromium", "success": True})
        return browser, attempts
    except Exception as chromium_exc:  # noqa: BLE001
        attempts.append(
            {
                "channel": "chromium",
                "success": False,
                "error_type": type(chromium_exc).__name__,
                "error_message": str(chromium_exc)[:400],
            }
        )
        raise


# JavaScript: gene-aware Shiny-ready predicate (corrections #1 and #2).
#
# Homepage sample art and dormant download anchors are present in the DropViz
# DOM at all times, so neither may count as a rendered gene result.
_SHINY_READY_JS = r"""
/* dropviz-script: shiny-ready */
(args) => {
  const gene = ((args && args.gene) || '').trim();
  const homepageMarkers = (args && args.homepageMarkers) || [];
  const body = document.body;
  if (!body) return {ready: false, reason: 'no_body'};

  const pageText = (body.innerText || '');
  const restoreErrorPresent = /RestoreContext\s+initialization/i.test(pageText)
    || /Session\s+\S+\s+not\s+found/i.test(pageText);

  const isHomepageAsset = (src) => {
    const value = (src || '').toLowerCase();
    return homepageMarkers.some(marker => value.indexOf(marker) !== -1);
  };

  const anchors = Array.from(
    document.querySelectorAll('a.shiny-download-link, a[id$=".dl"]')
  );
  const anchorEnabled = (a) => {
    const href = a.getAttribute('href') || '';
    const ariaDisabled = (a.getAttribute('aria-disabled') || '').toLowerCase();
    const className = (a.className || '').toString().toLowerCase();
    if (!href) return false;
    if (ariaDisabled === 'true') return false;
    if (className.indexOf('disabled') !== -1) return false;
    if (a.hasAttribute('disabled')) return false;
    return true;
  };
  const enabledDownloads = anchors.filter(anchorEnabled).map(a => a.id || null);

  const imgs = Array.from(document.images || []);
  const homepageAssets = imgs.filter(img => isHomepageAsset(img.currentSrc || img.src));
  const genePlotImages = imgs.filter(img => {
    const src = (img.currentSrc || img.src || '');
    if (isHomepageAsset(src)) return false;
    if (!(img.complete && img.naturalWidth > 0 && img.naturalHeight > 0)) return false;
    return /^blob:/i.test(src) || /^data:/i.test(src)
      || src.indexOf('/session/') !== -1
      || !!img.closest('.shiny-image-output, .shiny-plot-output');
  }).map(img => ({
    id: img.id || null,
    src: (img.currentSrc || img.src || '').slice(0, 300),
    width: img.naturalWidth,
    height: img.naturalHeight,
  }));

  let geneVisible = false;
  if (gene) {
    const escaped = gene.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    const pattern = new RegExp('\\b' + escaped + '\\b', 'i');
    const selectizeItems = Array.from(
      document.querySelectorAll('.selectize-input .item, .selectize-input')
    ).map(el => el.textContent || '').join(' ');
    const plotAlts = imgs.map(
      img => (img.alt || '') + ' ' + (img.currentSrc || img.src || '')
    ).join(' ');
    const outputText = Array.from(
      document.querySelectorAll('.shiny-bound-output, .tab-pane.active, .dataTables_wrapper')
    ).map(el => el.innerText || '').join(' ');
    geneVisible = pattern.test(selectizeItems) || pattern.test(plotAlts)
      || pattern.test(outputText);
  }

  const geneSpecificEvidence = enabledDownloads.length > 0
    || genePlotImages.length > 0
    || geneVisible;
  const homepageOnly = !geneSpecificEvidence && homepageAssets.length > 0;

  const evidence = {
    shinyOk: (typeof window.Shiny !== 'undefined'),
    restore_error_present: restoreErrorPresent,
    enabled_download_count: enabledDownloads.length,
    enabled_download_ids: enabledDownloads,
    dormant_download_count: anchors.length - enabledDownloads.length,
    gene_specific_plot_count: genePlotImages.length,
    gene_specific_plots: genePlotImages.slice(0, 20),
    homepage_asset_count: homepageAssets.length,
    requested_gene_visible: geneVisible,
    homepage_only: homepageOnly,
  };

  const fail = (reason) => Object.assign({ready: false, reason: reason}, evidence);

  if (restoreErrorPresent) return fail('restore_error');
  if (body.classList.contains('shiny-disconnected')) return fail('shiny_disconnected');
  const disconnectOverlay = document.getElementById('shiny-disconnected-overlay');
  if (disconnectOverlay && disconnectOverlay.offsetParent !== null) {
    return fail('shiny_disconnected');
  }
  if (body.classList.contains('shiny-busy')) return fail('shiny_busy');
  if (homepageOnly) return fail('homepage_only');
  if (!geneSpecificEvidence) return fail('no_gene_specific_output');

  return Object.assign({ready: true, reason: 'ok'}, evidence);
}
"""


# Reads the ranking already rendered in the cluster table, so a dropped CSV
# export does not also cost us the region needed by the local t-SNE.
_CLUSTER_TABLE_TOP_ROW_JS = """
() => {
  /* dropviz-script: cluster-table-top-row */
  const visible = (el) => {
    const rect = el.getBoundingClientRect();
    return rect.width > 0 && rect.height > 0;
  };
  for (const table of Array.from(document.querySelectorAll('table'))) {
    if (!visible(table)) continue;
    const headers = Array.from(table.querySelectorAll('thead th')).map(
      th => (th.textContent || '').trim()
    );
    const regionIndex = headers.findIndex(h => h.toLowerCase() === 'region');
    if (regionIndex < 0) continue;
    const row = table.querySelector('tbody tr');
    if (!row) continue;
    const cells = Array.from(row.querySelectorAll('td')).map(
      td => (td.textContent || '').trim()
    );
    if (cells.length <= regionIndex) continue;
    const cell = (name) => {
      const i = headers.findIndex(h => h.toLowerCase() === name);
      return i >= 0 && i < cells.length ? cells[i] : null;
    };
    return {
      headers: headers,
      cells: cells,
      region: cells[regionIndex],
      cluster: cell('cluster'),
      cell_class: cell('class'),
    };
  }
  return null;
}
"""


def _ready_js_args(gene_symbol: str) -> dict[str, Any]:
    """Arguments passed to :data:`_SHINY_READY_JS`."""
    return {
        "gene": (gene_symbol or "").strip(),
        "homepageMarkers": list(HOMEPAGE_ASSET_MARKERS),
    }


def evaluate_shiny_ready(page: Any, gene_symbol: str) -> dict[str, Any]:
    """Evaluate the gene-aware readiness predicate. Never raises."""
    try:
        result = page.evaluate(_SHINY_READY_JS, _ready_js_args(gene_symbol))
    except Exception as exc:  # noqa: BLE001
        return {"ready": False, "reason": "evaluate_failed", "error": str(exc)[:300]}
    if not isinstance(result, dict):
        return {"ready": False, "reason": "invalid_ready_result"}
    return result


def acceptance_from_evidence(
    evidence: dict[str, Any],
    *,
    state_failure: str | None = None,
    gene_query_submitted: bool = False,
) -> dict[str, Any]:
    """Map readiness evidence onto the Section 2c acceptance criteria."""
    restore_error = bool(evidence.get("restore_error_present"))
    gene_plots = int(evidence.get("gene_specific_plot_count") or 0)
    return {
        "requested_gene_visible": bool(evidence.get("requested_gene_visible")),
        "homepage_only": bool(evidence.get("homepage_only")),
        "restore_error_present": restore_error,
        "gene_query_submitted": bool(gene_query_submitted),
        "gene_specific_plot_count": gene_plots,
        "enabled_download_count": int(evidence.get("enabled_download_count") or 0),
        "usable_for_section_2c": bool(
            evidence.get("ready") and not restore_error and state_failure is None
        ),
    }


_PAGE_INVENTORY_JS = """
/* dropviz-script: page-inventory */
() => {
  const body = document.body;
  const text = (body && body.innerText) ? body.innerText.slice(0, 20000) : '';
  const html = (body && body.innerHTML) ? body.innerHTML.slice(0, 5000) : '';
  const downloads = Array.from(
    document.querySelectorAll('a.shiny-download-link, a[id$=".dl"], a[download]')
  ).map(a => ({
    id: a.id || null,
    href: a.getAttribute('href') || null,
    text: (a.innerText || a.textContent || '').trim().slice(0, 120),
    className: a.className || '',
  }));
  const images = Array.from(document.images || []).map(img => ({
    id: img.id || null,
    src: (img.currentSrc || img.src || '').slice(0, 300),
    complete: !!img.complete,
    naturalWidth: img.naturalWidth || 0,
    naturalHeight: img.naturalHeight || 0,
    alt: (img.alt || '').slice(0, 120),
  }));
  const classes = body ? Array.from(body.classList || []) : [];
  return {
    title: document.title || '',
    text,
    html_preview: html,
    downloads,
    images,
    body_classes: classes,
    url: location.href,
  };
}
"""

# Assets that ship with the DropViz marketing homepage. Their presence never
# demonstrates that a gene-specific Shiny result rendered.
HOMEPAGE_ASSET_MARKERS = (
    "tsne-sample.jpg",
    "scatter-sample.jpg",
    "rank-sample.jpg",
    "stanleycenter-web.png",
    "hms.png",
    "mccarrolllab.org",
)

# Shiny bookmark restore failure, e.g.
# "Error in RestoreContext initialization: Session <state id> not found".
RESTORE_ERROR_PATTERNS = (
    re.compile(r"RestoreContext\s+initialization", re.I),
    re.compile(r"Session\s+\S+\s+not\s+found", re.I),
)


def detect_restore_error(page_text: str) -> bool:
    """True when the page shows a Shiny bookmark restore failure."""
    text = page_text or ""
    return any(pattern.search(text) for pattern in RESTORE_ERROR_PATTERNS)


_EXPIRED_PATTERNS = (
    re.compile(r"state\s+(not\s+found|expired|unavailable)", re.I),
    re.compile(r"bookmark.*(not\s+found|expired|invalid)", re.I),
    re.compile(r"saved\s+state.*(not\s+found|expired|unavailable)", re.I),
    re.compile(r"unable to restore", re.I),
    re.compile(r"_state_id_.*(invalid|unknown|expired)", re.I),
)


def detect_state_failure(page_text: str) -> str | None:
    """Return a state-failure status string, or None if content looks valid."""
    text = page_text or ""
    lower = text.lower()
    # A failed bookmark restore silently falls back to the homepage, so this
    # must be checked before any generic heuristics.
    if detect_restore_error(text):
        return "state_not_found"
    if "shiny-disconnected" in lower and "reconnect" in lower:
        return "saved_state_unavailable"
    if re.search(r"state\s+expired|saved\s+state\s+expired|bookmark.*expired", text, re.I):
        return "state_expired"
    if re.search(
        r"state\s+not\s+found|bookmark.*not\s+found|unknown state|error:\s*bookmark",
        text,
        re.I,
    ):
        return "state_not_found"
    if re.search(
        r"saved\s+state.*(unavailable|not\s+found)|unable to restore|_state_id_.*(invalid|unknown)",
        text,
        re.I,
    ):
        return "saved_state_unavailable"
    for pattern in _EXPIRED_PATTERNS:
        match = pattern.search(text)
        if not match:
            continue
        matched = match.group(0).lower()
        if "expired" in matched:
            return "state_expired"
        if "not found" in matched or "unavailable" in matched:
            return "state_not_found" if "not found" in matched else "saved_state_unavailable"
        return "saved_state_unavailable"
    return None


class DropVizClient:
    """Playwright-backed DropViz acquisition client (never raises)."""

    def __init__(
        self,
        *,
        navigation_timeout_ms: int = 90_000,
        shiny_ready_timeout_ms: int = 15_000,
    ) -> None:
        self.navigation_timeout_ms = navigation_timeout_ms
        self.shiny_ready_timeout_ms = shiny_ready_timeout_ms

    def _launch_browser(self, pw: Any) -> tuple[Any, list[dict[str, Any]], str]:
        browser, attempts = _launch_playwright_chromium(pw)
        channel = "chromium"
        for attempt in reversed(attempts):
            if attempt.get("success"):
                channel = str(attempt.get("channel") or "chromium")
                break
        return browser, attempts, channel

    def inspect_saved_state(
        self,
        *,
        gene_symbol: str,
        state_url: str,
        output_dir: Path | str,
        page: Any | None = None,
        browser_channel: str | None = None,
        launch_attempts: list[dict[str, Any]] | None = None,
    ) -> ToolResult:
        """Open a saved-state URL, wait for Shiny-ready, inventory, screenshot.

        When ``page`` is provided (tests / shared browser), launch is skipped.
        """
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        endpoint = "inspect_saved_state"
        request_params: dict[str, Any] = {
            "gene_symbol": gene_symbol,
            "state_url": state_url,
            "output_dir": str(out_dir),
        }

        if not is_allowed_dropviz_url(state_url):
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol,
                request_url=state_url,
                request_params=request_params,
                success=False,
                error_type="host_not_allowed",
                error_message=f"URL host not in ALLOWED_HOSTS: {state_url}",
                data=_envelope(
                    status="host_not_allowed",
                    payload={
                        "gene_symbol": gene_symbol,
                        "state_url": state_url,
                        "view_type": VIEW_UNKNOWN,
                        "artifacts": [],
                        "downloads": [],
                        "extractions": [],
                        "acquisition_status": "failed",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                    },
                    audit={"launch_attempts": launch_attempts or []},
                ),
            )

        canonical = normalize_dropviz_url(state_url) or state_url

        owns_browser = page is None
        browser = None
        context = None
        attempts = list(launch_attempts or [])
        channel = browser_channel or "unknown"
        network_events: list[dict[str, Any]] = []
        temporary_hrefs: list[str] = []

        try:
            if owns_browser:
                try:
                    from playwright.sync_api import sync_playwright
                except Exception as exc:  # noqa: BLE001
                    return _tool_result(
                        endpoint_name=endpoint,
                        gene_symbol=gene_symbol,
                        request_url=canonical,
                        request_params=request_params,
                        success=False,
                        error_type="playwright_unavailable",
                        error_message=str(exc)[:500],
                        data=_envelope(
                            status="playwright_unavailable",
                            payload={
                                "gene_symbol": gene_symbol,
                                "state_url": state_url,
                                "view_type": VIEW_UNKNOWN,
                                "artifacts": [],
                                "downloads": [],
                                "extractions": [],
                                "acquisition_status": "failed",
                                "rank_extraction_status": "not_attempted",
                                "regional_quantitative_status": "unavailable",
                                "regional_evidence_type": "figure_derived",
                            },
                            audit={"launch_attempts": []},
                        ),
                    )
                try:
                    pw_cm = sync_playwright()
                    pw = pw_cm.__enter__()
                except Exception as exc:  # noqa: BLE001
                    return _tool_result(
                        endpoint_name=endpoint,
                        gene_symbol=gene_symbol,
                        request_url=canonical,
                        request_params=request_params,
                        success=False,
                        error_type="playwright_unavailable",
                        error_message=str(exc)[:500],
                        data=_envelope(
                            status="playwright_unavailable",
                            payload={
                                "gene_symbol": gene_symbol,
                                "state_url": state_url,
                                "view_type": VIEW_UNKNOWN,
                                "artifacts": [],
                                "downloads": [],
                                "extractions": [],
                                "acquisition_status": "failed",
                                "rank_extraction_status": "not_attempted",
                                "regional_quantitative_status": "unavailable",
                                "regional_evidence_type": "figure_derived",
                            },
                            audit={"launch_attempts": []},
                        ),
                    )
                try:
                    browser, attempts, channel = self._launch_browser(pw)
                except Exception as launch_exc:  # noqa: BLE001
                    try:
                        pw_cm.__exit__(None, None, None)
                    except Exception:  # noqa: BLE001
                        pass
                    return _tool_result(
                        endpoint_name=endpoint,
                        gene_symbol=gene_symbol,
                        request_url=canonical,
                        request_params=request_params,
                        success=False,
                        error_type="browser_launch_failed",
                        error_message=str(launch_exc)[:500],
                        data=_envelope(
                            status="browser_launch_failed",
                            payload={
                                "gene_symbol": gene_symbol,
                                "state_url": state_url,
                                "view_type": VIEW_UNKNOWN,
                                "artifacts": [],
                                "downloads": [],
                                "extractions": [],
                                "acquisition_status": "failed",
                                "rank_extraction_status": "not_attempted",
                                "regional_quantitative_status": "unavailable",
                                "regional_evidence_type": "figure_derived",
                            },
                            audit={"launch_attempts": attempts, "browser_channel": None},
                        ),
                    )
                context = browser.new_context(accept_downloads=True)
                page = context.new_page()
            else:
                pw_cm = None  # type: ignore[assignment]

            assert page is not None

            def _on_response(response: Any) -> None:
                try:
                    resp_url = str(getattr(response, "url", "") or "")
                    status = getattr(response, "status", None)
                    network_events.append(
                        {
                            "url": resp_url[:500],
                            "status": status,
                            "resource_type": getattr(response, "request", None)
                            and getattr(response.request, "resource_type", None),
                        }
                    )
                    # Track session download hrefs for audit only.
                    if "/session/" in resp_url or "download" in resp_url.lower():
                        if resp_url not in temporary_hrefs:
                            temporary_hrefs.append(resp_url)
                except Exception:  # noqa: BLE001
                    return

            try:
                page.on("response", _on_response)
            except Exception:  # noqa: BLE001
                pass

            try:
                page.goto(
                    canonical,
                    wait_until="domcontentloaded",
                    timeout=self.navigation_timeout_ms,
                )
            except Exception as goto_exc:  # noqa: BLE001
                final_url = str(getattr(page, "url", "") or "")
                # about:blank / empty means navigation never committed — not a redirect.
                if not final_url or final_url.startswith("about:"):
                    return _tool_result(
                        endpoint_name=endpoint,
                        gene_symbol=gene_symbol,
                        request_url=canonical,
                        request_params=request_params,
                        success=False,
                        error_type="site_unavailable",
                        error_message=str(goto_exc)[:500],
                        data=_envelope(
                            status="site_unavailable",
                            payload={
                                "gene_symbol": gene_symbol,
                                "state_url": state_url,
                                "view_type": VIEW_UNKNOWN,
                                "artifacts": [],
                                "downloads": [],
                                "extractions": [],
                                "acquisition_status": "failed",
                                "rank_extraction_status": "not_attempted",
                                "regional_quantitative_status": "unavailable",
                                "regional_evidence_type": "figure_derived",
                            },
                            audit={
                                "browser_channel": channel,
                                "launch_attempts": attempts,
                                "final_url": final_url or None,
                                "goto_error": str(goto_exc)[:400],
                            },
                        ),
                    )
                if not is_allowed_dropviz_url(final_url):
                    return _tool_result(
                        endpoint_name=endpoint,
                        gene_symbol=gene_symbol,
                        request_url=canonical,
                        request_params=request_params,
                        success=False,
                        error_type="redirect_not_allowed",
                        error_message=f"Navigation left allowlisted hosts: {final_url}",
                        data=_envelope(
                            status="redirect_not_allowed",
                            payload={
                                "gene_symbol": gene_symbol,
                                "state_url": state_url,
                                "view_type": VIEW_UNKNOWN,
                                "artifacts": [],
                                "downloads": [],
                                "extractions": [],
                                "acquisition_status": "failed",
                                "rank_extraction_status": "not_attempted",
                                "regional_quantitative_status": "unavailable",
                                "regional_evidence_type": "figure_derived",
                            },
                            audit={
                                "browser_channel": channel,
                                "launch_attempts": attempts,
                                "final_url": final_url,
                                "goto_error": str(goto_exc)[:400],
                            },
                        ),
                    )
                # Continue — some Shiny apps throw soft navigation errors.

            final_url = str(getattr(page, "url", canonical) or canonical)
            if not is_allowed_dropviz_url(final_url):
                return _tool_result(
                    endpoint_name=endpoint,
                    gene_symbol=gene_symbol,
                    request_url=canonical,
                    request_params=request_params,
                    success=False,
                    error_type="redirect_not_allowed",
                    error_message=f"Final URL host not allowlisted: {final_url}",
                    data=_envelope(
                        status="redirect_not_allowed",
                        payload={
                            "gene_symbol": gene_symbol,
                            "state_url": state_url,
                            "view_type": VIEW_UNKNOWN,
                            "artifacts": [],
                            "downloads": [],
                            "extractions": [],
                            "acquisition_status": "failed",
                            "rank_extraction_status": "not_attempted",
                            "regional_quantitative_status": "unavailable",
                            "regional_evidence_type": "figure_derived",
                        },
                        audit={
                            "browser_channel": channel,
                            "launch_attempts": attempts,
                            "final_url": final_url,
                        },
                    ),
                )

            shiny_ready: dict[str, Any] = {"ready": False, "reason": "not_checked"}
            try:
                # Poll gene-aware readiness: no restore error, not
                # disconnected/busy, and gene-specific output present.
                deadline_ms = min(self.shiny_ready_timeout_ms, 60_000)
                steps = max(1, int(deadline_ms / 250))
                for _ in range(steps):
                    shiny_ready = evaluate_shiny_ready(page, gene_symbol)
                    if shiny_ready.get("ready") or shiny_ready.get("restore_error_present"):
                        break
                    try:
                        page.wait_for_timeout(250)
                    except Exception:  # noqa: BLE001
                        break
            except Exception as wait_exc:  # noqa: BLE001
                shiny_ready = {
                    "ready": False,
                    "reason": "wait_failed",
                    "error": str(wait_exc)[:300],
                }

            try:
                inventory = page.evaluate(_PAGE_INVENTORY_JS)
            except Exception as inv_exc:  # noqa: BLE001
                inventory = {
                    "title": "",
                    "text": "",
                    "downloads": [],
                    "images": [],
                    "body_classes": [],
                    "url": final_url,
                    "error": str(inv_exc)[:300],
                }

            page_text = str(inventory.get("text") or "")
            state_failure = detect_state_failure(page_text)
            if state_failure is None and shiny_ready.get("restore_error_present"):
                state_failure = "state_not_found"
            download_ids = [
                d.get("id")
                for d in (inventory.get("downloads") or [])
                if d.get("id")
            ]
            enabled_download_ids = [
                str(x) for x in (shiny_ready.get("enabled_download_ids") or []) if x
            ]
            # Capture temporary hrefs from inventory into audit only.
            for d in inventory.get("downloads") or []:
                href = d.get("href")
                if href and href.startswith("http") and href not in temporary_hrefs:
                    temporary_hrefs.append(href)
                elif href and href.startswith("/") and "session" in href:
                    abs_href = urljoin(final_url, href)
                    if abs_href not in temporary_hrefs:
                        temporary_hrefs.append(abs_href)

            # Dormant homepage anchors must never drive classification: only an
            # actually rendered gene result can name a view.
            if state_failure or not shiny_ready.get("ready"):
                view_type = VIEW_UNKNOWN
            else:
                view_type = classify_view(
                    download_link_ids=enabled_download_ids,
                    page_text=page_text,
                )

            artifacts: list[dict[str, Any]] = []
            # Screenshots only when images are naturally sized (or page screenshot).
            screenshot_status = "skipped"
            screenshot_path = out_dir / "page_screenshot.png"
            images = inventory.get("images") or []
            ready_images = [
                img
                for img in images
                if img.get("complete")
                and int(img.get("naturalWidth") or 0) > 0
                and int(img.get("naturalHeight") or 0) > 0
            ]
            try:
                if ready_images or shiny_ready.get("ready"):
                    page.screenshot(path=str(screenshot_path), full_page=True)
                    if screenshot_path.is_file() and screenshot_path.stat().st_size > 0:
                        artifacts.append(
                            {
                                "kind": "screenshot",
                                "path": str(screenshot_path),
                                "sha256": sha256_file(screenshot_path),
                                "byte_size": screenshot_path.stat().st_size,
                            }
                        )
                        screenshot_status = "success"
                    else:
                        screenshot_status = "empty"
                else:
                    screenshot_status = "images_not_ready"
            except Exception as shot_exc:  # noqa: BLE001
                screenshot_status = "failed"
                logger.info("DropViz screenshot failed: %s", shot_exc)

            acceptance = acceptance_from_evidence(
                shiny_ready, state_failure=state_failure
            )
            manifest = {
                "gene_symbol": gene_symbol,
                "state_url": state_url,
                "canonical_url": canonical,
                "final_url": final_url,
                "view_type": view_type,
                "shiny_ready": shiny_ready,
                "state_failure": state_failure,
                "acceptance": acceptance,
                "download_link_ids": download_ids,
                "enabled_download_link_ids": enabled_download_ids,
                "image_inventory": images,
                "screenshot_status": screenshot_status,
                "title": inventory.get("title"),
            }
            manifest_path = out_dir / "manifest.json"
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            artifacts.append(
                {
                    "kind": "manifest",
                    "path": str(manifest_path),
                    "sha256": sha256_file(manifest_path),
                }
            )

            network_path = out_dir / "network.json"
            network_path.write_text(
                json.dumps({"events": network_events[:500]}, indent=2) + "\n",
                encoding="utf-8",
            )
            artifacts.append(
                {
                    "kind": "network_manifest",
                    "path": str(network_path),
                    "sha256": sha256_file(network_path),
                }
            )

            if state_failure:
                status = state_failure
                acquisition_status = state_failure
                success = False
            elif not shiny_ready.get("ready"):
                reason = str(shiny_ready.get("reason") or "not_ready")
                if reason == "shiny_disconnected":
                    status = "shiny_disconnected"
                elif reason == "shiny_busy":
                    status = "shiny_busy"
                elif reason == "restore_error":
                    status = "state_not_found"
                elif reason in {
                    "homepage_only",
                    "no_gene_specific_output",
                    "no_plot_or_download",
                    "images_not_ready",
                }:
                    status = "blank_or_unready"
                else:
                    status = "shiny_not_ready"
                acquisition_status = status
                success = False
            else:
                status = "success"
                acquisition_status = "success"
                success = True

            payload = {
                "gene_symbol": gene_symbol,
                "state_url": state_url,
                "view_type": view_type,
                "artifacts": artifacts,
                "downloads": [],
                "extractions": [],
                "acquisition_status": acquisition_status,
                "rank_extraction_status": "not_attempted",
                "regional_quantitative_status": "unavailable",
                "regional_evidence_type": "figure_derived",
                "download_link_ids": download_ids,
                "enabled_download_link_ids": enabled_download_ids,
                "shiny_ready": shiny_ready,
                "acceptance": acceptance,
                "screenshot_status": screenshot_status,
            }
            audit = {
                "browser_channel": channel,
                "launch_attempts": attempts,
                "final_url": final_url,
                "temporary_download_hrefs": temporary_hrefs,
                "network_manifest_path": str(network_path),
                "manifest_path": str(manifest_path),
            }
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol,
                request_url=canonical,
                request_params=request_params,
                success=success,
                status_code=None,
                error_type=None if success else status,
                error_message=None if success else f"inspect_saved_state: {status}",
                data=_envelope(status=status, payload=payload, audit=audit),
            )
        except Exception as exc:  # noqa: BLE001
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol,
                request_url=canonical,
                request_params=request_params,
                success=False,
                error_type="inspect_failed",
                error_message=str(exc)[:500],
                data=_envelope(
                    status="inspect_failed",
                    payload={
                        "gene_symbol": gene_symbol,
                        "state_url": state_url,
                        "view_type": VIEW_UNKNOWN,
                        "artifacts": [],
                        "downloads": [],
                        "extractions": [],
                        "acquisition_status": "failed",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                    },
                    audit={
                        "browser_channel": channel,
                        "launch_attempts": attempts,
                        "temporary_download_hrefs": temporary_hrefs,
                    },
                ),
            )
        finally:
            if owns_browser:
                try:
                    if context is not None:
                        context.close()
                except Exception:  # noqa: BLE001
                    pass
                try:
                    if browser is not None:
                        browser.close()
                except Exception:  # noqa: BLE001
                    pass
                try:
                    if pw_cm is not None:
                        pw_cm.__exit__(None, None, None)
                except Exception:  # noqa: BLE001
                    pass

    def download_shiny_export(
        self,
        *,
        page: Any,
        link_id: str,
        output_path: Path | str,
        gene_symbol: str = "",
        enable_timeout_ms: int = 30_000,
        expected_kind: str = "zip",
    ) -> ToolResult:
        """Click a Shiny downloadLink and validate its bytes. Never raises.

        ``expected_kind`` selects ZIP validation (plotDownload exports) or CSV
        passthrough (tableDownload exports).
        """
        endpoint = "download_shiny_export"
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Attribute selector is more reliable than CSS-escaping dotted Shiny ids.
        selector = f'a[id="{link_id}"]'
        request_url = str(getattr(page, "url", BASE_URL) or BASE_URL)
        request_params = {
            "link_id": link_id,
            "selector": selector,
            "css_selector": download_selector(link_id),
        }
        temporary_hrefs: list[str] = []

        try:
            locator = page.locator(selector)
            # Wait for Shiny to enable the download link when possible.
            try:
                page.wait_for_function(
                    """(id) => {
                      const a = document.getElementById(id);
                      if (!a) return false;
                      const disabled = a.classList.contains('disabled')
                        || a.getAttribute('aria-disabled') === 'true';
                      return !disabled;
                    }""",
                    arg=link_id,
                    timeout=enable_timeout_ms,
                )
            except Exception:  # noqa: BLE001
                pass

            # Inspect enablement before click for honest error typing.
            meta = {}
            try:
                meta = page.evaluate(
                    """/* dropviz-script: download-meta */
                    (id) => {
                      const a = document.getElementById(id);
                      if (!a) return {present: false};
                      return {
                        present: true,
                        disabled: a.classList.contains('disabled')
                          || a.getAttribute('aria-disabled') === 'true',
                        href: a.getAttribute('href') || '',
                      };
                    }""",
                    link_id,
                )
            except Exception:  # noqa: BLE001
                meta = {}

            if meta.get("present") and meta.get("disabled"):
                return _tool_result(
                    endpoint_name=endpoint,
                    gene_symbol=gene_symbol or "UNKNOWN",
                    request_url=request_url,
                    request_params=request_params,
                    success=False,
                    error_type="download_link_disabled",
                    error_message=(
                        f"Shiny download link remains disabled: {link_id}"
                    ),
                    data=_envelope(
                        status="download_link_disabled",
                        payload={
                            "gene_symbol": gene_symbol,
                            "link_id": link_id,
                            "path": str(out_path),
                            "artifacts": [],
                            "downloads": [],
                            "extractions": [],
                            "acquisition_status": "failed",
                            "rank_extraction_status": "not_attempted",
                            "regional_quantitative_status": "unavailable",
                            "regional_evidence_type": "figure_derived",
                            "link_meta": meta,
                        },
                        audit={"temporary_download_hrefs": temporary_hrefs},
                    ),
                )

            # Shiny renders downloadLink with target="_blank"; left in place the
            # download event fires on a popup page instead of this one.
            try:
                page.evaluate(
                    """/* dropviz-script: strip-target */
                    (id) => {
                      const a = document.getElementById(id);
                      if (a) a.removeAttribute('target');
                      return true;
                    }""",
                    link_id,
                )
            except Exception:  # noqa: BLE001
                pass

            capture_via = "download_event"
            try:
                with page.expect_download(timeout=60_000) as download_info:
                    try:
                        locator.click(timeout=30_000)
                    except Exception:
                        # Hidden-but-enabled links: force as last resort.
                        locator.click(timeout=15_000, force=True)
                download = download_info.value
                href = str(getattr(download, "url", "") or "")
                if href:
                    temporary_hrefs.append(href)
                download.save_as(str(out_path))
            except Exception as download_exc:  # noqa: BLE001
                # Same-session fallback: fetch the session href through the
                # browser context so cookies and the Shiny session still apply.
                session_href = str(meta.get("href") or "")
                if not session_href:
                    raise
                absolute = urljoin(str(getattr(page, "url", BASE_URL)), session_href)
                if not is_allowed_dropviz_url(absolute):
                    raise
                context = page.context
                response = context.request.get(absolute, timeout=120_000)
                if not response.ok:
                    # Shiny's session download URL is single-use, so a 404 here
                    # usually just confirms the aborted first attempt.
                    raise RuntimeError(
                        f"context request {response.status} after primary failure: "
                        f"{str(download_exc)[:200]}"
                    )
                out_path.write_bytes(response.body())
                capture_via = "context_request"
                if absolute not in temporary_hrefs:
                    temporary_hrefs.append(absolute)
            request_params["capture_via"] = capture_via
        except Exception as exc:  # noqa: BLE001
            # "canceled" plus a socket hang up means the Shiny server aborted the
            # export stream, which is distinct from a client-side failure.
            text = str(exc).lower()
            if "canceled" in text or "socket hang up" in text:
                failure_status = "export_canceled_by_server"
            else:
                failure_status = "download_failed"
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol or "UNKNOWN",
                request_url=request_url,
                request_params=request_params,
                success=False,
                error_type=failure_status,
                error_message=str(exc)[:500],
                data=_envelope(
                    status=failure_status,
                    payload={
                        "gene_symbol": gene_symbol,
                        "link_id": link_id,
                        "path": str(out_path),
                        "artifacts": [],
                        "downloads": [],
                        "extractions": [],
                        "acquisition_status": "failed",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                    },
                    audit={"temporary_download_hrefs": temporary_hrefs},
                ),
            )

        if not out_path.is_file() or out_path.stat().st_size == 0:
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol or "UNKNOWN",
                request_url=request_url,
                request_params=request_params,
                success=False,
                error_type="empty_download",
                error_message="Downloaded file is missing or empty",
                data=_envelope(
                    status="empty_download",
                    payload={
                        "gene_symbol": gene_symbol,
                        "link_id": link_id,
                        "path": str(out_path),
                        "artifacts": [],
                        "downloads": [],
                        "extractions": [],
                        "acquisition_status": "failed",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                    },
                    audit={"temporary_download_hrefs": temporary_hrefs},
                ),
            )

        content = out_path.read_bytes()
        if is_html_payload(content):
            try:
                out_path.unlink(missing_ok=True)
            except Exception:  # noqa: BLE001
                pass
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol or "UNKNOWN",
                request_url=request_url,
                request_params=request_params,
                success=False,
                error_type="html_download_rejected",
                error_message="Download body looks like HTML, not a ZIP export",
                data=_envelope(
                    status="html_download_rejected",
                    payload={
                        "gene_symbol": gene_symbol,
                        "link_id": link_id,
                        "artifacts": [],
                        "downloads": [],
                        "extractions": [],
                        "acquisition_status": "failed",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                    },
                    audit={"temporary_download_hrefs": temporary_hrefs},
                ),
            )

        if expected_kind == "csv":
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol or "UNKNOWN",
                request_url=request_url,
                request_params=request_params,
                success=True,
                data=_envelope(
                    status="success",
                    payload={
                        "gene_symbol": gene_symbol,
                        "link_id": link_id,
                        "path": str(out_path),
                        "sha256": sha256_bytes(content),
                        "byte_size": len(content),
                        "content_kind": "csv",
                        "artifacts": [
                            {
                                "kind": "download",
                                "path": str(out_path),
                                "sha256": sha256_bytes(content),
                                "byte_size": len(content),
                            }
                        ],
                        "downloads": [{"link_id": link_id, "path": str(out_path)}],
                        "extractions": [],
                        "acquisition_status": "success",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                        "stable_endpoint": False,
                    },
                    audit={
                        "temporary_download_hrefs": temporary_hrefs,
                        "session_href_is_stable_endpoint": False,
                    },
                ),
            )

        if not is_zip_bytes(content):
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol or "UNKNOWN",
                request_url=request_url,
                request_params=request_params,
                success=False,
                error_type="invalid_zip_magic",
                error_message="Download does not start with ZIP magic bytes",
                data=_envelope(
                    status="invalid_zip_magic",
                    payload={
                        "gene_symbol": gene_symbol,
                        "link_id": link_id,
                        "path": str(out_path),
                        "byte_size": len(content),
                        "artifacts": [],
                        "downloads": [],
                        "extractions": [],
                        "acquisition_status": "failed",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                    },
                    audit={"temporary_download_hrefs": temporary_hrefs},
                ),
            )

        inventory = inventory_zip_basenames(out_path)
        if not inventory.get("ok"):
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=gene_symbol or "UNKNOWN",
                request_url=request_url,
                request_params=request_params,
                success=False,
                error_type=str(inventory.get("error_type") or "invalid_zip"),
                error_message=str(inventory.get("error_message") or "ZIP inventory failed"),
                data=_envelope(
                    status=str(inventory.get("error_type") or "invalid_zip"),
                    payload={
                        "gene_symbol": gene_symbol,
                        "link_id": link_id,
                        "path": str(out_path),
                        "zip_inventory": inventory,
                        "artifacts": [],
                        "downloads": [],
                        "extractions": [],
                        "acquisition_status": "failed",
                        "rank_extraction_status": "not_attempted",
                        "regional_quantitative_status": "unavailable",
                        "regional_evidence_type": "figure_derived",
                    },
                    audit={"temporary_download_hrefs": temporary_hrefs},
                ),
            )

        digest = sha256_bytes(content)
        download_record = {
            "link_id": link_id,
            "path": str(out_path),
            "sha256": digest,
            "byte_size": len(content),
            "zip_inventory": inventory,
            # Explicit: session href is not a stable endpoint.
            "stable_endpoint": False,
        }
        return _tool_result(
            endpoint_name=endpoint,
            gene_symbol=gene_symbol or "UNKNOWN",
            request_url=request_url,
            request_params=request_params,
            success=True,
            status_code=None,
            data=_envelope(
                status="success",
                payload={
                    "gene_symbol": gene_symbol,
                    "link_id": link_id,
                    "artifacts": [],
                    "downloads": [download_record],
                    "extractions": [],
                    "acquisition_status": "success",
                    "rank_extraction_status": "not_attempted",
                    "regional_quantitative_status": "unavailable",
                    "regional_evidence_type": "figure_derived",
                    "path": str(out_path),
                    "sha256": digest,
                    "stable_endpoint": False,
                },
                audit={
                    "temporary_download_hrefs": temporary_hrefs,
                    "session_href_is_stable_endpoint": False,
                },
            ),
        )

    # ---- dynamic Query workflow -------------------------------------------------

    def _open_query_tab(self, page: Any, diagnostics: dict[str, Any]) -> bool:
        """Navigate the navbar to the Query panel."""
        for selector in (QUERY_TAB_SELECTOR, QUERY_TAB_FALLBACK_SELECTOR):
            try:
                locator = page.locator(selector).first
                if locator.count() == 0:
                    continue
                locator.click(timeout=15_000)
                self._settle(page, 1_500)
                if page.locator(GENE_SELECTIZE_INPUT_SELECTOR).count() > 0:
                    diagnostics["query_tab_selector"] = selector
                    return True
            except Exception as exc:  # noqa: BLE001
                diagnostics.setdefault("query_tab_errors", []).append(
                    {"selector": selector, "error": str(exc)[:200]}
                )
        # The Query panel may already be active without an explicit click.
        try:
            if page.locator(GENE_SELECTIZE_INPUT_SELECTOR).count() > 0:
                diagnostics["query_tab_selector"] = "already_active"
                return True
        except Exception:  # noqa: BLE001
            pass
        return False

    def _settle(self, page: Any, timeout_ms: int) -> None:
        try:
            page.wait_for_timeout(timeout_ms)
        except Exception:  # noqa: BLE001
            pass

    def _commit_selectize_value(
        self,
        page: Any,
        *,
        input_selector: str,
        select_id: str,
        value: str,
        diagnostics: dict[str, Any],
        key: str,
        allow_create: bool,
    ) -> bool:
        """Type into a selectize control and commit a matching option."""
        try:
            field = page.locator(input_selector).first
            field.click(timeout=15_000)
            field.type(value, delay=60)
        except Exception as exc:  # noqa: BLE001
            diagnostics[f"{key}_input_error"] = str(exc)[:300]
            return False

        self._settle(page, 2_000)

        # Prefer a real autocomplete match over selectize's create=TRUE path.
        try:
            options = page.locator(selectize_dropdown_option_selector(select_id))
            if options.count() == 0:
                options = page.locator(GENE_SELECTIZE_DROPDOWN_OPTION)
            count = min(options.count(), 25)
            diagnostics[f"{key}_option_count"] = count
            texts = [
                (options.nth(i).inner_text() or "").strip() for i in range(count)
            ]
            index = next(
                (i for i, t in enumerate(texts) if t.lower() == value.lower()),
                None,
            )
            if index is None:
                index = next(
                    (i for i, t in enumerate(texts) if value.lower() in t.lower()),
                    None,
                )
            if index is not None:
                options.nth(index).click(timeout=10_000)
                diagnostics[f"{key}_selected"] = texts[index]
            elif allow_create:
                page.keyboard.press("Enter")
                diagnostics[f"{key}_selected"] = "enter_fallback"
            else:
                diagnostics[f"{key}_selected"] = None
                diagnostics[f"{key}_options"] = texts[:25]
                return False
        except Exception as exc:  # noqa: BLE001
            diagnostics[f"{key}_autocomplete_error"] = str(exc)[:300]
            if not allow_create:
                return False
            try:
                page.keyboard.press("Enter")
            except Exception:  # noqa: BLE001
                return False

        self._settle(page, 1_000)
        for selector in (
            selectize_item_selector(select_id),
            GENE_SELECTED_ITEM_FALLBACK_SELECTOR,
        ):
            try:
                items = page.locator(selector)
                count = items.count()
                if count == 0:
                    continue
                texts = [
                    (items.nth(i).inner_text() or "").strip()
                    for i in range(min(count, 10))
                ]
                diagnostics[f"{key}_committed_items"] = texts
                if any(t.lower() == value.lower() for t in texts):
                    return True
            except Exception as exc:  # noqa: BLE001
                diagnostics.setdefault(f"{key}_commit_errors", []).append(
                    {"selector": selector, "error": str(exc)[:200]}
                )
        return False

    def _select_gene(
        self,
        page: Any,
        gene_symbol: str,
        diagnostics: dict[str, Any],
    ) -> bool:
        return self._commit_selectize_value(
            page,
            input_selector=GENE_SELECTIZE_INPUT_SELECTOR,
            select_id=GENE_SELECT_ID,
            value=gene_symbol,
            diagnostics=diagnostics,
            key="gene",
            allow_create=True,
        )

    def _select_region(
        self,
        page: Any,
        region: str,
        diagnostics: dict[str, Any],
    ) -> bool:
        return self._commit_selectize_value(
            page,
            input_selector=REGION_SELECTIZE_INPUT_SELECTOR,
            select_id=REGION_SELECT_ID,
            value=region,
            diagnostics=diagnostics,
            key="region",
            allow_create=False,
        )

    def read_top_cluster_row(self, page: Any) -> dict[str, Any] | None:
        """Read the top row of the rendered cluster table.

        The table export is the richer source, but the live server drops export
        connections intermittently. The same ranking is already on screen, so
        reading it keeps region-dependent views working when the CSV is lost.
        """
        try:
            row = page.evaluate(_CLUSTER_TABLE_TOP_ROW_JS)
        except Exception:  # noqa: BLE001
            return None
        if not isinstance(row, dict) or not row.get("region"):
            return None
        return row

    def _submit_query(self, page: Any, diagnostics: dict[str, Any]) -> bool:
        try:
            page.locator(QUERY_UPDATE_BUTTON_SELECTOR).first.click(timeout=15_000)
            diagnostics["query_submitted"] = True
            self._settle(page, 3_000)
            return True
        except Exception as exc:  # noqa: BLE001
            diagnostics["query_submit_error"] = str(exc)[:300]
            diagnostics["query_submitted"] = False
            return False

    def _activate_view(self, page: Any, view: dict[str, str]) -> dict[str, Any]:
        """Switch the main and inner tabsets to the requested view."""
        result: dict[str, Any] = {"main_tab": view["main_tab"], "sub_tab": view["sub_tab"]}
        for selector in (
            f'#mainpanel a[data-value="{view["main_tab"]}"]',
            f'#{view["sub_tab_container"]} a[data-value="{view["sub_tab"]}"]',
        ):
            try:
                locator = page.locator(selector).first
                if locator.count() == 0:
                    result.setdefault("missing_tabs", []).append(selector)
                    continue
                locator.click(timeout=15_000)
                self._settle(page, 1_500)
            except Exception as exc:  # noqa: BLE001
                result.setdefault("tab_errors", []).append(
                    {"selector": selector, "error": str(exc)[:200]}
                )
        return result

    def _wait_for_plot_image(
        self,
        page: Any,
        selector: str,
        timeout_ms: int,
    ) -> dict[str, Any]:
        """Wait until a Shiny-rendered image has real pixel dimensions."""
        deadline = max(1, int(timeout_ms / 500))
        info: dict[str, Any] = {"rendered": False, "reason": "not_found"}
        for _ in range(deadline):
            try:
                info = page.evaluate(
                    """/* dropviz-script: plot-image */
                    (sel) => {
                      const img = document.querySelector(sel);
                      if (!img) return {rendered: false, reason: 'not_found'};
                      const src = img.currentSrc || img.src || '';
                      const ready = !!(img.complete && img.naturalWidth > 0
                        && img.naturalHeight > 0);
                      return {
                        rendered: ready,
                        reason: ready ? 'ok' : 'not_painted',
                        src: src.slice(0, 300),
                        width: img.naturalWidth,
                        height: img.naturalHeight,
                      };
                    }""",
                    selector,
                )
            except Exception as exc:  # noqa: BLE001
                info = {"rendered": False, "reason": "evaluate_failed", "error": str(exc)[:200]}
            if isinstance(info, dict) and info.get("rendered"):
                return info
            self._settle(page, 500)
        return info if isinstance(info, dict) else {"rendered": False, "reason": "unknown"}

    def collect_dynamic_query(
        self,
        *,
        page: Any,
        mouse_gene_symbol: str,
        output_dir: Path | str,
        network_events: list[dict[str, Any]] | None = None,
        browser_channel: str | None = None,
        launch_attempts: list[dict[str, Any]] | None = None,
        skip_extraction: bool = False,
        plot_timeout_ms: int = 90_000,
    ) -> ToolResult:
        """Run the real DropViz Query workflow for one mouse gene. Never raises.

        Opens the Query panel, commits the gene through the selectize control,
        submits, then captures the rank, global t-SNE and local t-SNE views and
        downloads each export once Shiny enables its link.
        """
        out_dir = Path(output_dir)
        images_dir = out_dir / "images"
        downloads_dir = out_dir / "downloads"
        extracted_dir = out_dir / "extracted"
        for directory in (out_dir, images_dir, downloads_dir, extracted_dir):
            directory.mkdir(parents=True, exist_ok=True)

        endpoint = "collect_dynamic_query"
        diagnostics: dict[str, Any] = {"mouse_gene_symbol": mouse_gene_symbol}
        artifacts: list[dict[str, Any]] = []
        downloads: list[dict[str, Any]] = []
        extractions: list[dict[str, Any]] = []
        temporary_hrefs: list[str] = []
        events = network_events if network_events is not None else []

        def _finish(
            *,
            status: str,
            success: bool,
            gene_query_submitted: bool,
            evidence: dict[str, Any],
            error_type: str | None = None,
            error_message: str | None = None,
        ) -> ToolResult:
            acceptance = acceptance_from_evidence(
                evidence, gene_query_submitted=gene_query_submitted
            )
            acceptance["gene_specific_plot_count"] = max(
                acceptance.get("gene_specific_plot_count", 0),
                sum(1 for a in artifacts if a.get("kind") == "plot_image"),
            )
            manifest = {
                "gene_symbol": mouse_gene_symbol,
                "mode": "dynamic_genex",
                "base_url": BASE_URL,
                "final_url": str(getattr(page, "url", BASE_URL) or BASE_URL),
                "status": status,
                "acceptance": acceptance,
                "shiny_ready": evidence,
                "artifacts": artifacts,
                "downloads": downloads,
                "extractions": extractions,
            }
            manifest_path = out_dir / "manifest.json"
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            network_path = out_dir / "network.json"
            network_path.write_text(
                json.dumps({"events": events[:500]}, indent=2) + "\n",
                encoding="utf-8",
            )
            diagnostics_path = out_dir / "query_diagnostics.json"
            diagnostics_path.write_text(
                json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
            )

            captured_views = {
                a.get("view") for a in artifacts if a.get("kind") == "plot_image"
            }
            if len(captured_views) > 1:
                view_type = VIEW_MIXED
            elif captured_views == {"rank"}:
                view_type = VIEW_RANK
            elif captured_views == {"tsne_global"}:
                view_type = VIEW_GLOBAL_TSNE
            elif captured_views == {"tsne_local"}:
                view_type = VIEW_REGIONAL_TSNE
            else:
                view_type = VIEW_UNKNOWN

            rank_extraction_status = next(
                (
                    e.get("status", "not_attempted")
                    for e in extractions
                    if e.get("view") == "rank"
                ),
                None,
            )
            if rank_extraction_status is None:
                rank_download = next(
                    (d for d in downloads if d.get("view") == "rank"), None
                )
                if rank_download and not rank_download.get("success"):
                    # The rank table lives only inside the ZIP the server refused
                    # to stream, so extraction never had an input.
                    rank_extraction_status = "source_export_unavailable"
                else:
                    rank_extraction_status = "not_attempted"

            payload = {
                "gene_symbol": mouse_gene_symbol,
                "state_url": BASE_URL,
                "view_type": view_type,
                "artifacts": artifacts,
                "downloads": downloads,
                "extractions": extractions,
                "acquisition_status": (
                    "success" if any(d.get("success") for d in downloads) else status
                ),
                "rank_extraction_status": rank_extraction_status,
                "regional_quantitative_status": "unavailable",
                "regional_evidence_type": "figure_derived",
                "acceptance": acceptance,
                "dynamic_ui": True,
            }
            audit = {
                "browser_channel": browser_channel,
                "launch_attempts": launch_attempts or [],
                "final_url": manifest["final_url"],
                "temporary_download_hrefs": temporary_hrefs,
                "network_manifest_path": str(network_path),
                "manifest_path": str(manifest_path),
                "query_diagnostics_path": str(diagnostics_path),
                "session_href_is_stable_endpoint": False,
            }
            return _tool_result(
                endpoint_name=endpoint,
                gene_symbol=mouse_gene_symbol,
                request_url=BASE_URL,
                request_params={"mouse_gene_symbol": mouse_gene_symbol},
                success=success,
                status_code=None,
                error_type=error_type,
                error_message=error_message,
                data=_envelope(status=status, payload=payload, audit=audit),
            )

        try:
            page.goto(
                BASE_URL, wait_until="domcontentloaded", timeout=self.navigation_timeout_ms
            )
            self._settle(page, 3_000)

            if not self._open_query_tab(page, diagnostics):
                return _finish(
                    status="query_ui_unavailable",
                    success=False,
                    gene_query_submitted=False,
                    evidence=evaluate_shiny_ready(page, mouse_gene_symbol),
                    error_type="query_ui_unavailable",
                    error_message="Could not reach the DropViz Query gene input",
                )

            if not self._select_gene(page, mouse_gene_symbol, diagnostics):
                return _finish(
                    status="gene_selection_failed",
                    success=False,
                    gene_query_submitted=False,
                    evidence=evaluate_shiny_ready(page, mouse_gene_symbol),
                    error_type="gene_selection_failed",
                    error_message=f"Could not commit {mouse_gene_symbol} in the gene control",
                )

            submitted = self._submit_query(page, diagnostics)

            view_reports: list[dict[str, Any]] = []

            def _capture_view(view: dict[str, str]) -> None:
                report: dict[str, Any] = {"view": view["key"]}
                report.update(self._activate_view(page, view))
                image_info = self._wait_for_plot_image(
                    page, view["image_selector"], plot_timeout_ms
                )
                report["image"] = image_info

                if image_info.get("rendered"):
                    image_path = images_dir / view["image_name"]
                    try:
                        page.locator(view["image_selector"]).first.screenshot(
                            path=str(image_path)
                        )
                        if image_path.is_file() and image_path.stat().st_size > 0:
                            artifacts.append(
                                {
                                    "kind": "plot_image",
                                    "view": view["key"],
                                    "path": str(image_path),
                                    "sha256": sha256_file(image_path),
                                    "byte_size": image_path.stat().st_size,
                                    "width": image_info.get("width"),
                                    "height": image_info.get("height"),
                                }
                            )
                            report["screenshot_status"] = "success"
                        else:
                            report["screenshot_status"] = "empty"
                    except Exception as exc:  # noqa: BLE001
                        report["screenshot_status"] = "failed"
                        report["screenshot_error"] = str(exc)[:200]

                    download_path = downloads_dir / view["download_name"]
                    dl_result = self.download_shiny_export(
                        page=page,
                        link_id=view["download_id"],
                        output_path=download_path,
                        gene_symbol=mouse_gene_symbol,
                    )
                    dl_data = dl_result.data if isinstance(dl_result.data, dict) else {}
                    dl_payload = dl_data.get("payload") or {}
                    dl_audit = dl_data.get("audit") or {}
                    for href in dl_audit.get("temporary_download_hrefs") or []:
                        if href not in temporary_hrefs:
                            temporary_hrefs.append(href)
                    downloads.append(
                        {
                            "view": view["key"],
                            "link_id": view["download_id"],
                            "success": bool(dl_result.success),
                            "status": dl_data.get("status"),
                            "error_message": dl_result.error_message,
                            "capture_via": (dl_result.request_params or {}).get(
                                "capture_via"
                            ),
                            "path": dl_payload.get("path"),
                            "sha256": dl_payload.get("sha256"),
                            "byte_size": dl_payload.get("byte_size"),
                            "zip_members": dl_payload.get("zip_members"),
                        }
                    )
                    report["download_status"] = dl_data.get("status")

                    if dl_result.success and not skip_extraction:
                        target = extracted_dir / view["extract_dir"]
                        if view["key"] == "rank":
                            extraction = run_extract_dropviz_rank(
                                zip_or_rdata_path=download_path, output_dir=target
                            )
                        else:
                            extraction = run_inspect_dropviz_rdata(
                                zip_or_rdata_path=download_path, output_dir=target
                            )
                        extraction["view"] = view["key"]
                        extraction["output_dir"] = str(target)
                        extractions.append(extraction)
                        report["extraction_status"] = extraction.get("status")
                else:
                    report["download_status"] = "plot_unavailable"
                view_reports.append(report)

            for view in DYNAMIC_VIEW_PLAN:
                if not view.get("requires_region"):
                    _capture_view(view)

            table_reports: list[dict[str, Any]] = []
            rendered_top_rows: dict[str, dict[str, Any]] = {}
            for export in DYNAMIC_TABLE_EXPORTS:
                table_report: dict[str, Any] = {"export": export["key"]}
                table_report.update(self._activate_view(page, export))
                rendered = self.read_top_cluster_row(page)
                if rendered:
                    rendered_top_rows[export["key"]] = rendered
                    table_report["rendered_top_row"] = rendered
                csv_path = downloads_dir / export["download_name"]
                # The table export is the only structured quantitative artifact
                # the live server still streams, and it drops connections under
                # load, so a dropped export is worth one retry.
                for attempt in range(2):
                    csv_result = self.download_shiny_export(
                        page=page,
                        link_id=export["download_id"],
                        output_path=csv_path,
                        gene_symbol=mouse_gene_symbol,
                        expected_kind="csv",
                    )
                    if csv_result.success:
                        break
                    table_report["retried"] = True
                    self._settle(page, 2_000)
                if table_report.get("retried"):
                    table_report["attempts"] = attempt + 1
                csv_data = csv_result.data if isinstance(csv_result.data, dict) else {}
                csv_payload = csv_data.get("payload") or {}
                for href in (csv_data.get("audit") or {}).get(
                    "temporary_download_hrefs"
                ) or []:
                    if href not in temporary_hrefs:
                        temporary_hrefs.append(href)
                downloads.append(
                    {
                        "view": export["key"],
                        "link_id": export["download_id"],
                        "success": bool(csv_result.success),
                        "status": csv_data.get("status"),
                        "error_message": csv_result.error_message,
                        "content_kind": "csv",
                        "path": csv_payload.get("path"),
                        "sha256": csv_payload.get("sha256"),
                        "byte_size": csv_payload.get("byte_size"),
                    }
                )
                if csv_result.success:
                    summary = summarize_cluster_table_csv(csv_path, mouse_gene_symbol)
                    summary_path = (
                        extracted_dir / export["key"] / "cluster_table_summary.json"
                    )
                    summary_path.parent.mkdir(parents=True, exist_ok=True)
                    summary_path.write_text(
                        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
                    )
                    summary["view"] = export["key"]
                    summary["output_dir"] = str(summary_path.parent)
                    extractions.append(summary)
                table_report["status"] = csv_data.get("status")
                table_reports.append(table_report)
            diagnostics["table_exports"] = table_reports

            # The local/regional t-SNE renders a "choose a region" placeholder
            # until a region filter is committed, so pick the region of the
            # gene's top-ranked cluster and re-run the query.
            region_phase: dict[str, Any] = {"attempted": False}
            top_region = next(
                (
                    (e.get("top_clusters") or [{}])[0].get("region")
                    for e in extractions
                    if e.get("view") == "clusters_table" and e.get("ok")
                ),
                None,
            )
            region_source = "clusters_table_csv" if top_region else None
            if not top_region:
                rendered = rendered_top_rows.get("clusters_table") or rendered_top_rows.get(
                    "subclusters_table"
                )
                if rendered:
                    top_region = rendered.get("region")
                    region_source = "rendered_cluster_table"
            region_views = [v for v in DYNAMIC_VIEW_PLAN if v.get("requires_region")]
            if top_region and region_views:
                region_phase = {
                    "attempted": True,
                    "region": top_region,
                    "region_source": region_source,
                }
                self._activate_view(page, {"main_tab": "clusters"} | {
                    "sub_tab_container": "clusterpanel",
                    "sub_tab": "rank",
                })
                selected = self._select_region(page, str(top_region), diagnostics)
                region_phase["region_selected"] = selected
                if selected:
                    region_phase["resubmitted"] = self._submit_query(page, diagnostics)
                    for view in region_views:
                        _capture_view(view)
                else:
                    for view in region_views:
                        view_reports.append(
                            {
                                "view": view["key"],
                                "skipped": "region_selection_failed",
                                "download_status": "requires_region_selection",
                            }
                        )
            else:
                for view in region_views:
                    view_reports.append(
                        {
                            "view": view["key"],
                            "skipped": "no_region_available",
                            "download_status": "requires_region_selection",
                        }
                    )
            diagnostics["region_phase"] = region_phase
            diagnostics["views"] = view_reports

            try:
                page.screenshot(path=str(out_dir / "page_screenshot.png"), full_page=True)
                shot = out_dir / "page_screenshot.png"
                if shot.is_file() and shot.stat().st_size > 0:
                    artifacts.append(
                        {
                            "kind": "screenshot",
                            "path": str(shot),
                            "sha256": sha256_file(shot),
                            "byte_size": shot.stat().st_size,
                        }
                    )
            except Exception as exc:  # noqa: BLE001
                diagnostics["page_screenshot_error"] = str(exc)[:200]

            evidence = evaluate_shiny_ready(page, mouse_gene_symbol)
            plot_count = sum(1 for a in artifacts if a.get("kind") == "plot_image")
            expected = len(DYNAMIC_VIEW_PLAN) + len(DYNAMIC_TABLE_EXPORTS)
            download_ok = sum(1 for d in downloads if d.get("success"))

            if plot_count == 0:
                status = "plot_unavailable"
                success = False
            elif download_ok == expected and all(e.get("ok") for e in extractions):
                status = "success"
                success = True
            else:
                status = "partial_success"
                success = True

            return _finish(
                status=status,
                success=success,
                gene_query_submitted=submitted,
                evidence=evidence,
                error_type=None if success else status,
                error_message=None if success else f"collect_dynamic_query: {status}",
            )
        except Exception as exc:  # noqa: BLE001
            diagnostics["fatal_error"] = str(exc)[:500]
            return _finish(
                status="dynamic_ui_failed",
                success=False,
                gene_query_submitted=bool(diagnostics.get("query_submitted")),
                evidence={"ready": False, "reason": "exception"},
                error_type="dynamic_ui_failed",
                error_message=str(exc)[:500],
            )

    def drive_genex_dynamic_ui(
        self,
        *,
        page: Any,
        mouse_gene_symbol: str,
        output_dir: Path | str,
        network_events: list[dict[str, Any]] | None = None,
        browser_channel: str | None = None,
        launch_attempts: list[dict[str, Any]] | None = None,
        skip_extraction: bool = False,
    ) -> ToolResult:
        """Backwards-compatible alias for :meth:`collect_dynamic_query`."""
        return self.collect_dynamic_query(
            page=page,
            mouse_gene_symbol=mouse_gene_symbol,
            output_dir=output_dir,
            network_events=network_events,
            browser_channel=browser_channel,
            launch_attempts=launch_attempts,
            skip_extraction=skip_extraction,
        )


def _run_dynamic_query(
    *,
    client: "DropVizClient",
    page: Any,
    mouse_gene_symbol: str,
    out_dir: Path,
    channel: str,
    attempts: list[dict[str, Any]],
    skip_extraction: bool,
    audit: dict[str, Any],
    reason: str,
) -> dict[str, Any]:
    """Run the dynamic Query workflow and fold it into a state record."""
    audit["used_dynamic_fallback"] = True
    audit["fallback_reason"] = reason
    dyn_dir = out_dir / "dynamic_genex"
    result = client.collect_dynamic_query(
        page=page,
        mouse_gene_symbol=mouse_gene_symbol,
        output_dir=dyn_dir,
        browser_channel=channel,
        launch_attempts=attempts,
        skip_extraction=skip_extraction,
    )
    data = result.data if isinstance(result.data, dict) else {}
    payload = data.get("payload") or {}
    dyn_audit = data.get("audit") or {}
    for href in dyn_audit.get("temporary_download_hrefs") or []:
        if href not in audit.setdefault("temporary_download_hrefs", []):
            audit["temporary_download_hrefs"].append(href)
    return {
        "state_url": BASE_URL,
        "state_dir": str(dyn_dir),
        "status": data.get("status"),
        "view_type": payload.get("view_type"),
        "acquisition_status": payload.get("acquisition_status"),
        "rank_extraction_status": payload.get("rank_extraction_status", "not_attempted"),
        "regional_quantitative_status": payload.get(
            "regional_quantitative_status", "unavailable"
        ),
        "regional_evidence_type": payload.get("regional_evidence_type", "figure_derived"),
        "acceptance": payload.get("acceptance"),
        "artifacts": payload.get("artifacts") or [],
        "downloads": payload.get("downloads") or [],
        "extractions": payload.get("extractions") or [],
        "dynamic_ui": True,
    }


def _aggregate_acceptance(state_results: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Report gene-level acceptance for the attempt Section 2c would consume.

    Section 2c reads exactly one attempt, so the gene-level block mirrors that
    attempt rather than OR-ing across attempts: an expired saved state would
    otherwise stamp ``restore_error_present`` onto an otherwise clean dynamic
    collection. Attempts that were tried and rejected stay visible in
    ``rejected_attempts`` so the discarded evidence is still auditable.
    """
    attempts = [r for r in state_results if isinstance(r.get("acceptance"), dict)]
    rejected = [
        {
            "state_url": r.get("state_url"),
            "status": r.get("status"),
            "restore_error_present": bool(
                (r.get("acceptance") or {}).get("restore_error_present")
            ),
        }
        for r in attempts
        if not (r.get("acceptance") or {}).get("usable_for_section_2c")
    ]
    accepted = next(
        (r for r in attempts if (r["acceptance"]).get("usable_for_section_2c")),
        None,
    )
    if accepted is None:
        return {
            "requested_gene_visible": False,
            "homepage_only": True,
            "restore_error_present": any(a["restore_error_present"] for a in rejected),
            "gene_query_submitted": False,
            "gene_specific_plot_count": 0,
            "enabled_download_count": 0,
            "usable_for_section_2c": False,
            "accepted_attempt": None,
            "rejected_attempts": rejected,
        }
    return dict(accepted["acceptance"]) | {
        "accepted_attempt": accepted.get("state_url"),
        "rejected_attempts": rejected,
    }


def _merge_status(statuses: Sequence[str]) -> str:
    """Combine per-state statuses into an overall gene status."""
    if not statuses:
        return "no_states"
    if all(s == "success" for s in statuses):
        return "success"
    if any(s == "success" or s == "partial_success" for s in statuses):
        return "partial_success"
    # Prefer specific failure labels when unanimous.
    unique = set(statuses)
    if len(unique) == 1:
        return next(iter(unique))
    return "partial_success"


def collect_dropviz_gene(
    *,
    mouse_gene_symbol: str,
    output_dir: Path | str,
    saved_state_urls: Sequence[str] = (),
    client: DropVizClient | None = None,
    skip_extraction: bool = False,
) -> ToolResult:
    """Collect DropViz artifacts for one mouse gene. Never raises.

    Dynamic base-site fallback runs only when no URLs were supplied or every
    saved-state attempt fails with a state-availability status (correction #3).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    client = client or DropVizClient()
    endpoint = "collect_dropviz_gene"
    urls = [u for u in (saved_state_urls or []) if str(u).strip()]
    state_results: list[dict[str, Any]] = []
    overall_artifacts: list[dict[str, Any]] = []
    overall_downloads: list[dict[str, Any]] = []
    overall_extractions: list[dict[str, Any]] = []
    audit: dict[str, Any] = {
        "browser_channel": None,
        "launch_attempts": [],
        "final_url": None,
        "temporary_download_hrefs": [],
        "network_manifest_path": None,
        "used_dynamic_fallback": False,
        "fallback_reason": None,
    }

    def _process_state(
        *,
        state_url: str,
        state_dir: Path,
        page: Any | None = None,
        channel: str | None = None,
        attempts: list[dict[str, Any]] | None = None,
    ) -> dict[str, Any]:
        inspect = client.inspect_saved_state(
            gene_symbol=mouse_gene_symbol,
            state_url=state_url,
            output_dir=state_dir,
            page=page,
            browser_channel=channel,
            launch_attempts=attempts,
        )
        data = inspect.data if isinstance(inspect.data, dict) else {}
        payload = dict(data.get("payload") or {})
        state_audit = dict(data.get("audit") or {})
        status = str(data.get("status") or ("success" if inspect.success else "failed"))

        # Merge audit (temporary hrefs only here).
        for href in state_audit.get("temporary_download_hrefs") or []:
            if href not in audit["temporary_download_hrefs"]:
                audit["temporary_download_hrefs"].append(href)
        if state_audit.get("browser_channel"):
            audit["browser_channel"] = state_audit["browser_channel"]
        if state_audit.get("launch_attempts"):
            audit["launch_attempts"] = state_audit["launch_attempts"]
        if state_audit.get("final_url"):
            audit["final_url"] = state_audit["final_url"]
        if state_audit.get("network_manifest_path"):
            audit["network_manifest_path"] = state_audit["network_manifest_path"]

        downloads: list[dict[str, Any]] = []
        extractions: list[dict[str, Any]] = []
        acquisition_status = payload.get("acquisition_status") or status
        rank_extraction_status = "not_attempted"
        regional_quantitative_status = "unavailable"
        regional_evidence_type = "figure_derived"
        view_type = payload.get("view_type") or VIEW_UNKNOWN

        # Download exports when acquisition succeeded and we have a live page.
        if inspect.success and page is not None:
            link_ids = list(payload.get("download_link_ids") or [])
            rank_id = find_rank_download_id(link_ids)
            tsne_ids = find_tsne_download_ids(link_ids)
            targets: list[tuple[str, str]] = []
            if rank_id:
                targets.append((rank_id, "rank.zip"))
            for tid in tsne_ids:
                safe = tid.replace(".", "_")
                targets.append((tid, f"{safe}.zip"))

            for link_id, filename in targets:
                dl_path = state_dir / filename
                dl_result = client.download_shiny_export(
                    page=page,
                    link_id=link_id,
                    output_path=dl_path,
                    gene_symbol=mouse_gene_symbol,
                )
                dl_data = dl_result.data if isinstance(dl_result.data, dict) else {}
                dl_payload = dict(dl_data.get("payload") or {})
                dl_audit = dict(dl_data.get("audit") or {})
                for href in dl_audit.get("temporary_download_hrefs") or []:
                    if href not in audit["temporary_download_hrefs"]:
                        audit["temporary_download_hrefs"].append(href)
                if dl_result.success:
                    for rec in dl_payload.get("downloads") or []:
                        downloads.append(rec)
                        overall_downloads.append(rec)
                else:
                    downloads.append(
                        {
                            "link_id": link_id,
                            "success": False,
                            "error_type": dl_result.error_type,
                            "error_message": dl_result.error_message,
                        }
                    )

            if not skip_extraction:
                for rec in list(downloads):
                    path = rec.get("path")
                    if not path or not rec.get("sha256"):
                        continue
                    path_obj = Path(str(path))
                    name = path_obj.name.lower()
                    extract_dir = state_dir / f"extract_{path_obj.stem}"
                    if "rank" in name or (
                        rank_id and rec.get("link_id") == rank_id
                    ):
                        extracted = run_extract_dropviz_rank(
                            zip_or_rdata_path=path_obj,
                            output_dir=extract_dir,
                        )
                        rank_extraction_status = str(
                            extracted.get("status") or "extraction_failed"
                        )
                        extractions.append(
                            {
                                "kind": "rank",
                                "source_zip": str(path_obj),
                                "parent_sha256": rec.get("sha256"),
                                "result": extracted,
                                "api_run": None,
                            }
                        )
                        overall_extractions.append(extractions[-1])
                    else:
                        extracted = run_inspect_dropviz_rdata(
                            zip_or_rdata_path=path_obj,
                            output_dir=extract_dir,
                        )
                        regional = extracted.get("regional") or {}
                        if regional:
                            regional_quantitative_status = str(
                                regional.get("regional_quantitative_status")
                                or "unavailable"
                            )
                            regional_evidence_type = str(
                                regional.get("regional_evidence_type")
                                or "figure_derived"
                            )
                        extractions.append(
                            {
                                "kind": "rdata_inspect",
                                "source_zip": str(path_obj),
                                "parent_sha256": rec.get("sha256"),
                                "result": extracted,
                                "api_run": None,
                            }
                        )
                        overall_extractions.append(extractions[-1])

        # Overall state status: acquisition success + extraction partial.
        state_status = status
        if acquisition_status == "success":
            if rank_extraction_status in {
                "rscript_unavailable",
                "extraction_failed",
                "missing_clusters_top",
                "rank_validation_failed",
            }:
                state_status = "partial_success"
            elif downloads and any(d.get("success") is False for d in downloads):
                state_status = "partial_success"
            else:
                state_status = "success"

        for art in payload.get("artifacts") or []:
            overall_artifacts.append(art)

        record = {
            "state_url": state_url,
            "state_dir": str(state_dir),
            "status": state_status,
            "view_type": view_type,
            "acquisition_status": acquisition_status,
            "rank_extraction_status": rank_extraction_status,
            "regional_quantitative_status": regional_quantitative_status,
            "regional_evidence_type": regional_evidence_type,
            "acceptance": payload.get("acceptance"),
            "artifacts": payload.get("artifacts") or [],
            "downloads": downloads,
            "extractions": extractions,
            "inspect_success": inspect.success,
            "error_type": inspect.error_type,
        }
        return record

    # --- Browser session spanning all states ---
    try:
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name=endpoint,
            gene_symbol=mouse_gene_symbol,
            request_url=urls[0] if urls else BASE_URL,
            request_params={
                "mouse_gene_symbol": mouse_gene_symbol,
                "saved_state_urls": list(urls),
            },
            success=False,
            error_type="playwright_unavailable",
            error_message=str(exc)[:500],
            data=_envelope(
                status="playwright_unavailable",
                payload={
                    "gene_symbol": mouse_gene_symbol,
                    "state_url": None,
                    "view_type": VIEW_UNKNOWN,
                    "artifacts": [],
                    "downloads": [],
                    "extractions": [],
                    "acquisition_status": "failed",
                    "rank_extraction_status": "not_attempted",
                    "regional_quantitative_status": "unavailable",
                    "regional_evidence_type": "figure_derived",
                    "states": [],
                },
                audit=audit,
            ),
        )

    try:
        with sync_playwright() as pw:
            try:
                browser, attempts, channel = client._launch_browser(pw)
            except Exception as launch_exc:  # noqa: BLE001
                return _tool_result(
                    endpoint_name=endpoint,
                    gene_symbol=mouse_gene_symbol,
                    request_url=urls[0] if urls else BASE_URL,
                    request_params={
                        "mouse_gene_symbol": mouse_gene_symbol,
                        "saved_state_urls": list(urls),
                    },
                    success=False,
                    error_type="browser_launch_failed",
                    error_message=str(launch_exc)[:500],
                    data=_envelope(
                        status="browser_launch_failed",
                        payload={
                            "gene_symbol": mouse_gene_symbol,
                            "state_url": None,
                            "view_type": VIEW_UNKNOWN,
                            "artifacts": [],
                            "downloads": [],
                            "extractions": [],
                            "acquisition_status": "failed",
                            "rank_extraction_status": "not_attempted",
                            "regional_quantitative_status": "unavailable",
                            "regional_evidence_type": "figure_derived",
                            "states": [],
                        },
                        audit={**audit, "launch_attempts": attempts if "attempts" in dir() else []},
                    ),
                )
            audit["browser_channel"] = channel
            audit["launch_attempts"] = attempts
            context = browser.new_context(accept_downloads=True)
            page = context.new_page()

            if urls:
                for idx, state_url in enumerate(urls, start=1):
                    state_dir = out_dir / f"state_{idx}"
                    state_dir.mkdir(parents=True, exist_ok=True)
                    record = _process_state(
                        state_url=state_url,
                        state_dir=state_dir,
                        page=page,
                        channel=channel,
                        attempts=attempts,
                    )
                    state_results.append(record)

                # Dynamic fallback ONLY for state availability failures / no URLs.
                all_state_failures = state_results and all(
                    (r.get("acquisition_status") in _STATE_FAILURE_STATUSES)
                    or (r.get("status") in _STATE_FAILURE_STATUSES)
                    for r in state_results
                )
                if all_state_failures:
                    state_results.append(
                        _run_dynamic_query(
                            client=client,
                            page=page,
                            mouse_gene_symbol=mouse_gene_symbol,
                            out_dir=out_dir,
                            channel=channel,
                            attempts=attempts,
                            skip_extraction=skip_extraction,
                            audit=audit,
                            reason="all_saved_states_unavailable",
                        )
                    )
            else:
                state_results.append(
                    _run_dynamic_query(
                        client=client,
                        page=page,
                        mouse_gene_symbol=mouse_gene_symbol,
                        out_dir=out_dir,
                        channel=channel,
                        attempts=attempts,
                        skip_extraction=skip_extraction,
                        audit=audit,
                        reason="no_saved_state_urls",
                    )
                )

            try:
                context.close()
            except Exception:  # noqa: BLE001
                pass
            try:
                browser.close()
            except Exception:  # noqa: BLE001
                pass
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name=endpoint,
            gene_symbol=mouse_gene_symbol,
            request_url=urls[0] if urls else BASE_URL,
            request_params={
                "mouse_gene_symbol": mouse_gene_symbol,
                "saved_state_urls": list(urls),
            },
            success=False,
            error_type="collect_failed",
            error_message=str(exc)[:500],
            data=_envelope(
                status="collect_failed",
                payload={
                    "gene_symbol": mouse_gene_symbol,
                    "state_url": None,
                    "view_type": VIEW_UNKNOWN,
                    "artifacts": overall_artifacts,
                    "downloads": overall_downloads,
                    "extractions": overall_extractions,
                    "acquisition_status": "failed",
                    "rank_extraction_status": "not_attempted",
                    "regional_quantitative_status": "unavailable",
                    "regional_evidence_type": "figure_derived",
                    "states": state_results,
                },
                audit=audit,
            ),
        )

    # Dynamic-query records carry their own artifacts/downloads/extractions.
    for record in state_results:
        if not record.get("dynamic_ui"):
            continue
        overall_artifacts.extend(record.get("artifacts") or [])
        overall_downloads.extend(record.get("downloads") or [])
        overall_extractions.extend(record.get("extractions") or [])

    statuses = [str(r.get("status") or "failed") for r in state_results]
    overall = _merge_status(statuses)

    # Aggregate acquisition / extraction flags.
    acq_statuses = [str(r.get("acquisition_status") or "") for r in state_results]
    if any(a == "success" for a in acq_statuses):
        acquisition_status = "success"
    elif acq_statuses and all(a in _STATE_FAILURE_STATUSES for a in acq_statuses):
        acquisition_status = acq_statuses[0]
    else:
        acquisition_status = acq_statuses[0] if acq_statuses else "failed"

    rank_statuses = [
        str(r.get("rank_extraction_status") or "not_attempted") for r in state_results
    ]
    if any(s == "success" for s in rank_statuses):
        rank_extraction_status = "success"
    elif any(s == "rscript_unavailable" for s in rank_statuses):
        rank_extraction_status = "rscript_unavailable"
        if acquisition_status == "success" and overall == "success":
            overall = "partial_success"
    elif any(
        s
        in {
            "extraction_failed",
            "missing_clusters_top",
            "rank_validation_failed",
            "source_export_unavailable",
        }
        for s in rank_statuses
    ):
        rank_extraction_status = next(
            s
            for s in rank_statuses
            if s
            in {
                "extraction_failed",
                "missing_clusters_top",
                "rank_validation_failed",
                "source_export_unavailable",
                "rscript_unavailable",
            }
        )
        if acquisition_status == "success":
            overall = "partial_success"
    else:
        rank_extraction_status = "not_attempted"

    regional_statuses = [
        str(r.get("regional_quantitative_status") or "unavailable") for r in state_results
    ]
    regional_quantitative_status = (
        "available" if any(s == "available" for s in regional_statuses) else "unavailable"
    )
    regional_evidence_type = (
        "table_derived"
        if regional_quantitative_status == "available"
        else "figure_derived"
    )

    # If acquisition ok but R missing → partial_success (correction #1).
    if acquisition_status == "success" and rank_extraction_status == "rscript_unavailable":
        overall = "partial_success"

    # Prefer a view named by an attempt that actually rendered the gene.
    view_types = [
        r.get("view_type")
        for r in state_results
        if r.get("view_type") and r.get("view_type") != VIEW_UNKNOWN
    ] or [r.get("view_type") for r in state_results if r.get("view_type")]
    primary_view = view_types[0] if view_types else VIEW_UNKNOWN

    acceptance = _aggregate_acceptance(state_results)
    rank_download = next(
        (d for d in overall_downloads if d.get("view") == "rank"), None
    )
    summary = {
        "gene_symbol": mouse_gene_symbol,
        "status": overall,
        "acquisition_status": acquisition_status,
        "acceptance": acceptance,
        "rank_download_status": (
            rank_download.get("status") if rank_download else "not_attempted"
        ),
        "rank_extraction_status": rank_extraction_status,
        "regional_quantitative_status": regional_quantitative_status,
        "regional_evidence_type": regional_evidence_type,
        "states": state_results,
        "artifact_count": len(overall_artifacts),
        "download_count": len([d for d in overall_downloads if d.get("sha256")]),
        "extraction_count": len(overall_extractions),
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    overall_artifacts.append(
        {
            "kind": "gene_summary",
            "path": str(summary_path),
            "sha256": sha256_file(summary_path),
        }
    )

    success = overall in {"success", "partial_success"}
    return _tool_result(
        endpoint_name=endpoint,
        gene_symbol=mouse_gene_symbol,
        request_url=urls[0] if urls else BASE_URL,
        request_params={
            "mouse_gene_symbol": mouse_gene_symbol,
            "saved_state_urls": list(urls),
            "output_dir": str(out_dir),
        },
        success=success,
        status_code=None,
        error_type=None if success else overall,
        error_message=None if success else f"collect_dropviz_gene: {overall}",
        data=_envelope(
            status=overall,
            payload={
                "gene_symbol": mouse_gene_symbol,
                "state_url": urls[0] if urls else None,
                "view_type": primary_view,
                "artifacts": overall_artifacts,
                "downloads": overall_downloads,
                "extractions": overall_extractions,
                "acquisition_status": acquisition_status,
                "rank_extraction_status": rank_extraction_status,
                "regional_quantitative_status": regional_quantitative_status,
                "regional_evidence_type": regional_evidence_type,
                "states": state_results,
                "summary_path": str(summary_path),
            },
            audit=audit,
        ),
    )


__all__ = [
    "ALLOWED_HOSTS",
    "BASE_URL",
    "CANDIDATE_RANK_LINK_IDS",
    "CANDIDATE_TSNE_LINK_IDS",
    "DropVizClient",
    "RANK_SORT_POLICY",
    "SOURCE_NAME",
    "VIEW_GLOBAL_TSNE",
    "VIEW_MIXED",
    "VIEW_RANK",
    "VIEW_REGIONAL_TSNE",
    "VIEW_UNKNOWN",
    "DYNAMIC_TABLE_EXPORTS",
    "DYNAMIC_VIEW_PLAN",
    "HOMEPAGE_ASSET_MARKERS",
    "acceptance_from_evidence",
    "assess_regional_quantitative_fields",
    "classify_view",
    "collect_dropviz_gene",
    "css_escape_id",
    "derive_rank_outputs_from_raw_csv",
    "detect_restore_error",
    "detect_state_failure",
    "evaluate_shiny_ready",
    "summarize_cluster_table_csv",
    "download_selector",
    "find_rank_download_id",
    "find_tsne_download_ids",
    "inventory_zip_basenames",
    "is_allowed_dropviz_url",
    "is_html_payload",
    "is_zip_bytes",
    "normalize_dropviz_url",
    "redirect_is_allowed",
    "rscript_available",
    "run_extract_dropviz_rank",
    "run_inspect_dropviz_rdata",
    "sha256_bytes",
    "sha256_file",
    "validate_rank_row",
]
