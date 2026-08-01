"""Bundle-only Section 1e homologues / ortholog helpers.

Owns NCBI Datasets ortholog pagination, taxonomy scope filtering, OrthoDB
supporting lookup, and official NCBI Orthologs table capture for the section
bundle. Not wired into the full dossier workflow's generic Datasets client path
when Section 1e is selected.
"""

from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Sequence
from urllib.parse import quote, urlparse

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import (
    ApiRun,
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceCoverageResult,
    SourceStatus,
    SourceType,
)
from gene_dossier.section_1c import (
    _append_evidence,
    _persist_artifact_bytes,
    _persist_tool_result_json,
    _save_api_run_failure,
    _tool_result_to_api_run,
    _validate_nonblank_image,
)
from gene_dossier.source_ids import make_source_id, slugify
from gene_dossier.tools import ncbi_datasets, ncbi_taxonomy, orthodb
from gene_dossier.workflow import DossierState, WorkflowTransientContext

logger = logging.getLogger(__name__)

SECTION_HOMOLOGUES = "Homologues"
SUBSECTION_1E = "Homologues in model animals"

SUPPORTED_SECTION_1E_SCOPES: dict[int, dict[str, Any]] = {
    7776: {
        "label": "jawed vertebrates",
        "ncbi_scope_tax_id": 7776,
        "ncbi_taxon_name": "Gnathostomata",
    },
    32523: {
        "label": "tetrapods",
        "ncbi_scope_tax_id": 32523,
        "ncbi_taxon_name": "Tetrapoda",
    },
}

MODEL_SPECIES_PRIORITY: tuple[int, ...] = (
    10090,
    10116,
    7955,
    9031,
    8364,
    9823,
    9913,
    9544,
    9598,
    9615,
)

TABLE_OFFICIAL = "official_capture"
TABLE_COMPLETE_FALLBACK = "complete_api_fallback"
TABLE_PARTIAL_FALLBACK = "partial_api_fallback"
TABLE_UNAVAILABLE = "unavailable"

REQUIRED_HEADERS: tuple[str, ...] = (
    "Scientific name",
    "Symbol",
    "Length (aa)",
    "Architecture",
)
CAPTURE_ORIGIN_LIVE = "live_ncbi_datasets_dom"

_MIN_CAPTURE_WIDTH = 500
_MIN_CAPTURE_HEIGHT = 250
_MIN_CAPTURE_BYTES = 15_000
_NCBI_ALLOWED_HOSTS = frozenset({"www.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov"})


@dataclass(frozen=True)
class Section1eConfig:
    ortholog_scope_tax_id: int = 7776
    max_visible_rows: int = 20

    def __post_init__(self) -> None:
        if int(self.ortholog_scope_tax_id) not in SUPPORTED_SECTION_1E_SCOPES:
            raise ValueError(
                f"Unsupported Section 1e scope tax id: {self.ortholog_scope_tax_id}"
            )
        if int(self.max_visible_rows) < 1:
            raise ValueError("max_visible_rows must be >= 1")

    @property
    def ortholog_scope_label(self) -> str:
        return str(SUPPORTED_SECTION_1E_SCOPES[int(self.ortholog_scope_tax_id)]["label"])

    @property
    def ncbi_taxon_name(self) -> str:
        return str(
            SUPPORTED_SECTION_1E_SCOPES[int(self.ortholog_scope_tax_id)]["ncbi_taxon_name"]
        )


def ortholog_ncbi_url(*, entrez_gene_id: str | int) -> str:
    gene_id = str(entrez_gene_id).strip()
    if not gene_id.isdigit():
        raise ValueError("A numeric Entrez Gene ID is required")
    return f"https://www.ncbi.nlm.nih.gov/datasets/gene/{gene_id}/#orthologs"


def ortholog_ncbi_legacy_link(
    *,
    entrez_gene_id: str | int,
    scope_tax_id: int,
    gene_symbol: str,
) -> str:
    """Audit/reference-only legacy Gene Orthologs URL. Never a Playwright target."""
    gene_id = str(entrez_gene_id).strip()
    symbol = (gene_symbol or "").strip() or gene_id
    return (
        f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}/ortholog/"
        f"?scope={int(scope_tax_id)}&term={quote(symbol)}"
    )


def _normalize_header(text: str) -> str:
    return re.sub(r"\s+", " ", (text or "").strip()).lower()


def _headers_cover_required(visible_headers: Sequence[str]) -> bool:
    visible = {_normalize_header(h) for h in visible_headers if str(h).strip()}
    required = {_normalize_header(h) for h in REQUIRED_HEADERS}
    return required <= visible


def _parse_row_cells_by_header(
    headers: Sequence[str],
    cells: Sequence[str],
) -> dict[str, str]:
    header_index = {_normalize_header(h): i for i, h in enumerate(headers)}
    out: dict[str, str] = {}
    for label in REQUIRED_HEADERS:
        key = _normalize_header(label)
        idx = header_index.get(key)
        if idx is None or idx >= len(cells):
            out[key] = ""
        else:
            out[key] = str(cells[idx] or "").strip()
    return out


def _human_reference_row_detected(
    parsed_rows: Sequence[dict[str, Any]],
    *,
    resolved_entrez_gene_id: str,
) -> bool:
    human_gid = str(resolved_entrez_gene_id or "").strip()
    return any(
        str(row.get("gene_id") or "").strip() == human_gid
        or str(row.get("tax_id") or "").strip() == "9606"
        for row in parsed_rows
    )


def _expected_displayed_gene_count(
    *,
    scoped_ortholog_gene_count: int,
    human_reference_row_detected: bool,
) -> int:
    expected = int(scoped_ortholog_gene_count)
    if human_reference_row_detected:
        expected += 1
    return expected


def _count_consistency_passed(
    *,
    retrieval_complete: bool,
    displayed_gene_count: int | None,
    expected_displayed_count: int | None,
) -> bool | None:
    """Return True/False when retrieval is complete; None when hard claim is skipped.

    ``retrieval_complete`` here means both ortholog pagination and taxonomy
    membership resolution are complete enough to support a hard count claim.
    """
    if not retrieval_complete:
        return None
    if displayed_gene_count is None or expected_displayed_count is None:
        return False
    return int(displayed_gene_count) == int(expected_displayed_count)


def _count_claim_ready(
    *,
    pagination_complete: bool,
    taxonomy_complete: bool,
) -> bool:
    """Hard displayed==scoped(+human) checks only when both retrievals finished."""
    return bool(pagination_complete) and bool(taxonomy_complete)


def _official_capture_gate(
    *,
    capture_api_success: bool,
    capture_metadata: dict[str, Any] | None,
    captured_bytes: bytes | None,
    resolved_entrez_gene_id: str,
    configured_scope_tax_id: int,
    retrieval_complete: bool,
) -> bool:
    """Accept live Datasets DOM captures; hard count only when retrieval is complete.

    Incomplete taxonomy must not discard a successful live capture. When
    ``retrieval_complete`` is False the count check is skipped (returns None).
    """
    if not capture_api_success or not capture_metadata or captured_bytes is None:
        return False
    if capture_metadata.get("capture_origin") != CAPTURE_ORIGIN_LIVE:
        return False
    if str(capture_metadata.get("entrez_gene_id") or "").strip() != str(
        resolved_entrez_gene_id
    ).strip():
        return False
    try:
        selected_scope = int(capture_metadata.get("selected_scope_tax_id"))
    except (TypeError, ValueError):
        return False
    if selected_scope != int(configured_scope_tax_id):
        return False
    visible = list(capture_metadata.get("visible_headers") or [])
    if not _headers_cover_required(visible):
        return False
    digest = str(capture_metadata.get("sha256") or "").strip()
    if not digest or digest != hashlib.sha256(captured_bytes).hexdigest():
        return False
    count_ok = _count_consistency_passed(
        retrieval_complete=retrieval_complete,
        displayed_gene_count=capture_metadata.get("displayed_gene_count"),
        expected_displayed_count=capture_metadata.get("expected_displayed_count"),
    )
    if count_ok is False:
        return False
    return True

def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    raw_artifact_id: str | None = None,
    api_run_id: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=make_source_id(
            "NCBI Datasets", gene_symbol, AssertionType.gene_identity, key
        ),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_HOMOLOGUES,
        subsection=SUBSECTION_1E,
        source_name="NCBI Datasets",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type=fact_type,
        species="human",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.B,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def _group_values(gene_groups: Any) -> set[str]:
    out: set[str] = set()
    if isinstance(gene_groups, dict):
        gene_groups = [gene_groups]
    if not isinstance(gene_groups, list):
        return out
    for item in gene_groups:
        if isinstance(item, dict):
            for key in ("id", "gene_id", "name", "group_id", "identifier"):
                raw = item.get(key)
                if raw is not None:
                    out.add(str(raw).strip())
            for value in item.values():
                if isinstance(value, (str, int)):
                    out.add(str(value).strip())
        elif item is not None:
            out.add(str(item).strip())
    return {x for x in out if x}


def classify_membership(
    summary: dict[str, Any],
    *,
    query_gene_id: str,
) -> tuple[str, str | None]:
    """Return (membership_verified, reject_reason)."""
    groups = summary.get("gene_groups")
    values = _group_values(groups)
    qid = str(query_gene_id).strip()
    if groups not in (None, [], {}):
        if qid in values or any(qid in value for value in values):
            return "explicit", None
        return "rejected", "gene_groups_mismatch"
    query_ids = [str(x).strip() for x in (summary.get("query_gene_ids") or [])]
    if qid in query_ids:
        return "endpoint_implicit", None
    return "unverified", "missing_gene_groups_and_query_provenance"


def select_species_names(
    records: Sequence[dict[str, Any]],
    *,
    limit: int = 10,
) -> list[str]:
    """Deterministic distinct species names for prose."""
    by_tax: dict[int, dict[str, Any]] = {}
    for rec in records:
        try:
            tax_id = int(rec.get("tax_id"))
        except (TypeError, ValueError):
            continue
        if tax_id == 9606:
            continue
        by_tax.setdefault(tax_id, rec)
    ordered: list[dict[str, Any]] = []
    for tax_id in MODEL_SPECIES_PRIORITY:
        if tax_id in by_tax:
            ordered.append(by_tax.pop(tax_id))
    rest = sorted(
        by_tax.values(),
        key=lambda rec: str(
            rec.get("common_name") or rec.get("scientific_name") or rec.get("taxname") or ""
        ).lower(),
    )
    ordered.extend(rest)
    names: list[str] = []
    for rec in ordered:
        name = str(
            rec.get("common_name") or rec.get("scientific_name") or rec.get("taxname") or ""
        ).strip()
        if name and name not in names:
            names.append(name)
        if len(names) >= limit:
            break
    return names


def build_complete_narrative(
    *,
    gene_symbol: str,
    scope_label: str,
    ortholog_gene_count: int,
    species_names: Sequence[str],
) -> str:
    if species_names:
        if len(species_names) == 1:
            species_phrase = species_names[0]
        elif len(species_names) == 2:
            species_phrase = f"{species_names[0]} and {species_names[1]}"
        else:
            species_phrase = (
                ", ".join(species_names[:-1]) + f", and {species_names[-1]}"
            )
        including = f", including records from {species_phrase}"
    else:
        including = ""
    return (
        f"NCBI Orthologs provides comparative information for human {gene_symbol} "
        f"across {scope_label}. The current query identified "
        f"{ortholog_gene_count} orthologous genes{including}. "
        f"The NCBI Orthologs view for human {gene_symbol} is shown below, and "
        f"additional orthology information is available in OrthoDB."
    )


def build_incomplete_narrative(*, gene_symbol: str, scope_label: str) -> str:
    return (
        f"NCBI Orthologs provides comparative information for human {gene_symbol} "
        f"across {scope_label}. Complete structured ortholog retrieval was "
        f"unavailable for this run, so an exact ortholog count is not reported. "
        f"The NCBI Orthologs view for human {gene_symbol} is linked below when "
        f"available, and additional orthology information is available in OrthoDB."
    )


def _capture_passes_quality(meta: dict[str, Any], content: bytes) -> tuple[bool, str | None]:
    width = int(meta.get("width") or 0)
    height = int(meta.get("height") or 0)
    byte_size = int(meta.get("byte_size") or len(content))
    if width < _MIN_CAPTURE_WIDTH or height < _MIN_CAPTURE_HEIGHT:
        return False, f"capture dimensions too small ({width}x{height})"
    if byte_size < _MIN_CAPTURE_BYTES:
        return False, f"capture byte size too small ({byte_size})"
    return True, None


class _CaptureFailure(Exception):
    """Structured capture failure with a stable audit class."""

    def __init__(self, failure_class: str, message: str) -> None:
        super().__init__(message)
        self.failure_class = failure_class
        self.message = message


def _classify_capture_exception(exc: BaseException) -> tuple[str, str]:
    if isinstance(exc, _CaptureFailure):
        return exc.failure_class, exc.message
    text = f"{type(exc).__name__}: {exc}"
    lower = text.lower()
    if "timeout" in lower or "playwrighttimeouterror" in lower:
        return "navigation_failure", text
    if any(
        token in lower
        for token in (
            "403",
            "access forbidden",
            "access denied",
            "blocked",
            "consent",
            "captcha",
        )
    ):
        return "blocked_or_consent_page", text
    if "architecture" in lower:
        return "architecture_not_rendered", text
    if "consistency" in lower or "mismatch" in lower:
        return "dom_consistency_rejection", text
    if "quality" in lower or "byte size" in lower or "dimensions too small" in lower:
        return "image_quality_rejection", text
    if "path" in lower or "resolve" in lower:
        return "persisted_image_path_resolution_failure", text
    if "selector" in lower or "table" in lower or "header" in lower:
        return "table_selector_failure", text
    if "navigat" in lower or "host" in lower:
        return "navigation_failure", text
    return "capture_failed", text


def _dom_consistency(
    *,
    row_meta: list[dict[str, Any]],
    scoped_records: Sequence[dict[str, Any]],
    human_entrez_gene_id: str,
) -> dict[str, Any]:
    """Prefer stable IDs; common-name drift is warning-only."""
    allowed_gene_ids = {
        str(rec.get("gene_id")).strip()
        for rec in scoped_records
        if rec.get("gene_id") is not None
    }
    allowed_tax_ids = {
        str(rec.get("tax_id")).strip()
        for rec in scoped_records
        if rec.get("tax_id") is not None
    }
    # Human query row is present on the live NCBI page and allowed in the crop,
    # even though it is excluded from ortholog counts / API fallback tables.
    human_gid = str(human_entrez_gene_id or "").strip()
    if human_gid:
        allowed_gene_ids.add(human_gid)
    allowed_tax_ids.add("9606")

    gene_ok = True
    tax_ok = True
    name_warnings: list[str] = []
    compared_gene = 0
    compared_tax = 0
    for row in row_meta:
        gid = str(row.get("gene_id") or "").strip()
        tid = str(row.get("tax_id") or "").strip()
        if gid:
            compared_gene += 1
            if gid not in allowed_gene_ids:
                gene_ok = False
        if tid:
            compared_tax += 1
            if tid not in allowed_tax_ids:
                tax_ok = False
        # Scientific / common-name-only differences are warnings, never hard rejects.
        if not gid and not tid:
            label = str(row.get("scientific_name") or row.get("common_name") or "").strip()
            if label:
                name_warnings.append(label)

    if compared_gene or compared_tax:
        status = "pass" if gene_ok and tax_ok else "mismatch"
    else:
        status = "pass_with_name_warnings" if name_warnings else "pass"
    return {
        "status": status,
        "compared_gene_ids": compared_gene,
        "compared_tax_ids": compared_tax,
        "gene_ok": gene_ok,
        "tax_ok": tax_ok,
        "name_warnings": name_warnings[:20],
        "allowed_human_entrez_gene_id": human_gid or None,
    }


_ORTHOLOG_TABLE_FIND_JS = """
() => {
  const wanted = ["scientific name", "symbol", "length (aa)", "architecture"];
  const norm = (t) => (t || "").replace(/\\s+/g, " ").trim().toLowerCase();
  const headerNodes = Array.from(
    document.querySelectorAll("th, [role='columnheader']")
  );
  let architecture = null;
  for (const node of headerNodes) {
    if (norm(node.textContent) === "architecture") {
      architecture = node;
      break;
    }
  }
  if (!architecture) return null;

  const isContainer = (el) => {
    if (!el || el.nodeType !== 1) return false;
    const role = (el.getAttribute("role") || "").toLowerCase();
    const tag = el.tagName.toLowerCase();
    return tag === "table" || role === "table" || role === "grid";
  };
  let container = architecture.parentElement;
  while (container && !isContainer(container)) {
    container = container.parentElement;
  }
  if (!container) return null;

  const labels = Array.from(
    container.querySelectorAll("th, [role='columnheader']")
  ).map((n) => norm(n.textContent));
  const hasAll = wanted.every((w) => labels.includes(w));
  if (!hasAll) return null;

  if (!container.dataset.gdOrthologCapture) {
    container.dataset.gdOrthologCapture = "1";
  }
  return {
    selector: '[data-gd-ortholog-capture="1"]',
    tag: container.tagName.toLowerCase(),
    role: container.getAttribute("role") || null,
    header_labels: labels.slice(0, 20),
  };
}
"""


_ORTHOLOG_READY_JS = """
(sel) => {
  const root = document.querySelector(sel);
  if (!root) return { ready: false, reason: "container_missing" };
  const text = (root.innerText || "").toLowerCase();
  if (text.includes("loading") && !(root.querySelectorAll("tr, [role='row']").length > 3)) {
    return { ready: false, reason: "loading_skeleton" };
  }
  const rows = Array.from(root.querySelectorAll("tbody tr, [role='row']")).filter((row) => {
    const role = (row.getAttribute("role") || "").toLowerCase();
    const tag = row.tagName.toLowerCase();
    if (tag === "tr" && row.parentElement && row.parentElement.tagName.toLowerCase() === "thead") {
      return false;
    }
    if (role === "row" && row.querySelector("[role='columnheader']")) return false;
    const t = (row.innerText || "").trim();
    return t.length > 0;
  });
  if (rows.length < 2) return { ready: false, reason: "insufficient_rows", row_count: rows.length };

  const headers = Array.from(root.querySelectorAll("th, [role='columnheader']"));
  let archIndex = -1;
  headers.forEach((h, i) => {
    if ((h.textContent || "").replace(/\\s+/g, " ").trim().toLowerCase() === "architecture") {
      archIndex = i;
    }
  });
  let archRendered = false;
  const sample = rows.slice(0, Math.min(rows.length, 12));
  for (const row of sample) {
    let cell = null;
    if (archIndex >= 0) {
      const cells = row.querySelectorAll("td, [role='gridcell']");
      cell = cells[archIndex] || null;
    }
    const scope = cell || row;
    const candidates = scope.querySelectorAll(
      "div, span, svg, img, canvas, [class*='architect'], [class*='domain']"
    );
    for (const node of candidates) {
      const r = node.getBoundingClientRect();
      if (r && r.width >= 8 && r.height >= 4) {
        archRendered = true;
        break;
      }
    }
    if (archRendered) break;
    if (cell) {
      const r = cell.getBoundingClientRect();
      if (r && r.width >= 40 && r.height >= 8 && (cell.innerText || "").trim().length < 8) {
        archRendered = true;
        break;
      }
    }
  }
  return {
    ready: true,
    row_count: rows.length,
    architecture_rendered: archRendered,
  };
}
"""


_PARSE_TABLE_JS = """
(sel) => {
  const root = document.querySelector(sel);
  if (!root) return null;
  const norm = (t) => (t || "").replace(/\\s+/g, " ").trim();
  const headers = Array.from(
    root.querySelectorAll("th, [role='columnheader']")
  ).map((n) => norm(n.textContent));
  const rows = Array.from(root.querySelectorAll("tbody tr, [role='row']")).filter((row) => {
    const role = (row.getAttribute("role") || "").toLowerCase();
    const tag = row.tagName.toLowerCase();
    if (tag === "tr" && row.parentElement && row.parentElement.tagName.toLowerCase() === "thead") {
      return false;
    }
    if (role === "row" && row.querySelector("[role='columnheader']")) return false;
    return (row.innerText || "").trim().length > 0;
  });
  const parsed = rows.map((row) => {
    const cells = Array.from(row.querySelectorAll("td, [role='gridcell']")).map((c) =>
      norm(c.innerText)
    );
    const hrefs = Array.from(row.querySelectorAll("a")).map((a) => ({
      href: a.href || "",
      text: norm(a.textContent),
      gene: a.getAttribute("data-gene-id") || (a.dataset && a.dataset.geneId) || "",
      tax: a.getAttribute("data-tax-id") || (a.dataset && a.dataset.taxId) || "",
    }));
    let gene_id = null;
    let tax_id = null;
    for (const item of hrefs) {
      if (!gene_id && item.gene) gene_id = String(item.gene).trim();
      if (!tax_id && item.tax) tax_id = String(item.tax).trim();
      const mg = (item.href || "").match(/\\/gene\\/(\\d+)/);
      if (mg && !gene_id) gene_id = mg[1];
      const mt = (item.href || "").match(/taxonomy(?:\\/|id=)(\\d+)/i);
      if (mt && !tax_id) tax_id = mt[1];
    }
    return { cells, gene_id, tax_id, hrefs };
  });
  return { headers, rows: parsed };
}
"""


_READ_DISPLAYED_COUNT_JS = """
() => {
  const text = document.body ? document.body.innerText || "" : "";
  const patterns = [
    /([\\d,]+)\\s+Genes?/i,
    /([\\d,]+)\\s+orthologs?/i,
    /showing\\s+([\\d,]+)/i,
    /1-\\d+\\s+of\\s+([\\d,]+)/i,
  ];
  for (const re of patterns) {
    const m = text.match(re);
    if (m) {
      const n = parseInt(m[1].replace(/,/g, ""), 10);
      if (!Number.isNaN(n) && n > 0) return n;
    }
  }
  const aria = document.querySelector("[aria-label*='ortholog' i], [data-testid*='ortholog' i]");
  if (aria) {
    const m = (aria.textContent || "").match(/([\\d,]+)/);
    if (m) {
      const n = parseInt(m[1].replace(/,/g, ""), 10);
      if (!Number.isNaN(n) && n > 0) return n;
    }
  }
  return null;
}
"""


def _datasets_gene_route_ok(url: str, *, entrez_gene_id: str) -> bool:
    parsed = urlparse(url)
    host = (parsed.hostname or "").lower()
    if host not in _NCBI_ALLOWED_HOSTS:
        return False
    path = (parsed.path or "").rstrip("/")
    return path == f"/datasets/gene/{str(entrez_gene_id).strip()}"


def _orthologs_tab_selected(page: Any) -> bool:
    """True when Orthologs UI is active (aria/content), not merely URL fragment."""
    tab = page.locator(
        "[role='tab']:has-text('Orthologs'), button:has-text('Orthologs'), "
        "a:has-text('Orthologs')"
    ).first
    if tab.count() == 0:
        return False
    aria = (tab.get_attribute("aria-selected") or "").strip().lower()
    if aria == "true":
        return True
    # Datasets may use class/state without aria-selected; require Orthologs UI chrome.
    taxa = page.locator(
        "text=/Selected taxa/i, [placeholder*='taxa' i], "
        "input[aria-label*='taxa' i], input[placeholder*='taxonomy' i]"
    )
    table_hint = page.locator("th:has-text('Architecture'), [role='columnheader']:has-text('Architecture')")
    heading = page.locator("text=/Orthologs/i").first
    return taxa.count() > 0 and table_hint.count() > 0 and heading.count() > 0


def _activate_orthologs_tab(page: Any) -> None:
    if _orthologs_tab_selected(page):
        return
    tab = page.locator(
        "[role='tab']:has-text('Orthologs'), button:has-text('Orthologs'), "
        "a[href*='ortholog' i]:has-text('Orthologs')"
    ).first
    if tab.count() == 0:
        raise _CaptureFailure(
            "table_selector_failure",
            "Orthologs tab not found on NCBI Datasets gene page",
        )
    tab.click(timeout=15_000)
    for _ in range(40):
        if _orthologs_tab_selected(page):
            return
        page.wait_for_timeout(250)
    raise _CaptureFailure(
        "table_selector_failure",
        "Orthologs tab did not become selected / Orthologs UI not visible",
    )


def _read_displayed_gene_count(page: Any) -> int | None:
    raw = page.evaluate(_READ_DISPLAYED_COUNT_JS)
    if raw is None:
        return None
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _first_page_row_ids(page: Any, selector: str) -> list[str]:
    parsed = page.evaluate(_PARSE_TABLE_JS, selector) or {}
    ids: list[str] = []
    for row in parsed.get("rows") or []:
        gid = str(row.get("gene_id") or "").strip()
        tid = str(row.get("tax_id") or "").strip()
        ids.append(f"{gid}|{tid}")
    return ids


def _wait_loading_cycle(page: Any) -> None:
    loading = page.locator(
        "text=/loading/i, [aria-busy='true'], .MuiCircularProgress-root, "
        "[class*='skeleton' i], [class*='Spinner' i]"
    )
    try:
        loading.first.wait_for(state="visible", timeout=2_000)
    except Exception:
        pass
    try:
        loading.first.wait_for(state="hidden", timeout=30_000)
    except Exception:
        pass


def _stabilize_filtered_count(
    page: Any,
    *,
    table_selector: str | None = None,
) -> dict[str, Any]:
    """Wait until displayed count + first-page row IDs are stable across two polls."""
    old_count = _read_displayed_gene_count(page)
    _wait_loading_cycle(page)
    for _ in range(60):
        count = _read_displayed_gene_count(page)
        if count is not None and count > 0:
            break
        page.wait_for_timeout(250)
    else:
        count = _read_displayed_gene_count(page)

    row_ids: list[str] = []
    if table_selector:
        for _ in range(40):
            row_ids = _first_page_row_ids(page, table_selector)
            if row_ids:
                break
            page.wait_for_timeout(250)

    stable_count = None
    stable_ids: list[str] = []
    for _ in range(20):
        page.wait_for_timeout(400)
        c2 = _read_displayed_gene_count(page)
        ids2 = _first_page_row_ids(page, table_selector) if table_selector else []
        if (
            c2 is not None
            and c2 == count
            and (not table_selector or (ids2 and ids2 == row_ids))
        ):
            stable_count = c2
            stable_ids = ids2
            break
        count = c2
        row_ids = ids2
    return {
        "old_count": old_count,
        "displayed_gene_count": stable_count if stable_count is not None else count,
        "first_page_row_ids": stable_ids or row_ids,
        "stabilized": stable_count is not None,
    }


def _dismiss_blocking_overlays(page: Any) -> None:
    """Remove common NCBI/Qualtrics overlays that intercept clicks."""
    try:
        page.evaluate(
            """() => {
              for (const sel of [
                '.QSIWebResponsive',
                '.QSIWebResponsive-creative-container-fade',
                '#onetrust-banner-sdk',
                '.onetrust-pc-dark-filter',
              ]) {
                document.querySelectorAll(sel).forEach((el) => el.remove());
              }
            }"""
        )
    except Exception:
        pass
    for sel in (
        "button:has-text('No thanks')",
        "button:has-text('Close')",
        "#onetrust-accept-btn-handler",
        "[aria-label='Close']",
    ):
        btn = page.locator(sel).first
        if btn.count() == 0:
            continue
        try:
            btn.click(timeout=1_000)
        except Exception:
            continue


def _select_scope_taxa(
    page: Any,
    *,
    taxon_name: str,
    scope_tax_id: int,
) -> dict[str, Any]:
    """Type ncbi_taxon_name into Selected taxa and pick the exact autocomplete result."""
    info: dict[str, Any] = {"taxon_name": taxon_name, "scope_tax_id": scope_tax_id}
    # Clear existing chips when a clear/remove control is available.
    for _ in range(8):
        clear_btn = page.locator(
            "button[aria-label*='clear' i], button[aria-label*='remove' i], "
            "[data-testid*='clear' i], button[aria-label*='Delete' i]"
        ).first
        if clear_btn.count() == 0:
            break
        try:
            clear_btn.click(timeout=1_000)
            page.wait_for_timeout(200)
        except Exception:
            break

    input_box = page.locator(
        "input[placeholder*='taxonomic names' i], "
        "input[placeholder*='taxa' i], input[aria-label*='taxa' i], "
        "input[placeholder*='taxonomy' i], input[aria-label*='Selected taxa' i], "
        "#taxonomy_autocomplete, "
        "[role='combobox'] input"
    ).first
    if input_box.count() == 0:
        label = page.locator("text=/Selected taxa/i").first
        if label.count() == 0:
            raise _CaptureFailure(
                "table_selector_failure",
                "Selected taxa control not found",
            )
        label.click(timeout=5_000)
        input_box = page.locator(
            "input[placeholder*='taxonomic names' i], input[placeholder*='taxa' i], "
            "#taxonomy_autocomplete"
        ).first
        if input_box.count() == 0:
            raise _CaptureFailure(
                "table_selector_failure",
                "Selected taxa input not found after label click",
            )

    _dismiss_blocking_overlays(page)
    input_box.click(timeout=10_000, force=True)
    input_box.fill("")
    input_box.type(taxon_name, delay=40)
    page.wait_for_timeout(700)

    matched = False
    chosen_index = None
    for _ in range(30):
        options = page.locator("[role='option']")
        count = options.count()
        if count == 0:
            page.wait_for_timeout(250)
            continue
        # Prefer exact TaxID match when exposed in the option label.
        for i in range(min(count, 20)):
            text = (options.nth(i).inner_text() or "").strip()
            if f"TaxID: {scope_tax_id}" in text or f"TaxID:{scope_tax_id}" in text:
                chosen_index = i
                info["selected_option_text"] = text
                break
        if chosen_index is None:
            for i in range(min(count, 20)):
                text = (options.nth(i).inner_text() or "").strip()
                lower = text.lower()
                if lower.startswith(taxon_name.lower() + " ") or lower == taxon_name.lower():
                    chosen_index = i
                    info["selected_option_text"] = text
                    break
        if chosen_index is not None:
            _dismiss_blocking_overlays(page)
            options.nth(chosen_index).click(timeout=5_000, force=True)
            matched = True
            break
        page.wait_for_timeout(250)
    if not matched:
        raise _CaptureFailure(
            "table_selector_failure",
            f"Autocomplete option for taxon '{taxon_name}' not found",
        )

    page.wait_for_timeout(800)
    body = page.locator("body").inner_text(timeout=5_000)
    info["tax_id_visible"] = str(scope_tax_id) in body
    return info


def _configure_visible_columns(page: Any) -> list[str]:
    """Prefer official Select columns UI; keep required headers, hide others where possible."""
    keep = {_normalize_header(h) for h in REQUIRED_HEADERS}
    hide_candidates = (
        "Chromosome",
        "Isoforms",
        "Protein accession",
        "Genomic coordinates",
        "Orientation",
        "Action",
        "Common name",
    )
    _dismiss_blocking_overlays(page)
    btn = page.locator(
        "button:has-text('Select columns'), button:has-text('Columns'), "
        "[aria-label*='Select columns' i], [aria-label*='columns' i]"
    ).first
    if btn.count() > 0:
        btn.click(timeout=10_000, force=True)
        page.wait_for_timeout(400)
        menu = page.locator(
            "[role='menu'], [role='dialog'], [role='listbox'], .MuiPopover-root, "
            ".MuiMenu-paper"
        ).last
        for label in hide_candidates:
            item = menu.locator(
                f"label:has-text('{label}'), [role='menuitem']:has-text('{label}'), "
                f"span:has-text('{label}')"
            ).first
            if item.count() == 0:
                continue
            try:
                checkbox = item.locator("input[type='checkbox']").first
                if checkbox.count() == 0:
                    # Click the row/label; many MUI menus wrap the checkbox.
                    parent = item.locator("xpath=ancestor-or-self::*[contains(@class,'MuiFormControlLabel') or @role='menuitem'][1]")
                    target = parent if parent.count() else item
                    # Only toggle off if currently checked.
                    inp = target.locator("input[type='checkbox']").first
                    if inp.count() and inp.is_checked():
                        target.click(timeout=2_000, force=True)
                elif checkbox.is_checked():
                    checkbox.uncheck(timeout=2_000, force=True)
            except Exception:
                continue
        for label in REQUIRED_HEADERS:
            item = menu.locator(
                f"label:has-text('{label}'), [role='menuitem']:has-text('{label}')"
            ).first
            if item.count() == 0:
                continue
            try:
                checkbox = item.locator("input[type='checkbox']").first
                if checkbox.count() > 0 and not checkbox.is_checked():
                    checkbox.check(timeout=2_000, force=True)
            except Exception:
                continue
        page.keyboard.press("Escape")
        page.wait_for_timeout(400)

    headers = page.evaluate(
        """() => {
          const nodes = Array.from(document.querySelectorAll("th, [role='columnheader']"));
          return nodes.map((n) => (n.textContent || "").replace(/\\s+/g, " ").trim()).filter(Boolean);
        }"""
    )
    visible = [str(h) for h in (headers or []) if str(h).strip()]
    if not _headers_cover_required(visible):
        page.wait_for_timeout(800)
        headers = page.evaluate(
            """() => {
              const nodes = Array.from(document.querySelectorAll("th, [role='columnheader']"));
              return nodes.map((n) => (n.textContent || "").replace(/\\s+/g, " ").trim()).filter(Boolean);
            }"""
        )
        visible = [str(h) for h in (headers or []) if str(h).strip()]
    if not _headers_cover_required(visible):
        raise _CaptureFailure(
            "table_selector_failure",
            f"Required headers not visible; saw {visible[:20]}",
        )
    _ = keep
    return visible


def _set_rows_per_page(page: Any, rows: int = 20) -> None:
    _dismiss_blocking_overlays(page)
    control = page.get_by_text("Rows per page", exact=False).first
    if control.count() == 0:
        control = page.locator("[aria-label*='Rows per page' i]").first
    if control.count() == 0:
        return
    try:
        select = page.locator(
            "[aria-label*='Rows per page' i]"
        ).first
        if select.count() == 0:
            # MUI outlined select near the label.
            select = page.locator(
                "div:near(:text('Rows per page')) >> [role='combobox']"
            ).first
        if select.count() == 0:
            control.click(timeout=3_000, force=True)
            opt = page.locator(f"[role='option']:has-text('{rows}')").first
            if opt.count() > 0:
                opt.click(timeout=3_000, force=True)
            return
        tag = select.evaluate("el => el.tagName.toLowerCase()")
        if tag == "select":
            select.select_option(str(rows))
        else:
            select.click(timeout=3_000, force=True)
            opt = page.locator(f"[role='option']:has-text('{rows}')").first
            if opt.count() > 0:
                opt.click(timeout=3_000, force=True)
    except Exception:
        logger.debug("rows-per-page control not adjusted", exc_info=True)


def _extract_row_stable_ids_from_parsed(
    headers: Sequence[str],
    row: dict[str, Any],
) -> dict[str, Any]:
    cells = [str(c) for c in (row.get("cells") or [])]
    by_header = _parse_row_cells_by_header(headers, cells)
    return {
        "gene_id": str(row.get("gene_id") or "").strip() or None,
        "tax_id": str(row.get("tax_id") or "").strip() or None,
        "scientific_name": by_header.get("scientific name") or "",
        "symbol": by_header.get("symbol") or "",
        "length_aa": by_header.get("length (aa)") or "",
        "architecture": by_header.get("architecture") or "",
        "common_name": "",
        "text": " | ".join(cells)[:300],
        "cells_by_header": by_header,
    }


def _capture_ncbi_ortholog_table(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    entrez_gene_id: str,
    scope_tax_id: int,
    page_url: str,
    max_visible_rows: int,
    scoped_records: Sequence[dict[str, Any]],
    scoped_ortholog_gene_count: int | None,
    retrieval_complete: bool,
    parent_raw_artifact_ids: Sequence[str],
    settings: Settings,
    persist_db: bool,
    taxon_name: str | None = None,
    count_claim_ready: bool | None = None,
) -> tuple[ApiRun, dict[str, Any] | None, EvidenceRecord | None, dict[str, Any], bytes | None]:
    scope_meta = SUPPORTED_SECTION_1E_SCOPES.get(int(scope_tax_id), {})
    taxon = (taxon_name or str(scope_meta.get("ncbi_taxon_name") or "")).strip()
    # Hard displayed==scoped checks require complete pagination+taxonomy.
    # Prefer the explicit flag when callers distinguish taxonomy from pagination.
    hard_count_ready = (
        bool(retrieval_complete)
        if count_claim_ready is None
        else bool(count_claim_ready)
    )
    api = ApiRun(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="NCBI Gene",
        endpoint_name="capture_ncbi_ortholog_table",
        request_url=page_url,
        request_params={
            "entrez_gene_id": entrez_gene_id,
            "scope_tax_id": scope_tax_id,
            "max_visible_rows": max_visible_rows,
            "retrieval_method": "official_web_element_capture",
            "capture_target": "ncbi_datasets_gene_orthologs",
            "count_claim_ready": hard_count_ready,
        },
        success=False,
    )
    audit: dict[str, Any] = {
        "status": "unavailable",
        "source_page_url": page_url,
        "scope_tax_id": scope_tax_id,
        "attempts": [],
        "failure_class": None,
        "count_claim_ready": hard_count_ready,
    }
    try:
        from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # noqa: BLE001
        api.error_type = "playwright_unavailable"
        api.error_message = f"{type(exc).__name__}: {exc}"
        _save_api_run_failure(api, persist_db=persist_db)
        audit["reason"] = api.error_message
        audit["failure_class"] = "navigation_failure"
        return api, None, None, audit, None

    last_error: str | None = None
    last_failure_class = "capture_failed"
    for attempt in range(1, 3):
        attempt_info: dict[str, Any] = {"attempt": attempt}
        browser = None
        try:
            with sync_playwright() as pw:
                try:
                    browser = pw.chromium.launch(headless=True, channel="chrome")
                except Exception:
                    browser = pw.chromium.launch(headless=True)
                page = browser.new_page(
                    viewport={"width": 1400, "height": 1800},
                    user_agent=(
                        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                        "AppleWebKit/537.36 (KHTML, like Gecko) "
                        "Chrome/122.0.0.0 Safari/537.36"
                    ),
                )
                page.set_extra_http_headers(
                    {
                        "Accept-Language": "en-US,en;q=0.9",
                        "Upgrade-Insecure-Requests": "1",
                    }
                )
                try:
                    page.goto(page_url, wait_until="domcontentloaded", timeout=60_000)
                except PlaywrightTimeoutError as exc:
                    raise _CaptureFailure(
                        "navigation_failure", f"navigation timeout: {exc}"
                    ) from exc
                # Datasets gene pages are SPAs; allow hydrate before identity checks.
                page.wait_for_timeout(3_000)
                _dismiss_blocking_overlays(page)
                host = urlparse(page.url).hostname or ""
                if host not in _NCBI_ALLOWED_HOSTS:
                    raise _CaptureFailure(
                        "navigation_failure",
                        f"Unexpected NCBI host after navigation: {page.url}",
                    )
                body = page.locator("body").inner_text(timeout=15_000)
                lower = body.lower()
                title = (page.title() or "").lower()
                if (
                    "403" in title
                    or "access forbidden" in lower
                    or "do not have access to this page" in lower
                    or "access denied" in lower
                    or ("enable cookies" in lower and "captcha" in lower)
                    or not body.strip()
                ):
                    raise _CaptureFailure(
                        "blocked_or_consent_page",
                        f"NCBI Datasets page blocked ({page.title() or 'unknown title'})",
                    )
                if gene_symbol.lower() not in lower:
                    raise _CaptureFailure(
                        "blocked_or_consent_page",
                        "Expected gene symbol not found on page",
                    )
                # Datasets SPA often keeps Entrez only in the route; require either
                # body text or the Datasets gene path for the resolved ID.
                if (
                    entrez_gene_id not in body
                    and not _datasets_gene_route_ok(page.url, entrez_gene_id=entrez_gene_id)
                ):
                    # Brief wait for SPA hydrate then re-check.
                    page.wait_for_timeout(2_000)
                    body = page.locator("body").inner_text(timeout=15_000)
                    if (
                        entrez_gene_id not in body
                        and not _datasets_gene_route_ok(
                            page.url, entrez_gene_id=entrez_gene_id
                        )
                    ):
                        raise _CaptureFailure(
                            "blocked_or_consent_page",
                            "Expected Entrez Gene ID not found on page or URL",
                        )
                if not _datasets_gene_route_ok(page.url, entrez_gene_id=entrez_gene_id):
                    # Allow brief SPA settle then re-check route only (no #orthologs req).
                    page.wait_for_timeout(1_000)
                    if not _datasets_gene_route_ok(page.url, entrez_gene_id=entrez_gene_id):
                        raise _CaptureFailure(
                            "navigation_failure",
                            f"Not on Datasets gene route for {entrez_gene_id}: {page.url}",
                        )

                _activate_orthologs_tab(page)
                _dismiss_blocking_overlays(page)
                if not _datasets_gene_route_ok(page.url, entrez_gene_id=entrez_gene_id):
                    raise _CaptureFailure(
                        "navigation_failure",
                        f"Left Datasets gene route after Orthologs tab: {page.url}",
                    )

                taxa_info = _select_scope_taxa(
                    page, taxon_name=taxon, scope_tax_id=scope_tax_id
                )
                attempt_info["taxa"] = taxa_info

                # Locate table container (may appear after taxa filter).
                found = None
                for _ in range(40):
                    found = page.evaluate(_ORTHOLOG_TABLE_FIND_JS)
                    if found and found.get("selector"):
                        break
                    page.wait_for_timeout(500)
                # Stabilize count even if table finder still settling.
                stabilize = _stabilize_filtered_count(
                    page,
                    table_selector=(found or {}).get("selector"),
                )
                attempt_info["stabilize"] = stabilize
                if not found or not found.get("selector"):
                    for _ in range(20):
                        found = page.evaluate(_ORTHOLOG_TABLE_FIND_JS)
                        if found and found.get("selector"):
                            break
                        page.wait_for_timeout(500)
                if not found or not found.get("selector"):
                    raise _CaptureFailure(
                        "table_selector_failure",
                        "Scientific name/Symbol/Length (aa)/Architecture table not found",
                    )
                selector = str(found["selector"])
                attempt_info["container_selector"] = selector
                attempt_info["container_meta"] = found

                # Re-stabilize with table selector for row-id stability.
                stabilize = _stabilize_filtered_count(page, table_selector=selector)
                attempt_info["stabilize"] = stabilize
                displayed_gene_count = stabilize.get("displayed_gene_count")

                visible_headers = _configure_visible_columns(page)
                attempt_info["visible_headers"] = visible_headers
                _set_rows_per_page(page, 20)

                # Re-find table after column changes.
                found = page.evaluate(_ORTHOLOG_TABLE_FIND_JS)
                if not found or not found.get("selector"):
                    raise _CaptureFailure(
                        "table_selector_failure",
                        "Ortholog table missing after column configuration",
                    )
                selector = str(found["selector"])

                ready = None
                for _ in range(40):
                    ready = page.evaluate(_ORTHOLOG_READY_JS, selector)
                    if ready and ready.get("ready"):
                        break
                    page.wait_for_timeout(500)
                if not ready or not ready.get("ready"):
                    raise _CaptureFailure(
                        "table_selector_failure",
                        f"Ortholog table not ready: {(ready or {}).get('reason')}",
                    )
                rendered_row_count = int(ready.get("row_count") or 0)
                if rendered_row_count < 2:
                    raise _CaptureFailure(
                        "table_selector_failure",
                        "Ortholog table rows not rendered",
                    )
                if not ready.get("architecture_rendered"):
                    raise _CaptureFailure(
                        "architecture_not_rendered",
                        "Architecture graphics not rendered with nonzero dimensions",
                    )

                parsed_table = page.evaluate(_PARSE_TABLE_JS, selector) or {}
                headers = [str(h) for h in (parsed_table.get("headers") or [])]
                if not _headers_cover_required(headers):
                    raise _CaptureFailure(
                        "table_selector_failure",
                        f"Required headers missing after parse; saw {headers[:20]}",
                    )
                visible_headers = headers
                raw_rows = list(parsed_table.get("rows") or [])
                visible_n = min(max(1, max_visible_rows), max(0, len(raw_rows)))
                if visible_n < 1:
                    raise _CaptureFailure(
                        "table_selector_failure", "No data rows available to crop"
                    )
                row_meta = [
                    _extract_row_stable_ids_from_parsed(headers, raw_rows[index])
                    for index in range(visible_n)
                ]

                human_detected = _human_reference_row_detected(
                    row_meta, resolved_entrez_gene_id=entrez_gene_id
                )
                expected_displayed: int | None = None
                if scoped_ortholog_gene_count is not None:
                    expected_displayed = _expected_displayed_gene_count(
                        scoped_ortholog_gene_count=int(scoped_ortholog_gene_count),
                        human_reference_row_detected=human_detected,
                    )
                count_ok = _count_consistency_passed(
                    retrieval_complete=hard_count_ready,
                    displayed_gene_count=displayed_gene_count
                    if isinstance(displayed_gene_count, int)
                    else None,
                    expected_displayed_count=expected_displayed,
                )
                if count_ok is False:
                    raise _CaptureFailure(
                        "dom_consistency_rejection",
                        (
                            f"displayed_gene_count={displayed_gene_count} != "
                            f"expected_displayed_count={expected_displayed} "
                            f"(scoped={scoped_ortholog_gene_count}, "
                            f"human_row={human_detected})"
                        ),
                    )

                consistency = {"status": "skipped"}
                if hard_count_ready:
                    consistency = _dom_consistency(
                        row_meta=row_meta,
                        scoped_records=scoped_records,
                        human_entrez_gene_id=entrez_gene_id,
                    )
                    if consistency.get("status") == "mismatch":
                        raise _CaptureFailure(
                            "dom_consistency_rejection",
                            "capture DOM consistency mismatch vs scoped records",
                        )

                table = page.locator(selector).first
                table.wait_for(state="visible", timeout=10_000)
                header = table.locator(
                    "thead tr, tr:has(th), [role='row']:has([role='columnheader'])"
                ).first
                if header.count() == 0:
                    header = table.locator("tr, [role='row']").first
                header_box = header.bounding_box()

                skip_header = False
                if table.locator("tbody tr").count() > 0:
                    data_rows = table.locator("tbody tr")
                elif table.locator(
                    "[role='row']:not(:has([role='columnheader']))"
                ).count() > 0:
                    data_rows = table.locator(
                        "[role='row']:not(:has([role='columnheader']))"
                    )
                else:
                    data_rows = table.locator("tr")
                    skip_header = True

                data_row_count = data_rows.count()
                if skip_header and data_row_count > 0:
                    data_row_count = max(0, data_row_count - 1)
                visible_n = min(max(1, max_visible_rows), max(0, data_row_count), visible_n)
                start = 1 if skip_header else 0
                target_index = start + visible_n - 1
                bottom_box = data_rows.nth(target_index).bounding_box()
                table_box = table.bounding_box()
                if not header_box or not bottom_box or not table_box:
                    raise _CaptureFailure(
                        "table_selector_failure", "Unable to compute table crop box"
                    )
                clip = {
                    "x": max(0, table_box["x"]),
                    "y": max(0, header_box["y"]),
                    "width": max(1, table_box["width"]),
                    "height": max(
                        1, (bottom_box["y"] + bottom_box["height"]) - header_box["y"]
                    ),
                }

                png_bytes = page.screenshot(type="png", clip=clip)
                digest = hashlib.sha256(png_bytes).hexdigest()
                validation = _validate_nonblank_image(png_bytes)
                ok, reason = _capture_passes_quality(validation, png_bytes)
                capture_metadata = {
                    "capture_origin": CAPTURE_ORIGIN_LIVE,
                    "sha256": digest,
                    "entrez_gene_id": entrez_gene_id,
                    "selected_scope_tax_id": int(scope_tax_id),
                    "visible_headers": list(visible_headers),
                    "displayed_gene_count": displayed_gene_count,
                    "expected_displayed_count": expected_displayed,
                    "human_reference_row_detected": human_detected,
                    "count_consistency_passed": count_ok,
                    "retrieval_complete": retrieval_complete,
                    "count_claim_ready": hard_count_ready,
                }
                attempt_info.update(
                    {
                        "final_url": page.url,
                        "clip": clip,
                        "rendered_row_count": rendered_row_count,
                        "captured_row_count": visible_n,
                        "visible_row_count": visible_n,
                        "width": validation.get("width"),
                        "height": validation.get("height"),
                        "byte_size": validation.get("byte_size", len(png_bytes)),
                        "consistency": consistency,
                        "stable_ids": row_meta,
                        **capture_metadata,
                    }
                )
                if not ok:
                    raise _CaptureFailure(
                        "image_quality_rejection",
                        reason or "capture quality gate failed",
                    )

                artifact, meta = _persist_artifact_bytes(
                    dossier_run_id=dossier_run_id,
                    source_name="NCBI Gene",
                    content=png_bytes,
                    extension="png",
                    artifact_type="png",
                    filename_hint=f"ncbi-orthologs-{slugify(gene_symbol)}",
                    settings=settings,
                    api_run=api,
                    persist_db=persist_db,
                    notes={
                        "artifact_class": "derived_capture",
                        "artifact_origin": "ncbi_orthologs_table",
                        "artifact_role": "section_1e_ortholog_capture",
                        "retrieval_method": "official_web_element_capture",
                        "capture_origin": CAPTURE_ORIGIN_LIVE,
                        "source_page_url": page_url,
                        "final_url": page.url,
                        "scope_tax_id": scope_tax_id,
                        "selected_scope_tax_id": int(scope_tax_id),
                        "entrez_gene_id": entrez_gene_id,
                        "container_selector": selector,
                        "rendered_row_count": rendered_row_count,
                        "captured_row_count": visible_n,
                        "visible_row_count": visible_n,
                        "visible_headers": list(visible_headers),
                        "displayed_gene_count": displayed_gene_count,
                        "expected_displayed_count": expected_displayed,
                        "human_reference_row_detected": human_detected,
                        "count_consistency_passed": count_ok,
                        "clip": clip,
                        "stable_ids": row_meta,
                        "sha256": digest,
                        "parent_raw_artifact_ids": list(parent_raw_artifact_ids),
                        "captured_at": datetime.now(timezone.utc).isoformat(),
                        "consistency": consistency,
                    },
                    validate=_validate_nonblank_image,
                )
                api.success = True
                if persist_db:
                    from gene_dossier.db import save_api_run, session_scope

                    with session_scope() as session:
                        save_api_run(session, api)
                value = {
                    "status": "success",
                    "source_page_url": page_url,
                    "final_url": page.url,
                    "scope_tax_id": scope_tax_id,
                    "entrez_gene_id": entrez_gene_id,
                    "selected_scope_tax_id": int(scope_tax_id),
                    "relative_path": meta.get("relative_path"),
                    "media_type": meta.get("media_type") or "image/png",
                    "width": meta.get("width"),
                    "height": meta.get("height"),
                    "sha256": digest,
                    "byte_size": meta.get("byte_size"),
                    "artifact_class": "derived_capture",
                    "capture_origin": CAPTURE_ORIGIN_LIVE,
                    "container_selector": selector,
                    "rendered_row_count": rendered_row_count,
                    "captured_row_count": visible_n,
                    "visible_row_count": visible_n,
                    "visible_headers": list(visible_headers),
                    "displayed_gene_count": displayed_gene_count,
                    "expected_displayed_count": expected_displayed,
                    "human_reference_row_detected": human_detected,
                    "count_consistency_passed": count_ok,
                    "clip": clip,
                    "stable_ids": row_meta,
                    "consistency": consistency,
                    "figure_raw_artifact_id": artifact.id,
                    "presentation_item_key": f"orthologs-{slugify(gene_symbol)}",
                }
                rec = _record(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    fact_type="ortholog_table_capture",
                    key=f"capture-{entrez_gene_id}",
                    value=value,
                    display_text=(
                        f"{gene_symbol} official NCBI Orthologs table capture "
                        f"(scope {scope_tax_id})."
                    ),
                    raw_artifact_id=artifact.id,
                    api_run_id=api.id,
                )
                audit = {
                    "status": "success",
                    "failure_class": None,
                    "source_page_url": page_url,
                    "final_url": page.url,
                    "scope_tax_id": scope_tax_id,
                    "attempts": audit.get("attempts", []) + [attempt_info],
                    "container_selector": selector,
                    "rendered_row_count": rendered_row_count,
                    "captured_row_count": visible_n,
                    "visible_row_count": visible_n,
                    "clip": clip,
                    "stable_ids": row_meta,
                    "capture_origin": CAPTURE_ORIGIN_LIVE,
                    "visible_headers": list(visible_headers),
                    "displayed_gene_count": displayed_gene_count,
                    "expected_displayed_count": expected_displayed,
                    "human_reference_row_detected": human_detected,
                    "count_consistency_passed": count_ok,
                    "quality": {
                        "width": meta.get("width"),
                        "height": meta.get("height"),
                        "byte_size": meta.get("byte_size"),
                        "sha256": digest,
                    },
                    "consistency": consistency,
                    "capture_metadata": capture_metadata,
                }
                return api, meta, rec, audit, png_bytes
        except Exception as exc:  # noqa: BLE001
            failure_class, message = _classify_capture_exception(exc)
            last_failure_class = failure_class
            last_error = message
            attempt_info["error"] = message
            attempt_info["failure_class"] = failure_class
        finally:
            audit.setdefault("attempts", []).append(attempt_info)
            if browser is not None:
                try:
                    browser.close()
                except Exception:
                    pass

    api.error_type = last_failure_class
    api.error_message = last_error or "NCBI Orthologs table capture failed"
    _save_api_run_failure(api, persist_db=persist_db)
    audit["status"] = "unavailable"
    audit["reason"] = api.error_message
    audit["failure_class"] = last_failure_class
    return api, None, None, audit, None

def node_generate_section_1e_derived_artifacts(
    state: DossierState,
    *,
    settings: Settings | None = None,
    persist_db: bool = True,
    transient: WorkflowTransientContext | None = None,
    config: Section1eConfig | None = None,
    skip_table_capture: bool = False,
) -> DossierState:
    """Sole Section 1e network owner: Datasets pages, taxonomy, OrthoDB, capture."""
    if state.get("run_type") != "section_bundle" or "1e" not in (
        state.get("selected_section_keys") or []
    ):
        return state

    cfg = settings or get_settings()
    section_cfg = config or Section1eConfig()
    run_id = state["dossier_run_id"]
    gene = state["gene_symbol"]
    gene_ids = dict(state.get("gene_ids") or {})
    entrez = str(gene_ids.get("entrez_gene_id") or "").strip()
    evidence = list(state.get("evidence_records") or [])
    api_runs = list(state.get("api_runs") or [])
    raw_meta = list(state.get("raw_artifacts") or [])
    errors = list(state.get("errors") or [])
    coverage_extra = list(state.get("coverage") or [])

    audit: dict[str, Any] = {
        "config": {
            "ortholog_scope_tax_id": section_cfg.ortholog_scope_tax_id,
            "ortholog_scope_label": section_cfg.ortholog_scope_label,
            "max_visible_rows": section_cfg.max_visible_rows,
        },
        "pagination": {},
        "membership": {},
        "taxonomy": {},
        "orthodb": {},
        "viewer_capture": {},
        "table_status": TABLE_UNAVAILABLE,
        "network_owner": "section_1e",
    }

    if not entrez:
        audit["reason"] = "human entrez_gene_id unavailable"
        coverage_extra.append(
            SourceCoverageResult(
                dossier_run_id=run_id,
                source_name="NCBI Datasets",
                status=SourceStatus.missing,
                evidence_record_count=0,
                error_message="Human Entrez Gene ID unavailable for Section 1e",
                report_sections_supported=["Homologues in model animals"],
            )
        )
        return {
            **state,
            "evidence_records": evidence,
            "api_runs": api_runs,
            "raw_artifacts": raw_meta,
            "errors": errors,
            "coverage": coverage_extra,
            "section_1e_status": {
                "rendering_status": {"overall": "empty", "count_status": "failed"},
                "audit": audit,
            },
        }

    pages, pagination = ncbi_datasets.iter_ortholog_pages(
        entrez, gene_symbol=gene, settings=cfg
    )
    audit["pagination"] = pagination
    parent_raw_ids: list[str] = []
    summaries: list[dict[str, Any]] = []
    for page in pages:
        api, meta = _persist_tool_result_json(
            tr=page,
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            filename_hint=f"ncbi-datasets-orthologs-{slugify(entrez)}",
        )
        api_runs.append(api)
        if meta:
            raw_meta.append(meta)
            rid = str(meta.get("id") or "").strip()
            if rid:
                parent_raw_ids.append(rid)
        if page.success and page.data is not None:
            payload = page.data if isinstance(page.data, dict) else {}
            for report in ncbi_datasets.extract_reports(payload):
                summaries.append(
                    ncbi_datasets.summarize_ortholog_report(report, payload=payload)
                )

    retrieval_complete = bool(pagination.get("retrieval_complete"))
    # Deduplicate by (tax_id, gene_id).
    unique: dict[tuple[str, str], dict[str, Any]] = {}
    membership_counts = {
        "explicitly_verified_count": 0,
        "endpoint_implicit_count": 0,
        "rejected_group_mismatch_count": 0,
        "unverified_excluded_count": 0,
    }
    verified_records: list[dict[str, Any]] = []
    query_human: dict[str, Any] | None = None
    for summary in summaries:
        tax = str(summary.get("tax_id") or "").strip()
        gid = str(summary.get("gene_id") or "").strip()
        if not tax or not gid:
            continue
        key = (tax, gid)
        if key in unique:
            continue
        unique[key] = summary
        membership, reason = classify_membership(summary, query_gene_id=entrez)
        summary = {**summary, "membership_verified": membership, "membership_reason": reason}
        unique[key] = summary
        if gid == entrez and tax == "9606":
            query_human = summary
        if membership == "explicit":
            membership_counts["explicitly_verified_count"] += 1
            verified_records.append(summary)
        elif membership == "endpoint_implicit":
            membership_counts["endpoint_implicit_count"] += 1
            verified_records.append(summary)
        elif membership == "rejected":
            membership_counts["rejected_group_mismatch_count"] += 1
        else:
            membership_counts["unverified_excluded_count"] += 1
    audit["membership"] = {
        **membership_counts,
        "total_unique_records": len(unique),
        "query_human_present": query_human is not None,
    }

    tax_ids = sorted(
        {
            str(rec.get("tax_id")).strip()
            for rec in verified_records
            if rec.get("tax_id") is not None
        }
    )
    membership_map, tax_results, tax_audit = ncbi_taxonomy.resolve_taxonomy_memberships(
        tax_ids,
        scope_tax_id=section_cfg.ortholog_scope_tax_id,
        gene_symbol=gene,
        settings=cfg,
    )
    taxonomy_complete = bool(tax_audit.get("taxonomy_complete"))
    if "taxonomy_complete" not in tax_audit:
        taxonomy_complete = ncbi_taxonomy.taxonomy_retrieval_complete(
            requested_count=int(tax_audit.get("requested_count") or 0),
            resolved_count=int(tax_audit.get("resolved_count") or 0),
            unresolved_count=int(
                tax_audit.get("unresolved_count")
                if tax_audit.get("unresolved_count") is not None
                else tax_audit.get("unknown_count")
                or 0
            ),
            failed_request_count=int(tax_audit.get("failed_request_count") or 0),
        )
        tax_audit = {**tax_audit, "taxonomy_complete": taxonomy_complete}
    audit["taxonomy"] = tax_audit
    for tr in tax_results:
        api, meta = _persist_tool_result_json(
            tr=tr,
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            filename_hint="ncbi-taxonomy-lineage",
        )
        api_runs.append(api)
        if meta:
            raw_meta.append(meta)
            rid = str(meta.get("id") or "").strip()
            if rid:
                parent_raw_ids.append(rid)

    scoped_records: list[dict[str, Any]] = []
    unknown_tax_excluded = 0
    out_of_scope_excluded = 0
    for rec in verified_records:
        tid = str(rec.get("tax_id")).strip()
        member = membership_map.get(tid)
        if member is True:
            scoped_records.append(rec)
        elif member is False:
            out_of_scope_excluded += 1
        else:
            unknown_tax_excluded += 1
    audit["taxonomy"]["unknown_or_out_of_scope_excluded"] = (
        unknown_tax_excluded + out_of_scope_excluded
    )
    audit["taxonomy"]["unknown_tax_excluded"] = unknown_tax_excluded
    audit["taxonomy"]["out_of_scope_excluded"] = out_of_scope_excluded

    scoped_orthologs = [
        rec
        for rec in scoped_records
        if not (
            str(rec.get("gene_id")).strip() == entrez
            and str(rec.get("tax_id")).strip() == "9606"
        )
    ]
    scoped_species = {
        str(rec.get("tax_id")).strip()
        for rec in scoped_orthologs
        if rec.get("tax_id") is not None
    }
    count_claim_ready = _count_claim_ready(
        pagination_complete=retrieval_complete,
        taxonomy_complete=taxonomy_complete,
    )
    count_status = (
        "complete"
        if count_claim_ready and scoped_orthologs is not None
        else ("incomplete" if summaries else "empty")
    )
    if not count_claim_ready:
        count_status = "incomplete" if summaries else count_status
    scoped_ortholog_gene_count: int | None
    if count_claim_ready:
        scoped_ortholog_gene_count = len(scoped_orthologs)
    else:
        scoped_ortholog_gene_count = None
    audit["count_claim_ready"] = count_claim_ready
    audit["taxonomy_complete"] = taxonomy_complete
    audit["pagination_complete"] = retrieval_complete
    species_names = select_species_names(scoped_orthologs, limit=10)
    ncbi_url = ortholog_ncbi_url(entrez_gene_id=entrez)
    legacy_ncbi_url = ortholog_ncbi_legacy_link(
        entrez_gene_id=entrez,
        scope_tax_id=section_cfg.ortholog_scope_tax_id,
        gene_symbol=gene,
    )
    orthodb_url = orthodb.public_gene_url(entrez)
    audit["legacy_ncbi_ortholog_url"] = legacy_ncbi_url

    # OrthoDB supporting source (non-fatal).
    release_tr = orthodb.fetch_release_id(gene_symbol=gene, settings=cfg)
    release_api, release_meta = _persist_tool_result_json(
        tr=release_tr,
        dossier_run_id=run_id,
        gene_symbol=gene,
        settings=cfg,
        persist_db=persist_db,
        filename_hint="orthodb-release",
    )
    api_runs.append(release_api)
    if release_meta:
        raw_meta.append(release_meta)
    search_tr = orthodb.fetch_gene_search(entrez, gene_symbol=gene, settings=cfg)
    search_api, search_meta = _persist_tool_result_json(
        tr=search_tr,
        dossier_run_id=run_id,
        gene_symbol=gene,
        settings=cfg,
        persist_db=persist_db,
        filename_hint=f"orthodb-genesearch-{slugify(entrez)}",
    )
    api_runs.append(search_api)
    if search_meta:
        raw_meta.append(search_meta)
    orthodb_ok = False
    orthodb_diag: dict[str, Any] = {}
    if search_tr.success:
        orthodb_ok, orthodb_diag = orthodb.validate_gene_search_for_human(
            search_tr.data, entrez_gene_id=entrez
        )
    else:
        orthodb_diag = {
            "reason": search_tr.error_type or "request_failed",
            "message": search_tr.error_message,
        }
    audit["orthodb"] = {
        "release": release_tr.data if release_tr.success else None,
        "validated": orthodb_ok,
        "diagnostics": orthodb_diag,
        "public_url": orthodb_url,
    }
    if orthodb_ok:
        od_rec = EvidenceRecord(
            source_id=make_source_id(
                "OrthoDB", gene, AssertionType.gene_identity, f"orthodb-{entrez}"
            ),
            dossier_run_id=run_id,
            gene_symbol=gene,
            official_symbol=gene,
            section=SECTION_HOMOLOGUES,
            subsection=SUBSECTION_1E,
            source_name="OrthoDB",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="orthodb_gene_search",
            species="human",
            taxon_id=9606,
            evidence_grade=EvidenceGrade.C,
            value={
                "entrez_gene_id": entrez,
                "public_url": orthodb_url,
                "release": release_tr.data if release_tr.success else None,
                "diagnostics": orthodb_diag,
            },
            display_text=f"{gene} OrthoDB supporting orthology record validated.",
            api_run_id=search_api.id,
        )
        _append_evidence(evidence, od_rec, persist_db=persist_db)

    table_status = TABLE_UNAVAILABLE
    capture_rec = None
    capture_origin: str | None = None
    displayed_gene_count: int | None = None
    human_reference_row_detected: bool | None = None
    if not skip_table_capture:
        (
            cap_api,
            cap_meta,
            capture_rec,
            cap_audit,
            captured_bytes,
        ) = _capture_ncbi_ortholog_table(
            dossier_run_id=run_id,
            gene_symbol=gene,
            entrez_gene_id=entrez,
            scope_tax_id=section_cfg.ortholog_scope_tax_id,
            page_url=ncbi_url,
            max_visible_rows=section_cfg.max_visible_rows,
            scoped_records=scoped_records,
            scoped_ortholog_gene_count=scoped_ortholog_gene_count,
            retrieval_complete=retrieval_complete,
            count_claim_ready=count_claim_ready,
            parent_raw_artifact_ids=parent_raw_ids,
            settings=cfg,
            persist_db=persist_db,
            taxon_name=section_cfg.ncbi_taxon_name,
        )
        api_runs.append(cap_api)
        if cap_meta:
            raw_meta.append(cap_meta)
        audit["viewer_capture"] = cap_audit
        capture_metadata = None
        if isinstance(cap_audit, dict):
            capture_metadata = cap_audit.get("capture_metadata")
            if capture_metadata is None and capture_rec is not None:
                capture_metadata = {
                    "capture_origin": (capture_rec.value or {}).get("capture_origin"),
                    "sha256": (capture_rec.value or {}).get("sha256"),
                    "entrez_gene_id": (capture_rec.value or {}).get("entrez_gene_id"),
                    "selected_scope_tax_id": (capture_rec.value or {}).get(
                        "selected_scope_tax_id"
                    ),
                    "visible_headers": (capture_rec.value or {}).get("visible_headers"),
                    "displayed_gene_count": (capture_rec.value or {}).get(
                        "displayed_gene_count"
                    ),
                    "expected_displayed_count": (capture_rec.value or {}).get(
                        "expected_displayed_count"
                    ),
                }
            displayed_gene_count = cap_audit.get("displayed_gene_count")
            human_reference_row_detected = cap_audit.get(
                "human_reference_row_detected"
            )
            capture_origin = cap_audit.get("capture_origin")
        # Prefer official live capture whenever Playwright succeeded with live
        # provenance. Hard count gate applies only when count_claim_ready.
        capture_is_official = _official_capture_gate(
            capture_api_success=bool(cap_api.success),
            capture_metadata=capture_metadata if isinstance(capture_metadata, dict) else None,
            captured_bytes=captured_bytes,
            resolved_entrez_gene_id=entrez,
            configured_scope_tax_id=section_cfg.ortholog_scope_tax_id,
            retrieval_complete=count_claim_ready,
        )
        if capture_is_official and capture_rec is not None:
            _append_evidence(evidence, capture_rec, persist_db=persist_db)
            table_status = TABLE_OFFICIAL
        elif capture_rec is not None:
            # Live bytes existed but failed official gate — do not promote to official.
            audit["viewer_capture"] = {
                **(cap_audit or {}),
                "official_gate_rejected": True,
            }
            capture_rec = None
            capture_origin = None
    else:
        audit["viewer_capture"] = {"status": "skipped", "reason": "skip_table_capture"}

    if table_status != TABLE_OFFICIAL:
        if count_claim_ready and scoped_orthologs:
            table_status = TABLE_COMPLETE_FALLBACK
        elif scoped_orthologs:
            table_status = TABLE_PARTIAL_FALLBACK
        else:
            table_status = TABLE_UNAVAILABLE
    audit["table_status"] = table_status

    if count_claim_ready and scoped_ortholog_gene_count is not None:
        narrative = build_complete_narrative(
            gene_symbol=gene,
            scope_label=section_cfg.ortholog_scope_label,
            ortholog_gene_count=scoped_ortholog_gene_count,
            species_names=species_names,
        )
    else:
        narrative = build_incomplete_narrative(
            gene_symbol=gene, scope_label=section_cfg.ortholog_scope_label
        )

    summary_value = {
        "entrez_gene_id": entrez,
        "gene_symbol": gene,
        "scope_tax_id": section_cfg.ortholog_scope_tax_id,
        "scope_label": section_cfg.ortholog_scope_label,
        "retrieval_complete": retrieval_complete,
        "taxonomy_complete": taxonomy_complete,
        "count_claim_ready": count_claim_ready,
        "count_status": count_status,
        "total_unique_records": len(unique),
        "scoped_unique_records": len(scoped_records),
        "scoped_ortholog_gene_count": scoped_ortholog_gene_count,
        "scoped_species_count": len(scoped_species),
        "species_names": species_names,
        "ncbi_url": ncbi_url,
        "orthodb_url": orthodb_url,
        "table_status": table_status,
        "capture_origin": capture_origin,
        "displayed_gene_count": displayed_gene_count,
        "human_reference_row_detected": human_reference_row_detected,
        "narrative": narrative,
        "fallback_rows": [
            {
                "species": rec.get("common_name")
                or rec.get("scientific_name")
                or rec.get("taxname"),
                "gene": rec.get("symbol"),
                "description": rec.get("description"),
                "gene_id": rec.get("gene_id"),
                "tax_id": rec.get("tax_id"),
            }
            for rec in scoped_orthologs[: section_cfg.max_visible_rows]
        ],
        "query_human_excluded_from_count": True,
        "presentation_item_key": f"orthologs-{slugify(gene)}",
    }
    summary_rec = _record(
        dossier_run_id=run_id,
        gene_symbol=gene,
        fact_type="ortholog_collection_summary",
        key=f"summary-{entrez}-{section_cfg.ortholog_scope_tax_id}",
        value=summary_value,
        display_text=narrative,
        raw_artifact_id=parent_raw_ids[0] if parent_raw_ids else None,
    )
    _append_evidence(evidence, summary_rec, persist_db=persist_db)

    rendering_status = {
        "overall": "success"
        if table_status in {TABLE_OFFICIAL, TABLE_COMPLETE_FALLBACK}
        or (count_claim_ready and scoped_ortholog_gene_count is not None)
        else ("partial" if summaries else "empty"),
        "count_status": count_status,
        "table_status": table_status,
        "retrieval_complete": retrieval_complete,
        "taxonomy_complete": taxonomy_complete,
        "count_claim_ready": count_claim_ready,
        "scoped_ortholog_gene_count": scoped_ortholog_gene_count,
        "scoped_species_count": len(scoped_species),
    }
    _ = transient
    return {
        **state,
        "evidence_records": evidence,
        "api_runs": api_runs,
        "raw_artifacts": raw_meta,
        "errors": errors,
        "coverage": coverage_extra,
        "section_1e_status": {
            "config": {
                "ortholog_scope_tax_id": section_cfg.ortholog_scope_tax_id,
                "ortholog_scope_label": section_cfg.ortholog_scope_label,
                "max_visible_rows": section_cfg.max_visible_rows,
            },
            "rendering_status": rendering_status,
            "summary": summary_value,
            "audit": audit,
        },
    }


__all__ = [
    "SUPPORTED_SECTION_1E_SCOPES",
    "MODEL_SPECIES_PRIORITY",
    "TABLE_OFFICIAL",
    "TABLE_COMPLETE_FALLBACK",
    "TABLE_PARTIAL_FALLBACK",
    "TABLE_UNAVAILABLE",
    "REQUIRED_HEADERS",
    "CAPTURE_ORIGIN_LIVE",
    "Section1eConfig",
    "ortholog_ncbi_url",
    "ortholog_ncbi_legacy_link",
    "classify_membership",
    "select_species_names",
    "build_complete_narrative",
    "build_incomplete_narrative",
    "node_generate_section_1e_derived_artifacts",
]
