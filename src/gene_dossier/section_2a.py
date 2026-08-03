"""Bundle-only Section 2a tissue-specific expression (GTEx + HBT).

Owns GTEx V8 sample/median/tissue-metadata requests and the official Human
Brain Transcriptome whole-brain PDF. Not wired through the generic GTEx client
path when Section 2a is selected.
"""

from __future__ import annotations

import hashlib
import io
import logging
import math
import statistics
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Sequence

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceCoverageResult,
    SourceStatus,
    SourceType,
    ToolResult,
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
from gene_dossier.tools import gtex, hbt
from gene_dossier.workflow import DossierState, WorkflowTransientContext

logger = logging.getLogger(__name__)

SECTION_EXPRESSION = "Tissue and cell expression"
SUBSECTION_2A = "Tissue-specific information"

GTEX_PORTAL_GENE_URL = "https://www.gtexportal.org/home/gene/{gene}"
HBT_HOME_URL = hbt.HBT_HOME

EXPECTED_GTEX_V8_BRAIN_TISSUE_COUNT = 13
EXPECTED_GTEX_V8_BRAIN_TISSUES: tuple[str, ...] = (
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
)

MEDIAN_REL_TOL = 1e-3
MEDIAN_ABS_TOL = 0.01
JITTER_SEED_SALT = "gene-dossier-section-2a-v1"
PLOT_VERSION = "section_2a_gtex_violin_v1"

FALLBACK_COLOR = "#AAAAAA"
BRAIN_FALLBACK_COLOR = "#EEE8AA"
IQR_COLOR = "#4A4A4A"
MEDIAN_COLOR = "#FFFFFF"
POINT_ALPHA = 0.35


@dataclass(frozen=True)
class Section2aConfig:
    gtex_dataset_id: str = gtex.DEFAULT_DATASET
    gtex_genome_build: str = gtex.DEFAULT_GENOME_BUILD
    gtex_items_per_page: int = gtex.DEFAULT_ITEMS_PER_PAGE
    brain_tissue_prefix: str = "Brain_"
    hbt_base_url: str = hbt.HBT_BASE
    plot_dpi: int = 180
    hbt_raster_dpi: int = 180

    def __post_init__(self) -> None:
        if self.gtex_dataset_id != gtex.DEFAULT_DATASET:
            raise ValueError(
                f"Section 2a requires dataset {gtex.DEFAULT_DATASET!r}; "
                f"got {self.gtex_dataset_id!r}"
            )
        if self.gtex_genome_build != gtex.DEFAULT_GENOME_BUILD:
            raise ValueError(
                f"Section 2a requires genome build {gtex.DEFAULT_GENOME_BUILD!r}; "
                f"got {self.gtex_genome_build!r}"
            )
        if int(self.gtex_items_per_page) < 1:
            raise ValueError("gtex_items_per_page must be >= 1")


def gtex_intro_text(gene_symbol: str) -> str:
    gene = (gene_symbol or "").strip() or "this gene"
    return (
        "The Genotype-Tissue Expression (GTEx) project is a comprehensive public "
        "resource that provides tissue-specific gene expression and regulation data, "
        "including RNA-Seq, eQTLs, and histology images. The resource contains "
        "information collected from 53 non-diseased tissue sites across nearly 1000 "
        "individuals, utilizing primarily molecular assays such as whole genome/exome "
        "sequencing and RNA-Seq. Below is information pertaining to "
        f"{gene} expression across all tissues and specifically in the brain."
    )


def hbt_intro_text() -> str:
    return (
        "The Human Brain Transcriptome (HBT) is a public database hosted by the "
        "Department of Neurobiology at Yale University School of Medicine. It contains "
        "transcriptome data and associated metadata from over 1,340 tissue samples, "
        "specifically focusing on the developing and adult human brain."
    )


def hbt_link_text(gene_symbol: str) -> str:
    gene = (gene_symbol or "").strip() or "Gene"
    return (
        f"{gene} gene expression (Link) along entire development and adulthood in the "
        "cerebellar cortex (CBC), mediodorsal nucleus of the thalamus (MD), striatum "
        "(STR), amygdala (AMY), hippocampus (HIP) and 11 areas of neocortex (NCX) is "
        "shown below."
    )


def gtex_gene_url(gene_symbol: str) -> str:
    return GTEX_PORTAL_GENE_URL.format(gene=(gene_symbol or "").strip())


def hbt_pdf_url(gene_symbol: str) -> str:
    return hbt.whole_brain_pdf_url(gene_symbol)


def _normalize_entrez_gene_id(value: Any) -> str | None:
    """Normalize numeric Entrez Gene IDs; reject bools and non-digit text."""
    if value is None or isinstance(value, bool):
        return None
    text = str(value).strip()
    if not text or not text.isdigit():
        return None
    return str(int(text))


def _finite_nonnegative(values: Sequence[Any]) -> list[float] | None:
    out: list[float] = []
    for raw in values:
        try:
            number = float(raw)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(number) or number < 0:
            return None
        out.append(number)
    return out


def _percentile(sorted_vals: Sequence[float], q: float) -> float:
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    pos = (len(sorted_vals) - 1) * q
    low = int(math.floor(pos))
    high = int(math.ceil(pos))
    if low == high:
        return float(sorted_vals[low])
    weight = pos - low
    return float(sorted_vals[low] * (1.0 - weight) + sorted_vals[high] * weight)


def tissue_stats(sample_values: Sequence[float]) -> dict[str, Any]:
    values = [float(v) for v in sample_values]
    ordered = sorted(values)
    return {
        "sample_count": len(values),
        "median_tpm": float(statistics.median(values)) if values else None,
        "q1_tpm": _percentile(ordered, 0.25) if values else None,
        "q3_tpm": _percentile(ordered, 0.75) if values else None,
        "min_tpm": float(min(values)) if values else None,
        "max_tpm": float(max(values)) if values else None,
        "unit": "TPM",
    }


def _display_label_from_id(tissue_id: str) -> str:
    text = tissue_id.replace("_", " ")
    return text


def parse_sample_expression_rows(
    rows: Sequence[Any],
    *,
    expected_gencode_id: str,
    expected_gene_symbol: str,
    expected_dataset_id: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Validate sample-level rows; return (valid_tissues, diagnostics)."""
    valid: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    seen: set[str] = set()
    gene_u = expected_gene_symbol.strip().upper()
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            diagnostics.append({"index": index, "reason": "non_object_row"})
            continue
        tissue_id = str(row.get("tissueSiteDetailId") or "").strip()
        if not tissue_id:
            diagnostics.append({"index": index, "reason": "empty_tissue_id"})
            continue
        if tissue_id in seen:
            diagnostics.append(
                {
                    "index": index,
                    "tissue_site_detail_id": tissue_id,
                    "reason": "duplicate_tissue_id",
                }
            )
            continue
        gencode = str(row.get("gencodeId") or expected_gencode_id).strip()
        symbol = str(row.get("geneSymbol") or expected_gene_symbol).strip()
        dataset = str(row.get("datasetId") or expected_dataset_id).strip()
        unit = str(row.get("unit") or "TPM").strip()
        if gencode and gencode != expected_gencode_id:
            diagnostics.append(
                {
                    "tissue_site_detail_id": tissue_id,
                    "reason": "gencode_mismatch",
                    "gencode_id": gencode,
                }
            )
            continue
        if symbol and symbol.upper() != gene_u:
            diagnostics.append(
                {
                    "tissue_site_detail_id": tissue_id,
                    "reason": "gene_symbol_mismatch",
                    "gene_symbol": symbol,
                }
            )
            continue
        if dataset and dataset != expected_dataset_id:
            diagnostics.append(
                {
                    "tissue_site_detail_id": tissue_id,
                    "reason": "dataset_mismatch",
                    "dataset_id": dataset,
                }
            )
            continue
        if unit.upper() != "TPM":
            diagnostics.append(
                {
                    "tissue_site_detail_id": tissue_id,
                    "reason": "unit_not_tpm",
                    "unit": unit,
                }
            )
            continue
        raw_values = row.get("data")
        if not isinstance(raw_values, list):
            diagnostics.append(
                {
                    "tissue_site_detail_id": tissue_id,
                    "reason": "missing_sample_array",
                }
            )
            continue
        samples = _finite_nonnegative(raw_values)
        if samples is None:
            diagnostics.append(
                {
                    "tissue_site_detail_id": tissue_id,
                    "reason": "invalid_sample_values",
                }
            )
            continue
        seen.add(tissue_id)
        stats = tissue_stats(samples)
        valid.append(
            {
                "tissue_site_detail_id": tissue_id,
                "ontology_id": row.get("ontologyId"),
                "sample_values": samples,
                "unit": "TPM",
                "dataset_id": expected_dataset_id,
                "gencode_id": expected_gencode_id,
                "gene_symbol": expected_gene_symbol,
                "subset_group": row.get("subsetGroup"),
                **stats,
            }
        )
    return valid, diagnostics


def parse_median_rows(rows: Sequence[Any]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        tissue_id = str(row.get("tissueSiteDetailId") or "").strip()
        if not tissue_id:
            continue
        median = row.get("median")
        try:
            median_f = float(median)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(median_f):
            continue
        out[tissue_id] = {
            "tissue_site_detail_id": tissue_id,
            "ontology_id": row.get("ontologyId"),
            "median": median_f,
            "unit": row.get("unit") or "TPM",
            "dataset_id": row.get("datasetId"),
            "gencode_id": row.get("gencodeId"),
            "gene_symbol": row.get("geneSymbol"),
        }
    return out


def parse_tissue_metadata(rows: Sequence[Any]) -> dict[str, dict[str, Any]]:
    """Preserve metadata order via insertion order of the returned dict."""
    out: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        tissue_id = str(row.get("tissueSiteDetailId") or "").strip()
        if not tissue_id or tissue_id in out:
            continue
        color = row.get("colorHex")
        if isinstance(color, str) and color and not color.startswith("#"):
            color = f"#{color}"
        out[tissue_id] = {
            "tissue_site_detail_id": tissue_id,
            "tissue_site": row.get("tissueSite"),
            "tissue_site_detail": row.get("tissueSiteDetail"),
            "ontology_id": row.get("ontologyId"),
            "color_hex": color,
            "color_rgb": row.get("colorRgb"),
            "dataset_id": row.get("datasetId"),
        }
    return out


def validate_medians(
    tissues: Sequence[dict[str, Any]],
    api_medians: dict[str, dict[str, Any]],
    *,
    rel_tol: float = MEDIAN_REL_TOL,
    abs_tol: float = MEDIAN_ABS_TOL,
) -> dict[str, Any]:
    matched = 0
    mismatched = 0
    missing = 0
    max_abs_diff = 0.0
    mismatches: list[dict[str, Any]] = []
    for tissue in tissues:
        tissue_id = str(tissue["tissue_site_detail_id"])
        calc = tissue.get("median_tpm")
        api_row = api_medians.get(tissue_id)
        if api_row is None or calc is None:
            missing += 1
            continue
        api_median = float(api_row["median"])
        diff = abs(float(calc) - api_median)
        max_abs_diff = max(max_abs_diff, diff)
        ok = diff <= max(abs_tol, rel_tol * max(abs(api_median), abs(float(calc))))
        if ok:
            matched += 1
        else:
            mismatched += 1
            mismatches.append(
                {
                    "tissue_site_detail_id": tissue_id,
                    "calculated_median": calc,
                    "api_median": api_median,
                    "abs_diff": diff,
                }
            )
    systematic_fail = mismatched > 0 and mismatched >= max(3, matched // 5)
    return {
        "matched_tissue_count": matched,
        "mismatched_tissue_count": mismatched,
        "missing_median_count": missing,
        "maximum_absolute_difference": max_abs_diff,
        "mismatches": mismatches[:50],
        "rel_tol": rel_tol,
        "abs_tol": abs_tol,
        "systematic_mismatch": systematic_fail,
    }


def order_tissues(
    tissues: Sequence[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    by_id = {str(t["tissue_site_detail_id"]): t for t in tissues}
    ordered: list[dict[str, Any]] = []
    for tissue_id in metadata:
        if tissue_id in by_id:
            ordered.append(by_id.pop(tissue_id))
    for tissue_id in sorted(by_id):
        ordered.append(by_id[tissue_id])
    return ordered


def brain_subset(
    tissues: Sequence[dict[str, Any]],
    *,
    prefix: str = "Brain_",
) -> list[dict[str, Any]]:
    return [
        t
        for t in tissues
        if str(t.get("tissue_site_detail_id") or "").startswith(prefix)
    ]


def enrich_tissue_display(
    tissues: Sequence[dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
    *,
    color_fallback: str = FALLBACK_COLOR,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for tissue in tissues:
        tissue_id = str(tissue["tissue_site_detail_id"])
        meta = metadata.get(tissue_id) or {}
        label = (
            str(meta.get("tissue_site_detail") or "").strip()
            or _display_label_from_id(tissue_id)
        )
        color = str(meta.get("color_hex") or "").strip() or color_fallback
        row = dict(tissue)
        row["display_label"] = label
        row["color_hex"] = color
        row["ontology_id"] = row.get("ontology_id") or meta.get("ontology_id")
        out.append(row)
    return out


def _jitter_offsets(n: int, *, tissue_id: str) -> list[float]:
    seed_material = f"{JITTER_SEED_SALT}:{tissue_id}:{n}".encode("utf-8")
    seed = int(hashlib.sha256(seed_material).hexdigest()[:8], 16)
    # Simple LCG for deterministic jitter in [-0.18, 0.18]
    offsets: list[float] = []
    state = seed
    for _ in range(n):
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        offsets.append((state / 0x7FFFFFFF) * 0.36 - 0.18)
    return offsets


def render_gtex_violin_png(
    tissues: Sequence[dict[str, Any]],
    *,
    gene_symbol: str,
    gencode_id: str,
    dpi: int = 180,
    figsize: tuple[float, float] = (11.0, 5.2),
) -> bytes:
    """Deterministic matplotlib violin plot from sample TPM arrays."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if not tissues:
        raise ValueError("no tissues to plot")

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    positions = list(range(1, len(tissues) + 1))
    data = [list(t["sample_values"]) for t in tissues]
    parts = ax.violinplot(
        data,
        positions=positions,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        widths=0.85,
    )
    for body, tissue in zip(parts["bodies"], tissues):
        color = tissue.get("color_hex") or FALLBACK_COLOR
        body.set_facecolor(color)
        body.set_edgecolor("#333333")
        body.set_alpha(0.85)
        body.set_linewidth(0.4)

    for index, tissue in enumerate(tissues):
        pos = positions[index]
        values = list(tissue["sample_values"])
        q1 = float(tissue["q1_tpm"])
        q3 = float(tissue["q3_tpm"])
        med = float(tissue["median_tpm"])
        box_width = 0.18
        ax.add_patch(
            Rectangle(
                (pos - box_width / 2.0, q1),
                box_width,
                max(q3 - q1, 1e-9),
                facecolor=IQR_COLOR,
                edgecolor=IQR_COLOR,
                linewidth=0.6,
                zorder=3,
            )
        )
        ax.hlines(
            med,
            pos - box_width / 2.0,
            pos + box_width / 2.0,
            colors=MEDIAN_COLOR,
            linewidth=1.6,
            zorder=4,
        )
        offsets = _jitter_offsets(len(values), tissue_id=str(tissue["tissue_site_detail_id"]))
        xs = [pos + off for off in offsets]
        ax.scatter(
            xs,
            values,
            s=4,
            c=tissue.get("color_hex") or FALLBACK_COLOR,
            alpha=POINT_ALPHA,
            linewidths=0,
            zorder=2,
        )

    ax.set_ylabel("TPM")
    ax.set_title(
        f"Bulk tissue gene expression for {gene_symbol} ({gencode_id})",
        fontsize=10,
        pad=8,
    )
    labels = [str(t.get("display_label") or t["tissue_site_detail_id"]) for t in tissues]
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=6)
    ax.set_xlim(0.3, len(tissues) + 0.7)
    ymax = max((float(t.get("max_tpm") or 0.0) for t in tissues), default=1.0)
    ax.set_ylim(bottom=0.0, top=ymax * 1.08 if ymax > 0 else 1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor="white")
    plt.close(fig)
    return buf.getvalue()


def _evidence(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    confidence_notes: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=make_source_id(source_name, gene_symbol, AssertionType.expression, key),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_EXPRESSION,
        subsection=SUBSECTION_2A,
        source_name=source_name,
        source_type=SourceType.expression_database,
        assertion_type=AssertionType.expression,
        fact_type=fact_type,
        organism="Homo sapiens",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.B,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
        confidence_notes=confidence_notes,
    )


def _validate_pdf(content: bytes) -> dict[str, Any]:
    if not content.startswith(b"%PDF"):
        raise ValueError("PDF magic bytes missing")
    if len(content) < 500:
        raise ValueError("PDF too small")
    return {"media_type": "application/pdf", "byte_size": len(content)}


def _collect_paginated(
    *,
    fetch_page,
    dossier_run_id: str,
    gene_symbol: str,
    settings: Settings,
    persist_db: bool,
    filename_hint: str,
    api_runs: list[Any],
    raw_meta: list[dict[str, Any]],
) -> tuple[list[Any], list[ToolResult], dict[str, Any]]:
    rows: list[Any] = []
    results: list[ToolResult] = []
    audit: dict[str, Any] = {"pages": [], "status": "success"}
    for page_result in gtex.iter_paginated(fetch_page):
        results.append(page_result)
        api, meta = _persist_tool_result_json(
            tr=page_result,
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            settings=settings,
            persist_db=persist_db,
            filename_hint=filename_hint,
        )
        api_runs.append(api)
        page_info = {
            "success": page_result.success,
            "status_code": page_result.status_code,
            "request_url": page_result.request_url,
            "api_run_id": api.id,
            "raw_artifact_id": meta.get("id") if meta else None,
            "error_type": page_result.error_type,
            "error_message": page_result.error_message,
        }
        audit["pages"].append(page_info)
        if not page_result.success:
            audit["status"] = "failed"
            audit["error"] = page_result.error_message
            break
        if meta:
            raw_meta.append(meta)
        rows.extend(gtex._data_list(page_result.data))
    return rows, results, audit


def node_generate_section_2a_derived_artifacts(
    state: DossierState,
    *,
    settings: Settings | None = None,
    persist_db: bool = True,
    transient: WorkflowTransientContext | None = None,
    config: Section2aConfig | None = None,
) -> DossierState:
    """Sole Section 2a network owner: GTEx V8 + official HBT PDF."""
    _ = transient
    if state.get("run_type") != "section_bundle" or "2a" not in (
        state.get("selected_section_keys") or []
    ):
        return state

    cfg = settings or get_settings()
    section_cfg = config or Section2aConfig()
    run_id = state["dossier_run_id"]
    gene = state["gene_symbol"]
    gene_ids = dict(state.get("gene_ids") or {})
    evidence = list(state.get("evidence_records") or [])
    api_runs = list(state.get("api_runs") or [])
    raw_meta = list(state.get("raw_artifacts") or [])
    errors = list(state.get("errors") or [])
    coverage_extra = list(state.get("coverage") or [])

    audit: dict[str, Any] = {
        "config": {
            "gtex_dataset_id": section_cfg.gtex_dataset_id,
            "gtex_genome_build": section_cfg.gtex_genome_build,
            "gtex_items_per_page": section_cfg.gtex_items_per_page,
            "brain_tissue_prefix": section_cfg.brain_tissue_prefix,
            "plot_dpi": section_cfg.plot_dpi,
            "hbt_raster_dpi": section_cfg.hbt_raster_dpi,
            "plot_version": PLOT_VERSION,
        },
        "network_owner": "section_2a",
        "gtex": {},
        "hbt": {},
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }

    gtex_status = {
        "resolve": "pending",
        "sample_expression": "pending",
        "median_expression": "pending",
        "tissue_metadata": "pending",
        "all_tissues_figure": "pending",
        "brain_figure": "pending",
        "color_status": "official",
    }
    hbt_status = {"pdf": "pending", "figure": "pending"}

    # --- GTEx resolve ---
    resolve_tr = gtex.resolve_gene(
        gene,
        genome_build=section_cfg.gtex_genome_build,
        settings=cfg,
        require_exact_unambiguous=True,
    )
    resolve_api, resolve_meta = _persist_tool_result_json(
        tr=resolve_tr,
        dossier_run_id=run_id,
        gene_symbol=gene,
        settings=cfg,
        persist_db=persist_db,
        filename_hint=f"gtex-resolve-{slugify(gene)}",
    )
    api_runs.append(resolve_api)
    if resolve_meta:
        raw_meta.append(resolve_meta)
    audit["gtex"]["resolve"] = {
        "success": resolve_tr.success,
        "error_type": resolve_tr.error_type,
        "error_message": resolve_tr.error_message,
        "api_run_id": resolve_api.id,
        "raw_artifact_id": resolve_meta.get("id") if resolve_meta else None,
    }

    gencode_id: str | None = None
    resolve_data = resolve_tr.data if isinstance(resolve_tr.data, dict) else {}
    if resolve_tr.success:
        candidate_gencode_id = str(resolve_data.get("gencode_id") or "") or None
        identity_entrez = _normalize_entrez_gene_id(gene_ids.get("entrez_gene_id"))
        gtex_entrez = _normalize_entrez_gene_id(resolve_data.get("entrez_gene_id"))
        if identity_entrez is not None and gtex_entrez is not None:
            identity_status = (
                "match" if identity_entrez == gtex_entrez else "mismatch"
            )
        else:
            identity_status = "not_comparable"
        entrez_identity_mismatch = identity_status == "mismatch"
        audit["gtex"]["resolve"]["identity_check"] = {
            "status": identity_status,
            "identity_entrez": identity_entrez,
            "gtex_entrez": gtex_entrez,
        }

        if entrez_identity_mismatch:
            audit["gtex"]["resolve"]["candidate_gencode_id"] = candidate_gencode_id
            gencode_id = None
            gtex_status["resolve"] = "identity_mismatch"
            gtex_status["sample_expression"] = "skipped_identity_mismatch"
            gtex_status["median_expression"] = "skipped_identity_mismatch"
            gtex_status["tissue_metadata"] = "skipped_identity_mismatch"
            gtex_status["all_tissues_figure"] = "unavailable"
            gtex_status["brain_figure"] = "unavailable"
            gtex_status["color_status"] = "unavailable"
            audit["gtex"]["identity_gate"] = {
                "status": "failed",
                "reason": "entrez_gene_id_mismatch",
                "downstream_requests_skipped": [
                    "gene_expression",
                    "median_expression",
                    "tissue_site_detail",
                ],
            }
            coverage_extra.append(
                SourceCoverageResult(
                    dossier_run_id=run_id,
                    source_name=gtex.SOURCE_NAME,
                    status=SourceStatus.failed,
                    evidence_record_count=0,
                    error_message=(
                        f"GTEx Entrez Gene ID mismatch: dossier={identity_entrez}, "
                        f"GTEx={gtex_entrez}; downstream GTEx retrieval was skipped."
                    ),
                    report_sections_supported=[SUBSECTION_2A],
                )
            )
        else:
            gencode_id = candidate_gencode_id
            if gencode_id:
                gtex_status["resolve"] = "success"
                evidence_entrez = identity_entrez or gtex_entrez
                ref_rec = _evidence(
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    source_name=gtex.SOURCE_NAME,
                    fact_type="gtex_gene_reference",
                    key=f"{gencode_id}-reference",
                    value={
                        "gene_symbol": gene,
                        "entrez_gene_id": evidence_entrez,
                        "gencode_id": gencode_id,
                        "genome_build": resolve_data.get("genome_build")
                        or section_cfg.gtex_genome_build,
                        "chromosome": resolve_data.get("chromosome"),
                        "start": resolve_data.get("start"),
                        "end": resolve_data.get("end"),
                        "strand": resolve_data.get("strand"),
                        "gene_type": resolve_data.get("gene_type"),
                        "dataset_id": section_cfg.gtex_dataset_id,
                    },
                    display_text=f"{gene} GTEx GENCODE ID is {gencode_id}.",
                    api_run_id=resolve_api.id,
                    raw_artifact_id=resolve_meta.get("id") if resolve_meta else None,
                    confidence_notes="GTEx is human-only.",
                )
                _append_evidence(evidence, ref_rec, persist_db=persist_db)
            else:
                gtex_status["resolve"] = "failed"
    else:
        gtex_status["resolve"] = "failed"
        coverage_extra.append(
            SourceCoverageResult(
                dossier_run_id=run_id,
                source_name=gtex.SOURCE_NAME,
                status=SourceStatus.failed,
                evidence_record_count=0,
                error_message=resolve_tr.error_message or "GTEx gene resolve failed",
                report_sections_supported=[SUBSECTION_2A],
            )
        )

    tissues: list[dict[str, Any]] = []
    metadata: dict[str, dict[str, Any]] = {}
    median_validation: dict[str, Any] = {}
    sample_parent_raw_ids: list[str] = []
    sample_parent_api_ids: list[str] = []

    if gencode_id:
        # Sample expression (primary plot input)
        sample_rows, _sample_trs, sample_audit = _collect_paginated(
            fetch_page=lambda page: gtex.gene_expression(
                gencode_id,
                gene_symbol=gene,
                dataset_id=section_cfg.gtex_dataset_id,
                page=page,
                items_per_page=section_cfg.gtex_items_per_page,
                settings=cfg,
            ),
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            filename_hint=f"gtex-gene-expression-{slugify(gencode_id)}",
            api_runs=api_runs,
            raw_meta=raw_meta,
        )
        audit["gtex"]["sample_expression"] = sample_audit
        sample_parent_raw_ids = [
            p["raw_artifact_id"]
            for p in sample_audit.get("pages") or []
            if p.get("raw_artifact_id")
        ]
        sample_parent_api_ids = [
            p["api_run_id"] for p in sample_audit.get("pages") or [] if p.get("api_run_id")
        ]
        if sample_audit.get("status") == "success":
            tissues, sample_diags = parse_sample_expression_rows(
                sample_rows,
                expected_gencode_id=gencode_id,
                expected_gene_symbol=gene,
                expected_dataset_id=section_cfg.gtex_dataset_id,
            )
            audit["gtex"]["sample_expression"]["parse_diagnostics"] = sample_diags[:50]
            gtex_status["sample_expression"] = "success" if tissues else "failed"
            if not tissues:
                audit["gtex"]["sample_expression"]["error"] = "no valid tissue sample arrays"
        else:
            gtex_status["sample_expression"] = "failed"

        # API medians (validation)
        median_rows, _median_trs, median_audit = _collect_paginated(
            fetch_page=lambda page: gtex.median_expression(
                gencode_id,
                gene_symbol=gene,
                dataset_id=section_cfg.gtex_dataset_id,
                page=page,
                items_per_page=section_cfg.gtex_items_per_page,
                settings=cfg,
            ),
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            filename_hint=f"gtex-median-expression-{slugify(gencode_id)}",
            api_runs=api_runs,
            raw_meta=raw_meta,
        )
        audit["gtex"]["median_expression"] = median_audit
        api_medians: dict[str, dict[str, Any]] = {}
        if median_audit.get("status") == "success":
            api_medians = parse_median_rows(median_rows)
            gtex_status["median_expression"] = "success"
            if tissues:
                median_validation = validate_medians(tissues, api_medians)
                audit["gtex"]["median_validation"] = median_validation
                if median_validation.get("systematic_mismatch"):
                    gtex_status["all_tissues_figure"] = "failed"
                    gtex_status["brain_figure"] = "failed"
                    audit["gtex"]["figure_gate"] = "systematic_median_mismatch"
        else:
            gtex_status["median_expression"] = "unavailable"
            median_validation = {
                "matched_tissue_count": 0,
                "mismatched_tissue_count": 0,
                "missing_median_count": len(tissues),
                "api_median_unavailable": True,
            }
            audit["gtex"]["median_validation"] = median_validation

        # Tissue metadata once
        meta_rows, _meta_trs, meta_audit = _collect_paginated(
            fetch_page=lambda page: gtex.tissue_site_detail(
                dataset_id=section_cfg.gtex_dataset_id,
                page=page,
                items_per_page=section_cfg.gtex_items_per_page,
                gene_symbol=gene,
                settings=cfg,
            ),
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            filename_hint=f"gtex-tissue-site-detail-{section_cfg.gtex_dataset_id}",
            api_runs=api_runs,
            raw_meta=raw_meta,
        )
        audit["gtex"]["tissue_metadata"] = meta_audit
        if meta_audit.get("status") == "success":
            metadata = parse_tissue_metadata(meta_rows)
            gtex_status["tissue_metadata"] = "success"
            gtex_status["color_status"] = "official"
        else:
            gtex_status["tissue_metadata"] = "unavailable"
            gtex_status["color_status"] = "metadata_color_fallback"
            audit["gtex"]["tissue_metadata"]["color_status"] = "metadata_color_fallback"

    # Order + enrich tissues; emit compact summaries
    ordered_all: list[dict[str, Any]] = []
    ordered_brain: list[dict[str, Any]] = []
    if tissues and gtex_status.get("all_tissues_figure") != "failed":
        ordered_all = enrich_tissue_display(
            order_tissues(tissues, metadata),
            metadata,
            color_fallback=(
                FALLBACK_COLOR
                if gtex_status["color_status"] == "official"
                else FALLBACK_COLOR
            ),
        )
        ordered_brain = brain_subset(
            ordered_all, prefix=section_cfg.brain_tissue_prefix
        )
        for brain_row in ordered_brain:
            if gtex_status["color_status"] != "official" and not (
                metadata.get(str(brain_row["tissue_site_detail_id"])) or {}
            ).get("color_hex"):
                brain_row["color_hex"] = BRAIN_FALLBACK_COLOR

        audit["gtex"]["tissue_counts"] = {
            "total_tissue_count": len(ordered_all),
            "brain_tissue_count": len(ordered_brain),
            "total_sample_value_count": sum(
                int(t.get("sample_count") or 0) for t in ordered_all
            ),
            "brain_tissue_ids": [
                t["tissue_site_detail_id"] for t in ordered_brain
            ],
            "expected_v8_brain_count": EXPECTED_GTEX_V8_BRAIN_TISSUE_COUNT,
        }

        parent_api = sample_parent_api_ids[0] if sample_parent_api_ids else None
        parent_raw = sample_parent_raw_ids[0] if sample_parent_raw_ids else None
        for tissue in ordered_all:
            tissue_id = str(tissue["tissue_site_detail_id"])
            rec = _evidence(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name=gtex.SOURCE_NAME,
                fact_type="gtex_tissue_expression_summary",
                key=f"{gencode_id}-{tissue_id}",
                value={
                    "tissue_site_detail_id": tissue_id,
                    "display_label": tissue.get("display_label"),
                    "ontology_id": tissue.get("ontology_id"),
                    "sample_count": tissue.get("sample_count"),
                    "median_tpm": tissue.get("median_tpm"),
                    "q1_tpm": tissue.get("q1_tpm"),
                    "q3_tpm": tissue.get("q3_tpm"),
                    "min_tpm": tissue.get("min_tpm"),
                    "max_tpm": tissue.get("max_tpm"),
                    "unit": "TPM",
                    "color_hex": tissue.get("color_hex"),
                    "dataset_id": section_cfg.gtex_dataset_id,
                    "gencode_id": gencode_id,
                },
                display_text=(
                    f"{gene} GTEx median TPM in {tissue.get('display_label') or tissue_id} "
                    f"is {tissue.get('median_tpm')}."
                ),
                api_run_id=parent_api,
                raw_artifact_id=parent_raw,
                confidence_notes="GTEx is human-only. Sample arrays live in raw artifacts.",
            )
            _append_evidence(evidence, rec, persist_db=persist_db)

        collection = _evidence(
            dossier_run_id=run_id,
            gene_symbol=gene,
            source_name=gtex.SOURCE_NAME,
            fact_type="gtex_expression_collection_summary",
            key=f"{gencode_id}-collection",
            value={
                "total_tissue_count": len(ordered_all),
                "total_sample_value_count": sum(
                    int(t.get("sample_count") or 0) for t in ordered_all
                ),
                "brain_tissue_count": len(ordered_brain),
                "median_validation": {
                    k: median_validation.get(k)
                    for k in (
                        "matched_tissue_count",
                        "mismatched_tissue_count",
                        "missing_median_count",
                        "maximum_absolute_difference",
                        "api_median_unavailable",
                        "systematic_mismatch",
                    )
                    if k in median_validation
                    or median_validation.get("api_median_unavailable")
                },
                "dataset_id": section_cfg.gtex_dataset_id,
                "gencode_id": gencode_id,
                "unit": "TPM",
                "color_status": gtex_status["color_status"],
                "presentation_item_key": f"gtex-{(gene or '').lower()}",
            },
            display_text=(
                f"{gene} has GTEx expression in {len(ordered_all)} tissues "
                f"({len(ordered_brain)} brain)."
            ),
            api_run_id=parent_api,
            raw_artifact_id=parent_raw,
            confidence_notes="GTEx is human-only.",
        )
        _append_evidence(evidence, collection, persist_db=persist_db)

        # Figures from sample arrays only
        if ordered_all and not median_validation.get("systematic_mismatch"):
            try:
                all_png = render_gtex_violin_png(
                    ordered_all,
                    gene_symbol=gene,
                    gencode_id=gencode_id or "",
                    dpi=section_cfg.plot_dpi,
                    figsize=(11.2, 5.4),
                )
                all_api = _tool_result_to_api_run(
                    ToolResult(
                        source_name=gtex.SOURCE_NAME,
                        endpoint_name="gtex_all_tissues_violin",
                        success=True,
                        gene_symbol=gene,
                        request_url=gtex_gene_url(gene),
                        request_params={
                            "derivation_type": "gtex_violin_plot",
                            "subset": "all_tissues",
                        },
                        data={"tissue_count": len(ordered_all)},
                    ),
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                )
                all_artifact, all_meta = _persist_artifact_bytes(
                    dossier_run_id=run_id,
                    source_name=gtex.SOURCE_NAME,
                    content=all_png,
                    extension="png",
                    artifact_type="png",
                    filename_hint=f"gtex-all-tissues-violin-{slugify(gene)}",
                    settings=cfg,
                    api_run=all_api,
                    persist_db=persist_db,
                    notes={
                        "artifact_class": "derived",
                        "derivation_type": "gtex_violin_plot",
                        "subset": "all_tissues",
                        "plot_version": PLOT_VERSION,
                        "parent_raw_artifact_ids": sample_parent_raw_ids,
                    },
                    validate=_validate_nonblank_image,
                )
                api_runs.append(all_api)
                raw_meta.append(all_meta)
                fig_rec = _evidence(
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    source_name=gtex.SOURCE_NAME,
                    fact_type="gtex_all_tissues_figure",
                    key=f"{gencode_id}-all-tissues-figure",
                    value={
                        "artifact_class": "derived",
                        "derivation_type": "gtex_violin_plot",
                        "subset": "all_tissues",
                        "relative_path": all_meta.get("relative_path"),
                        "local_artifact_path": all_meta.get("relative_path"),
                        "media_type": "image/png",
                        "width": all_meta.get("width"),
                        "height": all_meta.get("height"),
                        "byte_size": all_meta.get("byte_size"),
                        "sha256": all_meta.get("content_hash") or all_artifact.content_hash,
                        "dataset_id": section_cfg.gtex_dataset_id,
                        "gencode_id": gencode_id,
                        "tissue_count": len(ordered_all),
                        "sample_count": sum(
                            int(t.get("sample_count") or 0) for t in ordered_all
                        ),
                        "parent_raw_artifact_ids": sample_parent_raw_ids,
                        "plot_version": PLOT_VERSION,
                        "presentation_item_key": f"gtex-{(gene or '').lower()}",
                    },
                    display_text=f"{gene} GTEx all-tissue expression violin plot.",
                    api_run_id=all_api.id,
                    raw_artifact_id=all_artifact.id,
                )
                _append_evidence(evidence, fig_rec, persist_db=persist_db)
                gtex_status["all_tissues_figure"] = "success"
                audit["gtex"]["all_tissues_figure"] = {
                    "relative_path": all_meta.get("relative_path"),
                    "sha256": all_meta.get("content_hash") or all_artifact.content_hash,
                    "width": all_meta.get("width"),
                    "height": all_meta.get("height"),
                    "byte_size": all_meta.get("byte_size"),
                }
            except Exception as exc:  # noqa: BLE001
                logger.exception("GTEx all-tissue violin failed for %s", gene)
                gtex_status["all_tissues_figure"] = "failed"
                errors.append(f"GTEx all-tissue figure failed: {exc}")
                audit["gtex"]["all_tissues_figure"] = {
                    "status": "failed",
                    "error": str(exc),
                }

        if ordered_brain and not median_validation.get("systematic_mismatch"):
            try:
                brain_png = render_gtex_violin_png(
                    ordered_brain,
                    gene_symbol=gene,
                    gencode_id=gencode_id or "",
                    dpi=section_cfg.plot_dpi,
                    figsize=(9.5, 5.0),
                )
                brain_api = _tool_result_to_api_run(
                    ToolResult(
                        source_name=gtex.SOURCE_NAME,
                        endpoint_name="gtex_brain_tissues_violin",
                        success=True,
                        gene_symbol=gene,
                        request_url=gtex_gene_url(gene),
                        request_params={
                            "derivation_type": "gtex_violin_plot",
                            "subset": "brain_tissues",
                        },
                        data={"tissue_count": len(ordered_brain)},
                    ),
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                )
                brain_artifact, brain_meta = _persist_artifact_bytes(
                    dossier_run_id=run_id,
                    source_name=gtex.SOURCE_NAME,
                    content=brain_png,
                    extension="png",
                    artifact_type="png",
                    filename_hint=f"gtex-brain-tissues-violin-{slugify(gene)}",
                    settings=cfg,
                    api_run=brain_api,
                    persist_db=persist_db,
                    notes={
                        "artifact_class": "derived",
                        "derivation_type": "gtex_violin_plot",
                        "subset": "brain_tissues",
                        "plot_version": PLOT_VERSION,
                        "parent_raw_artifact_ids": sample_parent_raw_ids,
                    },
                    validate=_validate_nonblank_image,
                )
                api_runs.append(brain_api)
                raw_meta.append(brain_meta)
                fig_rec = _evidence(
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    source_name=gtex.SOURCE_NAME,
                    fact_type="gtex_brain_tissues_figure",
                    key=f"{gencode_id}-brain-tissues-figure",
                    value={
                        "artifact_class": "derived",
                        "derivation_type": "gtex_violin_plot",
                        "subset": "brain_tissues",
                        "relative_path": brain_meta.get("relative_path"),
                        "local_artifact_path": brain_meta.get("relative_path"),
                        "media_type": "image/png",
                        "width": brain_meta.get("width"),
                        "height": brain_meta.get("height"),
                        "byte_size": brain_meta.get("byte_size"),
                        "sha256": brain_meta.get("content_hash")
                        or brain_artifact.content_hash,
                        "dataset_id": section_cfg.gtex_dataset_id,
                        "gencode_id": gencode_id,
                        "tissue_count": len(ordered_brain),
                        "sample_count": sum(
                            int(t.get("sample_count") or 0) for t in ordered_brain
                        ),
                        "parent_raw_artifact_ids": sample_parent_raw_ids,
                        "plot_version": PLOT_VERSION,
                        "presentation_item_key": f"gtex-{(gene or '').lower()}",
                    },
                    display_text=f"{gene} GTEx brain-tissue expression violin plot.",
                    api_run_id=brain_api.id,
                    raw_artifact_id=brain_artifact.id,
                )
                _append_evidence(evidence, fig_rec, persist_db=persist_db)
                gtex_status["brain_figure"] = "success"
                audit["gtex"]["brain_figure"] = {
                    "relative_path": brain_meta.get("relative_path"),
                    "sha256": brain_meta.get("content_hash") or brain_artifact.content_hash,
                    "width": brain_meta.get("width"),
                    "height": brain_meta.get("height"),
                    "byte_size": brain_meta.get("byte_size"),
                }
            except Exception as exc:  # noqa: BLE001
                logger.exception("GTEx brain violin failed for %s", gene)
                gtex_status["brain_figure"] = "failed"
                errors.append(f"GTEx brain figure failed: {exc}")
                audit["gtex"]["brain_figure"] = {"status": "failed", "error": str(exc)}
        elif gtex_status.get("sample_expression") != "success":
            gtex_status["all_tissues_figure"] = "unavailable"
            gtex_status["brain_figure"] = "unavailable"

    # --- HBT (attempt even if GTEx failed) ---
    hbt_tr = hbt.fetch_whole_brain_pdf(gene, settings=cfg)
    # Persist PDF bytes separately from ToolResult JSON (bytes not JSON-serializable)
    hbt_pdf_bytes = None
    hbt_data = dict(hbt_tr.data) if isinstance(hbt_tr.data, dict) else {}
    if hbt_tr.success:
        hbt_pdf_bytes = hbt_data.pop("pdf_bytes", None)
        hbt_tr = ToolResult(
            source_name=hbt_tr.source_name,
            endpoint_name=hbt_tr.endpoint_name,
            success=True,
            gene_symbol=hbt_tr.gene_symbol,
            request_url=hbt_tr.request_url,
            request_params=hbt_tr.request_params,
            status_code=hbt_tr.status_code,
            data={k: v for k, v in hbt_data.items() if k != "pdf_bytes"},
        )
    hbt_api = _tool_result_to_api_run(hbt_tr, dossier_run_id=run_id, gene_symbol=gene)
    api_runs.append(hbt_api)
    audit["hbt"]["fetch"] = {
        "success": hbt_tr.success,
        "error_type": hbt_tr.error_type,
        "error_message": hbt_tr.error_message,
        "request_url": hbt_tr.request_url,
        "api_run_id": hbt_api.id,
    }
    if not hbt_tr.success:
        _save_api_run_failure(hbt_api, persist_db=persist_db)

    if hbt_tr.success and isinstance(hbt_pdf_bytes, (bytes, bytearray)):
        try:
            pdf_artifact, pdf_meta = _persist_artifact_bytes(
                dossier_run_id=run_id,
                source_name=hbt.SOURCE_NAME,
                content=bytes(hbt_pdf_bytes),
                extension="pdf",
                artifact_type="pdf",
                filename_hint=f"hbt-whole-brain-{slugify(gene)}",
                settings=cfg,
                api_run=hbt_api,
                persist_db=persist_db,
                notes={
                    "artifact_class": "external_raw",
                    "artifact_origin": "hbt_whole_brain_pdf",
                    "source_url": hbt_pdf_url(gene),
                    "retrieval_method": "http_get",
                },
                validate=_validate_pdf,
            )
            raw_meta.append(pdf_meta)
            hbt_status["pdf"] = "success"
            pdf_rec = _evidence(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name=hbt.SOURCE_NAME,
                fact_type="hbt_whole_brain_pdf",
                key=f"{gene}-hbt-pdf",
                value={
                    "source_url": hbt_pdf_url(gene),
                    "relative_path": pdf_meta.get("relative_path"),
                    "media_type": "application/pdf",
                    "byte_size": pdf_meta.get("byte_size"),
                    "sha256": pdf_meta.get("content_hash") or pdf_artifact.content_hash,
                    "page_count": hbt_data.get("page_count"),
                    "selected_page_index": hbt_data.get("selected_page_index"),
                    "gene_text_found": hbt_data.get("gene_text_found"),
                },
                display_text=f"{gene} Human Brain Transcriptome whole-brain PDF.",
                api_run_id=hbt_api.id,
                raw_artifact_id=pdf_artifact.id,
            )
            _append_evidence(evidence, pdf_rec, persist_db=persist_db)
            audit["hbt"]["pdf"] = {
                "relative_path": pdf_meta.get("relative_path"),
                "sha256": pdf_meta.get("content_hash") or pdf_artifact.content_hash,
                "byte_size": pdf_meta.get("byte_size"),
                "source_url": hbt_pdf_url(gene),
            }

            page_index = int(hbt_data.get("selected_page_index") or 0)
            raster = hbt.rasterize_pdf_page(
                bytes(hbt_pdf_bytes),
                page_index=page_index,
                dpi=section_cfg.hbt_raster_dpi,
            )
            if raster is None:
                hbt_status["figure"] = "failed"
                audit["hbt"]["figure"] = {
                    "status": "failed",
                    "error": "pdf_rasterization_failed",
                }
            else:
                png_bytes, raster_meta = raster
                fig_api = _tool_result_to_api_run(
                    ToolResult(
                        source_name=hbt.SOURCE_NAME,
                        endpoint_name="hbt_whole_brain_raster",
                        success=True,
                        gene_symbol=gene,
                        request_url=hbt_pdf_url(gene),
                        request_params={
                            "derivation_type": "pdf_page_rasterization",
                            "page_index": page_index,
                            "dpi": section_cfg.hbt_raster_dpi,
                        },
                        data=raster_meta,
                    ),
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                )
                fig_artifact, fig_meta = _persist_artifact_bytes(
                    dossier_run_id=run_id,
                    source_name=hbt.SOURCE_NAME,
                    content=png_bytes,
                    extension="png",
                    artifact_type="png",
                    filename_hint=f"hbt-whole-brain-figure-{slugify(gene)}",
                    settings=cfg,
                    api_run=fig_api,
                    persist_db=persist_db,
                    notes={
                        "artifact_class": "official_source_render",
                        "derivation_type": "pdf_page_rasterization",
                        "source": "Human Brain Transcriptome",
                        "source_pdf_url": hbt_pdf_url(gene),
                        "source_pdf_sha256": pdf_artifact.content_hash,
                        "page_index": page_index,
                    },
                    validate=_validate_nonblank_image,
                )
                api_runs.append(fig_api)
                raw_meta.append(fig_meta)
                fig_rec = _evidence(
                    dossier_run_id=run_id,
                    gene_symbol=gene,
                    source_name=hbt.SOURCE_NAME,
                    fact_type="hbt_whole_brain_figure",
                    key=f"{gene}-hbt-figure",
                    value={
                        "artifact_class": "official_source_render",
                        "derivation_type": "pdf_page_rasterization",
                        "source": "Human Brain Transcriptome",
                        "source_url": hbt_pdf_url(gene),
                        "source_pdf_sha256": pdf_artifact.content_hash,
                        "relative_path": fig_meta.get("relative_path"),
                        "local_artifact_path": fig_meta.get("relative_path"),
                        "media_type": "image/png",
                        "width": fig_meta.get("width") or raster_meta.get("width"),
                        "height": fig_meta.get("height") or raster_meta.get("height"),
                        "byte_size": fig_meta.get("byte_size"),
                        "sha256": fig_meta.get("content_hash") or fig_artifact.content_hash,
                        "page_index": page_index,
                        "dpi": section_cfg.hbt_raster_dpi,
                        "presentation_item_key": f"hbt-{(gene or '').lower()}",
                    },
                    display_text=(
                        f"{gene} Human Brain Transcriptome developmental-expression figure."
                    ),
                    api_run_id=fig_api.id,
                    raw_artifact_id=fig_artifact.id,
                )
                _append_evidence(evidence, fig_rec, persist_db=persist_db)
                hbt_status["figure"] = "success"
                audit["hbt"]["figure"] = {
                    "relative_path": fig_meta.get("relative_path"),
                    "sha256": fig_meta.get("content_hash") or fig_artifact.content_hash,
                    "width": fig_meta.get("width") or raster_meta.get("width"),
                    "height": fig_meta.get("height") or raster_meta.get("height"),
                    "byte_size": fig_meta.get("byte_size"),
                    "page_index": page_index,
                }
        except Exception as exc:  # noqa: BLE001
            logger.exception("HBT persistence/raster failed for %s", gene)
            hbt_status["pdf"] = "failed"
            hbt_status["figure"] = "failed"
            errors.append(f"HBT figure failed: {exc}")
            audit["hbt"]["error"] = str(exc)
    else:
        hbt_status["pdf"] = "unavailable"
        hbt_status["figure"] = "unavailable"
        coverage_extra.append(
            SourceCoverageResult(
                dossier_run_id=run_id,
                source_name=hbt.SOURCE_NAME,
                status=SourceStatus.failed if not hbt_tr.success else SourceStatus.skipped,
                evidence_record_count=0,
                error_message=hbt_tr.error_message or "HBT PDF unavailable",
                report_sections_supported=[SUBSECTION_2A],
            )
        )

    overall = "complete"
    if (
        gtex_status.get("all_tissues_figure") != "success"
        and hbt_status.get("figure") != "success"
    ):
        overall = "empty" if gtex_status.get("resolve") == "failed" else "partial"
    elif (
        gtex_status.get("all_tissues_figure") != "success"
        or hbt_status.get("figure") != "success"
    ):
        overall = "partial"

    section_status = {
        "rendering_status": {
            "overall": overall,
            **{f"gtex_{k}": v for k, v in gtex_status.items()},
            **{f"hbt_{k}": v for k, v in hbt_status.items()},
        },
        "summary": {
            "gene_symbol": gene,
            "gencode_id": gencode_id,
            "dataset_id": section_cfg.gtex_dataset_id,
            "total_tissue_count": len(ordered_all),
            "brain_tissue_count": len(ordered_brain),
            "median_validation": median_validation,
            "gtex_portal_url": gtex_gene_url(gene),
            "hbt_pdf_url": hbt_pdf_url(gene),
            "hbt_home_url": HBT_HOME_URL,
            "color_status": gtex_status.get("color_status"),
            "presentation_item_key": f"tissue-{(gene or '').lower()}",
        },
        "audit": audit,
    }

    return {
        **state,
        "evidence_records": evidence,
        "api_runs": api_runs,
        "raw_artifacts": raw_meta,
        "errors": errors,
        "coverage": coverage_extra,
        "section_2a_status": section_status,
    }


__all__ = [
    "SECTION_EXPRESSION",
    "SUBSECTION_2A",
    "EXPECTED_GTEX_V8_BRAIN_TISSUE_COUNT",
    "EXPECTED_GTEX_V8_BRAIN_TISSUES",
    "Section2aConfig",
    "gtex_intro_text",
    "hbt_intro_text",
    "hbt_link_text",
    "gtex_gene_url",
    "hbt_pdf_url",
    "tissue_stats",
    "parse_sample_expression_rows",
    "parse_median_rows",
    "parse_tissue_metadata",
    "validate_medians",
    "order_tissues",
    "brain_subset",
    "enrich_tissue_display",
    "render_gtex_violin_png",
    "node_generate_section_2a_derived_artifacts",
]
