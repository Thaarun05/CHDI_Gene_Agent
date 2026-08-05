"""Bundle-only Section 2c single-cell / single-nucleus cell-type expression.

Three independent evidence branches, none of which may abort another:

- Allen Human M1 10x trimmed means + dendrogram (single-nucleus RNA sequencing).
- Allen mouse whole cortex and hippocampus 10x trimmed means + dendrogram.
- DropViz population-level expression derived from the published GSE116470
  metacell matrix (single-cell RNA sequencing using Drop-seq).

Dataset-level source files are shared by every gene, so they are read from the
protected accepted-source layout in :mod:`gene_dossier.section_2c_sources` and
re-acquired only under an explicit ``force_refresh``. Production Section 2c never
calls the live DropViz Shiny client, so ``dropviz_live_status`` defaults to
``not_attempted_optional``.

Every quantitative claim is scoped to the dataset that produced it, a value
greater than zero is reported as *nonzero aggregate expression* rather than a
biological detection, and the DropViz ranking is charted only when the matrix
value semantics establish an interpretable metric and unit.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import (
    ApiRun,
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
    _save_api_run_failure,
    _tool_result_to_api_run,
    _validate_nonblank_image,
)
from gene_dossier.section_2b import _resolve_mouse_symbol
from gene_dossier.section_2c_sources import (
    Section2cPaths,
    accept_gene_report,
    load_accepted_source,
    paths_for,
    write_json_atomic,
)
from gene_dossier.source_ids import make_source_id, slugify
from gene_dossier.tools import allen_celltypes as ac
from gene_dossier.tools import dropviz_geo as dg
from gene_dossier.workflow import DossierState, WorkflowTransientContext

logger = logging.getLogger(__name__)

SECTION_EXPRESSION = "Tissue and cell expression"
SUBSECTION_2C = "snRNA-Seq gene expression in cell type database"

CALCULATION_VERSION = "section_2c_v1"
CHART_VERSION = "dropviz_top_populations_v1"

GEO_SERIES_URL = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dg.GEO_ACCESSION}"

# Status vocabulary. Kept as constants so the section node, the presentation
# layer and the tests all agree on the exact strings.
STATUS_SUCCESS = "success"
STATUS_SOURCE_UNAVAILABLE = "source_unavailable"
STATUS_NOT_ATTEMPTED = "not_attempted_optional"
STATUS_NORMALIZATION_UNRESOLVED = dg.RANK_STATUS_NORMALIZATION_UNRESOLVED
STATUS_FIGURE_SUPPRESSED = "suppressed_normalization_unresolved"
STATUS_NOT_AVAILABLE = "not_available"
STATUS_PARTIAL = "partial"
STATUS_UNAVAILABLE = "unavailable"
STATUS_NOT_IN_SCOPE = "not_in_production_scope"

CONFIDENCE_INTERVAL_STATUS = dg.CONFIDENCE_INTERVAL_METHOD_UNRESOLVED
DROPVIZ_REGIONAL_TSNE_STATUS = STATUS_NOT_IN_SCOPE
DROPVIZ_DISPLAY_LABEL_STATUS = dg.LABEL_MAPPING_PARTIAL
DROPVIZ_HISTORICAL_APP_URL = "http://dropviz.org/"
FALLBACK_TOP_N = 5

SECTION_2C_VISUAL_COMPLETE_ROLES = frozenset(
    {
        "section_2c_human_scatter_figure",
        "section_2c_human_heatmap_figure",
        "section_2c_mouse_scatter_figure",
        "section_2c_mouse_heatmap_figure",
        "section_2c_dropviz_rank_figure",
    }
)

# DropViz assay wording is Drop-seq specific: it is not a single-nucleus assay,
# so the shared Allen terminology map cannot describe it.
DROPVIZ_ASSAY_TERMINOLOGY = "single-cell RNA sequencing using Drop-seq"
DROPVIZ_SAMPLING_SCOPE = (
    "broader adult mouse-brain population context, subject to the regions and "
    f"processing represented by {dg.GEO_ACCESSION}"
)

SECTION_2C_INTRO_TEXT = (
    "Single-cell and single-nucleus transcriptomic datasets provide complementary "
    "views of gene expression across brain cell populations. The Allen datasets "
    "characterize transcriptomic cell types in human motor cortex and mouse "
    "cortex/hippocampus, while DropViz provides adult mouse-brain population-level "
    "expression derived from Drop-seq data. Absolute expression values should be "
    "interpreted within each dataset and should not be compared directly across "
    "Human M1, Mouse CTX-HPF, and DropViz because the assays, aggregation "
    "procedures, and expression scales differ."
)

DROPVIZ_POPULATION_IDENTIFIER_NOTE = (
    "Population identifiers encode brain-region and cluster numbers; descriptive "
    "cell-class mapping was unavailable."
)

DROPVIZ_REGIONAL_TSNE_LIMITATION_NOTE = (
    "Historical DropViz regional t-SNE views were not included because the saved "
    "application states are no longer reproducible; production results use the "
    "archived GSE116470 matrix."
)

TRIMMED_MEAN_VALUE_LABEL = "Trimmed mean expression"
HUMAN_LABEL_COLUMN = "Human M1 transcriptomic cell type"
MOUSE_LABEL_COLUMN = "Mouse cortex/hippocampus cluster"
DROPVIZ_LABEL_COLUMN = "DropViz population"
DROPVIZ_FALLBACK_VALUE_LABEL = "Population-level expression"

THERAPEUTIC_CAVEAT = (
    "Aggregate trimmed-mean and population-level expression values report RNA "
    "abundance in the populations each dataset sampled; they do not establish "
    "protein abundance, functional dependence in any cell type, or therapeutic "
    "tractability, and the regions, species and processing represented by each "
    "dataset bound every statement above. These observations are cell-type "
    "localization context, not a target assessment."
)

# Cache keys for the five protected dataset-level sources.
_HUMAN_SOURCE_KEYS = (
    (ac.CACHE_KEY_HUMAN_TRIMMED_MEANS, "trimmed_means"),
    (ac.CACHE_KEY_HUMAN_TAXONOMY, "taxonomy"),
)
_MOUSE_SOURCE_KEYS = (
    (ac.CACHE_KEY_MOUSE_TRIMMED_MEANS, "trimmed_means"),
    (ac.CACHE_KEY_MOUSE_TAXONOMY, "taxonomy"),
)
GEO_SOURCE_KEY = "gse116470_metacells"

_FIGURE_VISUALIZATIONS = (ac.VISUALIZATION_SCATTER, ac.VISUALIZATION_HEATMAP)


@dataclass(frozen=True)
class Section2cConfig:
    """Section 2c policy knobs. Every default keeps the section offline-safe."""

    output_root: str | None = None
    force_refresh: bool = False
    top_n: int = dg.TOP_N_DEFAULT
    plot_dpi: int = 180
    attempt_allen_figures: bool = True
    figure_max_attempts: int = 2
    # Production Section 2c never reads data from the live DropViz Shiny server;
    # the live client stays an opt-in diagnostic.
    attempt_live_dropviz: bool = False
    documented_semantics: str | None = None
    documented_unit: str | None = None
    documentation_reference: str | None = None

    def __post_init__(self) -> None:
        if int(self.top_n) < 1:
            raise ValueError("top_n must be >= 1")
        if int(self.plot_dpi) < 72:
            raise ValueError("plot_dpi must be >= 72")
        if int(self.figure_max_attempts) < 1:
            raise ValueError("figure_max_attempts must be >= 1")
        allowed = {
            dg.VALUE_SEMANTICS_RAW_COUNTS,
            dg.VALUE_SEMANTICS_COUNT_COMPATIBLE,
            dg.VALUE_SEMANTICS_NORMALIZED,
            dg.VALUE_SEMANTICS_TRANSFORMED,
            dg.VALUE_SEMANTICS_UNRESOLVED,
        }
        if self.documented_semantics is not None:
            if self.documented_semantics not in allowed:
                raise ValueError(
                    f"Unsupported documented_semantics: {self.documented_semantics!r}"
                )
            if not (self.documentation_reference or "").strip():
                raise ValueError(
                    "documented_semantics requires a documentation_reference citation"
                )


# --------------------------------------------------------------------------------------
# Formatting helpers
# --------------------------------------------------------------------------------------
def format_value(value: Any, *, digits: int = 4) -> str:
    """Format one numeric cell for display, or return an explicit unavailability."""
    if value is None or isinstance(value, bool):
        return STATUS_NOT_AVAILABLE.replace("_", " ")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return STATUS_NOT_AVAILABLE.replace("_", " ")
    if not math.isfinite(number):
        return STATUS_NOT_AVAILABLE.replace("_", " ")
    text = f"{number:.{digits}f}".rstrip("0").rstrip(".")
    return text or "0"


def _pct(value: Any) -> str:
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return STATUS_NOT_AVAILABLE.replace("_", " ")


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
    organism: str = "Homo sapiens",
    taxon_id: int = 9606,
    confidence_notes: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=make_source_id(source_name, gene_symbol, AssertionType.expression, key),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_EXPRESSION,
        subsection=SUBSECTION_2C,
        source_name=source_name,
        source_type=SourceType.expression_database,
        assertion_type=AssertionType.expression,
        fact_type=fact_type,
        organism=organism,
        taxon_id=taxon_id,
        evidence_grade=EvidenceGrade.B,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
        confidence_notes=confidence_notes,
    )


def _validate_bytes(content: bytes, *, media_type: str) -> dict[str, Any]:
    if not content:
        raise ValueError("derived artifact is empty")
    return {"media_type": media_type, "byte_size": len(content)}


def _persist_derived(
    *,
    dossier_run_id: str,
    source_name: str,
    content: bytes,
    extension: str,
    media_type: str,
    filename_hint: str,
    artifact_role: str,
    parent_raw_artifact_ids: Sequence[str],
    settings: Settings,
    persist_db: bool,
    raw_meta: list,
    extra_notes: dict[str, Any] | None = None,
    validate: Any = None,
) -> tuple[Any, dict[str, Any]] | tuple[None, None]:
    """Persist one locally derived artifact with ``api_run=None`` and parent links.

    A local derivation is never given a fabricated ApiRun; provenance flows
    through ``parent_raw_artifact_ids`` instead.
    """
    try:
        artifact, meta = _persist_artifact_bytes(
            dossier_run_id=dossier_run_id,
            source_name=source_name,
            content=content,
            extension=extension,
            artifact_type=extension,
            filename_hint=filename_hint,
            settings=settings,
            api_run=None,
            persist_db=persist_db,
            notes={
                "artifact_class": "derived",
                "artifact_origin": "section_2c",
                "artifact_role": artifact_role,
                "retrieval_method": "local_derivation",
                "calculation_version": CALCULATION_VERSION,
                "parent_raw_artifact_ids": [str(i) for i in parent_raw_artifact_ids if i],
                **(extra_notes or {}),
            },
            validate=validate or (lambda b: _validate_bytes(b, media_type=media_type)),
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning("Section 2c derived artifact %s failed: %s", artifact_role, exc)
        return None, None
    raw_meta.append(meta)
    return artifact, meta


def _json_bytes(payload: Any) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False, default=str) + "\n"
    ).encode("utf-8")


def _csv_bytes(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> bytes:
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(list(headers))
    for row in rows:
        writer.writerow(["" if cell is None else cell for cell in row])
    return buf.getvalue().encode("utf-8")


# --------------------------------------------------------------------------------------
# Protected dataset-level source access
# --------------------------------------------------------------------------------------
@dataclass
class _SourcePayload:
    """One dataset-level source resolved for this run."""

    ok: bool
    source_key: str
    content: bytes | None = None
    sha256: str | None = None
    byte_size: int | None = None
    official_url: str | None = None
    origin: str = "accepted_pointer"
    raw_artifact_id: str | None = None
    error_type: str | None = None
    error_message: str | None = None

    def audit(self) -> dict[str, Any]:
        return {
            "source_key": self.source_key,
            "resolved": self.ok,
            "origin": self.origin,
            "sha256": self.sha256,
            "byte_size": self.byte_size,
            "official_url": self.official_url,
            "raw_artifact_id": self.raw_artifact_id,
            "error_type": self.error_type,
            "error_message": self.error_message,
        }


def _payload_from_accepted(record: dict[str, Any], source_key: str) -> _SourcePayload:
    artifact = Path(str(record.get("artifact_path") or ""))
    try:
        content = artifact.read_bytes()
    except OSError as exc:
        return _SourcePayload(
            ok=False,
            source_key=source_key,
            error_type="accepted_artifact_unreadable",
            error_message=str(exc),
        )
    return _SourcePayload(
        ok=True,
        source_key=source_key,
        content=content,
        sha256=str(record.get("sha256") or ""),
        byte_size=int(record.get("byte_size") or len(content)),
        official_url=str(record.get("official_url") or "") or None,
        origin="accepted_pointer",
    )


def _resolve_dataset_source(
    *,
    paths: Any,
    source_key: str,
    official_url: str | None,
    force_refresh: bool,
    downloader: Any | None,
    dossier_run_id: str,
    gene_symbol: str,
    settings: Settings,
    persist_db: bool,
    api_runs: list,
    raw_meta: list,
) -> _SourcePayload:
    """Reuse the accepted dataset source; re-download only under ``force_refresh``.

    Locally registered sources have no downloader, so a missing pointer for those
    is reported honestly rather than reconstructed.
    """
    if not force_refresh:
        record = load_accepted_source(
            paths, source_key=source_key, official_url=official_url
        )
        if record:
            return _payload_from_accepted(record, source_key)

    if downloader is None:
        return _SourcePayload(
            ok=False,
            source_key=source_key,
            official_url=official_url,
            origin="accepted_pointer_missing",
            error_type="accepted_source_missing",
            error_message=(
                f"no accepted dataset source for {source_key!r}; this source is "
                "registered locally and cannot be re-downloaded automatically"
            ),
        )

    tr: ToolResult = downloader()
    api = _tool_result_to_api_run(
        tr, dossier_run_id=dossier_run_id, gene_symbol=gene_symbol
    )
    data = tr.data if isinstance(tr.data, dict) else {}
    content = data.get("content") or b""
    if not tr.success or not content:
        _save_api_run_failure(api, persist_db=persist_db)
        api_runs.append(api)
        return _SourcePayload(
            ok=False,
            source_key=source_key,
            official_url=official_url,
            origin="live_download",
            error_type=tr.error_type or "download_failed",
            error_message=tr.error_message,
        )

    api_runs.append(api)
    meta: dict[str, Any] | None = None
    try:
        _artifact, meta = _persist_artifact_bytes(
            dossier_run_id=dossier_run_id,
            source_name=tr.source_name,
            content=content,
            extension="bin",
            artifact_type="bin",
            filename_hint=f"section-2c-{slugify(source_key)}",
            settings=settings,
            api_run=api,
            persist_db=persist_db,
            notes={
                "artifact_class": "external_raw",
                "artifact_origin": f"section_2c_{source_key}",
                "artifact_role": source_key,
                "source_url": tr.request_url,
                "retrieval_method": "http_bytes",
                "sha256": data.get("sha256"),
            },
            validate=lambda b: _validate_bytes(b, media_type="application/octet-stream"),
        )
        raw_meta.append(meta)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Section 2c refresh persist failed for %s: %s", source_key, exc)

    return _SourcePayload(
        ok=True,
        source_key=source_key,
        content=content,
        sha256=str(data.get("sha256") or "") or None,
        byte_size=len(content),
        official_url=official_url,
        origin="live_download",
        raw_artifact_id=(meta or {}).get("id"),
    )


# --------------------------------------------------------------------------------------
# Allen branch
# --------------------------------------------------------------------------------------
def analyze_allen_dataset(
    *,
    dataset: str,
    gene_symbol: str,
    trimmed_means: bytes,
    taxonomy_json: bytes,
    count_noun: str,
    top_n: int = dg.TOP_N_DEFAULT,
) -> dict[str, Any]:
    """Structured trimmed-means analysis for one Allen dataset (never raises)."""
    out: dict[str, Any] = {
        "ok": False,
        "dataset": dataset,
        "dataset_label": ac.DATASET_LABELS.get(dataset, dataset),
        "assay_terminology": ac.DATASET_ASSAY_TERMS.get(dataset),
        "sampling_scope": ac.DATASET_SAMPLING_SCOPE.get(dataset),
        "count_noun": count_noun,
        "error_type": None,
        "error_message": None,
    }
    try:
        extraction = ac.extract_gene_row(ac.text_lines(trimmed_means), gene_symbol)
    except Exception as exc:  # noqa: BLE001
        out["error_type"] = "trimmed_means_unreadable"
        out["error_message"] = str(exc)
        return out

    out["match_count"] = extraction.match_count
    out["gene_row_count"] = extraction.gene_row_count
    out["malformed_row_count"] = extraction.malformed_row_count
    if not extraction.ok:
        out["error_type"] = extraction.error_type
        out["error_message"] = extraction.error_message
        return out

    try:
        taxonomy = ac.parse_dendrogram(json.loads(taxonomy_json.decode("utf-8")))
    except Exception as exc:  # noqa: BLE001
        out["error_type"] = "taxonomy_unreadable"
        out["error_message"] = str(exc)
        return out

    reconciliation = ac.reconcile_taxonomy(
        taxonomy_leaves=taxonomy.leaf_order,
        matrix_labels=extraction.celltype_labels,
    )
    summary = ac.summarize_celltype_expression(
        gene_symbol=gene_symbol,
        source_symbol=extraction.source_symbol,
        celltype_labels=extraction.celltype_labels,
        values=extraction.values,
        taxonomy=taxonomy,
        count_noun=count_noun,
        top_n=top_n,
    )
    out.update(
        {
            "ok": True,
            "summary": summary,
            "taxonomy_reconciliation": reconciliation,
            "taxonomy_leaf_count": taxonomy.leaf_count,
            "named_internal_node_count": taxonomy.internal_label_count,
            "source_symbol": extraction.source_symbol,
        }
    )
    return out


def _format_list_with_and(items: Sequence[str]) -> str:
    """Join items with commas and 'and' before the final entry."""
    values = [str(item) for item in items if str(item).strip()]
    if not values:
        return ""
    if len(values) == 1:
        return values[0]
    if len(values) == 2:
        return f"{values[0]} and {values[1]}"
    return ", ".join(values[:-1]) + f", and {values[-1]}"


def _display_taxonomy_ancestor(label: str) -> str:
    """Presentation-only ancestor label normalization."""
    text = str(label or "").strip()
    if text == "CGE (most)":
        return "CGE-derived"
    return text


def _overlapping_taxonomy_ancestors_sentence(
    ancestors: Sequence[dict[str, Any]],
    *,
    top_count: int,
    entity_noun: str,
) -> str | None:
    """Build overlapping-ancestor prose without implying mutual exclusivity."""
    if not ancestors or top_count <= 0:
        return None
    named = [
        f"{_display_taxonomy_ancestor(str(item.get('ancestor') or ''))}"
        f" ({item.get('top_member_count')})"
        for item in list(ancestors)[:4]
        if str(item.get("ancestor") or "").strip()
    ]
    if not named:
        return None
    return (
        f"Among the {top_count} highest-ranked {entity_noun}, the following "
        f"overlapping taxonomy ancestors were represented: "
        f"{_format_list_with_and(named)}."
    )


def allen_branch_narrative(analysis: dict[str, Any], *, gene_symbol: str) -> str:
    """Dataset-scoped narrative for one Allen branch, derived only from the data."""
    label = analysis.get("dataset_label") or "Allen dataset"
    assay = analysis.get("assay_terminology") or "single-cell RNA sequencing"
    scope = analysis.get("sampling_scope") or "the sampled populations"
    summary = dict(analysis.get("summary") or {})
    noun = str(analysis.get("count_noun") or "celltype")
    noun_plural = "transcriptomic cell types" if noun == "celltype" else "transcriptomic clusters"
    entity_noun = "cell types" if noun == "celltype" else "clusters"
    valid = summary.get(f"valid_{noun}_count")
    nonzero = summary.get(f"nonzero_{noun}_count")
    dataset_phrase = (
        "Within the Human M1 dataset"
        if analysis.get("dataset") == ac.DATASET_HUMAN_M1
        else f"Within the {label} dataset"
    )

    parts = [
        f"The Allen {label} dataset profiles {assay} across {noun_plural}. "
        f"{dataset_phrase}, the gene showed nonzero aggregate expression across "
        f"{nonzero} of {valid} {noun_plural} "
        f"({_pct(summary.get('nonzero_percentage'))}% of the "
        f"{noun_plural} carrying a valid trimmed-mean value)."
    ]

    top = list(summary.get("top") or [])
    if top:
        leader = top[0]
        parts.append(
            f"The highest aggregate trimmed mean was {format_value(leader.get('value'))} in "
            f"{leader.get('label')}, against a median of "
            f"{format_value(summary.get('median'))} across {noun_plural} "
            f"(a maximum-to-median ratio of "
            f"{format_value(summary.get('max_to_median_ratio'), digits=2)})."
        )
    ancestors = list(summary.get("top_taxonomy_ancestors") or [])
    ancestor_sentence = _overlapping_taxonomy_ancestors_sentence(
        ancestors, top_count=len(top), entity_noun=entity_noun
    )
    if ancestor_sentence:
        parts.append(ancestor_sentence)

    reconciliation = dict(analysis.get("taxonomy_reconciliation") or {})
    missing = list(reconciliation.get("missing_expression_clusters") or [])
    if missing:
        parts.append(
            f"The published taxonomy carries {reconciliation.get('taxonomy_leaf_count')} "
            f"leaves against {reconciliation.get('expression_cluster_count')} expression "
            f"columns; {', '.join(missing)} has no expression column in the matrix and "
            "was recorded as missing rather than filled with a zero."
            if len(missing) == 1
            else (
                f"The published taxonomy carries {reconciliation.get('taxonomy_leaf_count')} "
                f"leaves against {reconciliation.get('expression_cluster_count')} expression "
                f"columns; {', '.join(missing)} have no expression column in the matrix and "
                "were recorded as missing rather than filled with zeros."
            )
        )

    parts.append(
        f"A value greater than zero is reported here as nonzero aggregate expression "
        f"rather than a detection, because the trimmed-mean matrix carries no "
        f"source-defined detection threshold. These values describe localization within "
        f"{scope} and do not extend to regions or cell types the dataset does not sample."
    )
    _ = gene_symbol
    return " ".join(parts)


def allen_unavailable_note(*, dataset: str, analysis: dict[str, Any] | None) -> str:
    """Honest one-sentence unavailability note for one Allen branch."""
    label = ac.DATASET_LABELS.get(dataset, dataset)
    reason = (analysis or {}).get("error_type") or STATUS_SOURCE_UNAVAILABLE
    detail = (analysis or {}).get("error_message")
    text = (
        f"The Allen {label} structured expression analysis is unavailable for this run "
        f"({reason})."
    )
    if detail:
        text = f"{text} {detail}."
    return text


# --------------------------------------------------------------------------------------
# DropViz branch
# --------------------------------------------------------------------------------------
def analyze_dropviz_matrix(
    *,
    content: bytes,
    mouse_gene_symbol: str,
    source_sha256: str | None,
    source_url: str | None,
    top_n: int = dg.TOP_N_DEFAULT,
    documented_semantics: str | None = None,
    documented_unit: str | None = None,
    documentation_reference: str | None = None,
) -> dict[str, Any]:
    """Scan, classify, and rank the DropViz GEO matrix for one gene (never raises)."""
    out: dict[str, Any] = {"ok": False, "error_type": None, "error_message": None}
    try:
        scan = dg.scan_matrix(dg.open_matrix_stream(content), target_gene=mouse_gene_symbol)
    except Exception as exc:  # noqa: BLE001
        out["error_type"] = "matrix_unreadable"
        out["error_message"] = str(exc)
        return out
    if not scan.ok:
        out["error_type"] = scan.error_type
        out["error_message"] = scan.error_message
        return out

    semantics = dg.classify_value_semantics(
        scan,
        documented_semantics=documented_semantics,
        documentation_reference=documentation_reference,
        documented_unit=documented_unit,
    )
    profile = dg.build_matrix_profile(
        scan,
        semantics=semantics,
        source_sha256=source_sha256,
        source_url=source_url,
    )
    value_semantics_status = str(semantics["value_semantics_status"])
    ranking = dg.build_ranking_records(
        scan,
        value_semantics_status=value_semantics_status,
        documented_unit=semantics.get("documented_unit") or documented_unit,
    )
    parsed_labels = [dg.parse_population_label(label) for label in scan.population_labels]
    label_status = dg.label_mapping_status(parsed_labels)

    out.update(
        {
            "ok": True,
            "scan_ok": True,
            "gene_row_count": scan.gene_row_count,
            "population_column_count": len(scan.population_labels),
            "target_row_match_count": scan.target_matches,
            "target_source_symbol": scan.target_source_symbol,
            "value_semantics_status": value_semantics_status,
            "value_semantics_basis": semantics.get("basis"),
            "documented_unit": semantics.get("documented_unit"),
            "matrix_profile": profile,
            "rank_status": ranking["status"],
            "rank_reason": ranking.get("reason"),
            "presentation": ranking["presentation"],
            "records": ranking["records"],
            "excluded": ranking["excluded"],
            "parsed_labels": parsed_labels,
            "label_mapping_status": label_status,
            "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
        }
    )

    if ranking["status"] == dg.RANK_STATUS_SUCCESS:
        ranked = dg.rank_records(ranking["records"])
        out["ranked_records"] = ranked
        out["ranking_summary"] = dg.summarize_ranking(
            ranking["records"],
            value_semantics_status=value_semantics_status,
            top_n=top_n,
        )
    else:
        out["ranked_records"] = []
        out["ranking_summary"] = {}
    return out


def dropviz_narrative(analysis: dict[str, Any], *, mouse_gene_symbol: str) -> str:
    """Deterministic DropViz narrative scoped to the published GSE116470 matrix."""
    summary = dict(analysis.get("ranking_summary") or {})
    presentation = dict(analysis.get("presentation") or {})
    unit = presentation.get("expression_unit") or ""
    unit_phrase = f" {unit}" if unit else ""
    parts = [
        f"DropViz population-level expression for {mouse_gene_symbol} was derived from "
        f"the published {dg.GEO_ACCESSION} metacell matrix "
        f"({DROPVIZ_ASSAY_TERMINOLOGY}), which reports "
        f"{analysis.get('gene_row_count')} gene rows across "
        f"{analysis.get('population_column_count')} population columns.",
        f"Matrix values were classified as {analysis.get('value_semantics_status')} "
        f"({analysis.get('value_semantics_basis')}), so population expression is "
        f"reported as {presentation.get('ranking_metric') or 'the documented source metric'}.",
        f"{mouse_gene_symbol} showed nonzero expression in "
        f"{summary.get('nonzero_population_count')} of "
        f"{summary.get('valid_population_count')} valid populations "
        f"({_pct(summary.get('nonzero_percentage'))}%).",
    ]
    top = list(summary.get("top_populations") or [])
    if top:
        leader = top[0]
        parts.append(
            f"The highest-ranked population was {leader.get('population_label')} at "
            f"{format_value(leader.get('ranking_value'), digits=2)}{unit_phrase}, against "
            f"a median of {format_value(summary.get('median'), digits=2)} "
            f"(a maximum-to-median ratio of "
            f"{format_value(summary.get('max_to_median_ratio'), digits=2)})."
        )
    normalized_share = summary.get("top_10_normalized_expression_share")
    raw_share = summary.get("top_10_raw_target_count_share")
    if normalized_share is not None:
        share_text = (
            f"The {len(top)} highest-ranked populations account for "
            f"{_pct(normalized_share)}% of the normalized population-level expression signal"
        )
        if raw_share is not None:
            share_text += (
                f" and {_pct(raw_share)}% of the raw target-gene counts summed across "
                "valid populations"
            )
        parts.append(share_text + ".")
    excluded = list(analysis.get("excluded") or [])
    if excluded:
        parts.append(
            f"{len(excluded)} population columns were excluded from the ranking and the "
            "reasons are recorded in the population audit."
        )
    parts.append(
        "No verified confidence-interval method is published for these population "
        "estimates, so the ranking figure shows point estimates without intervals."
    )
    parts.append(
        f"These results describe the {DROPVIZ_SAMPLING_SCOPE}, and do not describe the "
        "human brain or regions the series did not sample."
    )
    return " ".join(parts)


def dropviz_unavailable_note(analysis: dict[str, Any] | None) -> str:
    """Honest unavailability / unresolved note for the DropViz branch."""
    reason = (analysis or {}).get("rank_status") or (analysis or {}).get("error_type")
    detail = (analysis or {}).get("rank_reason") or (analysis or {}).get("error_message")
    if reason == dg.RANK_STATUS_NORMALIZATION_UNRESOLVED:
        text = (
            f"A DropViz population ranking was not published for this run because the "
            f"{dg.GEO_ACCESSION} matrix value semantics remain unresolved, so no "
            "interpretable expression unit could be established"
        )
    else:
        text = (
            f"The DropViz {dg.GEO_ACCESSION} population ranking is unavailable for this "
            f"run ({reason or STATUS_SOURCE_UNAVAILABLE})"
        )
    return f"{text}{f': {detail}' if detail else ''}."


def build_therapeutic_narrative(
    *,
    gene_symbol: str,
    mouse_symbol: str | None,
    human: dict[str, Any] | None,
    mouse: dict[str, Any] | None,
    dropviz: dict[str, Any] | None,
) -> str:
    """Cautious target-localization context, derived only from this run's numbers."""
    localization: list[str] = []
    breadth: list[str] = []

    def _branch_phrase(
        analysis: dict[str, Any] | None, *, noun: str, label_fallback: str
    ) -> None:
        if not analysis or not analysis.get("ok"):
            return
        summary = dict(analysis.get("summary") or {})
        top = list(summary.get("top") or [])
        label = analysis.get("dataset_label") or label_fallback
        if top:
            leaders = ", ".join(str(item.get("label")) for item in top[:3])
            localization.append(f"in {label} the highest values fall in {leaders}")
        valid = summary.get(f"valid_{noun}_count")
        nonzero = summary.get(f"nonzero_{noun}_count")
        if valid:
            breadth.append(
                f"{nonzero} of {valid} sampled {label} populations carry nonzero "
                "aggregate expression"
            )

    _branch_phrase(human, noun="celltype", label_fallback="Human M1 10x")
    _branch_phrase(mouse, noun="cluster", label_fallback="mouse cortex and hippocampus")

    if dropviz and dropviz.get("rank_status") == dg.RANK_STATUS_SUCCESS:
        summary = dict(dropviz.get("ranking_summary") or {})
        top = list(summary.get("top_populations") or [])
        if top:
            leaders = ", ".join(str(item.get("population_label")) for item in top[:3])
            localization.append(
                f"in the DropViz {dg.GEO_ACCESSION} adult mouse-brain populations the "
                f"highest-ranked populations are {leaders}"
            )
        if summary.get("valid_population_count"):
            breadth.append(
                f"{summary.get('nonzero_population_count')} of "
                f"{summary.get('valid_population_count')} valid DropViz populations "
                "carry nonzero expression"
            )

    if not localization and not breadth:
        return (
            "No cell-type localization context could be assembled for this run, so no "
            f"therapeutic interpretation is offered for {gene_symbol}. "
            + THERAPEUTIC_CAVEAT
        )

    parts = [
        "These expression patterns identify sampled cell populations in which target "
        "engagement could produce on-target effects, assuming the therapeutic reaches "
        "the relevant tissue and successfully modulates the target. Populations with "
        "higher aggregate RNA expression may have greater potential for on-target "
        "pharmacology, but RNA abundance alone does not establish delivery, protein "
        "abundance, target engagement, functional dependency, or therapeutic response.",
    ]
    gene_phrase = gene_symbol + (f" (mouse {mouse_symbol})" if mouse_symbol else "")
    if localization:
        parts.append(
            f"Expression localization for {gene_phrase}: "
            + "; ".join(localization)
            + "."
        )
    if breadth:
        parts.append(
            "Expression breadth: "
            + "; ".join(breadth)
            + ". This indicates whether RNA expression is distributed across many "
            "sampled populations or is more localized within each dataset."
        )
        parts.append(
            "Potential safety context: populations outside an intended target set that "
            "also express the target could experience on-target effects if a "
            "therapeutic reaches and engages those cells. The relative RNA-expression "
            "rankings and maximum-to-median ratios reported above provide quantitative "
            "context for comparing expression localization within each dataset."
        )
    parts.append(THERAPEUTIC_CAVEAT)
    return " ".join(parts)


# --------------------------------------------------------------------------------------
# Node
# --------------------------------------------------------------------------------------
def node_generate_section_2c_derived_artifacts(
    state: DossierState,
    *,
    settings: Settings | None = None,
    persist_db: bool = True,
    transient: WorkflowTransientContext | None = None,
    config: Section2cConfig | None = None,
) -> DossierState:
    """Sole Section 2c owner: Allen Human M1, Allen mouse CTX-HPF, DropViz GEO."""
    _ = transient
    if state.get("run_type") != "section_bundle" or "2c" not in (
        state.get("selected_section_keys") or []
    ):
        return state

    cfg = settings or get_settings()
    section_cfg = config or Section2cConfig()
    run_id = state["dossier_run_id"]
    gene = state["gene_symbol"]
    evidence = list(state.get("evidence_records") or [])
    api_runs = list(state.get("api_runs") or [])
    raw_meta = list(state.get("raw_artifacts") or [])
    errors = list(state.get("errors") or [])
    coverage_extra = list(state.get("coverage") or [])
    item_key = f"celltype-{slugify(gene) or 'gene'}"

    audit: dict[str, Any] = {
        "network_owner": "section_2c",
        "calculation_version": CALCULATION_VERSION,
        "config": {
            "force_refresh": section_cfg.force_refresh,
            "top_n": section_cfg.top_n,
            "plot_dpi": section_cfg.plot_dpi,
            "attempt_allen_figures": section_cfg.attempt_allen_figures,
            "figure_max_attempts": section_cfg.figure_max_attempts,
            "attempt_live_dropviz": section_cfg.attempt_live_dropviz,
            "documented_semantics": section_cfg.documented_semantics,
            "documented_unit": section_cfg.documented_unit,
            "documentation_reference": section_cfg.documentation_reference,
            "output_root": str(section_cfg.output_root or cfg.output_path),
        },
        "sources": {},
        "allen_human": {},
        "allen_mouse": {},
        "dropviz": {},
        "figures": {},
        "artifacts": {},
        "artifact_root": str(cfg.raw_data_path),
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }

    rendering: dict[str, Any] = {
        "overall": "pending",
        "scientific_status": "pending",
        "visual_status": STATUS_NOT_ATTEMPTED,
        "allen_human_analysis_status": STATUS_SOURCE_UNAVAILABLE,
        "allen_human_scatter_status": STATUS_NOT_ATTEMPTED,
        "allen_human_heatmap_status": STATUS_NOT_ATTEMPTED,
        "allen_mouse_analysis_status": STATUS_SOURCE_UNAVAILABLE,
        "allen_mouse_scatter_status": STATUS_NOT_ATTEMPTED,
        "allen_mouse_heatmap_status": STATUS_NOT_ATTEMPTED,
        "dropviz_background_status": STATUS_SOURCE_UNAVAILABLE,
        "dropviz_geo_matrix_status": STATUS_SOURCE_UNAVAILABLE,
        "dropviz_rank_status": STATUS_SOURCE_UNAVAILABLE,
        "dropviz_rank_figure_status": STATUS_NOT_ATTEMPTED,
        "dropviz_confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
        "dropviz_regional_tsne_status": DROPVIZ_REGIONAL_TSNE_STATUS,
        "dropviz_display_label_status": DROPVIZ_DISPLAY_LABEL_STATUS,
        "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
        "label_mapping_status": STATUS_NOT_AVAILABLE,
        "value_semantics_status": STATUS_NOT_AVAILABLE,
        # Production 2c never reads data from the live DropViz Shiny server.
        "dropviz_live_status": STATUS_NOT_ATTEMPTED,
        # Historical DropViz application URL is audit-only, never polished evidence.
        "dropviz_historical_app_url": DROPVIZ_HISTORICAL_APP_URL,
    }
    unresolved: list[dict[str, str]] = []

    mouse_symbol = _resolve_mouse_symbol(
        state,
        gene=gene,
        settings=cfg,
        persist_db=persist_db,
        api_runs=api_runs,
        raw_meta=raw_meta,
        audit=audit,
    )
    if not mouse_symbol:
        unresolved.append(
            {
                "field": "mouse_symbol",
                "reason": (
                    "the mouse ortholog symbol could not be resolved from NCBI, so the "
                    "mouse datasets were queried with the requested symbol as supplied"
                ),
            }
        )
    mouse_query = mouse_symbol or gene

    paths = paths_for(section_cfg.output_root or cfg.output_path)
    gene_attempt = paths.new_gene_attempt(slugify(gene) or gene)
    audit["gene_attempt_dir"] = str(gene_attempt)

    # ------------------------------------------------------------------
    # Dataset-level sources (accepted-pointer reuse; refresh only on request)
    # ------------------------------------------------------------------
    payloads: dict[str, _SourcePayload] = {}
    for cache_key, source_key in _HUMAN_SOURCE_KEYS:
        official_url = ac.human_m1_source_url(source_key)
        payloads[cache_key] = _resolve_dataset_source(
            paths=paths,
            source_key=cache_key,
            official_url=official_url,
            force_refresh=section_cfg.force_refresh,
            downloader=(
                (lambda sk=source_key: ac.download_human_m1_source(
                    sk, gene_symbol=gene, settings=cfg
                ))
                if section_cfg.force_refresh
                else None
            ),
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            api_runs=api_runs,
            raw_meta=raw_meta,
        )
    for cache_key, _source_key in _MOUSE_SOURCE_KEYS:
        payloads[cache_key] = _resolve_dataset_source(
            paths=paths,
            source_key=cache_key,
            official_url=None,
            force_refresh=section_cfg.force_refresh,
            downloader=None,
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            api_runs=api_runs,
            raw_meta=raw_meta,
        )
    payloads[GEO_SOURCE_KEY] = _resolve_dataset_source(
        paths=paths,
        source_key=GEO_SOURCE_KEY,
        official_url=dg.supplementary_url(),
        force_refresh=section_cfg.force_refresh,
        downloader=(
            (lambda: dg.download_metacell_matrix(gene_symbol=mouse_query, settings=cfg))
            if section_cfg.force_refresh
            else None
        ),
        dossier_run_id=run_id,
        gene_symbol=gene,
        settings=cfg,
        persist_db=persist_db,
        api_runs=api_runs,
        raw_meta=raw_meta,
    )
    audit["sources"] = {key: payload.audit() for key, payload in payloads.items()}

    def _parents(*keys: str) -> list[str]:
        out: list[str] = []
        for key in keys:
            aid = payloads[key].raw_artifact_id if key in payloads else None
            if aid:
                out.append(str(aid))
        return out

    # ------------------------------------------------------------------
    # Branch 1: Allen Human M1 (single-nucleus RNA sequencing)
    # ------------------------------------------------------------------
    human_analysis: dict[str, Any] | None = None
    human_tm = payloads[ac.CACHE_KEY_HUMAN_TRIMMED_MEANS]
    human_tax = payloads[ac.CACHE_KEY_HUMAN_TAXONOMY]
    if human_tm.ok and human_tax.ok:
        human_analysis = analyze_allen_dataset(
            dataset=ac.DATASET_HUMAN_M1,
            gene_symbol=gene,
            trimmed_means=human_tm.content or b"",
            taxonomy_json=human_tax.content or b"",
            count_noun="celltype",
            top_n=section_cfg.top_n,
        )
    else:
        human_analysis = {
            "ok": False,
            "dataset": ac.DATASET_HUMAN_M1,
            "dataset_label": ac.DATASET_LABELS[ac.DATASET_HUMAN_M1],
            "error_type": "dataset_source_unavailable",
            "error_message": (
                human_tm.error_message or human_tax.error_message or "accepted source missing"
            ),
        }
    rendering["allen_human_analysis_status"] = (
        STATUS_SUCCESS if human_analysis.get("ok") else STATUS_SOURCE_UNAVAILABLE
    )
    audit["allen_human"] = {
        k: v for k, v in human_analysis.items() if k not in {"summary", "records"}
    }
    audit["allen_human"]["taxonomy_reconciliation"] = human_analysis.get(
        "taxonomy_reconciliation"
    )

    human_artifact_id: str | None = None
    human_evidence_rec: EvidenceRecord | None = None
    if human_analysis.get("ok"):
        artifact_payload = ac.build_celltype_summary_artifact(
            dataset=ac.DATASET_HUMAN_M1,
            summary=human_analysis["summary"],
            reconciliation=human_analysis["taxonomy_reconciliation"],
            match_count=int(human_analysis.get("match_count") or 0),
            source_checksums={
                ac.CACHE_KEY_HUMAN_TRIMMED_MEANS: human_tm.sha256,
                ac.CACHE_KEY_HUMAN_TAXONOMY: human_tax.sha256,
            },
            source_urls={
                ac.CACHE_KEY_HUMAN_TRIMMED_MEANS: human_tm.official_url,
                ac.CACHE_KEY_HUMAN_TAXONOMY: human_tax.official_url,
            },
        )
        art, meta = _persist_derived(
            dossier_run_id=run_id,
            source_name=ac.SOURCE_NAME,
            content=_json_bytes(artifact_payload),
            extension="json",
            media_type="application/json",
            filename_hint="allen-human-celltype-summary",
            artifact_role="allen_human_celltype_summary.json",
            parent_raw_artifact_ids=_parents(
                ac.CACHE_KEY_HUMAN_TRIMMED_MEANS, ac.CACHE_KEY_HUMAN_TAXONOMY
            ),
            settings=cfg,
            persist_db=persist_db,
            raw_meta=raw_meta,
            extra_notes={"dataset": ac.DATASET_HUMAN_M1},
        )
        if art is not None and meta is not None:
            human_artifact_id = art.id
            audit["artifacts"]["allen_human_celltype_summary.json"] = meta.get(
                "relative_path"
            )
        summary = human_analysis["summary"]
        human_evidence_rec = _evidence(
            dossier_run_id=run_id,
            gene_symbol=gene,
            source_name=ac.SOURCE_NAME,
            fact_type="allen_human_celltype_summary",
            key="allen-human-m1",
            value={
                "dataset": ac.DATASET_HUMAN_M1,
                "dataset_label": human_analysis["dataset_label"],
                "assay_terminology": human_analysis["assay_terminology"],
                "sampling_scope": human_analysis["sampling_scope"],
                "source_symbol": summary.get("source_symbol"),
                "valid_celltype_count": summary.get("valid_celltype_count"),
                "nonzero_celltype_count": summary.get("nonzero_celltype_count"),
                "nonzero_percentage": summary.get("nonzero_percentage"),
                "maximum": summary.get("maximum"),
                "median": summary.get("median"),
                "max_to_median_ratio": summary.get("max_to_median_ratio"),
                "top": summary.get("top"),
                "taxonomy_reconciliation": human_analysis["taxonomy_reconciliation"],
                "unit": TRIMMED_MEAN_VALUE_LABEL,
                "artifact_class": "derived",
                "presentation_item_key": item_key,
            },
            display_text=(
                f"{gene} Allen Human M1 nonzero aggregate expression in "
                f"{summary.get('nonzero_celltype_count')} of "
                f"{summary.get('valid_celltype_count')} transcriptomic cell types."
            ),
            api_run_id=None,
            raw_artifact_id=human_artifact_id,
        )
        _append_evidence(evidence, human_evidence_rec, persist_db=persist_db)
    else:
        coverage_extra.append(
            SourceCoverageResult(
                dossier_run_id=run_id,
                source_name=ac.SOURCE_NAME,
                status=SourceStatus.failed,
                evidence_record_count=0,
                error_message=(
                    human_analysis.get("error_message")
                    or "Allen Human M1 cell-type analysis unavailable"
                ),
                report_sections_supported=[SUBSECTION_2C],
            )
        )
        unresolved.append(
            {
                "field": "allen_human_analysis",
                "reason": str(
                    human_analysis.get("error_message")
                    or human_analysis.get("error_type")
                    or STATUS_SOURCE_UNAVAILABLE
                ),
            }
        )

    # ------------------------------------------------------------------
    # Branch 2: Allen mouse whole cortex and hippocampus 10x
    # ------------------------------------------------------------------
    mouse_tm = payloads[ac.CACHE_KEY_MOUSE_TRIMMED_MEANS]
    mouse_tax = payloads[ac.CACHE_KEY_MOUSE_TAXONOMY]
    if mouse_tm.ok and mouse_tax.ok:
        mouse_analysis = analyze_allen_dataset(
            dataset=ac.DATASET_MOUSE_CTX_HPF,
            gene_symbol=mouse_query,
            trimmed_means=mouse_tm.content or b"",
            taxonomy_json=mouse_tax.content or b"",
            count_noun="cluster",
            top_n=section_cfg.top_n,
        )
    else:
        mouse_analysis = {
            "ok": False,
            "dataset": ac.DATASET_MOUSE_CTX_HPF,
            "dataset_label": ac.DATASET_LABELS[ac.DATASET_MOUSE_CTX_HPF],
            "error_type": "dataset_source_unavailable",
            "error_message": (
                mouse_tm.error_message or mouse_tax.error_message or "accepted source missing"
            ),
        }
    rendering["allen_mouse_analysis_status"] = (
        STATUS_SUCCESS if mouse_analysis.get("ok") else STATUS_SOURCE_UNAVAILABLE
    )
    audit["allen_mouse"] = {
        k: v for k, v in mouse_analysis.items() if k not in {"summary", "records"}
    }

    mouse_artifact_id: str | None = None
    mouse_evidence_rec: EvidenceRecord | None = None
    if mouse_analysis.get("ok"):
        artifact_payload = ac.build_celltype_summary_artifact(
            dataset=ac.DATASET_MOUSE_CTX_HPF,
            summary=mouse_analysis["summary"],
            reconciliation=mouse_analysis["taxonomy_reconciliation"],
            match_count=int(mouse_analysis.get("match_count") or 0),
            source_checksums={
                ac.CACHE_KEY_MOUSE_TRIMMED_MEANS: mouse_tm.sha256,
                ac.CACHE_KEY_MOUSE_TAXONOMY: mouse_tax.sha256,
            },
            source_urls={
                ac.CACHE_KEY_MOUSE_TRIMMED_MEANS: mouse_tm.official_url,
                ac.CACHE_KEY_MOUSE_TAXONOMY: mouse_tax.official_url,
            },
        )
        art, meta = _persist_derived(
            dossier_run_id=run_id,
            source_name=ac.SOURCE_NAME,
            content=_json_bytes(artifact_payload),
            extension="json",
            media_type="application/json",
            filename_hint="allen-mouse-celltype-summary",
            artifact_role="allen_mouse_celltype_summary.json",
            parent_raw_artifact_ids=_parents(
                ac.CACHE_KEY_MOUSE_TRIMMED_MEANS, ac.CACHE_KEY_MOUSE_TAXONOMY
            ),
            settings=cfg,
            persist_db=persist_db,
            raw_meta=raw_meta,
            extra_notes={"dataset": ac.DATASET_MOUSE_CTX_HPF},
        )
        if art is not None and meta is not None:
            mouse_artifact_id = art.id
            audit["artifacts"]["allen_mouse_celltype_summary.json"] = meta.get(
                "relative_path"
            )
        summary = mouse_analysis["summary"]
        mouse_evidence_rec = _evidence(
            dossier_run_id=run_id,
            gene_symbol=gene,
            source_name=ac.SOURCE_NAME,
            fact_type="allen_mouse_celltype_summary",
            key="allen-mouse-ctx-hpf",
            value={
                "dataset": ac.DATASET_MOUSE_CTX_HPF,
                "dataset_label": mouse_analysis["dataset_label"],
                "assay_terminology": mouse_analysis["assay_terminology"],
                "sampling_scope": mouse_analysis["sampling_scope"],
                "requested_symbol": mouse_query,
                "source_symbol": summary.get("source_symbol"),
                "valid_cluster_count": summary.get("valid_cluster_count"),
                "nonzero_cluster_count": summary.get("nonzero_cluster_count"),
                "nonzero_percentage": summary.get("nonzero_percentage"),
                "maximum": summary.get("maximum"),
                "median": summary.get("median"),
                "max_to_median_ratio": summary.get("max_to_median_ratio"),
                "top": summary.get("top"),
                "taxonomy_reconciliation": mouse_analysis["taxonomy_reconciliation"],
                "unit": TRIMMED_MEAN_VALUE_LABEL,
                "artifact_class": "derived",
                "presentation_item_key": item_key,
            },
            display_text=(
                f"{mouse_query} Allen mouse cortex/hippocampus nonzero aggregate "
                f"expression in {summary.get('nonzero_cluster_count')} of "
                f"{summary.get('valid_cluster_count')} clusters."
            ),
            api_run_id=None,
            raw_artifact_id=mouse_artifact_id,
            organism="Mus musculus",
            taxon_id=10090,
        )
        _append_evidence(evidence, mouse_evidence_rec, persist_db=persist_db)
    else:
        coverage_extra.append(
            SourceCoverageResult(
                dossier_run_id=run_id,
                source_name=ac.SOURCE_NAME,
                status=SourceStatus.failed,
                evidence_record_count=0,
                error_message=(
                    mouse_analysis.get("error_message")
                    or "Allen mouse cortex/hippocampus analysis unavailable"
                ),
                report_sections_supported=[SUBSECTION_2C],
            )
        )
        unresolved.append(
            {
                "field": "allen_mouse_analysis",
                "reason": str(
                    mouse_analysis.get("error_message")
                    or mouse_analysis.get("error_type")
                    or STATUS_SOURCE_UNAVAILABLE
                ),
            }
        )

    # ------------------------------------------------------------------
    # Branch 3: DropViz GSE116470 published metacell matrix
    # ------------------------------------------------------------------
    geo = payloads[GEO_SOURCE_KEY]
    dropviz_analysis: dict[str, Any]
    if geo.ok:
        rendering["dropviz_geo_matrix_status"] = STATUS_SUCCESS
        dropviz_analysis = analyze_dropviz_matrix(
            content=geo.content or b"",
            mouse_gene_symbol=mouse_query,
            source_sha256=geo.sha256,
            source_url=geo.official_url,
            top_n=section_cfg.top_n,
            documented_semantics=section_cfg.documented_semantics,
            documented_unit=section_cfg.documented_unit,
            documentation_reference=section_cfg.documentation_reference,
        )
        if not dropviz_analysis.get("ok"):
            rendering["dropviz_geo_matrix_status"] = STATUS_SOURCE_UNAVAILABLE
    else:
        dropviz_analysis = {
            "ok": False,
            "error_type": geo.error_type or "dataset_source_unavailable",
            "error_message": geo.error_message or "accepted GEO matrix missing",
        }

    dropviz_parent_ids = _parents(GEO_SOURCE_KEY)
    dropviz_artifact_ids: dict[str, str] = {}
    dropviz_rank_rec: EvidenceRecord | None = None
    figure_meta: dict[str, Any] | None = None

    if dropviz_analysis.get("ok"):
        rendering["value_semantics_status"] = str(
            dropviz_analysis.get("value_semantics_status") or STATUS_NOT_AVAILABLE
        )
        rendering["label_mapping_status"] = str(
            dropviz_analysis.get("label_mapping_status") or STATUS_NOT_AVAILABLE
        )
        rendering["dropviz_rank_status"] = str(dropviz_analysis.get("rank_status"))
        if rendering["value_semantics_status"] not in dg.COUNT_COMPATIBLE_SEMANTICS and (
            rendering["value_semantics_status"] != dg.VALUE_SEMANTICS_NORMALIZED
        ):
            unresolved.append(
                {
                    "field": "value_semantics_status",
                    "reason": (
                        f"GSE116470 matrix values classified as "
                        f"{rendering['value_semantics_status']}, so no documented "
                        "expression unit could be established for ranking"
                    ),
                }
            )

        art, meta = _persist_derived(
            dossier_run_id=run_id,
            source_name=dg.SOURCE_NAME,
            content=_json_bytes(dropviz_analysis["matrix_profile"]),
            extension="json",
            media_type="application/json",
            filename_hint="dropviz-geo-matrix-profile",
            artifact_role="dropviz_geo_matrix_profile.json",
            parent_raw_artifact_ids=dropviz_parent_ids,
            settings=cfg,
            persist_db=persist_db,
            raw_meta=raw_meta,
            extra_notes={"accession": dg.GEO_ACCESSION},
        )
        if art is not None and meta is not None:
            rendering["dropviz_background_status"] = STATUS_SUCCESS
            dropviz_artifact_ids["dropviz_geo_matrix_profile.json"] = art.id
            audit["artifacts"]["dropviz_geo_matrix_profile.json"] = meta.get(
                "relative_path"
            )
            profile_rec = _evidence(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name=dg.SOURCE_NAME,
                fact_type="dropviz_geo_matrix_profile",
                key="dropviz-geo-matrix-profile",
                value={
                    **{
                        k: v
                        for k, v in dropviz_analysis["matrix_profile"].items()
                        if k != "value_semantics_evidence"
                    },
                    "assay_terminology": DROPVIZ_ASSAY_TERMINOLOGY,
                    "sampling_scope": DROPVIZ_SAMPLING_SCOPE,
                    "artifact_class": "derived",
                    "presentation_item_key": item_key,
                },
                display_text=(
                    f"{dg.GEO_ACCESSION} metacell matrix profile: "
                    f"{dropviz_analysis.get('gene_row_count')} gene rows x "
                    f"{dropviz_analysis.get('population_column_count')} populations, "
                    f"value semantics {dropviz_analysis.get('value_semantics_status')}."
                ),
                api_run_id=None,
                raw_artifact_id=art.id,
                organism="Mus musculus",
                taxon_id=10090,
            )
            _append_evidence(evidence, profile_rec, persist_db=persist_db)

        # Per-population CSVs: source order, then descending rank.
        records = list(dropviz_analysis.get("records") or [])
        if records:
            headers = [
                "population_label",
                "target_source_value",
                "population_total",
                "transcripts_per_100k",
                "ranking_value",
                "ranking_metric",
                "expression_unit",
            ]
            raw_rows = [
                [
                    r.get("population_label"),
                    r.get("target_source_value"),
                    r.get("population_total"),
                    r.get("transcripts_per_100k"),
                    r.get("ranking_value"),
                    r.get("ranking_metric"),
                    r.get("expression_unit"),
                ]
                for r in records
            ]
            art, meta = _persist_derived(
                dossier_run_id=run_id,
                source_name=dg.SOURCE_NAME,
                content=_csv_bytes(headers, raw_rows),
                extension="csv",
                media_type="text/csv",
                filename_hint="dropviz-population-expression-raw",
                artifact_role="dropviz_population_expression_raw.csv",
                parent_raw_artifact_ids=dropviz_parent_ids,
                settings=cfg,
                persist_db=persist_db,
                raw_meta=raw_meta,
                extra_notes={"row_order": "source_column_order"},
            )
            if art is not None and meta is not None:
                dropviz_artifact_ids["dropviz_population_expression_raw.csv"] = art.id
                audit["artifacts"]["dropviz_population_expression_raw.csv"] = meta.get(
                    "relative_path"
                )

            ranked = list(dropviz_analysis.get("ranked_records") or [])
            ranked_rows = [
                [index + 1, *row]
                for index, row in enumerate(
                    [
                        [
                            r.get("population_label"),
                            r.get("target_source_value"),
                            r.get("population_total"),
                            r.get("transcripts_per_100k"),
                            r.get("ranking_value"),
                            r.get("ranking_metric"),
                            r.get("expression_unit"),
                        ]
                        for r in ranked
                    ]
                )
            ]
            art, meta = _persist_derived(
                dossier_run_id=run_id,
                source_name=dg.SOURCE_NAME,
                content=_csv_bytes(["rank", *headers], ranked_rows),
                extension="csv",
                media_type="text/csv",
                filename_hint="dropviz-population-expression-ranked",
                artifact_role="dropviz_population_expression_ranked.csv",
                parent_raw_artifact_ids=dropviz_parent_ids,
                settings=cfg,
                persist_db=persist_db,
                raw_meta=raw_meta,
                extra_notes={"row_order": "descending_ranking_value"},
            )
            if art is not None and meta is not None:
                dropviz_artifact_ids["dropviz_population_expression_ranked.csv"] = art.id
                audit["artifacts"]["dropviz_population_expression_ranked.csv"] = meta.get(
                    "relative_path"
                )

        # Population audit: every exclusion and every unresolved decision.
        population_audit = {
            "calculation_version": CALCULATION_VERSION,
            "accession": dg.GEO_ACCESSION,
            "source_sha256": geo.sha256,
            "requested_mouse_symbol": mouse_query,
            "target_row_match_count": dropviz_analysis.get("target_row_match_count"),
            "target_source_symbol": dropviz_analysis.get("target_source_symbol"),
            "value_semantics_status": dropviz_analysis.get("value_semantics_status"),
            "value_semantics_basis": dropviz_analysis.get("value_semantics_basis"),
            "value_semantics_evidence": (
                dropviz_analysis.get("matrix_profile") or {}
            ).get("value_semantics_evidence"),
            "rank_status": dropviz_analysis.get("rank_status"),
            "rank_reason": dropviz_analysis.get("rank_reason"),
            "label_mapping_status": dropviz_analysis.get("label_mapping_status"),
            "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
            "confidence_interval_note": (
                "No verified DropViz confidence-interval formula is vendored, so no "
                "interval is computed or drawn."
            ),
            "valid_population_count": (
                dropviz_analysis.get("ranking_summary") or {}
            ).get("valid_population_count"),
            "excluded_population_count": len(dropviz_analysis.get("excluded") or []),
            "excluded_populations": dropviz_analysis.get("excluded"),
            "parsed_population_labels": dropviz_analysis.get("parsed_labels"),
        }
        art, meta = _persist_derived(
            dossier_run_id=run_id,
            source_name=dg.SOURCE_NAME,
            content=_json_bytes(population_audit),
            extension="json",
            media_type="application/json",
            filename_hint="dropviz-population-audit",
            artifact_role="dropviz_population_audit.json",
            parent_raw_artifact_ids=dropviz_parent_ids,
            settings=cfg,
            persist_db=persist_db,
            raw_meta=raw_meta,
        )
        if art is not None and meta is not None:
            dropviz_artifact_ids["dropviz_population_audit.json"] = art.id
            audit["artifacts"]["dropviz_population_audit.json"] = meta.get("relative_path")

        if dropviz_analysis.get("rank_status") == dg.RANK_STATUS_SUCCESS:
            ranking_summary = dropviz_analysis["ranking_summary"]
            top_payload = {
                **ranking_summary,
                "accession": dg.GEO_ACCESSION,
                "source_sha256": geo.sha256,
                "requested_mouse_symbol": mouse_query,
                "label_mapping_status": dropviz_analysis.get("label_mapping_status"),
                "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
                "presentation": dropviz_analysis.get("presentation"),
            }
            art, meta = _persist_derived(
                dossier_run_id=run_id,
                source_name=dg.SOURCE_NAME,
                content=_json_bytes(top_payload),
                extension="json",
                media_type="application/json",
                filename_hint="dropviz-top-populations",
                artifact_role="dropviz_top_populations.json",
                parent_raw_artifact_ids=dropviz_parent_ids,
                settings=cfg,
                persist_db=persist_db,
                raw_meta=raw_meta,
            )
            top_artifact_id = None
            if art is not None and meta is not None:
                top_artifact_id = art.id
                dropviz_artifact_ids["dropviz_top_populations.json"] = art.id
                audit["artifacts"]["dropviz_top_populations.json"] = meta.get(
                    "relative_path"
                )
            dropviz_rank_rec = _evidence(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name=dg.SOURCE_NAME,
                fact_type="dropviz_top_populations",
                key="dropviz-top-populations",
                value={
                    "accession": dg.GEO_ACCESSION,
                    "requested_mouse_symbol": mouse_query,
                    "assay_terminology": DROPVIZ_ASSAY_TERMINOLOGY,
                    "sampling_scope": DROPVIZ_SAMPLING_SCOPE,
                    "value_semantics_status": dropviz_analysis.get(
                        "value_semantics_status"
                    ),
                    "ranking_metric": ranking_summary.get("ranking_metric"),
                    "expression_unit": ranking_summary.get("expression_unit"),
                    "valid_population_count": ranking_summary.get(
                        "valid_population_count"
                    ),
                    "nonzero_population_count": ranking_summary.get(
                        "nonzero_population_count"
                    ),
                    "nonzero_percentage": ranking_summary.get("nonzero_percentage"),
                    "maximum": ranking_summary.get("maximum"),
                    "median": ranking_summary.get("median"),
                    "max_to_median_ratio": ranking_summary.get("max_to_median_ratio"),
                    "top_populations": ranking_summary.get("top_populations"),
                    "top_10_normalized_expression_share": ranking_summary.get(
                        "top_10_normalized_expression_share"
                    ),
                    "top_10_raw_target_count_share": ranking_summary.get(
                        "top_10_raw_target_count_share"
                    ),
                    "label_mapping_status": dropviz_analysis.get("label_mapping_status"),
                    "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
                    "artifact_class": "derived",
                    "presentation_item_key": item_key,
                },
                display_text=(
                    f"{mouse_query} DropViz {dg.GEO_ACCESSION} nonzero expression in "
                    f"{ranking_summary.get('nonzero_population_count')} of "
                    f"{ranking_summary.get('valid_population_count')} populations."
                ),
                api_run_id=None,
                raw_artifact_id=top_artifact_id,
                organism="Mus musculus",
                taxon_id=10090,
            )
            _append_evidence(evidence, dropviz_rank_rec, persist_db=persist_db)

            # Chart only when an interpretable ranking metric AND unit both exist.
            presentation = dict(dropviz_analysis.get("presentation") or {})
            axis_label = presentation.get("axis_label")
            if presentation.get("chartable") and axis_label:
                try:
                    png = dg.render_top_populations_png(
                        mouse_gene_symbol=mouse_query,
                        top_populations=list(
                            ranking_summary.get("top_populations") or []
                        ),
                        axis_label=str(axis_label),
                        confidence_interval_status=CONFIDENCE_INTERVAL_STATUS,
                        dpi=section_cfg.plot_dpi,
                    )
                    fig_art, figure_meta = _persist_derived(
                        dossier_run_id=run_id,
                        source_name=dg.SOURCE_NAME,
                        content=png,
                        extension="png",
                        media_type="image/png",
                        filename_hint=f"{CHART_VERSION}-{slugify(mouse_query)}",
                        artifact_role=f"{CHART_VERSION}.png",
                        parent_raw_artifact_ids=[
                            *dropviz_parent_ids,
                            *(
                                [dropviz_artifact_ids["dropviz_top_populations.json"]]
                                if "dropviz_top_populations.json" in dropviz_artifact_ids
                                else []
                            ),
                        ],
                        settings=cfg,
                        persist_db=persist_db,
                        raw_meta=raw_meta,
                        extra_notes={
                            "chart_version": CHART_VERSION,
                            "derivation_type": "dropviz_population_point_chart",
                            "axis_label": str(axis_label),
                            "expression_unit": presentation.get("expression_unit"),
                            "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
                            "point_only": True,
                        },
                        validate=_validate_nonblank_image,
                    )
                    if fig_art is not None and figure_meta is not None:
                        rendering["dropviz_rank_figure_status"] = STATUS_SUCCESS
                        dropviz_artifact_ids[f"{CHART_VERSION}.png"] = fig_art.id
                        audit["artifacts"][f"{CHART_VERSION}.png"] = figure_meta.get(
                            "relative_path"
                        )
                        fig_rec = _evidence(
                            dossier_run_id=run_id,
                            gene_symbol=gene,
                            source_name=dg.SOURCE_NAME,
                            fact_type="dropviz_top_populations_figure",
                            key="dropviz-top-populations-figure",
                            value={
                                "relative_path": figure_meta.get("relative_path"),
                                "sha256": figure_meta.get("content_hash")
                                or fig_art.content_hash,
                                "width": figure_meta.get("width"),
                                "height": figure_meta.get("height"),
                                "byte_size": figure_meta.get("byte_size"),
                                "media_type": "image/png",
                                "artifact_class": "derived",
                                "derivation_type": "dropviz_population_point_chart",
                                "chart_version": CHART_VERSION,
                                "axis_label": str(axis_label),
                                "expression_unit": presentation.get("expression_unit"),
                                "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
                                "point_only": True,
                                "presentation_item_key": f"{item_key}-dropviz-rank",
                                "figure_raw_artifact_id": fig_art.id,
                            },
                            display_text=(
                                f"{mouse_query} DropViz highest-expressing populations "
                                f"({axis_label}), point estimates without intervals."
                            ),
                            api_run_id=None,
                            raw_artifact_id=fig_art.id,
                            organism="Mus musculus",
                            taxon_id=10090,
                        )
                        _append_evidence(evidence, fig_rec, persist_db=persist_db)
                    else:
                        rendering["dropviz_rank_figure_status"] = STATUS_SOURCE_UNAVAILABLE
                except Exception as exc:  # noqa: BLE001
                    logger.exception("Section 2c DropViz rank chart failed for %s", gene)
                    rendering["dropviz_rank_figure_status"] = STATUS_SOURCE_UNAVAILABLE
                    errors.append(f"Section 2c DropViz rank chart failed: {exc}")
            else:
                rendering["dropviz_rank_figure_status"] = STATUS_FIGURE_SUPPRESSED
        else:
            # No interpretable metric: emit no chart at all.
            rendering["dropviz_rank_figure_status"] = STATUS_FIGURE_SUPPRESSED
            unresolved.append(
                {
                    "field": "dropviz_rank_status",
                    "reason": str(
                        dropviz_analysis.get("rank_reason")
                        or dropviz_analysis.get("rank_status")
                    ),
                }
            )
    else:
        coverage_extra.append(
            SourceCoverageResult(
                dossier_run_id=run_id,
                source_name=dg.SOURCE_NAME,
                status=SourceStatus.failed,
                evidence_record_count=0,
                error_message=(
                    dropviz_analysis.get("error_message")
                    or f"DropViz {dg.GEO_ACCESSION} matrix unavailable"
                ),
                report_sections_supported=[SUBSECTION_2C],
            )
        )
        unresolved.append(
            {
                "field": "dropviz_geo_matrix",
                "reason": str(
                    dropviz_analysis.get("error_message")
                    or dropviz_analysis.get("error_type")
                    or STATUS_SOURCE_UNAVAILABLE
                ),
            }
        )

    audit["dropviz"] = {
        k: v
        for k, v in dropviz_analysis.items()
        if k not in {"records", "ranked_records", "parsed_labels", "matrix_profile"}
    }
    audit["dropviz"]["artifact_ids"] = dropviz_artifact_ids
    if section_cfg.attempt_live_dropviz:
        # Only an actually invoked-and-failed optional live client may report a
        # source failure; the default path never touches the Shiny server.
        rendering["dropviz_live_status"] = STATUS_SOURCE_UNAVAILABLE
        audit["dropviz"]["live_client"] = {
            "invoked": True,
            "status": STATUS_SOURCE_UNAVAILABLE,
            "note": (
                "The optional live DropViz Shiny client is not used as a data source by "
                "Section 2c; it was requested and reported unavailable."
            ),
        }
    else:
        audit["dropviz"]["live_client"] = {
            "invoked": False,
            "status": STATUS_NOT_ATTEMPTED,
            "note": (
                "Production Section 2c reads the published GEO matrix and never calls "
                "the live DropViz Shiny server."
            ),
        }

    # ------------------------------------------------------------------
    # Optional Allen explorer figures (never gate a structured analysis)
    # ------------------------------------------------------------------
    figures: dict[str, dict[str, Any] | None] = {}
    figure_notes: dict[str, list[str]] = {"human": [], "mouse": []}
    figure_specs = (
        ("human", ac.DATASET_HUMAN_M1, gene, human_analysis.get("ok")),
        ("mouse", ac.DATASET_MOUSE_CTX_HPF, mouse_query, mouse_analysis.get("ok")),
    )
    for species, dataset, symbol, analysis_ok in figure_specs:
        for visualization in _FIGURE_VISUALIZATIONS:
            status_key = f"allen_{species}_{visualization}_status"
            figures[f"{species}_{visualization}"] = None
            if not section_cfg.attempt_allen_figures or not analysis_ok:
                rendering[status_key] = STATUS_NOT_ATTEMPTED
                audit["figures"][f"{species}_{visualization}"] = {
                    "status": STATUS_NOT_ATTEMPTED,
                    "reason": (
                        "figure capture disabled by configuration"
                        if not section_cfg.attempt_allen_figures
                        else "structured analysis unavailable, so no figure was requested"
                    ),
                }
                continue
            captured, fig_record, note = _capture_allen_figure(
                dataset=dataset,
                gene_symbol=symbol,
                visualization=visualization,
                dossier_run_id=run_id,
                report_gene_symbol=gene,
                item_key=item_key,
                section_cfg=section_cfg,
                settings=cfg,
                persist_db=persist_db,
                api_runs=api_runs,
                raw_meta=raw_meta,
                evidence=evidence,
                audit=audit,
                parent_raw_artifact_ids=(
                    _parents(ac.CACHE_KEY_HUMAN_TRIMMED_MEANS)
                    if species == "human"
                    else _parents(ac.CACHE_KEY_MOUSE_TRIMMED_MEANS)
                ),
            )
            rendering[status_key] = (
                STATUS_SUCCESS if captured else STATUS_SOURCE_UNAVAILABLE
            )
            figures[f"{species}_{visualization}"] = fig_record
            if not captured:
                figure_notes[species].append(note)
                unresolved.append(
                    {
                        "field": status_key,
                        "reason": note,
                    }
                )

    # ------------------------------------------------------------------
    # Narratives + presentation summary
    # ------------------------------------------------------------------
    human_ok = bool(human_analysis.get("ok"))
    mouse_ok = bool(mouse_analysis.get("ok"))
    dropviz_ok = dropviz_analysis.get("rank_status") == dg.RANK_STATUS_SUCCESS

    dropviz_presentation = dict(dropviz_analysis.get("presentation") or {})
    dropviz_summary = dict(dropviz_analysis.get("ranking_summary") or {})
    axis_label = dropviz_presentation.get("axis_label")

    human_block = _branch_presentation(
        analysis=human_analysis,
        gene_symbol=gene,
        label_column=HUMAN_LABEL_COLUMN,
        count_noun="celltype",
        explorer_symbol=gene,
        figure_notes=figure_notes["human"],
        figure_statuses={
            visualization: rendering[f"allen_human_{visualization}_status"]
            for visualization in _FIGURE_VISUALIZATIONS
        },
        top_n=section_cfg.top_n,
    )
    mouse_block = _branch_presentation(
        analysis=mouse_analysis,
        gene_symbol=mouse_query,
        label_column=MOUSE_LABEL_COLUMN,
        count_noun="cluster",
        explorer_symbol=mouse_query,
        figure_notes=figure_notes["mouse"],
        figure_statuses={
            visualization: rendering[f"allen_mouse_{visualization}_status"]
            for visualization in _FIGURE_VISUALIZATIONS
        },
        top_n=section_cfg.top_n,
    )

    dropviz_block: dict[str, Any] = {
        "analysis_status": STATUS_SUCCESS if dropviz_ok else STATUS_SOURCE_UNAVAILABLE,
        "rank_status": str(
            dropviz_analysis.get("rank_status")
            or dropviz_analysis.get("error_type")
            or STATUS_SOURCE_UNAVAILABLE
        ),
        "accession": dg.GEO_ACCESSION,
        "source_sha256": geo.sha256,
        "source_url": geo.official_url,
        "assay_terminology": DROPVIZ_ASSAY_TERMINOLOGY,
        "sampling_scope": DROPVIZ_SAMPLING_SCOPE,
        "value_semantics_status": dropviz_analysis.get("value_semantics_status"),
        "ranking_metric": dropviz_presentation.get("ranking_metric"),
        "expression_unit": dropviz_presentation.get("expression_unit"),
        "axis_label": axis_label,
        "value_label": str(axis_label or DROPVIZ_FALLBACK_VALUE_LABEL),
        "label_column_header": DROPVIZ_LABEL_COLUMN,
        "label_mapping_status": dropviz_analysis.get("label_mapping_status"),
        "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
        "valid_population_count": dropviz_summary.get("valid_population_count"),
        "nonzero_population_count": dropviz_summary.get("nonzero_population_count"),
        "nonzero_percentage": dropviz_summary.get("nonzero_percentage"),
        "maximum": dropviz_summary.get("maximum"),
        "median": dropviz_summary.get("median"),
        "max_to_median_ratio": dropviz_summary.get("max_to_median_ratio"),
        "top_10_normalized_expression_share": dropviz_summary.get(
            "top_10_normalized_expression_share"
        ),
        "top_10_raw_target_count_share": dropviz_summary.get(
            "top_10_raw_target_count_share"
        ),
        "excluded_population_count": len(dropviz_analysis.get("excluded") or []),
        "target_row_match_count": dropviz_analysis.get("target_row_match_count"),
        "geo_url": GEO_SERIES_URL,
        "geo_link_label": f"NCBI GEO {dg.GEO_ACCESSION}",
        "geo_attribution": (
            f"Source data: NCBI GEO {dg.GEO_ACCESSION} processed metacell matrix."
        ),
        "population_identifier_note": DROPVIZ_POPULATION_IDENTIFIER_NOTE,
        "regional_tsne_limitation_note": DROPVIZ_REGIONAL_TSNE_LIMITATION_NOTE,
        "figure_status": rendering["dropviz_rank_figure_status"],
        "figure_status_notes": [],
        "regional_tsne_status": DROPVIZ_REGIONAL_TSNE_STATUS,
        "display_label_status": rendering.get("dropviz_display_label_status"),
        "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
        "live_status": rendering.get("dropviz_live_status"),
        "top_populations": [
            {
                "population_label": row.get("population_label"),
                "ranking_value": row.get("ranking_value"),
                "target_source_value": row.get("target_source_value"),
                "population_total": row.get("population_total"),
                "value_display": format_value(row.get("ranking_value"), digits=2),
            }
            for row in list(dropviz_summary.get("top_populations") or [])
        ],
    }
    if dropviz_ok:
        dropviz_block["narrative"] = dropviz_narrative(
            dropviz_analysis, mouse_gene_symbol=mouse_query
        )
        dropviz_block["status_note"] = None
    else:
        dropviz_block["narrative"] = dropviz_unavailable_note(dropviz_analysis)
        dropviz_block["status_note"] = dropviz_unavailable_note(dropviz_analysis)
    if rendering["dropviz_rank_figure_status"] == STATUS_FIGURE_SUPPRESSED:
        dropviz_block["figure_status_notes"].append(
            "No DropViz ranking figure was generated because no interpretable "
            "expression metric and unit were established for the matrix values."
        )
    elif rendering["dropviz_rank_figure_status"] == STATUS_SOURCE_UNAVAILABLE:
        dropviz_block["figure_status_notes"].append(
            "The DropViz ranking figure could not be rendered for this run."
        )

    therapeutic = build_therapeutic_narrative(
        gene_symbol=gene,
        mouse_symbol=mouse_symbol,
        human=human_analysis,
        mouse=mouse_analysis,
        dropviz=dropviz_analysis,
    )

    # Scientific and visual statuses are independent. Visual-complete acceptance
    # of the gene pointer happens only after report/PDF rendering.
    figure_keys = (
        "allen_human_scatter_status",
        "allen_human_heatmap_status",
        "allen_mouse_scatter_status",
        "allen_mouse_heatmap_status",
    )
    if human_ok and mouse_ok and dropviz_ok:
        scientific_status = STATUS_SUCCESS
    elif not any([human_ok, mouse_ok, dropviz_ok]):
        scientific_status = "empty"
    else:
        scientific_status = STATUS_PARTIAL

    attempted_figures = [
        key for key in figure_keys if rendering[key] != STATUS_NOT_ATTEMPTED
    ]
    if not attempted_figures and not section_cfg.attempt_allen_figures:
        visual_status = STATUS_NOT_ATTEMPTED
    else:
        success_count = sum(
            1 for key in figure_keys if rendering[key] == STATUS_SUCCESS
        )
        if success_count == 4:
            visual_status = STATUS_SUCCESS
        elif success_count == 0:
            visual_status = STATUS_UNAVAILABLE
        else:
            visual_status = STATUS_PARTIAL

    # overall remains a coarse compatibility rollup for existing callers.
    if scientific_status == STATUS_SUCCESS and visual_status in {
        STATUS_SUCCESS,
        STATUS_NOT_ATTEMPTED,
    }:
        overall = STATUS_SUCCESS
    elif scientific_status == "empty" and visual_status in {
        STATUS_UNAVAILABLE,
        STATUS_NOT_ATTEMPTED,
    }:
        overall = "empty"
    else:
        overall = STATUS_PARTIAL
    rendering["scientific_status"] = scientific_status
    rendering["visual_status"] = visual_status
    rendering["overall"] = overall
    rendering["dropviz_confidence_interval_status"] = CONFIDENCE_INTERVAL_STATUS
    rendering["dropviz_regional_tsne_status"] = DROPVIZ_REGIONAL_TSNE_STATUS
    rendering["dropviz_display_label_status"] = str(
        dropviz_analysis.get("label_mapping_status") or DROPVIZ_DISPLAY_LABEL_STATUS
    )

    summary_payload: dict[str, Any] = {
        "gene_symbol": gene,
        "mouse_symbol": mouse_symbol,
        "presentation_item_key": item_key,
        "intro_text": SECTION_2C_INTRO_TEXT,
        "overall_status": overall,
        "therapeutic_narrative": therapeutic,
        "unresolved_issues": unresolved,
        "human": human_block,
        "mouse": mouse_block,
        "dropviz": dropviz_block,
    }

    # ------------------------------------------------------------------
    # Compact evidence model (no raw matrices)
    # ------------------------------------------------------------------
    evidence_model = {
        "calculation_version": CALCULATION_VERSION,
        "section_key": "2c",
        "subsection": SUBSECTION_2C,
        "requested_gene_symbol": gene,
        "resolved_mouse_symbol": mouse_symbol,
        "queried_mouse_symbol": mouse_query,
        "overall_status": overall,
        "statuses": dict(rendering),
        "allen_human": _branch_evidence_model(
            analysis=human_analysis,
            block=human_block,
            count_noun="celltype",
            figures={
                "scatter": figures.get("human_scatter"),
                "heatmap": figures.get("human_heatmap"),
            },
            figure_statuses={
                "scatter": rendering["allen_human_scatter_status"],
                "heatmap": rendering["allen_human_heatmap_status"],
            },
            source_checksums={
                ac.CACHE_KEY_HUMAN_TRIMMED_MEANS: human_tm.sha256,
                ac.CACHE_KEY_HUMAN_TAXONOMY: human_tax.sha256,
            },
            artifact_relative_path=audit["artifacts"].get(
                "allen_human_celltype_summary.json"
            ),
        ),
        "allen_mouse": _branch_evidence_model(
            analysis=mouse_analysis,
            block=mouse_block,
            count_noun="cluster",
            figures={
                "scatter": figures.get("mouse_scatter"),
                "heatmap": figures.get("mouse_heatmap"),
            },
            figure_statuses={
                "scatter": rendering["allen_mouse_scatter_status"],
                "heatmap": rendering["allen_mouse_heatmap_status"],
            },
            source_checksums={
                ac.CACHE_KEY_MOUSE_TRIMMED_MEANS: mouse_tm.sha256,
                ac.CACHE_KEY_MOUSE_TAXONOMY: mouse_tax.sha256,
            },
            artifact_relative_path=audit["artifacts"].get(
                "allen_mouse_celltype_summary.json"
            ),
        ),
        "dropviz": {
            "accession": dg.GEO_ACCESSION,
            "source_sha256": geo.sha256,
            "source_url": geo.official_url,
            "assay_terminology": DROPVIZ_ASSAY_TERMINOLOGY,
            "sampling_scope": DROPVIZ_SAMPLING_SCOPE,
            "geo_matrix_status": rendering["dropviz_geo_matrix_status"],
            "background_status": rendering["dropviz_background_status"],
            "rank_status": rendering["dropviz_rank_status"],
            "rank_figure_status": rendering["dropviz_rank_figure_status"],
            "live_status": rendering["dropviz_live_status"],
            "value_semantics_status": rendering["value_semantics_status"],
            "label_mapping_status": rendering["label_mapping_status"],
            "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
            "gene_row_count": dropviz_analysis.get("gene_row_count"),
            "population_column_count": dropviz_analysis.get("population_column_count"),
            "target_row_match_count": dropviz_analysis.get("target_row_match_count"),
            "top_population_summary": {
                "valid_population_count": dropviz_block["valid_population_count"],
                "nonzero_population_count": dropviz_block["nonzero_population_count"],
                "nonzero_percentage": dropviz_block["nonzero_percentage"],
                "maximum": dropviz_block["maximum"],
                "median": dropviz_block["median"],
                "max_to_median_ratio": dropviz_block["max_to_median_ratio"],
                "ranking_metric": dropviz_block["ranking_metric"],
                "expression_unit": dropviz_block["expression_unit"],
                "top_populations": dropviz_block["top_populations"],
                "top_10_normalized_expression_share": dropviz_block[
                    "top_10_normalized_expression_share"
                ],
                "top_10_raw_target_count_share": dropviz_block[
                    "top_10_raw_target_count_share"
                ],
                "excluded_population_count": dropviz_block["excluded_population_count"],
            },
            "therapeutic_localization_summary": dropviz_block["narrative"],
            "rank_figure": (
                {
                    "relative_path": figure_meta.get("relative_path"),
                    "sha256": figure_meta.get("content_hash"),
                    "width": figure_meta.get("width"),
                    "height": figure_meta.get("height"),
                    "chart_version": CHART_VERSION,
                    "axis_label": axis_label,
                    "ranking_metric": dropviz_block["ranking_metric"],
                    "confidence_interval_status": CONFIDENCE_INTERVAL_STATUS,
                    "confidence_intervals_drawn": False,
                    "point_only": True,
                }
                if figure_meta
                else None
            ),
            "artifacts": {
                name: audit["artifacts"].get(name)
                for name in (
                    "dropviz_geo_matrix_profile.json",
                    "dropviz_population_expression_raw.csv",
                    "dropviz_population_expression_ranked.csv",
                    "dropviz_top_populations.json",
                    "dropviz_population_audit.json",
                    f"{CHART_VERSION}.png",
                )
            },
        },
        "intro_text": SECTION_2C_INTRO_TEXT,
        "therapeutic_narrative": therapeutic,
        "unresolved_issues": unresolved,
    }
    art, meta = _persist_derived(
        dossier_run_id=run_id,
        source_name=ac.SOURCE_NAME,
        content=_json_bytes(evidence_model),
        extension="json",
        media_type="application/json",
        filename_hint="section-2c-evidence",
        artifact_role="section_2c_evidence.json",
        parent_raw_artifact_ids=[
            i
            for i in [
                human_artifact_id,
                mouse_artifact_id,
                *dropviz_artifact_ids.values(),
            ]
            if i
        ],
        settings=cfg,
        persist_db=persist_db,
        raw_meta=raw_meta,
    )
    evidence_artifact_id = art.id if art is not None else None
    if meta is not None:
        audit["artifacts"]["section_2c_evidence.json"] = meta.get("relative_path")

    summary_rec = _evidence(
        dossier_run_id=run_id,
        gene_symbol=gene,
        source_name=ac.SOURCE_NAME,
        fact_type="section_2c_summary",
        key="section-2c-summary",
        value={
            "gene_symbol": gene,
            "mouse_symbol": mouse_symbol,
            "presentation_item_key": item_key,
            "overall_status": overall,
            "statuses": dict(rendering),
            "intro_text": SECTION_2C_INTRO_TEXT,
            "therapeutic_narrative": therapeutic,
            "unresolved_issues": unresolved,
            "artifact_class": "derived",
        },
        display_text=(
            f"{gene} single-cell / single-nucleus cell-type expression summary "
            f"({overall})."
        ),
        api_run_id=None,
        raw_artifact_id=evidence_artifact_id,
    )
    _append_evidence(evidence, summary_rec, persist_db=persist_db)

    section_status = {
        "rendering_status": rendering,
        "summary": summary_payload,
        "audit": audit,
    }

    # Immutable gene attempt + optional accepted pointer. Dataset matrices stay
    # under accepted/sources; this pointer only covers the gene-specific report.
    write_json_atomic(gene_attempt / "section_2c_status.json", section_status)
    write_json_atomic(
        gene_attempt / "manifest.json",
        {
            "gene_symbol": gene,
            "mouse_symbol": mouse_symbol,
            "overall_status": overall,
            "rendering_status": rendering,
            "attempt_dir": str(gene_attempt),
            "artifacts": dict(audit.get("artifacts") or {}),
        },
    )
    _maybe_accept_gene_report(
        paths,
        gene_symbol=gene,
        attempt_dir=gene_attempt,
        overall=overall,
        rendering=rendering,
        artifacts=dict(audit.get("artifacts") or {}),
    )

    return {
        **state,
        "evidence_records": evidence,
        "api_runs": api_runs,
        "raw_artifacts": raw_meta,
        "errors": errors,
        "coverage": coverage_extra,
        "section_2c_status": section_status,
    }


def scientific_core_accepted(rendering: dict[str, Any]) -> bool:
    """True when both Allen analyses and a DropViz ranking metric succeeded."""
    return (
        str(rendering.get("allen_human_analysis_status") or "") == STATUS_SUCCESS
        and str(rendering.get("allen_mouse_analysis_status") or "") == STATUS_SUCCESS
        and str(rendering.get("dropviz_rank_status") or "") == STATUS_SUCCESS
    )


def _maybe_accept_gene_report(
    paths: Section2cPaths,
    *,
    gene_symbol: str,
    attempt_dir: Path,
    overall: str,
    rendering: dict[str, Any],
    artifacts: dict[str, Any],
) -> Path | None:
    """Section 2c never pins visual-complete acceptance.

    The gene-attempt manifest is immutable; ``accepted/genes/<GENE>.json`` is
    written only after report/PDF rendering confirms embedded figure roles.
    """
    _ = (paths, gene_symbol, attempt_dir, overall, rendering, artifacts)
    return None


def evaluate_section_2c_visual_complete(
    *,
    rendering: dict[str, Any],
    embedded_figure_roles: set[str] | frozenset[str],
    pdf_render_status: str,
) -> dict[str, Any]:
    """Decide visual-complete acceptance after presentation + PDF rendering."""
    scientific = str(rendering.get("scientific_status") or "")
    visual = str(rendering.get("visual_status") or "")
    roles = set(embedded_figure_roles or ())
    missing = sorted(SECTION_2C_VISUAL_COMPLETE_ROLES - roles)
    visual_complete = (
        scientific == STATUS_SUCCESS
        and visual == STATUS_SUCCESS
        and pdf_render_status == STATUS_SUCCESS
        and SECTION_2C_VISUAL_COMPLETE_ROLES.issubset(roles)
    )
    return {
        "visual_complete": visual_complete,
        "scientific_status": scientific,
        "visual_status": visual,
        "pdf_render_status": pdf_render_status,
        "embedded_figure_roles": sorted(roles),
        "missing_figure_roles": missing,
        "required_roles": sorted(SECTION_2C_VISUAL_COMPLETE_ROLES),
    }


def accept_visual_complete_gene_report(
    paths: Section2cPaths,
    *,
    gene_symbol: str,
    attempt_dir: Path,
    rendering: dict[str, Any],
    artifacts: dict[str, Any] | None,
    evaluation: dict[str, Any],
    promote_existing: bool = False,
) -> Path | None:
    """Atomically pin a visual-complete gene pointer; never downgrade an accepted one."""
    if not evaluation.get("visual_complete"):
        return None
    pointer = paths.accepted_gene_pointer(gene_symbol)
    replaced_prior_visual_complete = False
    if pointer.is_file():
        try:
            existing = json.loads(pointer.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            existing = {}
        prior = dict(existing.get("acceptance") or {})
        if prior.get("section_2c_visual_complete") is True:
            if not promote_existing:
                # Keep the first accepted visual-complete pointer.
                return None
            replaced_prior_visual_complete = True
    return accept_gene_report(
        paths,
        gene_symbol=gene_symbol,
        attempt_dir=attempt_dir,
        acceptance={
            "overall_status": STATUS_SUCCESS,
            "scientific_status": rendering.get("scientific_status"),
            "visual_status": rendering.get("visual_status"),
            "scientific_core_accepted": True,
            "section_2c_visual_complete": True,
            "promotion_requested": bool(promote_existing),
            "replaced_prior_visual_complete": replaced_prior_visual_complete,
            "rendering_status": dict(rendering),
            "evaluation": dict(evaluation),
        },
        artifacts=artifacts or {},
    )


def _allen_short_intro(analysis: dict[str, Any], *, gene_symbol: str) -> str:
    """Brief dataset introduction used when explorer figures lead the layout."""
    label = analysis.get("dataset_label") or "Allen dataset"
    assay = analysis.get("assay_terminology") or "single-cell RNA sequencing"
    summary = dict(analysis.get("summary") or {})
    noun = str(analysis.get("count_noun") or "celltype")
    noun_plural = (
        "transcriptomic cell types" if noun == "celltype" else "transcriptomic clusters"
    )
    valid = summary.get(f"valid_{noun}_count")
    nonzero = summary.get(f"nonzero_{noun}_count")
    _ = gene_symbol
    return (
        f"The Allen {label} dataset profiles {assay} across {noun_plural}. "
        f"Within this dataset, the gene showed nonzero aggregate expression across "
        f"{nonzero} of {valid} {noun_plural} "
        f"({_pct(summary.get('nonzero_percentage'))}% of {noun_plural} with a valid "
        f"trimmed-mean value)."
    )


def _allen_scatter_interpretation(analysis: dict[str, Any], *, gene_symbol: str) -> str:
    summary = dict(analysis.get("summary") or {})
    top = list(summary.get("top") or [])
    label = analysis.get("dataset_label") or "Allen dataset"
    noun = str(analysis.get("count_noun") or "celltype")
    entity_noun = "cell types" if noun == "celltype" else "clusters"
    if not top:
        return (
            f"The paired cell-type and gene-expression scatter for {gene_symbol} "
            f"in the Allen {label} explorer places the gene in the sampled taxonomy."
        )
    leader = top[0]
    ancestors = list(summary.get("top_taxonomy_ancestors") or [])
    parts = [
        f"In the paired scatter view, the highest aggregate trimmed mean was "
        f"{format_value(leader.get('value'))} in {leader.get('label')}, against a "
        f"median of {format_value(summary.get('median'))} "
        f"(maximum-to-median ratio "
        f"{format_value(summary.get('max_to_median_ratio'), digits=2)})."
    ]
    ancestor_sentence = _overlapping_taxonomy_ancestors_sentence(
        ancestors, top_count=len(top), entity_noun=entity_noun
    )
    if ancestor_sentence:
        parts.append(ancestor_sentence)
    return " ".join(parts)


def _allen_heatmap_interpretation(analysis: dict[str, Any], *, gene_symbol: str) -> str:
    summary = dict(analysis.get("summary") or {})
    top = list(summary.get("top") or [])[:5]
    label = analysis.get("dataset_label") or "Allen dataset"
    if not top:
        return (
            f"The complete heatmap for {gene_symbol} in the Allen {label} explorer "
            f"shows the gene among comparison markers across the taxonomy."
        )
    leaders = ", ".join(
        f"{row.get('label')} ({format_value(row.get('value'))})" for row in top[:3]
    )
    return (
        f"The complete heatmap places {gene_symbol} in dendrogram context with "
        f"comparison genes. The leading trimmed-mean populations include {leaders}."
    )


def _branch_presentation(
    *,
    analysis: dict[str, Any],
    gene_symbol: str,
    label_column: str,
    count_noun: str,
    explorer_symbol: str,
    figure_notes: Sequence[str],
    figure_statuses: dict[str, str],
    top_n: int,
) -> dict[str, Any]:
    """Presentation-ready view of one Allen branch, success or not."""
    dataset = str(analysis.get("dataset") or "")
    summary = dict(analysis.get("summary") or {})
    ok = bool(analysis.get("ok"))
    figures_ok = all(
        str(figure_statuses.get(name) or "") == STATUS_SUCCESS
        for name in (ac.VISUALIZATION_SCATTER, ac.VISUALIZATION_HEATMAP)
    )
    source_page_url: str | None = None
    explorer_url: str | None = None
    if dataset in ac.DATASET_SOURCE_PAGES:
        try:
            source_page_url = ac.dataset_source_page(dataset)
        except KeyError:
            source_page_url = None
    if explorer_symbol and dataset in ac.EXPLORER_DATASET_SLUGS:
        try:
            explorer_url = ac.figure_url(
                dataset=dataset,
                gene_symbol=explorer_symbol,
                visualization=ac.VISUALIZATION_SCATTER,
            )
        except (KeyError, ValueError):
            explorer_url = None
    label = analysis.get("dataset_label") or dataset
    # Evidence keeps the full top_n; polished fallback uses at most five rows.
    evidence_top = list(summary.get("top") or [])[:top_n]
    display_top = [] if figures_ok else evidence_top[:FALLBACK_TOP_N]
    if ok and figures_ok:
        narrative = _allen_short_intro(analysis, gene_symbol=gene_symbol)
        scatter_interp = _allen_scatter_interpretation(
            analysis, gene_symbol=gene_symbol
        )
        heatmap_interp = _allen_heatmap_interpretation(
            analysis, gene_symbol=gene_symbol
        )
    elif ok:
        narrative = allen_branch_narrative(analysis, gene_symbol=gene_symbol)
        scatter_interp = None
        heatmap_interp = None
    else:
        narrative = allen_unavailable_note(dataset=dataset, analysis=analysis)
        scatter_interp = None
        heatmap_interp = None
    source_link_label = None
    if source_page_url:
        if dataset == ac.DATASET_HUMAN_M1:
            source_link_label = "Allen Brain Map: Human M1 10x"
        elif dataset == ac.DATASET_MOUSE_CTX_HPF:
            source_link_label = (
                "Allen Brain Map: Mouse Whole Cortex and Hippocampus 10x"
            )
        else:
            source_link_label = f"Allen Brain Map: {label}"
    explorer_link_label = (
        f"Open {explorer_symbol} in Transcriptomics Explorer"
        if explorer_url and explorer_symbol
        else None
    )
    return {
        "dataset": dataset,
        "dataset_label": label,
        "assay_terminology": analysis.get("assay_terminology"),
        "sampling_scope": analysis.get("sampling_scope"),
        "analysis_status": STATUS_SUCCESS if ok else STATUS_SOURCE_UNAVAILABLE,
        "figures_complete": figures_ok,
        "narrative": narrative,
        "scatter_interpretation": scatter_interp,
        "heatmap_interpretation": heatmap_interp,
        "status_note": (
            None if ok else allen_unavailable_note(dataset=dataset, analysis=analysis)
        ),
        "figure_status_notes": list(figure_notes),
        "figure_statuses": dict(figure_statuses),
        "value_label": TRIMMED_MEAN_VALUE_LABEL,
        "label_column_header": label_column,
        "source_page_url": source_page_url,
        "source_link_label": source_link_label,
        "explorer_url": explorer_url,
        "explorer_link_label": explorer_link_label,
        # Compatibility aliases used by older presentation helpers.
        "database_url": source_page_url,
        "database_link_label": source_link_label,
        "requested_symbol": explorer_symbol,
        "source_symbol": summary.get("source_symbol"),
        f"valid_{count_noun}_count": summary.get(f"valid_{count_noun}_count"),
        f"nonzero_{count_noun}_count": summary.get(f"nonzero_{count_noun}_count"),
        "nonzero_percentage": summary.get("nonzero_percentage"),
        "maximum": summary.get("maximum"),
        "median": summary.get("median"),
        "max_to_median_ratio": summary.get("max_to_median_ratio"),
        "taxonomy_reconciliation": analysis.get("taxonomy_reconciliation"),
        "top": [
            {
                "label": row.get("label"),
                "value": row.get("value"),
                "value_display": format_value(row.get("value")),
                "ancestors": row.get("ancestors"),
            }
            for row in display_top
        ],
        "top_evidence": [
            {
                "label": row.get("label"),
                "value": row.get("value"),
                "value_display": format_value(row.get("value")),
                "ancestors": row.get("ancestors"),
            }
            for row in evidence_top
        ],
    }


def _branch_evidence_model(
    *,
    analysis: dict[str, Any],
    block: dict[str, Any],
    count_noun: str,
    figures: dict[str, dict[str, Any] | None],
    figure_statuses: dict[str, str],
    source_checksums: dict[str, str | None],
    artifact_relative_path: str | None,
) -> dict[str, Any]:
    """Compact per-branch evidence entry (no raw matrices)."""
    return {
        "dataset": block.get("dataset"),
        "dataset_label": block.get("dataset_label"),
        "assay_terminology": block.get("assay_terminology"),
        "sampling_scope": block.get("sampling_scope"),
        "analysis_status": block.get("analysis_status"),
        "requested_symbol": block.get("requested_symbol"),
        "source_symbol": block.get("source_symbol"),
        "target_row_match_count": analysis.get("match_count"),
        "summary": {
            f"valid_{count_noun}_count": block.get(f"valid_{count_noun}_count"),
            f"nonzero_{count_noun}_count": block.get(f"nonzero_{count_noun}_count"),
            "nonzero_percentage": block.get("nonzero_percentage"),
            "maximum": block.get("maximum"),
            "median": block.get("median"),
            "max_to_median_ratio": block.get("max_to_median_ratio"),
            "top": block.get("top_evidence") or block.get("top"),
        },
        "taxonomy_reconciliation": analysis.get("taxonomy_reconciliation"),
        "figure_statuses": dict(figure_statuses),
        "figures": {name: value for name, value in figures.items()},
        "source_checksums": dict(source_checksums),
        "artifacts": {"celltype_summary_json": artifact_relative_path},
        "localization_summary": block.get("narrative"),
    }


def _capture_allen_figure(
    *,
    dataset: str,
    gene_symbol: str,
    visualization: str,
    dossier_run_id: str,
    report_gene_symbol: str,
    item_key: str,
    section_cfg: Section2cConfig,
    settings: Settings,
    persist_db: bool,
    api_runs: list,
    raw_meta: list,
    evidence: list,
    audit: dict[str, Any],
    parent_raw_artifact_ids: Sequence[str],
) -> tuple[bool, dict[str, Any] | None, str]:
    """Attempt one optional Allen explorer capture. Failure is reported, not raised."""
    species = "human" if dataset == ac.DATASET_HUMAN_M1 else "mouse"
    label = ac.DATASET_LABELS.get(dataset, dataset)
    audit_key = f"{species}_{visualization}"
    try:
        result = ac.capture_explorer_figure(
            dataset=dataset,
            gene_symbol=gene_symbol,
            visualization=visualization,
            max_attempts=section_cfg.figure_max_attempts,
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning("Allen %s %s capture raised: %s", label, visualization, exc)
        audit["figures"][audit_key] = {
            "status": STATUS_SOURCE_UNAVAILABLE,
            "error_type": "capture_exception",
            "error_message": str(exc)[:400],
        }
        return (
            False,
            None,
            f"The Allen {label} {visualization} figure is unavailable for this run "
            f"(capture_exception).",
        )

    # A browser navigation is a real network request, so it gets one ApiRun.
    api = ApiRun(
        dossier_run_id=dossier_run_id,
        gene_symbol=report_gene_symbol,
        source_name=ac.SOURCE_NAME,
        endpoint_name="capture_allen_celltypes_explorer_figure",
        request_url=str(result.get("source_page_url") or ""),
        request_params={
            "dataset": dataset,
            "visualization": visualization,
            "requested_gene": result.get("requested_gene"),
            "retrieval_method": ac.CAPTURE_RETRIEVAL_METHOD,
            "browser_channel": result.get("browser_channel"),
            "allowlisted_hosts": sorted(ac.ALLOWED_FIGURE_HOSTS),
        },
        success=bool(result.get("ok")),
        error_type=result.get("error_type"),
        error_message=result.get("error_message"),
    )
    audit["figures"][audit_key] = {
        "status": result.get("status"),
        "source_page_url": result.get("source_page_url"),
        "outer_page_url": result.get("outer_page_url") or result.get("final_url"),
        "final_url": result.get("final_url"),
        "viewer_frame_url": result.get("viewer_frame_url"),
        "dom_selector": result.get("dom_selector"),
        "capture_kind": result.get("capture_kind"),
        "candidate_count": result.get("candidate_count"),
        "accepted_candidate_score": result.get("accepted_candidate_score"),
        "plot_panel_count": result.get("plot_panel_count"),
        "includes_celltype_panel": result.get("includes_celltype_panel"),
        "includes_gene_expression_panel": result.get("includes_gene_expression_panel"),
        "includes_taxonomy_dendrogram": result.get("includes_taxonomy_dendrogram"),
        "validation": result.get("validation"),
        "browser_channel": result.get("browser_channel"),
        "launch_attempts": result.get("launch_attempts"),
        "attempts": result.get("attempts"),
        "error_type": result.get("error_type"),
        "error_message": result.get("error_message"),
        "api_run_id": api.id,
    }

    png = result.get("png")
    if not result.get("ok") or not png:
        _save_api_run_failure(api, persist_db=persist_db)
        api_runs.append(api)
        reason = result.get("error_type") or STATUS_SOURCE_UNAVAILABLE
        return (
            False,
            None,
            f"The Allen {label} {visualization} figure is unavailable for this run "
            f"({reason}).",
        )

    api_runs.append(api)
    try:
        artifact, meta = _persist_artifact_bytes(
            dossier_run_id=dossier_run_id,
            source_name=ac.SOURCE_NAME,
            content=png,
            extension="png",
            artifact_type="png",
            filename_hint=f"allen-{species}-{visualization}-{slugify(gene_symbol)}",
            settings=settings,
            api_run=api,
            persist_db=persist_db,
            notes={
                "artifact_class": "derived_capture",
                "artifact_origin": "allen_celltypes_explorer",
                "artifact_role": f"section_2c_{species}_{visualization}_figure",
                "retrieval_method": ac.CAPTURE_RETRIEVAL_METHOD,
                "dataset": dataset,
                "visualization": visualization,
                "source_page_url": result.get("source_page_url"),
                "final_url": result.get("final_url"),
                "dom_selector": result.get("dom_selector"),
                "requested_gene": result.get("requested_gene"),
                "parent_raw_artifact_ids": [str(i) for i in parent_raw_artifact_ids if i],
                "presentation_item_key": f"{item_key}-{species}-{visualization}",
            },
            validate=_validate_nonblank_image,
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning("Allen %s %s capture failed validation: %s", label, visualization, exc)
        audit["figures"][audit_key]["persist_error"] = str(exc)[:400]
        return (
            False,
            None,
            f"The Allen {label} {visualization} figure did not pass image validation "
            f"for this run.",
        )

    figure_value = {
        "dataset": dataset,
        "dataset_label": label,
        "visualization": visualization,
        "requested_gene": result.get("requested_gene"),
        "source_page_url": result.get("source_page_url"),
        "final_url": result.get("final_url"),
        "dom_selector": result.get("dom_selector"),
        "relative_path": meta.get("relative_path"),
        "sha256": meta.get("content_hash") or artifact.content_hash,
        "width": meta.get("width"),
        "height": meta.get("height"),
        "byte_size": meta.get("byte_size"),
        "media_type": "image/png",
        "artifact_class": "derived_capture",
        "retrieval_method": ac.CAPTURE_RETRIEVAL_METHOD,
        "presentation_item_key": f"{item_key}-{species}-{visualization}",
        "figure_raw_artifact_id": artifact.id,
        "figure_api_run_id": api.id,
    }
    rec = _evidence(
        dossier_run_id=dossier_run_id,
        gene_symbol=report_gene_symbol,
        source_name=ac.SOURCE_NAME,
        fact_type="allen_celltype_explorer_figure",
        key=f"{dataset}-{visualization}",
        value=figure_value,
        display_text=(
            f"Allen {label} {visualization} visualization for "
            f"{result.get('requested_gene')}."
        ),
        api_run_id=api.id,
        raw_artifact_id=artifact.id,
        organism="Homo sapiens" if species == "human" else "Mus musculus",
        taxon_id=9606 if species == "human" else 10090,
    )
    _append_evidence(evidence, rec, persist_db=persist_db)
    audit["figures"][audit_key]["relative_path"] = meta.get("relative_path")
    audit["figures"][audit_key]["sha256"] = meta.get("content_hash")
    return True, figure_value, ""


__all__ = [
    "CALCULATION_VERSION",
    "CHART_VERSION",
    "CONFIDENCE_INTERVAL_STATUS",
    "DROPVIZ_ASSAY_TERMINOLOGY",
    "DROPVIZ_POPULATION_IDENTIFIER_NOTE",
    "DROPVIZ_REGIONAL_TSNE_LIMITATION_NOTE",
    "DROPVIZ_SAMPLING_SCOPE",
    "GEO_SERIES_URL",
    "GEO_SOURCE_KEY",
    "SECTION_2C_INTRO_TEXT",
    "SECTION_EXPRESSION",
    "STATUS_FIGURE_SUPPRESSED",
    "STATUS_NORMALIZATION_UNRESOLVED",
    "STATUS_NOT_ATTEMPTED",
    "STATUS_NOT_AVAILABLE",
    "STATUS_SOURCE_UNAVAILABLE",
    "STATUS_SUCCESS",
    "SUBSECTION_2C",
    "THERAPEUTIC_CAVEAT",
    "TRIMMED_MEAN_VALUE_LABEL",
    "Section2cConfig",
    "allen_branch_narrative",
    "allen_unavailable_note",
    "analyze_allen_dataset",
    "analyze_dropviz_matrix",
    "build_therapeutic_narrative",
    "dropviz_narrative",
    "dropviz_unavailable_note",
    "format_value",
    "node_generate_section_2c_derived_artifacts",
    "SECTION_2C_VISUAL_COMPLETE_ROLES",
    "accept_visual_complete_gene_report",
    "evaluate_section_2c_visual_complete",
    "scientific_core_accepted",
]
