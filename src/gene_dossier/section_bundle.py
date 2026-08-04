"""Section-scoped dossier generation for curator review (Sections 1a–1e).

Builds a standalone Section 1 document without LLM synthesis or full-report
rendering. Provenance IDs live only in the audit JSON; polished outputs use
deterministic opaque evidence references.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from gene_dossier.config import Settings, get_settings
from gene_dossier.db import get_dossier_run, init_db, save_dossier_run, session_scope
from gene_dossier.models import DossierRun, EvidenceRecord, utcnow
from gene_dossier.rancho_report import (
    SECTION_1C_PDF_PAGE_BREAK,
    _escape,
    _rancho_css,
    _render_subsection,
    rasterize_pdf_pages_to_pngs,
    render_rancho_pdf,
    render_section_1c_subsection_segments,
    sanitize_polished_citation_tokens,
)
from gene_dossier.report_presentation import build_section_presentation
from gene_dossier.report_schema import (
    REPORT_SECTIONS,
    REPORT_STYLE,
    ReportContentBlock,
    ReportCover,
    ReportDocument,
    ReportMajorSection,
    ReportSubsection,
    infer_chromosome,
)
from gene_dossier.section_1c import node_generate_section_1c_derived_artifacts
from gene_dossier.section_1d import (
    evaluate_section_1d_reference_genes_acceptance,
    node_generate_section_1d_derived_artifacts,
)
from gene_dossier.section_1e import (
    Section1eConfig,
    node_generate_section_1e_derived_artifacts,
)
from gene_dossier.section_2a import (
    Section2aConfig,
    node_generate_section_2a_derived_artifacts,
)
from gene_dossier.section_2b import (
    Section2bConfig,
    node_generate_section_2b_derived_artifacts,
)
from gene_dossier.ucsc_figure import redact_api_key
from gene_dossier.workflow import (
    DossierState,
    WorkflowTransientContext,
    coverage_updates_from_state,
    node_call_source_clients,
    node_normalize_evidence,
    node_resolve_gene_identity,
    node_save_raw_artifacts,
)

logger = logging.getLogger(__name__)

SUPPORTED_SECTION_BUNDLE_KEYS = ("1a", "1b", "1c", "1d", "1e", "2a", "2b")
DEFAULT_SECTION_BUNDLE_KEYS = ("1a", "1b")

SECTION_SOURCE_DEPENDENCIES: dict[str, set[str]] = {
    "1a": set(),
    "1b": {"UCSC"},
    "1c": {"CDD", "PDBe"},
    "1d": {"AlphaFold"},
    "1e": {"NCBI Datasets", "OrthoDB"},
    "2a": {"GTEx"},
    "2b": {"Allen Brain Atlas", "BrainRNASeq"},
}

_OPAQUE_REF_BY_ROLE = {
    ("1a", "gene_aliases_table"): "ev-1a-gene-aliases-table",
    ("1b", "ucsc_conservation_figure"): "ev-1b-conservation-figure",
    ("1c", "section_1c_cdd_link"): "ev-1c-cdd-link",
    ("1c", "section_1c_domain_table"): "ev-1c-cdd-domain-table",
    ("1c", "section_1c_domain_architecture_figure"): "ev-1c-cdd-architecture",
    ("1c", "section_1c_pdb_table"): "ev-1c-pdbe-candidate-table",
    ("1c", "section_1c_pdb_assembly_figure"): "ev-1c-pdbe-assembly-figure",
    ("1c", "section_1c_pdb_domain_focus_figure"): "ev-1c-pdbe-domain-focus-figure",
    ("1c", "section_1c_image_attribution"): "ev-1c-image-attribution",
    ("1e", "section_1e_narrative"): "ev-1e-ortholog-summary",
    ("1e", "section_1e_fallback_table"): "ev-1e-ortholog-fallback-table",
    ("1e", "section_1e_attribution"): "ev-1e-ortholog-attribution",
    ("2a", "section_2a_gtex_intro"): "ev-2a-gtex-introduction",
    ("2a", "section_2a_gtex_all_tissues_link"): "ev-2a-gtex-all-tissues-summary",
    ("2a", "section_2a_gtex_all_tissues_figure"): "ev-2a-gtex-all-tissues-figure",
    ("2a", "section_2a_gtex_brain_link"): "ev-2a-gtex-brain-summary",
    ("2a", "section_2a_gtex_brain_figure"): "ev-2a-gtex-brain-figure",
    ("2a", "section_2a_hbt_intro"): "ev-2a-hbt-summary",
    ("2a", "section_2a_hbt_link"): "ev-2a-hbt-summary",
    ("2a", "section_2a_hbt_figure"): "ev-2a-hbt-figure",
    ("2b", "section_2b_intro"): "ev-2b-introduction",
    ("2b", "section_2b_summary_table"): "ev-2b-summary-table",
    ("2b", "section_2b_category_status"): "ev-2b-category-status",
    ("2b", "section_2b_celltype_intro"): "ev-2b-celltype-introduction",
    ("2b", "section_2b_source_link"): "ev-2b-source-link",
    ("2b", "section_2b_celltype_figure"): "ev-2b-celltype-figure",
}

_SECTION_1C_REF_SUFFIX_BY_ROLE = {
    "section_1c_domain_summary": "summary",
    "section_1c_domain_thumbnail": "thumbnail",
    "section_1c_feature_summary": "summary",
    "section_1c_feature_thumbnail": "thumbnail",
    "section_1c_pdb_link": "summary",
    "section_1c_pdb_official_image": "official-image",
}

_SECTION_1D_REF_SUFFIX_BY_ROLE = {
    "section_1d_species_link": "summary",
    "section_1d_human_structure_capture": "viewer-capture",
}

_SECTION_1E_REF_SUFFIX_BY_ROLE = {
    "section_1e_ortholog_capture": "viewer-capture",
}

# Status lines and the deterministic report-side legend are visible but are not
# source-backed structural evidence.
_SECTION_1D_NON_EVIDENCE_ROLES = frozenset(
    {"section_1d_species_status", "section_1d_confidence_legend"}
)
_SECTION_2A_NON_EVIDENCE_ROLES = frozenset({"section_2a_source_status"})
_SECTION_2B_NON_EVIDENCE_ROLES = frozenset({"section_2b_source_status"})

_SAFE_ITEM_KEY_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")

_RAW_ID_KEYS = frozenset(
    {
        "source_ids",
        "evidence_record_ids",
        "api_run_ids",
        "raw_artifact_ids",
        "source_id",
        "evidence_record_id",
        "api_run_id",
        "raw_artifact_id",
        "figure_api_run_id",
        "figure_raw_artifact_id",
    }
)


class SectionBundleError(ValueError):
    """Invalid section-bundle request or failed generation."""


@dataclass
class SectionBundleResult:
    """Outputs from a successful section-bundle run."""

    gene_symbol: str
    dossier_run_id: str
    selected_section_keys: list[str]
    output_dir: Path
    output_paths: dict[str, Path] = field(default_factory=dict)
    status: str = "completed"
    errors: list[str] = field(default_factory=list)


def validate_section_keys(section_keys: Iterable[str] | None) -> list[str]:
    """Accept only supported section keys; dedupe; canonicalize order."""
    raw = [str(k).strip().lower() for k in (section_keys or [])]
    if not raw:
        raise SectionBundleError(
            "At least one section key is required "
            f"({', '.join(SUPPORTED_SECTION_BUNDLE_KEYS)})."
        )
    normalized: list[str] = []
    for key in raw:
        if key in {"1.a", "a"}:
            key = "1a"
        elif key in {"1.b", "b"}:
            key = "1b"
        elif key in {"1.c", "c"}:
            key = "1c"
        elif key in {"1.d", "d"}:
            key = "1d"
        elif key in {"1.e", "e"}:
            key = "1e"
        elif key in {"2.a"}:
            key = "2a"
        elif key in {"2.b"}:
            key = "2b"
        if key not in SUPPORTED_SECTION_BUNDLE_KEYS:
            raise SectionBundleError(
                f"Unsupported section key {key!r}. Supported: "
                f"{', '.join(SUPPORTED_SECTION_BUNDLE_KEYS)}"
            )
        if key not in normalized:
            normalized.append(key)
    return [k for k in SUPPORTED_SECTION_BUNDLE_KEYS if k in normalized]


def sources_for_sections(section_keys: Sequence[str]) -> list[str]:
    """Non-identity sources required by the selected sections.

    Section 1d declares AlphaFold as a dependency for accounting, but the
    generic human-only AlphaFold client must **not** be invoked when 1d is
    selected — ``section_1d`` owns those network requests exclusively.

    Section 1e likewise owns NCBI Datasets ortholog/taxonomy and OrthoDB
    network requests; the generic Datasets ortholog client must not run when
    1e is selected.

    Section 2a owns GTEx (and HBT) network requests; the generic GTEx client
    must not run when 2a is selected.

    Section 2b owns Allen Brain Atlas and BrainRNASeq network requests; the
    generic clients must not run when 2b is selected.
    """
    needed: set[str] = set()
    for key in section_keys:
        needed |= SECTION_SOURCE_DEPENDENCIES.get(key, set())
    if "1d" in section_keys:
        needed.discard("AlphaFold")
    if "1e" in section_keys:
        needed.discard("NCBI Datasets")
        needed.discard("OrthoDB")
    if "2a" in section_keys:
        needed.discard("GTEx")
        needed.discard("HBT")
        needed.discard("Human Brain Transcriptome")
    if "2b" in section_keys:
        needed.discard("Allen Brain Atlas")
        needed.discard("BrainRNASeq")
        needed.discard("Allen Brain")
        needed.discard("Barres Lab")
    return sorted(needed)


def validate_safe_item_key(value: str | None) -> str:
    """Validate a safe polished-output grouping key."""
    key = (value or "").strip()
    if not key or not _SAFE_ITEM_KEY_RE.fullmatch(key):
        raise SectionBundleError(f"Unsafe or empty presentation_item_key: {value!r}")
    return key


def section_1c_evidence_ref(block: ReportContentBlock) -> str:
    """Dynamic Section 1c opaque ref from safe item key + role suffix."""
    role = str(block.presentation_role or "")
    suffix = _SECTION_1C_REF_SUFFIX_BY_ROLE.get(role)
    if not suffix:
        raise SectionBundleError(f"Section 1c role is not dynamically referenceable: {role!r}")
    key = validate_safe_item_key(block.presentation_item_key)
    return f"ev-1c-{key}-{suffix}"


def section_1d_evidence_ref(block: ReportContentBlock) -> str:
    """Dynamic Section 1d opaque ref from safe item key + role suffix."""
    role = str(block.presentation_role or "")
    suffix = _SECTION_1D_REF_SUFFIX_BY_ROLE.get(role)
    if not suffix:
        raise SectionBundleError(f"Section 1d role is not dynamically referenceable: {role!r}")
    key = validate_safe_item_key(block.presentation_item_key)
    return f"ev-1d-{key}-{suffix}"


def section_1e_evidence_ref(block: ReportContentBlock) -> str:
    """Dynamic Section 1e opaque ref from safe item key + role suffix."""
    role = str(block.presentation_role or "")
    suffix = _SECTION_1E_REF_SUFFIX_BY_ROLE.get(role)
    if not suffix:
        raise SectionBundleError(f"Section 1e role is not dynamically referenceable: {role!r}")
    key = validate_safe_item_key(block.presentation_item_key)
    return f"ev-1e-{key}-{suffix}"


def opaque_evidence_ref(section_key: str, block: ReportContentBlock, *, index: int) -> str:
    """Deterministic opaque ref from section + presentation role / kind."""
    role = block.presentation_role
    if section_key == "1c" and role in _SECTION_1C_REF_SUFFIX_BY_ROLE:
        return section_1c_evidence_ref(block)
    if section_key == "1d" and role in _SECTION_1D_REF_SUFFIX_BY_ROLE:
        return section_1d_evidence_ref(block)
    if section_key == "1e" and role in _SECTION_1E_REF_SUFFIX_BY_ROLE:
        return section_1e_evidence_ref(block)
    mapped = _OPAQUE_REF_BY_ROLE.get((section_key, role)) if role else None
    if mapped:
        return mapped
    if section_key == "1b" and block.kind == "link":
        return "ev-1b-transcript-link"
    if section_key == "1b" and block.kind == "narrative" and index == 0:
        return "ev-1b-summary"
    if section_key == "1a" and block.kind == "table":
        return "ev-1a-gene-aliases-table"
    if section_key == "1c" and block.kind == "narrative" and index == 0:
        return "ev-1c-introduction"
    return f"ev-{section_key}-block-{index + 1}"


def sanitize_credentials(value: Any) -> Any:
    """Recursively redact API-key material from nested JSON-compatible values."""
    if isinstance(value, str):
        return redact_api_key(value)
    if isinstance(value, dict):
        out: dict[str, Any] = {}
        for key, item in value.items():
            key_l = str(key).lower()
            if key_l in {"apikey", "api_key", "ucsc_browser_api_key"}:
                out[key] = "REDACTED"
            else:
                out[key] = sanitize_credentials(item)
        return out
    if isinstance(value, list):
        return [sanitize_credentials(item) for item in value]
    if isinstance(value, tuple):
        return [sanitize_credentials(item) for item in value]
    return value


def sanitize_polished_text(text: str | None) -> str:
    """Credential redaction plus removal of polished ``[source_id=...]`` tokens."""
    return sanitize_polished_citation_tokens(sanitize_credentials(text or "") or "")


def sanitize_secrets(value: Any) -> Any:
    """Backward-compatible alias for credential-only sanitization."""
    return sanitize_credentials(value)


def _strip_raw_ids(payload: Any) -> Any:
    if isinstance(payload, dict):
        return {
            key: _strip_raw_ids(value)
            for key, value in payload.items()
            if key not in _RAW_ID_KEYS
        }
    if isinstance(payload, list):
        return [_strip_raw_ids(item) for item in payload]
    return payload


def build_provenance_index(
    *,
    evidence_records: Sequence[EvidenceRecord],
    api_runs: Sequence[Any] | None = None,
    raw_artifacts: Sequence[Any] | None = None,
) -> dict[str, Any]:
    """Index evidence / API-run / artifact IDs for opaque evidence-ref mapping."""
    by_evidence_id: dict[str, EvidenceRecord] = {
        rec.id: rec for rec in evidence_records if rec.id
    }
    artifact_by_id: dict[str, Any] = {}
    artifact_by_path: dict[str, str] = {}
    for meta in raw_artifacts or []:
        if isinstance(meta, dict):
            aid = meta.get("id")
            path = meta.get("file_path") or meta.get("relative_path")
            if aid:
                artifact_by_id[str(aid)] = meta
            if aid and path:
                artifact_by_path[str(path)] = str(aid)
        else:
            aid = getattr(meta, "id", None)
            path = getattr(meta, "file_path", None)
            if aid:
                artifact_by_id[str(aid)] = meta
            if aid and path:
                artifact_by_path[str(path)] = str(aid)

    api_by_id: dict[str, Any] = {}
    for api in api_runs or []:
        if isinstance(api, dict):
            aid = api.get("id")
            if aid:
                api_by_id[str(aid)] = api
        else:
            aid = getattr(api, "id", None)
            if aid:
                api_by_id[str(aid)] = api

    return {
        "by_evidence_id": by_evidence_id,
        "artifact_by_id": artifact_by_id,
        "artifact_by_path": artifact_by_path,
        "api_by_id": api_by_id,
    }


def _ids_for_block(
    block: ReportContentBlock,
    *,
    provenance: dict[str, Any] | None,
) -> dict[str, list[str]]:
    evidence_ids = list(dict.fromkeys(block.evidence_record_ids or []))
    source_ids = list(dict.fromkeys(block.source_ids or []))
    artifact_ids: list[str] = []
    api_run_ids: list[str] = []
    if provenance:
        by_ev = provenance.get("by_evidence_id") or {}
        by_path = provenance.get("artifact_by_path") or {}
        for eid in evidence_ids:
            rec = by_ev.get(eid)
            if rec is None:
                continue
            if rec.raw_artifact_id:
                artifact_ids.append(rec.raw_artifact_id)
            if rec.api_run_id:
                api_run_ids.append(rec.api_run_id)
            if isinstance(rec.value, dict):
                for key in ("figure_raw_artifact_id", "raw_artifact_id"):
                    val = rec.value.get(key)
                    if val:
                        artifact_ids.append(str(val))
                for key in ("figure_api_run_id", "api_run_id"):
                    val = rec.value.get(key)
                    if val:
                        api_run_ids.append(str(val))
                rel = rec.value.get("relative_path") or rec.value.get(
                    "local_artifact_path"
                )
                if rel and str(rel) in by_path:
                    artifact_ids.append(by_path[str(rel)])
        if block.figure_path and str(block.figure_path) in by_path:
            artifact_ids.append(by_path[str(block.figure_path)])
    return {
        "evidence_record_ids": evidence_ids,
        "source_ids": source_ids,
        "raw_artifact_ids": list(dict.fromkeys(artifact_ids)),
        "api_run_ids": list(dict.fromkeys(api_run_ids)),
    }


def create_section_bundle_run(
    *,
    gene_symbol: str,
    selected_section_keys: Sequence[str],
    settings: Settings | None = None,
    persist_db: bool = True,
    dossier_run_id: str | None = None,
) -> tuple[DossierRun, DossierState]:
    """Create a ``section_bundle`` DossierRun and a complete compatible state."""
    gene = gene_symbol.strip()
    if not gene:
        raise SectionBundleError("gene_symbol is required")
    keys = validate_section_keys(selected_section_keys)
    _ = settings or get_settings()

    run = DossierRun(
        gene_symbol=gene,
        run_type="section_bundle",
        status="running",
        config={"selected_section_keys": list(keys)},
        notes="section_scoped_generation",
    )
    if dossier_run_id:
        run.id = dossier_run_id

    if persist_db:
        init_db()
        with session_scope() as session:
            save_dossier_run(session, run)

    state: DossierState = {
        "gene_symbol": gene,
        "dossier_run_id": run.id,
        "run_type": "section_bundle",
        "selected_section_keys": list(keys),
        "gene_ids": {},
        "tool_results": [],
        "api_runs": [],
        "raw_artifacts": [],
        "evidence_records": [],
        "coverage": [],
        "sections": [],
        "claims": [],
        "verification_results": [],
        "synthesis_notes": [],
        "output_paths": {},
        "errors": [],
        "status": "running",
    }
    return run, state


def finalize_section_bundle_run(
    *,
    dossier_run_id: str,
    status: str,
    selected_section_keys: Sequence[str],
    errors: Sequence[str] | None = None,
    persist_db: bool = True,
) -> None:
    """Mark the scoped run completed/failed without rewriting run_type."""
    if not persist_db:
        return
    init_db()
    with session_scope() as session:
        existing = get_dossier_run(session, dossier_run_id)
        if existing is None:
            raise SectionBundleError(
                f"Section-bundle DossierRun not found: {dossier_run_id}"
            )
        if existing.run_type != "section_bundle":
            raise SectionBundleError(
                f"Refusing to finalize non-section-bundle run: {dossier_run_id}"
            )
        existing.status = status
        if status in {"completed", "failed"}:
            existing.completed_at = utcnow()
        cfg = dict(existing.config or {})
        cfg["selected_section_keys"] = list(selected_section_keys)
        existing.config = cfg
        if errors:
            note = "; ".join(
                str(item) for item in sanitize_credentials(list(errors))
            )
            existing.notes = (
                f"{existing.notes}; {note}" if existing.notes else note
            )
        save_dossier_run(session, existing)


def _major_section_spec(number: int) -> Any:
    return next(spec for spec in REPORT_SECTIONS if spec.number == number)


def _major_section_1_spec() -> Any:
    return _major_section_spec(1)


def assign_opaque_refs(
    *,
    section_key: str,
    blocks: Sequence[ReportContentBlock],
    provenance: dict[str, Any] | None = None,
) -> tuple[list[ReportContentBlock], dict[str, dict[str, list[str]]]]:
    """Attach deterministic evidence_ref and build audit map entries.

    Blocks with ``presentation_role == section_1d_species_status`` remain visible
    but receive ``evidence_ref=None`` and are omitted from ``evidence_reference_map``.
    """
    polished: list[ReportContentBlock] = []
    ref_map: dict[str, dict[str, list[str]]] = {}
    for index, block in enumerate(blocks):
        if block.presentation_role in _SECTION_1D_NON_EVIDENCE_ROLES:
            polished.append(block.model_copy(update={"evidence_ref": None}))
            continue
        if block.presentation_role in _SECTION_2A_NON_EVIDENCE_ROLES:
            polished.append(block.model_copy(update={"evidence_ref": None}))
            continue
        if block.presentation_role in _SECTION_2B_NON_EVIDENCE_ROLES:
            polished.append(block.model_copy(update={"evidence_ref": None}))
            continue
        base_ref = opaque_evidence_ref(section_key, block, index=index)
        # One item can carry several official images (e.g. two Cadherin_C
        # structure thumbnails), so repeats take a deterministic ordinal rather
        # than colliding.
        ref = base_ref
        ordinal = 2
        while ref in ref_map:
            ref = f"{base_ref}-{ordinal}"
            ordinal += 1
        polished_block = block.model_copy(update={"evidence_ref": ref})
        polished.append(polished_block)
        ref_map[ref] = _ids_for_block(block, provenance=provenance)
    if len(ref_map) != len(set(ref_map)):
        raise SectionBundleError("Duplicate evidence refs generated")
    return polished, ref_map


def serialize_presentation_block(block: ReportContentBlock) -> dict[str, Any]:
    """Serialize one polished block without raw IDs (figure path allowed)."""
    payload: dict[str, Any] = {
        "kind": block.kind,
        "title": block.title,
        "text": sanitize_polished_text(block.text) if block.text else None,
        "presentation_role": block.presentation_role,
        "presentation_item_key": block.presentation_item_key,
        "presentation_page_break_before": bool(block.presentation_page_break_before),
        "evidence_ref": block.evidence_ref,
    }
    if block.table_headers:
        payload["table_headers"] = list(block.table_headers)
        payload["table_rows"] = [list(row) for row in block.table_rows]
    if block.kind == "figure":
        payload["figure_path"] = block.figure_path
        payload["figure_alt"] = block.figure_caption or block.text
        payload["media_type"] = "image/png"
    if block.kind == "link" and block.links:
        payload["links"] = [
            {
                "label": sanitize_polished_text(link.get("label")),
                "url": redact_api_key(link.get("url") or ""),
            }
            for link in block.links
        ]
    return sanitize_credentials(_strip_raw_ids(payload))


def build_section_bundle_document(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    section_keys: Sequence[str],
    evidence_records: Sequence[EvidenceRecord],
    api_runs: Sequence[Any] | None = None,
    raw_artifacts: Sequence[Any] | None = None,
    section_status_by_key: dict[str, Any] | None = None,
) -> tuple[ReportDocument, dict[str, Any], dict[str, Any]]:
    """Build a section-bundle document (Major 1 and/or Major 2) plus audit."""
    keys = validate_section_keys(section_keys)
    provenance = build_provenance_index(
        evidence_records=evidence_records,
        api_runs=api_runs,
        raw_artifacts=raw_artifacts,
    )
    status_by_key = section_status_by_key or {}

    evidence_reference_map: dict[str, dict[str, list[str]]] = {}
    diagnostics: list[dict[str, Any]] = []
    figure_notes: list[str] = []

    keys_by_major: dict[int, list[str]] = {}
    for section_key in keys:
        major_num = int(section_key[0])
        keys_by_major.setdefault(major_num, []).append(section_key)

    report_majors: list[ReportMajorSection] = []
    presentation_majors: list[dict[str, Any]] = []

    for major_num in sorted(keys_by_major):
        major_spec = _major_section_spec(major_num)
        sub_by_key = {sub.key: sub for sub in major_spec.subsections}
        presentation_subsections: list[dict[str, Any]] = []
        report_subsections: list[ReportSubsection] = []

        for section_key in keys_by_major[major_num]:
            letter = section_key[-1]
            sub_spec = sub_by_key[letter]
            result = build_section_presentation(
                section_key=section_key,
                gene_symbol=gene_symbol,
                evidence_records=evidence_records,
                section_status=status_by_key.get(section_key),
            )
            for diag in result.diagnostics:
                diagnostics.append(
                    {
                        "section_key": section_key,
                        "field": diag.field,
                        "reason": diag.reason,
                        "severity": diag.severity,
                    }
                )
                if diag.field == "figure_note":
                    figure_notes.append(diag.reason)
            polished_blocks, ref_map = assign_opaque_refs(
                section_key=section_key,
                blocks=result.blocks,
                provenance=provenance,
            )
            evidence_reference_map.update(ref_map)

            figure_meta: dict[str, Any] = {}
            for rec in evidence_records:
                if rec.fact_type == "ucsc_conservation_figure" and isinstance(
                    rec.value, dict
                ):
                    figure_meta = dict(rec.value)
                    break

            serialized_blocks: list[dict[str, Any]] = []
            for block in polished_blocks:
                item = serialize_presentation_block(block)
                if block.presentation_role == "ucsc_conservation_figure" and figure_meta:
                    item["figure_path"] = (
                        figure_meta.get("relative_path")
                        or figure_meta.get("local_artifact_path")
                        or item.get("figure_path")
                    )
                    item["media_type"] = figure_meta.get("media_type") or "image/png"
                    item["width"] = figure_meta.get("width")
                    item["height"] = figure_meta.get("height")
                    item["sha256"] = figure_meta.get("sha256") or figure_meta.get(
                        "content_hash"
                    )
                    item["byte_size"] = figure_meta.get("byte_size")
                serialized_blocks.append(item)

            presentation_subsections.append(
                {
                    "key": section_key,
                    "title": sub_spec.title,
                    "blocks": serialized_blocks,
                }
            )
            report_subsections.append(
                ReportSubsection(
                    key=letter,
                    title=sub_spec.title,
                    toc_title=sub_spec.toc_title,
                    presentation_blocks=polished_blocks,
                    status="populated" if polished_blocks else "empty",
                )
            )

        report_majors.append(
            ReportMajorSection(
                number=major_num,
                key=str(major_num),
                title=major_spec.title,
                toc_title=major_spec.toc_title,
                subsections=report_subsections,
                status="populated" if report_subsections else "empty",
                narrative_markdown=None,
                synthesis_status=None,
            )
        )
        presentation_majors.append(
            {
                "number": major_num,
                "title": major_spec.title,
                "subsections": presentation_subsections,
            }
        )

    cover = ReportCover(
        gene_symbol=gene_symbol,
        chromosome=infer_chromosome(list(evidence_records)),
        curator="Gene Dossier Platform",
    )
    document = ReportDocument(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        cover=cover,
        sections=report_majors,
    )

    # Backward-compatible presentation shape: keep ``major_section`` for single
    # Major 1 runs; add ``major_sections`` whenever multiple majors are present.
    presentation_raw: dict[str, Any] = {
        "document_type": "section_bundle",
        "gene_symbol": gene_symbol,
        "dossier_run_id": dossier_run_id,
        "selected_section_keys": list(keys),
        "major_sections": presentation_majors,
    }
    if len(presentation_majors) == 1:
        presentation_raw["major_section"] = presentation_majors[0]
    elif presentation_majors:
        # Prefer Major 1 for the legacy key when present.
        major1 = next((m for m in presentation_majors if m["number"] == 1), None)
        presentation_raw["major_section"] = major1 or presentation_majors[0]

    presentation = sanitize_credentials(presentation_raw)

    def _polish_strings(node: Any) -> Any:
        if isinstance(node, str):
            return sanitize_polished_text(node)
        if isinstance(node, dict):
            return {k: _polish_strings(v) for k, v in node.items()}
        if isinstance(node, list):
            return [_polish_strings(v) for v in node]
        return node

    presentation = _polish_strings(presentation)

    audit = {
        "document_type": "section_bundle_audit",
        "gene_symbol": gene_symbol,
        "dossier_run_id": dossier_run_id,
        "run_type": "section_bundle",
        "selected_section_keys": list(keys),
        "evidence_reference_map": evidence_reference_map,
        "diagnostics": diagnostics,
        "figure_notes": figure_notes,
        "section_1d_status": status_by_key.get("1d"),
        "section_1e_status": status_by_key.get("1e"),
        "section_2a_status": status_by_key.get("2a"),
        "section_2b_status": status_by_key.get("2b"),
    }
    return document, presentation, audit


def render_section_bundle_html(
    document: ReportDocument,
    *,
    show_header_logos: bool = True,
    include_page_chrome: bool = True,
    include_major_heading: bool = True,
) -> str:
    """Render section-bundle HTML (Major 1 and/or Major 2; no cover/TOC)."""
    if not document.sections:
        raise SectionBundleError("Section bundle document has no sections")

    from gene_dossier.rancho_report import (
        render_section_2a_subsection_segments,
        render_section_2b_subsection_segments,
    )

    body_parts: list[str] = []
    page_break = SECTION_1C_PDF_PAGE_BREAK

    for major_index, major in enumerate(document.sections):
        if major_index > 0:
            body_parts.append(page_break)

        if major.number == 1:
            body_parts.extend(
                _render_major_section_1_pages(
                    major, include_major_heading=include_major_heading
                )
            )
            continue

        # Major 2: 2a segments first (unchanged three-page split), then 2b on a
        # clean page after all 2a pages. Assembled 2a+2b does not repeat Major 2
        # heading on 2b pages. Focused 2b still gets the Major 2 heading.
        heading = f"{major.number}. {major.title}"
        first_page_parts: list[str] = [
            (
                f'<section id="section-{major.number}" '
                f'class="report-page section-bundle-body section-{major.number}-page">'
            ),
        ]
        if include_major_heading:
            first_page_parts.append(
                f'<h2 class="major-heading" style="color:{REPORT_STYLE.green_major};">'
                f"{_escape(heading)}</h2>"
            )
        continuation_2a: list[str] = []
        segments_2b: list[str] = []
        other_subs: list[str] = []
        has_2a = False
        has_2b = False
        for sub in major.subsections:
            if sub.key == "a" and any(
                str(b.presentation_role or "").startswith("section_2a_")
                for b in (sub.presentation_blocks or [])
            ):
                has_2a = True
                segments = render_section_2a_subsection_segments(sub)
                first_page_parts.append(segments[0])
                continuation_2a.extend(segments[1:])
            elif sub.key == "b" and any(
                str(b.presentation_role or "").startswith("section_2b_")
                for b in (sub.presentation_blocks or [])
            ):
                has_2b = True
                segments_2b = render_section_2b_subsection_segments(sub)
            else:
                other_subs.append(_render_subsection(sub))

        # Legacy non-2a/2b Major 2 content stays on page 1 (should be empty today).
        first_page_parts.extend(other_subs)
        first_page_parts.append("</section>")

        if has_2a:
            body_parts.extend(first_page_parts)
            for index, segment in enumerate(continuation_2a):
                body_parts.append(page_break)
                body_parts.append(
                    f'<section id="section-{major.number}-cont-{index + 2}" '
                    f'class="report-page section-bundle-body section-2a-continuation">'
                    f"{segment}</section>"
                )
            if has_2b and segments_2b:
                # Page break after all 2a pages, then 2b segment 0 (h3 once).
                body_parts.append(page_break)
                body_parts.append(
                    f'<section id="section-{major.number}-2b" '
                    f'class="report-page section-bundle-body section-2b-page">'
                    f"{segments_2b[0]}</section>"
                )
                for index, segment in enumerate(segments_2b[1:]):
                    body_parts.append(page_break)
                    body_parts.append(
                        f'<section id="section-{major.number}-2b-cont-{index + 2}" '
                        f'class="report-page section-bundle-body section-2b-continuation">'
                        f"{segment}</section>"
                    )
        elif has_2b and segments_2b:
            # Focused 2b (or Major 2 without 2a): heading already on first page.
            first_page_parts.insert(-1, segments_2b[0])
            body_parts.extend(first_page_parts)
            for index, segment in enumerate(segments_2b[1:]):
                body_parts.append(page_break)
                body_parts.append(
                    f'<section id="section-{major.number}-2b-cont-{index + 2}" '
                    f'class="report-page section-bundle-body section-2b-continuation">'
                    f"{segment}</section>"
                )
        else:
            body_parts.extend(first_page_parts)

    body = "\n".join(body_parts)

    from gene_dossier.rancho_report import _asset_data_uri, _img_tag

    header = ""
    footer = ""
    if include_page_chrome:
        rancho = _asset_data_uri("rancho_wordmark.png")
        chdi = _asset_data_uri("chdi_wordmark.png")
        rancho_header = _asset_data_uri("rancho_header_bar.png") or rancho
        rancho_footer = _asset_data_uri("rancho_footer.png")
        if show_header_logos:
            header_inner = (
                '<div class="page-header">'
                f'{_img_tag(rancho_header, cls="rancho", alt="Rancho BioSciences")}'
                f'{_img_tag(chdi, cls="chdi", alt="CHDI Foundation")}'
                "</div>"
            )
            footer_inner = (
                '<div class="page-footer">'
                f'{_img_tag(rancho_footer, cls="rancho", alt="Rancho BioSciences")}'
                f"<span>{_escape(REPORT_STYLE.footer_url)}</span>"
                "</div>"
            )
        else:
            header_inner = '<div class="page-header"></div>'
            footer_inner = (
                f'<div class="page-footer">'
                f"<span>{_escape(REPORT_STYLE.footer_url)}</span></div>"
            )
        header = f'<div class="report-page report-chrome">{header_inner}</div>'
        footer = f'<div class="report-page report-chrome">{footer_inner}</div>'

    major_nums = [sec.number for sec in document.sections]
    if major_nums == [1]:
        title_suffix = "Section 1"
    elif major_nums == [2]:
        title_suffix = "Section 2"
    else:
        title_suffix = "Section Bundle"
    title = _escape(f"{document.gene_symbol} — {title_suffix}")
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title}</title>
  <style>
{_rancho_css()}
  </style>
</head>
<body>
{header}
{body}
{footer}
</body>
</html>
"""


def _render_major_section_1_pages(
    major: ReportMajorSection,
    *,
    include_major_heading: bool,
) -> list[str]:
    """Preserve existing Section 1 segmentation (1c/1d/1e) exactly."""
    heading = f"{major.number}. {major.title}"
    body_parts = [
        (
            f'<section id="section-{major.number}" '
            f'class="report-page section-bundle-body">'
        ),
    ]
    if include_major_heading:
        body_parts.append(
            f'<h2 class="major-heading" style="color:{REPORT_STYLE.green_major};">'
            f"{_escape(heading)}</h2>"
        )
    continuation_segments: list[str] = []
    subsection_d = next((s for s in major.subsections if s.key == "d"), None)
    subsection_c = next((s for s in major.subsections if s.key == "c"), None)
    merge_1d_onto_pdb = bool(
        subsection_c
        and subsection_c.presentation_blocks
        and subsection_d
        and subsection_d.presentation_blocks
    )
    has_prior_to_1e = any(s.key != "e" for s in major.subsections)

    for sub in major.subsections:
        if sub.key == "d" and merge_1d_onto_pdb:
            continue
        if sub.key == "e":
            e_html = _render_subsection(sub)
            if has_prior_to_1e:
                continuation_segments.append(e_html)
            else:
                body_parts.append(e_html)
            continue
        if sub.key == "c" and sub.presentation_blocks:
            segments = render_section_1c_subsection_segments(sub)
            body_parts.extend(segments[:1])
            rest = list(segments[1:])
            if merge_1d_onto_pdb and subsection_d is not None:
                d_html = _render_subsection(subsection_d)
                if rest:
                    last_1c_segment_parts = [rest[-1], d_html]
                    rest[-1] = "\n".join(last_1c_segment_parts)
                else:
                    if body_parts:
                        body_parts[-1] = "\n".join([body_parts[-1], d_html])
            continuation_segments.extend(rest)
        else:
            body_parts.append(_render_subsection(sub))
    body_parts.append("</section>")
    out = list(body_parts)
    for index, segment in enumerate(continuation_segments):
        out.append(SECTION_1C_PDF_PAGE_BREAK)
        is_1e_page = 'class="report-subsection subsection-e"' in segment
        page_class = (
            "report-page section-bundle-body section-1e-page"
            if is_1e_page
            else "report-page section-bundle-body section-1c-continuation"
        )
        out.append(
            f'<section id="section-{major.number}-cont-{index + 2}" '
            f'class="{page_class}">'
            f"{segment}</section>"
        )
    return out


_BUNDLE_OUTPUT_NAMES = (
    "section_1.json",
    "section_1_audit.json",
    "section_1.html",
    "section_1.pdf",
    "section_1.png",
    "section_1_contact_sheet.png",
)


def _bundle_dir_is_populated(out: Path) -> bool:
    if not out.is_dir():
        return False
    for name in _BUNDLE_OUTPUT_NAMES:
        if (out / name).exists():
            return True
    if any(out.glob("section_1_page_*.png")):
        return True
    return False


def _write_json(path: Path, payload: Any, *, credentials_only: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cleaned = (
        sanitize_credentials(payload)
        if credentials_only
        else sanitize_credentials(payload)
    )
    path.write_text(
        json.dumps(cleaned, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _cleanup_attempt_outputs(paths: Sequence[Path]) -> None:
    for path in paths:
        try:
            if path.is_file():
                path.unlink(missing_ok=True)
        except OSError:
            logger.warning("Failed to remove incomplete bundle output %s", path)


def write_section_bundle_outputs(
    *,
    document: ReportDocument,
    presentation: dict[str, Any],
    audit: dict[str, Any],
    output_dir: str | Path,
    write_pdf: bool = True,
    dpi: int = 150,
    allow_rerender: bool = False,
    include_major_heading: bool = True,
) -> dict[str, Path]:
    """Write presentation/audit JSON, HTML, PDF, and all PDF page PNGs.

    Rejects an already-populated run output directory unless ``allow_rerender``.
    On failure, deletes only files that did not exist before this attempt.
    """
    out = Path(output_dir)
    if out.exists() and _bundle_dir_is_populated(out) and not allow_rerender:
        raise SectionBundleError(
            f"Section-bundle output directory already populated: {out}. "
            "Pass allow_rerender=True to overwrite."
        )
    out.mkdir(parents=True, exist_ok=True)

    presentation_path = out / "section_1.json"
    audit_path = out / "section_1_audit.json"
    html_path = out / "section_1.html"
    pdf_path = out / "section_1.pdf"

    existed_before = {
        presentation_path: presentation_path.exists(),
        audit_path: audit_path.exists(),
        html_path: html_path.exists(),
        pdf_path: pdf_path.exists(),
    }
    # Track prior PNG stems too.
    for prior in list(out.glob("section_1.png")) + list(out.glob("section_1_page_*.png")):
        existed_before[prior] = True

    newly_created: list[Path] = []
    paths: dict[str, Path] = {}

    def _mark_written(path: Path) -> None:
        if not existed_before.get(path, False):
            newly_created.append(path)

    try:
        _write_json(presentation_path, presentation, credentials_only=True)
        _mark_written(presentation_path)
        paths["section_1_json"] = presentation_path

        # Audit: credentials only — preserve [source_id=...] diagnostics.
        audit_path.write_text(
            json.dumps(
                sanitize_credentials(audit),
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )
        _mark_written(audit_path)
        paths["section_1_audit_json"] = audit_path

        html = render_section_bundle_html(
            document,
            include_page_chrome=True,
            include_major_heading=include_major_heading,
        )
        html = redact_api_key(html)
        html_path.write_text(html, encoding="utf-8")
        _mark_written(html_path)
        paths["section_1_html"] = html_path

        if write_pdf:
            section_one_keys = [
                sub.key
                for sec in document.sections
                if sec.number == 1
                for sub in sec.subsections
            ]
            stamp_first_page = section_one_keys == ["c"]
            # Focused 1c previews have no cover/earlier page; stamp page 1 so the
            # visual preview matches the Rancho body-page chrome without changing
            # assembled 1a/1b output.
            # Omit HTML chrome to avoid one-shot header/footer and stacked print padding.
            pdf_html = render_section_bundle_html(
                document,
                include_page_chrome=False,
                include_major_heading=include_major_heading,
            )
            pdf_html = redact_api_key(pdf_html)
            rendered = render_rancho_pdf(
                pdf_html,
                pdf_path,
                page_size="letter",
                stamp_page_chrome=True,
                stamp_cover=stamp_first_page,
            )
            if rendered is not None:
                _mark_written(Path(rendered))
                paths["section_1_pdf"] = Path(rendered)
                pngs = rasterize_pdf_pages_to_pngs(
                    rendered, out, stem="section_1", dpi=dpi
                )
                for index, png in enumerate(pngs):
                    if not existed_before.get(png, False):
                        newly_created.append(png)
                    if len(pngs) == 1:
                        paths["section_1_png"] = png
                    else:
                        paths[f"section_1_page_{index + 1}_png"] = png
        return paths
    except Exception:
        _cleanup_attempt_outputs(newly_created)
        raise


def run_section_bundle(
    gene_symbol: str,
    *,
    section_keys: Sequence[str] | None = None,
    output_dir: str | Path | None = None,
    settings: Settings | None = None,
    call_network: bool = True,
    persist_db: bool = True,
    write_pdf: bool = True,
    dpi: int = 150,
    dossier_run_id: str | None = None,
    allow_rerender: bool = False,
    preloaded_state: DossierState | None = None,
    acceptance_profile: str | None = None,
    section_1e_config: Section1eConfig | None = None,
    section_2a_config: Section2aConfig | None = None,
    section_2b_config: Section2bConfig | None = None,
) -> SectionBundleResult:
    """Execute identity (+ section-owned sources) and write a section bundle."""
    cfg = settings or get_settings()
    keys = validate_section_keys(section_keys or DEFAULT_SECTION_BUNDLE_KEYS)
    gene = gene_symbol.strip()
    run, state = create_section_bundle_run(
        gene_symbol=gene,
        selected_section_keys=keys,
        settings=cfg,
        persist_db=persist_db,
        dossier_run_id=dossier_run_id,
    )
    if preloaded_state is not None:
        # Offline / test path: merge preloaded tool/evidence into the new state.
        for key in (
            "gene_ids",
            "tool_results",
            "api_runs",
            "raw_artifacts",
            "evidence_records",
            "errors",
            "section_1d_status",
            "section_1e_status",
            "section_2a_status",
            "section_2b_status",
            "coverage",
        ):
            if key in preloaded_state and preloaded_state[key] is not None:
                state[key] = preloaded_state[key]  # type: ignore[literal-required]

    run_id = run.id
    gene_dir = Path(output_dir) if output_dir else (cfg.output_path / "section_validation" / gene)
    out_dir = gene_dir / run_id
    transient = WorkflowTransientContext()
    created_outputs: dict[str, Path] = {}
    errors: list[str] = []
    focused_1c = keys == ["1c"]
    focused_1d = keys == ["1d"]
    focused_1e = keys == ["1e"]
    focused_2a = keys == ["2a"]
    include_major_heading = not (focused_1c or focused_1d or focused_1e)

    try:
        if call_network:
            state = node_resolve_gene_identity(state, settings=cfg)
            needed = sources_for_sections(keys)
            if needed:
                state = node_call_source_clients(
                    state,
                    settings=cfg,
                    sources=needed,
                    call_network=True,
                    transient=transient,
                )
            state = node_save_raw_artifacts(
                state,
                settings=cfg,
                persist_db=persist_db,
                transient=transient,
            )
            state = node_normalize_evidence(state, persist_db=persist_db)
            if "1c" in keys:
                state = {
                    **state,
                    "run_type": "section_bundle",
                    "selected_section_keys": list(keys),
                    "acceptance_profile": acceptance_profile,
                }
                state = node_generate_section_1c_derived_artifacts(
                    state,
                    settings=cfg,
                    persist_db=persist_db,
                    transient=transient,
                )
            if "1d" in keys:
                state = {
                    **state,
                    "run_type": "section_bundle",
                    "selected_section_keys": list(keys),
                    "acceptance_profile": acceptance_profile,
                }
                state = node_generate_section_1d_derived_artifacts(
                    state,
                    settings=cfg,
                    persist_db=persist_db,
                    transient=transient,
                )
            if "1e" in keys:
                state = {
                    **state,
                    "run_type": "section_bundle",
                    "selected_section_keys": list(keys),
                    "acceptance_profile": acceptance_profile,
                }
                state = node_generate_section_1e_derived_artifacts(
                    state,
                    settings=cfg,
                    persist_db=persist_db,
                    transient=transient,
                    config=section_1e_config or Section1eConfig(),
                )
            if "2a" in keys:
                state = {
                    **state,
                    "run_type": "section_bundle",
                    "selected_section_keys": list(keys),
                    "acceptance_profile": acceptance_profile,
                }
                state = node_generate_section_2a_derived_artifacts(
                    state,
                    settings=cfg,
                    persist_db=persist_db,
                    transient=transient,
                    config=section_2a_config or Section2aConfig(),
                )
            if "2b" in keys:
                state = {
                    **state,
                    "run_type": "section_bundle",
                    "selected_section_keys": list(keys),
                    "acceptance_profile": acceptance_profile,
                }
                state = node_generate_section_2b_derived_artifacts(
                    state,
                    settings=cfg,
                    persist_db=persist_db,
                    transient=transient,
                    config=section_2b_config or Section2bConfig(),
                )

        evidence = list(state.get("evidence_records") or [])
        section_status_by_key: dict[str, Any] = {}
        if state.get("section_1d_status"):
            section_status_by_key["1d"] = state["section_1d_status"]
        if state.get("section_1e_status"):
            section_status_by_key["1e"] = state["section_1e_status"]
        if state.get("section_2a_status"):
            section_status_by_key["2a"] = state["section_2a_status"]
        if state.get("section_2b_status"):
            section_status_by_key["2b"] = state["section_2b_status"]
        document, presentation, audit = build_section_bundle_document(
            dossier_run_id=run_id,
            gene_symbol=gene,
            section_keys=keys,
            evidence_records=evidence,
            api_runs=list(state.get("api_runs") or []),
            raw_artifacts=list(state.get("raw_artifacts") or []),
            section_status_by_key=section_status_by_key,
        )
        coverage = coverage_updates_from_state(state)
        # Include Section 1d/1e coverage diagnostics emitted outside the generic client.
        for row in list(state.get("coverage") or []):
            if getattr(row, "source_name", None) in {
                "AlphaFold",
                "NCBI Datasets",
                "OrthoDB",
                "NCBI Gene",
                "GTEx",
                "Human Brain Transcriptome",
                "Allen Brain Atlas",
                "BrainRNASeq",
            }:
                coverage.append(row)
        audit["coverage"] = [
            {
                "source_name": row.source_name,
                "status": row.status.value if hasattr(row.status, "value") else str(row.status),
                "evidence_record_count": row.evidence_record_count,
                "error_message": row.error_message,
                "notes": getattr(row, "notes", None),
            }
            for row in coverage
        ]
        transcript_prov: dict[str, Any] = {}
        figure_prov: dict[str, Any] = {}
        for rec in evidence:
            if rec.fact_type == "ucsc_canonical_transcript" and isinstance(rec.value, dict):
                transcript_prov = {
                    "evidence_record_id": rec.id,
                    "source_id": rec.source_id,
                    "raw_artifact_id": rec.raw_artifact_id,
                    "api_run_id": rec.api_run_id,
                    "transcript_id": rec.value.get("transcript_id"),
                    "is_mane_select": rec.value.get("is_mane_select"),
                    "is_ensembl_canonical": rec.value.get("is_ensembl_canonical"),
                    "is_gencode_primary": rec.value.get("is_gencode_primary"),
                    "is_canonical_tier": rec.value.get("is_canonical_tier"),
                    "selection_reasons": rec.value.get("selection_reasons"),
                    "display_position": rec.value.get("display_position"),
                }
            if rec.fact_type == "ucsc_conservation_figure" and isinstance(rec.value, dict):
                figure_prov = {
                    "evidence_record_id": rec.id,
                    "source_id": rec.source_id,
                    "raw_artifact_id": rec.raw_artifact_id
                    or rec.value.get("figure_raw_artifact_id"),
                    "api_run_id": rec.api_run_id or rec.value.get("figure_api_run_id"),
                    "relative_path": rec.value.get("relative_path"),
                    "sha256": rec.value.get("sha256") or rec.value.get("content_hash"),
                    "width": rec.value.get("width"),
                    "height": rec.value.get("height"),
                    "media_type": rec.value.get("media_type"),
                    "track_preset_id": rec.value.get("track_preset_id"),
                    "retrieval_method": rec.value.get("retrieval_method"),
                    "source_note": rec.value.get("source_note"),
                    "caption": rec.value.get("caption"),
                }
        audit["transcript_selection_provenance"] = transcript_prov
        audit["figure_provenance"] = figure_prov
        if "1c" in keys:
            audit["section_1c"] = sanitize_credentials(state.get("section_1c") or {})
            if acceptance_profile:
                section_1c_audit = dict(audit.get("section_1c") or {})
                rendering = section_1c_audit.get("rendering_status") or {}
                reasons: list[str] = []
                if (
                    acceptance_profile == "section_1c_reference_genes"
                    and rendering.get("cdd_architecture") != "success"
                ):
                    reasons.append("official CDD architecture missing")
                section_1c_audit["acceptance_validation"] = {
                    "profile": acceptance_profile,
                    "status": "failed" if reasons else "success",
                    "reasons": reasons,
                }
                audit["section_1c"] = section_1c_audit
                if reasons:
                    errors.extend(reasons)
        if "1d" in keys and state.get("section_1d_status"):
            audit["section_1d"] = sanitize_credentials(
                (state.get("section_1d_status") or {}).get("audit")
                or state.get("section_1d_status")
            )
        if "1e" in keys and state.get("section_1e_status"):
            audit["section_1e"] = sanitize_credentials(
                (state.get("section_1e_status") or {}).get("audit")
                or state.get("section_1e_status")
            )
        if "2a" in keys and state.get("section_2a_status"):
            audit["section_2a"] = sanitize_credentials(
                (state.get("section_2a_status") or {}).get("audit")
                or state.get("section_2a_status")
            )
        if "2b" in keys and state.get("section_2b_status"):
            audit["section_2b"] = sanitize_credentials(
                (state.get("section_2b_status") or {}).get("audit")
                or state.get("section_2b_status")
            )
        audit["errors"] = list(state.get("errors") or [])
        if errors:
            audit["errors"] = list(dict.fromkeys([*audit["errors"], *errors]))
        audit = sanitize_credentials(audit)

        created_outputs = write_section_bundle_outputs(
            document=document,
            presentation=presentation,
            audit=audit,
            output_dir=out_dir,
            write_pdf=write_pdf,
            dpi=dpi,
            allow_rerender=allow_rerender,
            include_major_heading=include_major_heading,
        )

        if "1d" in keys and acceptance_profile == "section_1d_reference_genes":
            section_status = dict(state.get("section_1d_status") or {})
            d_blocks: list[ReportContentBlock] = []
            for sec in document.sections:
                for sub in sec.subsections:
                    if sub.key == "d":
                        d_blocks = list(sub.presentation_blocks or [])
            html_text = None
            html_path = created_outputs.get("section_1_html")
            if html_path and Path(html_path).is_file():
                html_text = Path(html_path).read_text(encoding="utf-8")
            pdf_path = created_outputs.get("section_1_pdf")
            reasons = evaluate_section_1d_reference_genes_acceptance(
                gene_symbol=gene,
                section_status=section_status,
                presentation_blocks=d_blocks,
                html=html_text,
                pdf_path=Path(pdf_path) if pdf_path else None,
                selected_section_keys=keys,
            )
            section_1d_audit = dict(audit.get("section_1d") or {})
            section_1d_audit["acceptance_validation"] = {
                "profile": acceptance_profile,
                "status": "failed" if reasons else "success",
                "reasons": reasons,
            }
            audit["section_1d"] = sanitize_credentials(section_1d_audit)
            if reasons:
                errors.extend(reasons)
                audit["errors"] = list(
                    dict.fromkeys([*(audit.get("errors") or []), *errors])
                )
            audit_path = created_outputs.get("section_1_audit_json")
            if audit_path:
                Path(audit_path).write_text(
                    json.dumps(
                        sanitize_credentials(audit),
                        indent=2,
                        sort_keys=True,
                        ensure_ascii=False,
                    )
                    + "\n",
                    encoding="utf-8",
                )

        finalize_section_bundle_run(
            dossier_run_id=run_id,
            status="completed",
            selected_section_keys=keys,
            persist_db=persist_db,
        )
        return SectionBundleResult(
            gene_symbol=gene,
            dossier_run_id=run_id,
            selected_section_keys=list(keys),
            output_dir=out_dir,
            output_paths=created_outputs,
            status="completed",
            errors=list(sanitize_credentials(errors)),
        )
    except Exception as exc:  # noqa: BLE001
        logger.exception("Section bundle failed for %s", gene)
        errors.append(str(exc))
        _cleanup_attempt_outputs(list(created_outputs.values()))
        try:
            finalize_section_bundle_run(
                dossier_run_id=run_id,
                status="failed",
                selected_section_keys=keys,
                errors=errors,
                persist_db=persist_db,
            )
        except SectionBundleError as finalize_exc:
            errors.append(str(finalize_exc))
        return SectionBundleResult(
            gene_symbol=gene,
            dossier_run_id=run_id,
            selected_section_keys=list(keys),
            output_dir=out_dir,
            output_paths={},
            status="failed",
            errors=list(sanitize_credentials(errors)),
        )
    finally:
        transient.clear_run(run_id)


__all__ = [
    "SUPPORTED_SECTION_BUNDLE_KEYS",
    "DEFAULT_SECTION_BUNDLE_KEYS",
    "SECTION_SOURCE_DEPENDENCIES",
    "SectionBundleError",
    "SectionBundleResult",
    "validate_section_keys",
    "sources_for_sections",
    "opaque_evidence_ref",
    "section_1c_evidence_ref",
    "validate_safe_item_key",
    "sanitize_credentials",
    "sanitize_polished_text",
    "sanitize_secrets",
    "build_provenance_index",
    "create_section_bundle_run",
    "finalize_section_bundle_run",
    "build_section_bundle_document",
    "render_section_bundle_html",
    "write_section_bundle_outputs",
    "run_section_bundle",
    "assign_opaque_refs",
    "serialize_presentation_block",
]
