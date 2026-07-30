"""Bundle-only Section 1d AlphaFold helpers.

Owns all Section 1d AlphaFold network requests (human / mouse / rat). Not wired
into the full dossier workflow's generic human-only AlphaFold client path.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence
from urllib.parse import urlparse

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
    _validate_nonblank_image,
    protein_seeds_by_species,
)
from gene_dossier.source_ids import make_source_id, slugify
from gene_dossier.tools import alphafold
from gene_dossier.workflow import DossierState, WorkflowTransientContext

logger = logging.getLogger(__name__)

SECTION_STRUCTURE = "AlphaFold / PDBe / CDD"
SUBSECTION_1D = "AlphaFold prediction"
ALPHAFOLD_HOST = "alphafold.ebi.ac.uk"

SPECIES_ORDER: tuple[tuple[str, int, str], ...] = (
    ("human", 9606, "Human"),
    ("rat", 10116, "Rat"),
    ("mouse", 10090, "Mouse"),
)

STATUS_ACCESSION_UNAVAILABLE = "accession_unavailable"
STATUS_MODEL_ABSENT = "model_absent"
STATUS_REQUEST_FAILED = "request_failed"
STATUS_SELECTED = "selected"
STATUS_CAPTURE_FAILED = "capture_failed"
STATUS_VISUALIZATION_UNAVAILABLE = "visualization_unavailable"

CAPTURE_FRESH = "fresh_official_capture"
CAPTURE_REUSED = "reused_official_capture"
CAPTURE_UNAVAILABLE = "capture_unavailable"

VISUALIZATION_UNAVAILABLE_TEXT = (
    "AlphaFold structure visualization temporarily unavailable"
)

_MIN_CAPTURE_WIDTH = 500
_MIN_CAPTURE_HEIGHT = 350
_MIN_CAPTURE_BYTES = 20_000
_ALPHAFOLD_ALLOWED_HOSTS = frozenset({ALPHAFOLD_HOST, f"www.{ALPHAFOLD_HOST}"})


@dataclass(frozen=True)
class SpeciesSlot:
    """Resolved identity + selection/status for one Section 1d species row."""

    species_key: str
    species_label: str
    taxon_id: int
    display_symbol: str
    accession: str | None
    status: str
    model_entity_id: str | None = None
    entry_url: str | None = None
    model_version: int | None = None
    message: str | None = None
    selection_diagnostics: tuple[dict[str, Any], ...] = ()


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
    taxon_id: int | None = None,
    organism: str | None = None,
    species: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=make_source_id(
            "AlphaFold", gene_symbol, AssertionType.protein_structure, key
        ),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_STRUCTURE,
        subsection=SUBSECTION_1D,
        source_name="AlphaFold",
        source_type=SourceType.structure_database,
        assertion_type=AssertionType.protein_structure,
        fact_type=fact_type,
        organism=organism,
        species=species,
        taxon_id=taxon_id,
        evidence_grade=EvidenceGrade.E,
        confidence_notes=(
            "Predicted structure from AlphaFold DB; not an experimental model."
        ),
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def _symbol_from_evidence(
    evidence_records: Sequence[EvidenceRecord],
    *,
    taxon_id: int,
    gene_ids: dict[str, Any],
    gene_symbol: str,
) -> str:
    """Prefer Section 1a identity symbols; fall back to gene_ids / request symbol."""
    if taxon_id == 9606:
        for key in ("official_symbol", "gene_symbol"):
            value = gene_ids.get(key)
            if value:
                return str(value).strip()
        return gene_symbol.strip()
    if taxon_id == 10090:
        value = gene_ids.get("mouse_symbol")
        if value:
            return str(value).strip()
    if taxon_id == 10116:
        value = gene_ids.get("rat_symbol")
        if value:
            return str(value).strip()

    for rec in evidence_records:
        if rec.taxon_id != taxon_id:
            continue
        value = rec.value if isinstance(rec.value, dict) else {}
        for key in (
            "nomenclaturesymbol",
            "gene_symbol",
            "symbol",
            "display_name",
            "gene_names",
        ):
            raw = value.get(key)
            if isinstance(raw, list) and raw:
                return str(raw[0]).strip()
            if raw:
                return str(raw).strip()
    # Ortholog casing heuristic when only the human symbol is known.
    base = gene_symbol.strip()
    if taxon_id in {10090, 10116} and base:
        return base[0] + base[1:].lower() if len(base) > 1 else base.upper()
    return base


def _taxon_of_record(rec: EvidenceRecord) -> int | None:
    if rec.taxon_id is not None:
        try:
            return int(rec.taxon_id)
        except (TypeError, ValueError):
            pass
    value = rec.value if isinstance(rec.value, dict) else {}
    for key in ("taxon_id", "tax_id", "organism_id"):
        raw = value.get(key)
        if raw is None:
            continue
        try:
            return int(raw)
        except (TypeError, ValueError):
            continue
    return None


def _accession_from_record_value(value: dict[str, Any]) -> str | None:
    for key in ("uniprot_accession", "uniprot_id", "primaryAccession", "accession"):
        raw = value.get(key)
        if isinstance(raw, list):
            for item in raw:
                text = str(item or "").strip()
                if text:
                    return text
        elif raw:
            text = str(raw).strip()
            if text:
                return text
    return None


def _accession_for_taxon(
    evidence_records: Sequence[EvidenceRecord],
    *,
    taxon_id: int,
    gene_ids: dict[str, Any],
) -> str | None:
    """Resolve UniProt accession using the same identity breadth as Section 1a."""
    # 1) Taxon-specific UniProt accession evidence.
    for rec in evidence_records:
        if rec.source_name != "UniProt" or rec.fact_type != "uniprot_accession":
            continue
        if _taxon_of_record(rec) != taxon_id:
            continue
        value = rec.value if isinstance(rec.value, dict) else {}
        accession = _accession_from_record_value(value)
        if accession:
            return accession

    # 2) Any taxon-specific ortholog/identity evidence carrying UniProt fields.
    for rec in evidence_records:
        if _taxon_of_record(rec) != taxon_id:
            continue
        value = rec.value if isinstance(rec.value, dict) else {}
        accession = _accession_from_record_value(value)
        if accession:
            return accession

    # 3) gene_ids species accession fallback.
    fallback = {
        9606: gene_ids.get("uniprot_accession"),
        10090: gene_ids.get("mouse_uniprot_accession"),
        10116: gene_ids.get("rat_uniprot_accession"),
    }.get(taxon_id)
    if fallback:
        text = str(fallback).strip()
        return text or None
    return None


def resolve_species_identity(
    *,
    gene_symbol: str,
    gene_ids: dict[str, Any],
    evidence_records: Sequence[EvidenceRecord],
) -> list[dict[str, Any]]:
    """Resolve Human/Rat/Mouse symbols and UniProt accessions from identity evidence."""
    seeds = {seed.taxon_id: seed for seed in protein_seeds_by_species(evidence_records)}
    out: list[dict[str, Any]] = []
    for species_key, taxon_id, label in SPECIES_ORDER:
        seed = seeds.get(taxon_id)
        accession = _accession_for_taxon(
            evidence_records, taxon_id=taxon_id, gene_ids=gene_ids
        )
        symbol = _symbol_from_evidence(
            evidence_records,
            taxon_id=taxon_id,
            gene_ids=gene_ids,
            gene_symbol=gene_symbol,
        )
        out.append(
            {
                "species_key": species_key,
                "species_label": label,
                "taxon_id": taxon_id,
                "display_symbol": symbol,
                "accession": accession,
                "seed": seed,
            }
        )
    return out


def unavailable_visible_text(
    *,
    species_label: str,
    display_symbol: str,
    status: str,
) -> str:
    """Plain-text fallback line for polished Section 1d (non-evidence)."""
    prefix = f"{species_label} {display_symbol}"
    if status == STATUS_ACCESSION_UNAVAILABLE:
        return f"{prefix}: UniProt accession not available"
    if status == STATUS_REQUEST_FAILED:
        return f"{prefix}: AlphaFold prediction temporarily unavailable"
    return f"{prefix}: AlphaFold prediction not available"


def _capture_passes_quality(meta: dict[str, Any], content: bytes) -> tuple[bool, str | None]:
    width = int(meta.get("width") or 0)
    height = int(meta.get("height") or 0)
    byte_size = int(meta.get("byte_size") or len(content))
    if width < _MIN_CAPTURE_WIDTH or height < _MIN_CAPTURE_HEIGHT:
        return False, f"capture dimensions too small ({width}x{height})"
    if byte_size < _MIN_CAPTURE_BYTES:
        return False, f"capture byte size too small ({byte_size})"
    return True, None


def _is_allowlisted_alphafold_entry_url(url: str | None, *, model_id: str) -> bool:
    if not url:
        return False
    parsed = urlparse(str(url))
    if parsed.hostname not in _ALPHAFOLD_ALLOWED_HOSTS:
        return False
    path = parsed.path or ""
    return f"/entry/{model_id}" in path


def _capture_cache_dir(settings: Settings) -> Path:
    path = Path(settings.raw_data_path) / "_section_1d_capture_cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_capture_cache_index(
    *,
    settings: Settings,
    model_id: str,
    payload: dict[str, Any],
) -> None:
    path = _capture_cache_dir(settings) / f"{slugify(model_id)}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _load_capture_cache_index(
    *,
    settings: Settings,
    model_id: str,
) -> dict[str, Any] | None:
    path = _capture_cache_dir(settings) / f"{slugify(model_id)}.json"
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return None
    return data if isinstance(data, dict) else None


def _validate_capture_candidate(
    *,
    value: dict[str, Any],
    model_id: str,
    accession: str,
    model_version: int | None,
    settings: Settings,
) -> tuple[dict[str, Any] | None, str | None]:
    """Return validated capture value metadata or a rejection reason."""
    if str(value.get("model_entity_id") or "") != model_id:
        return None, "model_entity_id_mismatch"
    if str(value.get("uniprot_accession") or "") != accession:
        return None, "accession_mismatch"
    if value.get("artifact_class") != "derived_capture":
        return None, "artifact_class_mismatch"
    if not _is_allowlisted_alphafold_entry_url(
        str(value.get("source_page_url") or ""), model_id=model_id
    ):
        return None, "source_page_url_not_entry"
    known_version = value.get("model_version")
    if (
        model_version is not None
        and known_version is not None
        and int(known_version) != int(model_version)
    ):
        return None, "model_version_mismatch"

    from gene_dossier.ucsc_figure import resolve_artifact_path, sha256_hex

    rel = value.get("relative_path") or value.get("local_artifact_path")
    if not rel:
        return None, "missing_relative_path"
    try:
        path = resolve_artifact_path(str(rel), root=settings.raw_data_path)
    except Exception:  # noqa: BLE001
        return None, "unresolvable_path"
    if not path.is_file():
        return None, "missing_managed_file"
    content = path.read_bytes()
    expected = str(value.get("sha256") or value.get("content_hash") or "")
    digest = sha256_hex(content)
    if expected and digest != expected:
        return None, "checksum_mismatch"
    try:
        validation = _validate_nonblank_image(content)
    except Exception as exc:  # noqa: BLE001
        return None, f"image_invalid:{exc}"
    ok, reason = _capture_passes_quality(validation, content)
    if not ok:
        return None, reason or "quality_gate_failed"
    enriched = {
        **value,
        "relative_path": str(rel),
        "sha256": digest,
        "width": validation.get("width") or value.get("width"),
        "height": validation.get("height") or value.get("height"),
        "byte_size": validation.get("byte_size") or len(content),
        "media_type": validation.get("media_type") or value.get("media_type") or "image/png",
        "artifact_class": "derived_capture",
    }
    return enriched, None


def find_reusable_official_capture(
    *,
    evidence_records: Sequence[EvidenceRecord],
    raw_artifacts: Sequence[Any] | None,
    model_id: str,
    accession: str,
    model_version: int | None,
    settings: Settings,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Find a previously validated official AlphaFold viewer capture."""
    diagnostics: list[dict[str, Any]] = []
    candidates: list[dict[str, Any]] = []

    for rec in evidence_records:
        if rec.source_name != "AlphaFold" or rec.fact_type != "alphafold_official_viewer_capture":
            continue
        value = rec.value if isinstance(rec.value, dict) else {}
        candidates.append(
            {
                **value,
                "figure_raw_artifact_id": value.get("figure_raw_artifact_id") or rec.raw_artifact_id,
                "reuse_source": "evidence_record",
                "reuse_evidence_record_id": rec.id,
            }
        )

    for meta in raw_artifacts or []:
        if not isinstance(meta, dict):
            continue
        if meta.get("source_name") != "AlphaFold":
            continue
        if meta.get("artifact_class") != "derived_capture":
            continue
        notes = meta.get("notes")
        if isinstance(notes, str):
            try:
                notes = json.loads(notes)
            except Exception:  # noqa: BLE001
                notes = {}
        if not isinstance(notes, dict):
            notes = {}
        payload = {
            "model_entity_id": notes.get("model_entity_id") or meta.get("model_entity_id"),
            "uniprot_accession": notes.get("uniprot_accession") or meta.get("uniprot_accession"),
            "model_version": notes.get("model_version") or meta.get("model_version"),
            "source_page_url": notes.get("source_page_url") or meta.get("source_page_url"),
            "relative_path": meta.get("relative_path") or meta.get("file_path"),
            "sha256": meta.get("expected_sha256") or meta.get("content_hash") or notes.get("sha256"),
            "width": meta.get("width") or notes.get("width"),
            "height": meta.get("height") or notes.get("height"),
            "byte_size": meta.get("byte_size") or notes.get("byte_size"),
            "artifact_class": "derived_capture",
            "figure_raw_artifact_id": meta.get("id"),
            "reuse_source": "raw_artifact_meta",
        }
        candidates.append(payload)

    cached = _load_capture_cache_index(settings=settings, model_id=model_id)
    if cached:
        candidates.append({**cached, "reuse_source": "capture_cache_index"})

    for candidate in candidates:
        validated, reason = _validate_capture_candidate(
            value=candidate,
            model_id=model_id,
            accession=accession,
            model_version=model_version,
            settings=settings,
        )
        if validated is None:
            diagnostics.append(
                {
                    "code": "reuse_rejected",
                    "reason": reason,
                    "reuse_source": candidate.get("reuse_source"),
                    "model_entity_id": candidate.get("model_entity_id"),
                }
            )
            continue
        diagnostics.append(
            {
                "code": "reuse_accepted",
                "reuse_source": candidate.get("reuse_source"),
                "model_entity_id": model_id,
                "relative_path": validated.get("relative_path"),
            }
        )
        return validated, diagnostics

    diagnostics.append(
        {
            "code": "reuse_unavailable",
            "model_entity_id": model_id,
            "message": "No validated reusable official viewer capture found",
        }
    )
    return None, diagnostics


def _capture_human_viewer(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    model_id: str,
    accession: str,
    taxon_id: int,
    model_version: int | None,
    parent_raw_artifact_ids: Sequence[str],
    settings: Settings,
    persist_db: bool,
    max_attempts: int = 3,
    headed: bool = False,
    channel: str | None = None,
    user_data_dir: str | Path | None = None,
) -> tuple[ApiRun, dict[str, Any] | None, EvidenceRecord | None, dict[str, Any]]:
    """Capture the official AlphaFold viewer container (or Mol* canvas fallback)."""
    page_url = alphafold.entry_url_for_model(model_id)
    api = ApiRun(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name="AlphaFold",
        endpoint_name="capture_alphafold_entry_viewer",
        request_url=page_url,
        request_params={
            "model_entity_id": model_id,
            "uniprot_accession": accession,
            "allowlisted_hosts": [ALPHAFOLD_HOST],
            "retrieval_method": "official_web_element_capture",
            "headed": headed,
            "channel": channel,
        },
        success=False,
    )
    audit: dict[str, Any] = {
        "status": "unavailable",
        "model_entity_id": model_id,
        "source_page_url": page_url,
        "attempts": [],
        "capture_mode": CAPTURE_UNAVAILABLE,
    }
    try:
        from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # noqa: BLE001
        api.error_type = "playwright_unavailable"
        api.error_message = f"{type(exc).__name__}: {exc}"
        _save_api_run_failure(api, persist_db=persist_db)
        audit["reason"] = api.error_message
        return api, None, None, audit

    last_error: str | None = None
    for attempt in range(1, max_attempts + 1):
        attempt_info: dict[str, Any] = {"attempt": attempt, "headed": headed}
        browser = None
        context = None
        try:
            with sync_playwright() as pw:
                launch_kwargs: dict[str, Any] = {"headless": not headed}
                if channel:
                    launch_kwargs["channel"] = channel
                if user_data_dir:
                    context = pw.chromium.launch_persistent_context(
                        user_data_dir=str(user_data_dir),
                        **launch_kwargs,
                        viewport={"width": 1400, "height": 1600},
                        user_agent=(
                            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                            "AppleWebKit/537.36 (KHTML, like Gecko) "
                            "Chrome/122.0.0.0 Safari/537.36"
                        ),
                    )
                    page = context.new_page()
                else:
                    browser = pw.chromium.launch(**launch_kwargs)
                    page = browser.new_page(
                        viewport={"width": 1400, "height": 1600},
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
                page.goto(page_url, wait_until="domcontentloaded", timeout=45_000)
                host = urlparse(page.url).hostname or ""
                if host not in _ALPHAFOLD_ALLOWED_HOSTS:
                    raise ValueError(f"Unexpected AlphaFold host after navigation: {page.url}")
                if not _is_allowlisted_alphafold_entry_url(page.url, model_id=model_id):
                    raise ValueError(
                        f"AlphaFold navigation left the selected entry page: {page.url}"
                    )
                body_text = page.locator("body").inner_text(timeout=8_000)
                lower = body_text.lower()
                if "access denied" in lower or "incorrectly blocked" in lower:
                    raise ValueError(
                        "AlphaFold entry page returned Access Denied to browser capture"
                    )
                for spinner in (
                    "text=/loading/i",
                    "[class*='spinner']",
                    "[class*='Loading']",
                ):
                    loc = page.locator(spinner)
                    if loc.count() > 0:
                        try:
                            loc.first.wait_for(state="hidden", timeout=15_000)
                        except Exception:
                            pass
                try:
                    page.wait_for_selector(
                        ".summary-molstar-container, .structure-container, canvas",
                        state="attached",
                        timeout=30_000,
                    )
                except Exception:
                    pass
                try:
                    page.wait_for_selector(
                        ".summary-molstar-container, canvas",
                        state="visible",
                        timeout=20_000,
                    )
                except Exception:
                    pass
                try:
                    page.wait_for_load_state("networkidle", timeout=12_000)
                except Exception:
                    pass
                page.wait_for_timeout(3_000)
                if not _is_allowlisted_alphafold_entry_url(page.url, model_id=model_id):
                    raise ValueError(
                        f"AlphaFold URL drifted away from selected entry: {page.url}"
                    )

                # Prefer the visible Mol* structure viewport; never capture PAE/heatmap.
                selectors = [
                    ".summary-molstar-container",
                    "[class*='summary-molstar']",
                    ".structure-container",
                    "[class*='structure-container']",
                    "[class*='molstarView']",
                    "[class*='molstar']",
                    "canvas",
                ]
                locator = None
                selector_used = None
                for selector in selectors:
                    candidate = page.locator(selector)
                    count = candidate.count()
                    if count <= 0:
                        continue
                    for index in range(count):
                        node = candidate.nth(index)
                        try:
                            class_name = (
                                node.get_attribute("class") or ""
                            ).lower()
                        except Exception:
                            class_name = ""
                        if any(
                            token in class_name
                            for token in ("pae", "heatmap", "axis-box", "img-box")
                        ):
                            continue
                        try:
                            if not node.is_visible():
                                continue
                        except Exception:
                            continue
                        try:
                            node.scroll_into_view_if_needed(timeout=5_000)
                        except Exception:
                            pass
                        box = node.bounding_box()
                        if not box:
                            continue
                        if box["width"] < 280 or box["height"] < 240:
                            continue
                        locator = node
                        selector_used = selector
                        break
                    if locator is not None:
                        break
                if locator is None or selector_used is None:
                    raise ValueError("AlphaFold viewer container/canvas not found")

                png = locator.screenshot(type="png")
                validation = _validate_nonblank_image(png)
                ok, reason = _capture_passes_quality(validation, png)
                attempt_info.update(
                    {
                        "selector": selector_used,
                        "width": validation.get("width"),
                        "height": validation.get("height"),
                        "byte_size": validation.get("byte_size", len(png)),
                        "pixel_variance": validation.get("pixel_variance"),
                        "final_url": page.url,
                    }
                )
                if not ok:
                    raise ValueError(reason or "capture quality gate failed")

                artifact, meta = _persist_artifact_bytes(
                    dossier_run_id=dossier_run_id,
                    source_name="AlphaFold",
                    content=png,
                    extension="png",
                    artifact_type="png",
                    filename_hint=f"alphafold-{slugify(accession)}-viewer",
                    settings=settings,
                    api_run=api,
                    persist_db=persist_db,
                    notes={
                        "artifact_class": "derived_capture",
                        "artifact_origin": "alphafold_official_entry_viewer",
                        "artifact_role": "section_1d_human_structure_capture",
                        "retrieval_method": "official_web_element_capture",
                        "source_page_url": page_url,
                        "model_entity_id": model_id,
                        "uniprot_accession": accession,
                        "taxon_id": taxon_id,
                        "model_version": model_version,
                        "dom_selector": selector_used,
                        "viewport": {"width": 1400, "height": 1600},
                        "parent_raw_artifact_ids": list(parent_raw_artifact_ids),
                        "presentation_item_key": (
                            f"alphafold-human-{accession.lower()}"
                        ),
                        "headed": headed,
                        "channel": channel,
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
                    "model_entity_id": model_id,
                    "uniprot_accession": accession,
                    "taxon_id": taxon_id,
                    "model_version": model_version,
                    "source_page_url": page_url,
                    "relative_path": meta.get("relative_path"),
                    "media_type": meta.get("media_type") or "image/png",
                    "width": meta.get("width"),
                    "height": meta.get("height"),
                    "sha256": meta.get("expected_sha256") or meta.get("content_hash"),
                    "byte_size": meta.get("byte_size"),
                    "artifact_class": "derived_capture",
                    "retrieval_method": "official_web_element_capture",
                    "dom_selector": selector_used,
                    "parent_raw_artifact_ids": list(parent_raw_artifact_ids),
                    "presentation_item_key": f"alphafold-human-{accession.lower()}",
                    "figure_raw_artifact_id": artifact.id,
                }
                _write_capture_cache_index(settings=settings, model_id=model_id, payload=value)
                rec = _record(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    fact_type="alphafold_official_viewer_capture",
                    key=f"viewer-{accession}",
                    value=value,
                    display_text=(
                        f"{gene_symbol} official AlphaFold viewer capture "
                        f"for {model_id}."
                    ),
                    raw_artifact_id=artifact.id,
                    api_run_id=api.id,
                    taxon_id=taxon_id,
                    species="human",
                )
                audit = {
                    "status": "success",
                    "model_entity_id": model_id,
                    "source_page_url": page_url,
                    "selector": selector_used,
                    "attempts": audit.get("attempts", []) + [attempt_info],
                    "quality": {
                        "width": meta.get("width"),
                        "height": meta.get("height"),
                        "byte_size": meta.get("byte_size"),
                    },
                    "capture_mode": CAPTURE_FRESH,
                }
                return api, meta, rec, audit
        except PlaywrightTimeoutError as exc:
            last_error = f"timeout: {exc}"
            attempt_info["error"] = last_error
        except Exception as exc:  # noqa: BLE001
            last_error = f"{type(exc).__name__}: {exc}"
            attempt_info["error"] = last_error
        finally:
            audit.setdefault("attempts", []).append(attempt_info)
            if context is not None:
                try:
                    context.close()
                except Exception:
                    pass
            if browser is not None:
                try:
                    browser.close()
                except Exception:
                    pass

    api.error_type = "capture_failed"
    api.error_message = last_error or "AlphaFold viewer capture failed"
    _save_api_run_failure(api, persist_db=persist_db)
    audit["status"] = "unavailable"
    audit["reason"] = api.error_message
    audit["capture_mode"] = CAPTURE_UNAVAILABLE
    return api, None, None, audit


def obtain_human_viewer_capture(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    model_id: str,
    accession: str,
    taxon_id: int,
    model_version: int | None,
    parent_raw_artifact_ids: Sequence[str],
    evidence_records: Sequence[EvidenceRecord],
    raw_artifacts: Sequence[Any] | None,
    settings: Settings,
    persist_db: bool,
    skip_viewer_capture: bool = False,
) -> tuple[list[ApiRun], list[dict[str, Any]], EvidenceRecord | None, dict[str, Any]]:
    """Reuse a validated prior official capture, else attempt a fresh capture."""
    reused, reuse_diags = find_reusable_official_capture(
        evidence_records=evidence_records,
        raw_artifacts=raw_artifacts,
        model_id=model_id,
        accession=accession,
        model_version=model_version,
        settings=settings,
    )
    if reused is not None:
        value = {
            **reused,
            "status": "success",
            "reuse": True,
            "capture_mode": CAPTURE_REUSED,
            "presentation_item_key": f"alphafold-human-{accession.lower()}",
        }
        rec = _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            fact_type="alphafold_official_viewer_capture",
            key=f"viewer-{accession}",
            value=value,
            display_text=(
                f"{gene_symbol} reused official AlphaFold viewer capture "
                f"for {model_id}."
            ),
            raw_artifact_id=str(reused.get("figure_raw_artifact_id") or "") or None,
            taxon_id=taxon_id,
            species="human",
        )
        audit = {
            "status": "success",
            "model_entity_id": model_id,
            "source_page_url": reused.get("source_page_url"),
            "capture_mode": CAPTURE_REUSED,
            "reuse_diagnostics": reuse_diags,
            "quality": {
                "width": reused.get("width"),
                "height": reused.get("height"),
                "byte_size": reused.get("byte_size"),
            },
        }
        return [], [], rec, audit

    if skip_viewer_capture:
        return [], [], None, {
            "status": "skipped",
            "reason": "skip_viewer_capture=True",
            "capture_mode": CAPTURE_UNAVAILABLE,
            "reuse_diagnostics": reuse_diags,
        }

    # Prefer installed Chrome when available; AFDB often blocks stock headless shell.
    api, meta, rec, audit = _capture_human_viewer(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        model_id=model_id,
        accession=accession,
        taxon_id=taxon_id,
        model_version=model_version,
        parent_raw_artifact_ids=parent_raw_artifact_ids,
        settings=settings,
        persist_db=persist_db,
        channel="chrome",
    )
    audit["reuse_diagnostics"] = reuse_diags
    metas = [meta] if meta else []
    return [api], metas, rec, audit


def evaluate_section_1d_reference_genes_acceptance(
    *,
    gene_symbol: str,
    section_status: dict[str, Any] | None,
    presentation_blocks: Sequence[Any] | None,
    html: str | None,
    pdf_path: Path | None,
    selected_section_keys: Sequence[str],
) -> list[str]:
    """Return failure reasons for the Section 1d visual acceptance profile."""
    reasons: list[str] = []
    status = section_status or {}
    slots = list(status.get("species_slots") or [])
    by_species = {str(s.get("species_key") or "").lower(): s for s in slots}
    for species_key in ("human", "rat", "mouse"):
        slot = by_species.get(species_key)
        if slot is None:
            reasons.append(f"{species_key} species slot missing")
            continue
        if slot.get("status") != STATUS_SELECTED:
            reasons.append(f"{species_key} model not selected ({slot.get('status')})")
        if not slot.get("entry_url"):
            reasons.append(f"{species_key} AlphaFold link missing")

    viewer = status.get("viewer_capture") or (status.get("audit") or {}).get(
        "viewer_capture"
    ) or {}
    if isinstance(status.get("audit"), dict) and not viewer:
        viewer = status["audit"].get("viewer_capture") or {}
    capture_mode = viewer.get("capture_mode")
    if viewer.get("status") != "success":
        reasons.append("human official viewer capture missing or failed")
    elif capture_mode not in {CAPTURE_FRESH, CAPTURE_REUSED}:
        reasons.append(f"unexpected capture mode: {capture_mode}")

    blocks = list(presentation_blocks or [])
    roles = [str(getattr(b, "presentation_role", None) or "") for b in blocks]
    if "section_1d_human_structure_capture" not in roles:
        reasons.append("human structure image block missing")
    if "section_1d_confidence_legend" not in roles:
        reasons.append("confidence legend block missing")
    if roles.count("section_1d_species_link") < 3:
        reasons.append("rat or mouse AlphaFold link missing")
    if any(
        str(getattr(b, "presentation_role", None) or "") == "section_1d_species_status"
        and VISUALIZATION_UNAVAILABLE_TEXT in str(getattr(b, "text", "") or "")
        for b in blocks
    ):
        reasons.append("visualization-unavailable status present")
    if any(
        str(getattr(b, "presentation_role", None) or "") == "section_1d_species_status"
        for b in blocks
    ):
        reasons.append("species unavailable/status line present")

    # Order: Human link → visual (capture+legend) → Rat → Mouse.
    link_indexes = [
        i for i, role in enumerate(roles) if role == "section_1d_species_link"
    ]
    capture_i = next(
        (i for i, role in enumerate(roles) if role == "section_1d_human_structure_capture"),
        None,
    )
    legend_i = next(
        (i for i, role in enumerate(roles) if role == "section_1d_confidence_legend"),
        None,
    )
    if len(link_indexes) >= 3 and capture_i is not None and legend_i is not None:
        human_i, rat_i, mouse_i = link_indexes[0], link_indexes[1], link_indexes[2]
        if not (human_i < capture_i < legend_i < rat_i < mouse_i):
            reasons.append("block order is not Human → visual → Rat → Mouse")

    text_blob = " ".join(str(getattr(b, "text", "") or "") for b in blocks)
    for mid in (
        str(by_species.get("human", {}).get("model_entity_id") or ""),
        str(by_species.get("rat", {}).get("model_entity_id") or ""),
        str(by_species.get("mouse", {}).get("model_entity_id") or ""),
    ):
        if mid and mid in text_blob:
            reasons.append(f"model id appears in visible prose: {mid}")

    html_text = html or ""
    if html_text:
        if "section-1d-visual-table" not in html_text:
            reasons.append("section-1d-visual-table missing from HTML")
        if "section-1d-confidence-legend" not in html_text:
            reasons.append("confidence key missing from HTML")
        if "pLDDT" not in html_text:
            reasons.append("explanatory blurb missing from HTML")
        if "section-1d-visual-row" in html_text and "section-1d-visual-table" not in html_text:
            reasons.append("flex visual-row used without visual table")
        # Standalone legend: legend appears before any structure image.
        legend_pos = html_text.find("section-1d-confidence-legend")
        img_pos = html_text.find("section-1d-human-structure-capture")
        if legend_pos >= 0 and (img_pos < 0 or legend_pos < img_pos):
            reasons.append("standalone confidence legend without preceding structure image")
        if "pae" in html_text.lower() and "pae image" in html_text.lower():
            reasons.append("PAE image present in polished HTML")

        if "1c" in selected_section_keys and "1d" in selected_section_keys:
            # 1d should appear after the PDB official image within the same page segment.
            pdb_pos = html_text.find("section-1c-pdb-official-image")
            d_pos = html_text.find("subsection-d")
            if pdb_pos < 0 or d_pos < 0 or d_pos < pdb_pos:
                reasons.append("assembled 1d does not follow 1c PDB image")

    if pdf_path is not None and Path(pdf_path).is_file():
        try:
            import fitz

            with fitz.open(str(pdf_path)) as doc:
                page_count = doc.page_count
            keys = list(selected_section_keys)
            if keys == ["1d"] and page_count != 1:
                reasons.append(
                    f"focused 1d output is {page_count} pages; expected one Letter page"
                )
        except Exception as exc:  # noqa: BLE001
            reasons.append(f"pdf page-count check failed: {exc}")

    _ = gene_symbol  # reserved for gene-specific diagnostics
    return reasons


def node_generate_section_1d_derived_artifacts(
    state: DossierState,
    *,
    settings: Settings | None = None,
    persist_db: bool = True,
    transient: WorkflowTransientContext | None = None,
    skip_viewer_capture: bool = False,
) -> DossierState:
    """Sole Section 1d network owner: fetch×species, select, capture, evidence."""
    if state.get("run_type") != "section_bundle" or "1d" not in (
        state.get("selected_section_keys") or []
    ):
        return state

    cfg = settings or get_settings()
    run_id = state["dossier_run_id"]
    gene = state["gene_symbol"]
    _ = transient  # reserved for future binary handoff; captures persist immediately
    evidence = list(state.get("evidence_records") or [])
    api_runs = list(state.get("api_runs") or [])
    raw_meta = list(state.get("raw_artifacts") or [])
    errors = list(state.get("errors") or [])
    coverage_extra = list(state.get("coverage") or [])

    identities = resolve_species_identity(
        gene_symbol=gene,
        gene_ids=dict(state.get("gene_ids") or {}),
        evidence_records=evidence,
    )
    slots: list[dict[str, Any]] = []
    audit: dict[str, Any] = {
        "species_slots": [],
        "selection_diagnostics": {},
        "viewer_capture": {},
        "forbidden_paths": {
            "pymol_invoked": False,
            "chimera_invoked": False,
            "local_coordinate_projection": False,
            "pae_image_used_as_structure": False,
            "llm_invoked": False,
        },
        "network_owner": "section_1d",
    }

    human_selected: dict[str, Any] | None = None
    human_parent_raw_ids: list[str] = []

    for identity in identities:
        species_key = identity["species_key"]
        label = identity["species_label"]
        taxon_id = identity["taxon_id"]
        symbol = identity["display_symbol"]
        accession = identity["accession"]
        slot: dict[str, Any] = {
            "species_key": species_key,
            "species_label": label,
            "taxon_id": taxon_id,
            "display_symbol": symbol,
            "accession": accession,
            "status": STATUS_ACCESSION_UNAVAILABLE,
            "model_entity_id": None,
            "entry_url": None,
            "model_version": None,
            "message": None,
            "presentation_item_key": f"alphafold-{species_key}-unavailable",
        }

        if not accession:
            slot["status"] = STATUS_ACCESSION_UNAVAILABLE
            slot["message"] = unavailable_visible_text(
                species_label=label, display_symbol=symbol, status=slot["status"]
            )
            coverage_extra.append(
                SourceCoverageResult(
                    dossier_run_id=run_id,
                    source_name="AlphaFold",
                    status=SourceStatus.skipped,
                    evidence_record_count=0,
                    error_message=(
                        f"{label}: UniProt accession not available for Section 1d"
                    ),
                    report_sections_supported=["AlphaFold protein structure prediction"],
                    notes=f"species={species_key}; no ApiRun created",
                )
            )
            slots.append(slot)
            audit["species_slots"].append(dict(slot))
            continue

        tr = alphafold.fetch_prediction(
            accession, gene_symbol=gene, settings=cfg
        )
        # Preserve response payload in transient for debugging when available.
        if transient is not None and tr.data is not None:
            transient.put(run_id, f"alphafold-{species_key}-{accession}", tr.data)

        api, meta = _persist_tool_result_json(
            tr=tr,
            dossier_run_id=run_id,
            gene_symbol=gene,
            settings=cfg,
            persist_db=persist_db,
            filename_hint=f"alphafold-{accession}-prediction",
        )
        # Failed requests still produce an ApiRun (no fabricated success artifact).
        if not tr.success:
            # _persist_tool_result_json already saved failure ApiRun without artifact.
            api_runs.append(api)
            slot["status"] = STATUS_REQUEST_FAILED
            slot["message"] = unavailable_visible_text(
                species_label=label, display_symbol=symbol, status=slot["status"]
            )
            coverage_extra.append(
                SourceCoverageResult(
                    dossier_run_id=run_id,
                    source_name="AlphaFold",
                    status=SourceStatus.failed,
                    evidence_record_count=0,
                    error_message=tr.error_message or tr.error_type,
                    report_sections_supported=["AlphaFold protein structure prediction"],
                    notes=f"species={species_key}; accession={accession}",
                )
            )
            slots.append(slot)
            audit["species_slots"].append(dict(slot))
            continue

        api_runs.append(api)
        if meta:
            raw_meta.append(meta)
            if meta.get("id"):
                if species_key == "human":
                    human_parent_raw_ids.append(str(meta["id"]))

        predictions = []
        if isinstance(tr.data, dict):
            predictions = tr.data.get("predictions") or []
        selected, sel_diags = alphafold.select_canonical_monomer_prediction(
            predictions,
            accession,
            expected_taxon_id=taxon_id,
        )
        audit["selection_diagnostics"][species_key] = sel_diags

        if selected is None:
            slot["status"] = STATUS_MODEL_ABSENT
            slot["message"] = unavailable_visible_text(
                species_label=label, display_symbol=symbol, status=slot["status"]
            )
            coverage_extra.append(
                SourceCoverageResult(
                    dossier_run_id=run_id,
                    source_name="AlphaFold",
                    status=SourceStatus.skipped,
                    evidence_record_count=0,
                    error_message=(
                        f"No hard-qualified canonical monomer for {accession}"
                    ),
                    report_sections_supported=["AlphaFold protein structure prediction"],
                    notes=f"species={species_key}",
                )
            )
            slots.append(slot)
            audit["species_slots"].append(dict(slot))
            continue

        summary = alphafold.summarize_prediction(selected)
        mid = summary.get("model_entity_id")
        entry = summary.get("entry_url") or (
            alphafold.entry_url_for_model(str(mid)) if mid else None
        )
        version = summary.get("latest_version")
        try:
            version_int = int(version) if version is not None else None
        except (TypeError, ValueError):
            version_int = None

        slot.update(
            {
                "status": STATUS_SELECTED,
                "model_entity_id": mid,
                "entry_url": entry,
                "model_version": version_int,
                "presentation_item_key": (
                    f"alphafold-{species_key}-{str(accession).lower()}"
                ),
                "message": None,
            }
        )
        value = {
            "availability_status": STATUS_SELECTED,
            "species_key": species_key,
            "species_label": label,
            "display_symbol": symbol,
            "uniprot_accession": accession,
            "taxon_id": taxon_id,
            "model_entity_id": mid,
            "entry_url": entry,
            "model_version": version_int,
            "global_metric_value": summary.get("global_metric_value"),
            "fraction_plddt_very_high": summary.get("fraction_plddt_very_high"),
            "fraction_plddt_confident": summary.get("fraction_plddt_confident"),
            "fraction_plddt_low": summary.get("fraction_plddt_low"),
            "fraction_plddt_very_low": summary.get("fraction_plddt_very_low"),
            "caveat": (
                "AlphaFold structure is a prediction, not an experimental model."
            ),
            "presentation_item_key": slot["presentation_item_key"],
            # Explicitly do not require paeImageUrl.
            "pae_image_url": summary.get("pae_image_url"),
        }
        rec = _record(
            dossier_run_id=run_id,
            gene_symbol=gene,
            fact_type="alphafold_species_prediction",
            key=f"{species_key}-{accession}",
            value=value,
            display_text=(
                f"{label} {symbol} AlphaFold prediction {mid} "
                f"(UniProt {accession})."
            ),
            raw_artifact_id=str(meta.get("id")) if meta and meta.get("id") else None,
            api_run_id=api.id,
            taxon_id=taxon_id,
            organism=summary.get("organism_scientific_name"),
            species=species_key,
        )
        _append_evidence(evidence, rec, persist_db=persist_db)
        if species_key == "human":
            human_selected = {
                "accession": accession,
                "model_entity_id": mid,
                "model_version": version_int,
                "taxon_id": taxon_id,
            }
        slots.append(slot)
        audit["species_slots"].append(dict(slot))

    if human_selected:
        cap_apis, cap_metas, cap_rec, cap_audit = obtain_human_viewer_capture(
            dossier_run_id=run_id,
            gene_symbol=gene,
            model_id=str(human_selected["model_entity_id"]),
            accession=str(human_selected["accession"]),
            taxon_id=int(human_selected["taxon_id"]),
            model_version=human_selected.get("model_version"),
            parent_raw_artifact_ids=human_parent_raw_ids,
            evidence_records=evidence,
            raw_artifacts=raw_meta,
            settings=cfg,
            persist_db=persist_db,
            skip_viewer_capture=skip_viewer_capture,
        )
        api_runs.extend(cap_apis)
        raw_meta.extend(cap_metas)
        if cap_rec is not None:
            _append_evidence(evidence, cap_rec, persist_db=persist_db)
        else:
            audit["viewer_capture_status"] = STATUS_CAPTURE_FAILED
            if not cap_audit.get("capture_mode"):
                cap_audit["capture_mode"] = CAPTURE_UNAVAILABLE
        audit["viewer_capture"] = cap_audit
        audit["capture_mode"] = cap_audit.get("capture_mode") or CAPTURE_UNAVAILABLE
    else:
        audit["viewer_capture"] = {
            "status": "unavailable",
            "capture_mode": CAPTURE_UNAVAILABLE,
            "reason": "human model not selected",
        }
        audit["capture_mode"] = CAPTURE_UNAVAILABLE

    selected_count = sum(1 for s in slots if s.get("status") == STATUS_SELECTED)
    capture_status = (audit.get("viewer_capture") or {}).get("status")
    rendering_status = {
        "human": next(
            (s["status"] for s in slots if s["species_key"] == "human"),
            STATUS_ACCESSION_UNAVAILABLE,
        ),
        "rat": next(
            (s["status"] for s in slots if s["species_key"] == "rat"),
            STATUS_ACCESSION_UNAVAILABLE,
        ),
        "mouse": next(
            (s["status"] for s in slots if s["species_key"] == "mouse"),
            STATUS_ACCESSION_UNAVAILABLE,
        ),
        "viewer_capture": capture_status,
        "capture_mode": audit.get("capture_mode") or CAPTURE_UNAVAILABLE,
        "visualization": (
            "success"
            if capture_status == "success"
            else STATUS_VISUALIZATION_UNAVAILABLE
        ),
        "overall": "success"
        if selected_count == 3 and capture_status == "success"
        else ("partial" if selected_count else "empty"),
    }

    return {
        **state,
        "evidence_records": evidence,
        "api_runs": api_runs,
        "raw_artifacts": raw_meta,
        "errors": errors,
        "coverage": coverage_extra,
        "section_1d_status": {
            "species_slots": slots,
            "rendering_status": rendering_status,
            "audit": audit,
        },
    }


__all__ = [
    "STATUS_ACCESSION_UNAVAILABLE",
    "STATUS_MODEL_ABSENT",
    "STATUS_REQUEST_FAILED",
    "STATUS_SELECTED",
    "STATUS_CAPTURE_FAILED",
    "STATUS_VISUALIZATION_UNAVAILABLE",
    "CAPTURE_FRESH",
    "CAPTURE_REUSED",
    "CAPTURE_UNAVAILABLE",
    "VISUALIZATION_UNAVAILABLE_TEXT",
    "SPECIES_ORDER",
    "SpeciesSlot",
    "resolve_species_identity",
    "unavailable_visible_text",
    "find_reusable_official_capture",
    "obtain_human_viewer_capture",
    "evaluate_section_1d_reference_genes_acceptance",
    "node_generate_section_1d_derived_artifacts",
]
