"""Source coverage reporting.

The platform must never silently omit a source. This module:

1. Builds a :class:`SourceCoverageResult` for every entry in the source registry
2. Applies runtime outcomes (success / failed / partial / ...) when provided
3. Writes markdown + JSON reports under ``data/outputs/``
4. Optionally persists coverage rows via ``db.save_source_coverage``

Coverage status vocabulary matches :class:`~gene_dossier.models.SourceStatus`.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from sqlmodel import Session

from .config import Settings, get_settings
from .db import save_source_coverage
from .models import SourceCoverageResult, SourceStatus
from .source_registry import SourceDefinition, get_all_sources


def _missing_required_keys(source: SourceDefinition, settings: Settings) -> list[str]:
    """Return required env keys that are absent for ``source``."""
    return [key for key in source.required_keys if not settings.has_key(key)]


def initial_status_for_source(
    source: SourceDefinition, settings: Settings | None = None
) -> SourceStatus:
    """Resolve the starting coverage status for a registry source.

    Priority:
    1. Required API keys missing → ``requires_key``
    2. Registry ``manual`` / ``deferred`` (and other non-key defaults) → keep as-is
    3. Registry ``requires_key`` but keys *are* present → ``not_implemented`` until run
    4. Client not implemented → ``not_implemented``
    """
    cfg = settings or get_settings()
    missing = _missing_required_keys(source, cfg)
    if missing:
        return SourceStatus.requires_key

    if source.default_status is SourceStatus.requires_key:
        # Keys are present (checked above); wait for the client to run.
        return SourceStatus.not_implemented

    if source.default_status != SourceStatus.not_implemented:
        return source.default_status

    return SourceStatus.not_implemented


def build_coverage_for_registry(
    dossier_run_id: str,
    *,
    settings: Settings | None = None,
    sources: Iterable[SourceDefinition] | None = None,
) -> list[SourceCoverageResult]:
    """Build a coverage row for every registered source (baseline before/during a run).

    Ensures the full source map appears in the report even if a source was never called.
    """
    cfg = settings or get_settings()
    results: list[SourceCoverageResult] = []
    for src in sources if sources is not None else get_all_sources():
        status = initial_status_for_source(src, cfg)
        missing = _missing_required_keys(src, cfg)
        notes = src.notes
        error_message = None
        if status is SourceStatus.requires_key and missing:
            error_message = f"missing required key(s): {', '.join(missing)}"
            notes = (notes + " | " if notes else "") + error_message

        results.append(
            SourceCoverageResult(
                dossier_run_id=dossier_run_id,
                source_name=src.name,
                status=status,
                evidence_record_count=0,
                error_message=error_message,
                report_sections_supported=list(src.report_sections),
                notes=notes,
            )
        )
    return results


def apply_coverage_updates(
    baseline: list[SourceCoverageResult],
    updates: Iterable[SourceCoverageResult],
) -> list[SourceCoverageResult]:
    """Merge runtime ``updates`` into ``baseline`` by ``source_name`` (case-insensitive).

    Sources present only in ``updates`` are appended. Order follows the baseline, then
    any extra update-only sources.
    """
    by_name = {r.source_name.lower(): r for r in baseline}
    extras: list[SourceCoverageResult] = []
    for upd in updates:
        key = upd.source_name.lower()
        if key in by_name:
            by_name[key] = upd
        else:
            extras.append(upd)

    ordered: list[SourceCoverageResult] = []
    seen: set[str] = set()
    for base in baseline:
        key = base.source_name.lower()
        ordered.append(by_name[key])
        seen.add(key)
    for extra in extras:
        if extra.source_name.lower() not in seen:
            ordered.append(extra)
    return ordered


def summarize_coverage(results: list[SourceCoverageResult]) -> dict[str, Any]:
    """Return counts by status and total source count."""
    counts = Counter(r.status.value for r in results)
    return {
        "total": len(results),
        "by_status": dict(sorted(counts.items())),
    }


def coverage_to_jsonable(results: list[SourceCoverageResult]) -> list[dict[str, Any]]:
    """Serialize coverage rows to JSON-friendly dicts."""
    return [r.model_dump(mode="json") for r in results]


def render_coverage_markdown(
    results: list[SourceCoverageResult],
    *,
    dossier_run_id: str,
    gene_symbol: str | None = None,
) -> str:
    """Render a human-readable markdown coverage report."""
    summary = summarize_coverage(results)
    title_gene = f" — {gene_symbol}" if gene_symbol else ""
    lines: list[str] = [
        f"# Source coverage report{title_gene}",
        "",
        f"- **dossier_run_id:** `{dossier_run_id}`",
        f"- **total sources:** {summary['total']}",
        f"- **by status:** {summary['by_status']}",
        "",
        "| Source | Status | Evidence | Artifact | Sections | Notes / error |",
        "| --- | --- | ---: | --- | --- | --- |",
    ]

    # Look up priority for display when possible.
    registry = {s.name.lower(): s for s in get_all_sources()}

    for r in results:
        src = registry.get(r.source_name.lower())
        priority = src.priority.value if src else "?"
        artifact = r.raw_artifact_path or "—"
        # Keep table cells compact.
        if len(artifact) > 48:
            artifact = "…" + artifact[-47:]
        sections = ", ".join(r.report_sections_supported) if r.report_sections_supported else "—"
        if len(sections) > 40:
            sections = sections[:37] + "…"
        note = r.error_message or r.notes or "—"
        if len(note) > 60:
            note = note[:57] + "…"
        lines.append(
            f"| {r.source_name} (P{priority}) | `{r.status.value}` | "
            f"{r.evidence_record_count} | `{artifact}` | {sections} | {note} |"
        )

    lines.extend(["", "## Status legend", ""])
    for status in SourceStatus:
        lines.append(f"- `{status.value}`")
    lines.append("")
    return "\n".join(lines)


def write_coverage_report(
    results: list[SourceCoverageResult],
    dossier_run_id: str,
    *,
    gene_symbol: str | None = None,
    output_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Write ``{run_id}_source_coverage.md`` and ``.json``; return their paths.

    Paths:
      ``{output_dir}/{dossier_run_id}_source_coverage.md``
      ``{output_dir}/{dossier_run_id}_source_coverage.json``
    """
    settings = get_settings()
    out = Path(output_dir) if output_dir is not None else settings.output_path
    out.mkdir(parents=True, exist_ok=True)

    md_path = out / f"{dossier_run_id}_source_coverage.md"
    json_path = out / f"{dossier_run_id}_source_coverage.json"

    md_path.write_text(
        render_coverage_markdown(
            results, dossier_run_id=dossier_run_id, gene_symbol=gene_symbol
        ),
        encoding="utf-8",
    )
    payload = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "summary": summarize_coverage(results),
        "sources": coverage_to_jsonable(results),
    }
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return {"markdown": md_path, "json": json_path}


def persist_coverage(session: Session, results: Iterable[SourceCoverageResult]) -> int:
    """Save coverage rows to the database; return the number of rows written."""
    count = 0
    for result in results:
        save_source_coverage(session, result)
        count += 1
    return count


def build_and_write_coverage(
    dossier_run_id: str,
    *,
    gene_symbol: str | None = None,
    updates: Iterable[SourceCoverageResult] | None = None,
    settings: Settings | None = None,
    output_dir: str | Path | None = None,
    session: Session | None = None,
) -> tuple[list[SourceCoverageResult], dict[str, Path]]:
    """Build registry baseline, apply optional updates, write reports, optionally persist.

    Returns ``(results, paths)`` where ``paths`` has ``markdown`` and ``json`` keys.
    """
    baseline = build_coverage_for_registry(dossier_run_id, settings=settings)
    results = apply_coverage_updates(baseline, updates or [])
    paths = write_coverage_report(
        results,
        dossier_run_id,
        gene_symbol=gene_symbol,
        output_dir=output_dir,
    )
    if session is not None:
        persist_coverage(session, results)
    return results, paths


__all__ = [
    "initial_status_for_source",
    "build_coverage_for_registry",
    "apply_coverage_updates",
    "summarize_coverage",
    "coverage_to_jsonable",
    "render_coverage_markdown",
    "write_coverage_report",
    "persist_coverage",
    "build_and_write_coverage",
]
