"""CHDI-style dossier rendering (presentation layer).

Assembles markdown (and optional JSON) reports from already-built
:class:`~gene_dossier.models.ReportSection`, :class:`~gene_dossier.models.Claim`,
:class:`~gene_dossier.models.VerificationResult`, and
:class:`~gene_dossier.models.SourceCoverageResult` objects.

No network I/O and no LLM calls — this module only formats what upstream layers
already produced.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import (
    Claim,
    ReportSection,
    SourceCoverageResult,
    SourceStatus,
    VerificationResult,
)
from gene_dossier.synthesis import CHDI_REPORT_SECTIONS, SynthesisResult

# Section statuses treated as having usable synthesized content.
_COMPLETED_SECTION_STATUSES = frozenset({"deterministic", "llm", "complete"})
_EMPTY_SECTION_STATUSES = frozenset({"empty"})
_DEFERRED_SECTION_STATUSES = frozenset({"deferred"})


def utc_now_iso() -> str:
    """UTC timestamp for report headers."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _status_value(status: SourceStatus | str) -> str:
    """Normalize a coverage status to its string value."""
    return status.value if isinstance(status, SourceStatus) else str(status)


def _is_success_status(status: SourceStatus | str) -> bool:
    """True if ``status`` represents successful coverage."""
    return _status_value(status) == SourceStatus.success.value


# --------------------------------------------------------------------------------------
# Classification helpers
# --------------------------------------------------------------------------------------
def classify_section_bucket(section: ReportSection) -> str:
    """Map a report section to a completion bucket name."""
    status = (section.status or "").strip().lower()
    if status in _COMPLETED_SECTION_STATUSES:
        return "completed"
    if status in _EMPTY_SECTION_STATUSES:
        return "missing"
    if status in _DEFERRED_SECTION_STATUSES:
        return "deferred"
    if status in {"partial", "draft"}:
        return "partial"
    # Unknown status with content → partial; without → missing.
    if (section.content_markdown or "").strip():
        return "partial"
    return "missing"


def summarize_sections(sections: Iterable[ReportSection]) -> dict[str, Any]:
    """Count sections by completion bucket and list names per bucket."""
    by_bucket: dict[str, list[str]] = defaultdict(list)
    for section in sections:
        bucket = classify_section_bucket(section)
        by_bucket[bucket].append(section.section_name)
    counts = {k: len(v) for k, v in sorted(by_bucket.items())}
    return {
        "total": sum(counts.values()),
        "by_bucket": counts,
        "names_by_bucket": {k: list(v) for k, v in sorted(by_bucket.items())},
    }


def summarize_coverage_by_status(
    coverage: Iterable[SourceCoverageResult] | None,
) -> dict[str, list[str]]:
    """Group source names by :class:`SourceStatus` value."""
    grouped: dict[str, list[str]] = defaultdict(list)
    if not coverage:
        return {}
    for row in coverage:
        grouped[_status_value(row.status)].append(row.source_name)
    return {k: list(v) for k, v in sorted(grouped.items())}


def summarize_verification(
    results: Iterable[VerificationResult] | None,
) -> dict[str, Any]:
    """Count verification verdicts and list non-pass items."""
    rows = list(results or [])
    counts = Counter(r.verdict for r in rows)
    flagged = [
        {
            "claim_id": r.claim_id,
            "verdict": r.verdict,
            "needs_human_review": r.needs_human_review,
            "reason": r.reason,
        }
        for r in rows
        if r.verdict != "pass" or r.needs_human_review
    ]
    return {
        "total": len(rows),
        "by_verdict": dict(sorted(counts.items())),
        "flagged_count": len(flagged),
        "flagged": flagged,
    }


# --------------------------------------------------------------------------------------
# Meta-section content
# --------------------------------------------------------------------------------------
def render_missing_sources_markdown(
    coverage: Iterable[SourceCoverageResult] | None,
) -> str:
    """Markdown for the Missing / deferred / manual sources section."""
    lines = [
        "## Missing / deferred / manual sources",
        "",
        "Sources that did not fully succeed in this run (never silently omitted):",
        "",
    ]
    coverage_list = list(coverage or [])
    rows = [r for r in coverage_list if not _is_success_status(r.status)]
    if not rows:
        lines.append(
            "_All registered sources reported `success`, or no coverage data was provided._"
        )
        lines.append("")
        return "\n".join(lines)

    preferred_order = [
        SourceStatus.failed,
        SourceStatus.requires_key,
        SourceStatus.manual,
        SourceStatus.deferred,
        SourceStatus.partial,
        SourceStatus.skipped,
        SourceStatus.not_implemented,
    ]
    preferred_values = {s.value for s in preferred_order}
    by_status: dict[str, list[SourceCoverageResult]] = defaultdict(list)
    for row in rows:
        by_status[_status_value(row.status)].append(row)

    for status in preferred_order:
        group = by_status.get(status.value) or []
        if not group:
            continue
        lines.append(f"### `{status.value}` ({len(group)})")
        lines.append("")
        for row in group:
            note = row.error_message or row.notes or ""
            note_bit = f" — {note}" if note else ""
            evidence = f"; evidence={row.evidence_record_count}"
            lines.append(f"- **{row.source_name}**{evidence}{note_bit}")
        lines.append("")

    for status_value, group in sorted(by_status.items()):
        if status_value in preferred_values:
            continue
        lines.append(f"### `{status_value}` ({len(group)})")
        lines.append("")
        for row in group:
            lines.append(f"- **{row.source_name}**")
        lines.append("")

    return "\n".join(lines)


def render_verification_warnings_markdown(
    verification_results: Iterable[VerificationResult] | None,
    claims: Iterable[Claim] | None = None,
) -> str:
    """Markdown for the Verification warnings section."""
    lines = [
        "## Verification warnings",
        "",
    ]
    results = list(verification_results or [])
    claim_by_id = {c.id: c for c in (claims or []) if c.id}

    flagged = [r for r in results if r.verdict != "pass" or r.needs_human_review]
    if not results:
        lines.append("_No verification results were provided for this run._")
        lines.append("")
        return "\n".join(lines)
    if not flagged:
        lines.append(
            f"_All {len(results)} claim(s) passed verification with no human-review flags._"
        )
        lines.append("")
        return "\n".join(lines)

    lines.append(
        f"{len(flagged)} of {len(results)} claim(s) need attention "
        "(fail / warning / human_review)."
    )
    lines.append("")
    lines.append("| Verdict | Review | Claim | Reason |")
    lines.append("| --- | --- | --- | --- |")

    for result in flagged:
        claim = claim_by_id.get(result.claim_id)
        text = (claim.claim_text if claim else result.claim_id) or result.claim_id
        text = " ".join(text.split())
        if len(text) > 80:
            text = text[:77] + "…"
        reason = (result.reason or "—").replace("|", "/")
        if len(reason) > 100:
            reason = reason[:97] + "…"
        review = "yes" if result.needs_human_review else "no"
        lines.append(f"| `{result.verdict}` | {review} | {text} | {reason} |")

    lines.append("")
    return "\n".join(lines)


def apply_meta_sections(
    sections: list[ReportSection],
    *,
    coverage: Iterable[SourceCoverageResult] | None = None,
    verification_results: Iterable[VerificationResult] | None = None,
    claims: Iterable[Claim] | None = None,
) -> list[ReportSection]:
    """Return sections with meta section bodies filled from coverage / verification."""
    out: list[ReportSection] = []
    for section in sections:
        name = section.section_name
        if name == "Missing / deferred / manual sources":
            out.append(
                section.model_copy(
                    update={
                        "content_markdown": render_missing_sources_markdown(coverage),
                        "status": "complete",
                    }
                )
            )
        elif name == "Verification warnings":
            out.append(
                section.model_copy(
                    update={
                        "content_markdown": render_verification_warnings_markdown(
                            verification_results, claims
                        ),
                        "status": "complete",
                    }
                )
            )
        else:
            out.append(section)
    return out


# --------------------------------------------------------------------------------------
# Full dossier markdown
# --------------------------------------------------------------------------------------
def _ordered_sections(sections: list[ReportSection]) -> list[ReportSection]:
    """Order by CHDI list; append any unknown section names at the end."""
    by_name: dict[str, ReportSection] = {}
    for section in sections:
        by_name[section.section_name] = section
    ordered: list[ReportSection] = []
    seen: set[str] = set()
    for name in CHDI_REPORT_SECTIONS:
        if name in by_name:
            ordered.append(by_name[name])
            seen.add(name)
    for section in sections:
        if section.section_name not in seen:
            ordered.append(section)
            seen.add(section.section_name)
    return ordered


def render_status_overview_markdown(
    *,
    sections: list[ReportSection],
    coverage: Iterable[SourceCoverageResult] | None = None,
    verification_results: Iterable[VerificationResult] | None = None,
) -> str:
    """Executive overview separating completed / missing / coverage buckets."""
    section_summary = summarize_sections(sections)
    coverage_groups = summarize_coverage_by_status(coverage)
    ver_summary = summarize_verification(verification_results)

    lines = [
        "## Status overview",
        "",
        "### Report sections",
        "",
    ]
    names_by_bucket = section_summary["names_by_bucket"]
    for bucket in ("completed", "partial", "missing", "deferred"):
        names = names_by_bucket.get(bucket) or []
        lines.append(f"- **{bucket}:** {len(names)}")
        if names and bucket != "completed":
            preview = ", ".join(names[:8])
            if len(names) > 8:
                preview += f", … (+{len(names) - 8})"
            lines.append(f"  - {preview}")
    lines.append("")

    lines.extend(["### Source coverage", ""])
    if not coverage_groups:
        lines.append("_No source coverage data provided._")
    else:
        for key in (
            "success",
            "partial",
            "failed",
            "requires_key",
            "manual",
            "deferred",
            "not_implemented",
            "skipped",
        ):
            names = coverage_groups.get(key) or []
            if not names and key not in coverage_groups:
                continue
            lines.append(f"- **{key}:** {len(names)}")
        for key, names in coverage_groups.items():
            if key in {
                "success",
                "partial",
                "failed",
                "requires_key",
                "manual",
                "deferred",
                "not_implemented",
                "skipped",
            }:
                continue
            lines.append(f"- **{key}:** {len(names)}")
    lines.append("")

    lines.extend(
        [
            "### Verification",
            "",
            f"- **claims verified:** {ver_summary['total']}",
            f"- **by verdict:** {ver_summary['by_verdict'] or '{}'}",
            f"- **flagged:** {ver_summary['flagged_count']}",
            "",
        ]
    )
    return "\n".join(lines)


def render_dossier_markdown(
    *,
    gene_symbol: str,
    dossier_run_id: str,
    sections: list[ReportSection],
    claims: Iterable[Claim] | None = None,
    coverage: Iterable[SourceCoverageResult] | None = None,
    verification_results: Iterable[VerificationResult] | None = None,
    synthesis_mode: str | None = None,
    synthesis_notes: Iterable[str] | None = None,
    generated_at: str | None = None,
) -> str:
    """Render the full CHDI-style gene dossier markdown document."""
    claim_list = list(claims or [])
    coverage_list = list(coverage or [])
    verification_list = list(verification_results or [])
    filled = apply_meta_sections(
        list(sections),
        coverage=coverage_list,
        verification_results=verification_list,
        claims=claim_list,
    )
    ordered = _ordered_sections(filled)
    stamp = generated_at or utc_now_iso()

    lines: list[str] = [
        f"# Gene dossier: {gene_symbol}",
        "",
        f"- **dossier_run_id:** `{dossier_run_id}`",
        f"- **gene_symbol:** `{gene_symbol}`",
        f"- **generated_at:** `{stamp}`",
    ]
    if synthesis_mode:
        lines.append(f"- **synthesis_mode:** `{synthesis_mode}`")
    lines.append(f"- **claims:** {len(claim_list)}")
    lines.append("")

    if synthesis_notes:
        lines.append("### Synthesis notes")
        lines.append("")
        for note in synthesis_notes:
            lines.append(f"- {note}")
        lines.append("")

    lines.append(
        render_status_overview_markdown(
            sections=ordered,
            coverage=coverage_list,
            verification_results=verification_list,
        )
    )

    lines.extend(["## Report sections", ""])
    for section in ordered:
        if section.section_name in {
            "Missing / deferred / manual sources",
            "Verification warnings",
        }:
            content = (section.content_markdown or "").rstrip() + "\n"
            lines.append(content)
            continue

        heading = f"## {section.section_name}"
        body = (section.content_markdown or "").strip()
        if body.startswith(f"**{section.section_name}"):
            lines.append(heading)
            lines.append("")
            lines.append(body)
        elif body:
            lines.append(heading)
            lines.append("")
            lines.append(body)
        else:
            lines.append(heading)
            lines.append("")
            lines.append(f"_No content for {section.section_name}._")
        lines.append("")
        if section.source_ids:
            cites = ", ".join(f"`{sid}`" for sid in section.source_ids[:40])
            if len(section.source_ids) > 40:
                cites += f", … (+{len(section.source_ids) - 40})"
            lines.append(f"_Section source_ids ({len(section.source_ids)}):_ {cites}")
            lines.append("")

    lines.extend(
        [
            "---",
            "",
            "_Provenance-first report: every claim must cite `source_id` values that "
            "resolve to evidence records. The LLM is never the source of truth._",
            "",
        ]
    )
    return "\n".join(lines)


def dossier_to_jsonable(
    *,
    gene_symbol: str,
    dossier_run_id: str,
    sections: list[ReportSection],
    claims: Iterable[Claim] | None = None,
    coverage: Iterable[SourceCoverageResult] | None = None,
    verification_results: Iterable[VerificationResult] | None = None,
    synthesis_mode: str | None = None,
    synthesis_notes: Iterable[str] | None = None,
    generated_at: str | None = None,
) -> dict[str, Any]:
    """JSON-serializable dossier payload (sidecar for the markdown report)."""
    claim_list = list(claims or [])
    coverage_list = list(coverage or [])
    verification_list = list(verification_results or [])
    filled = apply_meta_sections(
        list(sections),
        coverage=coverage_list,
        verification_results=verification_list,
        claims=claim_list,
    )
    ordered = _ordered_sections(filled)
    return {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "generated_at": generated_at or utc_now_iso(),
        "synthesis_mode": synthesis_mode,
        "synthesis_notes": list(synthesis_notes or []),
        "section_summary": summarize_sections(ordered),
        "coverage_by_status": summarize_coverage_by_status(coverage_list),
        "verification_summary": summarize_verification(verification_list),
        "sections": [s.model_dump(mode="json") for s in ordered],
        "claims": [c.model_dump(mode="json") for c in claim_list],
        "verification_results": [v.model_dump(mode="json") for v in verification_list],
        "coverage": [r.model_dump(mode="json") for r in coverage_list],
    }


# --------------------------------------------------------------------------------------
# Write helpers
# --------------------------------------------------------------------------------------
def write_dossier_report(
    *,
    gene_symbol: str,
    dossier_run_id: str,
    sections: list[ReportSection],
    claims: Iterable[Claim] | None = None,
    coverage: Iterable[SourceCoverageResult] | None = None,
    verification_results: Iterable[VerificationResult] | None = None,
    synthesis_mode: str | None = None,
    synthesis_notes: Iterable[str] | None = None,
    output_dir: str | Path | None = None,
    settings: Settings | None = None,
) -> dict[str, Path]:
    """Write ``{run_id}_report.md`` and ``{run_id}_report.json``; return paths.

    Paths:
      ``{output_dir}/{dossier_run_id}_report.md``
      ``{output_dir}/{dossier_run_id}_report.json``
    """
    cfg = settings or get_settings()
    out = Path(output_dir) if output_dir is not None else cfg.output_path
    out.mkdir(parents=True, exist_ok=True)

    claim_list = list(claims or [])
    coverage_list = list(coverage or [])
    verification_list = list(verification_results or [])

    stamp = utc_now_iso()
    md = render_dossier_markdown(
        gene_symbol=gene_symbol,
        dossier_run_id=dossier_run_id,
        sections=sections,
        claims=claim_list,
        coverage=coverage_list,
        verification_results=verification_list,
        synthesis_mode=synthesis_mode,
        synthesis_notes=synthesis_notes,
        generated_at=stamp,
    )
    payload = dossier_to_jsonable(
        gene_symbol=gene_symbol,
        dossier_run_id=dossier_run_id,
        sections=sections,
        claims=claim_list,
        coverage=coverage_list,
        verification_results=verification_list,
        synthesis_mode=synthesis_mode,
        synthesis_notes=synthesis_notes,
        generated_at=stamp,
    )

    md_path = out / f"{dossier_run_id}_report.md"
    json_path = out / f"{dossier_run_id}_report.json"
    md_path.write_text(md, encoding="utf-8")
    json_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return {"markdown": md_path, "json": json_path}


def render_synthesis_result(
    synthesis: SynthesisResult,
    *,
    coverage: Iterable[SourceCoverageResult] | None = None,
    verification_results: Iterable[VerificationResult] | None = None,
    output_dir: str | Path | None = None,
    settings: Settings | None = None,
    write: bool = True,
) -> tuple[str, dict[str, Path] | None]:
    """Render (and optionally write) a dossier from a :class:`SynthesisResult`."""
    coverage_list = list(coverage or [])
    verification_list = list(verification_results or [])
    claim_list = list(synthesis.claims or [])
    md = render_dossier_markdown(
        gene_symbol=synthesis.gene_symbol,
        dossier_run_id=synthesis.dossier_run_id,
        sections=synthesis.sections,
        claims=claim_list,
        coverage=coverage_list,
        verification_results=verification_list,
        synthesis_mode=synthesis.mode,
        synthesis_notes=synthesis.notes,
    )
    paths: dict[str, Path] | None = None
    if write:
        paths = write_dossier_report(
            gene_symbol=synthesis.gene_symbol,
            dossier_run_id=synthesis.dossier_run_id,
            sections=synthesis.sections,
            claims=claim_list,
            coverage=coverage_list,
            verification_results=verification_list,
            synthesis_mode=synthesis.mode,
            synthesis_notes=synthesis.notes,
            output_dir=output_dir,
            settings=settings,
        )
    return md, paths


__all__ = [
    "classify_section_bucket",
    "summarize_sections",
    "summarize_coverage_by_status",
    "summarize_verification",
    "render_missing_sources_markdown",
    "render_verification_warnings_markdown",
    "apply_meta_sections",
    "render_status_overview_markdown",
    "render_dossier_markdown",
    "dossier_to_jsonable",
    "write_dossier_report",
    "render_synthesis_result",
    "utc_now_iso",
]
