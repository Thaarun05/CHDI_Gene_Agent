#!/usr/bin/env python3
"""Smoke-test registered biomedical source clients.

Calls each dispatched client for a gene (default: SREBF2), soft-fails per
source, and writes a compact markdown + JSON smoke report. This is **not** a
full dossier pass — it only exercises live (or key-gated) client calls and
records success/failure without normalization, synthesis, or Rancho rendering.

Usage::

    python scripts/run_source_smoke_tests.py
    python scripts/run_source_smoke_tests.py --gene SREBF2 --priority A
    python scripts/run_source_smoke_tests.py --sources "NCBI Gene,Ensembl,UniProt"

Missing API keys degrade to ``requires_key``. Network/HTTP failures become
``failed``. The script always exits 0 unless an unexpected fatal error occurs
outside a source call (so CI can archive the report even when sources fail).
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import Settings, get_settings  # noqa: E402
from gene_dossier.models import ToolResult  # noqa: E402
from gene_dossier.source_registry import (  # noqa: E402
    SourceDefinition,
    SourcePriority,
    get_all_sources,
    get_source,
)
from gene_dossier.workflow import (  # noqa: E402
    CLIENT_DISPATCH,
    IDENTITY_SOURCES,
    _safe_call_client,
    extract_gene_ids_from_tool_result,
)

DEFAULT_GENE = "SREBF2"
LOGGER = logging.getLogger("run_source_smoke_tests")

_SECRET_QUERY_KEYS = {
    "api_key",
    "apikey",
    "access_key",
    "accesskey",
    "token",
    "auth",
    "authorization",
    "key",
}


def _redact_url(url: str | None) -> str:
    """Return a truncated URL with credential-bearing query params redacted."""
    if not url:
        return ""
    try:
        parts = urlsplit(url)
        secret_keys = {k.lower() for k in _SECRET_QUERY_KEYS}
        query: list[tuple[str, str]] = []
        for key, value in parse_qsl(parts.query, keep_blank_values=True):
            if key.lower() in secret_keys:
                query.append((key, "REDACTED"))
            else:
                query.append((key, value))
        redacted = urlunsplit(
            (parts.scheme, parts.netloc, parts.path, urlencode(query), parts.fragment)
        )
        return redacted[:240]
    except Exception:  # noqa: BLE001
        return "[unparseable_url_redacted]"


def _row_from_call(
    src: SourceDefinition,
    *,
    status: str,
    result: ToolResult | None,
    elapsed: float,
    missing: list[str],
) -> dict[str, Any]:
    return {
        "source_name": src.name,
        "priority": src.priority.value,
        "status": status,
        "elapsed_seconds": round(elapsed, 3),
        "endpoint_name": result.endpoint_name if result else None,
        "status_code": result.status_code if result else None,
        "success": bool(result.success) if result else False,
        "error_type": result.error_type
        if result
        else ("not_implemented" if src.name not in CLIENT_DISPATCH else None),
        "error_message": (
            result.error_message
            if result
            else (
                f"missing required key(s): {', '.join(missing)}"
                if missing
                else "no client dispatch"
            )
        ),
        "request_url": _redact_url(result.request_url if result else ""),
        "missing_keys": missing,
    }


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Smoke-test Gene Dossier source clients (default gene: SREBF2). "
            "Writes a markdown + JSON smoke report under the output directory."
        )
    )
    parser.add_argument(
        "--gene",
        default=DEFAULT_GENE,
        help=f"Gene symbol (default: {DEFAULT_GENE})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for smoke reports (default: settings.OUTPUT_DIR)",
    )
    parser.add_argument(
        "--sources",
        default=None,
        help='Comma-separated source names (e.g. "NCBI Gene,Ensembl,GTEx")',
    )
    parser.add_argument(
        "--priority",
        choices=["A", "B", "C"],
        default=None,
        help="Limit to registry priority A, B, or C",
    )
    parser.add_argument(
        "--identity-only",
        action="store_true",
        help="Only call NCBI Gene, Ensembl, and UniProt",
    )
    parser.add_argument(
        "--skip-identity-bootstrap",
        action="store_true",
        help="Do not call identity sources first to seed gene_ids",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable DEBUG logging",
    )
    return parser.parse_args(argv)


def _split_csv(raw: str | None) -> list[str] | None:
    if not raw:
        return None
    return [part.strip() for part in raw.split(",") if part.strip()]


def _db_kind(settings: Settings) -> str:
    url = settings.database_url or ""
    if url.startswith(("postgresql://", "postgresql+psycopg://")):
        return "postgres"
    if url.startswith("sqlite:"):
        return "sqlite"
    return "other"


def _select_sources(
    *,
    names: list[str] | None,
    priority: str | None,
    identity_only: bool,
) -> list[SourceDefinition]:
    if identity_only:
        selected = []
        for name in IDENTITY_SOURCES:
            src = get_source(name)
            if src is not None:
                selected.append(src)
        return selected

    all_sources = get_all_sources()
    if names:
        wanted = {n.lower() for n in names}
        selected = [s for s in all_sources if s.name.lower() in wanted]
        missing = wanted - {s.name.lower() for s in selected}
        if missing:
            LOGGER.warning("Unknown source name(s) ignored: %s", ", ".join(sorted(missing)))
        return selected

    if priority:
        p = SourcePriority(priority.upper())
        return [s for s in all_sources if s.priority == p]

    # Default: only sources with a workflow client dispatch entry.
    return [s for s in all_sources if s.name in CLIENT_DISPATCH]


def _missing_keys(src: SourceDefinition, settings: Settings) -> list[str]:
    return [k for k in src.required_keys if not settings.has_key(k)]


def _status_for_result(
    result: ToolResult | None,
    *,
    missing_keys: list[str],
    no_client: bool,
) -> str:
    if no_client:
        return "not_implemented"
    if missing_keys:
        return "requires_key"
    if result is None:
        return "failed"
    if result.error_type == "requires_key":
        return "requires_key"
    if result.error_type == "missing_identifier":
        return "skipped"
    if result.success:
        return "success"
    return "failed"


def _call_one(
    src: SourceDefinition,
    *,
    gene_symbol: str,
    gene_ids: dict[str, Any],
    settings: Settings,
) -> tuple[str, ToolResult | None, float, list[str]]:
    """Return (status, tool_result|None, elapsed_sec, missing_keys)."""
    missing = _missing_keys(src, settings)
    fn = CLIENT_DISPATCH.get(src.name)
    if fn is None:
        return "not_implemented", None, 0.0, missing
    if missing:
        stub = ToolResult(
            source_name=src.name,
            endpoint_name="smoke_skip",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="requires_key",
            error_message=f"missing required key(s): {', '.join(missing)}",
        )
        return "requires_key", stub, 0.0, missing

    started = time.perf_counter()
    result = _safe_call_client(
        src.name, fn, gene_symbol=gene_symbol, gene_ids=gene_ids, settings=settings
    )
    elapsed = time.perf_counter() - started
    status = _status_for_result(result, missing_keys=[], no_client=False)
    return status, result, elapsed, missing


def _render_markdown(
    *,
    gene_symbol: str,
    rows: list[dict[str, Any]],
    gene_ids: dict[str, Any],
) -> str:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    lines = [
        f"# Source smoke report — {gene_symbol}",
        "",
        f"- **generated_at:** {datetime.now(timezone.utc).isoformat()}",
        f"- **total sources tested:** {len(rows)}",
        f"- **by status:** {counts}",
        "",
        "## Resolved identifiers",
        "",
    ]
    if gene_ids:
        for key, value in sorted(gene_ids.items()):
            lines.append(f"- `{key}`: `{value}`")
    else:
        lines.append("- _(none)_")
    lines.extend(
        [
            "",
            "| Source | Priority | Status | Seconds | Endpoint | HTTP | Error |",
            "| --- | --- | --- | ---: | --- | ---: | --- |",
        ]
    )
    for row in rows:
        err = row.get("error_message") or "—"
        if len(err) > 80:
            err = err[:77] + "…"
        lines.append(
            f"| {row['source_name']} | {row['priority']} | `{row['status']}` | "
            f"{row['elapsed_seconds']:.2f} | `{row.get('endpoint_name') or '—'}` | "
            f"{row.get('status_code') if row.get('status_code') is not None else '—'} | "
            f"{err} |"
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    settings = get_settings()
    settings.ensure_dirs()
    gene = (args.gene or DEFAULT_GENE).strip()
    out_dir = Path(args.output_dir) if args.output_dir else settings.output_path
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = _select_sources(
        names=_split_csv(args.sources),
        priority=args.priority,
        identity_only=args.identity_only,
    )
    if not selected:
        LOGGER.error("No sources selected to smoke-test.")
        return 1

    LOGGER.info("Starting source smoke tests for %s (%d sources)", gene, len(selected))
    LOGGER.info("database=%s", _db_kind(settings))
    LOGGER.info("outputs=%s", out_dir)

    gene_ids: dict[str, Any] = {}
    bootstrap_rows: dict[str, dict[str, Any]] = {}
    # Bootstrap identity so ID-dependent clients (Reactome, Open Targets, …) can run.
    # Results are reused for selected identity sources to avoid a second live call.
    if not args.skip_identity_bootstrap:
        for name in IDENTITY_SOURCES:
            src = get_source(name)
            if src is None or src.name not in CLIENT_DISPATCH:
                continue
            status, result, elapsed, missing = _call_one(
                src, gene_symbol=gene, gene_ids=gene_ids, settings=settings
            )
            LOGGER.info(
                "identity bootstrap %s -> %s (%.2fs)", name, status, elapsed
            )
            if result is not None:
                gene_ids = extract_gene_ids_from_tool_result(result, gene_ids)
            bootstrap_rows[src.name] = _row_from_call(
                src,
                status=status,
                result=result,
                elapsed=elapsed,
                missing=missing,
            )

    rows: list[dict[str, Any]] = []
    for src in selected:
        if src.name in bootstrap_rows:
            row = bootstrap_rows[src.name]
            rows.append(row)
            LOGGER.info(
                "%-20s %-12s %.2fs %s (reused bootstrap)",
                src.name,
                row["status"],
                row["elapsed_seconds"],
                row.get("error_type") or "",
            )
            continue
        status, result, elapsed, missing = _call_one(
            src, gene_symbol=gene, gene_ids=gene_ids, settings=settings
        )
        if result is not None and result.success:
            gene_ids = extract_gene_ids_from_tool_result(result, gene_ids)
        row = _row_from_call(
            src,
            status=status,
            result=result,
            elapsed=elapsed,
            missing=missing,
        )
        rows.append(row)
        LOGGER.info(
            "%-20s %-12s %.2fs %s",
            src.name,
            status,
            elapsed,
            row.get("error_type") or "",
        )

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    stem = f"{gene.lower()}_source_smoke_{stamp}"
    md_path = out_dir / f"{stem}.md"
    json_path = out_dir / f"{stem}.json"
    md_path.write_text(
        _render_markdown(gene_symbol=gene, rows=rows, gene_ids=gene_ids),
        encoding="utf-8",
    )
    payload = {
        "gene_symbol": gene,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "gene_ids": gene_ids,
        "sources": rows,
        "summary": {
            "total": len(rows),
            "by_status": {
                status: sum(1 for r in rows if r["status"] == status)
                for status in sorted({r["status"] for r in rows})
            },
        },
    }
    json_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    summary = payload["summary"]["by_status"]
    print()
    print("=" * 72)
    print(f"Source smoke tests — {gene}")
    print("=" * 72)
    print(f"total:    {len(rows)}")
    print(f"by status:{summary}")
    if gene_ids:
        print("gene_ids:", ", ".join(f"{k}={v}" for k, v in sorted(gene_ids.items())))
    print(f"markdown: {md_path}")
    print(f"json:     {json_path}")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
