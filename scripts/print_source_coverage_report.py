#!/usr/bin/env python3
"""Print (and optionally write) a source coverage report.

The platform must never silently omit a source. This script prints coverage for:

1. A persisted dossier run (``--run-id``) loaded from the provenance DB
2. An existing coverage JSON file (``--from-json``)
3. The live registry baseline for the current environment (default)

Usage::

    python scripts/print_source_coverage_report.py
    python scripts/print_source_coverage_report.py --run-id <dossier_run_id>
    python scripts/print_source_coverage_report.py --from-json data/outputs/run_source_coverage.json
    python scripts/print_source_coverage_report.py --write --gene SREBF2

Never logs raw ``DATABASE_URL`` credentials.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import Settings, get_settings  # noqa: E402
from gene_dossier.coverage import (  # noqa: E402
    build_coverage_for_registry,
    coverage_to_jsonable,
    render_coverage_markdown,
    summarize_coverage,
    write_coverage_report,
)
from gene_dossier.db import (  # noqa: E402
    get_dossier_run,
    init_db,
    list_source_coverage,
    session_scope,
)
from gene_dossier.models import SourceCoverageResult, SourceStatus  # noqa: E402

LOGGER = logging.getLogger("print_source_coverage_report")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Print a Gene Dossier source coverage report from the registry baseline, "
            "a persisted dossier run, or an existing coverage JSON file."
        )
    )
    parser.add_argument(
        "--run-id",
        default=None,
        help="Load coverage rows for this dossier_run_id from the database",
    )
    parser.add_argument(
        "--from-json",
        type=Path,
        default=None,
        help="Load coverage from an existing *_source_coverage.json file",
    )
    parser.add_argument(
        "--gene",
        default=None,
        help="Optional gene symbol shown in the report header",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory used when --write is set (default: settings.OUTPUT_DIR)",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Also write markdown + JSON coverage files under the output directory",
    )
    parser.add_argument(
        "--json-stdout",
        action="store_true",
        help="Print JSON to stdout instead of markdown",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable DEBUG logging",
    )
    return parser.parse_args(argv)


def _db_kind(settings: Settings) -> str:
    url = settings.database_url or ""
    if url.startswith(("postgresql://", "postgresql+psycopg://")):
        return "postgres"
    if url.startswith("sqlite:"):
        return "sqlite"
    return "other"


def _results_from_json(path: Path) -> tuple[str, str | None, list[SourceCoverageResult]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    dossier_run_id = str(payload.get("dossier_run_id") or path.stem)
    gene_symbol = payload.get("gene_symbol")
    rows_raw = payload.get("sources") or payload.get("coverage") or []
    results: list[SourceCoverageResult] = []
    for row in rows_raw:
        if not isinstance(row, dict):
            continue
        status_raw = row.get("status") or SourceStatus.not_implemented.value
        try:
            status = SourceStatus(status_raw)
        except ValueError:
            status = SourceStatus.not_implemented
        results.append(
            SourceCoverageResult(
                dossier_run_id=str(row.get("dossier_run_id") or dossier_run_id),
                source_name=str(row.get("source_name") or "unknown"),
                status=status,
                raw_artifact_path=row.get("raw_artifact_path"),
                evidence_record_count=int(row.get("evidence_record_count") or 0),
                error_message=row.get("error_message"),
                report_sections_supported=list(row.get("report_sections_supported") or []),
                notes=row.get("notes"),
            )
        )
    return dossier_run_id, gene_symbol, results


def _results_from_db(
    dossier_run_id: str,
) -> tuple[str, str | None, list[SourceCoverageResult]]:
    init_db()
    with session_scope() as session:
        run = get_dossier_run(session, dossier_run_id)
        if run is None:
            raise SystemExit(f"Dossier run not found: {dossier_run_id}")
        rows = list_source_coverage(session, dossier_run_id)
        gene_symbol = run.gene_symbol
    if not rows:
        LOGGER.warning(
            "No coverage rows persisted for run %s; falling back to registry baseline",
            dossier_run_id,
        )
        rows = build_coverage_for_registry(dossier_run_id)
    return dossier_run_id, gene_symbol, rows


def _print_summary(results: list[SourceCoverageResult]) -> None:
    summary = summarize_coverage(results)
    print(f"total sources: {summary['total']}")
    print(f"by status:     {summary['by_status']}")


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.run_id and args.from_json:
        LOGGER.error("Use only one of --run-id or --from-json")
        return 2

    settings = get_settings()
    settings.ensure_dirs()
    LOGGER.info("database=%s", _db_kind(settings))

    gene_symbol = args.gene
    if args.from_json is not None:
        path = Path(args.from_json)
        if not path.is_file():
            LOGGER.error("Coverage JSON not found: %s", path)
            return 1
        dossier_run_id, file_gene, results = _results_from_json(path)
        gene_symbol = gene_symbol or file_gene
        LOGGER.info("Loaded coverage from %s (%d sources)", path, len(results))
    elif args.run_id:
        dossier_run_id, run_gene, results = _results_from_db(args.run_id)
        gene_symbol = gene_symbol or run_gene
        LOGGER.info(
            "Loaded coverage for run %s (%d sources)", dossier_run_id, len(results)
        )
    else:
        dossier_run_id = "registry-baseline"
        results = build_coverage_for_registry(dossier_run_id, settings=settings)
        LOGGER.info("Built registry baseline coverage (%d sources)", len(results))

    if args.json_stdout:
        payload: dict[str, Any] = {
            "dossier_run_id": dossier_run_id,
            "gene_symbol": gene_symbol,
            "summary": summarize_coverage(results),
            "sources": coverage_to_jsonable(results),
        }
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        print()
        print("=" * 72)
        title = f"Source coverage — {gene_symbol}" if gene_symbol else "Source coverage"
        print(title)
        print(f"dossier_run_id: {dossier_run_id}")
        print("=" * 72)
        _print_summary(results)
        print()
        print(
            render_coverage_markdown(
                results, dossier_run_id=dossier_run_id, gene_symbol=gene_symbol
            )
        )

    if args.write:
        out_dir = Path(args.output_dir) if args.output_dir else settings.output_path
        paths = write_coverage_report(
            results,
            dossier_run_id,
            gene_symbol=gene_symbol,
            output_dir=out_dir,
        )
        # Keep write notices off the JSON stdout stream.
        sink = sys.stderr if args.json_stdout else sys.stdout
        print(file=sink)
        print(f"wrote markdown: {paths['markdown']}", file=sink)
        print(f"wrote json:     {paths['json']}", file=sink)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
