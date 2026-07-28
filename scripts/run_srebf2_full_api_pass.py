#!/usr/bin/env python3
"""Run a full Gene Dossier API pass for SREBF2 (or another gene).

Primary deliverable entrypoint from ``IMPLEMENTATION_PLAN.md``::

    python scripts/run_srebf2_full_api_pass.py

What it does:
  1. Creates a dossier run
  2. Resolves gene identity (NCBI Gene → Ensembl → UniProt)
  3. Calls every registered source client (soft-fail per source)
  4. Saves raw artifacts + optional DB provenance rows
  5. Normalizes evidence records
  6. Synthesizes sections (deterministic unless ``--allow-llm``)
  7. Verifies claims
  8. Writes coverage, debug markdown, and polished Rancho report

All API / LLM keys are optional. Missing keys degrade gracefully
(``requires_key`` / deterministic synthesis) rather than aborting the run.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Allow ``python scripts/run_srebf2_full_api_pass.py`` without an editable install.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import get_settings  # noqa: E402
from gene_dossier.report_presentation import build_section_presentation  # noqa: E402
from gene_dossier.workflow import run_gene_dossier_full_api_pass  # noqa: E402

DEFAULT_GENE = "SREBF2"
LOGGER = logging.getLogger("run_srebf2_full_api_pass")


def _llm_summary(settings) -> tuple[str | None, str | None]:
    provider = settings.default_llm_provider or None
    model = (
        settings.default_llm_model
        or settings.nvidia_nim_model
        or None
    )
    return provider, model


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the Gene Dossier full API pass (default gene: SREBF2). "
            "Writes coverage, debug markdown, and Rancho HTML/JSON/(PDF)."
        )
    )
    parser.add_argument(
        "--gene",
        default=DEFAULT_GENE,
        help=f"Gene symbol to dossier (default: {DEFAULT_GENE})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: settings.OUTPUT_DIR / data/outputs)",
    )
    parser.add_argument(
        "--run-id",
        default=None,
        help="Optional fixed dossier_run_id (otherwise auto-generated)",
    )
    parser.add_argument(
        "--sources",
        default=None,
        help=(
            "Optional comma-separated source name filter for non-identity clients "
            '(e.g. "GTEx,STRING,Reactome"). Identity sources always run first.'
        ),
    )
    parser.add_argument(
        "--allow-llm",
        action="store_true",
        help="Allow LangChain LLM synthesis when API keys are configured "
        "(default: deterministic-only).",
    )
    parser.add_argument(
        "--no-rancho",
        action="store_true",
        help="Skip polished Rancho HTML/PDF report (debug markdown still written).",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Write Rancho HTML/JSON but skip PDF export.",
    )
    parser.add_argument(
        "--no-db",
        action="store_true",
        help="Do not persist dossier/api/evidence/coverage rows to the database.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable DEBUG logging.",
    )
    return parser.parse_args(argv)


def _split_sources(raw: str | None) -> list[str] | None:
    if not raw:
        return None
    items = [part.strip() for part in raw.split(",")]
    return [item for item in items if item]


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    settings = get_settings()
    settings.ensure_dirs()
    gene = (args.gene or DEFAULT_GENE).strip()
    sources = _split_sources(args.sources)
    actual_output_dir = args.output_dir or settings.output_path
    db_url = settings.database_url or ""
    db_kind = (
        "postgres"
        if db_url.startswith(("postgresql://", "postgresql+psycopg://"))
        else "sqlite"
    )

    LOGGER.info("Starting full API pass for %s", gene)
    LOGGER.info("database=%s", db_kind)
    LOGGER.info("raw=%s  outputs=%s", settings.raw_data_path, actual_output_dir)
    if settings.has_llm() and args.allow_llm:
        LOGGER.info("LLM keys present; LLM synthesis enabled (--allow-llm).")
    elif settings.has_llm():
        LOGGER.info("LLM keys present but unused (pass --allow-llm to enable).")
    else:
        LOGGER.info("No LLM keys; using deterministic synthesis.")

    result = run_gene_dossier_full_api_pass(
        gene,
        settings=settings,
        output_dir=args.output_dir,
        dossier_run_id=args.run_id,
        sources=sources,
        call_network=True,
        force_deterministic=not args.allow_llm,
        write_rancho=not args.no_rancho,
        write_pdf=not args.no_pdf and not args.no_rancho,
        persist_db=not args.no_db,
    )

    print()
    print("=" * 72)
    print(f"Gene Dossier full API pass — {result.gene_symbol}")
    print("=" * 72)
    print(f"status:          {result.status}")
    print(f"dossier_run_id:  {result.dossier_run_id}")
    print(f"evidence:        {len(result.evidence_records)}")
    print(f"claims:          {len(result.claims)}")
    print(f"verification:    {len(result.verification_results)}")
    print(f"coverage rows:   {len(result.coverage)}")
    print(f"synthesis mode:  {result.synthesis_mode or 'n/a'}")
    provider, model = _llm_summary(settings)
    if result.synthesis_mode == "llm":
        print(f"llm provider:    {provider or 'configured'}")
        print(f"llm model:       {model or 'default'}")

    sec1a = build_section_presentation(
        section_key="1a",
        gene_symbol=result.gene_symbol,
        evidence_records=result.evidence_records,
    )
    sec1b = build_section_presentation(
        section_key="1b",
        gene_symbol=result.gene_symbol,
        evidence_records=result.evidence_records,
    )
    ucsc_facts = {
        rec.fact_type
        for rec in result.evidence_records
        if rec.source_name == "UCSC" and rec.fact_type.startswith("ucsc_")
    }
    figure_recs = [
        rec
        for rec in result.evidence_records
        if rec.fact_type == "ucsc_conservation_figure"
    ]
    figure_status = "missing"
    if figure_recs:
        value = figure_recs[-1].value or {}
        figure_status = str(value.get("retrieval_method") or "present")

    print(f"section 1a:      {len(sec1a.blocks)} presentation blocks")
    print(f"section 1b:      {len(sec1b.blocks)} presentation blocks")
    print(
        "ucsc facts:      "
        + ", ".join(
            t
            for t in (
                "ucsc_gene_locus",
                "ucsc_transcript_inventory",
                "ucsc_canonical_transcript",
                "ucsc_conservation_figure",
            )
            if t in ucsc_facts
        )
    )
    print(f"ucsc figure:     {figure_status}")

    if result.gene_ids:
        print()
        print("Resolved identifiers:")
        for key in (
            "official_symbol",
            "entrez_gene_id",
            "ensembl_id",
            "uniprot_accession",
            "chromosome",
            "gtex_gencode_id",
            "mouse_entrez_id",
            "mgi_id",
        ):
            if key in result.gene_ids:
                print(f"  {key}: {result.gene_ids[key]}")

    if result.output_paths:
        print()
        print("Outputs:")
        for key in sorted(result.output_paths):
            print(f"  {key}: {result.output_paths[key]}")

    if result.errors:
        print()
        print(f"Source-level errors / notes ({len(result.errors)}):")
        for err in result.errors[:25]:
            print(f"  - {err}")
        if len(result.errors) > 25:
            print(f"  … {len(result.errors) - 25} more")

    if result.synthesis_notes:
        print()
        print("Synthesis notes:")
        for note in result.synthesis_notes:
            print(f"  - {note}")

    print()
    if result.status == "completed":
        LOGGER.info("Pass completed successfully.")
        return 0

    LOGGER.error("Pass finished with status=%s", result.status)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
