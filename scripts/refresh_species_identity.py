#!/usr/bin/env python3
"""Refresh Human/Mouse/Rat identity evidence for an existing dossier run.

Fetches only NCBI Gene / Ensembl / UniProt species identity (no LLM, no other
sources), normalizes EvidenceRecords, persists them, and optionally regenerates
the Section 1a preview.

Example::

    PYTHONPATH=src .venv/bin/python scripts/refresh_species_identity.py \\
      --gene SREBF2 \\
      --dossier-run-id 9923bf6d326246a4a9abb6d56053c898 \\
      --render-preview
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import get_settings  # noqa: E402
from gene_dossier.db import (  # noqa: E402
    list_evidence_for_run,
    save_api_run,
    save_evidence_record,
    save_raw_artifact,
    session_scope,
)
from gene_dossier.models import ApiRun  # noqa: E402
from gene_dossier.raw_store import RawStore  # noqa: E402
from gene_dossier.species_identity import fetch_species_identity_results  # noqa: E402
from gene_dossier.workflow import (  # noqa: E402
    extract_gene_ids_from_tool_result,
    normalize_tool_result,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Refresh species identity evidence (NCBI/Ensembl/UniProt)."
    )
    parser.add_argument("--gene", required=True)
    parser.add_argument("--dossier-run-id", required=True)
    parser.add_argument(
        "--skip-human",
        action="store_true",
        help="Skip human taxon (9606); refresh mouse/rat only",
    )
    parser.add_argument(
        "--render-preview",
        action="store_true",
        help="Regenerate Section 1a HTML preview after refresh",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_REPO_ROOT / "data" / "outputs" / "section_previews",
    )
    args = parser.parse_args(argv)

    gene = args.gene.strip()
    run_id = args.dossier_run_id.strip()
    settings = get_settings()
    skip = {9606} if args.skip_human else set()

    print(f"Fetching species identity for {gene} (skip_taxons={sorted(skip) or 'none'})...")
    results = fetch_species_identity_results(
        gene, settings=settings, skip_taxons=skip
    )

    store = RawStore(base_dir=settings.raw_data_path)
    gene_ids: dict = {}
    new_evidence = []
    for result in results:
        gene_ids = extract_gene_ids_from_tool_result(result, gene_ids)
        status = "OK" if result.success else f"FAIL:{result.error_type}"
        print(f"  {result.source_name}/{result.endpoint_name}: {status}")
        if (
            result.source_name == "NCBI Gene"
            and isinstance(result.data, dict)
            and result.data.get("selection_warnings")
        ):
            print(f"    warnings={result.data.get('selection_warnings')}")

        api_run = ApiRun(
            dossier_run_id=run_id,
            gene_symbol=gene,
            source_name=result.source_name,
            endpoint_name=result.endpoint_name,
            request_url=result.request_url,
            request_params=dict(result.request_params or {}),
            status_code=result.status_code,
            success=result.success,
            error_type=result.error_type,
            error_message=result.error_message,
        )
        artifact = None
        if result.data is not None:
            artifact = store.save_json(
                run_id,
                result.source_name,
                result.data,
                api_run_id=api_run.id,
                original_url=result.request_url or None,
                filename_hint=result.endpoint_name,
            )
            api_run.raw_artifact_id = artifact.id
            result.raw_artifact_id = artifact.id

        with session_scope() as session:
            save_api_run(session, api_run)
            if artifact is not None:
                save_raw_artifact(session, artifact)

        if not result.success:
            continue
        batch = normalize_tool_result(
            result,
            dossier_run_id=run_id,
            api_run_id=api_run.id,
            raw_artifact_id=result.raw_artifact_id,
        )
        with session_scope() as session:
            existing_by_source = {
                e.source_id: e
                for e in list_evidence_for_run(session, run_id)
                if e.source_id
            }
            for rec in batch:
                prior = existing_by_source.get(rec.source_id)
                if prior is not None:
                    # Upsert richer identity fields onto the prior row id.
                    rec = rec.model_copy(update={"id": prior.id})
                save_evidence_record(session, rec)
                existing_by_source[rec.source_id] = rec
                new_evidence.append(rec)

    print(f"Persisted {len(new_evidence)} new evidence records.")
    print(f"gene_ids keys: {sorted(gene_ids)}")

    if args.render_preview:
        from gene_dossier.rancho_report import render_rancho_section_fragment
        from gene_dossier.report_presentation import build_section_presentation
        from gene_dossier.report_schema import build_report_document

        with session_scope() as session:
            evidence = list_evidence_for_run(session, run_id)
        presentation = build_section_presentation(
            section_key="1a", gene_symbol=gene, evidence_records=evidence
        )
        for d in presentation.diagnostics:
            print(f"  diagnostic[{d.severity}]: {d.field}: {d.reason}")
        doc = build_report_document(
            gene_symbol=gene,
            dossier_run_id=run_id,
            evidence_records=evidence,
        )
        html = render_rancho_section_fragment(
            document=doc,
            section_number=1,
            subsection_key="a",
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        out = args.output_dir / f"{gene}_1a_gene_aliases.html"
        out.write_text(html, encoding="utf-8")
        print(f"Wrote preview: {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
