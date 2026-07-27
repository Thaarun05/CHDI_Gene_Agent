#!/usr/bin/env python3
"""Deduplicate stale species-identity EvidenceRecords for one dossier run.

Example::

    PYTHONPATH=src .venv/bin/python scripts/dedupe_species_identity_evidence.py \\
      --gene SREBF2 \\
      --dossier-run-id 9923bf6d326246a4a9abb6d56053c898
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.db import (  # noqa: E402
    delete_evidence_record,
    get_dossier_run,
    list_evidence_for_run,
    session_scope,
)
from gene_dossier.identity_hygiene import (  # noqa: E402
    assert_dossier_gene_matches,
    canonical_primary_identifier,
    dedupe_species_identity_records,
    is_species_identity_record,
    resolve_taxon_id,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Remove stale duplicate species-identity EvidenceRecords."
    )
    parser.add_argument("--gene", required=True, help="Dossier query gene symbol")
    parser.add_argument("--dossier-run-id", required=True)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report removals without deleting rows",
    )
    args = parser.parse_args(argv)

    query = args.gene.strip()
    run_id = args.dossier_run_id.strip()

    with session_scope() as session:
        run = get_dossier_run(session, run_id)
        if run is None:
            raise SystemExit(f"No dossier run found for dossier_run_id={run_id!r}.")
        try:
            assert_dossier_gene_matches(run.gene_symbol, query)
        except ValueError as exc:
            raise SystemExit(f"{exc} Aborting.") from exc

        evidence = list_evidence_for_run(session, run_id)
        before_identity = [r for r in evidence if is_species_identity_record(r)]
        result = dedupe_species_identity_records(evidence, query_symbol=query)
        retained_identity = [
            r for r in result.retained if is_species_identity_record(r)
        ]

        print(f"Identity records before: {len(before_identity)}")
        print(f"Identity records retained: {len(retained_identity)}")
        print(f"Stale identity removed: {len(result.removed)}")
        for rec in result.removed:
            print(
                f"  REMOVE id={rec.id} source={rec.source_name} "
                f"fact={rec.fact_type} tax={resolve_taxon_id(rec)} "
                f"primary={canonical_primary_identifier(rec)} "
                f"gene_symbol={rec.gene_symbol!r} "
                f"source_id={rec.source_id}"
            )

        if args.dry_run:
            print("Dry run: no deletions performed.")
            return 0

        deleted = 0
        for rec in result.removed:
            if delete_evidence_record(session, rec.id):
                deleted += 1
        print(f"Deleted evidence rows: {deleted}")
        print("RawArtifacts and ApiRuns were not modified.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
