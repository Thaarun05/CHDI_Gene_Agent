#!/usr/bin/env python
"""Acquire and pin the CTD chemical–gene interactions bulk gzip for Section 6a.

Downloads::

    https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz

using the same download/persistence path as gene-run bootstrap so a successful
DB-backed acquire always produces one real ApiRun and one raw gzip RawArtifact,
then pins ``accepted/sources/ctd_chem_gene_ixns.json`` with those ORIGINAL ids.

``--no-db`` may download/validate for debugging and keep local attempt artifacts,
but never creates or replaces the accepted source pointer.

Usage::

    PYTHONPATH=src .venv/bin/python scripts/acquire_ctd_chem_gene_ixns.py \\
        --output-root data/outputs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_SRC = _REPO / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import get_settings  # noqa: E402
from gene_dossier.db import init_db, save_dossier_run, session_scope  # noqa: E402
from gene_dossier.models import DossierRun  # noqa: E402
from gene_dossier.section_6a import (  # noqa: E402
    PARSER_VERSION,
    download_persist_and_pin_ctd_bulk,
)
from gene_dossier.section_6a_sources import (  # noqa: E402
    OFFICIAL_URL,
    SOURCE_KEY,
    load_accepted_source,
    paths_for,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Acquire CTD chem-gene bulk for Section 6a")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Base output root (default: settings.output_path)",
    )
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="Re-download even when an accepted pointer exists",
    )
    parser.add_argument(
        "--no-db",
        action="store_true",
        help=(
            "Debug download only: skip DB persistence and do NOT create or replace "
            "accepted/sources/ctd_chem_gene_ixns.json (accepted pointers require "
            "resolvable ApiRun + RawArtifact provenance)"
        ),
    )
    args = parser.parse_args(argv)
    cfg = get_settings()
    paths = paths_for(args.output_root or cfg.output_path)
    persist_db = not args.no_db
    # Accepted CTD pointers must cite ORIGINAL resolvable ApiRun/RawArtifact IDs.
    promote = persist_db

    if not args.force_refresh:
        cached = load_accepted_source(paths, source_key=SOURCE_KEY, official_url=OFFICIAL_URL)
        if cached:
            print(
                f"[ctd] cache hit sha256={cached.get('sha256')} "
                f"bytes={cached.get('byte_size')} attempt={cached.get('source_attempt_id')} "
                f"api_run_id={cached.get('api_run_id')} "
                f"raw_artifact_id={cached.get('raw_artifact_id')}"
            )
            return 0

    run = DossierRun(
        gene_symbol="CTD_BULK",
        run_type="ctd_source_acquire",
        status="running",
        notes=f"section_6a CTD bulk acquire ({PARSER_VERSION})",
        config={"source_key": SOURCE_KEY, "official_url": OFFICIAL_URL},
    )
    if persist_db:
        init_db()
        with session_scope() as session:
            save_dossier_run(session, run)

    print(f"[ctd] downloading {OFFICIAL_URL}")
    payload = download_persist_and_pin_ctd_bulk(
        paths=paths,
        dossier_run_id=run.id,
        settings=cfg,
        gene_symbol="CTD_BULK",
        persist_db=persist_db,
        promote=promote,
        origin="acquire_script",
    )
    if not payload.ok:
        print(f"[ctd] FAILED {payload.error_type}: {payload.error_message}")
        if persist_db:
            run.status = "failed"
            with session_scope() as session:
                save_dossier_run(session, run)
        return 1

    if persist_db:
        run.status = "completed"
        with session_scope() as session:
            save_dossier_run(session, run)

    if promote:
        print(
            f"[ctd] accepted sha256={payload.sha256} bytes={payload.byte_size} "
            f"attempt={payload.source_attempt_id} "
            f"api_run_id={payload.api_run_id} raw_artifact_id={payload.raw_artifact_id}"
        )
    else:
        print(
            f"[ctd] downloaded (not pinned; --no-db) sha256={payload.sha256} "
            f"bytes={payload.byte_size} attempt={payload.source_attempt_id}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
