#!/usr/bin/env python3
"""One-time official AlphaFold viewer capture utility.

Use when AFDB blocks headless Chromium. Persists a managed derived_capture that
Section 1d can reuse. Never renders coordinates locally.

Example::

    PYTHONPATH=src .venv/bin/python scripts/capture_alphafold_viewer.py \\
      --model-entity-id AF-Q12772-F1 \\
      --accession Q12772 \\
      --headed
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import tempfile
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import get_settings  # noqa: E402
from gene_dossier.section_1d import (  # noqa: E402
    CAPTURE_FRESH,
    _capture_human_viewer,
    _write_capture_cache_index,
)

LOGGER = logging.getLogger("capture_alphafold_viewer")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Capture the official AlphaFold entry viewer for reuse."
    )
    parser.add_argument("--model-entity-id", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--taxon-id", type=int, default=9606)
    parser.add_argument("--model-version", type=int, default=None)
    parser.add_argument("--dossier-run-id", default="manual-alphafold-capture")
    parser.add_argument("--gene-symbol", default="MANUAL")
    parser.add_argument(
        "--headed",
        action="store_true",
        help="Launch a headed browser (recommended when AFDB blocks headless).",
    )
    parser.add_argument(
        "--channel",
        default=None,
        help='Optional Playwright browser channel, e.g. "chrome".',
    )
    parser.add_argument(
        "--user-data-dir",
        type=Path,
        default=None,
        help="Optional persistent Chromium profile directory.",
    )
    parser.add_argument(
        "--persist-profile",
        action="store_true",
        help="Use a temporary persistent profile under the system temp dir.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    settings = get_settings()
    user_data_dir = args.user_data_dir
    temp_profile: tempfile.TemporaryDirectory[str] | None = None
    if args.persist_profile and user_data_dir is None:
        temp_profile = tempfile.TemporaryDirectory(prefix="afdb-capture-profile-")
        user_data_dir = Path(temp_profile.name)

    try:
        api, meta, rec, audit = _capture_human_viewer(
            dossier_run_id=args.dossier_run_id,
            gene_symbol=args.gene_symbol,
            model_id=args.model_entity_id,
            accession=args.accession,
            taxon_id=args.taxon_id,
            model_version=args.model_version,
            parent_raw_artifact_ids=[],
            settings=settings,
            persist_db=False,
            headed=args.headed,
            channel=args.channel,
            user_data_dir=user_data_dir,
        )
    finally:
        if temp_profile is not None:
            temp_profile.cleanup()

    payload = {
        "api_run_id": api.id,
        "success": api.success,
        "error": api.error_message,
        "audit": audit,
        "meta": meta,
        "evidence_value": rec.value if rec is not None else None,
        "capture_mode": audit.get("capture_mode") or CAPTURE_FRESH,
    }
    print(json.dumps(payload, indent=2, sort_keys=True, default=str))
    if not api.success or rec is None:
        return 1

    # Ensure the reusable cache index points at the persisted artifact.
    value = rec.value if isinstance(rec.value, dict) else {}
    if value.get("model_entity_id"):
        _write_capture_cache_index(
            settings=settings,
            model_id=str(value["model_entity_id"]),
            payload=value,
        )
    LOGGER.info(
        "Captured %s → %s (sha256=%s)",
        args.model_entity_id,
        value.get("relative_path"),
        value.get("sha256"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
