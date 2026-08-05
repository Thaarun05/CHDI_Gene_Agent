#!/usr/bin/env python3
"""Offline smoke: run Section 2c then render its focused bundle HTML.

Uses the already-accepted dataset-level sources under
``data/outputs/section_2c/accepted/sources`` and skips the optional Allen
Cell Types Explorer browser captures unless ``--figures`` is passed.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import get_settings  # noqa: E402
from gene_dossier.section_2c import (  # noqa: E402
    Section2cConfig,
    node_generate_section_2c_derived_artifacts,
)
from gene_dossier.section_bundle import (  # noqa: E402
    build_section_bundle_document,
    render_section_bundle_html,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--mouse-symbol", required=True)
    parser.add_argument("--figures", action="store_true")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)

    tmp = Path(tempfile.mkdtemp())
    # Figure resolution reads the process-wide artifact root, so managed
    # figure bytes must land in the configured raw-data directory.
    settings = get_settings()
    run_id = f"smoke-2c-{args.gene.lower()}"
    state = {
        "run_type": "section_bundle",
        "selected_section_keys": ["2c"],
        "dossier_run_id": run_id,
        "gene_symbol": args.gene,
        "gene_ids": {"mouse_symbol": args.mouse_symbol},
        "evidence_records": [],
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "coverage": [],
    }
    state = node_generate_section_2c_derived_artifacts(
        state,
        settings=settings,
        persist_db=False,
        config=Section2cConfig(
            output_root="data/outputs",
            attempt_allen_figures=args.figures,
        ),
    )
    status = state["section_2c_status"]
    document, presentation, audit = build_section_bundle_document(
        dossier_run_id=run_id,
        gene_symbol=args.gene,
        section_keys=["2c"],
        evidence_records=list(state.get("evidence_records") or []),
        api_runs=list(state.get("api_runs") or []),
        raw_artifacts=list(state.get("raw_artifacts") or []),
        section_status_by_key={"2c": status},
    )
    html = render_section_bundle_html(document)

    out_dir = args.out or (tmp / "render")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "section_2.html").write_text(html, encoding="utf-8")
    (out_dir / "presentation.json").write_text(
        json.dumps(presentation, indent=2), encoding="utf-8"
    )
    (out_dir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(json.dumps(status["rendering_status"], indent=2))
    print(f"pages={html.count('class=\"report-page section-bundle-body')}")
    print(f"page_breaks={html.count('page-break')}")
    print(f"figures={html.count('<img')}")
    print(f"blocks={len(presentation['major_sections'][0]['subsections'][0]['blocks'])}")
    print(f"out={out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
