#!/usr/bin/env python3
"""Smoke-render Section 4a for one gene via the section bundle."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gene_dossier.config import get_settings
from gene_dossier.section_4a import Section4aConfig
from gene_dossier.section_bundle import run_section_bundle


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--harmonizome-max-curated-display", type=int, default=14)
    parser.add_argument("--harmonizome-max-predicted-display", type=int, default=25)
    parser.add_argument("--promote-section-4a-accepted", action="store_true")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args(argv)
    settings = get_settings()
    result = run_section_bundle(
        args.gene,
        section_keys=["4a"],
        output_dir=args.output_dir,
        settings=settings,
        write_pdf=not args.no_pdf,
        promote_section_4a_accepted=args.promote_section_4a_accepted,
        section_4a_config=Section4aConfig(
            max_displayed_curated_associations=args.harmonizome_max_curated_display,
            max_displayed_predicted_associations=args.harmonizome_max_predicted_display,
        ),
    )
    print(f"status={result.status}")
    print(f"output_dir={result.output_dir}")
    for key, path in sorted(result.output_paths.items()):
        print(f"  {key}: {path}")
    return 0 if result.status == "completed" and not result.errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
