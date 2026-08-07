#!/usr/bin/env python3
"""Smoke-render Section 5a for one gene via the section bundle."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gene_dossier.config import get_settings
from gene_dossier.section_5a import Section5aConfig
from gene_dossier.section_bundle import run_section_bundle


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--string-required-score", type=int, default=400)
    parser.add_argument("--no-string-network-figure", action="store_true")
    parser.add_argument("--promote-section-5a-accepted", action="store_true")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args(argv)
    settings = get_settings()
    result = run_section_bundle(
        args.gene,
        section_keys=["5a"],
        output_dir=args.output_dir,
        settings=settings,
        write_pdf=not args.no_pdf,
        promote_section_5a_accepted=args.promote_section_5a_accepted,
        section_5a_config=Section5aConfig(
            required_score=args.string_required_score,
            attempt_network_figure=not args.no_string_network_figure,
        ),
    )
    print(f"status={result.status}")
    print(f"output_dir={result.output_dir}")
    for key, path in sorted(result.output_paths.items()):
        print(f"  {key}: {path}")
    return 0 if result.status == "completed" and not result.errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
