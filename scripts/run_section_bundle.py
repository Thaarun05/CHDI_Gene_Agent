#!/usr/bin/env python3
"""Section-scoped generation for Sections 1a / 1b / opt-in 1c.

Example::

    PYTHONPATH=src .venv/bin/python scripts/run_section_bundle.py \\
      --gene SREBF2 --sections 1a 1b \\
      --output-dir data/outputs/section_validation/SREBF2 -v
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import get_settings  # noqa: E402
from gene_dossier.section_bundle import (  # noqa: E402
    DEFAULT_SECTION_BUNDLE_KEYS,
    SectionBundleError,
    run_section_bundle,
    validate_section_keys,
)

LOGGER = logging.getLogger("run_section_bundle")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a standalone Section 1 bundle (1a Gene Aliases / "
            "1b UCSC conservation / opt-in 1c Known structure) without LLM "
            "synthesis or full-report rendering."
        )
    )
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g. SREBF2)")
    parser.add_argument(
        "--sections",
        nargs="+",
        default=list(DEFAULT_SECTION_BUNDLE_KEYS),
        help="Section keys to include (1a, 1b, and/or opt-in 1c). Default: 1a 1b",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Base output directory; files are written under "
            "<output-dir>/<dossier_run_id>/ (default: "
            "data/outputs/section_validation/<GENE>)"
        ),
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Skip PDF/PNG generation (HTML + JSON still written)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="PNG rasterization DPI (default: 150)",
    )
    parser.add_argument(
        "--acceptance-profile",
        default=None,
        choices=["section_1c_reference_genes"],
        help=(
            "Optional acceptance-validation profile. When it fails, outputs are "
            "preserved but the CLI exits nonzero."
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose logging",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    try:
        keys = validate_section_keys(args.sections)
    except SectionBundleError as exc:
        LOGGER.error("%s", exc)
        return 2

    settings = get_settings()
    result = run_section_bundle(
        args.gene,
        section_keys=keys,
        output_dir=args.output_dir,
        settings=settings,
        write_pdf=not args.no_pdf,
        dpi=args.dpi,
        acceptance_profile=args.acceptance_profile,
    )

    print(f"status={result.status}")
    print(f"gene={result.gene_symbol}")
    print(f"dossier_run_id={result.dossier_run_id}")
    print(f"sections={','.join(result.selected_section_keys)}")
    print(f"output_dir={result.output_dir}")
    for key in sorted(result.output_paths):
        print(f"  {key}: {result.output_paths[key]}")
    if result.errors:
        for err in result.errors:
            print(f"error: {err}")
        return 1
    return 0 if result.status == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
