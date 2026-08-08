#!/usr/bin/env python3
"""Section-scoped dossier generation for Sections 1a–5a.

Example::

    PYTHONPATH=src .venv/bin/python scripts/run_section_bundle.py \\
      --gene SREBF2 \\
      --output-dir data/outputs/section_validation/SREBF2 -v

Omit ``--sections`` to run the default bundle (1a through 4a). Pass
``--sections`` explicitly to override (e.g. ``--sections 5a``).
Section 5a (STRING) is supported/opt-in and not in the default bundle;
structured network scope is fixed at add_nodes=30.
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
from gene_dossier.section_1e import (  # noqa: E402
    SUPPORTED_SECTION_1E_SCOPES,
    Section1eConfig,
)
from gene_dossier.section_2a import Section2aConfig  # noqa: E402
from gene_dossier.section_2b import Section2bConfig  # noqa: E402
from gene_dossier.section_2c import Section2cConfig  # noqa: E402
from gene_dossier.section_3a import Section3aConfig  # noqa: E402
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
            "Generate a standalone section bundle for deterministic sections "
            "1a–4a by default, with opt-in 5a (STRING PPI). "
            "Default sections: 1a 1b 1c 1d 1e 2a 2b 2c 3a 4a."
        )
    )
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g. SREBF2)")
    parser.add_argument(
        "--sections",
        nargs="+",
        default=list(DEFAULT_SECTION_BUNDLE_KEYS),
        help=(
            "Section keys to include (1a–4a default; opt-in 5a). "
            "Default: 1a 1b 1c 1d 1e 2a 2b 2c 3a 4a. Explicit --sections "
            "overrides the default."
        ),
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
        "--ortholog-scope-tax-id",
        type=int,
        default=7776,
        help=(
            "NCBI taxonomy scope for Section 1e orthologs "
            f"(supported: {', '.join(str(k) for k in sorted(SUPPORTED_SECTION_1E_SCOPES))}; "
            "default: 7776 jawed vertebrates)"
        ),
    )
    parser.add_argument(
        "--max-visible-rows",
        type=int,
        default=20,
        help=(
            "Max NCBI Orthologs table rows to crop for Section 1e capture "
            "(default: 20)"
        ),
    )
    parser.add_argument(
        "--allen-probe-id",
        type=int,
        default=None,
        help=(
            "Optional explicit Allen HumanMA Agilent probe ID for Section 2b "
            "(ignored unless 2b is selected)"
        ),
    )
    parser.add_argument(
        "--no-allen-celltype-figures",
        action="store_true",
        help=(
            "Skip the optional Allen Cell Types Explorer browser captures for "
            "Section 2c (ignored unless 2c is selected)"
        ),
    )
    parser.add_argument(
        "--refresh-section-2c-sources",
        action="store_true",
        help=(
            "Re-download Section 2c dataset-level sources instead of reusing "
            "accepted artifacts (ignored unless 2c is selected)"
        ),
    )
    parser.add_argument(
        "--geo-max-candidates",
        type=int,
        default=500,
        help=(
            "Max GEO Profiles candidate UIDs to enrich for Section 3a "
            "(default: 500; ignored unless 3a is selected)"
        ),
    )
    parser.add_argument(
        "--geo-max-selected",
        type=int,
        default=6,
        help=(
            "Max polished GEO Profiles to select for Section 3a "
            "(default: 6; ignored unless 3a is selected)"
        ),
    )
    parser.add_argument(
        "--no-geo-profile-figures",
        action="store_true",
        help=(
            "Skip GEO Profile chart acquisition for Section 3a "
            "(ignored unless 3a is selected)"
        ),
    )
    parser.add_argument(
        "--refresh-section-3a",
        action="store_true",
        help=(
            "Force a fresh Section 3a attempt directory (ignored unless 3a is selected)"
        ),
    )
    parser.add_argument(
        "--acceptance-profile",
        default=None,
        choices=["section_1c_reference_genes", "section_1d_reference_genes"],
        help=(
            "Optional acceptance-validation profile. When it fails, outputs are "
            "preserved but the CLI exits nonzero."
        ),
    )
    parser.add_argument(
        "--promote-section-2c-accepted",
        action="store_true",
        help=(
            "Replace an existing accepted visual-complete Section 2c gene pointer only "
            "when the newly rendered report also passes the complete scientific, visual, "
            "figure-role, and PDF acceptance checks."
        ),
    )
    parser.add_argument(
        "--promote-section-3a-visual-accepted",
        action="store_true",
        help=(
            "Optionally pin/replace a visual-complete Section 3a accepted pointer after "
            "scientific-complete acceptance when charts and PDF checks pass."
        ),
    )
    parser.add_argument(
        "--harmonizome-max-curated-display",
        type=int,
        default=14,
        help=(
            "Max curated Harmonizome associations to display for Section 4a "
            "(default: 14; ignored unless 4a is selected)"
        ),
    )
    parser.add_argument(
        "--harmonizome-max-predicted-display",
        type=int,
        default=25,
        help=(
            "Max predicted Harmonizome associations to display for Section 4a "
            "(default: 25; ignored unless 4a is selected)"
        ),
    )
    parser.add_argument(
        "--promote-section-4a-accepted",
        action="store_true",
        help=(
            "Replace an existing successful Section 4a accepted pointer when the new "
            "attempt is also complete (default: keep prior successful pointer)."
        ),
    )
    parser.add_argument(
        "--string-required-score",
        type=int,
        default=400,
        help=(
            "STRING required_score (0-1000) for Section 5a structured network, "
            "official PNG, and get_link (default: 400; ignored unless 5a is selected). "
            "Structured add_nodes is fixed at 30."
        ),
    )
    parser.add_argument(
        "--no-string-network-figure",
        action="store_true",
        help=(
            "Skip the official STRING high-resolution network PNG for Section 5a "
            "(ignored unless 5a is selected)"
        ),
    )
    parser.add_argument(
        "--promote-section-5a-accepted",
        action="store_true",
        help=(
            "Replace an existing successful Section 5a accepted pointer when the new "
            "attempt is also complete (default: keep prior successful pointer)."
        ),
    )
    parser.add_argument(
        "--no-biogrid-network-figure",
        action="store_true",
        help=(
            "Skip the official BioGRID Network Viewer cy.png capture for Section 5b "
            "(ignored unless 5b is selected)"
        ),
    )
    parser.add_argument(
        "--promote-section-5b-accepted",
        action="store_true",
        help=(
            "Replace an existing successful Section 5b accepted pointer when the new "
            "attempt is also complete (default: keep prior successful pointer)."
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

    section_1e_config = None
    if "1e" in keys:
        try:
            section_1e_config = Section1eConfig(
                ortholog_scope_tax_id=args.ortholog_scope_tax_id,
                max_visible_rows=args.max_visible_rows,
            )
        except ValueError as exc:
            LOGGER.error("%s", exc)
            return 2

    section_2a_config = Section2aConfig() if "2a" in keys else None

    section_2b_config = None
    if "2b" in keys:
        try:
            section_2b_config = Section2bConfig(allen_probe_id=args.allen_probe_id)
        except ValueError as exc:
            LOGGER.error("%s", exc)
            return 2

    section_2c_config = None
    if "2c" in keys:
        try:
            section_2c_config = Section2cConfig(
                attempt_allen_figures=not args.no_allen_celltype_figures,
                force_refresh=args.refresh_section_2c_sources,
            )
        except ValueError as exc:
            LOGGER.error("%s", exc)
            return 2

    section_3a_config = None
    if "3a" in keys:
        try:
            section_3a_config = Section3aConfig(
                force_refresh=args.refresh_section_3a,
                max_discovery_profiles=args.geo_max_candidates,
                max_selected_profiles=args.geo_max_selected,
                attempt_figures=not args.no_geo_profile_figures,
            )
        except ValueError as exc:
            LOGGER.error("%s", exc)
            return 2

    section_4a_config = None
    if "4a" in keys:
        try:
            from gene_dossier.section_4a import Section4aConfig

            section_4a_config = Section4aConfig(
                max_displayed_curated_associations=args.harmonizome_max_curated_display,
                max_displayed_predicted_associations=args.harmonizome_max_predicted_display,
            )
        except ValueError as exc:
            LOGGER.error("%s", exc)
            return 2

    section_5a_config = None
    if "5a" in keys:
        try:
            from gene_dossier.section_5a import Section5aConfig

            section_5a_config = Section5aConfig(
                required_score=args.string_required_score,
                attempt_network_figure=not args.no_string_network_figure,
            )
        except ValueError as exc:
            LOGGER.error("%s", exc)
            return 2

    section_5b_config = None
    if "5b" in keys:
        from gene_dossier.section_5b import Section5bConfig

        section_5b_config = Section5bConfig(
            attempt_network_figure=not args.no_biogrid_network_figure,
        )

    settings = get_settings()
    result = run_section_bundle(
        args.gene,
        section_keys=keys,
        output_dir=args.output_dir,
        settings=settings,
        write_pdf=not args.no_pdf,
        dpi=args.dpi,
        acceptance_profile=args.acceptance_profile,
        promote_section_2c_accepted=args.promote_section_2c_accepted,
        promote_section_3a_visual_accepted=args.promote_section_3a_visual_accepted,
        promote_section_4a_accepted=args.promote_section_4a_accepted,
        promote_section_5a_accepted=args.promote_section_5a_accepted,
        promote_section_5b_accepted=args.promote_section_5b_accepted,
        section_1e_config=section_1e_config,
        section_2a_config=section_2a_config,
        section_2b_config=section_2b_config,
        section_2c_config=section_2c_config,
        section_3a_config=section_3a_config,
        section_4a_config=section_4a_config,
        section_5a_config=section_5a_config,
        section_5b_config=section_5b_config,
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
