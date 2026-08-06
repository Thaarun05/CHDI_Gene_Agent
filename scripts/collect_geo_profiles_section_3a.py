#!/usr/bin/env python3
"""Collect Section 3a GEO Profiles candidates for a gene (no polished report)."""

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
from gene_dossier.section_3a_sources import paths_for, write_json_atomic  # noqa: E402
from gene_dossier.tools import geo_profiles as gp  # noqa: E402

LOGGER = logging.getLogger("collect_geo_profiles_section_3a")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Collect and rank GEO Profiles for Section 3a screening."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Output root (default: settings.output_path/section_3a)",
    )
    parser.add_argument("--max-candidates", type=int, default=500)
    parser.add_argument("--max-selected", type=int, default=6)
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Skip profile HTML/chart acquisition",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Create a new attempt directory (always true for this CLI)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    _ = args.refresh

    settings = get_settings()
    gene = args.gene.strip()
    collected = gp.collect_section_3a_profiles(
        gene,
        max_discovery_profiles=args.max_candidates,
        max_selected_profiles=args.max_selected,
        attempt_figures=not args.no_figures,
        settings=settings,
    )
    paths = paths_for(args.output_root or settings.output_path)
    attempt = paths.new_gene_attempt(gene)
    # Strip image bytes before JSON persistence.
    candidates = []
    for cand in list(collected.get("candidates") or []):
        row = {k: v for k, v in cand.items() if k not in {"graph_image_bytes", "gds_metadata"}}
        row["gds_metadata"] = {
            k: v
            for k, v in dict(cand.get("gds_metadata") or {}).items()
            if k != "samples"
        }
        candidates.append(row)
    selected = []
    for cand in list(collected.get("selected_profiles") or []):
        selected.append(
            {k: v for k, v in cand.items() if k not in {"graph_image_bytes"}}
        )
    write_json_atomic(
        attempt / "summary.json",
        {
            "gene_symbol": gene,
            "exact_profile_count": collected.get("exact_profile_count"),
            "neural_profile_count": collected.get("neural_profile_count"),
            "subset_effect_profile_count": collected.get("subset_effect_profile_count"),
            "candidate_union_count": collected.get("candidate_union_count"),
            "selected_profile_count": collected.get("selected_profile_count"),
            "scientific_status": collected.get("scientific_status"),
            "visual_status": collected.get("visual_status"),
            "max_chart_candidates": collected.get("max_chart_candidates"),
            "attempt_figures": collected.get("attempt_figures"),
        },
    )
    write_json_atomic(attempt / "candidates.json", {"candidates": candidates})
    write_json_atomic(attempt / "selected_profiles.json", {"selected_profiles": selected})

    print(f"gene={gene}")
    print(f"attempt_dir={attempt}")
    print(f"exact_profile_count={collected.get('exact_profile_count')}")
    print(f"neural_profile_count={collected.get('neural_profile_count')}")
    print(f"subset_effect_profile_count={collected.get('subset_effect_profile_count')}")
    print(f"candidate_union_count={collected.get('candidate_union_count')}")
    print(f"selected_profile_count={collected.get('selected_profile_count')}")
    print(f"scientific_status={collected.get('scientific_status')}")
    print(f"visual_status={collected.get('visual_status')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
