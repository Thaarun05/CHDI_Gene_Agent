#!/usr/bin/env python3
"""Collect DropViz Section 2c acquisition artifacts for one mouse gene.

Validation saved-state URLs live **only** in this CLI / test config — never in
``gene_dossier.tools.dropviz``.

Example::

    PYTHONPATH=src .venv/bin/python scripts/collect_dropviz_section_2c.py \\
      --mouse-gene-symbol Srebf2 \\
      --state-url 'http://dropviz.org/?_state_id_=05190dfa61f331d8' \\
      --state-url 'http://dropviz.org/?_state_id_=5c4cbc26b012914c' \\
      --output-dir data/outputs/section_2c_collection/Srebf2 -v
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.tools.dropviz import collect_dropviz_gene  # noqa: E402

LOGGER = logging.getLogger("collect_dropviz_section_2c")

# Validation-only defaults (CLI/config — not production client constants).
VALIDATION_STATE_URLS: dict[str, tuple[str, ...]] = {
    "Srebf2": (
        "http://dropviz.org/?_state_id_=05190dfa61f331d8",
        "http://dropviz.org/?_state_id_=5c4cbc26b012914c",
    ),
    "Cdh10": (
        "http://dropviz.org/?_state_id_=719bad29fe0f17fc",
        "http://dropviz.org/?_state_id_=3814b82b16caaf25",
    ),
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Collect DropViz Section 2c data for a mouse gene symbol."
    )
    parser.add_argument(
        "--mouse-gene-symbol",
        required=True,
        help="Mouse gene symbol (e.g. Srebf2, Cdh10).",
    )
    parser.add_argument(
        "--state-url",
        action="append",
        default=[],
        help="Saved-state DropViz URL (repeatable). "
        "If omitted and symbol has validation defaults, those URLs are used.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: data/outputs/section_2c_collection/<symbol>).",
    )
    parser.add_argument(
        "--no-default-states",
        action="store_true",
        help="Do not inject validation default state URLs when --state-url is omitted.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    symbol = args.mouse_gene_symbol.strip()
    urls = list(args.state_url or [])
    if not urls and not args.no_default_states:
        urls = list(VALIDATION_STATE_URLS.get(symbol, ()))
        if urls:
            LOGGER.info("Using %d validation default state URL(s) for %s", len(urls), symbol)

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = (
            _REPO_ROOT / "data" / "outputs" / "section_2c_collection" / symbol
        )

    LOGGER.info("Collecting DropViz for %s → %s", symbol, output_dir)
    result = collect_dropviz_gene(
        mouse_gene_symbol=symbol,
        output_dir=output_dir,
        saved_state_urls=urls,
    )

    data = result.data if isinstance(result.data, dict) else {}
    print(
        json.dumps(
            {
                "success": result.success,
                "error_type": result.error_type,
                "error_message": result.error_message,
                "status": data.get("status"),
                "payload_summary": {
                    "gene_symbol": (data.get("payload") or {}).get("gene_symbol"),
                    "acquisition_status": (data.get("payload") or {}).get(
                        "acquisition_status"
                    ),
                    "rank_extraction_status": (data.get("payload") or {}).get(
                        "rank_extraction_status"
                    ),
                    "regional_quantitative_status": (data.get("payload") or {}).get(
                        "regional_quantitative_status"
                    ),
                    "view_type": (data.get("payload") or {}).get("view_type"),
                    "summary_path": (data.get("payload") or {}).get("summary_path"),
                    "state_count": len((data.get("payload") or {}).get("states") or []),
                },
                "audit_keys": sorted((data.get("audit") or {}).keys()),
            },
            indent=2,
        )
    )
    return 0 if result.success else 1


if __name__ == "__main__":
    raise SystemExit(main())
