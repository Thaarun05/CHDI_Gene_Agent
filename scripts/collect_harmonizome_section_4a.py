#!/usr/bin/env python3
"""Collect Section 4a Harmonizome associations for one gene (diagnostic CLI)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gene_dossier.config import get_settings
from gene_dossier.tools.harmonizome_section4a import collect_section_4a_harmonizome


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--harmonizome-max-curated-display", type=int, default=14)
    parser.add_argument("--harmonizome-max-predicted-display", type=int, default=25)
    parser.add_argument("--timeout", type=float, default=None)
    args = parser.parse_args(argv)
    settings = get_settings()
    if args.timeout is not None:
        object.__setattr__(settings, "http_timeout_seconds", float(args.timeout))
    result = collect_section_4a_harmonizome(
        args.gene,
        max_displayed_curated=args.harmonizome_max_curated_display,
        max_displayed_predicted=args.harmonizome_max_predicted_display,
        settings=settings,
    )
    # Drop ToolResult object for JSON dump
    printable = {k: v for k, v in result.items() if k != "tool_result"}
    tr = result.get("tool_result")
    if tr is not None:
        printable["tool_result"] = {
            "success": tr.success,
            "status_code": tr.status_code,
            "request_url": tr.request_url,
            "error_type": tr.error_type,
            "error_message": tr.error_message,
        }
    if printable.get("collection"):
        coll = dict(printable["collection"])
        coll.pop("curated_records", None)
        coll.pop("predicted_records", None)
        printable["collection"] = coll
    print(json.dumps(printable, indent=2, default=str))
    return 0 if result.get("scientific_status") == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
