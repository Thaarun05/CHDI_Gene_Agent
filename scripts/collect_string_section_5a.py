#!/usr/bin/env python3
"""Collect Section 5a STRING network for one gene (diagnostic CLI)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gene_dossier.config import get_settings
from gene_dossier.tools import string_db as sd
from gene_dossier.section_5a import Section5aConfig, canonicalize_network


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--string-required-score", type=int, default=400)
    parser.add_argument("--timeout", type=float, default=None)
    args = parser.parse_args(argv)
    settings = get_settings()
    if args.timeout is not None:
        object.__setattr__(settings, "http_timeout_seconds", float(args.timeout))
    cfg = Section5aConfig(required_score=args.string_required_score)

    resolved = sd.resolve_string_identifier(
        args.gene,
        species=cfg.species_taxon_id,
        settings=settings,
        caller_identity=sd.SECTION_5A_CALLER_IDENTITY,
    )
    out: dict = {
        "resolve_success": resolved.success,
        "resolve_error": resolved.error_type,
        "resolve_data": {
            k: v
            for k, v in dict(resolved.data or {}).items()
            if k not in {"raw", "_string_meta"}
        }
        if isinstance(resolved.data, dict)
        else resolved.data,
    }
    if not resolved.success:
        print(json.dumps(out, indent=2, default=str))
        return 1
    string_id = (resolved.data or {}).get("string_id")
    network = sd.fetch_network(
        str(string_id),
        gene_symbol=args.gene,
        species=cfg.species_taxon_id,
        add_nodes=cfg.add_nodes,
        required_score=cfg.required_score,
        network_type=cfg.network_type,
        settings=settings,
        caller_identity=sd.SECTION_5A_CALLER_IDENTITY,
    )
    rows = sd.extract_network_rows(network.data) if network.success else []
    canon = canonicalize_network(
        rows,
        query_string_id=str(string_id),
        species_taxon_id=cfg.species_taxon_id,
        required_score=cfg.required_score,
    )
    out.update(
        {
            "network_success": network.success,
            "network_error": network.error_type,
            "request_params": network.request_params,
            "stats": canon.get("stats"),
            "warning_count": len(canon.get("warnings") or []),
        }
    )
    print(json.dumps(out, indent=2, default=str))
    return 0 if network.success else 1


if __name__ == "__main__":
    raise SystemExit(main())
