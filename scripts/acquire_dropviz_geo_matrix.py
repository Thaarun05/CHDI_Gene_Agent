#!/usr/bin/env python
"""Acquire and profile the DropViz GEO GSE116470 metacell matrix.

Dataset-level acquisition: downloads once into an immutable attempt directory,
profiles the actual matrix, classifies its value semantics, and only then pins
an accepted source pointer. Re-running reuses the accepted artifact unless
``--force-refresh`` is given.

The profile is the scientific gate for Section 2c: the DropViz chart must not be
built until the matrix scale and units are established here.

Usage::

    PYTHONPATH=src python scripts/acquire_dropviz_geo_matrix.py \
        --output-root data/outputs --probe-gene Srebf2
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from gene_dossier.section_2c_sources import (
    accept_source,
    load_accepted_source,
    paths_for,
    sha256_bytes,
    write_json_atomic,
)
from gene_dossier.tools import dropviz_geo as dg

SOURCE_KEY = "gse116470_metacells"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="data/outputs")
    parser.add_argument(
        "--probe-gene",
        default=None,
        help="Optional symbol used only to confirm exact-match behaviour in the profile.",
    )
    parser.add_argument("--force-refresh", action="store_true")
    parser.add_argument(
        "--local-file",
        default=None,
        help="Use a local .csv.gz instead of downloading (offline development).",
    )
    args = parser.parse_args()

    paths = paths_for(args.output_root)
    paths.ensure()
    official_url = dg.supplementary_url()

    if not args.force_refresh:
        cached = load_accepted_source(
            paths, source_key=SOURCE_KEY, official_url=official_url
        )
        if cached:
            print(f"cache hit: {cached['artifact_path']} ({cached['byte_size']} bytes)")
            print(f"sha256={cached['sha256']}")
            print(f"value_semantics_status={cached.get('value_semantics_status')}")
            return 0

    attempt = paths.new_source_attempt(SOURCE_KEY)
    print(f"attempt dir: {attempt}")

    started = time.time()
    if args.local_file:
        content = Path(args.local_file).read_bytes()
        download_meta = {
            "retrieval_method": "local_file",
            "local_file": str(args.local_file),
            "resolved_url": None,
        }
        check = dg.validate_gzip_payload(content)
        if not check["ok"]:
            write_json_atomic(
                attempt / "manifest.json",
                {"status": "failed", "error_type": check["error_type"], **download_meta},
            )
            print(f"FAILED validation: {check['error_type']}: {check['error_message']}")
            return 2
        status_code = None
    else:
        print(f"downloading {official_url} ...")
        tr = dg.download_metacell_matrix(
            gene_symbol=args.probe_gene or "", settings=None
        )
        data = tr.data if isinstance(tr.data, dict) else {}
        if not tr.success:
            write_json_atomic(
                attempt / "manifest.json",
                {
                    "status": "failed",
                    "error_type": tr.error_type,
                    "error_message": tr.error_message,
                    "status_code": tr.status_code,
                    "request_url": tr.request_url,
                    "attempts": data.get("attempts"),
                },
            )
            print(f"FAILED download: {tr.error_type}: {tr.error_message}")
            return 2
        content = data.get("content") or b""
        status_code = tr.status_code
        download_meta = {
            "retrieval_method": "httpx_direct",
            "resolved_url": data.get("resolved_url"),
            "content_type": data.get("content_type"),
            "response_headers": data.get("response_headers"),
            "attempts": data.get("attempts"),
        }

    digest = sha256_bytes(content)
    artifact_path = attempt / dg.METACELL_FILENAME
    artifact_path.write_bytes(content)
    print(
        f"downloaded {len(content)} bytes in {time.time() - started:.1f}s sha256={digest}"
    )

    print("scanning matrix (single streaming pass) ...")
    scan_started = time.time()
    scan = dg.scan_matrix(
        dg.open_matrix_stream(content), target_gene=args.probe_gene
    )
    if not scan.ok:
        write_json_atomic(
            attempt / "manifest.json",
            {
                "status": "failed",
                "error_type": scan.error_type,
                "error_message": scan.error_message,
                "sha256": digest,
                **download_meta,
            },
        )
        print(f"FAILED scan: {scan.error_type}: {scan.error_message}")
        return 3

    semantics = dg.classify_value_semantics(
        scan,
        documentation_reference=(
            "NCBI GEO GSE116470 series record and supplementary file description"
        ),
    )
    profile = dg.build_matrix_profile(
        scan,
        semantics=semantics,
        source_sha256=digest,
        source_url=official_url,
    )
    write_json_atomic(attempt / "dropviz_geo_matrix_profile.json", profile)
    print(f"scan completed in {time.time() - scan_started:.1f}s")

    manifest = {
        "status": "validated",
        "source_key": SOURCE_KEY,
        "accession": dg.GEO_ACCESSION,
        "supplementary_filename": dg.METACELL_FILENAME,
        "official_url": official_url,
        "status_code": status_code,
        "sha256": digest,
        "byte_size": len(content),
        "artifact_path": str(artifact_path),
        "probe_gene": args.probe_gene,
        "probe_gene_match_count": scan.target_matches,
        "calculation_version": dg.CALCULATION_VERSION,
        **download_meta,
    }
    write_json_atomic(attempt / "manifest.json", manifest)

    accept_source(
        paths,
        source_key=SOURCE_KEY,
        attempt_dir=attempt,
        artifact_path=artifact_path,
        official_url=official_url,
        sha256=digest,
        byte_size=len(content),
        validation={
            "gene_row_count": profile["gene_row_count"],
            "population_column_count": profile["population_column_count"],
            "shape_matches_expected": profile["shape_matches_expected"],
        },
        extra={
            "value_semantics_status": semantics["value_semantics_status"],
            "value_semantics_basis": semantics["basis"],
            "matrix_profile_path": str(attempt / "dropviz_geo_matrix_profile.json"),
        },
    )

    print("--- matrix profile ---")
    print(
        json.dumps(
            {
                k: profile[k]
                for k in (
                    "gene_row_count",
                    "population_column_count",
                    "shape_matches_expected",
                    "duplicate_gene_symbols",
                    "malformed_row_count",
                    "value_semantics_status",
                    "value_semantics_basis",
                )
            },
            indent=2,
        )
    )
    ev = profile["value_semantics_evidence"]
    print(
        json.dumps(
            {
                "integer_fraction": ev["integer_fraction"],
                "minimum": ev["minimum"],
                "maximum": ev["maximum"],
                "negative_value_count": ev["negative_value_count"],
                "column_sum_distribution": ev["column_sum_distribution"],
            },
            indent=2,
        )
    )
    if args.probe_gene:
        print(f"probe gene {args.probe_gene!r} exact matches: {scan.target_matches}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
