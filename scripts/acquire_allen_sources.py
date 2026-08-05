#!/usr/bin/env python
"""Acquire and pin the Allen Cell Types dataset-level sources for Section 2c.

Human M1 trimmed means and taxonomy are downloaded from the official Allen
bucket; the mouse cortex/hippocampus bundle is registered from local files under
canonical names. Each validated source is written into an immutable attempt
directory and pinned by an accepted pointer, so later gene runs reuse it and
perform no repeated download unless ``--force-refresh`` is given.

Usage::

    PYTHONPATH=src python scripts/acquire_allen_sources.py \
        --output-root data/outputs \
        --mouse-trimmed-means response.xls --mouse-taxonomy response.json
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
    write_json_atomic,
)
from gene_dossier.tools import allen_celltypes as ac

HUMAN_ATTEMPT_KEY = "allen_human_m1"
MOUSE_ATTEMPT_KEY = "mouse_ctx_hpf"


def acquire_human(paths, *, force_refresh: bool) -> int:
    failures = 0
    attempt = None
    for source_key, cache_key in ac.HUMAN_CACHE_KEYS.items():
        official_url = ac.human_m1_source_url(source_key)
        if not force_refresh:
            cached = load_accepted_source(
                paths, source_key=cache_key, official_url=official_url
            )
            if cached:
                print(f"[human/{source_key}] cache hit sha256={cached['sha256']}")
                continue

        if attempt is None:
            attempt = paths.new_source_attempt(HUMAN_ATTEMPT_KEY)
            print(f"[human] attempt dir: {attempt}")

        print(f"[human/{source_key}] downloading {official_url}")
        started = time.time()
        tr = ac.download_human_m1_source(source_key)
        data = tr.data if isinstance(tr.data, dict) else {}
        if not tr.success:
            failures += 1
            write_json_atomic(
                attempt / f"{cache_key}_failure.json",
                {
                    "status": "failed",
                    "source_key": cache_key,
                    "error_type": tr.error_type,
                    "error_message": tr.error_message,
                    "status_code": tr.status_code,
                    "request_url": tr.request_url,
                },
            )
            print(f"[human/{source_key}] FAILED {tr.error_type}: {tr.error_message}")
            continue

        content = data.get("content") or b""
        filename = Path(official_url).name
        artifact_path = attempt / filename
        artifact_path.write_bytes(content)
        print(
            f"[human/{source_key}] {len(content)} bytes in {time.time() - started:.1f}s "
            f"sha256={data.get('sha256')}"
        )

        validation = {"byte_size": len(content)}
        if source_key == "trimmed_means":
            labels = ac.matrix_celltype_labels(ac.text_lines(content))
            validation["celltype_column_count"] = len(labels)
            validation["first_columns"] = labels[:5]
        else:
            taxonomy = ac.parse_dendrogram(json.loads(content.decode("utf-8")))
            validation["taxonomy_leaf_count"] = taxonomy.leaf_count
            validation["named_internal_node_count"] = taxonomy.internal_label_count

        accept_source(
            paths,
            source_key=cache_key,
            attempt_dir=attempt,
            artifact_path=artifact_path,
            official_url=official_url,
            sha256=str(data.get("sha256")),
            byte_size=len(content),
            validation=validation,
            extra={"dataset": ac.DATASET_HUMAN_M1, "retrieval_method": "httpx_direct"},
        )
        print(f"[human/{source_key}] accepted: {json.dumps(validation)[:200]}")
    return failures


def acquire_mouse(paths, *, trimmed_means: str, taxonomy_path: str, force_refresh: bool) -> int:
    specs = (
        (ac.CACHE_KEY_MOUSE_TRIMMED_MEANS, trimmed_means, ac.MOUSE_TRIMMED_MEANS_FILENAME),
        (ac.CACHE_KEY_MOUSE_TAXONOMY, taxonomy_path, ac.MOUSE_TAXONOMY_FILENAME),
    )
    failures = 0
    attempt = None
    for cache_key, source_path, canonical_name in specs:
        if not force_refresh:
            cached = load_accepted_source(paths, source_key=cache_key)
            if cached:
                print(f"[mouse/{cache_key}] cache hit sha256={cached['sha256']}")
                continue

        if attempt is None:
            attempt = paths.new_source_attempt(MOUSE_ATTEMPT_KEY)
            print(f"[mouse] attempt dir: {attempt}")

        registered = ac.register_local_source(source_path, canonical_name=canonical_name)
        if not registered["ok"]:
            failures += 1
            print(f"[mouse/{cache_key}] FAILED {registered['error_type']}")
            continue

        content = registered["content"]
        artifact_path = attempt / canonical_name
        artifact_path.write_bytes(content)

        validation = {"byte_size": registered["byte_size"]}
        if cache_key == ac.CACHE_KEY_MOUSE_TRIMMED_MEANS:
            labels = ac.matrix_celltype_labels(ac.text_lines(content))
            validation["expression_cluster_count"] = len(labels)
        else:
            parsed = ac.parse_dendrogram(json.loads(content.decode("utf-8")))
            validation["taxonomy_leaf_count"] = parsed.leaf_count
            validation["named_internal_node_count"] = parsed.internal_label_count

        accept_source(
            paths,
            source_key=cache_key,
            attempt_dir=attempt,
            artifact_path=artifact_path,
            official_url="",
            sha256=registered["sha256"],
            byte_size=registered["byte_size"],
            validation=validation,
            extra={
                "dataset": ac.DATASET_MOUSE_CTX_HPF,
                "retrieval_method": "local_registration",
                "original_filename": registered["original_filename"],
                "original_path": registered["original_path"],
            },
        )
        print(
            f"[mouse/{cache_key}] accepted from {registered['original_filename']} "
            f"-> {canonical_name} {json.dumps(validation)[:160]}"
        )
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="data/outputs")
    parser.add_argument("--mouse-trimmed-means", default="response.xls")
    parser.add_argument("--mouse-taxonomy", default="response.json")
    parser.add_argument("--force-refresh", action="store_true")
    parser.add_argument("--skip-human", action="store_true")
    parser.add_argument("--skip-mouse", action="store_true")
    args = parser.parse_args()

    paths = paths_for(args.output_root)
    paths.ensure()

    failures = 0
    if not args.skip_human:
        failures += acquire_human(paths, force_refresh=args.force_refresh)
    if not args.skip_mouse:
        failures += acquire_mouse(
            paths,
            trimmed_means=args.mouse_trimmed_means,
            taxonomy_path=args.mouse_taxonomy,
            force_refresh=args.force_refresh,
        )

    print(f"done with {failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
