#!/usr/bin/env python3
"""Refresh UCSC conservation evidence for an existing dossier run (Section 1b).

Does not rerun unrelated sources or invoke an LLM.

Approved SREBF2 import path::

    PYTHONPATH=src .venv/bin/python scripts/refresh_ucsc_conservation.py \\
      --gene SREBF2 \\
      --dossier-run-id 9923bf6d326246a4a9abb6d56053c898 \\
      --search-json tests/fixtures/ucsc/srebf2_search_relevant.json \\
      --track-json tests/fixtures/ucsc/srebf2_known_gene_region.json \\
      --figure-file tests/fixtures/ucsc/srebf2_comprehensive_conservation.png
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.config import get_settings  # noqa: E402
from gene_dossier.db import (  # noqa: E402
    get_dossier_run,
    get_evidence_by_source_id,
    list_evidence_for_run,
    save_api_run,
    save_evidence_record,
    save_raw_artifact,
    session_scope,
)
from gene_dossier.models import ApiRun, RawArtifact  # noqa: E402
from gene_dossier.normalize.ucsc_conservation import (  # noqa: E402
    build_conservation_evidence,
)
from gene_dossier.raw_store import RawStore, compute_hash  # noqa: E402
from gene_dossier.ucsc_figure import (  # noqa: E402
    stage_and_commit_figure,
    validate_image_bytes,
)
from gene_dossier.ucsc_parse import parse_search_response  # noqa: E402


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"Expected JSON object in {path}")
    return data


def _upsert_evidence(session, record) -> str:
    existing = get_evidence_by_source_id(
        session, record.source_id, dossier_run_id=record.dossier_run_id
    )
    if existing is not None:
        record.id = existing.id
        record.created_at = existing.created_at
    save_evidence_record(session, record)
    return record.id


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Refresh UCSC conservation evidence for Section 1b."
    )
    parser.add_argument("--gene", required=True)
    parser.add_argument("--dossier-run-id", required=True)
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--search-json", type=Path)
    parser.add_argument("--track-json", type=Path)
    parser.add_argument("--figure-file", type=Path)
    parser.add_argument("--fetch-search", action="store_true")
    parser.add_argument("--fetch-track", action="store_true")
    parser.add_argument("--fetch-figure", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--render-preview",
        action="store_true",
        help="Regenerate Section 1b preview after refresh",
    )
    args = parser.parse_args(argv)

    gene = args.gene.strip()
    run_id = args.dossier_run_id.strip()
    genome = args.genome.strip() or "hg38"
    settings = get_settings()

    with session_scope() as session:
        run = get_dossier_run(session, run_id)
        if run is None:
            raise SystemExit(f"No dossier run {run_id!r}")
        if (run.gene_symbol or "").strip().upper() != gene.upper():
            raise SystemExit(
                f"Gene mismatch: run gene={run.gene_symbol!r} requested={gene!r}"
            )

    # --- Validate / load inputs before any DB writes ---
    search_payload: dict[str, Any] | None = None
    track_payload: dict[str, Any] | None = None
    figure_bytes: bytes | None = None

    if args.search_json:
        search_payload = _load_json(args.search_json)
    if args.track_json:
        track_payload = _load_json(args.track_json)
    if args.figure_file:
        figure_bytes = args.figure_file.read_bytes()
        validated, err = validate_image_bytes(figure_bytes)
        if err or validated is None:
            raise SystemExit(f"Invalid figure file: {err.message if err else 'unknown'}")

    if args.fetch_search or args.fetch_track or args.fetch_figure:
        from gene_dossier.tools import ucsc as ucsc_client

        if args.fetch_search or (args.fetch_track and search_payload is None):
            searched = ucsc_client.search(gene, genome=genome, settings=settings)
            if not searched.success or not isinstance(searched.data, dict):
                raise SystemExit(f"UCSC search failed: {searched.error_message}")
            search_payload = searched.data
        if args.fetch_track:
            inv = parse_search_response(
                search_payload, gene_symbol=gene, genome=genome
            )
            display = inv.selected_display_interval
            if display is None:
                raise SystemExit("Could not resolve display locus from search")
            track_res = ucsc_client.get_track_data(
                display.chrom,
                display.api_start_0_based,
                display.api_end_exclusive,
                gene_symbol=gene,
                genome=genome,
                settings=settings,
            )
            if not track_res.success or not isinstance(track_res.data, dict):
                raise SystemExit(f"UCSC track failed: {track_res.error_message}")
            track_payload = track_res.data
        if args.fetch_figure:
            print(
                "Live --fetch-figure is implemented only as a non-blocking optional "
                "path; use --figure-file for the approved SREBF2 deliverable."
            )
            if figure_bytes is None:
                print("WARNING: no figure imported (live fetch deferred/partial).")

    if search_payload is None and track_payload is None:
        raise SystemExit("Provide --search-json/--track-json and/or fetch flags")

    # Pre-build evidence (validation) before writing.
    figure_value: dict[str, Any] | None = None
    staged_rel: str | None = None
    staged_abs: Path | None = None
    staged_hash: str | None = None
    if figure_bytes is not None:
        # Stage to temp managed path before DB transaction; commit path after.
        # For dry-run, validate only.
        if not args.dry_run:
            staged_abs, staged_rel, staged_hash = stage_and_commit_figure(
                dossier_run_id=run_id,
                content=figure_bytes,
                extension=args.figure_file.suffix.lstrip(".") or "png",
                settings=settings,
            )
        else:
            validated, _ = validate_image_bytes(figure_bytes)
            assert validated is not None
            staged_hash = validated.sha256
            staged_rel = f"{run_id}/ucsc/figures/{validated.sha256}.png"

        # Fill figure_value after we know locus/transcript (below). Temporary stub.
        figure_value = {
            "sha256": staged_hash,
            "relative_path": staged_rel,
            "local_artifact_path": staged_rel,
            "media_type": "image/png",
            "retrieval_method": "attached_validated_ucsc_render",
            "origin_endpoint": "hgRenderTracks",
            "origin_generation_confirmed_by_user": True,
            "api_key_used_to_generate": True,
            "api_key_used": True,
            "api_key_persisted": False,
            "byte_size": len(figure_bytes),
        }
        validated, _ = validate_image_bytes(figure_bytes)
        if validated:
            figure_value["width"] = validated.width
            figure_value["height"] = validated.height
            figure_value["media_type"] = validated.media_type

    # Peek selection to enrich figure metadata
    records, diagnostics = build_conservation_evidence(
        dossier_run_id=run_id,
        gene_symbol=gene,
        genome=genome,
        search_payload=search_payload,
        track_payload=track_payload,
        figure_value=None,
    )
    locus = next((r for r in records if r.fact_type == "ucsc_gene_locus"), None)
    transcript = next(
        (r for r in records if r.fact_type == "ucsc_canonical_transcript"), None
    )
    if figure_value is not None:
        if locus and isinstance(locus.value, dict):
            figure_value["genome"] = locus.value.get("genome")
            figure_value["display_position"] = locus.value.get("display_position")
        if transcript and isinstance(transcript.value, dict):
            figure_value["selected_transcript"] = transcript.value.get("transcript_id")
            figure_value.setdefault(
                "display_position", transcript.value.get("display_position")
            )
        records, diagnostics = build_conservation_evidence(
            dossier_run_id=run_id,
            gene_symbol=gene,
            genome=genome,
            search_payload=search_payload,
            track_payload=track_payload,
            figure_value=figure_value,
        )

    print(f"Prepared {len(records)} UCSC EvidenceRecords")
    for d in diagnostics:
        print(f"  diagnostic[{d.get('severity')}]: {d.get('code')}: {d.get('message')}")
    for rec in records:
        print(f"  fact={rec.fact_type} source_id={rec.source_id}")

    if args.dry_run:
        print("Dry-run complete; no database writes.")
        return 0

    store = RawStore(base_dir=settings.raw_data_path)
    created_ids: list[str] = []

    with session_scope() as session:
        search_api = None
        search_art = None
        if search_payload is not None:
            search_api = ApiRun(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name="UCSC",
                endpoint_name="search_import",
                method="IMPORT",
                request_url="import://ucsc/search",
                request_params={"genome": genome, "search": gene},
                status_code=None,
                success=True,
            )
            search_art = store.save_json(
                run_id,
                "UCSC",
                search_payload,
                api_run_id=search_api.id,
                filename_hint="search",
                notes="UCSC search import",
            )
            # Store portable relative path when under raw root.
            try:
                from gene_dossier.ucsc_figure import relative_to_artifact_root

                search_art.file_path = relative_to_artifact_root(Path(search_art.file_path))
            except ValueError:
                pass
            search_api.raw_artifact_id = search_art.id
            save_api_run(session, search_api)
            save_raw_artifact(session, search_art)

        track_api = None
        track_art = None
        if track_payload is not None:
            track_api = ApiRun(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name="UCSC",
                endpoint_name="knownGene_import",
                method="IMPORT",
                request_url="import://ucsc/knownGene",
                request_params={"genome": genome, "track": "knownGene"},
                status_code=None,
                success=True,
            )
            track_art = store.save_json(
                run_id,
                "UCSC",
                track_payload,
                api_run_id=track_api.id,
                filename_hint="knownGene",
                notes="UCSC knownGene import",
            )
            try:
                from gene_dossier.ucsc_figure import relative_to_artifact_root

                track_art.file_path = relative_to_artifact_root(Path(track_art.file_path))
            except ValueError:
                pass
            track_api.raw_artifact_id = track_art.id
            save_api_run(session, track_api)
            save_raw_artifact(session, track_art)

        figure_api = None
        figure_art = None
        if figure_value is not None and staged_abs is not None and staged_rel is not None:
            figure_api = ApiRun(
                dossier_run_id=run_id,
                gene_symbol=gene,
                source_name="UCSC",
                endpoint_name="attached_figure_import",
                method="IMPORT",
                request_url="import://ucsc/attached_figure",
                request_params={
                    "retrieval_method": "attached_validated_ucsc_render",
                    "origin_endpoint": "hgRenderTracks",
                    "relative_path": staged_rel,
                    "sha256": staged_hash,
                },
                status_code=None,
                success=True,
            )
            figure_art = RawArtifact(
                dossier_run_id=run_id,
                api_run_id=figure_api.id,
                source_name="UCSC",
                artifact_type="image",
                file_path=staged_rel,
                original_url=None,
                content_hash=staged_hash or compute_hash(figure_bytes or b""),
                notes="attached_validated_ucsc_render",
            )
            figure_api.raw_artifact_id = figure_art.id
            save_api_run(session, figure_api)
            save_raw_artifact(session, figure_art)

        # Rebuild with provenance IDs
        records, diagnostics = build_conservation_evidence(
            dossier_run_id=run_id,
            gene_symbol=gene,
            genome=genome,
            search_payload=search_payload,
            track_payload=track_payload,
            figure_value=figure_value,
            search_api_run_id=search_api.id if search_api else None,
            search_raw_artifact_id=search_art.id if search_art else None,
            track_api_run_id=track_api.id if track_api else None,
            track_raw_artifact_id=track_art.id if track_art else None,
            figure_api_run_id=figure_api.id if figure_api else None,
            figure_raw_artifact_id=figure_art.id if figure_art else None,
        )
        for rec in records:
            eid = _upsert_evidence(session, rec)
            created_ids.append(eid)

        # Supersede prior conservation-figure evidence when the SHA changes.
        figure_rec = next(
            (r for r in records if r.fact_type == "ucsc_conservation_figure"), None
        )
        if figure_rec is not None:
            from gene_dossier.db import delete_evidence_record

            for existing in list_evidence_for_run(session, run_id):
                if (
                    existing.source_name == "UCSC"
                    and existing.fact_type == "ucsc_conservation_figure"
                    and existing.source_id != figure_rec.source_id
                ):
                    delete_evidence_record(session, existing.id)
                    print(
                        f"Superseded stale figure EvidenceRecord {existing.id} "
                        f"({existing.source_id})"
                    )

        ucsc_count = sum(
            1
            for r in list_evidence_for_run(session, run_id)
            if r.source_name == "UCSC"
            and r.fact_type.startswith("ucsc_")
        )
        print(f"UCSC conservation EvidenceRecords for run: {ucsc_count}")
        print("Created/updated EvidenceRecord IDs:")
        for eid in created_ids:
            print(f"  {eid}")

    if args.render_preview:
        import subprocess

        cmd = [
            sys.executable,
            str(_REPO_ROOT / "scripts" / "render_rancho_section_preview.py"),
            "--gene",
            gene,
            "--section",
            "1b",
            "--dossier-run-id",
            run_id,
        ]
        return subprocess.call(cmd, cwd=str(_REPO_ROOT), env={**dict(**{k: v for k, v in __import__("os").environ.items()}), "PYTHONPATH": str(_SRC)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
