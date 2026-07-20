"""Tests for source coverage reporting (no network required)."""

from __future__ import annotations

import json
from pathlib import Path

from gene_dossier import models as m
from gene_dossier.config import Settings
from gene_dossier.coverage import (
    apply_coverage_updates,
    build_and_write_coverage,
    build_coverage_for_registry,
    initial_status_for_source,
    summarize_coverage,
    write_coverage_report,
)
from gene_dossier.db import get_engine, init_db, list_source_coverage, save_dossier_run, session_scope
from gene_dossier.source_registry import get_source


def _bare_settings(**overrides) -> Settings:
    """Settings with no API keys and no .env file."""
    base = {
        "ncbi_api_key": None,
        "biogrid_accesskey": None,
        "omim_api_key": None,
        "serpapi_api_key": None,
        "openai_api_key": None,
        "anthropic_api_key": None,
    }
    base.update(overrides)
    return Settings(_env_file=None, **base)


def test_initial_status_requires_key_when_missing():
    bare = _bare_settings()
    assert initial_status_for_source(get_source("BioGRID"), bare) is m.SourceStatus.requires_key
    assert initial_status_for_source(get_source("OMIM"), bare) is m.SourceStatus.requires_key
    assert initial_status_for_source(get_source("Patents"), bare) is m.SourceStatus.requires_key


def test_initial_status_not_implemented_when_key_present():
    with_key = _bare_settings(biogrid_accesskey="present")
    assert (
        initial_status_for_source(get_source("BioGRID"), with_key)
        is m.SourceStatus.not_implemented
    )


def test_initial_status_manual_and_deferred():
    bare = _bare_settings()
    assert initial_status_for_source(get_source("Antibodies"), bare) is m.SourceStatus.manual
    assert initial_status_for_source(get_source("HDinHD"), bare) is m.SourceStatus.deferred
    assert (
        initial_status_for_source(get_source("NCBI Gene"), bare)
        is m.SourceStatus.not_implemented
    )


def test_build_coverage_includes_all_registry_sources():
    results = build_coverage_for_registry("run1", settings=_bare_settings())
    assert len(results) == 32
    summary = summarize_coverage(results)
    assert summary["total"] == 32
    assert summary["by_status"]["requires_key"] == 3
    assert summary["by_status"]["manual"] == 1
    assert summary["by_status"]["deferred"] == 6
    assert summary["by_status"]["not_implemented"] == 22

    biogrid = next(r for r in results if r.source_name == "BioGRID")
    assert biogrid.error_message is not None
    assert "BIOGRID_ACCESSKEY" in biogrid.error_message


def test_apply_coverage_updates_merges_by_name():
    baseline = build_coverage_for_registry("run1", settings=_bare_settings())
    upd = m.SourceCoverageResult(
        dossier_run_id="run1",
        source_name="ncbi gene",  # case-insensitive match
        status=m.SourceStatus.success,
        raw_artifact_path="/tmp/ncbi.json",
        evidence_record_count=2,
        report_sections_supported=["General gene information"],
    )
    merged = apply_coverage_updates(baseline, [upd])
    assert len(merged) == 32
    ncbi = next(r for r in merged if r.source_name.lower() == "ncbi gene")
    assert ncbi.status is m.SourceStatus.success
    assert ncbi.evidence_record_count == 2
    # Unrelated source unchanged.
    assert next(r for r in merged if r.source_name == "Antibodies").status is m.SourceStatus.manual


def test_write_coverage_report_md_and_json(tmp_path: Path):
    results = build_coverage_for_registry("run1", settings=_bare_settings())
    upd = m.SourceCoverageResult(
        dossier_run_id="run1",
        source_name="NCBI Gene",
        status=m.SourceStatus.success,
        raw_artifact_path="/tmp/ncbi.json",
        evidence_record_count=1,
        report_sections_supported=["General gene information"],
    )
    results = apply_coverage_updates(results, [upd])
    paths = write_coverage_report(
        results, "run1", gene_symbol="SREBF2", output_dir=tmp_path
    )
    assert paths["markdown"].name == "run1_source_coverage.md"
    assert paths["json"].name == "run1_source_coverage.json"
    md = paths["markdown"].read_text(encoding="utf-8")
    assert "Source coverage report — SREBF2" in md
    assert "NCBI Gene" in md
    assert "`success`" in md
    assert "BioGRID" in md

    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert payload["dossier_run_id"] == "run1"
    assert payload["gene_symbol"] == "SREBF2"
    assert payload["summary"]["total"] == 32
    assert len(payload["sources"]) == 32


def test_build_and_write_coverage_persists_rows(tmp_path: Path):
    engine = get_engine("sqlite://")
    init_db(engine)
    run = m.DossierRun(id="run-persist", gene_symbol="SREBF2")
    upd = m.SourceCoverageResult(
        dossier_run_id="run-persist",
        source_name="NCBI Gene",
        status=m.SourceStatus.success,
        evidence_record_count=3,
        report_sections_supported=["General gene information"],
    )
    with session_scope(engine) as session:
        save_dossier_run(session, run)
        results, paths = build_and_write_coverage(
            "run-persist",
            gene_symbol="SREBF2",
            updates=[upd],
            settings=_bare_settings(),
            output_dir=tmp_path,
            session=session,
        )
        assert len(results) == 32
        assert paths["markdown"].exists()
        rows = list_source_coverage(session, "run-persist")
        assert len(rows) == 32
        ncbi = next(r for r in rows if r.source_name == "NCBI Gene")
        assert ncbi.status is m.SourceStatus.success
        assert ncbi.evidence_record_count == 3
