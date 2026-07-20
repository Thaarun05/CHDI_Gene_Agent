"""Tests for the provenance database layer (in-memory SQLite; no network)."""

from __future__ import annotations

import pytest
from sqlalchemy import inspect

from gene_dossier import models as m
from gene_dossier.db import (
    SourceCoverageResultRow,
    get_engine,
    get_evidence_by_source_id,
    init_db,
    is_postgres_url,
    is_sqlite_url,
    list_evidence_for_run,
    list_source_coverage,
    normalize_database_url,
    save_api_run,
    save_dossier_run,
    save_evidence_record,
    save_raw_artifact,
    save_source_coverage,
    session_scope,
)


EXPECTED_TABLES = {
    "dossier_runs",
    "api_runs",
    "raw_artifacts",
    "evidence_records",
    "report_sections",
    "claims",
    "verification_results",
    "source_coverage_results",
}


@pytest.fixture
def engine():
    eng = get_engine("sqlite://")
    init_db(eng)
    return eng


def test_url_helpers_sqlite_and_postgres():
    assert is_sqlite_url("sqlite:///data/gene_dossier.db")
    assert is_sqlite_url("sqlite://")
    assert not is_sqlite_url("postgresql+psycopg://u:p@h/db")

    assert is_postgres_url("postgresql://u:p@h/db")
    assert is_postgres_url("postgres://u:p@h/db")
    assert is_postgres_url("postgresql+psycopg://u:p@h/db")
    assert not is_postgres_url("sqlite:///x.db")


def test_normalize_database_url_uses_psycopg():
    assert (
        normalize_database_url("postgresql://u:p@host/db")
        == "postgresql+psycopg://u:p@host/db"
    )
    assert (
        normalize_database_url("postgres://u:p@host/db")
        == "postgresql+psycopg://u:p@host/db"
    )
    assert normalize_database_url("sqlite:///data/x.db") == "sqlite:///data/x.db"
    assert (
        normalize_database_url("postgresql+psycopg://u:p@host/db")
        == "postgresql+psycopg://u:p@host/db"
    )


def test_init_db_creates_all_tables(engine):
    names = set(inspect(engine).get_table_names())
    assert EXPECTED_TABLES <= names
    assert SourceCoverageResultRow.__tablename__ == "source_coverage_results"


def test_dossier_run_and_evidence_round_trip(engine):
    run = m.DossierRun(gene_symbol="SREBF2", notes="unit")
    api = m.ApiRun(
        dossier_run_id=run.id,
        gene_symbol="SREBF2",
        source_name="NCBI Gene",
        endpoint_name="esearch",
        request_url="https://example.com",
        request_params={"db": "gene"},
        success=True,
        status_code=200,
    )
    art = m.RawArtifact(
        dossier_run_id=run.id,
        api_run_id=api.id,
        source_name="NCBI Gene",
        artifact_type="json",
        file_path="/tmp/fake.json",
        content_hash="abc123",
    )
    ev = m.EvidenceRecord(
        source_id="ncbi-gene:srebf2:gene_identity:6721",
        dossier_run_id=run.id,
        gene_symbol="SREBF2",
        section="General gene information",
        source_name="NCBI Gene",
        source_type=m.SourceType.curated_database,
        assertion_type=m.AssertionType.gene_identity,
        fact_type="entrez_id",
        evidence_grade=m.EvidenceGrade.C,
        value={"entrez_id": "6721"},
        display_text="SREBF2 Entrez Gene ID is 6721.",
        api_run_id=api.id,
        raw_artifact_id=art.id,
    )

    with session_scope(engine) as session:
        save_dossier_run(session, run)
        save_api_run(session, api)
        save_raw_artifact(session, art)
        save_evidence_record(session, ev)

    with session_scope(engine) as session:
        rows = list_evidence_for_run(session, run.id)
        assert len(rows) == 1
        assert rows[0].evidence_grade is m.EvidenceGrade.C
        assert rows[0].source_type is m.SourceType.curated_database
        assert rows[0].value["entrez_id"] == "6721"
        found = get_evidence_by_source_id(session, ev.source_id, run.id)
        assert found is not None
        assert found.display_text.startswith("SREBF2")


def test_source_coverage_round_trip_and_upsert(engine):
    run = m.DossierRun(gene_symbol="SREBF2")
    cov = m.SourceCoverageResult(
        dossier_run_id=run.id,
        source_name="NCBI Gene",
        status=m.SourceStatus.success,
        raw_artifact_path="/tmp/x.json",
        evidence_record_count=1,
        report_sections_supported=["General gene information"],
    )
    missing_key = m.SourceCoverageResult(
        dossier_run_id=run.id,
        source_name="BioGRID",
        status=m.SourceStatus.requires_key,
        error_message="BIOGRID_ACCESSKEY missing",
    )

    with session_scope(engine) as session:
        save_dossier_run(session, run)
        save_source_coverage(session, cov)
        save_source_coverage(session, missing_key)

    with session_scope(engine) as session:
        rows = list_source_coverage(session, run.id)
        assert len(rows) == 2
        by_name = {r.source_name: r for r in rows}
        assert by_name["NCBI Gene"].status is m.SourceStatus.success
        assert by_name["BioGRID"].status is m.SourceStatus.requires_key
        # Upsert same run+source updates count without duplicating.
        save_source_coverage(
            session, cov.model_copy(update={"evidence_record_count": 3})
        )

    with session_scope(engine) as session:
        rows = list_source_coverage(session, run.id)
        assert len(rows) == 2
        ncbi = next(r for r in rows if r.source_name == "NCBI Gene")
        assert ncbi.evidence_record_count == 3


def test_raw_artifact_row_stores_metadata_not_payload(engine):
    """DB stores path/hash only — never the raw response body."""
    run = m.DossierRun(gene_symbol="SREBF2")
    art = m.RawArtifact(
        dossier_run_id=run.id,
        source_name="UniProt",
        artifact_type="json",
        file_path="/data/raw/run/uniprot/x.json",
        content_hash="deadbeef",
        original_url="https://rest.uniprot.org/example",
    )
    dump = art.model_dump()
    assert "file_path" in dump and "content_hash" in dump
    assert "data" not in dump and "body" not in dump and "content" not in dump

    with session_scope(engine) as session:
        save_dossier_run(session, run)
        save_raw_artifact(session, art)
