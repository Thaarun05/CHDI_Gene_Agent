"""Provenance database (SQLModel) — SQLite or Supabase Postgres via DATABASE_URL.

Persists dossier runs, API call metadata, raw-artifact *pointers*, evidence records,
report sections, claims, verification results, and source coverage rows. Domain models
in ``models.py`` stay storage-agnostic; this module is the table layer.

Truth hierarchy:

- Raw artifact *bytes* live on disk (``raw_store``); later Supabase Storage / S3.
- This database stores structured provenance (Postgres preferred for shared/dev;
  SQLite for local/offline tests). Never store large raw API response bodies here.
- Chroma indexes evidence for semantic search only — it is not the source of truth.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
import sqlite3
from typing import Any, Optional
from urllib.parse import urlparse

from sqlalchemy import Column, JSON, UniqueConstraint, create_engine, event
from sqlmodel import Field, Session, SQLModel, select

from . import models as domain
from .config import PROJECT_ROOT, get_settings
from .models import new_id
from .source_ids import slugify


# --------------------------------------------------------------------------------------
# Table models
# --------------------------------------------------------------------------------------
class DossierRunRow(SQLModel, table=True):
    """Persisted dossier generation run."""

    __tablename__ = "dossier_runs"

    id: str = Field(primary_key=True)
    gene_symbol: str = Field(index=True)
    official_symbol: Optional[str] = None
    run_type: str = "full_dossier"
    focus: Optional[str] = None
    status: str = Field(default="created", index=True)
    started_at: datetime
    completed_at: Optional[datetime] = None
    config: Optional[dict[str, Any]] = Field(default=None, sa_column=Column(JSON))
    notes: Optional[str] = None


class ApiRunRow(SQLModel, table=True):
    """Persisted metadata for one API call."""

    __tablename__ = "api_runs"

    id: str = Field(primary_key=True)
    dossier_run_id: str = Field(index=True, foreign_key="dossier_runs.id")
    gene_symbol: str = Field(index=True)
    source_name: str = Field(index=True)
    endpoint_name: str
    method: str = "GET"
    request_url: str
    request_params: dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    status_code: Optional[int] = None
    success: bool = False
    error_type: Optional[str] = None
    error_message: Optional[str] = None
    retrieved_at: datetime
    raw_artifact_id: Optional[str] = Field(default=None, index=True)


class RawArtifactRow(SQLModel, table=True):
    """Persisted metadata pointer to a raw artifact (path + hash — not the bytes)."""

    __tablename__ = "raw_artifacts"

    id: str = Field(primary_key=True)
    dossier_run_id: str = Field(index=True, foreign_key="dossier_runs.id")
    api_run_id: Optional[str] = Field(default=None, index=True)
    source_name: str = Field(index=True)
    artifact_type: str
    file_path: str
    original_url: Optional[str] = None
    content_hash: str = Field(index=True)
    captured_at: datetime
    notes: Optional[str] = None


class EvidenceRecordRow(SQLModel, table=True):
    """Persisted normalized evidence fact."""

    __tablename__ = "evidence_records"

    id: str = Field(primary_key=True)
    source_id: str = Field(index=True)
    dossier_run_id: str = Field(index=True, foreign_key="dossier_runs.id")
    gene_symbol: str = Field(index=True)
    official_symbol: Optional[str] = None
    section: str = Field(index=True)
    subsection: Optional[str] = None
    source_name: str = Field(index=True)
    source_type: str
    assertion_type: str = Field(index=True)
    fact_type: str
    organism: Optional[str] = None
    species: Optional[str] = None
    taxon_id: Optional[int] = None
    evidence_grade: str = Field(index=True)
    manual_review_required: bool = False
    confidence_notes: Optional[str] = None
    value: dict[str, Any] = Field(default_factory=dict, sa_column=Column(JSON))
    display_text: str
    api_run_id: Optional[str] = None
    raw_artifact_id: Optional[str] = Field(default=None, index=True)
    raw_response_pointer: Optional[str] = None
    created_at: datetime


class ReportSectionRow(SQLModel, table=True):
    """Persisted report section with cited source_ids."""

    __tablename__ = "report_sections"

    id: str = Field(primary_key=True)
    dossier_run_id: str = Field(index=True, foreign_key="dossier_runs.id")
    section_name: str
    subsection_name: Optional[str] = None
    content_markdown: str = ""
    source_ids: list[str] = Field(default_factory=list, sa_column=Column(JSON))
    figure_ids: list[str] = Field(default_factory=list, sa_column=Column(JSON))
    table_ids: list[str] = Field(default_factory=list, sa_column=Column(JSON))
    status: str = "draft"
    created_at: datetime
    updated_at: datetime


class ClaimRow(SQLModel, table=True):
    """Persisted factual claim citing source_ids."""

    __tablename__ = "claims"

    id: str = Field(primary_key=True)
    dossier_run_id: str = Field(index=True, foreign_key="dossier_runs.id")
    section_id: Optional[str] = None
    claim_text: str
    source_ids: list[str] = Field(default_factory=list, sa_column=Column(JSON))
    evidence_grade: Optional[str] = None
    claim_type: Optional[str] = None
    created_at: datetime


class VerificationResultRow(SQLModel, table=True):
    """Persisted verification outcome for one claim."""

    __tablename__ = "verification_results"

    id: str = Field(primary_key=True)
    claim_id: str = Field(index=True, foreign_key="claims.id")
    source_id_presence_passed: bool
    source_exists_passed: bool
    semantic_support: str = "pass"
    causal_language_check: str = "pass"
    evidence_strength_check: str = "pass"
    verdict: str = "pass"
    reason: Optional[str] = None
    needs_human_review: bool = False
    created_at: datetime


class SourceCoverageResultRow(SQLModel, table=True):
    """Persisted per-source coverage line for a dossier run."""

    __tablename__ = "source_coverage_results"
    __table_args__ = (
        UniqueConstraint("dossier_run_id", "source_name", name="uq_coverage_run_source"),
    )

    id: str = Field(primary_key=True)
    dossier_run_id: str = Field(index=True, foreign_key="dossier_runs.id")
    source_name: str = Field(index=True)
    status: str = Field(index=True)
    raw_artifact_path: Optional[str] = None
    evidence_record_count: int = 0
    error_message: Optional[str] = None
    report_sections_supported: list[str] = Field(default_factory=list, sa_column=Column(JSON))
    notes: Optional[str] = None


class GeneratedReportRow(SQLModel, table=True):
    """Durable metadata pointer for one generated dossier report."""

    __tablename__ = "generated_reports"
    __table_args__ = (
        UniqueConstraint("dossier_run_id", name="uq_generated_report_dossier_run"),
    )

    id: str = Field(primary_key=True)
    dossier_run_id: str = Field(index=True, foreign_key="dossier_runs.id")
    gene_symbol: str = Field(index=True)
    title: str
    status: str = Field(index=True)
    created_at: datetime = Field(index=True)
    sections: list[str] = Field(default_factory=list, sa_column=Column(JSON))
    html_path: str
    pdf_path: Optional[str] = None
    job_id: Optional[str] = Field(default=None, index=True)


# --------------------------------------------------------------------------------------
# Engine / session
# --------------------------------------------------------------------------------------
def is_sqlite_url(database_url: str) -> bool:
    """Return True if ``database_url`` targets SQLite."""
    return database_url.strip().lower().startswith("sqlite:")


def is_postgres_url(database_url: str) -> bool:
    """Return True if ``database_url`` targets Postgres (incl. Supabase)."""
    scheme = urlparse(database_url.strip()).scheme.lower()
    return scheme in {"postgresql", "postgres", "postgresql+psycopg", "postgresql+psycopg2"}


def normalize_database_url(database_url: str) -> str:
    """Normalize URL schemes for SQLAlchemy.

    - ``postgres://`` → ``postgresql+psycopg://``
    - ``postgresql://`` → ``postgresql+psycopg://`` (explicit psycopg3 driver)
    - SQLite URLs are returned unchanged
    """
    url = database_url.strip()
    lower = url.lower()
    if lower.startswith("postgres://"):
        return "postgresql+psycopg://" + url.split("://", 1)[1]
    if lower.startswith("postgresql://"):
        return "postgresql+psycopg://" + url.split("://", 1)[1]
    return url


def _sqlite_path_from_url(database_url: str) -> Path | None:
    """Return the filesystem path for a ``sqlite:///...`` URL, else None."""
    if not database_url.startswith("sqlite:///"):
        return None
    raw = database_url.removeprefix("sqlite:///")
    if raw in {":memory:", ""}:
        return None
    # Absolute sqlite URL: sqlite:////abs/path
    if database_url.startswith("sqlite:////"):
        return Path("/" + database_url.removeprefix("sqlite:////"))
    path = Path(raw)
    return path if path.is_absolute() else (PROJECT_ROOT / path)


def get_engine(database_url: str | None = None):
    """Create a SQLAlchemy engine for ``database_url`` (defaults to settings).

    - **SQLite**: creates parent dirs for file DBs, ``check_same_thread=False``,
      and enables ``PRAGMA foreign_keys=ON``.
    - **Postgres / Supabase**: normal engine using the ``psycopg`` driver
      (``postgresql+psycopg://...``). No SQLite-specific options.
    """
    raw_url = database_url or get_settings().database_url
    url = normalize_database_url(raw_url)

    if is_sqlite_url(url):
        path = _sqlite_path_from_url(url)
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
        engine = create_engine(
            url,
            echo=False,
            connect_args={"check_same_thread": False},
        )

        @event.listens_for(engine, "connect")
        def _set_sqlite_pragma(dbapi_connection, _connection_record) -> None:  # noqa: ANN001
            cursor = dbapi_connection.cursor()
            cursor.execute("PRAGMA foreign_keys=ON")
            cursor.close()

        return engine

    if is_postgres_url(url):
        # Pool settings suitable for Supabase session / transaction poolers.
        return create_engine(url, echo=False, pool_pre_ping=True)

    # Fallback for any other SQLAlchemy-supported URL.
    return create_engine(url, echo=False)


def get_read_only_engine(database_url: str | None = None):
    """Open an existing SQLite provenance database without permitting writes.

    The immutable SQLite URI prevents table creation, migrations, run creation,
    and evidence persistence. Read-only scientific-agent evaluation is local-only;
    other database backends are rejected instead of receiving weaker guarantees.
    """
    raw_url = database_url or get_settings().database_url
    url = normalize_database_url(raw_url)
    if not is_sqlite_url(url):
        raise ValueError("read-only scientific-agent mode requires SQLite")
    path = _sqlite_path_from_url(url)
    if path is None:
        raise ValueError("read-only scientific-agent mode requires a file database")
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError("read-only provenance database does not exist")
    sqlite_uri = f"{resolved.as_uri()}?mode=ro&immutable=1"
    return create_engine(
        "sqlite://",
        echo=False,
        creator=lambda: sqlite3.connect(
            sqlite_uri,
            uri=True,
            check_same_thread=False,
        ),
    )


def init_db(engine=None) -> None:
    """Create all tables if they do not exist."""
    eng = engine or get_engine()
    SQLModel.metadata.create_all(eng)


@contextmanager
def session_scope(engine=None) -> Iterator[Session]:
    """Yield a session that commits on success and rolls back on error."""
    eng = engine or get_engine()
    with Session(eng) as session:
        try:
            yield session
            session.commit()
        except Exception:
            session.rollback()
            raise


# --------------------------------------------------------------------------------------
# Domain <-> row converters
# --------------------------------------------------------------------------------------
def dossier_run_to_row(run: domain.DossierRun) -> DossierRunRow:
    """Convert a domain :class:`DossierRun` to a table row."""
    return DossierRunRow(**run.model_dump())


def dossier_run_from_row(row: DossierRunRow) -> domain.DossierRun:
    """Convert a table row to a domain :class:`DossierRun`."""
    return domain.DossierRun.model_validate(row.model_dump())


def api_run_to_row(run: domain.ApiRun) -> ApiRunRow:
    """Convert a domain :class:`ApiRun` to a table row."""
    return ApiRunRow(**run.model_dump())


def api_run_from_row(row: ApiRunRow) -> domain.ApiRun:
    """Convert a table row to a domain :class:`ApiRun`."""
    return domain.ApiRun.model_validate(row.model_dump())


def raw_artifact_to_row(artifact: domain.RawArtifact) -> RawArtifactRow:
    """Convert a domain :class:`RawArtifact` to a table row (metadata only)."""
    return RawArtifactRow(**artifact.model_dump())


def raw_artifact_from_row(row: RawArtifactRow) -> domain.RawArtifact:
    """Convert a table row to a domain :class:`RawArtifact`."""
    return domain.RawArtifact.model_validate(row.model_dump())


def evidence_to_row(evidence: domain.EvidenceRecord) -> EvidenceRecordRow:
    """Convert a domain :class:`EvidenceRecord` to a table row (enums -> strings)."""
    data = evidence.model_dump()
    data["source_type"] = evidence.source_type.value
    data["assertion_type"] = evidence.assertion_type.value
    data["evidence_grade"] = evidence.evidence_grade.value
    return EvidenceRecordRow(**data)


def evidence_from_row(row: EvidenceRecordRow) -> domain.EvidenceRecord:
    """Convert a table row to a domain :class:`EvidenceRecord`."""
    return domain.EvidenceRecord.model_validate(row.model_dump())


def claim_to_row(claim: domain.Claim) -> ClaimRow:
    """Convert a domain :class:`Claim` to a table row."""
    data = claim.model_dump()
    if claim.evidence_grade is not None:
        data["evidence_grade"] = claim.evidence_grade.value
    return ClaimRow(**data)


def claim_from_row(row: ClaimRow) -> domain.Claim:
    """Convert a table row to a domain :class:`Claim`."""
    return domain.Claim.model_validate(row.model_dump())


def coverage_row_id(dossier_run_id: str, source_name: str) -> str:
    """Stable primary key for one coverage row per run + source."""
    return f"{dossier_run_id}:{slugify(source_name, allow_underscore=True) or new_id()}"


def coverage_to_row(result: domain.SourceCoverageResult) -> SourceCoverageResultRow:
    """Convert a domain :class:`SourceCoverageResult` to a table row."""
    return SourceCoverageResultRow(
        id=coverage_row_id(result.dossier_run_id, result.source_name),
        dossier_run_id=result.dossier_run_id,
        source_name=result.source_name,
        status=result.status.value,
        raw_artifact_path=result.raw_artifact_path,
        evidence_record_count=result.evidence_record_count,
        error_message=result.error_message,
        report_sections_supported=list(result.report_sections_supported),
        notes=result.notes,
    )


def coverage_from_row(row: SourceCoverageResultRow) -> domain.SourceCoverageResult:
    """Convert a table row to a domain :class:`SourceCoverageResult`."""
    data = row.model_dump()
    data.pop("id", None)
    return domain.SourceCoverageResult.model_validate(data)


# --------------------------------------------------------------------------------------
# Convenience writers / readers
# --------------------------------------------------------------------------------------
def save_dossier_run(session: Session, run: domain.DossierRun) -> domain.DossierRun:
    """Insert or replace a dossier run and return the domain object."""
    session.merge(dossier_run_to_row(run))
    session.flush()
    return run


def get_dossier_run(session: Session, dossier_run_id: str) -> domain.DossierRun | None:
    """Load a dossier run by id, or None if missing."""
    row = session.get(DossierRunRow, dossier_run_id)
    return dossier_run_from_row(row) if row else None


def save_api_run(session: Session, run: domain.ApiRun) -> domain.ApiRun:
    """Insert or replace an API run record."""
    session.merge(api_run_to_row(run))
    session.flush()
    return run


def save_raw_artifact(session: Session, artifact: domain.RawArtifact) -> domain.RawArtifact:
    """Insert or replace a raw artifact *metadata* pointer (not the file bytes)."""
    session.merge(raw_artifact_to_row(artifact))
    session.flush()
    return artifact


def save_evidence_record(
    session: Session, evidence: domain.EvidenceRecord
) -> domain.EvidenceRecord:
    """Insert or replace an evidence record."""
    session.merge(evidence_to_row(evidence))
    session.flush()
    return evidence


def delete_evidence_record(session: Session, evidence_id: str) -> bool:
    """Delete one evidence record by primary key. Returns True if a row was removed.

    Does not delete related RawArtifact or ApiRun rows.
    """
    row = session.get(EvidenceRecordRow, evidence_id)
    if row is None:
        return False
    session.delete(row)
    session.flush()
    return True


def list_evidence_for_run(
    session: Session, dossier_run_id: str
) -> list[domain.EvidenceRecord]:
    """Return all evidence records for a dossier run."""
    rows = session.exec(
        select(EvidenceRecordRow).where(EvidenceRecordRow.dossier_run_id == dossier_run_id)
    ).all()
    return [evidence_from_row(r) for r in rows]


def get_evidence_by_source_id(
    session: Session, source_id: str, dossier_run_id: str | None = None
) -> domain.EvidenceRecord | None:
    """Return the first evidence record matching ``source_id`` (optionally scoped to a run)."""
    stmt = select(EvidenceRecordRow).where(EvidenceRecordRow.source_id == source_id)
    if dossier_run_id is not None:
        stmt = stmt.where(EvidenceRecordRow.dossier_run_id == dossier_run_id)
    row = session.exec(stmt).first()
    return evidence_from_row(row) if row else None


def save_source_coverage(
    session: Session, result: domain.SourceCoverageResult
) -> domain.SourceCoverageResult:
    """Insert or replace a source coverage row for a run."""
    session.merge(coverage_to_row(result))
    session.flush()
    return result


def list_source_coverage(
    session: Session, dossier_run_id: str
) -> list[domain.SourceCoverageResult]:
    """Return all source coverage rows for a dossier run."""
    rows = session.exec(
        select(SourceCoverageResultRow).where(
            SourceCoverageResultRow.dossier_run_id == dossier_run_id
        )
    ).all()
    return [coverage_from_row(r) for r in rows]


def canonical_generated_report_id(dossier_run_id: str) -> str:
    """Return the only valid generated-report id for a dossier run."""
    return f"report-{dossier_run_id}"


def _copy_generated_report(row: GeneratedReportRow) -> GeneratedReportRow:
    return GeneratedReportRow.model_validate(row.model_dump())


def save_generated_report(
    session: Session, report: GeneratedReportRow
) -> GeneratedReportRow:
    """Idempotently upsert generated-report metadata for one dossier run."""
    expected_id = canonical_generated_report_id(report.dossier_run_id)
    if report.id != expected_id:
        raise ValueError(
            f"Generated report id must be {expected_id!r} for dossier run "
            f"{report.dossier_run_id!r}."
        )
    existing = session.exec(
        select(GeneratedReportRow).where(
            GeneratedReportRow.dossier_run_id == report.dossier_run_id
        )
    ).first()
    if existing is not None and existing.id != report.id:
        raise ValueError(
            f"Dossier run {report.dossier_run_id!r} already belongs to "
            f"generated report {existing.id!r}."
        )
    merged = session.merge(report)
    session.flush()
    session.refresh(merged)
    return _copy_generated_report(merged)


def get_generated_report(
    session: Session, report_id: str
) -> GeneratedReportRow | None:
    """Load one generated report by its exact canonical id."""
    row = session.get(GeneratedReportRow, report_id)
    return _copy_generated_report(row) if row else None


def get_generated_report_for_run(
    session: Session, dossier_run_id: str
) -> GeneratedReportRow | None:
    """Load generated-report metadata for one exact dossier run."""
    row = session.exec(
        select(GeneratedReportRow).where(
            GeneratedReportRow.dossier_run_id == dossier_run_id
        )
    ).first()
    return _copy_generated_report(row) if row else None


def list_generated_reports(session: Session) -> list[GeneratedReportRow]:
    """Return generated reports newest first."""
    rows = session.exec(
        select(GeneratedReportRow).order_by(GeneratedReportRow.created_at.desc())
    ).all()
    return [_copy_generated_report(row) for row in rows]


__all__ = [
    "DossierRunRow",
    "ApiRunRow",
    "RawArtifactRow",
    "EvidenceRecordRow",
    "ReportSectionRow",
    "ClaimRow",
    "VerificationResultRow",
    "SourceCoverageResultRow",
    "GeneratedReportRow",
    "is_sqlite_url",
    "is_postgres_url",
    "normalize_database_url",
    "get_engine",
    "get_read_only_engine",
    "init_db",
    "session_scope",
    "dossier_run_to_row",
    "dossier_run_from_row",
    "api_run_to_row",
    "api_run_from_row",
    "raw_artifact_to_row",
    "raw_artifact_from_row",
    "evidence_to_row",
    "evidence_from_row",
    "claim_to_row",
    "claim_from_row",
    "coverage_row_id",
    "coverage_to_row",
    "coverage_from_row",
    "save_dossier_run",
    "get_dossier_run",
    "save_api_run",
    "save_raw_artifact",
    "save_evidence_record",
    "delete_evidence_record",
    "list_evidence_for_run",
    "get_evidence_by_source_id",
    "save_source_coverage",
    "list_source_coverage",
    "canonical_generated_report_id",
    "save_generated_report",
    "get_generated_report",
    "get_generated_report_for_run",
    "list_generated_reports",
]
