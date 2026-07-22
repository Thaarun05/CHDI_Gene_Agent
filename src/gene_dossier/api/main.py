"""FastAPI application for the Gene Dossier Platform.

Exposes a thin HTTP surface over the existing provenance-first pipeline:

- health / version
- source registry inspection
- dossier run start (optional background) via LangGraph workflow
- read-back of runs, evidence, and coverage from the DB
- keyword evidence search (Chroma optional via retrieval; keyword always works)

The LLM is never the source of truth. Chroma is never the source of truth.
Live biomedical API calls happen only when a dossier run is started with
``call_network=True``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Literal

from fastapi import BackgroundTasks, FastAPI, HTTPException, Query
from pydantic import BaseModel, Field

from gene_dossier import __version__
from gene_dossier.config import Settings, get_settings
from gene_dossier.db import (
    get_dossier_run,
    init_db,
    list_evidence_for_run,
    list_source_coverage,
    session_scope,
)
from gene_dossier.models import new_id
from gene_dossier.retrieval import search_evidence_keyword
from gene_dossier.source_registry import get_all_sources, get_source, registry_summary
from gene_dossier.workflow import DossierPassResult, run_gene_dossier_full_api_pass

logger = logging.getLogger(__name__)

app = FastAPI(
    title="Gene Dossier Platform",
    description=(
        "Provenance-first CHDI-style gene dossiers. Every claim cites "
        "``source_id``s that resolve to EvidenceRecords backed by raw artifacts."
    ),
    version=__version__,
)


# --------------------------------------------------------------------------------------
# Request / response models
# --------------------------------------------------------------------------------------
class HealthResponse(BaseModel):
    status: Literal["ok"] = "ok"
    version: str
    database: Literal["sqlite", "postgres", "other"]


class SourceOut(BaseModel):
    name: str
    priority: str
    source_type: str
    client_module: str | None = None
    normalizer_module: str | None = None
    required_keys: list[str] = Field(default_factory=list)
    optional_keys: list[str] = Field(default_factory=list)
    report_sections: list[str] = Field(default_factory=list)
    default_status: str
    notes: str | None = None


class DossierRunRequest(BaseModel):
    """Start a full API pass for one gene."""

    gene_symbol: str = Field(..., min_length=1, examples=["SREBF2"])
    run_id: str | None = Field(
        default=None, description="Optional fixed dossier_run_id"
    )
    sources: list[str] | None = Field(
        default=None,
        description="Optional non-identity source filter (identity sources always run).",
    )
    call_network: bool = Field(
        default=True,
        description="If false, skip live API calls (offline / test mode).",
    )
    allow_llm: bool = Field(
        default=False,
        description="If true, allow LangChain synthesis when LLM keys are present.",
    )
    write_rancho: bool = True
    write_pdf: bool = True
    persist_db: bool = True
    wait: bool = Field(
        default=False,
        description=(
            "If true, run synchronously and return the completed result. "
            "If false (default), accept the run and execute in a background task."
        ),
    )
    output_dir: str | None = None


class DossierRunAccepted(BaseModel):
    dossier_run_id: str
    gene_symbol: str
    status: Literal["accepted", "completed", "failed"]
    message: str
    output_paths: dict[str, str] = Field(default_factory=dict)
    evidence_count: int | None = None
    claim_count: int | None = None
    errors: list[str] = Field(default_factory=list)
    synthesis_notes: list[str] = Field(default_factory=list)


class DossierRunOut(BaseModel):
    id: str
    gene_symbol: str
    official_symbol: str | None = None
    run_type: str
    status: str
    started_at: str | None = None
    completed_at: str | None = None
    notes: str | None = None


class EvidenceSearchRequest(BaseModel):
    query: str = ""
    gene_symbol: str | None = None
    section: str | None = None
    source_name: str | None = None
    evidence_grade: str | None = None
    assertion_type: str | None = None
    limit: int = Field(default=20, ge=1, le=200)


class RetrievalHitOut(BaseModel):
    source_id: str
    score: float
    method: str
    gene_symbol: str
    section: str
    source_name: str
    evidence_grade: str
    assertion_type: str
    display_text: str


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def _db_kind(settings: Settings | None = None) -> Literal["sqlite", "postgres", "other"]:
    url = (settings or get_settings()).database_url or ""
    if url.startswith(("postgresql://", "postgresql+psycopg://")):
        return "postgres"
    if url.startswith("sqlite:"):
        return "sqlite"
    return "other"


def _paths_as_str(result: DossierPassResult) -> dict[str, str]:
    return {key: str(path) for key, path in result.output_paths.items()}


def _execute_dossier_pass(body: DossierRunRequest, dossier_run_id: str) -> DossierPassResult:
    settings = get_settings()
    settings.ensure_dirs()
    out = Path(body.output_dir) if body.output_dir else None
    return run_gene_dossier_full_api_pass(
        body.gene_symbol.strip(),
        settings=settings,
        output_dir=out,
        dossier_run_id=dossier_run_id,
        sources=body.sources,
        call_network=body.call_network,
        force_deterministic=not body.allow_llm,
        write_rancho=body.write_rancho,
        write_pdf=body.write_pdf and body.write_rancho,
        persist_db=body.persist_db,
    )


def _background_dossier_pass(body: DossierRunRequest, dossier_run_id: str) -> None:
    try:
        result = _execute_dossier_pass(body, dossier_run_id)
        logger.info(
            "Background dossier pass finished gene=%s run_id=%s status=%s evidence=%d",
            result.gene_symbol,
            result.dossier_run_id,
            result.status,
            len(result.evidence_records),
        )
    except Exception as exc:  # noqa: BLE001
        logger.exception(
            "Background dossier pass failed run_id=%s: %s", dossier_run_id, exc
        )


# --------------------------------------------------------------------------------------
# Routes
# --------------------------------------------------------------------------------------
@app.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    """Liveness probe; does not expose DATABASE_URL credentials."""
    return HealthResponse(
        status="ok",
        version=__version__,
        database=_db_kind(),
    )


@app.get("/version")
def version() -> dict[str, str]:
    return {"version": __version__}


@app.get("/sources", response_model=list[SourceOut])
def list_sources() -> list[SourceOut]:
    """Return the full registered source map (Priority A/B/C)."""
    return [
        SourceOut(
            name=s.name,
            priority=s.priority.value,
            source_type=s.source_type.value,
            client_module=s.client_module,
            normalizer_module=s.normalizer_module,
            required_keys=list(s.required_keys),
            optional_keys=list(s.optional_keys),
            report_sections=list(s.report_sections),
            default_status=s.default_status.value,
            notes=s.notes,
        )
        for s in get_all_sources()
    ]


@app.get("/sources/summary")
def sources_summary() -> dict[str, Any]:
    return registry_summary()


@app.get("/sources/{name}", response_model=SourceOut)
def get_source_detail(name: str) -> SourceOut:
    src = get_source(name)
    if src is None:
        raise HTTPException(status_code=404, detail=f"Unknown source: {name}")
    return SourceOut(
        name=src.name,
        priority=src.priority.value,
        source_type=src.source_type.value,
        client_module=src.client_module,
        normalizer_module=src.normalizer_module,
        required_keys=list(src.required_keys),
        optional_keys=list(src.optional_keys),
        report_sections=list(src.report_sections),
        default_status=src.default_status.value,
        notes=src.notes,
    )


@app.post("/dossier/runs", response_model=DossierRunAccepted)
def start_dossier_run(
    body: DossierRunRequest,
    background_tasks: BackgroundTasks,
) -> DossierRunAccepted:
    """Start a full dossier API pass for ``gene_symbol``.

    Default behavior accepts the run (``status=accepted``) and executes the
    LangGraph pass in a background task. Set ``wait=true`` for a synchronous
    completion response (useful for offline tests).
    """
    gene = body.gene_symbol.strip()
    if not gene:
        raise HTTPException(status_code=422, detail="gene_symbol is required")

    if not body.wait and not body.persist_db:
        raise HTTPException(
            status_code=422,
            detail=(
                "persist_db=false is only supported with wait=true because "
                "background runs cannot be polled without DB persistence."
            ),
        )

    dossier_run_id = (body.run_id or new_id()).strip()
    if body.persist_db:
        init_db()

    if body.wait:
        result = _execute_dossier_pass(body, dossier_run_id)
        status: Literal["accepted", "completed", "failed"] = (
            "completed" if result.status == "completed" else "failed"
        )
        return DossierRunAccepted(
            dossier_run_id=result.dossier_run_id,
            gene_symbol=result.gene_symbol,
            status=status,
            message=f"Dossier pass finished with status={result.status}",
            output_paths=_paths_as_str(result),
            evidence_count=len(result.evidence_records),
            claim_count=len(result.claims),
            errors=list(result.errors),
            synthesis_notes=list(result.synthesis_notes),
        )

    background_tasks.add_task(_background_dossier_pass, body, dossier_run_id)
    return DossierRunAccepted(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene,
        status="accepted",
        message=(
            "Dossier pass accepted and scheduled in a background task. "
            "Poll GET /dossier/runs/{dossier_run_id} for status."
        ),
    )


@app.get("/dossier/runs/{dossier_run_id}", response_model=DossierRunOut)
def get_run(dossier_run_id: str) -> DossierRunOut:
    """Load a persisted dossier run (requires persist_db during the pass)."""
    init_db()
    with session_scope() as session:
        run = get_dossier_run(session, dossier_run_id)
    if run is None:
        raise HTTPException(status_code=404, detail=f"Run not found: {dossier_run_id}")
    return DossierRunOut(
        id=run.id,
        gene_symbol=run.gene_symbol,
        official_symbol=run.official_symbol,
        run_type=run.run_type,
        status=run.status,
        started_at=run.started_at.isoformat() if run.started_at else None,
        completed_at=run.completed_at.isoformat() if run.completed_at else None,
        notes=run.notes,
    )


@app.get("/dossier/runs/{dossier_run_id}/evidence")
def get_run_evidence(
    dossier_run_id: str,
    limit: int = Query(default=200, ge=1, le=2000),
) -> dict[str, Any]:
    """Return evidence records for a run (provenance DB; not Chroma)."""
    init_db()
    with session_scope() as session:
        run = get_dossier_run(session, dossier_run_id)
        if run is None:
            raise HTTPException(
                status_code=404, detail=f"Run not found: {dossier_run_id}"
            )
        records = list_evidence_for_run(session, dossier_run_id)
    sliced = records[:limit]
    return {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": run.gene_symbol,
        "count": len(records),
        "returned": len(sliced),
        "evidence": [r.model_dump(mode="json") for r in sliced],
    }


@app.get("/dossier/runs/{dossier_run_id}/coverage")
def get_run_coverage(dossier_run_id: str) -> dict[str, Any]:
    """Return source coverage rows for a run."""
    init_db()
    with session_scope() as session:
        run = get_dossier_run(session, dossier_run_id)
        if run is None:
            raise HTTPException(
                status_code=404, detail=f"Run not found: {dossier_run_id}"
            )
        rows = list_source_coverage(session, dossier_run_id)
    return {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": run.gene_symbol,
        "count": len(rows),
        "coverage": [r.model_dump(mode="json") for r in rows],
    }


@app.post(
    "/dossier/runs/{dossier_run_id}/search",
    response_model=list[RetrievalHitOut],
)
def search_run_evidence(
    dossier_run_id: str,
    body: EvidenceSearchRequest,
) -> list[RetrievalHitOut]:
    """Keyword/metadata search over a run's evidence (always available; no Chroma)."""
    init_db()
    with session_scope() as session:
        run = get_dossier_run(session, dossier_run_id)
        if run is None:
            raise HTTPException(
                status_code=404, detail=f"Run not found: {dossier_run_id}"
            )
        records = list_evidence_for_run(session, dossier_run_id)

    hits = search_evidence_keyword(
        records,
        body.query,
        gene_symbol=body.gene_symbol,
        section=body.section,
        source_name=body.source_name,
        evidence_grade=body.evidence_grade,
        assertion_type=body.assertion_type,
        limit=body.limit,
    )
    return [
        RetrievalHitOut(
            source_id=h.source_id,
            score=h.score,
            method=h.method,
            gene_symbol=h.record.gene_symbol,
            section=h.record.section,
            source_name=h.record.source_name,
            evidence_grade=str(
                getattr(h.record.evidence_grade, "value", h.record.evidence_grade)
            ),
            assertion_type=str(
                getattr(h.record.assertion_type, "value", h.record.assertion_type)
            ),
            display_text=h.record.display_text,
        )
        for h in hits
    ]


def create_app() -> FastAPI:
    """Factory for ASGI servers / tests."""
    return app


__all__ = ["app", "create_app"]
