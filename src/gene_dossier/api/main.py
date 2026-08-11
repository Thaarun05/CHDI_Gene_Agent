"""FastAPI application for the Gene Dossier Platform.

Exposes a thin HTTP surface over the existing provenance-first pipeline:

- health / version
- source registry inspection
- dossier run start (optional background) via LangGraph workflow
- read-back of runs, evidence, and coverage from the DB
- real frontend API integration (/api/genes, /api/jobs, /api/reports, /api/ask, /api/compare)
"""

from __future__ import annotations

import logging
import re
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from fastapi import APIRouter, BackgroundTasks, FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse
from pydantic import BaseModel, Field

from gene_dossier import __version__
from gene_dossier.config import PROJECT_ROOT, Settings, get_settings
from gene_dossier.db import (
    DossierRunRow,
    EvidenceRecordRow,
    RawArtifactRow,
    evidence_from_row,
    get_dossier_run,
    init_db,
    list_evidence_for_run,
    list_source_coverage,
    session_scope,
)
from gene_dossier.models import EvidenceRecord, new_id
from gene_dossier.retrieval import (
    ChromaEvidenceIndex,
    RetrievalHit,
    build_local_minilm_embedding_function,
    search_evidence_keyword,
    vector_id_for_record,
)
from gene_dossier.section_bundle import SUPPORTED_SECTION_BUNDLE_KEYS, run_section_bundle
from gene_dossier.source_registry import get_all_sources, get_source, registry_summary
from gene_dossier.workflow import DossierPassResult, run_gene_dossier_full_api_pass
from sqlmodel import select

logger = logging.getLogger(__name__)

app = FastAPI(
    title="Gene Dossier Platform",
    description=(
        "Provenance-first CHDI-style gene dossiers. Every claim cites "
        "``source_id``s that resolve to EvidenceRecords backed by raw artifacts."
    ),
    version=__version__,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
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


class GeneOut(BaseModel):
    symbol: str
    name: str
    organism: str
    entrezGeneId: str
    uniprotAccession: str
    summary: str


class EvidenceCoverageRowOut(BaseModel):
    category: str
    status: str
    detail: str | None = None


EvidenceUniverseName = Literal[
    "accepted_demo",
    "explicit_run",
    "accepted_demo_with_tool_overlay",
]


class EvidenceUniverseOut(BaseModel):
    baseEvidenceRunId: str | None = None
    toolRunIds: list[str] = Field(default_factory=list)
    dossierRunIds: list[str] = Field(default_factory=list)
    evidenceUniverse: EvidenceUniverseName


class EvidenceCoverageResponseOut(EvidenceUniverseOut):
    geneSymbol: str
    rows: list[EvidenceCoverageRowOut]


class EvidenceRecordFrontendOut(BaseModel):
    id: str
    geneSymbol: str
    sourceName: str
    evidenceType: str
    factType: str
    # evidenceClass is only set when a scientific evidence classification
    # (e.g. direct_target_evidence, indirect_pathway_effect) is present in
    # the record's value dict.  EvidenceGrade (A/B/C) is a separate concept.
    evidenceClass: str | None = None
    evidenceGrade: str | None = None
    section: str
    subsection: str | None = None
    sourceIdentifier: str | None = None
    retrievedAt: str
    displayText: str
    status: str | None = "Available"
    apiRunId: str | None = None
    rawArtifactId: str | None = None
    sourceUrl: str | None = None


class EvidenceListResponseOut(EvidenceUniverseOut):
    geneSymbol: str
    records: list[EvidenceRecordFrontendOut]


class WorkflowJobStageOut(BaseModel):
    id: str
    label: str
    status: str


class WorkflowJobOut(BaseModel):
    id: str
    geneSymbol: str
    jobType: str
    status: str
    stages: list[WorkflowJobStageOut]
    createdAt: str
    completedAt: str | None = None
    artifactIds: list[str] | None = None
    dossierRunId: str | None = None
    sectionKeys: list[str] | None = None
    errors: list[str] | None = None


# --------------------------------------------------------------------------------------
# In-memory job store (lives for server lifetime; keyed by frontend job_id)
# --------------------------------------------------------------------------------------
# Structure: {job_id: {"job": WorkflowJobOut, "dossier_run_id": str, "output_paths": dict}}
_JOB_STORE: dict[str, dict[str, Any]] = {}
_JOB_STORE_LOCK = threading.Lock()

# Default HD dossier demo section keys (do NOT modify DEFAULT_SECTION_BUNDLE_KEYS)
_HD_DOSSIER_DEFAULT_SECTIONS: list[str] = [
    "1a", "1b", "1c", "1d", "1e",
    "2a", "2b", "2c",
    "3a",
    "4a",
    "5a", "5b",
    "6a",
    "7a",
]

# Friday demo provenance registry (keyed by gene symbol).
# base_evidence_run_id is the default evidence universe; report_run_id is only
# the accepted rendered dossier artifact.
DEMO_GENE_REGISTRY: dict[str, dict[str, Any]] = {
    "SREBF2": {
        "report_id": "rep-srebf2",
        "html_path": (
            PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "SREBF2_full_1a7a"
            / "407e1a4293c6424e8b6b830a1f0a7c60"
            / "section_1.html"
        ),
        "pdf_path": (
            PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "SREBF2_full_1a7a"
            / "407e1a4293c6424e8b6b830a1f0a7c60"
            / "section_1.pdf"
        ),
        "base_evidence_run_id": "407e1a4293c6424e8b6b830a1f0a7c60",
        "report_run_id": "cb9030ab81dc42db80b81dd15d48e653",
        "name": "Sterol regulatory element binding transcription factor 2",
        "entrez_gene_id": "6721",
        "uniprot_accession": "Q12772",
        "description": (
            "SREBF2 (Sterol Regulatory Element Binding Transcription Factor 2) is a basic "
            "helix-loop-helix-leucine zipper (bHLH-Zip) transcription factor that serves as "
            "the master regulator of cellular cholesterol homeostasis."
        ),
    },
    "CDH10": {
        "report_id": "rep-cdh10",
        "html_path": (
            PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "CDH10_full_1a7a"
            / "d94f392f4a3941d5a59f697f58d18234"
            / "section_1.html"
        ),
        "pdf_path": (
            PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "CDH10_full_1a7a"
            / "d94f392f4a3941d5a59f697f58d18234"
            / "section_1.pdf"
        ),
        "base_evidence_run_id": "d94f392f4a3941d5a59f697f58d18234",
        "report_run_id": "ae97cb43e4d94732b72ef86cecc3f40d",
        "name": "Cadherin 10",
        "entrez_gene_id": "1008",
        "uniprot_accession": "Q9Y6N8",
        "description": (
            "CDH10 (Cadherin 10) is a type II classical cadherin belonging to the cadherin "
            "superfamily. It encodes a calcium-dependent cell-cell adhesion glycoprotein that "
            "is predominantly expressed in brain tissue."
        ),
    },
}

_DEMO_REPORT_REGISTRY = DEMO_GENE_REGISTRY


class ReportArtifactOut(BaseModel):
    id: str
    geneSymbol: str
    title: str
    status: str
    createdAt: str
    sections: list[str]
    htmlUrl: str | None = None
    pdfUrl: str | None = None


class ArtifactOut(BaseModel):
    id: str
    path: str
    artifactType: str
    exists: bool


class JobArtifactsResponseOut(BaseModel):
    jobId: str
    dossierRunId: str | None = None
    report: ReportArtifactOut | None = None
    supplementaryArtifacts: list[ArtifactOut] = Field(default_factory=list)


class CitationOut(BaseModel):
    id: str
    label: str
    evidenceRecordId: str
    sourceName: str


class AskRequest(BaseModel):
    question: str
    gene_symbol: str = "SREBF2"
    dossier_run_id: str | None = None
    refresh_if_available: bool = False
    tool_run_ids: list[str] = Field(default_factory=list)


class AskResponseOut(BaseModel):
    status: str
    question: str
    geneSymbol: str
    summary: str
    retrievalMethod: str  # 'semantic' | 'keyword' | 'abstain'
    generationMethod: str  # 'grounded_llm' | 'deterministic' | 'abstain'
    embeddingBackend: str  # 'local_minilm' | 'real' | 'hash_test_fallback' | 'unavailable'
    baseEvidenceRunId: str | None = None
    toolRunIds: list[str] = Field(default_factory=list)
    dossierRunIds: list[str] = Field(default_factory=list)
    evidenceUniverse: EvidenceUniverseName
    evidenceBlocks: list[dict[str, Any]]
    limitations: list[str]
    citations: list[CitationOut]
    evidenceUsedCount: int
    sourcesCount: int
    sourcesUsed: list[str]
    toolsInvokedCount: int
    toolActivity: list[dict[str, Any]]
    agentActivity: list[str]


class CompareRequest(BaseModel):
    genes: list[str] = Field(default_factory=lambda: ["SREBF2", "CDH10"])
    tool_run_ids: dict[str, list[str]] = Field(default_factory=dict)
    dossier_run_ids: dict[str, list[str]] = Field(default_factory=dict)


class ComparisonCellOut(BaseModel):
    status: str
    summary: str
    evidenceCount: int
    evidenceRecordIds: list[str]


class ComparisonResponseOut(BaseModel):
    genes: list[str]
    dimensions: list[str]
    matrix: list[dict[str, Any]]
    narrative: str
    evidenceUniverses: dict[str, EvidenceUniverseOut]


# --------------------------------------------------------------------------------------
# Helper functions
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


def _map_evidence_row(row: EvidenceRecordRow, art_map: dict[str, str]) -> EvidenceRecordFrontendOut:
    source_url = None
    if row.raw_artifact_id:
        source_url = art_map.get(row.raw_artifact_id)
    if not source_url and isinstance(row.value, dict):
        source_url = (
            row.value.get("url")
            or row.value.get("source_url")
            or row.value.get("link")
        )

    # evidenceGrade is the provenance quality tier (A/B/C)
    grade_val = str(getattr(row.evidence_grade, "value", row.evidence_grade)) if row.evidence_grade else None

    # evidenceClass is a *scientific* classification distinct from grade.
    # Only expose it when the record's value dict contains an explicit
    # classification key (e.g. evidence_class, assertion_class).
    ev_class: str | None = None
    if isinstance(row.value, dict):
        ev_class = (
            row.value.get("evidence_class")
            or row.value.get("assertion_class")
            or row.value.get("evidence_classification")
        )

    return EvidenceRecordFrontendOut(
        id=row.id,
        geneSymbol=row.gene_symbol,
        sourceName=row.source_name,
        evidenceType=row.assertion_type or row.source_type,
        factType=row.fact_type,
        evidenceClass=ev_class,
        evidenceGrade=grade_val,
        section=row.section,
        subsection=row.subsection,
        sourceIdentifier=row.source_id,
        retrievedAt=row.created_at.isoformat() if row.created_at else datetime.now().isoformat(),
        displayText=row.display_text,
        status="Available",
        apiRunId=row.api_run_id,
        rawArtifactId=row.raw_artifact_id,
        sourceUrl=source_url,
    )


def _unique_ids(ids: list[str] | tuple[str, ...] | set[str] | None) -> list[str]:
    return [rid for rid in dict.fromkeys(str(rid).strip() for rid in (ids or [])) if rid]


def resolve_evidence_universe(
    gene_symbol: str,
    requested_run_id: str | None = None,
    tool_run_ids: list[str] | None = None,
) -> EvidenceUniverseOut:
    """Resolve the explicit run universe for frontend-facing evidence endpoints."""
    gene = gene_symbol.strip().upper()
    tool_ids = _unique_ids(tool_run_ids)
    if requested_run_id:
        base_id = requested_run_id.strip()
        dossier_run_ids = _unique_ids([base_id, *tool_ids])
        return EvidenceUniverseOut(
            baseEvidenceRunId=base_id,
            toolRunIds=[rid for rid in tool_ids if rid != base_id],
            dossierRunIds=dossier_run_ids,
            evidenceUniverse="explicit_run",
        )

    demo = DEMO_GENE_REGISTRY.get(gene)
    if not demo:
        raise HTTPException(status_code=404, detail=f"No accepted demo evidence baseline for gene: {gene_symbol}")

    base_id = str(demo["base_evidence_run_id"])
    overlay_ids = [rid for rid in tool_ids if rid != base_id]
    return EvidenceUniverseOut(
        baseEvidenceRunId=base_id,
        toolRunIds=overlay_ids,
        dossierRunIds=_unique_ids([base_id, *overlay_ids]),
        evidenceUniverse="accepted_demo_with_tool_overlay" if overlay_ids else "accepted_demo",
    )


def _load_evidence_rows_for_universe(
    gene_symbol: str,
    universe: EvidenceUniverseOut,
    *,
    limit: int | None = None,
    section: str | None = None,
    source: str | None = None,
    assertion_type: str | None = None,
) -> list[EvidenceRecordRow]:
    """Load only records in the selected provenance universe."""
    sym = gene_symbol.strip().upper()
    if not universe.dossierRunIds:
        return []
    with session_scope() as session:
        stmt = select(EvidenceRecordRow).where(
            EvidenceRecordRow.gene_symbol == sym,
            EvidenceRecordRow.dossier_run_id.in_(universe.dossierRunIds),
        )
        if section:
            stmt = stmt.where(EvidenceRecordRow.section.contains(section))
        if source:
            stmt = stmt.where(EvidenceRecordRow.source_name == source)
        if assertion_type:
            stmt = stmt.where(EvidenceRecordRow.assertion_type == assertion_type)
        if limit is not None:
            stmt = stmt.limit(limit)
        return list(session.exec(stmt).all())


def _load_domain_records_for_universe(
    gene_symbol: str,
    universe: EvidenceUniverseOut,
    *,
    limit: int | None = None,
    section: str | None = None,
    source: str | None = None,
    assertion_type: str | None = None,
) -> list[EvidenceRecord]:
    sym = gene_symbol.strip().upper()
    if not universe.dossierRunIds:
        return []
    with session_scope() as session:
        stmt = select(EvidenceRecordRow).where(
            EvidenceRecordRow.gene_symbol == sym,
            EvidenceRecordRow.dossier_run_id.in_(universe.dossierRunIds),
        )
        if section:
            stmt = stmt.where(EvidenceRecordRow.section.contains(section))
        if source:
            stmt = stmt.where(EvidenceRecordRow.source_name == source)
        if assertion_type:
            stmt = stmt.where(EvidenceRecordRow.assertion_type == assertion_type)
        if limit is not None:
            stmt = stmt.limit(limit)
        rows = session.exec(stmt).all()
        return [evidence_from_row(r) for r in rows]


def _load_frontend_evidence_for_universe(
    gene_symbol: str,
    universe: EvidenceUniverseOut,
    *,
    limit: int,
    section: str | None = None,
    source: str | None = None,
    assertion_type: str | None = None,
) -> list[EvidenceRecordFrontendOut]:
    sym = gene_symbol.strip().upper()
    if not universe.dossierRunIds:
        return []
    with session_scope() as session:
        stmt = select(EvidenceRecordRow).where(
            EvidenceRecordRow.gene_symbol == sym,
            EvidenceRecordRow.dossier_run_id.in_(universe.dossierRunIds),
        )
        if section:
            stmt = stmt.where(EvidenceRecordRow.section.contains(section))
        if source:
            stmt = stmt.where(EvidenceRecordRow.source_name == source)
        if assertion_type:
            stmt = stmt.where(EvidenceRecordRow.assertion_type == assertion_type)
        rows = session.exec(stmt.limit(limit * 2)).all()

        art_ids = _unique_ids([r.raw_artifact_id for r in rows if r.raw_artifact_id])
        art_map: dict[str, str] = {}
        if art_ids:
            arts = session.exec(select(RawArtifactRow).where(RawArtifactRow.id.in_(art_ids))).all()
            art_map = {a.id: a.original_url for a in arts if a.original_url}

        seen: set[str] = set()
        out: list[EvidenceRecordFrontendOut] = []
        for r in rows:
            key = f"{r.dossier_run_id}:{r.id}"
            if key in seen:
                continue
            seen.add(key)
            out.append(_map_evidence_row(r, art_map))
            if len(out) >= limit:
                break
        return out


def _artifact_map_for_rows(rows: list[EvidenceRecordRow]) -> dict[str, str]:
    art_ids = _unique_ids([r.raw_artifact_id for r in rows if r.raw_artifact_id])
    if not art_ids:
        return {}
    with session_scope() as session:
        arts = session.exec(select(RawArtifactRow).where(RawArtifactRow.id.in_(art_ids))).all()
        return {a.id: a.original_url for a in arts if a.original_url}


def _row_ref(row: EvidenceRecordRow) -> dict[str, Any]:
    assertion_type = getattr(getattr(row, "assertion_type", None), "value", getattr(row, "assertion_type", ""))
    return {
        "id": getattr(row, "id", ""),
        "section": getattr(row, "section", "") or "",
        "subsection": getattr(row, "subsection", "") or "",
        "source_name": getattr(row, "source_name", "") or "",
        "source_type": getattr(row, "source_type", "") or "",
        "assertion_type": str(assertion_type or ""),
        "fact_type": getattr(row, "fact_type", "") or "",
        "source_id": getattr(row, "source_id", "") or "",
        "display_text": getattr(row, "display_text", "") or "",
        "dossier_run_id": getattr(row, "dossier_run_id", "") or "",
    }


def _record_ref(record: EvidenceRecord) -> dict[str, Any]:
    assertion_type = getattr(record.assertion_type, "value", record.assertion_type)
    source_type = getattr(record.source_type, "value", record.source_type)
    return {
        "id": record.id,
        "section": record.section or "",
        "subsection": record.subsection or "",
        "source_name": record.source_name or "",
        "source_type": str(source_type or ""),
        "assertion_type": str(assertion_type or ""),
        "fact_type": record.fact_type or "",
        "source_id": record.source_id or "",
        "display_text": record.display_text or "",
        "dossier_run_id": record.dossier_run_id or "",
    }


def _text_fields(ref: dict[str, Any]) -> str:
    return " ".join(str(ref.get(k) or "") for k in (
        "section", "subsection", "source_name", "source_type",
        "assertion_type", "fact_type", "source_id", "display_text",
    )).lower()


def _evidence_category(ref: dict[str, Any]) -> str:
    text = _text_fields(ref)
    source = str(ref.get("source_name") or "").lower()
    assertion = str(ref.get("assertion_type") or "").lower()

    if source == "ctd" or "section_6a" in text or "ctd" in text or "chemical perturbation" in text:
        return "chemical_perturbation"
    if (
        assertion == "chemical_tool"
        or "section_7a" in text
        or "chemical tool" in text
        or source in {"chembl", "pubchem", "drugbank", "pubtator3"}
        or "ncats" in source
        or "bioactivity" in text
        or "small molecule" in text
    ):
        return "chemical_tool"
    if assertion == "ppi" or "protein-protein" in text or "ppi" in text or source in {"string", "biogrid"}:
        return "ppi"
    if "harmonizome" in source or "transcription factor" in text or assertion == "transcription_factor_association":
        return "tf"
    if "geo profiles" in source or source == "geo" or "section_3a" in text or "geo perturbation" in text:
        return "geo"
    if assertion in {"expression", "cell_type_expression"} or any(
        token in text for token in ("expression", "gtex", "allen brain atlas", "brainrnaseq", "human brain transcriptome")
    ):
        return "expression"
    if any(token in text for token in ("structure", "domain", "alphafold", "pdbe", "ucsc", "cdd")):
        return "structure"
    if assertion in {"gene_identity", "orthology"} or any(
        token in text for token in ("general gene", "identity", "ncbi gene", "uniprot", "ensembl")
    ):
        return "identity"
    return "other"


_CATEGORY_LABELS: dict[str, str] = {
    "identity": "Gene Identity",
    "structure": "Genomic & Domain Structure",
    "expression": "Expression Profile",
    "geo": "GEO Perturbations",
    "tf": "Transcription Factors",
    "ppi": "Protein Interactions",
    "chemical_perturbation": "Chemical Perturbations",
    "chemical_tool": "Chemical Tools",
}

_CATEGORY_DETAILS: dict[str, tuple[str, str]] = {
    "identity": ("Identity records: {count}", "Identity records not found"),
    "structure": ("Structural/domain records: {count}", "Structural records not found"),
    "expression": ("Expression records: {count}", "Expression records not found"),
    "geo": ("GEO perturbation records: {count}", "No GEO perturbation records"),
    "tf": ("TF association records: {count}", "No TF association records"),
    "ppi": ("PPI records: {count}", "No PPI records"),
    "chemical_perturbation": ("CTD chemical-gene interaction records ({count})", "No CTD perturbation records"),
    "chemical_tool": ("Chemical tool/bioactivity records ({count})", "No chemical tool records"),
}


def _coverage_rows_from_refs(refs: list[dict[str, Any]]) -> list[EvidenceCoverageRowOut]:
    counts: dict[str, int] = {key: 0 for key in _CATEGORY_LABELS}
    for ref in refs:
        category = _evidence_category(ref)
        if category in counts:
            counts[category] += 1

    rows: list[EvidenceCoverageRowOut] = []
    for key, label in _CATEGORY_LABELS.items():
        count = counts[key]
        detail_ok, detail_missing = _CATEGORY_DETAILS[key]
        rows.append(
            EvidenceCoverageRowOut(
                category=label,
                status="Available" if count > 0 else "Not available",
                detail=detail_ok.format(count=count) if count > 0 else detail_missing,
            )
        )
    return rows


def _filter_records_by_category(records: list[EvidenceRecord], category: str | None) -> list[EvidenceRecord]:
    if not category:
        return records
    return [record for record in records if _evidence_category(_record_ref(record)) == category]


def _find_report_files(gene_symbol: str) -> tuple[str | None, Path | None, Path | None]:
    """Resolve report HTML/PDF for a gene using the demo registry first, then filesystem scan."""
    gene = gene_symbol.upper()

    # 1. Check known-good demo registry
    demo = _DEMO_REPORT_REGISTRY.get(gene)
    if demo:
        html_p = demo["html_path"]
        pdf_p = demo["pdf_path"]
        return (
            demo["report_id"],
            html_p if html_p.exists() else None,
            pdf_p if pdf_p.exists() else None,
        )

    # 2. Filesystem scan over data/outputs
    outputs_dir = PROJECT_ROOT / "data" / "outputs"
    if outputs_dir.exists():
        for json_file in sorted(outputs_dir.glob("*_report.json"), reverse=True):
            run_id = json_file.name.split("_")[0]
            html_p = outputs_dir / f"{run_id}_rancho_report.html"
            pdf_p = outputs_dir / f"{run_id}_rancho_report.pdf"
            if html_p.exists():
                try:
                    text = json_file.read_text(encoding="utf-8", errors="replace")
                    if gene in text:
                        return run_id, html_p, pdf_p if pdf_p.exists() else None
                except Exception:
                    pass

    return None, None, None


# --------------------------------------------------------------------------------------
# Shared endpoint handlers
# --------------------------------------------------------------------------------------
def _extract_gene_identity(rows: list[EvidenceRecordRow]) -> dict[str, str]:
    """Parse canonical identity fields from stored EvidenceRecords.

    Returns a dict with keys: name, entrez_gene_id, uniprot_accession, description.
    Only fills keys that are backed by actual stored evidence.
    """
    identity: dict[str, str] = {}
    for r in rows:
        val = r.value if isinstance(r.value, dict) else {}
        # NCBI Gene row
        if r.source_name in {"NCBI Gene", "NCBI Datasets"} and r.fact_type in {
            "entrez_gene_id",
            "ncbi_gene",
            "gene_summary",
        }:
            if "entrez_gene_id" not in identity and val.get("entrez_gene_id"):
                identity["entrez_gene_id"] = str(val["entrez_gene_id"])
            if "name" not in identity and val.get("nomenclaturesymbol"):
                identity["name"] = str(val["nomenclaturesymbol"])
            if "description" not in identity and val.get("description"):
                identity["description"] = str(val["description"])
        # UniProt row
        if r.source_name == "UniProt" and r.fact_type in {
            "uniprot_record",
            "uniprot_accession",
            "protein_annotation",
        }:
            if "uniprot_accession" not in identity and val.get("uniprot_accession"):
                identity["uniprot_accession"] = str(val["uniprot_accession"])
            if "uniprot_accession" not in identity and val.get("primaryAccession"):
                identity["uniprot_accession"] = str(val["primaryAccession"])
            if "protein_name" not in identity and val.get("protein_name"):
                identity["protein_name"] = str(val["protein_name"])
        # Ensembl row
        if r.source_name == "Ensembl" and val.get("description") and "description" not in identity:
            identity["description"] = str(val["description"])
    return identity


def handle_get_gene(symbol: str) -> GeneOut:
    sym = symbol.strip().upper()
    init_db()

    # Known demo fixture fallbacks (used only if canonical DB records are unavailable)
    _DEMO_FIXTURES: dict[str, dict[str, str]] = {
        "SREBF2": {
            "name": "Sterol regulatory element binding transcription factor 2",
            "entrez_gene_id": "6721",
            "uniprot_accession": "Q12772",
            "description": (
                "SREBF2 (Sterol Regulatory Element Binding Transcription Factor 2) is a basic "
                "helix-loop-helix-leucine zipper (bHLH-Zip) transcription factor that serves as "
                "the master regulator of cellular cholesterol homeostasis."
            ),
        },
        "CDH10": {
            "name": "Cadherin 10",
            "entrez_gene_id": "1008",
            "uniprot_accession": "Q9Y6N8",
            "description": (
                "CDH10 (Cadherin 10) is a type II classical cadherin belonging to the cadherin "
                "superfamily. It encodes a calcium-dependent cell-cell adhesion glycoprotein that "
                "is predominantly expressed in brain tissue."
            ),
        },
    }

    with session_scope() as session:
        recs = session.exec(
            select(EvidenceRecordRow).where(EvidenceRecordRow.gene_symbol == sym)
        ).all()
        if not recs:
            if sym not in _DEMO_FIXTURES:
                raise HTTPException(status_code=404, detail=f"Gene not found: {symbol}")
        identity = _extract_gene_identity(list(recs))
        total = len(recs)

    fixture = _DEMO_FIXTURES.get(sym, {})

    # Prefer DB-derived values; fall back to fixtures (log if using fixture)
    entrez = identity.get("entrez_gene_id") or fixture.get("entrez_gene_id") or "Not available"
    uniprot = identity.get("uniprot_accession") or fixture.get("uniprot_accession") or "Not available"
    name = (
        identity.get("protein_name")
        or identity.get("name")
        or fixture.get("name")
        or f"{sym} gene product"
    )
    description = (
        identity.get("description")
        or fixture.get("description")
        or f"Biological evidence dossier for {sym} ({total} evidence records)."
    )

    if fixture and not identity:
        logger.info(
            "handle_get_gene: using fixture fallback for %s (no canonical DB identity records)", sym
        )

    return GeneOut(
        symbol=sym,
        name=name,
        organism="Homo sapiens",
        entrezGeneId=entrez,
        uniprotAccession=uniprot,
        summary=description,
    )


def handle_get_gene_coverage(symbol: str) -> EvidenceCoverageResponseOut:
    """Derive coverage from the accepted baseline evidence universe by default."""
    sym = symbol.strip().upper()
    init_db()
    universe = resolve_evidence_universe(sym)

    # Reduce ORM rows to plain values while the session is open. SQLAlchemy
    # expires instances on commit, so reading row attributes after the session
    # closes can trigger DetachedInstanceError.
    with session_scope() as session:
        stmt = select(EvidenceRecordRow).where(
            EvidenceRecordRow.gene_symbol == sym,
            EvidenceRecordRow.dossier_run_id.in_(universe.dossierRunIds),
        )
        rows = session.exec(stmt).all()
        record_refs = [_row_ref(r) for r in rows]

    return EvidenceCoverageResponseOut(
        geneSymbol=sym,
        **universe.model_dump(),
        rows=_coverage_rows_from_refs(record_refs),
    )


def handle_list_gene_evidence(
    symbol: str,
    limit: int = 200,
    section: str | None = None,
    source: str | None = None,
    type: str | None = None,
) -> EvidenceListResponseOut:
    sym = symbol.strip().upper()
    init_db()
    universe = resolve_evidence_universe(sym)
    records = _load_frontend_evidence_for_universe(
        sym,
        universe,
        limit=limit,
        section=section,
        source=source,
        assertion_type=type,
    )
    return EvidenceListResponseOut(
        geneSymbol=sym,
        **universe.model_dump(),
        records=records,
    )


def handle_get_evidence_record(evidence_id: str) -> EvidenceRecordFrontendOut:
    init_db()
    with session_scope() as session:
        row = session.get(EvidenceRecordRow, evidence_id)
        if not row:
            row = session.exec(
                select(EvidenceRecordRow).where(EvidenceRecordRow.source_id == evidence_id)
            ).first()
        if not row:
            raise HTTPException(status_code=404, detail=f"EvidenceRecord not found: {evidence_id}")

        art_map: dict[str, str] = {}
        if row.raw_artifact_id:
            art = session.get(RawArtifactRow, row.raw_artifact_id)
            if art and art.original_url:
                art_map[art.id] = art.original_url

        return _map_evidence_row(row, art_map)


def _report_timestamp_from_file(path: Path) -> str:
    """Return ISO 8601 UTC timestamp from a file's mtime, or now() if unavailable."""
    try:
        mtime = path.stat().st_mtime
        return datetime.fromtimestamp(mtime, tz=timezone.utc).isoformat()
    except Exception:
        return datetime.now(tz=timezone.utc).isoformat()


def handle_list_reports() -> list[ReportArtifactOut]:
    """Build report list from the demo registry using actual filesystem timestamps."""
    sections = [
        "Section 1: Identity & Structure",
        "Section 2: Expression Profile",
        "Section 3: GEO Perturbations",
        "Section 4: Transcription Factors",
        "Section 5: Protein Interactions",
        "Section 6: Chemical Perturbations",
        "Section 7: Chemical Tools",
    ]

    out: list[ReportArtifactOut] = []
    for gene_sym, demo in _DEMO_REPORT_REGISTRY.items():
        html_p: Path = demo["html_path"]
        pdf_p: Path = demo["pdf_path"]
        report_id = demo["report_id"]
        # Use html file mtime as authoritative timestamp
        created_at = _report_timestamp_from_file(html_p) if html_p.exists() else _report_timestamp_from_file(pdf_p)
        out.append(
            ReportArtifactOut(
                id=report_id,
                geneSymbol=gene_sym,
                title=f"HD Gene Dossier \u2014 {gene_sym}",
                status="Completed" if html_p.exists() else "Unavailable",
                createdAt=created_at,
                sections=sections,
                htmlUrl=f"/api/reports/{report_id}/html" if html_p.exists() else None,
                pdfUrl=f"/api/reports/{report_id}/pdf" if pdf_p.exists() else None,
            )
        )
    return out


def handle_get_report(report_id: str) -> ReportArtifactOut:
    reports = handle_list_reports()
    for r in reports:
        if r.id.lower() == report_id.lower() or r.geneSymbol.lower() in report_id.lower():
            return r

    raise HTTPException(status_code=404, detail=f"Report not found: {report_id}")


def handle_get_report_html(report_id: str):
    rid = report_id.lower()
    if rid in {"rep-srebf2", "srebf2"}:
        run_id, html_p, _ = _find_report_files("SREBF2")
        if html_p and html_p.exists():
            return FileResponse(html_p, media_type="text/html")
    if rid in {"rep-cdh10", "cdh10"}:
        run_id, html_p, _ = _find_report_files("CDH10")
        if html_p and html_p.exists():
            return FileResponse(html_p, media_type="text/html")

    outputs_dir = get_settings().project_root / "data" / "outputs"
    html_p = outputs_dir / f"{report_id}_rancho_report.html"
    if html_p.exists():
        return FileResponse(html_p, media_type="text/html")

    raise HTTPException(status_code=404, detail=f"Report HTML not found: {report_id}")


def handle_get_report_pdf(report_id: str):
    rid = report_id.lower()
    project_root = PROJECT_ROOT
    if rid in {"rep-srebf2", "srebf2"}:
        _, _, pdf_p = _find_report_files("SREBF2")
        if pdf_p and pdf_p.exists():
            return FileResponse(pdf_p, media_type="application/pdf", filename="SREBF2_HD_Dossier.pdf")
    if rid in {"rep-cdh10", "cdh10"}:
        _, _, pdf_p = _find_report_files("CDH10")
        if pdf_p and pdf_p.exists():
            return FileResponse(pdf_p, media_type="application/pdf", filename="CDH10_HD_Dossier.pdf")

    outputs_dir = project_root / "data" / "outputs"
    pdf_p = outputs_dir / f"{report_id}_rancho_report.pdf"
    if pdf_p.exists():
        return FileResponse(pdf_p, media_type="application/pdf", filename=f"{report_id}.pdf")

    raise HTTPException(status_code=404, detail=f"Report PDF not found: {report_id}")


_PHARMACOLOGY_QUERY_TERMS = {
    "bioactivity",
    "chemical",
    "compound",
    "drug",
    "inhibit",
    "inhibitor",
    "manipulate",
    "manipulated",
    "modulate",
    "modulator",
    "pharmacologic",
    "pharmacological",
    "pharmacologically",
    "small molecule",
    "small-molecule",
    "tool",
}

_CHEMICAL_EVIDENCE_SOURCES = {
    "chembl",
    "ctd",
    "drugbank",
    "ncats",
    "pubchem",
    "pubmed",
    "pubtator3",
}

_CHEMICAL_EVIDENCE_TERMS = {
    "6a",
    "7a",
    "bioactivity",
    "chemical",
    "chembl",
    "compound",
    "ctd",
    "drug",
    "inhibitor",
    "pharmacolog",
    "pubchem",
    "section_6a",
    "section_7a",
    "small molecule",
    "small-molecule",
    "tool",
}


def _has_any_text_term(text: str, terms: set[str]) -> bool:
    lower = text.lower()
    return any(term in lower for term in terms)


def _is_pharmacology_question(question: str) -> bool:
    return _has_any_text_term(question, _PHARMACOLOGY_QUERY_TERMS)


def _is_chemical_evidence_record(record: Any) -> bool:
    source_name = (getattr(record, "source_name", "") or "").lower()
    if source_name in _CHEMICAL_EVIDENCE_SOURCES:
        return True

    assertion_type = getattr(getattr(record, "assertion_type", None), "value", None)
    assertion_type = assertion_type or getattr(record, "assertion_type", "")
    haystack = " ".join(
        [
            getattr(record, "section", "") or "",
            getattr(record, "subsection", "") or "",
            source_name,
            getattr(record, "fact_type", "") or "",
            getattr(record, "source_id", "") or "",
            str(assertion_type or ""),
        ]
    )
    return _has_any_text_term(haystack, _CHEMICAL_EVIDENCE_TERMS)


_FRIDAY_RAG_COLLECTION = "friday_demo_minilm_l6_v2_v1"
_SEMANTIC_EMBEDDING_BACKENDS = {"local_minilm", "real"}

_CONTROLLED_TOOL_REGISTRY: dict[str, dict[str, Any]] = {
    "get_identity": {"category": "identity", "sectionKeys": ["1a"], "label": "identity"},
    "get_expression": {"category": "expression", "sectionKeys": ["2a", "2b", "2c"], "label": "expression"},
    "get_geo": {"category": "geo", "sectionKeys": ["3a"], "label": "GEO"},
    "get_tf": {"category": "tf", "sectionKeys": ["4a"], "label": "transcription factors"},
    "get_ppi": {"category": "ppi", "sectionKeys": ["5a", "5b"], "label": "protein interactions"},
    "get_chemical_perturbations": {
        "category": "chemical_perturbation",
        "sectionKeys": ["6a"],
        "label": "chemical perturbations",
    },
    "get_chemical_tools": {"category": "chemical_tool", "sectionKeys": ["7a"], "label": "chemical tools"},
    "get_full_dossier": {
        "category": None,
        "sectionKeys": list(_HD_DOSSIER_DEFAULT_SECTIONS),
        "label": "full dossier",
    },
}


def _infer_required_category(question: str) -> str | None:
    q = question.lower()
    if any(term in q for term in ("interact", "interaction", "ppi", "protein partner", "binds")):
        return "ppi"
    if any(term in q for term in ("ctd", "chemical perturbation", "chemical-gene", "chemical gene")):
        return "chemical_perturbation"
    if _is_pharmacology_question(question):
        return "chemical_tool"
    if any(term in q for term in ("expression", "expressed", "tissue", "cell type", "brain region")):
        return "expression"
    if any(term in q for term in ("geo", "differential", "knockdown", "perturbation")):
        return "geo"
    if any(term in q for term in ("transcription factor", "tf", "regulates")):
        return "tf"
    if any(term in q for term in ("identity", "entrez", "uniprot", "alias")):
        return "identity"
    return None


def _tool_for_category(category: str | None) -> tuple[str, dict[str, Any]] | None:
    if category is None:
        return None
    for name, spec in _CONTROLLED_TOOL_REGISTRY.items():
        if spec["category"] == category:
            return name, spec
    return None


def _category_scope_label(category: str | None) -> str:
    if category == "chemical_tool":
        return "chemical/tool evidence"
    if category == "chemical_perturbation":
        return "6a chemical perturbation evidence"
    if category == "ppi":
        return "PPI evidence"
    if category and category in _CATEGORY_LABELS:
        return _CATEGORY_LABELS[category].lower()
    return "all evidence"


def _category_exists(records: list[EvidenceRecord], category: str | None) -> bool:
    if category is None:
        return bool(records)
    return bool(_filter_records_by_category(records, category))


def _persistent_chroma_index() -> ChromaEvidenceIndex:
    return ChromaEvidenceIndex(
        collection_name=_FRIDAY_RAG_COLLECTION,
        embedding_function=build_local_minilm_embedding_function(),
        ephemeral=False,
        allow_hash_fallback=False,
        allow_external_embedding_provider=False,
    )


def _index_records(records: list[EvidenceRecord]) -> int:
    if not records:
        return 0
    try:
        index = _persistent_chroma_index()
        return index.upsert_evidence(records) if index.available else 0
    except Exception as exc:  # noqa: BLE001
        logger.warning("ask: persistent Chroma indexing failed: %s", exc)
        return 0


def _retrieve_grounded_hits(
    *,
    question: str,
    gene: str,
    universe: EvidenceUniverseOut,
    records: list[EvidenceRecord],
    category: str | None,
    limit: int,
) -> tuple[list[RetrievalHit], str, list[str], str]:
    retrieval_records = _filter_records_by_category(records, category)
    if category is None:
        retrieval_records = records
    activities = [
        f"Retrieval scope: {_category_scope_label(category)}",
        f"Evidence universe run IDs: {', '.join(universe.dossierRunIds)}",
    ]

    hits: list[RetrievalHit] = []
    retrieval_method = "abstain"
    embedding_backend = "unavailable"
    if retrieval_records:
        try:
            index = _persistent_chroma_index()
            embedding_backend = getattr(index.status, "embedding_backend", "unavailable")
            if index.available and embedding_backend in _SEMANTIC_EMBEDDING_BACKENDS:
                indexed = index.upsert_evidence(records)
                if indexed > 0:
                    lookup = {vector_id_for_record(r): r for r in retrieval_records}
                    semantic_hits = index.query(
                        question,
                        gene_symbol=gene,
                        dossier_run_ids=universe.dossierRunIds,
                        limit=max(50, limit * 4),
                        record_lookup=lookup,
                    )
                    if semantic_hits:
                        hits = semantic_hits[:limit]
                        retrieval_method = "semantic"
                    activities.append(f"Semantic retrieval attempted first; indexed/upserted {indexed} records")
                else:
                    embedding_backend = "unavailable"
                    activities.append("Semantic retrieval unavailable; real embedding upsert failed")
            else:
                if embedding_backend == "hash_test_fallback":
                    activities.append("Semantic retrieval skipped; hash embeddings are test fallback only")
                else:
                    activities.append("Semantic retrieval unavailable; real embedding provider is not available")
        except Exception as exc:  # noqa: BLE001
            embedding_backend = "unavailable"
            activities.append(f"Semantic retrieval failed; falling back to keyword ({type(exc).__name__})")
            logger.warning("ask: semantic retrieval failed: %s", exc)

    if len(hits) < 2 and retrieval_records:
        keyword_hits = search_evidence_keyword(
            retrieval_records,
            question,
            gene_symbol=gene,
            limit=limit,
        )
        if not keyword_hits and category:
            category_terms = {
                "ppi": "protein interaction ppi string biogrid interactor partner",
                "chemical_tool": "chemical tool bioactivity compound drug inhibitor chembl pubchem",
                "chemical_perturbation": "chemical perturbation ctd chemical-gene interaction",
                "expression": "expression tissue gtex allen brain region",
                "geo": "geo perturbation differential expression",
                "tf": "transcription factor association harmonizome",
                "identity": "gene identity ncbi uniprot ensembl",
            }.get(category, "")
            keyword_hits = search_evidence_keyword(
                retrieval_records,
                f"{question} {category_terms}",
                gene_symbol=gene,
                limit=limit,
            )
        merged: dict[str, RetrievalHit] = {h.source_id: h for h in hits}
        for h in keyword_hits:
            vector_id = vector_id_for_record(h.record)
            if vector_id in merged:
                continue
            merged[vector_id] = RetrievalHit(
                record=h.record,
                score=h.score,
                method="keyword" if retrieval_method == "abstain" else h.method,
                source_id=vector_id,
            )
        hits = sorted(merged.values(), key=lambda h: (-h.score, h.source_id))[:limit]
        if hits and retrieval_method == "abstain":
            retrieval_method = "keyword"
            activities.append("Keyword retrieval used as fallback")
        elif keyword_hits:
            activities.append("Keyword retrieval augmented sub-threshold semantic hits")

    return hits, retrieval_method, activities, embedding_backend


def _has_sufficient_evidence(
    *,
    category: str | None,
    records: list[EvidenceRecord],
    hits: list[RetrievalHit],
    min_hits: int,
) -> bool:
    if not _category_exists(records, category):
        return False
    if category is None:
        return len(hits) >= min_hits
    relevant_hits = [
        h for h in hits if _evidence_category(_record_ref(h.record)) == category
    ]
    return len(relevant_hits) >= min_hits


def _include_tool_overlay_hits(
    *,
    hits: list[RetrievalHit],
    records: list[EvidenceRecord],
    category: str | None,
    tool_run_ids: list[str],
    limit: int,
) -> list[RetrievalHit]:
    overlay_ids = set(tool_run_ids)
    if not overlay_ids:
        return hits
    existing_record_ids = {h.record.id for h in hits}
    overlay_records = [
        record
        for record in _filter_records_by_category(records, category)
        if record.dossier_run_id in overlay_ids and record.id not in existing_record_ids
    ]
    overlay_hits = [
        RetrievalHit(
            record=record,
            score=1.0,
            method="tool_overlay",
            source_id=vector_id_for_record(record),
        )
        for record in overlay_records
    ]
    return [*overlay_hits, *hits][:limit] if overlay_hits else hits


def _execute_controlled_tool(gene: str, tool_name: str, spec: dict[str, Any]) -> dict[str, Any]:
    section_keys = list(spec["sectionKeys"])
    activity: dict[str, Any] = {
        "toolName": tool_name,
        "sectionKeys": section_keys,
        "status": "started",
        "dossierRunId": None,
        "evidenceRecordsPersisted": 0,
    }
    try:
        result = run_section_bundle(
            gene,
            section_keys=section_keys,
            settings=get_settings(),
            persist_db=True,
            write_pdf=False,
        )
        activity["status"] = result.status
        activity["dossierRunId"] = result.dossier_run_id
        activity["outputPaths"] = {k: str(v) for k, v in result.output_paths.items()}
        activity["errors"] = list(result.errors)

        run_universe = EvidenceUniverseOut(
            baseEvidenceRunId=result.dossier_run_id,
            toolRunIds=[],
            dossierRunIds=[result.dossier_run_id],
            evidenceUniverse="explicit_run",
        )
        tool_records = _load_domain_records_for_universe(gene, run_universe)
        activity["evidenceRecordsPersisted"] = len(tool_records)
        activity["indexedRecords"] = _index_records(tool_records)
        return activity
    except Exception as exc:  # noqa: BLE001
        logger.exception("Controlled tool execution failed: gene=%s tool=%s", gene, tool_name)
        activity["status"] = "failed"
        activity["error"] = str(exc)
        return activity


def _deterministic_grounded_summary(
    *,
    gene: str,
    hits: list[RetrievalHit],
    retrieval_method: str,
    category: str | None,
) -> str:
    source_set = sorted({h.record.source_name for h in hits if h.record.source_name})
    snippets = [h.record.display_text for h in hits[:3] if h.record.display_text]
    lines = [
        f"Retrieved {len(hits)} {gene} evidence records via {retrieval_method} search within {_category_scope_label(category)}.",
    ]
    if source_set:
        lines.append(f"Sources represented: {', '.join(source_set)}.")
    if snippets:
        lines.append("Representative retrieved evidence:")
        lines.extend(f"[{i}] {snippet[:260]}" for i, snippet in enumerate(snippets, start=1))
    return " ".join(lines)


_HEX_EVIDENCE_ID_RE = re.compile(r"\b[a-f0-9]{32}\b", re.IGNORECASE)
_EVIDENCE_RECORD_LABEL_RE = re.compile(
    r"\bEvidenceRecord(?:\s+ID)?[:#\s]+([A-Za-z0-9_.:-]+)",
    re.IGNORECASE,
)


def _extract_cited_evidence_ids(text: str) -> list[str]:
    """Extract literal EvidenceRecord identifiers cited by generated text."""
    candidates: list[str] = []
    candidates.extend(match.group(1).strip(".,;()[]{}") for match in _EVIDENCE_RECORD_LABEL_RE.finditer(text or ""))
    candidates.extend(match.group(0) for match in _HEX_EVIDENCE_ID_RE.finditer(text or ""))
    return _unique_ids(candidates)


def _llm_citations_are_supported(text: str, allowed_ids: set[str]) -> bool:
    cited_ids = _extract_cited_evidence_ids(text)
    if not cited_ids:
        return False
    return all(cited_id in allowed_ids for cited_id in cited_ids)


def _try_grounded_llm_summary(
    *,
    question: str,
    gene: str,
    hits: list[RetrievalHit],
) -> tuple[str | None, str]:
    settings = get_settings()
    if not settings.has_llm() or not hits:
        return None, "deterministic"
    try:
        from gene_dossier.synthesis import build_chat_model_candidates
    except Exception as exc:  # noqa: BLE001
        logger.warning("ask: LangChain synthesis unavailable: %s", exc)
        return None, "deterministic"

    evidence_lines = []
    allowed_ids = {h.record.id for h in hits}
    for h in hits[:8]:
        evidence_lines.append(
            f"EvidenceRecord {h.record.id} | source={h.record.source_name} | text={h.record.display_text}"
        )
    prompt = (
        "Answer only from the supplied EvidenceRecords. Do not add outside biology. "
        "Cite EvidenceRecord IDs literally. State limitations when evidence is thin. "
        f"Gene: {gene}\nQuestion: {question}\n\n" + "\n".join(evidence_lines)
    )
    for candidate in build_chat_model_candidates(settings):
        try:
            response = candidate.model.invoke(prompt)
            text = getattr(response, "content", response)
            if isinstance(text, list):
                text = " ".join(str(part) for part in text)
            text = str(text).strip()
            if text and _llm_citations_are_supported(text, allowed_ids):
                return text, "grounded_llm"
        except Exception as exc:  # noqa: BLE001
            logger.warning("ask: grounded LLM synthesis failed via %s: %s", candidate.provider, exc)
    return None, "deterministic"


def handle_ask_question(body: AskRequest) -> AskResponseOut:
    """Answer a question grounded strictly in the selected EvidenceRecord universe."""
    MIN_EVIDENCE_FOR_ANSWER = 2
    RETRIEVE_LIMIT = 15
    MAX_BLOCKS = 10

    gene = (body.gene_symbol or "SREBF2").strip().upper()
    q = body.question.strip()
    init_db()

    universe = resolve_evidence_universe(
        gene,
        requested_run_id=body.dossier_run_id,
        tool_run_ids=body.tool_run_ids,
    )
    records = _load_domain_records_for_universe(gene, universe)
    category = _infer_required_category(q)
    tool_activity: list[dict[str, Any]] = []
    agent_activity = [
        f"Resolved target gene symbol: {gene}",
        f"Planner category: {category or 'general'}",
        f"Initial evidence universe: {universe.evidenceUniverse}",
    ]

    hits, retrieval_method, retrieval_activity, embedding_backend = _retrieve_grounded_hits(
        question=q,
        gene=gene,
        universe=universe,
        records=records,
        category=category,
        limit=RETRIEVE_LIMIT,
    )
    hits = _include_tool_overlay_hits(
        hits=hits,
        records=records,
        category=category,
        tool_run_ids=universe.toolRunIds,
        limit=RETRIEVE_LIMIT,
    )
    agent_activity.extend(retrieval_activity)

    sufficient = _has_sufficient_evidence(
        category=category,
        records=records,
        hits=hits,
        min_hits=MIN_EVIDENCE_FOR_ANSWER,
    )
    planned_tool = _tool_for_category(category)
    should_run_tool = bool(planned_tool and (body.refresh_if_available or not sufficient))
    if should_run_tool and planned_tool:
        tool_name, spec = planned_tool
        agent_activity.append(
            f"Selected controlled tool: {tool_name} ({','.join(spec['sectionKeys'])})"
        )
        activity = _execute_controlled_tool(gene, tool_name, spec)
        tool_activity.append(activity)
        new_run_id = activity.get("dossierRunId")
        if new_run_id and activity.get("status") != "failed":
            universe = resolve_evidence_universe(
                gene,
                requested_run_id=body.dossier_run_id,
                tool_run_ids=[*body.tool_run_ids, str(new_run_id)],
            )
            records = _load_domain_records_for_universe(gene, universe)
            hits, retrieval_method, retrieval_activity, embedding_backend = _retrieve_grounded_hits(
                question=q,
                gene=gene,
                universe=universe,
                records=records,
                category=category,
                limit=RETRIEVE_LIMIT,
            )
            hits = _include_tool_overlay_hits(
                hits=hits,
                records=records,
                category=category,
                tool_run_ids=universe.toolRunIds,
                limit=RETRIEVE_LIMIT,
            )
            agent_activity.extend(retrieval_activity)
            sufficient = _has_sufficient_evidence(
                category=category,
                records=records,
                hits=hits,
                min_hits=MIN_EVIDENCE_FOR_ANSWER,
            )
    elif body.refresh_if_available and not planned_tool:
        agent_activity.append("Refresh requested but no allowlisted tool matches the question category")

    if not sufficient:
        return AskResponseOut(
            status="abstain",
            question=q,
            geneSymbol=gene,
            summary=(
                f"Insufficient evidence to answer: the selected universe lacks the required "
                f"{_category_scope_label(category)} with at least {MIN_EVIDENCE_FOR_ANSWER} relevant retrieval hits."
            ),
            retrievalMethod="abstain",
            generationMethod="abstain",
            embeddingBackend=embedding_backend,
            **universe.model_dump(),
            evidenceBlocks=[],
            limitations=[
                "No answer was generated beyond the selected provenance-backed evidence universe.",
                "A matching controlled deterministic workflow may be required if no tool was already invoked.",
            ],
            citations=[],
            evidenceUsedCount=0,
            sourcesCount=0,
            sourcesUsed=[],
            toolsInvokedCount=len(tool_activity),
            toolActivity=tool_activity,
            agentActivity=[
                *agent_activity,
                "Abstained: category-aware sufficiency was not met",
            ],
        )

    # ---- Build citations and evidence blocks from retrieved records only ----
    citations: list[CitationOut] = []
    items: list[dict[str, Any]] = []

    for idx, hit in enumerate(hits[:MAX_BLOCKS]):
        cid = f"cit-{idx + 1}"
        rec = hit.record
        citations.append(
            CitationOut(
                id=cid,
                label=f"[{idx + 1}]",
                evidenceRecordId=rec.id,
                sourceName=rec.source_name,
            )
        )
        items.append({
            "text": rec.display_text,
            "citationIds": [cid],
        })

    source_set = {h.record.source_name for h in hits if h.record.source_name}
    llm_summary, generation_method = _try_grounded_llm_summary(
        question=q,
        gene=gene,
        hits=hits[:MAX_BLOCKS],
    )
    if llm_summary:
        summary = llm_summary
    else:
        summary = _deterministic_grounded_summary(
            gene=gene,
            hits=hits[:MAX_BLOCKS],
            retrieval_method=retrieval_method,
            category=category,
        )

    return AskResponseOut(
        status="answered",
        question=q,
        geneSymbol=gene,
        summary=summary,
        retrievalMethod=retrieval_method,
        generationMethod=generation_method,
        embeddingBackend=embedding_backend,
        **universe.model_dump(),
        evidenceBlocks=[{
            "sourceGroup": f"{gene} Evidence Records",
            "items": items,
        }],
        limitations=[
            "Evidence is derived strictly from normalized provenance records in the dossier database.",
            "Retrieval is scoped to the explicit evidence universe reported with this response.",
            "Live biological validation recommended before clinical translation.",
            f"Retrieval method: {retrieval_method}.",
        ],
        citations=citations,
        evidenceUsedCount=len(items),
        sourcesCount=len(source_set),
        sourcesUsed=sorted(source_set),
        toolsInvokedCount=len(tool_activity),
        toolActivity=tool_activity,
        agentActivity=[
            *agent_activity,
            f"Retrieved {len(hits)} hits via {retrieval_method} retrieval",
            f"Generation method: {generation_method}",
        ],
    )


def handle_compare_genes(body: CompareRequest) -> ComparisonResponseOut:
    """Compare genes using only each gene's selected evidence universe."""
    genes = [g.strip().upper() for g in (body.genes or ["SREBF2", "CDH10"])]
    dimensions = [
        "Gene Identity",
        "Expression",
        "GEO Perturbations",
        "Protein Interactions",
        "Chemical Perturbations",
        "Chemical Tools",
    ]

    init_db()
    evidence_universes: dict[str, EvidenceUniverseOut] = {}
    with session_scope() as session:
        gene_recs: dict[str, list[dict[str, Any]]] = {}
        for g in genes:
            explicit_ids = _unique_ids((body.dossier_run_ids or {}).get(g, []))
            requested_base = explicit_ids[0] if explicit_ids else None
            explicit_tool_ids = [*explicit_ids[1:], *((body.tool_run_ids or {}).get(g, []))]
            universe = resolve_evidence_universe(
                g,
                requested_run_id=requested_base,
                tool_run_ids=explicit_tool_ids,
            )
            evidence_universes[g] = universe
            rows = session.exec(
                select(EvidenceRecordRow).where(
                    EvidenceRecordRow.gene_symbol == g,
                    EvidenceRecordRow.dossier_run_id.in_(universe.dossierRunIds),
                )
            ).all()
            gene_recs[g] = [_row_ref(r) for r in rows]

    def _cell(matched: list[dict[str, Any]], detail_available: str, detail_unavailable: str) -> ComparisonCellOut:
        count = len(matched)
        ids = [r["id"] for r in matched[:5]]  # up to 5 real IDs
        return ComparisonCellOut(
            status="Available" if count > 0 else "Not available",
            summary=detail_available.format(count=count) if count > 0 else detail_unavailable,
            evidenceCount=count,
            evidenceRecordIds=ids,
        )

    matrix: list[dict[str, Any]] = []
    # Track per-gene dimension statuses for narrative
    gene_dim_counts: dict[str, dict[str, int]] = {g: {} for g in genes}

    for dim in dimensions:
        cells: dict[str, ComparisonCellOut] = {}
        for g in genes:
            if dim == "Gene Identity":
                matched = [r for r in gene_recs.get(g, []) if _evidence_category(r) == "identity"]
                cells[g] = _cell(matched,
                    "Gene identity records: {count}",
                    "No identity records",
                )
            elif dim == "Expression":
                matched = [r for r in gene_recs.get(g, []) if _evidence_category(r) == "expression"]
                cells[g] = _cell(matched,
                    "Expression records: {count}",
                    "No expression records",
                )
            elif dim == "GEO Perturbations":
                matched = [r for r in gene_recs.get(g, []) if _evidence_category(r) == "geo"]
                cells[g] = _cell(matched,
                    "GEO perturbation records: {count}",
                    "No GEO perturbation records",
                )
            elif dim == "Protein Interactions":
                matched = [r for r in gene_recs.get(g, []) if _evidence_category(r) == "ppi"]
                cells[g] = _cell(matched,
                    "PPI records: {count}",
                    "No PPI records",
                )
            elif dim == "Chemical Perturbations":
                matched = [r for r in gene_recs.get(g, []) if _evidence_category(r) == "chemical_perturbation"]
                cells[g] = _cell(matched,
                    "CTD chemical-gene interaction records: {count}",
                    "No CTD perturbation records",
                )
            elif dim == "Chemical Tools":
                matched = [r for r in gene_recs.get(g, []) if _evidence_category(r) == "chemical_tool"]
                cells[g] = _cell(matched,
                    "Chemical tool records: {count}",
                    "No chemical tool records",
                )
            else:
                matched = []
                cells[g] = ComparisonCellOut(
                    status="Not available", summary="Dimension not mapped",
                    evidenceCount=0, evidenceRecordIds=[],
                )

            gene_dim_counts[g][dim] = cells[g].evidenceCount

        matrix.append({"dimension": dim, "cells": cells})

    # ---- Deterministic narrative from actual matrix data ----
    narrative_parts: list[str] = []
    for g in genes:
        counts = gene_dim_counts[g]
        available = [dim for dim, cnt in counts.items() if cnt > 0]
        unavailable = [dim for dim, cnt in counts.items() if cnt == 0]
        total = sum(counts.values())
        line = f"{g}: {total} evidence records across {len(available)} dimension(s)"
        if available:
            line += f" ({', '.join(available)})"
        if unavailable:
            line += f"; no records for: {', '.join(unavailable)}"
        narrative_parts.append(line)
    narrative = "  |  ".join(narrative_parts)

    return ComparisonResponseOut(
        genes=genes,
        dimensions=dimensions,
        matrix=matrix,
        narrative=narrative,
        evidenceUniverses=evidence_universes,
    )


# --------------------------------------------------------------------------------------
# API Router for /api (and mounted directly at root for proxy compatibility)
# --------------------------------------------------------------------------------------
api_router = APIRouter(prefix="/api", tags=["api"])


@api_router.get("/genes/{symbol}", response_model=GeneOut)
@app.get("/genes/{symbol}", response_model=GeneOut)
def get_gene_endpoint(symbol: str) -> GeneOut:
    return handle_get_gene(symbol)


@api_router.get("/genes/{symbol}/coverage", response_model=EvidenceCoverageResponseOut)
@app.get("/genes/{symbol}/coverage", response_model=EvidenceCoverageResponseOut)
def get_gene_coverage_endpoint(symbol: str) -> EvidenceCoverageResponseOut:
    return handle_get_gene_coverage(symbol)


@api_router.get("/genes/{symbol}/evidence", response_model=EvidenceListResponseOut)
@app.get("/genes/{symbol}/evidence", response_model=EvidenceListResponseOut)
def list_gene_evidence_endpoint(
    symbol: str,
    limit: int = Query(default=200, ge=1, le=2000),
    section: str | None = None,
    source: str | None = None,
    type: str | None = None,
) -> EvidenceListResponseOut:
    return handle_list_gene_evidence(symbol, limit=limit, section=section, source=source, type=type)


@api_router.get("/evidence/{evidence_record_id}", response_model=EvidenceRecordFrontendOut)
@app.get("/evidence/{evidence_record_id}", response_model=EvidenceRecordFrontendOut)
def get_evidence_record_endpoint(evidence_record_id: str) -> EvidenceRecordFrontendOut:
    return handle_get_evidence_record(evidence_record_id)


@api_router.get("/evidence", response_model=EvidenceListResponseOut)
@app.get("/evidence", response_model=EvidenceListResponseOut)
def list_all_evidence_endpoint(
    gene: str | None = None,
    source: str | None = None,
    type: str | None = None,
    section: str | None = None,
    limit: int = Query(default=200, ge=1, le=2000),
) -> EvidenceListResponseOut:
    symbol = gene or "SREBF2"
    return handle_list_gene_evidence(symbol, limit=limit, section=section, source=source, type=type)


def _run_section_bundle_job(job_id: str, gene: str, section_keys: list[str]) -> None:
    """Execute run_section_bundle in a background thread and update _JOB_STORE."""
    try:
        with _JOB_STORE_LOCK:
            entry = _JOB_STORE.get(job_id)
            if entry:
                entry["job"].status = "Running"

        result = run_section_bundle(
            gene,
            section_keys=section_keys,
            settings=get_settings(),
            persist_db=True,
            write_pdf=False,  # PDF generation is slow; skip for background jobs
        )

        final_status = result.status  # 'completed' | 'partial' | 'failed'
        artifact_ids: list[str] = list(result.output_paths.keys())
        output_paths_str = {k: str(v) for k, v in result.output_paths.items()}

        # Build stages from section_keys
        stages = [
            WorkflowJobStageOut(id=f"s-{k}", label=f"Section {k}", status="Complete")
            for k in section_keys
        ]

        completed_at = datetime.now(tz=timezone.utc).isoformat()
        with _JOB_STORE_LOCK:
            if job_id in _JOB_STORE:
                job_out = _JOB_STORE[job_id]["job"]
                job_out.status = final_status.capitalize() if final_status else "Completed"
                job_out.stages = stages
                job_out.completedAt = completed_at
                job_out.artifactIds = artifact_ids
                job_out.dossierRunId = result.dossier_run_id
                job_out.errors = list(result.errors) if result.errors else []
                _JOB_STORE[job_id]["dossier_run_id"] = result.dossier_run_id
                _JOB_STORE[job_id]["output_paths"] = output_paths_str

        logger.info(
            "job %s finished gene=%s status=%s run_id=%s",
            job_id, gene, final_status, result.dossier_run_id,
        )

    except Exception as exc:
        logger.exception("job %s failed: %s", job_id, exc)
        with _JOB_STORE_LOCK:
            if job_id in _JOB_STORE:
                _JOB_STORE[job_id]["job"].status = "Failed"
                _JOB_STORE[job_id]["job"].errors = [str(exc)]


@api_router.post("/jobs", response_model=WorkflowJobOut)
@app.post("/jobs", response_model=WorkflowJobOut)
def create_job_endpoint(body: dict[str, Any]) -> WorkflowJobOut:
    """Create and enqueue a real section bundle job.

    Body fields:
      gene_symbol: str           (required)
      sections: list[str] | None (optional; defaults to HD dossier demo set)
      use_existing_accepted: bool (if true, return known-good demo run without calling live APIs)
    """
    gene_raw = body.get("gene_symbol") or body.get("gene") or ""
    gene = str(gene_raw).strip().upper()
    if not gene:
        raise HTTPException(status_code=422, detail="gene_symbol is required")

    use_existing = bool(body.get("use_existing_accepted", False))
    requested_sections = body.get("sections") or body.get("section_keys")

    section_keys: list[str]
    if requested_sections:
        # Validate against supported keys
        bad = [k for k in requested_sections if k not in SUPPORTED_SECTION_BUNDLE_KEYS]
        if bad:
            raise HTTPException(
                status_code=422,
                detail=f"Unsupported section keys: {bad}. Supported: {sorted(SUPPORTED_SECTION_BUNDLE_KEYS)}",
            )
        section_keys = [str(k) for k in requested_sections]
    else:
        section_keys = list(_HD_DOSSIER_DEFAULT_SECTIONS)

    now_str = datetime.now(tz=timezone.utc).isoformat()
    job_id = f"job-{gene.lower()}-{new_id()[:8]}"

    # ---- use_existing_accepted: return known-good demo run immediately ----
    if use_existing:
        demo = _DEMO_REPORT_REGISTRY.get(gene)
        if demo:
            stages = [
                WorkflowJobStageOut(id=f"s-{k}", label=f"Section {k}", status="Complete")
                for k in section_keys
            ]
            job_out = WorkflowJobOut(
                id=job_id,
                geneSymbol=gene,
                jobType="hd_dossier",
                status="Completed",
                stages=stages,
                createdAt=now_str,
                completedAt=now_str,
                artifactIds=[demo["report_id"]],
                dossierRunId=demo["report_run_id"],
                sectionKeys=section_keys,
                errors=[],
            )
            with _JOB_STORE_LOCK:
                _JOB_STORE[job_id] = {
                    "job": job_out,
                    "dossier_run_id": demo["report_run_id"],
                    "output_paths": {
                        "html": str(demo["html_path"]),
                        "pdf": str(demo["pdf_path"]),
                    },
                }
            return job_out
        else:
            raise HTTPException(
                status_code=404,
                detail=f"No accepted demo run registered for gene: {gene}",
            )

    # ---- Real execution: create queued job, launch background thread ----
    init_db()
    pending_stages = [
        WorkflowJobStageOut(id=f"s-{k}", label=f"Section {k}", status="Queued")
        for k in section_keys
    ]
    job_out = WorkflowJobOut(
        id=job_id,
        geneSymbol=gene,
        jobType="hd_dossier",
        status="Queued",
        stages=pending_stages,
        createdAt=now_str,
        completedAt=None,
        artifactIds=None,
        dossierRunId=None,
        sectionKeys=section_keys,
        errors=[],
    )
    with _JOB_STORE_LOCK:
        _JOB_STORE[job_id] = {
            "job": job_out,
            "dossier_run_id": None,
            "output_paths": {},
        }

    thread = threading.Thread(
        target=_run_section_bundle_job,
        args=(job_id, gene, section_keys),
        daemon=True,
        name=f"dossier-job-{job_id}",
    )
    thread.start()

    logger.info(
        "Job %s enqueued gene=%s sections=%s thread=%s",
        job_id, gene, section_keys, thread.name,
    )
    return job_out


@api_router.get("/jobs/{job_id}", response_model=WorkflowJobOut)
@app.get("/jobs/{job_id}", response_model=WorkflowJobOut)
def get_job_endpoint(job_id: str) -> WorkflowJobOut:
    with _JOB_STORE_LOCK:
        entry = _JOB_STORE.get(job_id)
    if entry:
        return entry["job"]
    raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")


@api_router.get("/jobs/{job_id}/artifacts", response_model=JobArtifactsResponseOut)
@app.get("/jobs/{job_id}/artifacts", response_model=JobArtifactsResponseOut)
def get_job_artifacts_endpoint(job_id: str) -> JobArtifactsResponseOut:
    with _JOB_STORE_LOCK:
        entry = _JOB_STORE.get(job_id)
    if not entry:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    output_paths = entry.get("output_paths", {})
    dossier_run_id = entry.get("dossier_run_id")
    job = entry.get("job")
    report: ReportArtifactOut | None = None
    if job and job.artifactIds:
        try:
            report = handle_get_report(job.artifactIds[0])
        except Exception:
            report = None

    supplementary: list[ArtifactOut] = []
    for key, path_str in output_paths.items():
        if report and key in {"html", "pdf"}:
            continue
        p = Path(path_str)
        supplementary.append(
            ArtifactOut(
                id=key,
                path=path_str,
                artifactType=key,
                exists=p.exists(),
            )
        )
    return JobArtifactsResponseOut(
        jobId=job_id,
        dossierRunId=dossier_run_id,
        report=report,
        supplementaryArtifacts=supplementary,
    )


@api_router.get("/reports", response_model=list[ReportArtifactOut])
@app.get("/reports", response_model=list[ReportArtifactOut])
def list_reports_endpoint() -> list[ReportArtifactOut]:
    return handle_list_reports()


@api_router.get("/reports/{report_id}", response_model=ReportArtifactOut)
@app.get("/reports/{report_id}", response_model=ReportArtifactOut)
def get_report_endpoint(report_id: str) -> ReportArtifactOut:
    return handle_get_report(report_id)


@api_router.get("/reports/{report_id}/html")
@app.get("/reports/{report_id}/html")
def get_report_html_endpoint(report_id: str):
    return handle_get_report_html(report_id)


@api_router.get("/reports/{report_id}/pdf")
@app.get("/reports/{report_id}/pdf")
def get_report_pdf_endpoint(report_id: str):
    return handle_get_report_pdf(report_id)


@api_router.get("/history")
@app.get("/history")
def list_history_endpoint():
    init_db()
    with session_scope() as session:
        runs = session.exec(
            select(DossierRunRow).order_by(DossierRunRow.started_at.desc()).limit(20)
        ).all()
        return [
            {
                "id": r.id,
                "geneLabel": r.gene_symbol,
                "workflow": r.run_type,
                "status": "Completed" if r.status in {"completed", "running"} else "Failed",
                "createdAt": r.started_at.isoformat() if r.started_at else datetime.now().isoformat(),
            }
            for r in runs
        ]


@api_router.get("/recent")
@app.get("/recent")
def list_recent_endpoint():
    return [
        {"id": "rw-1", "label": "SREBF2 Target Dossier", "href": "/genes/SREBF2"},
        {"id": "rw-2", "label": "CDH10 Target Dossier", "href": "/genes/CDH10"},
        {"id": "rw-3", "label": "SREBF2 HD Report", "href": "/reports/rep-srebf2"},
    ]


@api_router.post("/ask", response_model=AskResponseOut)
@app.post("/ask", response_model=AskResponseOut)
def ask_question_endpoint(body: AskRequest) -> AskResponseOut:
    return handle_ask_question(body)


@api_router.post("/compare", response_model=ComparisonResponseOut)
@app.post("/compare", response_model=ComparisonResponseOut)
def compare_genes_endpoint(body: CompareRequest) -> ComparisonResponseOut:
    return handle_compare_genes(body)


app.include_router(api_router)


# --------------------------------------------------------------------------------------
# Existing standalone routes
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
