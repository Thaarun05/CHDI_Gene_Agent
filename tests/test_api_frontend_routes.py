"""Regression tests for frontend-facing API handlers."""

from __future__ import annotations

import json
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

from gene_dossier.api import main as api_main
from gene_dossier.agent.evidence import public_evidence_reference, public_run_reference
from gene_dossier.agent.models import (
    AgentEvidenceUniverse,
    AnswerMode,
    AnswerStatus,
    CapabilityId,
    EvidenceNeed,
    EvidenceRequirement,
    EvidenceRequirementAssessment,
    PlannerMethod,
    RequirementStatus,
    ScientificAgentResult,
    ScientificEntities,
    ScientificIntent,
    ScientificQuestionPlan,
)
from gene_dossier import synthesis
from gene_dossier.config import Settings
from gene_dossier.db import (
    get_engine,
    get_generated_report,
    init_db,
    save_dossier_run,
    save_evidence_record,
    session_scope,
)
from gene_dossier.models import (
    AssertionType,
    DossierRun,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.section_bundle import SectionBundleResult


class _Result:
    def __init__(self, rows: list[object]) -> None:
        self._rows = rows

    def all(self) -> list[object]:
        return self._rows


class _Session:
    def __init__(self, rows: list[object]) -> None:
        self._rows = rows

    def exec(self, _stmt: object) -> _Result:
        return _Result(self._rows)


class _DetachedAwareRow:
    def __init__(self, state: dict[str, bool], section: str, source_name: str) -> None:
        self._state = state
        self._section = section
        self._source_name = source_name

    def _assert_attached(self) -> None:
        if self._state["closed"]:
            raise AssertionError("row attribute accessed after session closed")

    @property
    def section(self) -> str:
        self._assert_attached()
        return self._section

    @property
    def source_name(self) -> str:
        self._assert_attached()
        return self._source_name


def test_gene_coverage_extracts_fields_before_session_closes(monkeypatch) -> None:
    state = {"closed": True}
    rows = [
        _DetachedAwareRow(state, "General gene information", "NCBI Gene"),
        _DetachedAwareRow(state, "CTD chemical perturbation", "CTD"),
        _DetachedAwareRow(state, "Chemical tool bioactivity", "ChEMBL"),
    ]

    @contextmanager
    def fake_session_scope():
        state["closed"] = False
        try:
            yield _Session(rows)
        finally:
            state["closed"] = True

    monkeypatch.setattr(api_main, "init_db", lambda: None)
    monkeypatch.setattr(api_main, "session_scope", fake_session_scope)

    coverage = api_main.handle_get_gene_coverage("SREBF2")
    by_category = {row.category: row for row in coverage.rows}

    assert by_category["Gene Identity"].status == "Available"
    assert by_category["Chemical Perturbations"].detail == (
        "CTD chemical-gene interaction records (1)"
    )
    assert by_category["Chemical Tools"].detail == "Chemical tool/bioactivity records (1)"


def _evidence_record(
    *,
    section: str,
    source_name: str,
    gene_symbol: str = "SREBF2",
    assertion_type: AssertionType = AssertionType.chemical_interaction,
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=f"{source_name.lower()}:{gene_symbol.lower()}:test",
        dossier_run_id="run-test",
        gene_symbol=gene_symbol,
        section=section,
        source_name=source_name,
        source_type=SourceType.curated_database,
        assertion_type=assertion_type,
        fact_type="test_fact",
        evidence_grade=EvidenceGrade.C,
        display_text="Stored test evidence.",
    )


def _run_fake_section_bundle_job(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    run_id: str,
    status: str,
    include_html: bool,
    include_pdf: bool = False,
    registration_failure: bool = False,
    gene: str = "CDH10",
) -> dict[str, object]:
    output_root = tmp_path / "outputs"
    output_dir = output_root / "section_validation" / gene / run_id
    output_dir.mkdir(parents=True)
    presentation_path = output_dir / "section_1.json"
    audit_path = output_dir / "section_1_audit.json"
    html_path = output_dir / "section_1.html"
    pdf_path = output_dir / "section_1.pdf"
    presentation_path.write_text("{}\n", encoding="utf-8")
    audit_path.write_text("{}\n", encoding="utf-8")
    if include_html:
        html_path.write_text("<html><body>Fresh report</body></html>", encoding="utf-8")
    if include_pdf:
        pdf_path.write_bytes(b"%PDF-1.4\n% exact run\n")

    output_paths = {
        "section_1_json": presentation_path,
        "section_1_audit_json": audit_path,
    }
    if include_html:
        output_paths["section_1_html"] = html_path
    if include_pdf:
        output_paths["section_1_pdf"] = pdf_path

    settings = Settings(
        database_url=f"sqlite:///{tmp_path / 'api-test.db'}",
        output_dir=output_root,
    )
    engine = get_engine(settings.database_url)
    init_db(engine)
    run = DossierRun(gene_symbol=gene, run_type="section_bundle", status=status)
    run.id = run_id
    with session_scope(engine) as session:
        save_dossier_run(session, run)
    monkeypatch.setattr(api_main, "get_settings", lambda: settings)
    monkeypatch.setattr(api_main, "init_db", lambda: init_db(engine))
    monkeypatch.setattr(api_main, "session_scope", lambda: session_scope(engine))

    result = SectionBundleResult(
        gene_symbol=gene,
        dossier_run_id=run_id,
        selected_section_keys=["1a", "5a", "5b"],
        output_dir=output_dir,
        output_paths=output_paths,
        status=status,
    )
    job_id = f"job-{run_id}"
    report_id = f"report-{run_id}"
    job = api_main.WorkflowJobOut(
        id=job_id,
        geneSymbol=gene,
        jobType="hd_dossier",
        status="Queued",
        stages=[],
        createdAt="2026-08-11T00:00:00+00:00",
        artifactIds=None,
        sectionKeys=["1a", "5a", "5b"],
        errors=[],
    )
    with api_main._JOB_STORE_LOCK:
        api_main._JOB_STORE[job_id] = {
            "job": job,
            "dossier_run_id": None,
            "output_paths": {},
        }
        api_main._GENERATED_REPORT_STORE.pop(report_id, None)

    run_kwargs: dict[str, object] = {}

    def fake_run(*_args, **kwargs):
        run_kwargs.update(kwargs)
        return result

    monkeypatch.setattr(api_main, "run_section_bundle", fake_run)
    if registration_failure:
        monkeypatch.setattr(
            api_main,
            "_persist_generated_report",
            lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("metadata boom")),
        )
    api_main._run_section_bundle_job(job_id, gene, ["1a", "5a", "5b"])
    return {
        "job_id": job_id,
        "report_id": report_id,
        "run_id": run_id,
        "html_path": html_path,
        "pdf_path": pdf_path if include_pdf else None,
        "run_kwargs": run_kwargs,
        "engine": engine,
    }


def _clean_fake_section_bundle_job(info: dict[str, object]) -> None:
    with api_main._JOB_STORE_LOCK:
        api_main._JOB_STORE.pop(str(info["job_id"]), None)
        api_main._GENERATED_REPORT_STORE.pop(str(info["report_id"]), None)


@pytest.fixture
def completed_fresh_job(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    info = _run_fake_section_bundle_job(
        monkeypatch,
        tmp_path,
        run_id="fresh-run-test-123",
        status="completed",
        include_html=True,
    )
    try:
        yield info
    finally:
        _clean_fake_section_bundle_job(info)


def test_fresh_job_uses_generated_report_id(completed_fresh_job) -> None:
    job = api_main.get_job_endpoint(str(completed_fresh_job["job_id"]))

    assert job.status == "Completed"
    assert job.dossierRunId == "fresh-run-test-123"
    assert job.artifactIds == ["report-fresh-run-test-123"]
    assert not any(key.startswith("section_1_") for key in job.artifactIds)
    assert completed_fresh_job["run_kwargs"]["write_pdf"] is True


def test_fresh_job_artifacts_returns_generated_report(completed_fresh_job) -> None:
    artifacts = api_main.get_job_artifacts_endpoint(str(completed_fresh_job["job_id"]))

    assert artifacts.dossierRunId == "fresh-run-test-123"
    assert artifacts.report is not None
    assert artifacts.report.id == "report-fresh-run-test-123"
    assert artifacts.report.geneSymbol == "CDH10"
    assert artifacts.report.htmlUrl.endswith("/api/reports/report-fresh-run-test-123/html")
    assert artifacts.report.pdfUrl is None


def test_generated_report_metadata_and_exact_html_resolve(completed_fresh_job) -> None:
    report_id = str(completed_fresh_job["report_id"])
    report = api_main.handle_get_report(report_id)
    response = api_main.handle_get_report_html(report_id)

    assert report.id == report_id
    assert report.geneSymbol == "CDH10"
    assert report.sections == ["Section 1a", "Section 5a", "Section 5b"]
    assert response.path == completed_fresh_job["html_path"]
    assert response.media_type == "text/html"


def test_generated_report_survives_runtime_cache_loss(completed_fresh_job) -> None:
    report_id = str(completed_fresh_job["report_id"])
    with api_main._JOB_STORE_LOCK:
        api_main._GENERATED_REPORT_STORE.clear()

    report = api_main.handle_get_report(report_id)
    response = api_main.handle_get_report_html(report_id)

    assert report.id == report_id
    assert report.reportOrigin == "generated"
    assert report.dossierRunId == completed_fresh_job["run_id"]
    assert response.path == completed_fresh_job["html_path"]
    with session_scope(completed_fresh_job["engine"]) as session:
        row = get_generated_report(session, report_id)
        assert row is not None
        assert not Path(row.html_path).is_absolute()
        assert row.html_path.endswith(f"/{completed_fresh_job['run_id']}/section_1.html")


def test_missing_generated_html_is_truthfully_unavailable(completed_fresh_job) -> None:
    report_id = str(completed_fresh_job["report_id"])
    Path(completed_fresh_job["html_path"]).unlink()
    with api_main._JOB_STORE_LOCK:
        api_main._GENERATED_REPORT_STORE.clear()

    report = api_main.handle_get_report(report_id)
    assert report.htmlUrl is None
    assert report.pdfUrl is None
    with pytest.raises(api_main.HTTPException) as exc_info:
        api_main.handle_get_report_html(report_id)
    assert exc_info.value.status_code == 404
    with session_scope(completed_fresh_job["engine"]) as session:
        assert get_generated_report(session, report_id) is not None


def test_generated_artifact_validation_rejects_unsafe_paths(
    completed_fresh_job,
    tmp_path: Path,
) -> None:
    run_id = str(completed_fresh_job["run_id"])
    root = api_main.get_settings().output_path

    with pytest.raises(ValueError, match="must be relative"):
        api_main._validate_generated_artifact(
            completed_fresh_job["html_path"],
            dossier_run_id=run_id,
            artifact_type="html",
            stored_relative=True,
        )

    sibling = root / "section_validation" / "CDH10" / f"{run_id}-other"
    sibling.mkdir(parents=True)
    sibling_html = sibling / "section_1.html"
    sibling_html.write_text("wrong run", encoding="utf-8")
    with pytest.raises(ValueError, match="exact dossier run"):
        api_main._validate_generated_artifact(
            sibling_html.relative_to(root),
            dossier_run_id=run_id,
            artifact_type="html",
            stored_relative=True,
        )

    outside = tmp_path / "outside" / "section_1.html"
    outside.parent.mkdir()
    outside.write_text("outside", encoding="utf-8")
    symlink_dir = root / "section_validation" / "CDH10" / run_id / "linked"
    symlink_dir.mkdir(parents=True)
    symlink = symlink_dir / "section_1.html"
    symlink.symlink_to(outside)
    with pytest.raises(ValueError, match="outside the output root"):
        api_main._validate_generated_artifact(
            symlink.relative_to(root),
            dossier_run_id=run_id,
            artifact_type="html",
            stored_relative=True,
        )

    wrong_type = root / "section_validation" / "CDH10" / run_id / "report.html"
    wrong_type.write_text("wrong name", encoding="utf-8")
    with pytest.raises(ValueError, match="section_1.html"):
        api_main._validate_generated_artifact(
            wrong_type.relative_to(root),
            dossier_run_id=run_id,
            artifact_type="html",
            stored_relative=True,
        )


def test_generated_report_id_never_falls_back_to_accepted(completed_fresh_job) -> None:
    with pytest.raises(api_main.HTTPException) as exc_info:
        api_main.handle_get_report("report-SREBF2-not-a-run")
    assert exc_info.value.status_code == 404


def test_reports_and_recent_include_durable_generated_report(completed_fresh_job) -> None:
    reports = api_main.handle_list_reports()
    by_id = {report.id: report for report in reports}
    assert by_id["rep-srebf2"].reportOrigin == "accepted"
    generated = by_id[str(completed_fresh_job["report_id"])]
    assert generated.reportOrigin == "generated"
    assert generated.dossierRunId == completed_fresh_job["run_id"]
    recent = api_main.list_recent_endpoint()
    assert any(item["href"] == f"/reports/{completed_fresh_job['report_id']}" for item in recent)
    assert generated.id != "rep-srebf2"


def test_recovery_report_resolves_durably_with_exact_artifact(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    run_id = "3b599459b3da4af8afcee4bf5891ad6d"
    report_id = f"report-{run_id}"
    output_root = api_main.PROJECT_ROOT / "data" / "outputs"
    exact_html = output_root / "section_validation" / "SREBF2" / run_id / "section_1.html"
    assert exact_html.is_file()
    settings = Settings(
        database_url=f"sqlite:///{tmp_path / 'recovery-test.db'}",
        output_dir=output_root,
    )
    engine = get_engine(settings.database_url)
    init_db(engine)
    run = DossierRun(gene_symbol="SREBF2", run_type="section_bundle", status="completed")
    run.id = run_id
    with session_scope(engine) as session:
        save_dossier_run(session, run)
    monkeypatch.setattr(api_main, "get_settings", lambda: settings)
    monkeypatch.setattr(api_main, "init_db", lambda: init_db(engine))
    monkeypatch.setattr(api_main, "session_scope", lambda: session_scope(engine))

    api_main._persist_generated_report(
        dossier_run_id=run_id,
        gene_symbol="SREBF2",
        status="Completed",
        created_at=run.started_at,
        sections=["Section 1a", "Section 7a"],
        html_path=exact_html,
        pdf_path=None,
        job_id="job-srebf2-7c461679",
    )
    with api_main._JOB_STORE_LOCK:
        api_main._GENERATED_REPORT_STORE.clear()

    report = api_main.handle_get_report(report_id)
    html_response = api_main.handle_get_report_html(report_id)
    reports = {item.id: item for item in api_main.handle_list_reports()}
    recent = api_main.list_recent_endpoint()

    assert report.id == report_id
    assert report.dossierRunId == run_id
    assert report.reportOrigin == "generated"
    assert report.htmlUrl == f"/api/reports/{report_id}/html"
    assert report.pdfUrl is None
    assert html_response.path == exact_html
    assert reports["rep-srebf2"].reportOrigin == "accepted"
    assert reports[report_id].reportOrigin == "generated"
    assert any(item["href"] == f"/reports/{report_id}" for item in recent)
    with pytest.raises(api_main.HTTPException):
        api_main.handle_get_report_pdf(report_id)


def test_generated_report_pdf_is_unavailable(completed_fresh_job) -> None:
    with pytest.raises(api_main.HTTPException) as exc_info:
        api_main.handle_get_report_pdf(str(completed_fresh_job["report_id"]))

    assert exc_info.value.status_code == 404


def test_future_generated_pdf_uses_exact_run_artifact(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    info = _run_fake_section_bundle_job(
        monkeypatch,
        tmp_path,
        run_id="future-pdf-run-123",
        status="completed",
        include_html=True,
        include_pdf=True,
    )
    try:
        report = api_main.handle_get_report(str(info["report_id"]))
        html_response = api_main.handle_get_report_html(str(info["report_id"]))
        pdf_response = api_main.handle_get_report_pdf(str(info["report_id"]))
        assert report.htmlUrl == "/api/reports/report-future-pdf-run-123/html"
        assert report.pdfUrl == "/api/reports/report-future-pdf-run-123/pdf"
        assert html_response.path == info["html_path"]
        assert pdf_response.path == info["pdf_path"]
        assert pdf_response.media_type == "application/pdf"
    finally:
        _clean_fake_section_bundle_job(info)


def test_generated_report_registration_failure_preserves_scientific_run(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    info = _run_fake_section_bundle_job(
        monkeypatch,
        tmp_path,
        run_id="registration-failure-run",
        status="completed",
        include_html=True,
        registration_failure=True,
    )
    try:
        job = api_main.get_job_endpoint(str(info["job_id"]))
        assert job.status == "Completed"
        assert job.artifactIds == []
        assert job.errors and "metadata boom" in job.errors[0]
        assert Path(info["html_path"]).is_file()
    finally:
        _clean_fake_section_bundle_job(info)


def test_generated_primary_html_is_not_supplementary(completed_fresh_job) -> None:
    artifacts = api_main.get_job_artifacts_endpoint(str(completed_fresh_job["job_id"]))
    supplementary_ids = {artifact.id for artifact in artifacts.supplementaryArtifacts}

    assert "section_1_html" not in supplementary_ids
    assert supplementary_ids == {"section_1_json", "section_1_audit_json"}


def test_partial_fresh_job_preserves_status_and_report(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    info = _run_fake_section_bundle_job(
        monkeypatch,
        tmp_path,
        run_id="fresh-partial-test-123",
        status="partial",
        include_html=True,
    )
    try:
        job = api_main.get_job_endpoint(str(info["job_id"]))
        artifacts = api_main.get_job_artifacts_endpoint(str(info["job_id"]))

        assert job.status == "Partial"
        assert job.artifactIds == ["report-fresh-partial-test-123"]
        assert artifacts.report is not None
        assert artifacts.report.status == "Partial"
    finally:
        _clean_fake_section_bundle_job(info)


def test_failed_fresh_job_does_not_advertise_report(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    info = _run_fake_section_bundle_job(
        monkeypatch,
        tmp_path,
        run_id="fresh-failed-test-123",
        status="failed",
        include_html=True,
    )
    try:
        job = api_main.get_job_endpoint(str(info["job_id"]))
        artifacts = api_main.get_job_artifacts_endpoint(str(info["job_id"]))

        assert job.status == "Failed"
        assert job.artifactIds == []
        assert artifacts.report is None
        with pytest.raises(api_main.HTTPException):
            api_main.handle_get_report(str(info["report_id"]))
    finally:
        _clean_fake_section_bundle_job(info)


def test_fresh_job_without_html_does_not_advertise_report(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    info = _run_fake_section_bundle_job(
        monkeypatch,
        tmp_path,
        run_id="fresh-no-html-test-123",
        status="completed",
        include_html=False,
    )
    try:
        job = api_main.get_job_endpoint(str(info["job_id"]))
        artifacts = api_main.get_job_artifacts_endpoint(str(info["job_id"]))

        assert job.status == "Completed"
        assert job.artifactIds == []
        assert artifacts.report is None
    finally:
        _clean_fake_section_bundle_job(info)


def test_pharmacology_question_scopes_to_chemical_evidence() -> None:
    assert api_main._is_pharmacology_question(
        "What evidence suggests SREBF2 can be pharmacologically manipulated?"
    )
    assert api_main._is_chemical_evidence_record(
        _evidence_record(section="Section 7a chemical tools", source_name="ChEMBL")
    )
    assert not api_main._is_chemical_evidence_record(
        _evidence_record(
            section="Expression Profile",
            source_name="Allen Brain Atlas",
            assertion_type=AssertionType.expression,
        )
    )


def test_demo_universe_defaults_to_baseline_only() -> None:
    universe = api_main.resolve_evidence_universe("CDH10")

    assert universe.baseEvidenceRunId == "d94f392f4a3941d5a59f697f58d18234"
    assert universe.toolRunIds == []
    assert universe.dossierRunIds == ["d94f392f4a3941d5a59f697f58d18234"]
    assert universe.evidenceUniverse == "accepted_demo"


def test_demo_universe_preserves_request_local_tool_overlay() -> None:
    universe = api_main.resolve_evidence_universe("CDH10", tool_run_ids=["new-tool-run"])

    assert universe.baseEvidenceRunId == "d94f392f4a3941d5a59f697f58d18234"
    assert universe.toolRunIds == ["new-tool-run"]
    assert universe.dossierRunIds == [
        "d94f392f4a3941d5a59f697f58d18234",
        "new-tool-run",
    ]
    assert universe.evidenceUniverse == "accepted_demo_with_tool_overlay"


def test_category_aware_sufficiency_requires_required_category() -> None:
    expression_hit = api_main.RetrievalHit(
        record=_evidence_record(
            section="Expression Profile",
            source_name="GTEx",
            assertion_type=AssertionType.expression,
        ),
        score=0.9,
        method="semantic",
        source_id="run-test:expr-1",
    )
    ppi_record = _evidence_record(
        section="Protein-protein interaction partners",
        source_name="STRING",
        assertion_type=AssertionType.ppi,
    )

    assert not api_main._has_sufficient_evidence(
        category="ppi",
        records=[ppi_record],
        hits=[expression_hit, expression_hit],
        min_hits=2,
    )


def test_controlled_tool_registry_selects_ppi_sections() -> None:
    category = api_main._infer_required_category("What proteins interact with CDH10?")
    planned = api_main._tool_for_category(category)

    assert category == "ppi"
    assert planned is not None
    tool_name, spec = planned
    assert tool_name == "get_ppi"
    assert spec["sectionKeys"] == ["5a", "5b"]


def test_ask_stored_ppi_rag_does_not_invoke_tool(monkeypatch) -> None:
    universe = api_main.resolve_evidence_universe("CDH10")
    records = [
        _evidence_record(
            section="Protein-protein interaction partners",
            source_name="STRING",
            gene_symbol="CDH10",
            assertion_type=AssertionType.ppi,
        ),
        _evidence_record(
            section="Protein-protein interaction partners",
            source_name="BioGRID",
            gene_symbol="CDH10",
            assertion_type=AssertionType.ppi,
        ),
    ]
    hits = [
        api_main.RetrievalHit(
            record=record,
            score=0.9,
            method="semantic",
            source_id=f"{record.dossier_run_id}:{record.id}",
        )
        for record in records
    ]

    monkeypatch.setattr(api_main, "init_db", lambda: None)
    monkeypatch.setattr(
        api_main, "_load_domain_records_for_universe", lambda *_args, **_kwargs: records
    )
    monkeypatch.setattr(
        api_main,
        "_retrieve_grounded_hits",
        lambda **_kwargs: (hits, "semantic", ["Semantic retrieval attempted first"], "real"),
    )
    monkeypatch.setattr(
        api_main,
        "_execute_controlled_tool",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("tool should not run")),
    )
    monkeypatch.setattr(
        api_main, "_try_grounded_llm_summary", lambda **_kwargs: (None, "deterministic")
    )

    response = api_main._handle_legacy_ask_question(
        api_main.AskRequest(question="What proteins interact with CDH10?", gene_symbol="CDH10")
    )

    assert response.retrievalMethod == "semantic"
    assert response.embeddingBackend == "real"
    assert response.toolsInvokedCount == 0
    assert response.baseEvidenceRunId == universe.baseEvidenceRunId
    assert response.dossierRunIds == universe.dossierRunIds


def test_ask_refresh_ppi_invokes_overlay_tool(monkeypatch) -> None:
    baseline_universe = api_main.resolve_evidence_universe("CDH10")
    baseline_records = [
        _evidence_record(
            section="Protein-protein interaction partners",
            source_name="STRING",
            gene_symbol="CDH10",
            assertion_type=AssertionType.ppi,
        ),
        _evidence_record(
            section="Protein-protein interaction partners",
            source_name="BioGRID",
            gene_symbol="CDH10",
            assertion_type=AssertionType.ppi,
        ),
    ]
    tool_records = [
        rec.model_copy(update={"dossier_run_id": "tool-run-ppi"}) for rec in baseline_records
    ]

    def fake_load(_gene, universe, **_kwargs):
        if "tool-run-ppi" in universe.dossierRunIds:
            return [*baseline_records, *tool_records]
        return baseline_records

    def fake_retrieve(**kwargs):
        records = kwargs["records"]
        hits = [
            api_main.RetrievalHit(
                record=record,
                score=0.9,
                method="semantic",
                source_id=f"{record.dossier_run_id}:{record.id}",
            )
            for record in records[:2]
        ]
        return hits, "semantic", ["Semantic retrieval attempted first"], "real"

    monkeypatch.setattr(api_main, "init_db", lambda: None)
    monkeypatch.setattr(api_main, "_load_domain_records_for_universe", fake_load)
    monkeypatch.setattr(api_main, "_retrieve_grounded_hits", fake_retrieve)
    monkeypatch.setattr(
        api_main,
        "_execute_controlled_tool",
        lambda _gene, tool_name, spec: {
            "toolName": tool_name,
            "sectionKeys": list(spec["sectionKeys"]),
            "status": "completed",
            "dossierRunId": "tool-run-ppi",
            "evidenceRecordsPersisted": 2,
            "indexedRecords": 2,
        },
    )
    monkeypatch.setattr(
        api_main, "_try_grounded_llm_summary", lambda **_kwargs: (None, "deterministic")
    )

    response = api_main._handle_legacy_ask_question(
        api_main.AskRequest(
            question="What proteins interact with CDH10?",
            gene_symbol="CDH10",
            refresh_if_available=True,
        )
    )

    assert response.toolsInvokedCount == 1
    assert response.embeddingBackend == "real"
    assert response.toolRunIds == ["tool-run-ppi"]
    assert response.dossierRunIds == [*baseline_universe.dossierRunIds, "tool-run-ppi"]
    assert response.evidenceUniverse == "accepted_demo_with_tool_overlay"
    assert response.toolActivity[0]["toolName"] == "get_ppi"
    assert response.toolActivity[0]["sectionKeys"] == ["5a", "5b"]


class _FakeIndexStatus:
    def __init__(self, embedding_backend: str) -> None:
        self.embedding_backend = embedding_backend


class _FakeHashIndex:
    available = True
    status = _FakeIndexStatus("hash_test_fallback")

    def upsert_evidence(self, _records):
        raise AssertionError("hash-backed vector retrieval should not run in demo path")


class _FakeRealIndex:
    available = True
    status = _FakeIndexStatus("real")

    def __init__(self, hits):
        self._hits = hits

    def upsert_evidence(self, _records):
        return 1

    def query(self, *_args, **_kwargs):
        return self._hits


class _FakeLocalMiniLMIndex(_FakeRealIndex):
    status = _FakeIndexStatus("local_minilm")


def test_hash_embedding_backend_forces_keyword_label(monkeypatch) -> None:
    record = _evidence_record(
        section="Protein-protein interaction partners",
        source_name="STRING",
        assertion_type=AssertionType.ppi,
    )
    monkeypatch.setattr(api_main, "_persistent_chroma_index", lambda: _FakeHashIndex())

    hits, method, _activity, backend = api_main._retrieve_grounded_hits(
        question="protein interaction partners",
        gene="SREBF2",
        universe=api_main.resolve_evidence_universe("SREBF2"),
        records=[record],
        category="ppi",
        limit=5,
    )

    assert method == "keyword"
    assert backend == "hash_test_fallback"
    assert hits


def test_real_embedding_backend_can_label_semantic(monkeypatch) -> None:
    record = _evidence_record(
        section="Protein-protein interaction partners",
        source_name="STRING",
        assertion_type=AssertionType.ppi,
    )
    semantic_hit = api_main.RetrievalHit(
        record=record,
        score=0.9,
        method="semantic",
        source_id=f"{record.dossier_run_id}:{record.id}",
    )
    monkeypatch.setattr(
        api_main, "_persistent_chroma_index", lambda: _FakeRealIndex([semantic_hit])
    )

    hits, method, _activity, backend = api_main._retrieve_grounded_hits(
        question="protein interaction partners",
        gene="SREBF2",
        universe=api_main.resolve_evidence_universe("SREBF2"),
        records=[record],
        category="ppi",
        limit=5,
    )

    assert method == "semantic"
    assert backend == "real"
    assert hits == [semantic_hit]


def test_local_minilm_backend_can_label_semantic(monkeypatch) -> None:
    record = _evidence_record(
        section="Protein-protein interaction partners",
        source_name="STRING",
        assertion_type=AssertionType.ppi,
    )
    semantic_hit = api_main.RetrievalHit(
        record=record,
        score=0.9,
        method="semantic",
        source_id=f"{record.dossier_run_id}:{record.id}",
    )
    monkeypatch.setattr(
        api_main,
        "_persistent_chroma_index",
        lambda: _FakeLocalMiniLMIndex([semantic_hit]),
    )

    hits, method, _activity, backend = api_main._retrieve_grounded_hits(
        question="protein interaction partners",
        gene="SREBF2",
        universe=api_main.resolve_evidence_universe("SREBF2"),
        records=[record],
        category="ppi",
        limit=5,
    )

    assert method == "semantic"
    assert backend == "local_minilm"
    assert hits == [semantic_hit]


def test_llm_output_with_valid_ordinal_citation_is_accepted(monkeypatch) -> None:
    valid = "a" * 32
    record = _evidence_record(
        section="Chemical tools",
        source_name="ChEMBL",
        assertion_type=AssertionType.chemical_tool,
    ).model_copy(update={"id": valid})
    hit = api_main.RetrievalHit(
        record=record,
        score=0.9,
        method="semantic",
        source_id=f"{record.dossier_run_id}:{record.id}",
    )

    class _Settings:
        def has_llm(self):
            return True

    class _Model:
        def invoke(self, _prompt):
            return "Grounded answer supported by stored chemical-tool evidence. [[1]]"

    monkeypatch.setattr(api_main, "get_settings", lambda: _Settings())
    monkeypatch.setattr(
        synthesis,
        "build_chat_model_candidates",
        lambda _settings: [SimpleNamespace(provider="fake", model=_Model())],
    )

    summary, method = api_main._try_grounded_llm_summary(
        question="Can this be manipulated?",
        gene="SREBF2",
        hits=[hit],
    )

    assert summary == "Grounded answer supported by stored chemical-tool evidence. [[1]]"
    assert method == "grounded_llm"


@pytest.mark.parametrize(
    "response_template",
    [
        "Grounded answer without a citation marker.",
        "Grounded answer with an out-of-range marker. [[2]]",
        "Grounded answer exposes EvidenceRecord {valid}. [[1]]",
        "Grounded answer exposes invented ID {invented}. [[1]]",
    ],
)
def test_invalid_llm_citation_output_falls_back(monkeypatch, response_template: str) -> None:
    valid = "a" * 32
    invented = "b" * 32
    record = _evidence_record(
        section="Chemical tools",
        source_name="ChEMBL",
        assertion_type=AssertionType.chemical_tool,
    ).model_copy(update={"id": valid})
    hit = api_main.RetrievalHit(
        record=record,
        score=0.9,
        method="semantic",
        source_id=f"{record.dossier_run_id}:{record.id}",
    )

    class _Settings:
        def has_llm(self):
            return True

    class _Model:
        def invoke(self, _prompt):
            return response_template.format(valid=valid, invented=invented)

    monkeypatch.setattr(api_main, "get_settings", lambda: _Settings())
    monkeypatch.setattr(
        synthesis,
        "build_chat_model_candidates",
        lambda _settings: [SimpleNamespace(provider="fake", model=_Model())],
    )

    summary, method = api_main._try_grounded_llm_summary(
        question="Can this be manipulated?",
        gene="SREBF2",
        hits=[hit],
    )

    assert summary is None
    assert method == "deterministic"


def test_deterministic_summary_uses_ordinals_and_redacts_record_ids() -> None:
    evidence_id = "a" * 32
    record = _evidence_record(
        section="Chemical tools",
        source_name="ChEMBL",
        assertion_type=AssertionType.chemical_tool,
    ).model_copy(
        update={
            "id": evidence_id,
            "display_text": f"EvidenceRecord {evidence_id} supports the retrieved claim.",
        }
    )
    hit = api_main.RetrievalHit(
        record=record,
        score=0.9,
        method="semantic",
        source_id=f"{record.dossier_run_id}:{record.id}",
    )

    summary = api_main._deterministic_grounded_summary(
        gene="SREBF2",
        hits=[hit],
        retrieval_method="semantic",
        category="chemical_tool",
    )

    assert "[[1]]" in summary
    assert evidence_id not in summary
    assert "[1]" not in summary.replace("[[1]]", "")


def test_accepted_job_preserves_selected_section_keys() -> None:
    job = api_main.create_job_endpoint(
        {
            "gene_symbol": "SREBF2",
            "use_existing_accepted": True,
            "sections": ["5a", "5b"],
        }
    )
    try:
        assert job.status == "Completed"
        assert job.sectionKeys == ["5a", "5b"]
        assert job.dossierRunId == "cb9030ab81dc42db80b81dd15d48e653"
        assert job.artifactIds == ["rep-srebf2"]
    finally:
        with api_main._JOB_STORE_LOCK:
            api_main._JOB_STORE.pop(job.id, None)


def test_accepted_cdh10_job_keeps_registered_report_id() -> None:
    job = api_main.create_job_endpoint(
        {
            "gene_symbol": "CDH10",
            "use_existing_accepted": True,
            "sections": ["5a", "5b"],
        }
    )
    try:
        assert job.status == "Completed"
        assert job.dossierRunId == "ae97cb43e4d94732b72ef86cecc3f40d"
        assert job.artifactIds == ["rep-cdh10"]
    finally:
        with api_main._JOB_STORE_LOCK:
            api_main._JOB_STORE.pop(job.id, None)


def test_generate_page_has_no_accepted_report_fallback() -> None:
    frontend_root = api_main.PROJECT_ROOT / "frontend" / "src"
    source = (frontend_root / "pages" / "GenerateDossierPage.tsx").read_text(encoding="utf-8")
    reports_source = (frontend_root / "pages" / "ReportsPage.tsx").read_text(encoding="utf-8")
    viewer_source = (frontend_root / "components" / "ReportViewer.tsx").read_text(encoding="utf-8")

    assert "artifact?.id ?? 'rep-srebf2'" not in source
    assert "Report artifact is not available." in source
    assert "const initialJobId = params.get('job')" in source
    assert "getJob(initialJobId)" in source
    assert "This job session has expired." in source
    assert "next.delete('job')" in source
    assert (
        "startDossierJob"
        not in source.split("useEffect(() => {", 1)[1].split(
            "}, [initialJobId, params, setParams])", 1
        )[0]
    )
    assert "PDF not generated" in source
    assert "href={report.pdfUrl || '#'}" not in reports_source
    assert "href={report.pdfUrl || '#'}" not in viewer_source
    assert "report.reportOrigin === 'generated'" in reports_source


def test_ask_frontend_preserves_gene_activity_and_ordinal_contracts() -> None:
    frontend_root = api_main.PROJECT_ROOT / "frontend" / "src"
    ask_source = (frontend_root / "pages" / "AskPage.tsx").read_text(encoding="utf-8")
    composer_source = (frontend_root / "components" / "SearchComposer.tsx").read_text(
        encoding="utf-8"
    )
    activity_source = (frontend_root / "components" / "AgentActivity.tsx").read_text(
        encoding="utf-8"
    )
    home_source = (frontend_root / "pages" / "HomePage.tsx").read_text(encoding="utf-8")
    client_source = (frontend_root / "api" / "client.ts").read_text(encoding="utf-8")

    assert "const [contextGene, setContextGene]" in ask_source
    assert "const initialQuestion = params.get('q') ?? ''" in ask_source
    assert "useState<string | null>(initialContext)" in ask_source
    assert "!initialQuestion.trim()" in ask_source
    assert "void submit(initialQuestion, initialContext, initialMode)" in ask_source
    assert "What evidence suggests SREBF2 can be pharmacologically manipulated?" not in ask_source
    assert "askEvidenceQuestion(question.trim(), selectedContext" in ask_source
    assert "requestGeneration.current" in ask_source
    assert "setResponse(null)" in ask_source
    assert "new URLSearchParams(current)" in ask_source
    assert "else next.delete('gene')" in ask_source
    assert "response.citations[ordinal - 1]" in ask_source
    assert "<AgentActivity steps={response.agentActivity} />" in ask_source
    assert "onSelectGene={selectContextGene}" in ask_source
    assert "Context Gene:" in (frontend_root / "components" / "SearchComposer.tsx").read_text(
        encoding="utf-8"
    )
    assert "evidenceSelection: 'accepted_or_latest_generated'" in ask_source
    assert "response.comparisonMatrix" in ask_source
    assert "response.evidenceGaps" in ask_source
    assert "response.answerSections" in ask_source
    assert "response.recommendations" in ask_source
    assert "result.failures" in ask_source
    assert 'label="Evidence" evidenceReference={item.public_evidence_ref}' in ask_source
    assert "label={`[${index + 1}]`}" not in ask_source
    assert "cells: Record<string, AgentComparisonCell>" in (
        frontend_root / "api" / "types.ts"
    ).read_text(encoding="utf-8")
    assert "onSelectGene?.(null)" in composer_source
    assert "No context gene" in composer_source
    assert "Stored Evidence Only" in composer_source
    assert "Math.min(Math.max(textarea.scrollHeight, 92), 240)" in composer_source
    assert "textarea.scrollHeight > 240 ? 'auto' : 'hidden'" in composer_source
    assert "onClick={() => onChange(suggestion)}" in composer_source
    assert "onClick={() => onSelectGene?.(suggestion)}" not in composer_source
    assert "const isActive = loading && i === steps.length - 1" in activity_source
    assert "useState<string | null>(null)" in home_source
    assert "if (contextGene) params.set('gene', contextGene)" in home_source
    assert "context_gene: contextGene" in client_source
    assert "geneSymbol = 'SREBF2'" not in client_source
    assert "gene === 'SREBF2'" in client_source
    assert "Mock mode has no qualifying persisted evidence" in client_source
    assert "const autoSubmitted = useRef(false)" in ask_source
    assert "autoSubmitted.current = true" in ask_source
    assert "const next = new URLSearchParams(current)" in ask_source
    assert "const REQUEST_TIMEOUT_MS = 120_000" in ask_source
    assert "backend_unavailable" in ask_source
    assert "provider_failure" in ask_source
    assert "cell.citationOrdinals.slice(0, 3)" in ask_source
    assert "View all {cell.evidenceCount}" in ask_source
    ordered_components = [
        "<StructuredAnswer",
        "<RequirementSummary",
        "<ComparisonDecisionSummary",
        "<HdModifierMatrix",
        "<EvidenceGapSummary",
        "<RecommendationSummary",
        "<Limitations",
        "<TechnicalDiagnostics",
    ]
    positions = [ask_source.index(component) for component in ordered_components]
    assert positions == sorted(positions)
    assert "const isLast = i === steps.length - 1" not in activity_source


def test_general_ask_adapter_exposes_authoritative_per_gene_contract(monkeypatch) -> None:
    record = _evidence_record(
        section="Protein-protein interaction partners",
        source_name="STRING",
        gene_symbol="CDH10",
        assertion_type=AssertionType.ppi,
    )
    requirement = EvidenceRequirement(
        id="ppi",
        label="Protein interactions",
        description="Qualifying CDH10 protein interaction evidence.",
        genes=["CDH10"],
        evidence_need=EvidenceNeed.protein_interaction,
        capability_ids=[CapabilityId.ppi],
        required=True,
        minimum_support=1,
        rationale="The question asks for interacting proteins.",
    )
    plan = ScientificQuestionPlan(
        intent=ScientificIntent.single_gene_question,
        entities=ScientificEntities(genes=["CDH10"]),
        primary_gene="CDH10",
        objective="Identify CDH10 interaction evidence.",
        analysis_lens="general",
        answer_mode=AnswerMode.fact,
        evidence_requirements=[requirement],
        requires_multi_gene=False,
        planner_method=PlannerMethod.deterministic_fallback,
    )
    result = ScientificAgentResult(
        status=AnswerStatus.answered,
        question="What proteins interact with CDH10?",
        context_gene="SREBF2",
        plan=plan,
        evidence_universes={
            "CDH10": AgentEvidenceUniverse(
                gene_symbol="CDH10",
                base_evidence_run_id="d94f",
                explicit_run_ids=["d94f"],
                dossier_run_ids=["d94f"],
                evidence_universe="accepted_demo",
            )
        },
        assessments=[
            EvidenceRequirementAssessment(
                requirement_id="ppi",
                gene_symbol="CDH10",
                evidence_need=EvidenceNeed.protein_interaction,
                required=True,
                minimum_support=1,
                status=RequirementStatus.sufficient,
                qualifying_count=1,
                evidence_record_ids=[record.id],
                contributing_capability_ids=[CapabilityId.ppi],
                detail="Threshold met.",
            )
        ],
        selected_records=[record],
        summary="Stored evidence supports a CDH10 interaction. [[1]]",
        retrieval_method="semantic",
        generation_method="hybrid",
        embedding_backend="local_minilm",
        agent_activity=["Completed"],
        metadata={
            "grounding": {
                "requestedSlotCount": 3,
                "acceptedSlotCount": 2,
                "fallbackSlotCount": 1,
                "diagnosticCounts": {"missing_slot": 1},
            }
        },
    )

    class FakeService:
        def __init__(self, **_kwargs):
            pass

        def execute(self, _request):
            return result

    monkeypatch.setattr(api_main, "ScientificAgentService", FakeService)
    response = api_main.handle_ask_question(
        api_main.AskRequest(
            question="What proteins interact with CDH10?",
            gene_symbol="SREBF2",
            context_gene="SREBF2",
        )
    )

    assert response.geneSymbol == "CDH10"
    assert response.contextGene == "SREBF2"
    assert response.evidenceUniverses["CDH10"].dossierRunIds == []
    assert response.evidenceUniverseRefs["CDH10"]["dossierRunRefs"] == [
        public_run_reference("d94f")
    ]
    assert response.evidenceRequirements[0]["required"] is True
    assert response.requirementAssessments[0]["status"] == "sufficient"
    assert response.citations[0].publicEvidenceRef == public_evidence_reference(record)
    assert response.citationRegistry == []
    assert response.evidenceCategories == []
    assert response.recommendations == []
    assert response.failures == []
    assert response.generationMethod == "hybrid"
    assert response.metadata["grounding"]["fallbackSlotCount"] == 1
    serialized = json.dumps(response.model_dump(mode="json"), sort_keys=True)
    assert record.id not in serialized
    assert "d94f" not in serialized
    assert response.citations[0].publicEvidenceRef in serialized


def test_ask_adapter_recursively_blocks_contextual_private_identifier_leak(monkeypatch) -> None:
    private_id = "private-contextual-record-id"

    class FakeService:
        def __init__(self, **_kwargs):
            pass

        def execute(self, request):
            return ScientificAgentResult(
                status=AnswerStatus.insufficient_evidence,
                question=request.question,
                context_gene=None,
                summary=f"Unsafe contextual payload {private_id}",
                private_identifiers={private_id},
            )

    monkeypatch.setattr(api_main, "ScientificAgentService", FakeService)

    with pytest.raises(api_main.HTTPException, match="serialized safely"):
        api_main.handle_ask_question(api_main.AskRequest(question="What evidence is available?"))


def test_public_evidence_detail_lookup_returns_public_safe_record(
    monkeypatch,
    tmp_path: Path,
) -> None:
    settings = Settings(
        database_url=f"sqlite:///{tmp_path / 'public-evidence.db'}",
        raw_data_dir=tmp_path / "raw",
        output_dir=tmp_path / "outputs",
        index_dir=tmp_path / "indexes",
    )
    engine = get_engine(settings.database_url)
    init_db(engine)
    record = _evidence_record(
        section="Protein interactions",
        source_name="STRING",
        gene_symbol="CDH10",
        assertion_type=AssertionType.ppi,
    )
    with session_scope(engine) as session:
        save_dossier_run(
            session,
            DossierRun(
                id="run-test",
                gene_symbol="CDH10",
                status="completed",
                run_type="offline_public_detail_test",
            ),
        )
        save_evidence_record(session, record)

    @contextmanager
    def disposable_session_scope():
        with session_scope(engine) as session:
            yield session

    monkeypatch.setattr(api_main, "init_db", lambda: None)
    monkeypatch.setattr(api_main, "session_scope", disposable_session_scope)
    public_ref = public_evidence_reference(record)

    detail = api_main.handle_get_evidence_record(public_ref)

    assert detail.id == public_ref
    assert detail.geneSymbol == "CDH10"
    assert detail.sourceName == "STRING"
    assert detail.sourceIdentifier == record.source_id
    assert detail.apiRunId is None
    assert detail.rawArtifactId is None


def test_ask_adapter_preserves_context_gene_omitted_null_and_explicit(monkeypatch) -> None:
    captured: list[object] = []
    plan = ScientificQuestionPlan(
        intent=ScientificIntent.single_gene_question,
        entities=ScientificEntities(genes=["MSH3"]),
        primary_gene="MSH3",
        objective="offline context test",
        analysis_lens="general",
        answer_mode=AnswerMode.synthesis,
        evidence_requirements=[],
        requires_multi_gene=False,
        planner_method=PlannerMethod.deterministic_fallback,
    )

    class FakeService:
        def __init__(self, **_kwargs):
            pass

        def execute(self, request):
            captured.append(request)
            return ScientificAgentResult(
                status=AnswerStatus.insufficient_evidence,
                question=request.question,
                context_gene=request.context_gene,
                plan=plan,
                evidence_universes={},
                summary="No evidence.",
                agent_activity=["done"],
            )

    monkeypatch.setattr(api_main, "ScientificAgentService", FakeService)

    omitted = api_main.handle_ask_question(
        api_main.AskRequest(question="What evidence links MSH3 to HD?")
    )
    explicit_null = api_main.handle_ask_question(
        api_main.AskRequest.model_validate(
            {"question": "What evidence links MSH3 to HD?", "context_gene": None}
        )
    )
    explicit_msh3 = api_main.handle_ask_question(
        api_main.AskRequest(question="What evidence links MSH3 to HD?", context_gene="MSH3")
    )
    explicit_srebf2 = api_main.handle_ask_question(
        api_main.AskRequest(question="What evidence links SREBF2 to HD?", context_gene="SREBF2")
    )

    assert [request.context_gene for request in captured] == [None, None, "MSH3", "SREBF2"]
    assert omitted.contextGene is None
    assert explicit_null.contextGene is None
    assert explicit_msh3.contextGene == "MSH3"
    assert explicit_srebf2.contextGene == "SREBF2"


def test_explicit_research_mode_supersedes_legacy_acquisition_flag(monkeypatch) -> None:
    captured: list[object] = []
    plan = ScientificQuestionPlan(
        intent=ScientificIntent.single_gene_question,
        entities=ScientificEntities(genes=["MSH3"]),
        primary_gene="MSH3",
        objective="offline research-mode precedence",
        analysis_lens="general",
        answer_mode=AnswerMode.synthesis,
        evidence_requirements=[],
        requires_multi_gene=False,
        planner_method=PlannerMethod.deterministic_fallback,
    )

    class FakeService:
        def __init__(self, **_kwargs):
            pass

        def execute(self, request):
            captured.append(request)
            return ScientificAgentResult(
                status=AnswerStatus.insufficient_evidence,
                question=request.question,
                context_gene=request.context_gene,
                plan=plan,
                summary="No evidence.",
            )

    monkeypatch.setattr(api_main, "ScientificAgentService", FakeService)
    cases = [
        api_main.AskRequest(
            question="MSH3 evidence",
            research_mode="auto",
            allow_tool_acquisition=False,
        ),
        api_main.AskRequest(
            question="MSH3 evidence",
            research_mode="deep_research",
            allow_tool_acquisition=False,
        ),
        api_main.AskRequest(
            question="MSH3 evidence",
            research_mode="stored_only",
            allow_tool_acquisition=True,
        ),
        api_main.AskRequest(question="MSH3 evidence", allow_tool_acquisition=False),
    ]
    for body in cases:
        api_main.handle_ask_question(body)

    assert [request.research_mode.value for request in captured] == [
        "auto",
        "deep_research",
        "stored_only",
        "auto",
    ]
    assert [request.allow_tool_acquisition for request in captured] == [
        True,
        True,
        False,
        False,
    ]


def test_accepted_report_html_resolves_validated_full_artifacts() -> None:
    expected = {
        "SREBF2": (
            api_main.PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "SREBF2_full_1a7a"
            / "407e1a4293c6424e8b6b830a1f0a7c60"
            / "section_1.html"
        ),
        "CDH10": (
            api_main.PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "CDH10_full_1a7a"
            / "d94f392f4a3941d5a59f697f58d18234"
            / "section_1.html"
        ),
    }

    for gene, expected_path in expected.items():
        report_id = api_main.DEMO_GENE_REGISTRY[gene]["report_id"]
        resolved_id, resolved_html, _ = api_main._find_report_files(gene)
        response = api_main.handle_get_report_html(report_id)

        assert expected_path.exists()
        assert api_main.DEMO_GENE_REGISTRY[gene]["html_path"] == expected_path
        assert resolved_id == report_id
        assert resolved_html == expected_path
        assert response.path == expected_path
        assert response.media_type == "text/html"
        assert response.headers["cache-control"] == "no-store, no-cache, must-revalidate"
        assert response.headers["pragma"] == "no-cache"


def test_accepted_report_pdf_resolves_beside_validated_html() -> None:
    expected = {
        "SREBF2": (
            api_main.PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "SREBF2_full_1a7a"
            / "407e1a4293c6424e8b6b830a1f0a7c60"
            / "section_1.pdf",
            api_main.PROJECT_ROOT / "SREBF2_report" / "SREBF2_report.pdf",
        ),
        "CDH10": (
            api_main.PROJECT_ROOT
            / "data"
            / "outputs"
            / "section_validation"
            / "CDH10_full_1a7a"
            / "d94f392f4a3941d5a59f697f58d18234"
            / "section_1.pdf",
            api_main.PROJECT_ROOT / "CDH10 report" / "CDH10_report.pdf",
        ),
    }

    for gene, (expected_pdf, old_pdf) in expected.items():
        report_id = api_main.DEMO_GENE_REGISTRY[gene]["report_id"]
        _, resolved_html, resolved_pdf = api_main._find_report_files(gene)
        response = api_main.handle_get_report_pdf(report_id)

        assert expected_pdf.exists()
        assert api_main.DEMO_GENE_REGISTRY[gene]["pdf_path"] == expected_pdf
        assert resolved_html is not None
        assert expected_pdf.parent == resolved_html.parent
        assert resolved_pdf == expected_pdf
        assert resolved_pdf != old_pdf
        assert response.path == expected_pdf
        assert response.media_type == "application/pdf"
        assert response.headers["cache-control"] == "no-store, no-cache, must-revalidate"
        assert response.headers["pragma"] == "no-cache"


def test_friday_baseline_chemical_perturbations_are_zero() -> None:
    srebf2 = api_main.handle_get_gene_coverage("SREBF2")
    cdh10 = api_main.handle_get_gene_coverage("CDH10")

    srebf2_by_category = {row.category: row for row in srebf2.rows}
    cdh10_by_category = {row.category: row for row in cdh10.rows}

    assert srebf2.baseEvidenceRunId == "407e1a4293c6424e8b6b830a1f0a7c60"
    assert cdh10.baseEvidenceRunId == "d94f392f4a3941d5a59f697f58d18234"
    assert srebf2_by_category["Chemical Perturbations"].status == "Not available"
    assert cdh10_by_category["Chemical Perturbations"].status == "Not available"


def test_compare_defaults_to_baseline_universes_only() -> None:
    response = api_main.handle_compare_genes(api_main.CompareRequest(genes=["SREBF2", "CDH10"]))
    chemical_row = next(
        row for row in response.matrix if row["dimension"] == "Chemical Perturbations"
    )

    assert response.evidenceUniverses["SREBF2"].dossierRunIds == [
        "407e1a4293c6424e8b6b830a1f0a7c60"
    ]
    assert response.evidenceUniverses["CDH10"].dossierRunIds == ["d94f392f4a3941d5a59f697f58d18234"]
    assert chemical_row["cells"]["SREBF2"].evidenceCount == 0
    assert chemical_row["cells"]["CDH10"].evidenceCount == 0
