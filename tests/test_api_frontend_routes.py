"""Regression tests for frontend-facing API handlers."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

from gene_dossier.api import main as api_main
from gene_dossier import synthesis
from gene_dossier.models import AssertionType, EvidenceGrade, EvidenceRecord, SourceType
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
    gene: str = "CDH10",
) -> dict[str, object]:
    output_dir = tmp_path / run_id
    output_dir.mkdir()
    presentation_path = output_dir / "section_1.json"
    audit_path = output_dir / "section_1_audit.json"
    html_path = output_dir / "section_1.html"
    presentation_path.write_text("{}\n", encoding="utf-8")
    audit_path.write_text("{}\n", encoding="utf-8")
    if include_html:
        html_path.write_text("<html><body>Fresh report</body></html>", encoding="utf-8")

    output_paths = {
        "section_1_json": presentation_path,
        "section_1_audit_json": audit_path,
    }
    if include_html:
        output_paths["section_1_html"] = html_path

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

    monkeypatch.setattr(
        api_main,
        "run_section_bundle",
        lambda *_args, **_kwargs: result,
    )
    api_main._run_section_bundle_job(job_id, gene, ["1a", "5a", "5b"])
    return {
        "job_id": job_id,
        "report_id": report_id,
        "run_id": run_id,
        "html_path": html_path,
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


def test_fresh_job_artifacts_returns_generated_report(completed_fresh_job) -> None:
    artifacts = api_main.get_job_artifacts_endpoint(str(completed_fresh_job["job_id"]))

    assert artifacts.dossierRunId == "fresh-run-test-123"
    assert artifacts.report is not None
    assert artifacts.report.id == "report-fresh-run-test-123"
    assert artifacts.report.geneSymbol == "CDH10"
    assert artifacts.report.htmlUrl.endswith(
        "/api/reports/report-fresh-run-test-123/html"
    )
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


def test_generated_report_pdf_is_unavailable(completed_fresh_job) -> None:
    with pytest.raises(api_main.HTTPException) as exc_info:
        api_main.handle_get_report_pdf(str(completed_fresh_job["report_id"]))

    assert exc_info.value.status_code == 404


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
    monkeypatch.setattr(api_main, "_load_domain_records_for_universe", lambda *_args, **_kwargs: records)
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
    monkeypatch.setattr(api_main, "_try_grounded_llm_summary", lambda **_kwargs: (None, "deterministic"))

    response = api_main.handle_ask_question(
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
        rec.model_copy(update={"dossier_run_id": "tool-run-ppi"})
        for rec in baseline_records
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
    monkeypatch.setattr(api_main, "_try_grounded_llm_summary", lambda **_kwargs: (None, "deterministic"))

    response = api_main.handle_ask_question(
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
    monkeypatch.setattr(api_main, "_persistent_chroma_index", lambda: _FakeRealIndex([semantic_hit]))

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
    source = (
        api_main.PROJECT_ROOT / "frontend" / "src" / "pages" / "GenerateDossierPage.tsx"
    ).read_text(encoding="utf-8")

    assert "artifact?.id ?? 'rep-srebf2'" not in source
    assert "Report artifact is not available." in source


def test_ask_frontend_preserves_gene_activity_and_ordinal_contracts() -> None:
    frontend_root = api_main.PROJECT_ROOT / "frontend" / "src"
    ask_source = (frontend_root / "pages" / "AskPage.tsx").read_text(encoding="utf-8")
    composer_source = (frontend_root / "components" / "SearchComposer.tsx").read_text(
        encoding="utf-8"
    )
    activity_source = (frontend_root / "components" / "AgentActivity.tsx").read_text(
        encoding="utf-8"
    )

    assert "const [selectedGene, setSelectedGene]" in ask_source
    assert "const initialQ = params.get('q') ?? ''" in ask_source
    assert "useState<AskGene>(initialGene)" in ask_source
    assert "if (initialQ.trim())" in ask_source
    assert "void submit(initialQ)" in ask_source
    assert "What evidence suggests SREBF2 can be pharmacologically manipulated?" not in ask_source
    assert "const requestGene = selectedGene" in ask_source
    assert "askEvidenceQuestion(q, requestGene" in ask_source
    assert "requestGeneration.current" in ask_source
    assert "setResponse(null)" in ask_source
    assert "new URLSearchParams(current)" in ask_source
    assert "next.set('gene', gene)" in ask_source
    assert "response.citations[ordinal - 1]" in ask_source
    assert "<AgentActivity steps={response.agentActivity} />" in ask_source
    assert "onSelectGene={selectGene}" in ask_source
    assert "onSelectGene(g)" in composer_source
    assert "const isActive = loading && i === steps.length - 1" in activity_source
    assert "const isLast = i === steps.length - 1" not in activity_source


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
    response = api_main.handle_compare_genes(
        api_main.CompareRequest(genes=["SREBF2", "CDH10"])
    )
    chemical_row = next(
        row for row in response.matrix if row["dimension"] == "Chemical Perturbations"
    )

    assert response.evidenceUniverses["SREBF2"].dossierRunIds == [
        "407e1a4293c6424e8b6b830a1f0a7c60"
    ]
    assert response.evidenceUniverses["CDH10"].dossierRunIds == [
        "d94f392f4a3941d5a59f697f58d18234"
    ]
    assert chemical_row["cells"]["SREBF2"].evidenceCount == 0
    assert chemical_row["cells"]["CDH10"].evidenceCount == 0
