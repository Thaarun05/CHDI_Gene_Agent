"""Regression tests for frontend-facing API handlers."""

from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace

from gene_dossier.api import main as api_main
from gene_dossier import synthesis
from gene_dossier.models import AssertionType, EvidenceGrade, EvidenceRecord, SourceType


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


def test_llm_output_with_invented_evidence_id_falls_back(monkeypatch) -> None:
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
            return f"Grounded answer cites EvidenceRecord {valid} and EvidenceRecord {invented}."

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


def test_accepted_job_preserves_selected_section_keys() -> None:
    job = api_main.create_job_endpoint(
        {
            "gene_symbol": "SREBF2",
            "use_existing_accepted": True,
            "sections": ["5a", "5b"],
        }
    )

    assert job.status == "Completed"
    assert job.sectionKeys == ["5a", "5b"]
    assert job.dossierRunId == "cb9030ab81dc42db80b81dd15d48e653"


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
