"""Tests for LangGraph dossier workflow orchestration (offline / no network)."""

from __future__ import annotations

from pathlib import Path

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    ReportSection,
    SourceStatus,
    SourceType,
    ToolResult,
)
from gene_dossier.source_ids import make_source_id
from gene_dossier.workflow import (
    _client_opentargets,
    _client_reactome,
    _coverage_updates_from_state,
    extract_gene_ids_from_tool_result,
    node_index_evidence_in_chroma,
    node_render_outputs,
    node_resolve_gene_identity,
    run_gene_dossier_full_api_pass,
)


def _ncbi_tool_result() -> ToolResult:
    return ToolResult(
        source_name="NCBI Gene",
        endpoint_name="lookup_gene",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/ncbi",
        data={
            "gene_symbol": "SREBF2",
            "selected_gene_id": "6721",
            "selection_method": "exact_symbol",
            "selected_summary": {
                "nomenclaturesymbol": "SREBF2",
                "chromosome": "22",
                "uid": "6721",
            },
            "candidate_ids": ["6721"],
        },
    )


def _ensembl_tool_result() -> ToolResult:
    return ToolResult(
        source_name="Ensembl",
        endpoint_name="lookup_symbol",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/ensembl",
        data={
            "ensembl_id": "ENSG00000198911",
            "gene_symbol": "SREBF2",
        },
    )


def _uniprot_tool_result() -> ToolResult:
    return ToolResult(
        source_name="UniProt",
        endpoint_name="search_reviewed",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/uniprot",
        data={
            "selected_accession": "Q12772",
            "gene_symbol": "SREBF2",
        },
    )


def test_client_reactome_requires_uniprot_accession():
    settings = get_settings()
    missing = _client_reactome(
        gene_symbol="SREBF2", gene_ids={}, settings=settings
    )
    assert missing.success is False
    assert missing.error_type == "missing_identifier"

    # With accession present, do not claim missing_identifier (network not required here).
    # We only assert the guard path; live Reactome call is out of scope for unit tests.


def test_client_opentargets_requires_ensembl_id():
    settings = get_settings()
    missing = _client_opentargets(
        gene_symbol="SREBF2", gene_ids={}, settings=settings
    )
    assert missing.success is False
    assert missing.error_type == "missing_identifier"
    assert "ensembl_id" in (missing.error_message or "")


def test_extract_gene_ids_from_identity_tool_results():
    gene_ids: dict = {}
    for tr in (_ncbi_tool_result(), _ensembl_tool_result(), _uniprot_tool_result()):
        gene_ids = extract_gene_ids_from_tool_result(tr, gene_ids)
    assert gene_ids["entrez_gene_id"] == "6721"
    assert gene_ids["chromosome"] == "22"
    assert gene_ids["ensembl_id"] == "ENSG00000198911"
    assert gene_ids["uniprot_accession"] == "Q12772"
    assert gene_ids["official_symbol"] == "SREBF2"


def test_extract_gene_ids_from_gtex_derives_ensembl_id():
    """GTEx GENCODE versioned IDs should backfill bare ensembl_id when missing."""
    tr = ToolResult(
        source_name="GTEx",
        endpoint_name="median_expression",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/gtex",
        data={"gencode_id": "ENSG00000198911.11"},
    )
    gene_ids = extract_gene_ids_from_tool_result(tr, {})
    assert gene_ids["gtex_gencode_id"] == "ENSG00000198911.11"
    assert gene_ids["ensembl_id"] == "ENSG00000198911"


def test_preloaded_identity_results_populate_gene_ids_offline():
    """call_network=False must still harvest IDs from preloaded identity ToolResults."""
    settings = Settings()
    state = {
        "gene_symbol": "SREBF2",
        "dossier_run_id": "offline-id-1",
        "gene_ids": {},
        "tool_results": [
            _ncbi_tool_result(),
            _ensembl_tool_result(),
            _uniprot_tool_result(),
        ],
        "errors": [],
    }
    updated = node_resolve_gene_identity(
        state, settings=settings, call_network=False  # type: ignore[arg-type]
    )
    ids = updated["gene_ids"]
    assert ids["entrez_gene_id"] == "6721"
    assert ids["ensembl_id"] == "ENSG00000198911"
    assert ids["uniprot_accession"] == "Q12772"
    assert ids["chromosome"] == "22"
    assert updated.get("official_symbol") == "SREBF2"


def test_offline_srebf2_smoke_persist_db_false(tmp_path: Path):
    result = run_gene_dossier_full_api_pass(
        "SREBF2",
        call_network=False,
        preloaded_tool_results=[
            _ncbi_tool_result(),
            _ensembl_tool_result(),
            _uniprot_tool_result(),
        ],
        output_dir=tmp_path,
        write_pdf=False,
        persist_db=False,
        force_deterministic=True,
        dossier_run_id="wf-offline-srebf2",
    )
    assert result.status == "completed"
    assert result.gene_ids.get("entrez_gene_id") == "6721"
    assert result.gene_ids.get("ensembl_id") == "ENSG00000198911"
    assert result.gene_ids.get("uniprot_accession") == "Q12772"
    assert result.output_paths.get("rancho_html")
    assert result.output_paths.get("coverage_markdown")
    assert result.output_paths.get("debug_markdown")
    assert result.evidence_records  # NCBI identity normalized
    assert any("Chroma indexing" in note for note in result.synthesis_notes)


def test_node_index_evidence_in_chroma_soft_fails_without_crash():
    """Empty evidence and indexing path must never abort the workflow."""
    settings = Settings()
    empty = node_index_evidence_in_chroma(
        {
            "gene_symbol": "SREBF2",
            "dossier_run_id": "chroma-empty",
            "evidence_records": [],
            "synthesis_notes": [],
        },
        settings=settings,
    )
    assert any(
        "No evidence records available for Chroma indexing." in n
        for n in empty["synthesis_notes"]
    )

    from gene_dossier.models import (
        AssertionType,
        EvidenceGrade,
        EvidenceRecord,
        SourceType,
    )
    from gene_dossier.source_ids import make_source_id

    record = EvidenceRecord(
        source_id=make_source_id(
            "NCBI Gene", "SREBF2", AssertionType.gene_identity, "6721"
        ),
        dossier_run_id="chroma-one",
        gene_symbol="SREBF2",
        section="General",
        source_name="NCBI Gene",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="entrez_gene_id",
        evidence_grade=EvidenceGrade.A,
        display_text="SREBF2 Entrez Gene ID is 6721.",
    )
    populated = node_index_evidence_in_chroma(
        {
            "gene_symbol": "SREBF2",
            "dossier_run_id": "chroma-one",
            "evidence_records": [record],
            "synthesis_notes": [],
        },
        settings=settings,
    )
    notes = populated["synthesis_notes"]
    assert notes
    assert any("Chroma indexing" in n for n in notes)
    # Soft-fail contract: either complete, unavailable, or failed — never raise.
    assert any(
        n.startswith("Chroma indexing complete:")
        or n.startswith("Chroma indexing unavailable:")
        or n.startswith("Chroma indexing failed:")
        for n in notes
    )


def test_coverage_access_forbidden_classified_as_deferred():
    """BrainRNASeq-style access_forbidden should be deferred, not failed."""
    result = ToolResult(
        source_name="BrainRNASeq",
        endpoint_name="fetch_gene_expression",
        success=False,
        gene_symbol="SREBF2",
        request_url=(
            "https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-124.csv"
        ),
        status_code=403,
        data={
            "species": "human",
            "content_type": "text/html; charset=UTF-8",
            "raw_text_preview": "Just a moment...",
        },
        error_type="access_forbidden",
        error_message="HTTP 403 Forbidden",
    )
    updates = _coverage_updates_from_state(
        {
            "dossier_run_id": "cov-access-forbidden",
            "gene_symbol": "SREBF2",
            "tool_results": [result],
            "evidence_records": [],
            "raw_artifacts": [],
        }
    )
    assert len(updates) == 1
    row = updates[0]
    assert row.source_name == "BrainRNASeq"
    assert row.status == SourceStatus.deferred
    assert row.evidence_record_count == 0
    assert row.error_message == "HTTP 403 Forbidden"


def test_node_render_outputs_forwards_sections_to_rancho(tmp_path: Path, monkeypatch):
    """Finalized state['sections'] must reach Rancho builder unchanged."""
    settings = Settings()
    source_id = make_source_id(
        "NCBI Gene", "SREBF2", AssertionType.gene_identity, "6721"
    )
    evidence_records = [
        EvidenceRecord(
            source_id=source_id,
            dossier_run_id="wf-rancho-wire",
            gene_symbol="SREBF2",
            section="General",
            source_name="NCBI Gene",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="entrez_gene_id",
            evidence_grade=EvidenceGrade.A,
            display_text="SREBF2 Entrez Gene ID is 6721.",
        )
    ]
    sections = [
        ReportSection(
            dossier_run_id="wf-rancho-wire",
            section_name="General",
            content_markdown="Finalized synthesis narrative.",
            source_ids=[source_id],
            status="complete",
        )
    ]
    state = {
        "gene_symbol": "SREBF2",
        "dossier_run_id": "wf-rancho-wire",
        "gene_ids": {"chromosome": "22"},
        "evidence_records": evidence_records,
        "sections": sections,
        "coverage": [],
        "claims": [],
        "verification_results": [],
        "errors": [],
        "output_paths": {},
    }

    captured: dict = {}

    def fake_build_and_write_coverage(*args, **kwargs):
        return [], {"markdown": tmp_path / "coverage.md", "json": tmp_path / "coverage.json"}

    def fake_write_dossier_report(**kwargs):
        return {
            "markdown": tmp_path / "debug.md",
            "json": tmp_path / "debug.json",
        }

    def fake_rancho_builder(**kwargs):
        captured.clear()
        captured.update(kwargs)
        return object(), {"html": tmp_path / "rancho.html", "json": tmp_path / "rancho.json"}

    monkeypatch.setattr(
        "gene_dossier.workflow.build_and_write_coverage",
        fake_build_and_write_coverage,
    )
    monkeypatch.setattr(
        "gene_dossier.workflow.write_dossier_report",
        fake_write_dossier_report,
    )
    monkeypatch.setattr(
        "gene_dossier.workflow.build_and_write_rancho_report",
        fake_rancho_builder,
    )

    node_render_outputs(
        state,  # type: ignore[arg-type]
        settings=settings,
        output_dir=tmp_path,
        write_rancho=True,
        write_pdf=False,
        persist_db=False,
    )
    assert captured["report_sections"] is sections
    assert captured["evidence_records"] is evidence_records
    assert captured["write_pdf"] is False

    state["sections"] = []
    node_render_outputs(
        state,  # type: ignore[arg-type]
        settings=settings,
        output_dir=tmp_path,
        write_rancho=True,
        write_pdf=False,
        persist_db=False,
    )
    assert captured["report_sections"] is None
