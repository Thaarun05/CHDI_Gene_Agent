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
    WorkflowTransientContext,
    _client_opentargets,
    _client_pubmed,
    _client_reactome,
    _coverage_updates_from_state,
    extract_gene_ids_from_tool_result,
    node_call_source_clients,
    node_save_raw_artifacts,
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
                "nomenclaturename": "sterol regulatory element binding transcription factor 2",
                "otheraliases": "SREBP2, SREBP-2",
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
            "entries": [
                {
                    "primaryAccession": "Q12772",
                    "organism_id": 9606,
                    "gene_synonyms": ["SREBP2", "bHLHd2"],
                }
            ],
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
    assert gene_ids["full_name"] == "sterol regulatory element binding transcription factor 2"
    assert gene_ids["pubmed_aliases"] == ["SREBP2", "SREBP-2", "bHLHd2"]


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


def test_pubmed_client_uses_identity_aliases_not_planner_values(monkeypatch):
    from gene_dossier.tools import pubmed

    captured: dict[str, object] = {}

    def fake_search_hd_literature(gene_symbol, **kwargs):
        captured.update(kwargs)
        return ToolResult(
            source_name="PubMed",
            endpoint_name="search_hd_literature",
            success=True,
            gene_symbol=gene_symbol,
            request_url="https://example.test/pubmed",
            data={"pmids": []},
        )

    monkeypatch.setattr(pubmed, "search_hd_literature", fake_search_hd_literature)

    _client_pubmed(
        gene_symbol="SREBF2",
        gene_ids={
            "pubmed_aliases": ["SREBP2"],
            "full_name": "sterol regulatory element binding transcription factor 2",
            "planner_aliases": ["MADEUP"],
        },
        settings=get_settings(),
    )

    assert captured["aliases"] == ["SREBP2"]
    assert captured["full_name"] == "sterol regulatory element binding transcription factor 2"


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


def test_ucsc_live_figure_stays_out_of_serializable_state(tmp_path: Path, monkeypatch):
    from gene_dossier.tools.ucsc import UCSCLiveFigurePayload, UCSCWorkflowExecution
    from gene_dossier.ucsc_figure import resolve_artifact_path, sha256_hex
    from gene_dossier.normalize.ucsc_conservation import normalize_ucsc_conservation

    transient = WorkflowTransientContext()
    png = (Path(__file__).parent / "fixtures" / "ucsc" / "srebf2_comprehensive_conservation.png").read_bytes()
    expected_sha = "3d165b72c20d11a0c921d16bf2cd17418a5169c2d0cec0537e297de5be0e3d6a"
    execution = UCSCWorkflowExecution(
        tool_result=ToolResult(
            source_name="UCSC",
            endpoint_name="fetch_gene_region",
            success=True,
            gene_symbol="SREBF2",
            request_url="https://api.genome.ucsc.edu/getData/track",
            request_params={"genome": "hg38", "track": "knownGene"},
            data={
                "gene_symbol": "SREBF2",
                "genome": "hg38",
                "search": {},
                "track_data": {},
                "figure": {
                    "status": "ok",
                    "display_position": "chr22:41833105-41907305",
                    "selected_transcript": "ENST00000361204.9",
                },
            },
        ),
        live_figure=UCSCLiveFigurePayload(
            content=png,
            media_type="image/png",
            width=1436,
            height=1192,
            sha256=expected_sha,
            byte_size=len(png),
            request_chain=(
                {
                    "endpoint_name": "hgRenderTracks",
                    "request_url": "https://genome.ucsc.edu/cgi-bin/hgRenderTracks",
                    "request_params": {"db": "hg38", "position": "chr22:41833105-41907305"},
                    "status_code": 200,
                    "success": True,
                    "error_type": None,
                    "error_message": None,
                    "parent_request_index": None,
                },
            ),
            image_request_index=0,
            track_preset_id="ucsc_section_1b_comprehensive_v1",
            track_preset_version=1,
            track_params={"db": "hg38", "position": "chr22:41833105-41907305", "pix": "1400"},
            genome="hg38",
            display_position="chr22:41833105-41907305",
            selected_transcript="ENST00000361204.9",
        ),
    )

    def _fake_fetch(*args, **kwargs):
        return execution

    import gene_dossier.tools.ucsc as ucsc

    monkeypatch.setattr(ucsc, "fetch_gene_region_execution", _fake_fetch)
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")
    state = node_call_source_clients(
        {
            "gene_symbol": "SREBF2",
            "dossier_run_id": "wf-ucsc-live",
            "gene_ids": {},
            "tool_results": [],
            "errors": [],
        },
        settings=settings,
        call_network=True,
        sources=["UCSC"],
        transient=transient,
    )
    tool_result = state["tool_results"][0]
    token = tool_result.data["_transient_ucsc_figure_token"]
    assert token.startswith("wf-ucsc-live:")
    assert len(token.split(":")[1]) == 32
    assert transient.live_figures[token].content == png
    assert "content" not in str(tool_result.data)

    saved = node_save_raw_artifacts(
        state,
        settings=settings,
        persist_db=False,
        transient=transient,
    )
    persisted = saved["tool_results"][0]
    assert "_transient_ucsc_figure_token" not in persisted.data
    assert transient.live_figures == {}
    rel = persisted.data["figure"]["local_artifact_path"]
    assert rel.endswith(".png")
    final_path = resolve_artifact_path(rel, root=settings.raw_data_path)
    assert final_path.is_file()
    assert sha256_hex(final_path.read_bytes()) == expected_sha
    assert any(a.get("artifact_type") == "image" for a in saved["raw_artifacts"])
    assert persisted.data["figure"]["figure_api_run_id"]
    assert persisted.data["figure"]["figure_raw_artifact_id"]

    # Four UCSC fact types after normalization, with figure artifact ID wired.
    import json

    search = json.loads(
        (Path(__file__).parent / "fixtures" / "ucsc" / "srebf2_search_relevant.json").read_text()
    )
    track = json.loads(
        (Path(__file__).parent / "fixtures" / "ucsc" / "srebf2_known_gene_region.json").read_text()
    )
    persisted.data["search"] = search
    persisted.data["track_data"] = track
    records = normalize_ucsc_conservation(
        persisted,
        dossier_run_id="wf-ucsc-live",
        api_run_id="json-api",
        raw_artifact_id="json-art",
    )
    types = {r.fact_type for r in records}
    assert types == {
        "ucsc_gene_locus",
        "ucsc_transcript_inventory",
        "ucsc_canonical_transcript",
        "ucsc_conservation_figure",
    }
    fig = next(r for r in records if r.fact_type == "ucsc_conservation_figure")
    assert fig.raw_artifact_id == persisted.data["figure"]["figure_raw_artifact_id"]
    assert fig.api_run_id == persisted.data["figure"]["figure_api_run_id"]


def test_empty_request_chain_is_rejected(tmp_path: Path):
    from gene_dossier.tools.ucsc import UCSCLiveFigurePayload
    from gene_dossier.raw_store import RawStore
    from gene_dossier.workflow import _persist_ucsc_result_with_live_figure

    png = (Path(__file__).parent / "fixtures" / "ucsc" / "srebf2_comprehensive_conservation.png").read_bytes()
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")
    transient = WorkflowTransientContext()
    token = "run-empty:abc"
    transient.put_figure(
        token,
        UCSCLiveFigurePayload(
            content=png,
            media_type="image/png",
            width=100,
            height=100,
            sha256="x",
            byte_size=len(png),
            request_chain=(),
            image_request_index=0,
            track_preset_id="ucsc_section_1b_comprehensive_v1",
            track_preset_version=1,
            track_params={},
            genome="hg38",
            display_position="chr22:1-2",
            selected_transcript=None,
        ),
    )
    result = ToolResult(
        source_name="UCSC",
        endpoint_name="fetch_gene_region",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test",
        data={"_transient_ucsc_figure_token": token, "figure": {"status": "ok"}},
    )
    try:
        _persist_ucsc_result_with_live_figure(
            result=result,
            dossier_run_id="run-empty",
            gene_symbol="SREBF2",
            settings=settings,
            store=RawStore(base_dir=settings.raw_data_path),
            transient=transient,
            persist_db=False,
        )
        assert False, "expected RuntimeError"
    except RuntimeError as exc:
        assert "missing HTTP request provenance" in str(exc)
    assert transient.live_figures == {}


def test_failed_final_move_cleans_new_json_and_image(tmp_path: Path, monkeypatch):
    from gene_dossier.tools.ucsc import UCSCLiveFigurePayload
    from gene_dossier.raw_store import RawStore
    from gene_dossier.workflow import _persist_ucsc_result_with_live_figure

    png = (Path(__file__).parent / "fixtures" / "ucsc" / "srebf2_comprehensive_conservation.png").read_bytes()
    expected_sha = "3d165b72c20d11a0c921d16bf2cd17418a5169c2d0cec0537e297de5be0e3d6a"
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")
    store = RawStore(base_dir=settings.raw_data_path)
    transient = WorkflowTransientContext()
    token = "run-fail:abcdef0123456789abcdef0123456789"
    transient.put_figure(
        token,
        UCSCLiveFigurePayload(
            content=png,
            media_type="image/png",
            width=80,
            height=40,
            sha256=expected_sha,
            byte_size=len(png),
            request_chain=(
                {
                    "endpoint_name": "hgRenderTracks",
                    "request_url": "https://genome.ucsc.edu/cgi-bin/hgRenderTracks",
                    "request_params": {"db": "hg38"},
                    "status_code": 200,
                    "success": True,
                    "error_type": None,
                    "error_message": None,
                    "parent_request_index": None,
                },
            ),
            image_request_index=0,
            track_preset_id="ucsc_section_1b_comprehensive_v1",
            track_preset_version=1,
            track_params={"pix": "1400"},
            genome="hg38",
            display_position="chr22:1-2",
            selected_transcript=None,
        ),
    )

    original_replace = Path.replace

    def _boom(self, target):  # noqa: ANN001
        raise OSError("simulated final move failure")

    monkeypatch.setattr(Path, "replace", _boom)
    result = ToolResult(
        source_name="UCSC",
        endpoint_name="fetch_gene_region",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test",
        data={
            "_transient_ucsc_figure_token": token,
            "gene_symbol": "SREBF2",
            "figure": {"status": "ok"},
        },
    )
    try:
        _persist_ucsc_result_with_live_figure(
            result=result,
            dossier_run_id="run-fail",
            gene_symbol="SREBF2",
            settings=settings,
            store=store,
            transient=transient,
            persist_db=False,
        )
        assert False, "expected OSError"
    except OSError:
        pass
    finally:
        monkeypatch.setattr(Path, "replace", original_replace)

    # Newly written JSON and any temp/final image from this attempt must be gone.
    ucsc_dir = settings.raw_data_path / "run-fail" / "ucsc"
    remaining = list(ucsc_dir.rglob("*")) if ucsc_dir.exists() else []
    remaining_files = [p for p in remaining if p.is_file()]
    assert remaining_files == []
    assert transient.live_figures == {}


def test_concurrent_same_gene_transient_isolation():
    from gene_dossier.tools.ucsc import UCSCLiveFigurePayload

    ctx = WorkflowTransientContext()
    png = b"\x89PNG\r\n\x1a\n" + b"0" * 100
    a = UCSCLiveFigurePayload(
        content=png,
        media_type="image/png",
        width=10,
        height=10,
        sha256="aaa",
        byte_size=len(png),
        request_chain=({"endpoint_name": "hgRenderTracks", "success": True, "status_code": 200},),
        image_request_index=0,
        track_preset_id="p",
        track_preset_version=1,
        track_params={},
        genome="hg38",
        display_position="chr1:1-2",
        selected_transcript=None,
    )
    b = UCSCLiveFigurePayload(
        content=png + b"1",
        media_type="image/png",
        width=10,
        height=10,
        sha256="bbb",
        byte_size=len(png) + 1,
        request_chain=({"endpoint_name": "hgRenderTracks", "success": True, "status_code": 200},),
        image_request_index=0,
        track_preset_id="p",
        track_preset_version=1,
        track_params={},
        genome="hg38",
        display_position="chr1:1-2",
        selected_transcript=None,
    )
    tok_a = "runA:11111111111111111111111111111111"
    tok_b = "runB:22222222222222222222222222222222"
    ctx.put_figure(tok_a, a)
    ctx.put_figure(tok_b, b)
    assert ctx.pop_figure(tok_a).sha256 == "aaa"
    assert ctx.pop_figure(tok_a) is None
    assert ctx.pop_figure(tok_b).sha256 == "bbb"
    ctx.put_figure("runA:zzzz", a)
    ctx.clear_run("runA")
    assert all(not k.startswith("runA:") for k in ctx.live_figures)


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
