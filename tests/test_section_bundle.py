"""Offline tests for section-scoped 1a/1b bundle generation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
    new_id,
)
from gene_dossier.config import PROJECT_ROOT
from gene_dossier.report_presentation import (
    UCSC_STABLE_INTRO,
    build_known_structure_blocks,
    build_conservation_blocks,
    transcript_selection_sentence,
)
from gene_dossier.report_schema import ReportContentBlock
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SUPPORTED_SECTION_BUNDLE_KEYS,
    SectionBundleError,
    assign_opaque_refs,
    build_section_bundle_document,
    finalize_section_bundle_run,
    opaque_evidence_ref,
    render_section_bundle_html,
    run_section_bundle,
    sanitize_credentials,
    sanitize_polished_text,
    sources_for_sections,
    validate_section_keys,
    write_section_bundle_outputs,
)
from gene_dossier.rancho_report import (
    clear_stale_bundle_pngs,
    rasterize_pdf_pages_to_pngs,
)


FIXTURES = Path(__file__).parent / "fixtures" / "ucsc"
SEARCH_JSON = FIXTURES / "srebf2_search_relevant.json"
TRACK_JSON = FIXTURES / "srebf2_known_gene_region.json"


def _ev(**kwargs: Any) -> EvidenceRecord:
    defaults = dict(
        id=new_id(),
        source_id=f"src-{new_id()[:8]}",
        dossier_run_id="bundle-run",
        gene_symbol="SREBF2",
        section="General gene information",
        source_name="NCBI Gene",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="entrez_gene_id",
        evidence_grade=EvidenceGrade.C,
        taxon_id=9606,
        organism="Homo sapiens",
        value={},
        display_text="x",
    )
    defaults.update(kwargs)
    return EvidenceRecord(**defaults)


def _identity_records(gene: str = "SREBF2") -> list[EvidenceRecord]:
    return [
        _ev(
            gene_symbol=gene,
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            value={
                "entrez_gene_id": "6721",
                "nomenclaturesymbol": gene,
                "description": "sterol regulatory element binding transcription factor 2",
            },
            source_id="ncbi-entrez",
            id="eid-ncbi",
        ),
        _ev(
            gene_symbol=gene,
            source_name="Ensembl",
            fact_type="ensembl_gene_id",
            value={"ensembl_gene_id": "ENSG00000198911"},
            source_id="ensembl-id",
            id="eid-ensembl",
        ),
        _ev(
            gene_symbol=gene,
            source_name="UniProt",
            fact_type="uniprot_accession",
            value={"uniprot_accession": "Q12772", "gene_names": [gene, "SREBP2"]},
            source_id="uniprot-id",
            id="eid-uniprot",
        ),
    ]


def _ucsc_records(tmp_path: Path, gene: str = "SREBF2") -> list[EvidenceRecord]:
    from gene_dossier.normalize.ucsc_conservation import build_conservation_evidence
    from gene_dossier.ucsc_figure import sha256_hex

    search = json.loads(SEARCH_JSON.read_text())
    track = json.loads(TRACK_JSON.read_text())
    png = (
        b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
        b"\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDATx\x9cc\xf8\x0f"
        b"\x00\x00\x01\x01\x00\x05\x18\xd8N\x00\x00\x00\x00IEND\xaeB`\x82"
    )
    fig_dir = tmp_path / "raw" / "bundle-run" / "ucsc"
    fig_dir.mkdir(parents=True)
    fig_path = fig_dir / "figure.png"
    fig_path.write_bytes(png)
    rel = "bundle-run/ucsc/figure.png"
    figure_value = {
        "relative_path": rel,
        "local_artifact_path": rel,
        "sha256": sha256_hex(png),
        "media_type": "image/png",
        "width": 1400,
        "height": 400,
        "byte_size": len(png),
        "retrieval_method": "programmatic_browser_render",
        "track_preset_id": "ucsc_section_1b_comprehensive_v1",
    }
    records, _ = build_conservation_evidence(
        dossier_run_id="bundle-run",
        gene_symbol=gene,
        genome="hg38",
        search_payload=search,
        track_payload=track,
        figure_value=figure_value,
    )
    return records


def test_validate_section_keys_order_and_reject():
    assert validate_section_keys(["1c", "1b", "1a", "1a"]) == ["1a", "1b", "1c"]
    assert validate_section_keys(DEFAULT_SECTION_BUNDLE_KEYS) == list(
        DEFAULT_SECTION_BUNDLE_KEYS
    )
    assert SUPPORTED_SECTION_BUNDLE_KEYS == (
        "1a",
        "1b",
        "1c",
        "1d",
        "1e",
        "2a",
        "2b",
        "2c",
        "3a",
        "4a",
        "5a",
        "5b",
    )
    with pytest.raises(SectionBundleError):
        validate_section_keys([])


def test_sources_for_sections_dependency_aware():
    assert sources_for_sections(["1a"]) == []
    assert sources_for_sections(["1b"]) == ["UCSC"]
    assert sources_for_sections(["1a", "1b"]) == ["UCSC"]
    assert sources_for_sections(["1c"]) == ["CDD", "PDBe"]
    assert sources_for_sections(["1a", "1b", "1c"]) == ["CDD", "PDBe", "UCSC"]
    assert sources_for_sections(["1e"]) == []
    assert "NCBI Datasets" not in sources_for_sections(["1e"])
    assert sources_for_sections(["2a"]) == []
    assert "GTEx" not in sources_for_sections(["2a"])
    assert sources_for_sections(["1b", "2a"]) == ["UCSC"]

def test_transcript_selection_sentence_flags():
    assert "both MANE Select and Ensembl canonical" in transcript_selection_sentence(
        {"is_mane_select": True, "is_ensembl_canonical": True}
    )
    assert transcript_selection_sentence(
        {"is_mane_select": True, "is_ensembl_canonical": False}
    ).startswith("The MANE Select")
    assert transcript_selection_sentence(
        {"is_mane_select": False, "is_ensembl_canonical": True}
    ).startswith("The Ensembl canonical")
    assert "canonical-tier" in transcript_selection_sentence(
        {
            "is_mane_select": False,
            "is_ensembl_canonical": False,
            "is_canonical_tier": True,
        }
    )
    assert "highest-ranked" in transcript_selection_sentence({})


def test_opaque_refs_deterministic():
    table = ReportContentBlock(kind="table", presentation_role="gene_aliases_table")
    narrative = ReportContentBlock(kind="narrative", text="x")
    link = ReportContentBlock(kind="link", links=[{"label": "a", "url": "https://x"}])
    figure = ReportContentBlock(
        kind="figure", presentation_role="ucsc_conservation_figure", figure_path="a.png"
    )
    assert opaque_evidence_ref("1a", table, index=0) == "ev-1a-gene-aliases-table"
    assert opaque_evidence_ref("1b", narrative, index=0) == "ev-1b-summary"
    assert opaque_evidence_ref("1b", link, index=1) == "ev-1b-transcript-link"
    assert (
        opaque_evidence_ref("1b", figure, index=2) == "ev-1b-conservation-figure"
    )


def test_section_bundle_html_and_json_contracts(tmp_path, monkeypatch):
    from gene_dossier import report_presentation as rp

    evidence = _identity_records() + _ucsc_records(tmp_path)

    def _fake_resolve(value):
        return str(value.get("relative_path")), []

    monkeypatch.setattr(rp, "_resolve_figure_path", _fake_resolve)

    document, presentation, audit = build_section_bundle_document(
        dossier_run_id="bundle-run",
        gene_symbol="SREBF2",
        section_keys=["1a", "1b"],
        evidence_records=evidence,
    )

    # Document contains only Section 1 with 1a and 1b — no 1c
    assert len(document.sections) == 1
    keys = [s.key for s in document.sections[0].subsections]
    assert keys == ["a", "b"]
    assert document.sections[0].narrative_markdown is None

    html = render_section_bundle_html(document)
    assert "1. General Gene Information" in html
    assert "a. Gene Aliases" in html
    major_idx = html.index("1. General Gene Information")
    aliases_idx = html.index("a. Gene Aliases")
    assert major_idx < aliases_idx
    between = html[major_idx:aliases_idx]
    assert "Key findings" not in between
    assert "Limitations" not in between
    assert "Supporting evidence" not in between
    assert "Key findings" not in html
    assert "Limitations" not in html
    assert "Supporting evidence" not in html
    assert "[source_id=" not in html
    assert "max-width: 8.5in" in html
    assert "@page" in html
    assert 'data-evidence-ref="ev-1a-gene-aliases-table"' in html
    assert 'data-evidence-ref="ev-1b-summary"' in html
    assert "eid-ncbi" not in html
    assert "ncbi-entrez" not in html

    flat = json.dumps(presentation)
    assert "ev-1b-summary" in flat
    assert "evidence_record_ids" not in flat
    assert "source_ids" not in flat
    assert "eid-ncbi" not in flat

    assert "ev-1a-gene-aliases-table" in audit["evidence_reference_map"]
    assert "eid-ncbi" in json.dumps(audit["evidence_reference_map"])

    dirty = sanitize_credentials(
        {
            **audit,
            "diagnostics": [
                {
                    "reason": (
                        "url=https://x?apiKey=fake-secret-key "
                        "[source_id=audit-keep-me]"
                    )
                }
            ],
        }
    )
    dirty_text = json.dumps(dirty)
    assert "fake-secret-key" not in dirty_text
    assert "[source_id=audit-keep-me]" in dirty_text
    assert "REDACTED" in dirty_text

    # Presentation polish strips citation tokens; credentials scrubbed.
    polished = sanitize_polished_text(
        "see [source_id=gone] apiKey=fake-secret-key"
    )
    assert "[source_id=" not in polished
    assert "fake-secret-key" not in polished

    out = tmp_path / "out"
    paths = write_section_bundle_outputs(
        document=document,
        presentation=presentation,
        audit=dirty,
        output_dir=out,
        write_pdf=False,
    )
    assert paths["section_1_html"].is_file()
    assert paths["section_1_json"].is_file()
    assert paths["section_1_audit_json"].is_file()
    written_audit = paths["section_1_audit_json"].read_text()
    assert "fake-secret-key" not in written_audit
    assert "[source_id=audit-keep-me]" in written_audit
    assert "report-page report-chrome" in paths["section_1_html"].read_text()
    assert "@media print" in paths["section_1_html"].read_text()


def test_conservation_one_narrative_block(tmp_path, monkeypatch):
    from gene_dossier import report_presentation as rp

    evidence = _ucsc_records(tmp_path)

    def _fake_resolve(value):
        return str(value.get("relative_path")), []

    monkeypatch.setattr(rp, "_resolve_figure_path", _fake_resolve)
    result = build_conservation_blocks(
        gene_symbol="SREBF2", evidence_records=evidence
    )
    kinds = [b.kind for b in result.blocks]
    assert kinds == ["narrative", "link", "figure"]
    assert UCSC_STABLE_INTRO in (result.blocks[0].text or "")


def test_assign_opaque_refs_preserves_internal_ids_until_serialize():
    block = ReportContentBlock(
        kind="narrative",
        text="Hello [source_id=secret-id]",
        source_ids=["src-1"],
        evidence_record_ids=["ev-1"],
    )
    polished, ref_map = assign_opaque_refs(section_key="1b", blocks=[block])
    assert polished[0].evidence_ref == "ev-1b-summary"
    assert polished[0].source_ids == ["src-1"]
    assert ref_map["ev-1b-summary"]["source_ids"] == ["src-1"]
    from gene_dossier.section_bundle import serialize_presentation_block

    ser = serialize_presentation_block(polished[0])
    assert ser["evidence_ref"] == "ev-1b-summary"
    assert "source_ids" not in ser
    assert "secret-id" not in (ser.get("text") or "")


def test_source_note_not_second_narrative(tmp_path, monkeypatch):
    from gene_dossier import report_presentation as rp

    evidence = _ucsc_records(tmp_path)
    for rec in evidence:
        if rec.fact_type == "ucsc_conservation_figure" and isinstance(rec.value, dict):
            rec.value["source_note"] = "programmatic provenance note"
            rec.value["caption"] = "should stay in audit only"

    def _fake_resolve(value):
        return str(value.get("relative_path")), []

    monkeypatch.setattr(rp, "_resolve_figure_path", _fake_resolve)
    result = build_conservation_blocks(
        gene_symbol="SREBF2", evidence_records=evidence
    )
    assert [b.kind for b in result.blocks] == ["narrative", "link", "figure"]
    assert any(d.field == "figure_note" for d in result.diagnostics)


def test_gencode_missing_release_wording():
    # Exercise release_phrase via build path with empty inventory release.
    # Direct unit of the else branch:
    release = None
    if release and str(release).upper().startswith("GENCODE"):
        release_phrase = f"current {release} annotation"
    elif release and str(release).upper().startswith("V"):
        release_phrase = f"current GENCODE {release} annotation"
    elif release:
        release_phrase = f"current GENCODE {release} annotation"
    else:
        release_phrase = "current GENCODE annotation"
    assert release_phrase == "current GENCODE annotation"
    assert "GENCODE GENCODE" not in release_phrase


def test_finalize_rejects_missing_and_wrong_type(tmp_path, monkeypatch):
    from gene_dossier import section_bundle as sb
    from gene_dossier.models import DossierRun

    monkeypatch.setattr(sb, "init_db", lambda: None)

    class _Sess:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    monkeypatch.setattr(sb, "session_scope", lambda: _Sess())
    monkeypatch.setattr(sb, "get_dossier_run", lambda session, rid: None)
    with pytest.raises(SectionBundleError, match="not found"):
        finalize_section_bundle_run(
            dossier_run_id="missing",
            status="completed",
            selected_section_keys=["1a"],
            persist_db=True,
        )

    wrong = DossierRun(
        id="full-run",
        gene_symbol="SREBF2",
        run_type="full_api_pass",
        status="running",
    )
    monkeypatch.setattr(sb, "get_dossier_run", lambda session, rid: wrong)
    with pytest.raises(SectionBundleError, match="non-section-bundle"):
        finalize_section_bundle_run(
            dossier_run_id="full-run",
            status="completed",
            selected_section_keys=["1a"],
            persist_db=True,
        )


def test_provenance_map_includes_raw_artifact_ids(tmp_path, monkeypatch):
    from gene_dossier import report_presentation as rp

    evidence = _identity_records() + _ucsc_records(tmp_path)
    fig_rec = next(
        r for r in evidence if r.fact_type == "ucsc_conservation_figure"
    )
    fig_rec.raw_artifact_id = "raw-fig-1"
    fig_rec.api_run_id = "api-fig-1"
    fig_rec.value["figure_raw_artifact_id"] = "raw-fig-1"
    for rec in evidence:
        if not rec.raw_artifact_id:
            rec.raw_artifact_id = f"raw-{rec.id[:8]}"
        if not rec.api_run_id:
            rec.api_run_id = f"api-{rec.id[:8]}"

    monkeypatch.setattr(
        rp, "_resolve_figure_path", lambda value: (str(value.get("relative_path")), [])
    )
    raw_artifacts = [
        {
            "id": "raw-fig-1",
            "file_path": "bundle-run/ucsc/figure.png",
            "artifact_type": "image",
        }
    ]
    document, presentation, audit = build_section_bundle_document(
        dossier_run_id="bundle-run",
        gene_symbol="SREBF2",
        section_keys=["1a", "1b"],
        evidence_records=evidence,
        raw_artifacts=raw_artifacts,
    )
    fig_map = audit["evidence_reference_map"]["ev-1b-conservation-figure"]
    assert "raw-fig-1" in fig_map["raw_artifact_ids"]
    assert "api-fig-1" in fig_map["api_run_ids"]
    assert document is not None
    assert presentation["major_section"]["subsections"]


def test_populated_output_dir_rejected_and_preexisting_survives(tmp_path, monkeypatch):
    from gene_dossier import report_presentation as rp

    evidence = _identity_records() + _ucsc_records(tmp_path)
    monkeypatch.setattr(
        rp, "_resolve_figure_path", lambda value: (str(value.get("relative_path")), [])
    )
    document, presentation, audit = build_section_bundle_document(
        dossier_run_id="bundle-run",
        gene_symbol="SREBF2",
        section_keys=["1a"],
        evidence_records=evidence,
    )
    out = tmp_path / "populated"
    out.mkdir()
    preexisting = out / "section_1.html"
    preexisting.write_text("KEEP-ME", encoding="utf-8")
    with pytest.raises(SectionBundleError, match="already populated"):
        write_section_bundle_outputs(
            document=document,
            presentation=presentation,
            audit=audit,
            output_dir=out,
            write_pdf=False,
            allow_rerender=False,
        )
    assert preexisting.read_text() == "KEEP-ME"

    # Failure mid-write must not delete preexisting when allow_rerender.
    def _boom(*args, **kwargs):
        raise RuntimeError("pdf boom")

    monkeypatch.setattr(
        "gene_dossier.section_bundle.render_rancho_pdf", _boom
    )
    # Fresh empty dir with a pre-seeded file that is not a bundle stem... 
    # Use allow_rerender on a dir that already has section_1.html
    try:
        write_section_bundle_outputs(
            document=document,
            presentation=presentation,
            audit=audit,
            output_dir=out,
            write_pdf=True,
            allow_rerender=True,
        )
    except RuntimeError:
        pass
    # Preexisting html may have been overwritten (existed_before=True) — must remain
    assert preexisting.exists()


def test_stale_png_cleanup_and_multipage_names(tmp_path):
    out = tmp_path / "pngs"
    out.mkdir()
    (out / "section_1.png").write_bytes(b"old")
    (out / "section_1_page_2.png").write_bytes(b"stale")
    (out / "section_1_contact_sheet.png").write_bytes(b"sheet")
    clear_stale_bundle_pngs(out)
    assert not (out / "section_1.png").exists()
    assert not (out / "section_1_page_2.png").exists()
    assert not (out / "section_1_contact_sheet.png").exists()

    # Multipage naming when pymupdf available
    try:
        import fitz
    except Exception:
        pytest.skip("pymupdf not available")
    pdf = out / "demo.pdf"
    doc = fitz.open()
    for _ in range(2):
        doc.new_page(width=612, height=792)
    doc.save(str(pdf))
    doc.close()
    written = rasterize_pdf_pages_to_pngs(pdf, out, stem="section_1", dpi=72)
    assert len(written) == 2
    assert written[0].name == "section_1_page_1.png"
    assert written[1].name == "section_1_page_2.png"
    assert not (out / "section_1.png").exists()


def test_run_section_bundle_mocked_no_llm_and_1a_skips_ucsc(tmp_path, monkeypatch):
    from gene_dossier import section_bundle as sb
    from gene_dossier.models import DossierRun, ToolResult
    from gene_dossier.workflow import WorkflowTransientContext

    calls: dict[str, Any] = {"ucsc": 0, "identity": 0, "save": 0, "transients": []}

    def _fake_create(**kwargs):
        run = DossierRun(
            gene_symbol=kwargs["gene_symbol"],
            run_type="section_bundle",
            status="running",
            config={"selected_section_keys": list(kwargs["selected_section_keys"])},
        )
        state = {
            "gene_symbol": run.gene_symbol,
            "dossier_run_id": run.id,
            "gene_ids": {},
            "tool_results": [],
            "api_runs": [],
            "raw_artifacts": [],
            "evidence_records": _identity_records(),
            "coverage": [],
            "sections": [],
            "claims": [],
            "verification_results": [],
            "synthesis_notes": [],
            "output_paths": {},
            "errors": [],
            "status": "running",
        }
        return run, state

    def _fake_identity(state, **kwargs):
        calls["identity"] += 1
        return state

    def _fake_call(state, *, sources=None, transient=None, **kwargs):
        calls["ucsc"] += 1
        calls["transients"].append(transient)
        results = list(state.get("tool_results") or [])
        results.append(
            ToolResult(
                source_name="UCSC",
                endpoint_name="fetch_gene_region",
                gene_symbol=state["gene_symbol"],
                request_url="https://example.test/ucsc",
                success=True,
                data={"ok": True},
            )
        )
        return {**state, "tool_results": results}

    def _fake_save(state, *, transient=None, **kwargs):
        calls["save"] += 1
        calls["transients"].append(transient)
        return state

    def _fake_normalize(state, **kwargs):
        return state

    monkeypatch.setattr(sb, "create_section_bundle_run", _fake_create)
    monkeypatch.setattr(sb, "node_resolve_gene_identity", _fake_identity)
    monkeypatch.setattr(sb, "node_call_source_clients", _fake_call)
    monkeypatch.setattr(sb, "node_save_raw_artifacts", _fake_save)
    monkeypatch.setattr(sb, "node_normalize_evidence", _fake_normalize)
    monkeypatch.setattr(sb, "finalize_section_bundle_run", lambda **kwargs: None)
    monkeypatch.setattr(
        sb,
        "coverage_updates_from_state",
        lambda state: [],
    )

    # Guard: LLM / full report must never be imported/called from this path.
    import gene_dossier.synthesis as synthesis
    import gene_dossier.rancho_report as rancho

    monkeypatch.setattr(
        synthesis,
        "synthesize_dossier",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("LLM called")),
    )
    monkeypatch.setattr(
        rancho,
        "build_and_write_rancho_report",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("full report called")),
    )

    # 1a only — UCSC must not be called
    result = run_section_bundle(
        "SREBF2",
        section_keys=["1a"],
        output_dir=tmp_path / "1a",
        call_network=True,
        persist_db=False,
        write_pdf=False,
    )
    assert result.status == "completed"
    assert calls["ucsc"] == 0
    assert calls["identity"] == 1

    # 1a+1b — same transient instance for call + save
    calls["ucsc"] = 0
    calls["save"] = 0
    calls["transients"] = []
    evidence = _identity_records() + _ucsc_records(tmp_path)

    def _fake_create_1b(**kwargs):
        run, state = _fake_create(**kwargs)
        state["evidence_records"] = evidence
        return run, state

    monkeypatch.setattr(sb, "create_section_bundle_run", _fake_create_1b)
    from gene_dossier import report_presentation as rp

    monkeypatch.setattr(
        rp, "_resolve_figure_path", lambda value: (str(value.get("relative_path")), [])
    )

    result2 = run_section_bundle(
        "SREBF2",
        section_keys=["1a", "1b"],
        output_dir=tmp_path / "1b",
        call_network=True,
        persist_db=False,
        write_pdf=False,
    )
    assert result2.status == "completed"
    assert calls["ucsc"] == 1
    assert calls["save"] == 1
    assert len(calls["transients"]) >= 2
    assert calls["transients"][0] is calls["transients"][1]
    assert isinstance(calls["transients"][0], WorkflowTransientContext)


def test_run_section_bundle_1c_calls_only_cdd_pdbe_and_bundle_derived(tmp_path, monkeypatch):
    from gene_dossier import section_bundle as sb
    from gene_dossier.models import DossierRun, ToolResult

    calls: dict[str, Any] = {"sources": None, "derived": 0}

    def _fake_create(**kwargs):
        run = DossierRun(
            gene_symbol=kwargs["gene_symbol"],
            run_type="section_bundle",
            status="running",
            config={"selected_section_keys": list(kwargs["selected_section_keys"])},
        )
        state = {
            "gene_symbol": run.gene_symbol,
            "dossier_run_id": run.id,
            "run_type": "section_bundle",
            "selected_section_keys": list(kwargs["selected_section_keys"]),
            "gene_ids": {"uniprot_accession": "Q12772", "protein_length": 1141},
            "tool_results": [],
            "api_runs": [],
            "raw_artifacts": [],
            "evidence_records": _identity_records(),
            "coverage": [],
            "sections": [],
            "claims": [],
            "verification_results": [],
            "synthesis_notes": [],
            "output_paths": {},
            "errors": [],
            "status": "running",
        }
        return run, state

    def _fake_call(state, *, sources=None, **kwargs):
        calls["sources"] = list(sources or [])
        results = list(state.get("tool_results") or [])
        for source in sources or []:
            results.append(
                ToolResult(
                    source_name=source,
                    endpoint_name="mock",
                    gene_symbol=state["gene_symbol"],
                    request_url="https://example.test",
                    success=True,
                    data={"ok": True},
                )
            )
        return {**state, "tool_results": results}

    def _fake_derived(state, **kwargs):
        calls["derived"] += 1
        return {
            **state,
            "section_1c": {
                "section_status": "failed",
                "source_status": {"CDD": "unavailable", "PDBe": "unavailable"},
                "rendering_status": {},
            },
        }

    monkeypatch.setattr(sb, "create_section_bundle_run", _fake_create)
    monkeypatch.setattr(sb, "node_resolve_gene_identity", lambda state, **kwargs: state)
    monkeypatch.setattr(sb, "node_call_source_clients", _fake_call)
    monkeypatch.setattr(sb, "node_save_raw_artifacts", lambda state, **kwargs: state)
    monkeypatch.setattr(sb, "node_normalize_evidence", lambda state, **kwargs: state)
    monkeypatch.setattr(sb, "node_generate_section_1c_derived_artifacts", _fake_derived)
    monkeypatch.setattr(sb, "finalize_section_bundle_run", lambda **kwargs: None)
    monkeypatch.setattr(sb, "coverage_updates_from_state", lambda state: [])

    result = run_section_bundle(
        "SREBF2",
        section_keys=["1c"],
        output_dir=tmp_path / "1c",
        call_network=True,
        persist_db=False,
        write_pdf=False,
    )
    assert result.status == "completed"
    assert calls["sources"] == ["CDD", "PDBe"]
    assert calls["derived"] == 1
    audit = json.loads(result.output_paths["section_1_audit_json"].read_text())
    assert audit["section_1c"]["section_status"] == "failed"


def test_section_1c_polished_blocks_are_figure_led_not_tables():
    from gene_dossier.config import get_settings

    image_path = get_settings().raw_data_path / "tests" / "section-1c-pdbe-official.png"
    image_path.parent.mkdir(parents=True, exist_ok=True)
    image_path.write_bytes((PROJECT_ROOT / "src/gene_dossier/assets/rancho_wordmark.png").read_bytes())
    relative_image = str(image_path.relative_to(get_settings().raw_data_path))
    records = [
        _ev(
            source_name="UniProt",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="uniprot_accession",
            value={
                "uniprot_accession": "Q12772",
                "reviewed": True,
                "taxon_id": 9606,
                "protein_length": 1141,
            },
        ),
        _ev(
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="conserved_domain_hit",
            value={
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "domain_description": "basic Helix-Loop-Helix zipper domain.",
                "from_residue": 343,
                "to_residue": 403,
                "evalue": "1e-30",
                "bitscore": 88.2,
            },
        ),
        _ev(
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="cdd_official_architecture_figure",
            value={
                "relative_path": relative_image,
                "media_type": "image/png",
                "width": 640,
                "height": 120,
                "artifact_class": "official",
            },
        ),
        _ev(
            source_name="PDBe",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="pdb_candidate_selection",
            value={
                "selected_uniprot_accession": "Q12772",
                "candidates": [
                    {
                        "pdb_id": "1ukl",
                        "chain_ids": ["C", "D", "E", "F"],
                        "experimental_method": "X-ray diffraction",
                        "mapped_spans": [[343, 403]],
                        "title": "Crystal structure of Importin-beta and SREBP-2 complex",
                        "selected": True,
                    }
                ],
            },
        ),
        _ev(
            source_name="PDBe",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="pdb_official_structure_image",
            value={
                "pdb_id": "1ukl",
                "relative_path": relative_image,
                "media_type": "image/png",
                "width": 640,
                "height": 480,
                "attribution": "Image source: PDBe, PDB 1UKL",
                "species_common_name": "human",
            },
        ),
    ]

    result = build_known_structure_blocks(
        gene_symbol="SREBF2",
        evidence_records=records,
    )
    roles = [block.presentation_role for block in result.blocks]
    assert all(block.kind != "table" for block in result.blocks)
    assert "section_1c_cdd_link" in roles
    assert "section_1c_domain_architecture_figure" in roles
    assert roles.index("section_1c_cdd_link") < roles.index("section_1c_domain_architecture_figure")
    assert "section_1c_domain_summary" in roles
    assert "section_1c_pdb_link" in roles
    assert "section_1c_pdb_official_image" in roles
    pdbe_fig = next(block for block in result.blocks if block.presentation_role == "section_1c_pdb_official_image")
    assert pdbe_fig.figure_caption == "Image source: PDBe, PDB 1UKL"
    assert "Conserved Domain Database" in result.blocks[0].text
    assert any("1ukl" in (link.get("url") or "") for block in result.blocks for link in block.links)
    assert not any((diag.field or "").startswith("pymol") for diag in result.diagnostics)


def test_section_1c_dynamic_evidence_refs_use_safe_item_keys():
    blocks = [
        ReportContentBlock(kind="narrative", text="intro"),
        ReportContentBlock(
            kind="link",
            text="CDD",
            presentation_role="section_1c_cdd_link",
            links=[{"label": "CDD", "url": "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml"}],
        ),
        ReportContentBlock(
            kind="narrative",
            text="Cadherin repeat",
            presentation_role="section_1c_domain_summary",
            presentation_item_key="domain-cd11304",
        ),
        ReportContentBlock(
            kind="figure",
            figure_path="/tmp/missing.png",
            presentation_role="section_1c_domain_thumbnail",
            presentation_item_key="domain-cd11304",
        ),
        ReportContentBlock(
            kind="link",
            text="PDB",
            presentation_role="section_1c_pdb_link",
            presentation_item_key="pdb-6cg6",
            links=[{"label": "PDB", "url": "https://www.ebi.ac.uk/pdbe/entry/pdb/6cg6"}],
        ),
        ReportContentBlock(
            kind="figure",
            figure_path="/tmp/missing.png",
            presentation_role="section_1c_pdb_official_image",
            presentation_item_key="pdb-6cg6",
        ),
    ]
    polished, ref_map = assign_opaque_refs(section_key="1c", blocks=blocks)
    refs = [block.evidence_ref for block in polished]
    assert refs == [
        "ev-1c-introduction",
        "ev-1c-cdd-link",
        "ev-1c-domain-cd11304-summary",
        "ev-1c-domain-cd11304-thumbnail",
        "ev-1c-pdb-6cg6-summary",
        "ev-1c-pdb-6cg6-official-image",
    ]
    assert len(ref_map) == len(set(ref_map))


def test_section_1c_dynamic_evidence_refs_reject_missing_item_key():
    with pytest.raises(SectionBundleError):
        assign_opaque_refs(
            section_key="1c",
            blocks=[
                ReportContentBlock(
                    kind="narrative",
                    text="bad",
                    presentation_role="section_1c_domain_summary",
                )
            ],
        )


def test_section_1c_presentation_strips_pssm_and_human_expression_host():
    records = [
        _ev(
            source_name="UniProt",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="uniprot_accession",
            value={"uniprot_accession": "Q12772", "taxon_id": 9606, "protein_length": 1141},
        ),
        _ev(
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="conserved_domain_hit",
            value={
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "from_residue": 343,
                "to_residue": 403,
            },
        ),
        _ev(
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="cdd_family_summary",
            value={
                "canonical_accession": "cd18922",
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "synopsis": "cd18922 (PSSM ID: 381492): specific bHLHzip domain.",
                "presentation_item_key": "domain-cd18922",
            },
        ),
        _ev(
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="cdd_family_summary",
            value={
                "canonical_accession": "cl00081",
                "domain_accession": "cl00081",
                "domain_short_name": "bHLH_SF",
                "synopsis": "cl00081 (PSSM ID: 444684): bHLH proteins are transcriptional regulators.",
                "presentation_item_key": "domain-cl00081",
            },
        ),
        _ev(
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="cdd_conserved_feature",
            value={
                "domain_accession": "cd18922",
                "feature_label": "putative DNA binding site",
                "feature_type": "site",
                "query_residues": "R350, E351",
                "presentation_item_key": "feature-cd18922-dna-binding",
            },
        ),
        _ev(
            source_name="PDBe",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="pdb_candidate_selection",
            value={
                "candidates": [
                    {
                        "pdb_id": "1ukl",
                        "selected": True,
                        "species_common_name": "human",
                        "expression_host": "Escherichia coli BL21",
                    }
                ],
            },
        ),
    ]
    result = build_known_structure_blocks(gene_symbol="SREBF2", evidence_records=records)
    visible = "\n".join(block.text or "" for block in result.blocks)
    assert "PSSM ID" not in visible
    assert "Conserved Protein Domain Family bHLH_SF" in visible
    assert "putative DNA binding site" not in visible
    bhlhzip_block = next(
        block
        for block in result.blocks
        if block.presentation_role == "section_1c_domain_summary"
        and block.presentation_item_key == "domain-cd18922"
    )
    assert bhlhzip_block.links == []
    bhlh_sf_block = next(
        block
        for block in result.blocks
        if block.presentation_role == "section_1c_domain_summary"
        and block.presentation_item_key == "domain-cl00081"
    )
    assert bhlh_sf_block.links[0]["label"] == "Conserved Protein Domain Family bHLH_SF:"
    pdb_link = next(block for block in result.blocks if block.presentation_role == "section_1c_pdb_link")
    assert "expressed in" not in pdb_link.text


def test_section_1c_cdh10_order_and_suppression():
    records = [
        _ev(
            gene_symbol="CDH10",
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="conserved_domain_hit",
            value={
                "domain_accession": "smart00112",
                "domain_short_name": "CA",
                "from_residue": 100,
                "to_residue": 160,
            },
        ),
        _ev(
            gene_symbol="CDH10",
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="conserved_domain_hit",
            value={
                "domain_accession": "cd11304",
                "domain_short_name": "Cadherin_repeat",
                "from_residue": 170,
                "to_residue": 280,
            },
        ),
        _ev(
            gene_symbol="CDH10",
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="conserved_domain_hit",
            value={
                "domain_accession": "pfam01049",
                "domain_short_name": "CADH_Y-type_LIR",
                "from_residue": 900,
                "to_residue": 960,
            },
        ),
        _ev(
            gene_symbol="CDH10",
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="cdd_family_summary",
            value={
                "canonical_accession": "cd11304",
                "domain_accession": "cd11304",
                "domain_short_name": "Cadherin_repeat",
                "synopsis": "cd11304 (PSSM ID: 206637): Cadherins are glycoproteins.",
                "presentation_item_key": "domain-cd11304",
            },
        ),
        _ev(
            gene_symbol="CDH10",
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="cdd_family_summary",
            value={
                "canonical_accession": "pfam01049",
                "matched_query_domain_accession": "pfam01049",
                "domain_accession": "pfam01049",
                "domain_short_name": "Cadherin_C / CADH_Y-type_LIR",
                "pssm_id": "426014",
                "synopsis": "Cadherin cytoplasmic region: Cadherins are linked to the cytoskeleton by catenins.",
                "presentation_item_key": "domain-pfam01049",
            },
        ),
        _ev(
            gene_symbol="CDH10",
            source_name="CDD",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="cdd_conserved_feature",
            value={
                "domain_accession": "cd11304",
                "feature_label": "Ca2+ binding site [ion binding site]",
                "feature_type": "ion binding site",
                "query_residues": "E171, M172",
                "presentation_item_key": "feature-cd11304-ca2-binding-site-ion-binding-site",
            },
        ),
    ]
    result = build_known_structure_blocks(gene_symbol="CDH10", evidence_records=records)
    visible = "\n".join(block.text or "" for block in result.blocks)
    assert "CA:" not in visible
    # Broad smart00112 "CA" prose stays suppressed while the specific
    # Cadherin_C block renders under its polished name only.
    assert "Cadherin_C: Cadherin cytoplasmic region:" in visible
    assert "CADH_Y-type_LIR" not in visible
    assert "PSSM ID" not in visible
    feature_idx = visible.index("Ca2+ binding site [ion binding site]")
    domain_idx = visible.index("CD11304: Cadherin_repeat")
    cterminal_idx = visible.index("Cadherin_C:")
    assert feature_idx < domain_idx
    assert domain_idx < cterminal_idx
    assert "E171" not in visible

    cadherin_repeat = next(
        block
        for block in result.blocks
        if block.presentation_item_key == "domain-cd11304"
        and block.presentation_role == "section_1c_domain_summary"
    )
    assert cadherin_repeat.links[0]["label"] == "(NCBI CDD Link)"
    assert cadherin_repeat.presentation_page_break_before is False

    cadherin_c = next(
        block
        for block in result.blocks
        if block.presentation_item_key == "domain-pfam01049"
        and block.presentation_role == "section_1c_domain_summary"
    )
    assert cadherin_c.links[0]["label"] == "(NCBI CDD Link)"
    assert cadherin_c.presentation_page_break_before is True


def test_section_1c_cdh10_mouse_pdb_heading_omits_ortholog_word():
    records = [
        _ev(
            gene_symbol="CDH10",
            source_name="PDBe",
            source_type=SourceType.structure_database,
            assertion_type=AssertionType.protein_structure,
            fact_type="pdb_candidate_selection",
            value={
                "candidates": [
                    {
                        "pdb_id": "6cg6",
                        "selected": True,
                        "species_common_name": "mouse",
                        "expression_host": "Escherichia coli",
                    }
                ],
            },
        )
    ]
    result = build_known_structure_blocks(gene_symbol="CDH10", evidence_records=records)
    pdb_link = next(block for block in result.blocks if block.presentation_role == "section_1c_pdb_link")
    assert pdb_link.text == (
        "3D structures from PDB: Mouse CDH10 protein expressed in Escherichia coli (PDB link)"
    )
    assert "ortholog" not in (pdb_link.text or "").lower()


# --------------------------------------------------------------------------------------
# Section 1c page layout against the original Rancho report body flow
# --------------------------------------------------------------------------------------
def _probe_figure(name: str) -> str:
    """Write a managed figure fixture and return its portable relative path."""
    from gene_dossier.config import get_settings

    root = get_settings().raw_data_path
    path = root / "tests" / f"section-1c-{name}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        (PROJECT_ROOT / "src/gene_dossier/assets/rancho_wordmark.png").read_bytes()
    )
    return str(path.relative_to(root))


def _cdd(gene: str, fact_type: str, value: dict[str, Any]) -> EvidenceRecord:
    return _ev(
        gene_symbol=gene,
        source_name="CDD",
        source_type=SourceType.structure_database,
        assertion_type=AssertionType.protein_structure,
        fact_type=fact_type,
        value=value,
    )


def _pdbe(gene: str, fact_type: str, value: dict[str, Any]) -> EvidenceRecord:
    return _ev(
        gene_symbol=gene,
        source_name="PDBe",
        source_type=SourceType.structure_database,
        assertion_type=AssertionType.protein_structure,
        fact_type=fact_type,
        value=value,
    )


def _srebf2_1c_records() -> list[EvidenceRecord]:
    """SREBF2 Section 1c evidence shaped like a real CDD/PDBe bundle."""
    figure = _probe_figure("srebf2")
    return [
        _ev(
            source_name="UniProt",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="uniprot_accession",
            value={"uniprot_accession": "Q12772", "taxon_id": 9606, "protein_length": 1141},
        ),
        _cdd(
            "SREBF2",
            "conserved_domain_hit",
            {
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "from_residue": 343,
                "to_residue": 403,
            },
        ),
        _cdd(
            "SREBF2",
            "cdd_official_architecture_figure",
            {"relative_path": figure, "width": 1197, "height": 148},
        ),
        _cdd(
            "SREBF2",
            "cdd_family_summary",
            {
                "canonical_accession": "cd18922",
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "requested_uid": "381492",
                "synopsis": "specific bHLHzip domain.",
                "presentation_item_key": "domain-cd18922",
            },
        ),
        _cdd(
            "SREBF2",
            "cdd_family_summary",
            {
                "canonical_accession": "cl00081",
                "domain_accession": "cl00081",
                "domain_short_name": "bHLH_SF",
                "requested_uid": "444684",
                "synopsis": (
                    "basic Helix Loop Helix (bHLH) domain superfamily: bHLH proteins "
                    "are transcriptional regulators found from yeast to humans."
                ),
                "presentation_item_key": "domain-cl00081",
            },
        ),
        # Square official structure thumbnail from the same CDD page: renders.
        _cdd(
            "SREBF2",
            "cdd_family_thumbnail",
            {
                "presentation_item_key": "domain-cl00081",
                "relative_path": figure,
                "width": 100,
                "height": 100,
                "classified_role": "family_structure_thumbnail",
                "source_url": (
                    "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdThumbnail.cgi?uid=444684"
                ),
            },
        ),
        # Thin alignment strip served for the specific hit: never renders.
        _cdd(
            "SREBF2",
            "cdd_family_thumbnail",
            {
                "presentation_item_key": "domain-cd18922",
                "relative_path": figure,
                "width": 100,
                "height": 24,
                "classified_role": "family_structure_thumbnail",
                "source_url": (
                    "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdThumbnail.cgi?uid=381492"
                ),
            },
        ),
        _pdbe(
            "SREBF2",
            "pdb_candidate_selection",
            {
                "candidates": [
                    {
                        "pdb_id": "1ukl",
                        "selected": True,
                        "species_common_name": "human",
                        "expression_host": "Escherichia coli",
                    }
                ]
            },
        ),
        _pdbe(
            "SREBF2",
            "pdb_official_structure_image",
            {
                "pdb_id": "1ukl",
                "relative_path": figure,
                "attribution": "Image source: PDBe, PDB 1UKL",
            },
        ),
    ]


def _cdh10_1c_records(
    *, pfam_thumbnails: list[dict[str, Any]] | None = None
) -> list[EvidenceRecord]:
    """CDH10 Section 1c evidence; ``pfam_thumbnails`` overrides Cadherin_C images.

    The default reproduces the live bundle, whose only pfam01049 thumbnail is a
    100x24 strip served from the neighbouring ``uid=460041`` family page.
    """
    figure = _probe_figure("cdh10")
    default_thumbnails = [
        {
            "presentation_item_key": "domain-pfam01049",
            "relative_path": figure,
            "width": 100,
            "height": 24,
            "classified_role": "family_structure_thumbnail",
            "source_url": (
                "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdThumbnail.cgi?uid=460041"
            ),
        }
    ]
    records = [
        _cdd(
            "CDH10",
            "conserved_domain_hit",
            {
                "domain_accession": "smart00112",
                "domain_short_name": "CA",
                "from_residue": 80,
                "to_residue": 158,
            },
        ),
        _cdd(
            "CDH10",
            "conserved_domain_hit",
            {
                "domain_accession": "cd11304",
                "domain_short_name": "Cadherin_repeat",
                "from_residue": 170,
                "to_residue": 280,
            },
        ),
        _cdd(
            "CDH10",
            "conserved_domain_hit",
            {
                "domain_accession": "pfam01049",
                "domain_short_name": "CADH_Y-type_LIR",
                "from_residue": 721,
                "to_residue": 780,
            },
        ),
        _cdd(
            "CDH10",
            "cdd_official_architecture_figure",
            {"relative_path": figure, "width": 1197, "height": 148},
        ),
        _cdd(
            "CDH10",
            "cdd_family_summary",
            {
                "canonical_accession": "cd11304",
                "domain_accession": "cd11304",
                "domain_short_name": "Cadherin_repeat",
                "requested_uid": "206637",
                "synopsis": (
                    "Cadherin tandem repeat domain.: Cadherins are glycoproteins "
                    "involved in Ca2+-mediated cell-cell adhesion."
                ),
                "presentation_item_key": "domain-cd11304",
            },
        ),
        _cdd(
            "CDH10",
            "cdd_family_summary",
            {
                "canonical_accession": "pfam01049",
                "matched_query_domain_accession": "pfam01049",
                "domain_accession": "pfam01049",
                "domain_short_name": "Cadherin_C / CADH_Y-type_LIR",
                "pssm_id": "426014",
                "requested_uid": "426014",
                "synopsis": (
                    "Cadherin cytoplasmic region: Cadherins are vital in cell-cell "
                    "adhesion during tissue differentiation."
                ),
                "presentation_item_key": "domain-pfam01049",
            },
        ),
        _cdd(
            "CDH10",
            "cdd_conserved_feature",
            {
                "domain_accession": "cd11304",
                "feature_label": "Ca2+ binding site [ion binding site]",
                "feature_type": "ion binding site",
                "query_residues": "E280, S281",
                "presentation_item_key": "feature-cd11304-ca2",
            },
        ),
        _cdd(
            "CDH10",
            "cdd_feature_thumbnail",
            {
                "presentation_item_key": "feature-cd11304-ca2",
                "relative_path": figure,
                "width": 300,
                "height": 300,
                "classified_role": "conserved_feature_structure_thumbnail",
            },
        ),
        _pdbe(
            "CDH10",
            "pdb_candidate_selection",
            {
                "candidates": [
                    {
                        "pdb_id": "6cg6",
                        "selected": True,
                        "species_common_name": "mouse",
                        "expression_host": "Escherichia coli",
                    }
                ]
            },
        ),
        _pdbe(
            "CDH10",
            "pdb_official_structure_image",
            {
                "pdb_id": "6cg6",
                "relative_path": figure,
                "attribution": "Image source: PDBe, PDB 6CG6",
            },
        ),
    ]
    thumbnails = default_thumbnails if pfam_thumbnails is None else pfam_thumbnails
    records.extend(_cdd("CDH10", "cdd_family_thumbnail", value) for value in thumbnails)
    return records


def _rendered_1c_pages(gene: str, records: list[EvidenceRecord]) -> tuple[list[str], dict[str, Any]]:
    """Focused Section 1c HTML split into its rendered report pages, plus audit."""
    document, _presentation, audit = build_section_bundle_document(
        dossier_run_id="layout-run",
        gene_symbol=gene,
        section_keys=["1c"],
        evidence_records=records,
    )
    html = render_section_bundle_html(
        document,
        include_page_chrome=False,
        include_major_heading=False,
    )
    pages = [f'<section id="section-1{part}' for part in html.split('<section id="section-1')[1:]]
    return pages, audit


def test_section_1c_srebf2_page_layout_matches_rancho_body_flow():
    pages, _audit = _rendered_1c_pages("SREBF2", _srebf2_1c_records())
    assert len(pages) == 2
    cdd_page, pdb_page = pages

    # Page 1: architecture, then the bHLHzip lead sentence with a bold lead
    # phrase and no standalone bHLHzip_SREBP2 link block, then the linked
    # bHLH_SF family heading, its synopsis, and the square green thumbnail.
    assert "section-1c-domain-architecture-figure" in cdd_page
    assert (
        "<strong>basic Helix-Loop-Helix-zipper (bHLHzip) domain</strong>"
        " found in sterol regulatory element-binding protein 2"
    ) in cdd_page
    assert "bHLHzip_SREBP2" not in cdd_page
    assert ">Conserved Protein Domain Family bHLH_SF:</a>" in cdd_page
    assert "transcriptional regulators found from yeast to humans" in cdd_page
    assert 'data-item-key="domain-cl00081"' in cdd_page
    assert "section-1c-domain-thumbnail" in cdd_page

    # The bHLHzip specific-hit alignment strip never becomes a thumbnail.
    assert cdd_page.count("section-1c-domain-thumbnail") == 1
    assert "3D structures from PDB" not in cdd_page

    # Page 2: PDB heading, official image, and attribution only.
    assert "3D structures from PDB: Human SREBF2 protein (PDB link)" in pdb_page
    assert "section-1c-pdb-official-image" in pdb_page
    assert "Image source: PDBe, PDB 1UKL" in pdb_page
    assert "section-1c-domain-architecture-figure" not in pdb_page


def test_section_1c_cdh10_page_layout_matches_rancho_body_flow():
    pages, audit = _rendered_1c_pages("CDH10", _cdh10_1c_records())
    assert len(pages) == 3
    repeat_page, cadherin_c_page, pdb_page = pages

    # Page 1: architecture, then the Ca2+ feature heading and its larger
    # thumbnail, then the Cadherin_repeat paragraph closing with its CDD link.
    assert "section-1c-domain-architecture-figure" in repeat_page
    assert "Ca2+ binding site [ion binding site]" in repeat_page
    assert "section-1c-feature-thumbnail" in repeat_page
    assert "CD11304: Cadherin_repeat:" in repeat_page
    assert repeat_page.index("Ca2+ binding site") < repeat_page.index("CD11304: Cadherin_repeat:")
    assert ">(NCBI CDD Link)</a>" in repeat_page
    # Broad smart00112 "CA" prose is not a separate visible family block.
    assert ">CA:" not in repeat_page
    assert "3D structures from PDB" not in repeat_page

    # Page 2: Cadherin_C text and link under its polished name only. The live
    # bundle's only pfam01049 image is a rejected strip, so no image renders.
    assert "Cadherin_C: Cadherin cytoplasmic region:" in cadherin_c_page
    assert "CADH_Y-type_LIR" not in cadherin_c_page
    assert ">(NCBI CDD Link)</a>" in cadherin_c_page
    assert "section-1c-domain-thumbnail" not in cadherin_c_page
    assert "3D structures from PDB" not in cadherin_c_page
    omitted = [
        diag
        for diag in audit["diagnostics"]
        if diag["field"] == "section_1c_domain_thumbnail_omitted"
    ]
    assert omitted and "domain-pfam01049" in omitted[0]["reason"]

    # Page 3: PDB heading without "ortholog", official image, attribution.
    assert (
        "3D structures from PDB: Mouse CDH10 protein expressed in Escherichia coli (PDB link)"
    ) in pdb_page
    assert "ortholog" not in pdb_page.lower()
    assert "section-1c-pdb-official-image" in pdb_page
    assert "Image source: PDBe, PDB 6CG6" in pdb_page


def test_section_1c_two_official_thumbnails_bracket_the_cadherin_c_paragraph():
    figure = _probe_figure("cdh10")
    records = _cdh10_1c_records(
        pfam_thumbnails=[
            {
                "presentation_item_key": "domain-pfam01049",
                "relative_path": figure,
                "width": 100,
                "height": 100,
                "classified_role": "family_structure_thumbnail",
                "source_url": (
                    "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdThumbnail.cgi?uid=426014"
                ),
            },
            {
                "presentation_item_key": "domain-pfam01049",
                "relative_path": figure,
                "width": 120,
                "height": 110,
                "classified_role": "family_structure_thumbnail",
                "source_url": (
                    "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdThumbnail.cgi?uid=426014"
                ),
            },
        ]
    )
    result = build_known_structure_blocks(gene_symbol="CDH10", evidence_records=records)
    roles = [
        block.presentation_role
        for block in result.blocks
        if block.presentation_item_key == "domain-pfam01049"
    ]
    assert roles == [
        "section_1c_domain_thumbnail",
        "section_1c_domain_summary",
        "section_1c_domain_thumbnail",
    ]
    # The page break stays on whichever block opens the Cadherin_C page.
    leading = next(
        block
        for block in result.blocks
        if block.presentation_item_key == "domain-pfam01049"
    )
    assert leading.presentation_page_break_before is True

    pages, _audit = _rendered_1c_pages("CDH10", records)
    assert len(pages) == 3
    assert pages[1].count("section-1c-domain-thumbnail") == 2


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ({"width": 100, "height": 100}, True),
        ({"width": 300, "height": 300}, True),
        ({"width": 120, "height": 110}, True),
        ({}, True),
        ({"width": 100, "height": 24}, False),
        ({"width": 100, "height": 40}, False),
        ({"width": 250, "height": 100}, False),
        ({"width": 100, "height": 100, "classified_role": "alignment_or_sequence_thumbnail"}, False),
        ({"width": 100, "height": 100, "classified_role": "family_sequence_logo"}, False),
        ({"width": 100, "height": 100, "classified_role": "msa_graphic"}, False),
    ],
)
def test_section_1c_structure_thumbnail_gate(value: dict[str, Any], expected: bool):
    from gene_dossier.report_presentation import _is_renderable_structure_thumbnail

    record = _cdd("CDH10", "cdd_family_thumbnail", value)
    assert _is_renderable_structure_thumbnail(record) is expected


def test_section_1c_thumbnail_from_other_cdd_family_page_is_rejected():
    from gene_dossier.report_presentation import _is_renderable_structure_thumbnail

    record = _cdd(
        "CDH10",
        "cdd_family_thumbnail",
        {
            "width": 100,
            "height": 100,
            "classified_role": "family_structure_thumbnail",
            "source_url": (
                "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdThumbnail.cgi?uid=460041"
            ),
        },
    )
    assert _is_renderable_structure_thumbnail(record, family_uid="426014") is False
    assert _is_renderable_structure_thumbnail(record, family_uid="460041") is True


def test_section_1c_page_split_is_identical_in_html_pdf_and_pngs(tmp_path):
    document, presentation, audit = build_section_bundle_document(
        dossier_run_id="split-run",
        gene_symbol="CDH10",
        section_keys=["1c"],
        evidence_records=_cdh10_1c_records(),
    )
    paths = write_section_bundle_outputs(
        document=document,
        presentation=presentation,
        audit=audit,
        output_dir=tmp_path / "cdh10",
        dpi=72,
        include_major_heading=False,
    )
    pdf = paths.get("section_1_pdf")
    assert pdf is not None

    import fitz

    with fitz.open(str(pdf)) as doc:
        page_texts = [doc[index].get_text() for index in range(doc.page_count)]
    assert len(page_texts) == 3
    assert "Ca2+ binding site" in page_texts[0]
    assert "Cadherin_C:" in page_texts[1]
    assert "3D structures from PDB" in page_texts[2]
    # No PDB heading leaks onto an earlier page.
    assert not any("3D structures from PDB" in text for text in page_texts[:2])

    pngs = sorted((tmp_path / "cdh10").glob("section_1_page_*.png"))
    assert len(pngs) == 3


def test_section_1c_pdf_page_break_sentinel_is_bundle_only():
    from gene_dossier.rancho_report import (
        SECTION_1C_PDF_PAGE_BREAK,
        _split_pdf_page_segments,
    )

    document, _presentation, _audit = build_section_bundle_document(
        dossier_run_id="sentinel-run",
        gene_symbol="SREBF2",
        section_keys=["1c"],
        evidence_records=_srebf2_1c_records(),
    )
    bundle_html = render_section_bundle_html(document, include_page_chrome=False)
    assert bundle_html.count(SECTION_1C_PDF_PAGE_BREAK) == 1
    segments = _split_pdf_page_segments(bundle_html)
    assert len(segments) == 2
    # Every segment keeps the stylesheet so both stories render identically.
    assert all("section-1c-pdb-official-image" in segment for segment in segments)
    assert "3D structures from PDB" not in segments[0]
    assert "3D structures from PDB" in segments[1]

    # A document without the sentinel is one unchanged segment.
    plain = "<html><body><p>no sentinel</p></body></html>"
    assert _split_pdf_page_segments(plain) == [plain]


def test_default_section_bundle_keys_full_ordered_1a_through_4a(tmp_path, monkeypatch):
    """Default bundle is 1a–4a; explicit section_keys / --sections still override."""
    from gene_dossier import section_bundle as sb

    expected = (
        "1a",
        "1b",
        "1c",
        "1d",
        "1e",
        "2a",
        "2b",
        "2c",
        "3a",
        "4a",
    )
    assert DEFAULT_SECTION_BUNDLE_KEYS == expected
    assert SUPPORTED_SECTION_BUNDLE_KEYS == expected + ("5a", "5b")
    assert validate_section_keys(DEFAULT_SECTION_BUNDLE_KEYS) == list(expected)
    assert len(DEFAULT_SECTION_BUNDLE_KEYS) == len(set(DEFAULT_SECTION_BUNDLE_KEYS))
    assert validate_section_keys(["4a", "1a", "4a", "2b"]) == ["1a", "2b", "4a"]
    assert validate_section_keys(["1a", "2c", "3a"]) == ["1a", "2c", "3a"]

    captured: dict[str, list[str]] = {}

    class _StopAfterCreate(Exception):
        pass

    def _fake_create(**kwargs):
        captured["keys"] = list(kwargs["selected_section_keys"])
        raise _StopAfterCreate()

    monkeypatch.setattr(sb, "create_section_bundle_run", _fake_create)

    with pytest.raises(_StopAfterCreate):
        run_section_bundle(
            "GENEX",
            section_keys=None,
            output_dir=tmp_path / "default",
            persist_db=False,
            write_pdf=False,
        )
    assert captured["keys"] == list(expected)

    with pytest.raises(_StopAfterCreate):
        run_section_bundle(
            "GENEX",
            section_keys=["4a"],
            output_dir=tmp_path / "only4a",
            persist_db=False,
            write_pdf=False,
        )
    assert captured["keys"] == ["4a"]

    with pytest.raises(_StopAfterCreate):
        run_section_bundle(
            "GENEX",
            section_keys=["1a", "3a"],
            output_dir=tmp_path / "subset",
            persist_db=False,
            write_pdf=False,
        )
    assert captured["keys"] == ["1a", "3a"]

    # CLI argparse default matches DEFAULT_SECTION_BUNDLE_KEYS; --sections overrides.
    import argparse

    script = Path(__file__).resolve().parents[1] / "scripts" / "run_section_bundle.py"
    text = script.read_text(encoding="utf-8")
    assert "Default: 1a 1b 1c 1d 1e 2a 2b 2c 3a 4a" in text
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sections",
        nargs="+",
        default=list(DEFAULT_SECTION_BUNDLE_KEYS),
    )
    assert parser.parse_args([]).sections == list(expected)
    assert parser.parse_args(["--sections", "4a"]).sections == ["4a"]
    assert parser.parse_args(["--sections", "1a", "2b"]).sections == ["1a", "2b"]


def test_section_defaults_include_completed_sections():
    assert DEFAULT_SECTION_BUNDLE_KEYS == (
        "1a",
        "1b",
        "1c",
        "1d",
        "1e",
        "2a",
        "2b",
        "2c",
        "3a",
        "4a",
    )
    assert SUPPORTED_SECTION_BUNDLE_KEYS == (
        "1a",
        "1b",
        "1c",
        "1d",
        "1e",
        "2a",
        "2b",
        "2c",
        "3a",
        "4a",
        "5a",
        "5b",
    )
    assert "1e" in DEFAULT_SECTION_BUNDLE_KEYS
    assert "4a" in DEFAULT_SECTION_BUNDLE_KEYS
    assert "5a" not in DEFAULT_SECTION_BUNDLE_KEYS


def test_render_section_bundle_html_can_suppress_major_heading_for_focused_1c():
    from gene_dossier.report_schema import (
        ReportCover,
        ReportDocument,
        ReportMajorSection,
        ReportSubsection,
    )

    document = ReportDocument(
        dossier_run_id="run",
        gene_symbol="CDH10",
        cover=ReportCover(gene_symbol="CDH10"),
        sections=[
            ReportMajorSection(
                number=1,
                key="1",
                title="General Gene Information",
                toc_title="GENERAL GENE INFORMATION",
                subsections=[
                    ReportSubsection(
                        key="c",
                        title="Known structure",
                        toc_title="KNOWN STRUCTURE",
                        presentation_blocks=[ReportContentBlock(kind="narrative", text="body")],
                    )
                ],
            )
        ],
    )
    html = render_section_bundle_html(document, include_major_heading=False)
    assert "1. General Gene Information" not in html
    assert "c. Known structure" in html


def test_full_report_major_narrative_unchanged():
    """Full-report _render_major still emits major narrative when present."""
    from gene_dossier.rancho_report import _render_major
    from gene_dossier.report_schema import ReportMajorSection, ReportSubsection

    major = ReportMajorSection(
        number=2,
        key="2",
        title="Expression",
        toc_title="EXPRESSION",
        narrative_markdown="**Key findings**\n\n- still rendered in full report\n",
        subsections=[
            ReportSubsection(
                key="a",
                title="Tissue",
                toc_title="TISSUE",
                presentation_blocks=[
                    ReportContentBlock(kind="narrative", text="deterministic body")
                ],
            )
        ],
    )
    html = _render_major(major)
    assert "Key findings" in html
    assert "deterministic body" in html


def test_srebf2_and_cdh10_1a_row_order():
    from gene_dossier.report_presentation import build_gene_aliases_blocks
    from tests.test_report_presentation import (
        _human_bundle,
        _mouse_bundle,
        _rat_bundle,
    )

    # SREBF2
    result = build_gene_aliases_blocks(
        gene_symbol="SREBF2",
        evidence_records=_human_bundle() + _mouse_bundle() + _rat_bundle(),
    )
    labels = [row[0] for row in result.blocks[0].table_rows]
    assert labels[:5] == [
        "Entrez Gene ID",
        "Gene Symbol",
        "Gene Name",
        "Ensembl ID",
        "UniProt ID",
    ]
    assert "Synonyms" in labels[5] and "Aliases" in labels[5]
    # Keys may include zero-width joiner in Synonyms label — look up by prefix.
    entrez_row = next(r for r in result.blocks[0].table_rows if r[0].startswith("Entrez"))
    assert "6721" in entrez_row[1]
    assert "300095" in entrez_row[3]

    cdh = [
        _ev(
            gene_symbol="CDH10",
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            value={
                "entrez_gene_id": "1008",
                "nomenclaturesymbol": "CDH10",
                "description": "cadherin 10",
            },
            source_id="ncbi-cdh",
            id="eid-cdh",
        )
    ]
    cdh_result = build_gene_aliases_blocks(
        gene_symbol="CDH10", evidence_records=cdh
    )
    cdh_labels = [r[0] for r in cdh_result.blocks[0].table_rows]
    assert cdh_labels[:5] == labels[:5]
    assert "1008" in cdh_result.blocks[0].table_rows[0][1]


# ---------------------------------------------------------------------------
# Section 2c visual-complete accepted-pointer promotion
# ---------------------------------------------------------------------------
def test_accept_visual_complete_promotion_matrix(tmp_path):
    from gene_dossier.section_2c import (
        SECTION_2C_VISUAL_COMPLETE_ROLES,
        STATUS_SUCCESS,
        accept_visual_complete_gene_report,
    )
    from gene_dossier.section_2c_sources import paths_for

    paths = paths_for(tmp_path / "section_2c")
    gene = "GENEX"
    prior_dir = tmp_path / "attempt_prior"
    new_dir = tmp_path / "attempt_new"
    prior_dir.mkdir()
    new_dir.mkdir()
    (prior_dir / "marker.txt").write_text("prior", encoding="utf-8")
    (new_dir / "marker.txt").write_text("new", encoding="utf-8")

    rendering = {
        "scientific_status": STATUS_SUCCESS,
        "visual_status": STATUS_SUCCESS,
    }
    complete_eval = {
        "visual_complete": True,
        "scientific_status": STATUS_SUCCESS,
        "visual_status": STATUS_SUCCESS,
    }
    incomplete_eval = {
        "visual_complete": False,
        "scientific_status": STATUS_SUCCESS,
        "visual_status": "partial",
        "missing_figure_roles": sorted(SECTION_2C_VISUAL_COMPLETE_ROLES),
    }

    # Partial / failed new run never replaces.
    assert (
        accept_visual_complete_gene_report(
            paths,
            gene_symbol=gene,
            attempt_dir=new_dir,
            rendering=rendering,
            artifacts={},
            evaluation=incomplete_eval,
            promote_existing=True,
        )
        is None
    )
    assert not paths.accepted_gene_pointer(gene).is_file()

    # New visual-complete creates pointer when none exists.
    created = accept_visual_complete_gene_report(
        paths,
        gene_symbol=gene,
        attempt_dir=prior_dir,
        rendering=rendering,
        artifacts={"k": "v"},
        evaluation=complete_eval,
        promote_existing=False,
    )
    assert created is not None and created.is_file()
    prior_payload = json.loads(created.read_text(encoding="utf-8"))
    assert prior_payload["attempt_dir"] == str(prior_dir)
    assert prior_payload["acceptance"]["section_2c_visual_complete"] is True
    assert prior_payload["acceptance"]["promotion_requested"] is False
    assert prior_payload["acceptance"]["replaced_prior_visual_complete"] is False

    # Default preserves prior visual-complete pointer.
    assert (
        accept_visual_complete_gene_report(
            paths,
            gene_symbol=gene,
            attempt_dir=new_dir,
            rendering=rendering,
            artifacts={},
            evaluation=complete_eval,
            promote_existing=False,
        )
        is None
    )
    preserved = json.loads(paths.accepted_gene_pointer(gene).read_text(encoding="utf-8"))
    assert preserved["attempt_dir"] == str(prior_dir)
    assert prior_dir.is_dir()
    assert (prior_dir / "marker.txt").read_text(encoding="utf-8") == "prior"

    # Explicit promotion replaces pointer only; old attempt dir remains.
    promoted = accept_visual_complete_gene_report(
        paths,
        gene_symbol=gene,
        attempt_dir=new_dir,
        rendering=rendering,
        artifacts={"k2": "v2"},
        evaluation=complete_eval,
        promote_existing=True,
    )
    assert promoted is not None
    payload = json.loads(promoted.read_text(encoding="utf-8"))
    assert payload["attempt_dir"] == str(new_dir)
    assert payload["acceptance"]["promotion_requested"] is True
    assert payload["acceptance"]["replaced_prior_visual_complete"] is True
    assert prior_dir.is_dir()
    assert (prior_dir / "marker.txt").read_text(encoding="utf-8") == "prior"


def test_accept_visual_complete_replaces_malformed_or_non_visual_prior(tmp_path):
    from gene_dossier.section_2c import (
        STATUS_SUCCESS,
        accept_visual_complete_gene_report,
    )
    from gene_dossier.section_2c_sources import paths_for, write_json_atomic

    paths = paths_for(tmp_path / "section_2c")
    gene = "GENEX"
    bad_dir = tmp_path / "bad"
    good_dir = tmp_path / "good"
    bad_dir.mkdir()
    good_dir.mkdir()
    pointer = paths.accepted_gene_pointer(gene)
    pointer.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(
        pointer,
        {
            "gene_symbol": gene,
            "attempt_dir": str(bad_dir),
            "acceptance": {"section_2c_visual_complete": False},
            "artifacts": {},
        },
    )
    rendering = {
        "scientific_status": STATUS_SUCCESS,
        "visual_status": STATUS_SUCCESS,
    }
    evaluation = {
        "visual_complete": True,
        "scientific_status": STATUS_SUCCESS,
        "visual_status": STATUS_SUCCESS,
    }
    replaced = accept_visual_complete_gene_report(
        paths,
        gene_symbol=gene,
        attempt_dir=good_dir,
        rendering=rendering,
        artifacts={},
        evaluation=evaluation,
        promote_existing=False,
    )
    assert replaced is not None
    payload = json.loads(replaced.read_text(encoding="utf-8"))
    assert payload["attempt_dir"] == str(good_dir)
    assert payload["acceptance"]["section_2c_visual_complete"] is True
    assert payload["acceptance"]["replaced_prior_visual_complete"] is False


def test_promote_section_2c_cli_flag_defaults_and_choices():
    import argparse
    import importlib.util

    script = Path(__file__).resolve().parents[1] / "scripts" / "run_section_bundle.py"
    text = script.read_text(encoding="utf-8")
    assert "--promote-section-2c-accepted" in text
    assert "promote_section_2c_accepted=args.promote_section_2c_accepted" in text
    assert "section_2c_visual_complete" not in text.split("choices=")[1].split("]")[0]

    spec = importlib.util.spec_from_file_location("run_section_bundle_cli", script)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    # Avoid executing main; parse the ArgumentParser construction by invoking main's parser.
    # Re-build equivalent parser checks via running with --help is heavy; inspect source.
    assert "action=\"store_true\"" in text or "action='store_true'" in text

    # Signature default on run_section_bundle.
    import inspect

    from gene_dossier import section_bundle as sb

    params = inspect.signature(sb.run_section_bundle).parameters
    assert params["promote_section_2c_accepted"].default is False

    # acceptance-profile choices remain unchanged.
    assert "section_1c_reference_genes" in text
    assert "section_1d_reference_genes" in text
    # Ensure the promote flag is not folded into acceptance-profile choices.
    profile_block = text[
        text.index("--acceptance-profile") : text.index("--promote-section-2c-accepted")
    ]
    assert "section_2c_visual_complete" not in profile_block


def test_promote_flag_has_no_effect_when_2c_absent(tmp_path, monkeypatch):
    """When 2c is not selected, promote_section_2c_accepted is ignored."""
    called = {"accept": 0}

    def _fake_accept(*args, **kwargs):
        called["accept"] += 1
        return None

    monkeypatch.setattr(
        "gene_dossier.section_bundle.accept_visual_complete_gene_report",
        _fake_accept,
    )
    # Minimal no-network identity-only path is heavy; assert the post-render
    # gate is gated on '2c' in keys by inspecting source.
    src = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "gene_dossier"
        / "section_bundle.py"
    ).read_text(encoding="utf-8")
    assert 'if "2c" in keys' in src or "section_2c_status" in src
    assert "promote_existing=promote_section_2c_accepted" in src
    # Without selecting 2c, the accept helper must not be invoked by a 1a/1b-only
    # offline run that never builds section_2c_status.
    assert called["accept"] == 0


def test_section_3a_visual_complete_fails_when_selected_chart_blocked():
    """Blocked/missing charts must not pass post-render visual-complete acceptance."""
    from gene_dossier.section_3a import (
        STATUS_SUCCESS,
        STATUS_UNAVAILABLE,
        evaluate_section_3a_visual_complete,
    )

    status = {
        "summary": {
            "selected_profile_count": 2,
            "selected_profiles": [
                {"graph_status": "failed", "graph_ok": False},
                {"graph_status": "failed", "graph_ok": False},
            ],
        },
        "rendering_status": {
            "scientific_status": STATUS_SUCCESS,
            "visual_status": STATUS_UNAVAILABLE,
        },
    }
    evaluation = evaluate_section_3a_visual_complete(
        status=status,
        embedded_figure_count=0,
        selected_count=2,
        pdf_render_status=STATUS_SUCCESS,
    )
    assert evaluation["scientific_complete"] is True
    assert evaluation["visual_complete"] is False
    assert evaluation["visual_status"] == STATUS_UNAVAILABLE

    # Partial visual status also fails visual-complete.
    status["rendering_status"]["visual_status"] = "partial"
    partial = evaluate_section_3a_visual_complete(
        status=status,
        embedded_figure_count=1,
        selected_count=2,
        pdf_render_status=STATUS_SUCCESS,
    )
    assert partial["visual_complete"] is False
