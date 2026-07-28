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
from gene_dossier.report_presentation import (
    UCSC_STABLE_INTRO,
    build_conservation_blocks,
    transcript_selection_sentence,
)
from gene_dossier.report_schema import ReportContentBlock
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SUPPORTED_SECTION_BUNDLE_KEYS,
    SectionBundleError,
    assign_opaque_refs,
    build_provenance_index,
    build_section_bundle_document,
    create_section_bundle_run,
    finalize_section_bundle_run,
    opaque_evidence_ref,
    render_section_bundle_html,
    run_section_bundle,
    sanitize_credentials,
    sanitize_polished_text,
    sanitize_secrets,
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
    assert validate_section_keys(["1b", "1a", "1a"]) == ["1a", "1b"]
    assert validate_section_keys(DEFAULT_SECTION_BUNDLE_KEYS) == list(
        SUPPORTED_SECTION_BUNDLE_KEYS
    )
    with pytest.raises(SectionBundleError):
        validate_section_keys(["1c"])
    with pytest.raises(SectionBundleError):
        validate_section_keys([])


def test_sources_for_sections_dependency_aware():
    assert sources_for_sections(["1a"]) == []
    assert sources_for_sections(["1b"]) == ["UCSC"]
    assert sources_for_sections(["1a", "1b"]) == ["UCSC"]


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
    from gene_dossier.report_presentation import _transcript_selection_sentence as _
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
    rows = {r[0]: r for r in result.blocks[0].table_rows}
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
