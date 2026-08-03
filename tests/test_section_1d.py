"""Offline tests for Section 1d AlphaFold selection and presentation."""

from __future__ import annotations

from pathlib import Path

import pytest

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.report_presentation import build_alphafold_blocks, build_section_presentation
from gene_dossier.report_schema import (
    ReportContentBlock,
    ReportCover,
    ReportDocument,
    ReportMajorSection,
    ReportSubsection,
)
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SUPPORTED_SECTION_BUNDLE_KEYS,
    assign_opaque_refs,
    build_section_bundle_document,
    render_section_bundle_html,
    sources_for_sections,
    validate_section_keys,
)
from gene_dossier.source_ids import make_source_id
from gene_dossier.tools import alphafold


def _pred(**overrides):
    base = {
        "modelEntityId": "AF-Q12772-F1",
        "uniprotAccession": "Q12772",
        "entityType": "protein",
        "isUniProt": True,
        "isComplex": False,
        "taxId": 9606,
        "isUniProtReviewed": True,
        "isUniProtReferenceProteome": True,
        "providerId": "GDM",
        "toolUsed": "AlphaFold Monomer v2.0",
        "latestVersion": 4,
        "sequenceStart": 1,
        "sequenceEnd": 100,
        "sequence": "M" * 20,
    }
    base.update(overrides)
    return base


def test_select_prefers_model_entity_id_and_rejects_isoform():
    predictions = [
        _pred(modelEntityId="AF-Q12772-2-F1", uniprotAccession="Q12772-2", latestVersion=9),
        _pred(modelEntityId="AF-Q12772-F1", uniprotAccession="Q12772", latestVersion=4),
    ]
    selected, diags = alphafold.select_canonical_monomer_prediction(
        predictions, "Q12772", expected_taxon_id=9606
    )
    assert selected is not None
    assert alphafold.model_entity_id(selected) == "AF-Q12772-F1"
    assert any(d.get("code") == "rejected_accession_mismatch" for d in diags)


def test_select_new_schema_fields_without_legacy_aliases():
    predictions = [
        {
            "modelEntityId": "AF-Q9Y6N8-F1",
            "uniprotAccession": "Q9Y6N8",
            "entityType": "protein",
            "isUniProt": True,
            "isComplex": False,
            "taxId": 9606,
            "isUniProtReviewed": True,
            "isUniProtReferenceProteome": True,
            "providerId": "GDM",
            "toolUsed": "AlphaFold Monomer",
            "latestVersion": 4,
            "sequenceStart": 1,
            "sequenceEnd": 50,
            "sequence": "ACDEFGHIKL",
            # paeImageUrl intentionally absent
        }
    ]
    selected, _ = alphafold.select_canonical_monomer_prediction(
        predictions, "Q9Y6N8", expected_taxon_id=9606
    )
    assert selected is not None
    summary = alphafold.summarize_prediction(selected)
    assert summary["model_entity_id"] == "AF-Q9Y6N8-F1"
    assert summary["entry_url"] == "https://alphafold.ebi.ac.uk/entry/AF-Q9Y6N8-F1"
    assert summary["pae_image_url"] is None
    assert summary["is_uniprot_reviewed"] is True


def test_select_mixed_old_and_new_fields_new_wins():
    predictions = [
        {
            "modelEntityId": "AF-Q12772-F1",
            "entryId": "LEGACY-SHOULD-LOSE",
            "uniprotAccession": "Q12772",
            "entityType": "protein",
            "isUniProt": True,
            "isComplex": False,
            "taxId": 9606,
            "isUniProtReviewed": True,
            "isReviewed": False,
            "isUniProtReferenceProteome": True,
            "isReferenceProteome": False,
            "sequenceStart": 1,
            "uniprotStart": 99,
            "sequenceEnd": 10,
            "uniprotEnd": 999,
            "sequence": "NEWSEQ",
            "uniprotSequence": "OLDSEQ",
            "latestVersion": 4,
            "providerId": "GDM",
            "toolUsed": "AlphaFold Monomer",
        }
    ]
    selected, _ = alphafold.select_canonical_monomer_prediction(
        predictions, "Q12772", expected_taxon_id=9606
    )
    assert alphafold.model_entity_id(selected) == "AF-Q12772-F1"
    assert alphafold.is_uniprot_reviewed(selected) is True
    assert alphafold.sequence_start(selected) == 1
    assert alphafold.sequence_value(selected) == "NEWSEQ"


def test_select_rejects_complex_taxon_mismatch_and_non_uniprot():
    predictions = [
        _pred(isComplex=True, modelEntityId="AF-Q12772-F1"),
        _pred(isUniProt=False, modelEntityId="AF-Q12772-F1"),
        _pred(taxId=10090, modelEntityId="AF-Q12772-F1"),
        _pred(entityType="complex", modelEntityId="AF-Q12772-F1"),
    ]
    selected, diags = alphafold.select_canonical_monomer_prediction(
        predictions, "Q12772", expected_taxon_id=9606
    )
    assert selected is None
    codes = {d.get("code") for d in diags}
    assert "no_hard_qualified_candidate" in codes


def test_select_does_not_prefer_higher_version_isoform_over_f1():
    predictions = [
        _pred(
            modelEntityId="AF-Q12772-F1",
            uniprotAccession="Q12772",
            latestVersion=3,
        ),
        _pred(
            modelEntityId="AF-OTHER-F1",
            uniprotAccession="Q12772",
            latestVersion=99,
            isUniProtReviewed=True,
        ),
    ]
    selected, _ = alphafold.select_canonical_monomer_prediction(
        predictions, "Q12772", expected_taxon_id=9606
    )
    assert alphafold.model_entity_id(selected) == "AF-Q12772-F1"


def _prediction_record(species: str, accession: str, symbol: str) -> EvidenceRecord:
    mid = f"AF-{accession}-F1"
    return EvidenceRecord(
        source_id=make_source_id(
            "AlphaFold", "GENE", AssertionType.protein_structure, f"{species}-{accession}"
        ),
        dossier_run_id="run-1d",
        gene_symbol="GENE",
        section="AlphaFold / PDBe / CDD",
        subsection="AlphaFold prediction",
        source_name="AlphaFold",
        source_type=SourceType.structure_database,
        assertion_type=AssertionType.protein_structure,
        fact_type="alphafold_species_prediction",
        species=species,
        evidence_grade=EvidenceGrade.E,
        value={
            "availability_status": "selected",
            "species_key": species,
            "species_label": species.title(),
            "display_symbol": symbol,
            "uniprot_accession": accession,
            "model_entity_id": mid,
            "entry_url": f"https://alphafold.ebi.ac.uk/entry/{mid}",
            "presentation_item_key": f"alphafold-{species}-{accession.lower()}",
        },
        display_text=f"{species} prediction",
    )


def test_build_alphafold_blocks_status_lines_have_no_evidence_refs():
    status = {
        "species_slots": [
            {
                "species_key": "human",
                "species_label": "Human",
                "display_symbol": "SREBF2",
                "accession": "Q12772",
                "status": "selected",
                "entry_url": "https://alphafold.ebi.ac.uk/entry/AF-Q12772-F1",
                "presentation_item_key": "alphafold-human-q12772",
            },
            {
                "species_key": "rat",
                "species_label": "Rat",
                "display_symbol": "Srebf2",
                "accession": None,
                "status": "accession_unavailable",
                "message": "Rat Srebf2: UniProt accession not available",
                "presentation_item_key": "alphafold-rat-unavailable",
            },
            {
                "species_key": "mouse",
                "species_label": "Mouse",
                "display_symbol": "Srebf2",
                "accession": "Q3U1N2",
                "status": "model_absent",
                "message": "Mouse Srebf2: AlphaFold prediction not available",
                "presentation_item_key": "alphafold-mouse-unavailable",
            },
        ]
    }
    records = [_prediction_record("human", "Q12772", "SREBF2")]
    result = build_alphafold_blocks(
        gene_symbol="SREBF2",
        evidence_records=records,
        section_status=status,
    )
    polished, ref_map = assign_opaque_refs(section_key="1d", blocks=result.blocks)
    roles = [b.presentation_role for b in polished]
    assert "section_1d_species_link" in roles
    assert "section_1d_confidence_legend" not in roles
    assert roles.count("section_1d_species_status") == 3
    viz = next(
        b
        for b in polished
        if b.presentation_role == "section_1d_species_status"
        and "visualization temporarily unavailable" in (b.text or "")
    )
    assert viz.evidence_ref is None
    for block in polished:
        if block.presentation_role == "section_1d_species_status":
            assert block.evidence_ref is None
            assert "unavailable" in (block.presentation_item_key or "")
    assert all("unavailable" not in ref for ref in ref_map)
    assert any(ref.startswith("ev-1d-alphafold-human-q12772-summary") for ref in ref_map)


def test_sources_for_sections_excludes_generic_alphafold_when_1d_selected():
    assert sources_for_sections(["1d"]) == []
    assert "AlphaFold" not in sources_for_sections(["1a", "1b", "1c", "1d"])
    assert SUPPORTED_SECTION_BUNDLE_KEYS == ("1a", "1b", "1c", "1d", "1e", "2a")
    assert DEFAULT_SECTION_BUNDLE_KEYS == ("1a", "1b")
    assert validate_section_keys(["1d", "1a"]) == ["1a", "1d"]
    assert validate_section_keys(["1e", "1a"]) == ["1a", "1e"]
    assert sources_for_sections(["1e"]) == []
    assert "NCBI Datasets" not in sources_for_sections(["1a", "1e"])


def test_assembled_1d_merges_after_final_1c_pdb_segment(tmp_path):
    # Minimal 1c PDB page-break + 1d blocks.
    pdb_link = ReportContentBlock(
        kind="link",
        text="PDB",
        presentation_role="section_1c_pdb_link",
        presentation_item_key="pdb-1abc",
        presentation_page_break_before=True,
        links=[{"label": "PDB", "url": "https://www.ebi.ac.uk/pdbe/entry/pdb/1abc"}],
    )
    pdb_img = ReportContentBlock(
        kind="narrative",
        text="pdb image placeholder",
        presentation_role="section_1c_pdb_official_image",
        presentation_item_key="pdb-1abc",
    )
    domain = ReportContentBlock(
        kind="narrative",
        text="domain page",
        presentation_role="section_1c_domain_summary",
        presentation_item_key="cd-test",
    )
    human_link = ReportContentBlock(
        kind="link",
        text="Human SREBF2: ",
        presentation_role="section_1d_species_link",
        presentation_item_key="alphafold-human-q12772",
        links=[
            {
                "label": "AlphaFold Protein Structure Link",
                "url": "https://alphafold.ebi.ac.uk/entry/AF-Q12772-F1",
            }
        ],
    )
    document = ReportDocument(
        dossier_run_id="run",
        gene_symbol="SREBF2",
        cover=ReportCover(gene_symbol="SREBF2"),
        sections=[
            ReportMajorSection(
                number=1,
                key="1",
                title="Gene information",
                toc_title="GENE INFORMATION",
                subsections=[
                    ReportSubsection(
                        key="c",
                        title="Known structure",
                        toc_title="KNOWN STRUCTURE",
                        presentation_blocks=[domain, pdb_link, pdb_img],
                        status="populated",
                    ),
                    ReportSubsection(
                        key="d",
                        title="AlphaFold protein structure prediction",
                        toc_title="ALPHAFOLD PROTEIN STRUCTURE PREDICTION",
                        presentation_blocks=[human_link],
                        status="populated",
                    ),
                ],
                status="populated",
            )
        ],
    )
    html = render_section_bundle_html(document, include_page_chrome=False)
    assert "AlphaFold protein structure prediction" in html
    assert "subsection-d" in html
    # 1d must appear in a continuation page together with PDB content, not before it.
    pdb_pos = html.find("section-1c-pdb")
    if pdb_pos < 0:
        pdb_pos = html.find("section_1c_pdb") if False else html.find("PDB")
    d_pos = html.find("subsection-d")
    assert d_pos > pdb_pos
    assert html.find("subsection-d") > html.find("section-1c-continuation") or (
        "section-1c-continuation" in html and "subsection-d" in html.split("section-1c-continuation")[-1]
    )


def test_build_section_presentation_routes_1d():
    result = build_section_presentation(
        section_key="1d",
        gene_symbol="SREBF2",
        evidence_records=[_prediction_record("human", "Q12772", "SREBF2")],
        section_status={
            "species_slots": [
                {
                    "species_key": "human",
                    "species_label": "Human",
                    "display_symbol": "SREBF2",
                    "status": "selected",
                    "entry_url": "https://alphafold.ebi.ac.uk/entry/AF-Q12772-F1",
                    "presentation_item_key": "alphafold-human-q12772",
                }
            ]
        },
    )
    assert result.blocks
    assert result.blocks[0].presentation_role == "section_1d_species_link"
    assert result.blocks[1].presentation_role == "section_1d_species_status"
    assert "visualization temporarily unavailable" in (result.blocks[1].text or "")
    assert all(b.presentation_role != "section_1d_confidence_legend" for b in result.blocks)


def _make_capture_png_bytes() -> bytes:
    from io import BytesIO
    import random

    from PIL import Image, ImageDraw

    img = Image.new("RGB", (800, 600), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    random.seed(7)
    colors = [(0, 83, 214), (101, 203, 243), (255, 219, 19), (255, 125, 69), (30, 30, 30)]
    for _ in range(2500):
        x1, y1 = random.randint(0, 799), random.randint(0, 599)
        x2, y2 = random.randint(0, 799), random.randint(0, 599)
        draw.line((x1, y1, x2, y2), fill=random.choice(colors), width=2)
    buf = BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def _write_managed_capture(
    tmp_path: Path,
    *,
    model_id: str = "AF-Q12772-F1",
    accession: str = "Q12772",
    model_version: int = 4,
    mutate_checksum: bool = False,
) -> tuple[object, EvidenceRecord, dict]:
    from gene_dossier.config import Settings
    from gene_dossier.ucsc_figure import sha256_hex

    root = tmp_path / "raw"
    root.mkdir(parents=True, exist_ok=True)
    rel = f"tests/alphafold/{accession.lower()}-viewer.png"
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    content = _make_capture_png_bytes()
    path.write_bytes(content)
    digest = sha256_hex(content)
    if mutate_checksum:
        digest = "0" * 64
    value = {
        "status": "success",
        "model_entity_id": model_id,
        "uniprot_accession": accession,
        "taxon_id": 9606,
        "model_version": model_version,
        "source_page_url": f"https://alphafold.ebi.ac.uk/entry/{model_id}",
        "relative_path": rel,
        "media_type": "image/png",
        "width": 800,
        "height": 600,
        "sha256": digest,
        "byte_size": len(content),
        "artifact_class": "derived_capture",
        "retrieval_method": "official_web_element_capture",
        "presentation_item_key": f"alphafold-human-{accession.lower()}",
        "figure_raw_artifact_id": "raw-capture-1",
    }
    rec = EvidenceRecord(
        source_id=make_source_id(
            "AlphaFold", "SREBF2", AssertionType.protein_structure, f"viewer-{accession}"
        ),
        dossier_run_id="run-1d",
        gene_symbol="SREBF2",
        section="AlphaFold / PDBe / CDD",
        subsection="AlphaFold prediction",
        source_name="AlphaFold",
        source_type=SourceType.structure_database,
        assertion_type=AssertionType.protein_structure,
        fact_type="alphafold_official_viewer_capture",
        species="human",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.E,
        value=value,
        display_text="capture",
        raw_artifact_id="raw-capture-1",
    )
    settings = Settings(raw_data_dir=root)
    return settings, rec, value


def _three_species_status() -> dict:
    return {
        "species_slots": [
            {
                "species_key": "human",
                "species_label": "Human",
                "display_symbol": "SREBF2",
                "accession": "Q12772",
                "status": "selected",
                "model_entity_id": "AF-Q12772-F1",
                "entry_url": "https://alphafold.ebi.ac.uk/entry/AF-Q12772-F1",
                "presentation_item_key": "alphafold-human-q12772",
            },
            {
                "species_key": "rat",
                "species_label": "Rat",
                "display_symbol": "Srebf2",
                "accession": "Q3T1I5",
                "status": "selected",
                "model_entity_id": "AF-Q3T1I5-F1",
                "entry_url": "https://alphafold.ebi.ac.uk/entry/AF-Q3T1I5-F1",
                "presentation_item_key": "alphafold-rat-q3t1i5",
            },
            {
                "species_key": "mouse",
                "species_label": "Mouse",
                "display_symbol": "Srebf2",
                "accession": "Q3U1N2",
                "status": "selected",
                "model_entity_id": "AF-Q3U1N2-F1",
                "entry_url": "https://alphafold.ebi.ac.uk/entry/AF-Q3U1N2-F1",
                "presentation_item_key": "alphafold-mouse-q3u1n2",
            },
        ],
        "audit": {
            "viewer_capture": {
                "status": "success",
                "capture_mode": "reused_official_capture",
            }
        },
        "viewer_capture": {
            "status": "success",
            "capture_mode": "reused_official_capture",
        },
    }


def test_capture_plus_legend_render_as_visual_table(tmp_path, monkeypatch):
    from gene_dossier.config import get_settings
    from gene_dossier.rancho_report import _render_section_1d_blocks
    from gene_dossier.ucsc_figure import sha256_hex

    settings, capture, value = _write_managed_capture(tmp_path)
    monkeypatch.setenv("RAW_DATA_DIR", str(settings.raw_data_dir))
    get_settings.cache_clear()

    records = [
        _prediction_record("human", "Q12772", "SREBF2"),
        _prediction_record("rat", "Q3T1I5", "Srebf2"),
        _prediction_record("mouse", "Q3U1N2", "Srebf2"),
        capture,
    ]
    result = build_alphafold_blocks(
        gene_symbol="SREBF2",
        evidence_records=records,
        section_status=_three_species_status(),
    )
    roles = [b.presentation_role for b in result.blocks]
    assert roles[0] == "section_1d_species_link"
    assert roles[1] == "section_1d_human_structure_capture"
    assert roles[2] == "section_1d_confidence_legend"
    assert roles[3] == "section_1d_species_link"
    assert roles[4] == "section_1d_species_link"
    assert "section_1d_species_status" not in roles

    html = _render_section_1d_blocks(list(result.blocks))
    assert "section-1d-visual-table" in html
    assert html.count("<td") == 3
    assert html.index("section-1d-human-structure-capture") < html.index(
        "section-1d-confidence-legend"
    )
    assert html.index("section-1d-confidence-legend") < html.index('class="section-1d-blurb"')
    assert html.index("</table>") < html.index("alphafold.ebi.ac.uk/entry/AF-Q3T1I5-F1")
    assert html.index("AF-Q3T1I5-F1") < html.index("AF-Q3U1N2-F1")
    assert "section-1d-visual-row" not in html
    assert html.count("<td") == 3
    assert "section-1d-visual-table" in html
    assert sha256_hex((settings.raw_data_path / value["relative_path"]).read_bytes()) == value[
        "sha256"
    ]
    get_settings.cache_clear()


def test_capture_absent_skips_legend_and_emits_unavailable_status():
    records = [
        _prediction_record("human", "Q12772", "SREBF2"),
        _prediction_record("rat", "Q3T1I5", "Srebf2"),
        _prediction_record("mouse", "Q3U1N2", "Srebf2"),
    ]
    result = build_alphafold_blocks(
        gene_symbol="SREBF2",
        evidence_records=records,
        section_status=_three_species_status(),
    )
    roles = [b.presentation_role for b in result.blocks]
    assert "section_1d_confidence_legend" not in roles
    assert "section_1d_human_structure_capture" not in roles
    status = next(
        b for b in result.blocks if b.presentation_role == "section_1d_species_status"
    )
    assert "visualization temporarily unavailable" in (status.text or "")
    polished, _ = assign_opaque_refs(section_key="1d", blocks=result.blocks)
    assert all(
        b.evidence_ref is None
        for b in polished
        if b.presentation_role == "section_1d_species_status"
    )


def test_find_reusable_official_capture_accepts_and_rejects(tmp_path, monkeypatch):
    from gene_dossier.config import get_settings
    from gene_dossier.section_1d import find_reusable_official_capture
    from gene_dossier.ucsc_figure import sha256_hex

    settings, capture, value = _write_managed_capture(tmp_path)
    monkeypatch.setenv("RAW_DATA_DIR", str(settings.raw_data_dir))
    get_settings.cache_clear()

    reused, diags = find_reusable_official_capture(
        evidence_records=[capture],
        raw_artifacts=None,
        model_id="AF-Q12772-F1",
        accession="Q12772",
        model_version=4,
        settings=settings,
    )
    assert reused is not None
    assert reused["relative_path"] == value["relative_path"]
    assert any(d.get("code") == "reuse_accepted" for d in diags)

    bad_checksum = capture.model_copy(
        update={"value": {**value, "sha256": "1" * 64}}
    )
    reused, diags = find_reusable_official_capture(
        evidence_records=[bad_checksum],
        raw_artifacts=None,
        model_id="AF-Q12772-F1",
        accession="Q12772",
        model_version=4,
        settings=settings,
    )
    assert reused is None
    assert any(d.get("reason") == "checksum_mismatch" for d in diags)

    wrong_model = capture.model_copy(
        update={"value": {**value, "model_entity_id": "AF-OTHER-F1"}}
    )
    reused, diags = find_reusable_official_capture(
        evidence_records=[wrong_model],
        raw_artifacts=None,
        model_id="AF-Q12772-F1",
        accession="Q12772",
        model_version=4,
        settings=settings,
    )
    assert reused is None
    assert any(d.get("reason") == "model_entity_id_mismatch" for d in diags)

    wrong_version = capture.model_copy(update={"value": {**value, "model_version": 99}})
    reused, diags = find_reusable_official_capture(
        evidence_records=[wrong_version],
        raw_artifacts=None,
        model_id="AF-Q12772-F1",
        accession="Q12772",
        model_version=4,
        settings=settings,
    )
    assert reused is None
    assert any(d.get("reason") == "model_version_mismatch" for d in diags)

    missing = capture.model_copy(
        update={"value": {**value, "relative_path": "tests/alphafold/missing.png"}}
    )
    reused, diags = find_reusable_official_capture(
        evidence_records=[missing],
        raw_artifacts=None,
        model_id="AF-Q12772-F1",
        accession="Q12772",
        model_version=4,
        settings=settings,
    )
    assert reused is None
    assert any(d.get("reason") == "missing_managed_file" for d in diags)

    # Provenance retained on accepted reuse.
    reused, _ = find_reusable_official_capture(
        evidence_records=[capture],
        raw_artifacts=None,
        model_id="AF-Q12772-F1",
        accession="Q12772",
        model_version=4,
        settings=settings,
    )
    assert reused is not None
    assert reused.get("figure_raw_artifact_id") == "raw-capture-1"
    assert reused.get("sha256") == sha256_hex(
        (settings.raw_data_path / value["relative_path"]).read_bytes()
    )
    get_settings.cache_clear()


def test_cdh10_rat_resolves_from_taxon_identity_evidence():
    from gene_dossier.section_1d import resolve_species_identity
    import inspect
    import gene_dossier.section_1d as s1d

    records = [
        EvidenceRecord(
            source_id=make_source_id(
                "HGNC", "CDH10", AssertionType.gene_identity, "human"
            ),
            dossier_run_id="run",
            gene_symbol="CDH10",
            section="Gene aliases",
            source_name="HGNC",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="gene_symbol",
            taxon_id=9606,
            evidence_grade=EvidenceGrade.A,
            value={"gene_symbol": "CDH10", "uniprot_accession": "Q9Y6N8"},
            display_text="human",
        ),
        EvidenceRecord(
            source_id=make_source_id(
                "RGD", "CDH10", AssertionType.gene_identity, "rat"
            ),
            dossier_run_id="run",
            gene_symbol="CDH10",
            section="Gene aliases",
            source_name="RGD",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ortholog",
            taxon_id=10116,
            evidence_grade=EvidenceGrade.A,
            value={"gene_symbol": "Cdh10", "uniprot_accession": "F1LR98"},
            display_text="rat",
        ),
        EvidenceRecord(
            source_id=make_source_id(
                "MGI", "CDH10", AssertionType.gene_identity, "mouse"
            ),
            dossier_run_id="run",
            gene_symbol="CDH10",
            section="Gene aliases",
            source_name="MGI",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ortholog",
            taxon_id=10090,
            evidence_grade=EvidenceGrade.A,
            value={"gene_symbol": "Cdh10", "uniprot_id": "P70408"},
            display_text="mouse",
        ),
    ]
    resolved = resolve_species_identity(
        gene_symbol="CDH10",
        gene_ids={},
        evidence_records=records,
    )
    assert [r["species_key"] for r in resolved] == ["human", "rat", "mouse"]
    by_key = {r["species_key"]: r for r in resolved}
    assert by_key["human"]["accession"] == "Q9Y6N8"
    assert by_key["rat"]["accession"] == "F1LR98"
    assert by_key["mouse"]["accession"] == "P70408"
    source = inspect.getsource(s1d)
    assert "F1LR98" not in source
    assert 'gene_symbol == "CDH10"' not in source
    assert '== "CDH10"' not in source


def test_section_1d_reference_genes_acceptance_profile(tmp_path, monkeypatch):
    from gene_dossier.config import get_settings
    from gene_dossier.rancho_report import _render_section_1d_blocks
    from gene_dossier.section_1d import evaluate_section_1d_reference_genes_acceptance

    settings, capture, _ = _write_managed_capture(tmp_path)
    monkeypatch.setenv("RAW_DATA_DIR", str(settings.raw_data_dir))
    get_settings.cache_clear()

    status = _three_species_status()
    records = [
        _prediction_record("human", "Q12772", "SREBF2"),
        _prediction_record("rat", "Q3T1I5", "Srebf2"),
        _prediction_record("mouse", "Q3U1N2", "Srebf2"),
        capture,
    ]
    result = build_alphafold_blocks(
        gene_symbol="SREBF2",
        evidence_records=records,
        section_status=status,
    )
    html = _render_section_1d_blocks(list(result.blocks))
    reasons = evaluate_section_1d_reference_genes_acceptance(
        gene_symbol="SREBF2",
        section_status=status,
        presentation_blocks=result.blocks,
        html=html,
        pdf_path=None,
        selected_section_keys=["1d"],
    )
    assert reasons == []

    # Fail without capture.
    no_capture = build_alphafold_blocks(
        gene_symbol="SREBF2",
        evidence_records=records[:-1],
        section_status={
            **status,
            "viewer_capture": {"status": "unavailable", "capture_mode": "capture_unavailable"},
            "audit": {
                "viewer_capture": {
                    "status": "unavailable",
                    "capture_mode": "capture_unavailable",
                }
            },
        },
    )
    reasons = evaluate_section_1d_reference_genes_acceptance(
        gene_symbol="SREBF2",
        section_status={
            **status,
            "viewer_capture": {"status": "unavailable", "capture_mode": "capture_unavailable"},
        },
        presentation_blocks=no_capture.blocks,
        html=_render_section_1d_blocks(list(no_capture.blocks)),
        pdf_path=None,
        selected_section_keys=["1d"],
    )
    assert any("viewer capture" in r for r in reasons)

    # Fail when rat/mouse link missing.
    truncated = _three_species_status()
    truncated["species_slots"] = truncated["species_slots"][:1]
    reasons = evaluate_section_1d_reference_genes_acceptance(
        gene_symbol="SREBF2",
        section_status=truncated,
        presentation_blocks=result.blocks[:2],
        html=html,
        pdf_path=None,
        selected_section_keys=["1d"],
    )
    assert any("missing" in r for r in reasons)
    get_settings.cache_clear()


def test_no_custom_coordinate_rendering_helpers_and_defaults():
    import gene_dossier.section_1d as s1d

    assert not hasattr(s1d, "render_pymol_png")
    assert not hasattr(s1d, "render_mmcif_projection_png")
    assert DEFAULT_SECTION_BUNDLE_KEYS == ("1a", "1b")
    assert SUPPORTED_SECTION_BUNDLE_KEYS == ("1a", "1b", "1c", "1d", "1e", "2a")
    assert "1d" not in DEFAULT_SECTION_BUNDLE_KEYS
    assert "1e" not in DEFAULT_SECTION_BUNDLE_KEYS
