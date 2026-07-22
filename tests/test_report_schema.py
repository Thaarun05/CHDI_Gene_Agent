"""Tests for the Rancho/CHDI report layout schema (no network / no LLM)."""

from __future__ import annotations

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.report_schema import (
    BACK_MATTER_TITLES,
    COMPILED_RELEVANT_DATABASES,
    REPORT_SECTIONS,
    REPORT_STYLE,
    ReportCover,
    build_report_document,
    cover_lines,
    infer_chromosome,
    iter_toc_entries,
    resolve_report_slot,
)
from gene_dossier.source_ids import make_source_id


def _ev(
    *,
    source_name: str,
    assertion_type: AssertionType,
    key: str,
    display_text: str,
    section: str = "General",
    subsection: str | None = None,
    fact_type: str = "fact",
    grade: EvidenceGrade = EvidenceGrade.C,
    source_type: SourceType = SourceType.curated_database,
    value: dict | None = None,
    gene_symbol: str = "SREBF2",
    dossier_run_id: str = "run-schema",
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=make_source_id(source_name, gene_symbol, assertion_type, key),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        section=section,
        subsection=subsection,
        source_name=source_name,
        source_type=source_type,
        assertion_type=assertion_type,
        fact_type=fact_type,
        evidence_grade=grade,
        display_text=display_text,
        value=value or {},
    )


# --------------------------------------------------------------------------------------
# Layout / style tokens
# --------------------------------------------------------------------------------------
def test_report_has_exactly_15_major_sections():
    assert len(REPORT_SECTIONS) == 15
    assert [s.number for s in REPORT_SECTIONS] == list(range(1, 16))
    assert [s.key for s in REPORT_SECTIONS] == [str(i) for i in range(1, 16)]


def test_section_titles_match_reference_wording():
    titles = [s.title for s in REPORT_SECTIONS]
    assert titles[0] == "General Gene Information"
    assert titles[3] == "Transcription factors that drive the gene’s expression"
    assert "’" in titles[3]  # typographic apostrophe
    assert titles[14] == "NIH and ERC grants"


def test_back_matter_and_compiled_databases():
    assert BACK_MATTER_TITLES == (
        "References",
        "Compiled List of Relevant Databases",
    )
    assert len(COMPILED_RELEVANT_DATABASES) >= 20
    names = [n for n, _ in COMPILED_RELEVANT_DATABASES]
    assert "GTEx Portal" in names
    assert "STRING" in names
    assert REPORT_STYLE.green_major.startswith("#")
    assert REPORT_STYLE.orange_sub.startswith("#")
    assert REPORT_STYLE.table_header_bg.startswith("#")


def test_iter_toc_entries_includes_majors_subs_and_back_matter():
    entries = iter_toc_entries()
    majors = [e for e in entries if e["kind"] == "major"]
    subs = [e for e in entries if e["kind"] == "subsection"]
    back = [e for e in entries if e["kind"] == "back_matter"]
    assert len(majors) == 15
    assert len(subs) > 0
    assert {e["display_title"] for e in back} == set(BACK_MATTER_TITLES)
    # Section 1a slot key used by TOC page-number maps
    assert any(e.get("slot") == "1a" for e in subs)


# --------------------------------------------------------------------------------------
# Slot mapping
# --------------------------------------------------------------------------------------
def test_resolve_report_slot_maps_core_srebf2_sources():
    cases = [
        (_ev(source_name="NCBI Gene", assertion_type=AssertionType.gene_identity,
             key="6721", display_text="Entrez 6721"), "1a"),
        (_ev(source_name="Ensembl", assertion_type=AssertionType.gene_identity,
             key="ENSG00000198911", display_text="ENSG00000198911"), "1a"),
        (_ev(source_name="UCSC", assertion_type=AssertionType.gene_identity,
             key="chr22", display_text="chr22:41833105-41907305"), "1b"),
        (_ev(source_name="UniProt", assertion_type=AssertionType.protein_function,
             key="Q12772", display_text="SREBP2"), "1c"),
        (_ev(source_name="AlphaFold", assertion_type=AssertionType.protein_structure,
             key="Q12772", display_text="AF model"), "1d"),
        (_ev(source_name="GTEx", assertion_type=AssertionType.expression,
             key="median", display_text="Brain high", section="Tissue expression"), "2a"),
        (_ev(source_name="Allen Brain Atlas", assertion_type=AssertionType.expression,
             key="1051154", display_text="HBA probe"), "2b"),
        (_ev(source_name="DropViz", assertion_type=AssertionType.cell_type_expression,
             key="neuron", display_text="Neuronal"), "2c"),
        (_ev(source_name="GEO", assertion_type=AssertionType.perturbation,
             key="GDS123", display_text="GEO hit"), "3a"),
        (_ev(source_name="Harmonizome",
             assertion_type=AssertionType.transcription_factor_association,
             key="TF", display_text="TF hit"), "4a"),
        (_ev(source_name="STRING", assertion_type=AssertionType.ppi,
             key="SCAP", display_text="SREBF2-SCAP"), "5a"),
        (_ev(source_name="BioGRID", assertion_type=AssertionType.ppi,
             key="SCAP", display_text="BioGRID hit"), "5b"),
        (_ev(source_name="CTD", assertion_type=AssertionType.chemical_interaction,
             key="statin", display_text="CTD hit"), "6a"),
        (_ev(source_name="ChEMBL", assertion_type=AssertionType.chemical_interaction,
             key="CHEMBL1", display_text="tool"), "7a"),
        (_ev(source_name="GTEx", assertion_type=AssertionType.expression,
             key="eqtl1", display_text="eQTL", section="eQTL brain",
             fact_type="eqtl"), "8a"),
        (_ev(source_name="ClinVar", assertion_type=AssertionType.variant_association,
             key="rs1", display_text="variant"), "9a"),
        (_ev(source_name="OMIM", assertion_type=AssertionType.disease_association,
             key="600481", display_text="MIM 600481"), "9b"),
        (_ev(source_name="Reactome", assertion_type=AssertionType.pathway_membership,
             key="R-HSA-1", display_text="pathway"), "10a"),
        (_ev(source_name="MGI", assertion_type=AssertionType.knockout_phenotype,
             key="MGI:107585", display_text="knockout"), "11b"),
        (_ev(source_name="PubMed", assertion_type=AssertionType.literature_summary,
             key="PMID:1", display_text="lab paper",
             source_type=SourceType.literature), "12"),
        (_ev(source_name="Antibodies", assertion_type=AssertionType.literature_summary,
             key="ab1", display_text="Abcam antibody",
             source_type=SourceType.commercial_source), "13"),
        (_ev(source_name="Patents", assertion_type=AssertionType.patent_claim,
             key="US1", display_text="patent",
             source_type=SourceType.patent_database), "14"),
        (_ev(source_name="NIH RePORTER", assertion_type=AssertionType.grant_project,
             key="R01", display_text="NIH grant",
             source_type=SourceType.grant_database), "15a"),
        (_ev(source_name="ERC", assertion_type=AssertionType.grant_project,
             key="ERC1", display_text="ERC grant",
             source_type=SourceType.grant_database), "15b"),
    ]
    for record, expected in cases:
        slot = resolve_report_slot(record)
        assert slot is not None, record.source_name
        assert slot.slot_id == expected, (record.source_name, slot.slot_id)


def test_unmapped_record_returns_none():
    record = _ev(
        source_name="Unknown Widget DB",
        assertion_type=AssertionType.visual_observation,
        key="x",
        display_text="noise",
        section="misc",
    )
    assert resolve_report_slot(record) is None


def test_infer_chromosome_from_display_text_and_value():
    from_text = _ev(
        source_name="UCSC",
        assertion_type=AssertionType.gene_identity,
        key="loc",
        display_text="Located on chr22:41833105-41907305",
    )
    assert infer_chromosome([from_text]) == "22"
    from_value = _ev(
        source_name="NCBI Gene",
        assertion_type=AssertionType.gene_identity,
        key="6721",
        display_text="gene",
        value={"chromosome": "22"},
    )
    assert infer_chromosome([from_value]) == "22"


# --------------------------------------------------------------------------------------
# build_report_document
# --------------------------------------------------------------------------------------
def test_build_report_document_buckets_evidence_and_keeps_provenance():
    records = [
        _ev(
            source_name="NCBI Gene",
            assertion_type=AssertionType.gene_identity,
            key="6721",
            display_text="SREBF2 Entrez Gene ID is 6721.",
            grade=EvidenceGrade.A,
        ),
        _ev(
            source_name="STRING",
            assertion_type=AssertionType.ppi,
            key="SCAP",
            display_text="SREBF2 interacts with SCAP.",
        ),
        _ev(
            source_name="Mystery Source",
            assertion_type=AssertionType.visual_observation,
            key="zzz",
            display_text="unmapped",
            section="misc",
        ),
    ]
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=records,
        curator="Test Curator",
        report_date="2026-07-22",
        chromosome="22",
        references=["Smith et al. 2020"],
    )

    assert doc.gene_symbol == "SREBF2"
    assert doc.cover.gene_line == "SREBF2 (CHR22)"
    assert doc.cover.curator == "Test Curator"
    assert len(doc.sections) == 15
    assert doc.references == ["Smith et al. 2020"]
    assert len(doc.compiled_databases) == len(COMPILED_RELEVANT_DATABASES)

    sec1 = doc.sections[0]
    aliases = next(s for s in sec1.subsections if s.key == "a")
    assert aliases.status == "populated"
    assert aliases.blocks[0].text.startswith("SREBF2 Entrez")
    assert aliases.source_ids
    assert aliases.blocks[0].source_ids == aliases.source_ids

    sec5 = doc.sections[4]
    string_sub = next(s for s in sec5.subsections if s.key == "a")
    assert string_sub.status == "populated"
    assert "SCAP" in (string_sub.blocks[0].text or "")

    assert len(doc.unmapped_source_ids) == 1
    assert doc.unmapped_source_ids[0] == records[2].source_id


def test_cover_lines_order():
    cover = ReportCover(
        gene_symbol="SREBF2",
        chromosome="22",
        curator="Ada",
        report_date="July 22, 2026",
    )
    lines = cover_lines(cover)
    assert lines[0] == "SREBF2 (CHR22)"
    assert lines[1] == "Gene Report"
    assert "Prepared for the CHDI Foundation" in lines
    assert any("Ada" in line for line in lines)
    assert any("July 22, 2026" in line for line in lines)
