"""Tests for the Rancho/CHDI report layout schema (no network / no LLM)."""

from __future__ import annotations

import pytest

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    ReportSection,
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
    resolve_report_slot_for_section,
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


def _report_section(
    *,
    section_name: str,
    content_markdown: str = "**Section**\n\nNarrative.",
    source_ids: list[str] | None = None,
    status: str = "llm",
    dossier_run_id: str = "run-schema",
) -> ReportSection:
    return ReportSection(
        dossier_run_id=dossier_run_id,
        section_name=section_name,
        content_markdown=content_markdown,
        source_ids=list(source_ids or []),
        status=status,
    )


CANONICAL_SYNTHESIS_SECTION_NAMES: list[str] = [
    "General gene information",
    "Gene aliases and identifiers",
    "Conservation / orthologs",
    "Known structure / domains",
    "AlphaFold / PDBe / CDD",
    "Homologues",
    "Tissue and cell expression",
    "GEO perturbations",
    "Transcription factors",
    "Protein-protein interactions",
    "CTD perturbations",
    "Chemical tools",
    "eQTLs",
    "ClinVar / OMIM / Open Targets / SNPs",
    "Pathways",
    "Knockouts / model phenotypes",
    "Major labs / literature",
    "Antibodies",
    "Patents",
    "NIH/ERC grants",
    "Missing / deferred / manual sources",
    "Verification warnings",
]

CANONICAL_CONTENT_SLOT_CASES: list[tuple[str, str]] = [
    ("General gene information", "1"),
    ("Gene aliases and identifiers", "1a"),
    ("Conservation / orthologs", "1b"),
    ("Known structure / domains", "1c"),
    ("AlphaFold / PDBe / CDD", "1d"),
    ("Homologues", "1e"),
    ("Tissue and cell expression", "2"),
    ("GEO perturbations", "3a"),
    ("Transcription factors", "4a"),
    ("Protein-protein interactions", "5"),
    ("CTD perturbations", "6a"),
    ("Chemical tools", "7"),
    ("eQTLs", "8a"),
    ("ClinVar / OMIM / Open Targets / SNPs", "9"),
    ("Pathways", "10a"),
    ("Knockouts / model phenotypes", "11"),
    ("Major labs / literature", "12"),
    ("Antibodies", "13"),
    ("Patents", "14"),
    ("NIH/ERC grants", "15"),
]


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


# --------------------------------------------------------------------------------------
# Synthesized ReportSection → Rancho slots
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("section_name,expected_slot", CANONICAL_CONTENT_SLOT_CASES)
def test_resolve_report_slot_for_section_canonical_names(
    section_name: str, expected_slot: str
):
    slot = resolve_report_slot_for_section(section_name)
    assert slot is not None
    assert slot.slot_id == expected_slot


def test_major_level_narrative_tissue_expression():
    markdown = (
        "**Tissue and cell expression** (SREBF2)\n\n"
        "Expression narrative from synthesis.\n\n"
        "**Key findings**\n\n- Cortex reported.\n"
    )
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[],
        report_sections=[
            _report_section(
                section_name="Tissue and cell expression",
                content_markdown=markdown,
                source_ids=["sid-synth-expr"],
                status="llm",
            )
        ],
    )
    major2 = next(s for s in doc.sections if s.key == "2")
    assert major2.narrative_markdown is not None
    assert "Expression narrative from synthesis." in major2.narrative_markdown
    assert major2.synthesis_status == "llm"
    for sub in major2.subsections:
        assert sub.narrative_markdown is None
        assert sub.synthesis_status is None
    assert doc.unmapped_report_sections == []


def test_subsection_narrative_pathways():
    markdown = "**Pathways**\n\nPathway narrative.\n"
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[],
        report_sections=[
            _report_section(
                section_name="Pathways",
                content_markdown=markdown,
                source_ids=["sid-synth-path"],
                status="llm",
            )
        ],
    )
    major10 = next(s for s in doc.sections if s.key == "10")
    sub_a = next(s for s in major10.subsections if s.key == "a")
    assert sub_a.narrative_markdown is not None
    assert "Pathway narrative." in sub_a.narrative_markdown
    assert sub_a.synthesis_status == "llm"
    assert major10.narrative_markdown is None
    assert major10.synthesis_status is None


def test_evidence_preserved_with_major_level_synthesis():
    evidence = _ev(
        source_name="GTEx",
        assertion_type=AssertionType.expression,
        key="median",
        display_text="GTEx median expression in cortex.",
        section="Tissue expression",
    )
    synth_ids = ["sid-synth-a", "sid-synth-b"]
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[evidence],
        report_sections=[
            _report_section(
                section_name="Tissue and cell expression",
                content_markdown="**Tissue**\n\nSynthesized tissue narrative.\n",
                source_ids=synth_ids,
                status="llm",
            )
        ],
    )
    major2 = next(s for s in doc.sections if s.key == "2")
    sub_2a = next(s for s in major2.subsections if s.key == "a")

    assert major2.narrative_markdown is not None
    assert "Synthesized tissue narrative." in major2.narrative_markdown
    assert major2.synthesis_status == "llm"
    for sid in synth_ids:
        assert sid in major2.source_ids

    assert sub_2a.status == "populated"
    assert any((b.text or "").startswith("GTEx median") for b in sub_2a.blocks)
    assert evidence.source_id in sub_2a.source_ids
    assert evidence.id in sub_2a.evidence_record_ids
    for sid in synth_ids:
        assert sid not in sub_2a.source_ids
        assert sid not in sub_2a.evidence_record_ids
        assert sid not in major2.evidence_record_ids


def test_source_id_order_preserving_dedupe():
    evidence = _ev(
        source_name="Reactome",
        assertion_type=AssertionType.pathway_membership,
        key="Q12772:pathway",
        display_text="Cholesterol biosynthesis pathway.",
        section="Pathways",
    )
    # Distinct IDs with intentional overlap/repeat against the evidence source_id.
    wiki_id = "wikipathways:SREBF2:pathway"
    synth_ids = [
        evidence.source_id,
        wiki_id,
        wiki_id,
        "sid-synth-extra",
    ]
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[evidence],
        report_sections=[
            _report_section(
                section_name="Pathways",
                content_markdown="**Pathways**\n\nPathway prose.\n",
                source_ids=synth_ids,
                status="llm",
            )
        ],
    )
    major10 = next(s for s in doc.sections if s.key == "10")
    sub_a = next(s for s in major10.subsections if s.key == "a")

    expected = [evidence.source_id, wiki_id, "sid-synth-extra"]
    assert sub_a.source_ids == expected
    assert major10.source_ids == expected
    assert sub_a.evidence_record_ids == [evidence.id]
    assert wiki_id not in sub_a.evidence_record_ids
    assert "sid-synth-extra" not in major10.evidence_record_ids


def test_report_sections_omitted_equals_empty_list():
    records = [
        _ev(
            source_name="NCBI Gene",
            assertion_type=AssertionType.gene_identity,
            key="6721",
            display_text="SREBF2 Entrez Gene ID is 6721.",
            grade=EvidenceGrade.A,
        )
    ]
    doc_omitted = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=records,
        chromosome="22",
    )
    doc_empty = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=records,
        report_sections=[],
        chromosome="22",
    )
    assert doc_omitted.model_dump(mode="json") == doc_empty.model_dump(mode="json")
    sec1 = doc_omitted.sections[0]
    aliases = next(s for s in sec1.subsections if s.key == "a")
    assert aliases.status == "populated"
    assert aliases.blocks
    assert aliases.narrative_markdown is None
    assert aliases.synthesis_status is None


def test_meta_sections_retained():
    sections = [
        _report_section(
            section_name="Missing / deferred / manual sources",
            content_markdown="_Deferred meta note._",
            status="deferred",
            source_ids=["sid-meta-1"],
        ),
        _report_section(
            section_name="Verification warnings",
            content_markdown="_Verification meta note._",
            status="deferred",
            source_ids=["sid-meta-2"],
        ),
    ]
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[],
        report_sections=sections,
    )
    assert len(doc.unmapped_report_sections) == 2
    by_name = {u.section_name: u for u in doc.unmapped_report_sections}
    for name, md, sid in [
        ("Missing / deferred / manual sources", "_Deferred meta note._", "sid-meta-1"),
        ("Verification warnings", "_Verification meta note._", "sid-meta-2"),
    ]:
        entry = by_name[name]
        assert entry.reason == "meta_section"
        assert entry.attempted_slot is None
        assert entry.content_markdown == md
        assert entry.status == "deferred"
        assert entry.source_ids == [sid]
    for major in doc.sections:
        assert major.narrative_markdown is None
        assert major.synthesis_status is None
        for sub in major.subsections:
            assert sub.narrative_markdown is None
            assert sub.synthesis_status is None


def test_unknown_section_retained():
    unknown = _report_section(
        section_name="Completely Unknown Section",
        content_markdown="Should not land on Rancho body.",
        source_ids=["sid-unknown"],
        status="llm",
    )
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[],
        report_sections=[unknown],
    )
    assert len(doc.unmapped_report_sections) == 1
    entry = doc.unmapped_report_sections[0]
    assert entry.reason == "unmapped_name"
    assert entry.section_name == "Completely Unknown Section"
    assert entry.content_markdown == "Should not land on Rancho body."
    assert entry.source_ids == ["sid-unknown"]
    assert entry.status == "llm"
    for major in doc.sections:
        assert major.narrative_markdown is None
        for sub in major.subsections:
            assert sub.narrative_markdown is None


def test_slot_collision_first_writer_wins():
    first = _report_section(
        section_name="Tissue and cell expression",
        content_markdown="First narrative wins.",
        source_ids=["sid-first"],
        status="llm",
    )
    second = _report_section(
        section_name="Tissue and cell expression",
        content_markdown="Second narrative must not overwrite.",
        source_ids=["sid-second"],
        status="deterministic",
    )
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[],
        report_sections=[first, second],
    )
    major2 = next(s for s in doc.sections if s.key == "2")
    assert major2.narrative_markdown is not None
    assert "First narrative wins." in major2.narrative_markdown
    assert "Second narrative" not in major2.narrative_markdown
    assert major2.synthesis_status == "llm"
    assert "sid-first" in major2.source_ids
    assert "sid-second" not in major2.source_ids

    conflicts = [u for u in doc.unmapped_report_sections if u.reason == "slot_conflict"]
    assert len(conflicts) == 1
    conflict = conflicts[0]
    assert conflict.attempted_slot == "2"
    assert conflict.content_markdown == "Second narrative must not overwrite."
    assert conflict.source_ids == ["sid-second"]
    assert conflict.status == "deterministic"


@pytest.mark.parametrize("status", ["empty", "deferred"])
def test_blank_narrative_remains_truthful(status: str):
    evidence = _ev(
        source_name="GTEx",
        assertion_type=AssertionType.expression,
        key=f"median-blank-{status}",
        display_text="GTEx still present.",
        section="Tissue expression",
    )
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[evidence],
        report_sections=[
            _report_section(
                section_name="Tissue and cell expression",
                content_markdown="   ",
                source_ids=["sid-empty"],
                status=status,
            )
        ],
    )
    major2 = next(s for s in doc.sections if s.key == "2")
    sub_2a = next(s for s in major2.subsections if s.key == "a")
    assert major2.narrative_markdown is None
    assert major2.synthesis_status == status
    assert sub_2a.status == "populated"
    assert any("GTEx still present." in (b.text or "") for b in sub_2a.blocks)


def test_all_canonical_chdi_sections_together():
    assert len(CANONICAL_SYNTHESIS_SECTION_NAMES) == 22
    report_sections = [
        _report_section(
            section_name=name,
            content_markdown=f"**{name}**\n\nProse for {name}.\n",
            source_ids=[f"sid-{i}"],
            status="llm",
        )
        for i, name in enumerate(CANONICAL_SYNTHESIS_SECTION_NAMES)
    ]
    doc = build_report_document(
        dossier_run_id="run-schema",
        gene_symbol="SREBF2",
        evidence_records=[],
        report_sections=report_sections,
    )

    mapped = 0
    for major in doc.sections:
        if major.synthesis_status is not None:
            mapped += 1
        for sub in major.subsections:
            if sub.synthesis_status is not None:
                mapped += 1
    assert mapped == 20

    reasons = [u.reason for u in doc.unmapped_report_sections]
    assert reasons.count("meta_section") == 2
    assert reasons.count("slot_conflict") == 0
    assert reasons.count("unmapped_name") == 0
    assert len(doc.unmapped_report_sections) == 2

