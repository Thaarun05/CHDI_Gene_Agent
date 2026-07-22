"""End-to-end tests for the Rancho/CHDI HTML+PDF renderer (no network / no LLM)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.rancho_report import (
    build_and_write_rancho_report,
    render_rancho_html,
    write_rancho_report,
)
from gene_dossier.report_schema import (
    REPORT_SECTIONS,
    REPORT_STYLE,
    ReportContentBlock,
    ReportCover,
    ReportDocument,
    ReportMajorSection,
    ReportSubsection,
    build_report_document,
)
from gene_dossier.source_ids import make_source_id


# --------------------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------------------
def _fake_report_document() -> ReportDocument:
    """Minimal hand-built ReportDocument (not built from evidence)."""
    sections: list[ReportMajorSection] = []
    for spec in REPORT_SECTIONS:
        subs = [
            ReportSubsection(
                key=sub.key,
                title=sub.title,
                toc_title=sub.toc_title,
                status="empty",
            )
            for sub in spec.subsections
        ]
        if spec.number == 1 and subs:
            subs[0].status = "populated"
            subs[0].blocks = [
                ReportContentBlock(
                    kind="narrative",
                    text="Fake gene aliases: SREBF2, SREBP2.",
                    source_ids=["src:fake:aliases"],
                ),
                ReportContentBlock(
                    kind="table",
                    title="Identifier table",
                    table_headers=["Source", "ID"],
                    table_rows=[["NCBI Gene", "6721"], ["Ensembl", "ENSG00000198911"]],
                    source_ids=["src:fake:ids"],
                ),
            ]
            subs[0].source_ids = ["src:fake:aliases", "src:fake:ids"]
        sections.append(
            ReportMajorSection(
                number=spec.number,
                key=spec.key,
                title=spec.title,
                toc_title=spec.toc_title,
                subsections=subs,
                status="partial" if spec.number == 1 else "empty",
                source_ids=subs[0].source_ids if spec.number == 1 and subs else [],
            )
        )
    return ReportDocument(
        dossier_run_id="fake-run",
        gene_symbol="FAKEG",
        cover=ReportCover(
            gene_symbol="FAKEG",
            chromosome="1",
            curator="Test Curator",
            report_date="2026-07-22",
        ),
        sections=sections,
        references=["Doe J et al. Fake Journal. 2024."],
        compiled_databases=[
            {"name": "GeneCards", "url": "https://www.genecards.org/"},
        ],
        unmapped_source_ids=[],
    )


def _srebf2_evidence() -> list[EvidenceRecord]:
    """Validated SREBF2 anchors shaped as EvidenceRecords (offline, no network)."""
    run_id = "srebf2-e2e"
    gene = "SREBF2"

    def ev(
        source_name: str,
        assertion: AssertionType,
        key: str,
        text: str,
        *,
        section: str = "General Gene Information",
        subsection: str | None = None,
        fact_type: str = "fact",
        grade: EvidenceGrade = EvidenceGrade.C,
        source_type: SourceType = SourceType.curated_database,
        value: dict | None = None,
    ) -> EvidenceRecord:
        return EvidenceRecord(
            source_id=make_source_id(source_name, gene, assertion, key),
            dossier_run_id=run_id,
            gene_symbol=gene,
            official_symbol=gene,
            section=section,
            subsection=subsection,
            source_name=source_name,
            source_type=source_type,
            assertion_type=assertion,
            fact_type=fact_type,
            evidence_grade=grade,
            display_text=text,
            value=value or {},
            organism="Homo sapiens",
            taxon_id=9606,
        )

    return [
        ev(
            "NCBI Gene",
            AssertionType.gene_identity,
            "6721",
            "SREBF2 Entrez Gene ID is 6721 (chr22).",
            grade=EvidenceGrade.A,
            value={"chromosome": "22", "entrez_gene_id": "6721"},
        ),
        ev(
            "Ensembl",
            AssertionType.gene_identity,
            "ENSG00000198911",
            "Ensembl gene ID for SREBF2 is ENSG00000198911.",
            grade=EvidenceGrade.A,
        ),
        ev(
            "UCSC",
            AssertionType.gene_identity,
            "chr22:41833105-41907305",
            "UCSC hg38 region for SREBF2: chr22:41833105-41907305 "
            "(canonical transcript ENST00000361204.9).",
        ),
        ev(
            "UniProt",
            AssertionType.protein_function,
            "Q12772",
            "UniProt Q12772: Sterol regulatory element-binding protein 2 (SREBP-2).",
            section="Known structure",
        ),
        ev(
            "AlphaFold",
            AssertionType.protein_structure,
            "Q12772",
            "AlphaFold structure prediction available for UniProt Q12772.",
            section="AlphaFold",
        ),
        ev(
            "NCBI Datasets",
            AssertionType.orthology,
            "20788",
            "Mouse ortholog of SREBF2 is Entrez 20788 (MGI:107585).",
            section="Homologues",
        ),
        ev(
            "GTEx",
            AssertionType.expression,
            "ENSG00000198911.11",
            "GTEx median expression for SREBF2 (ENSG00000198911.11) is high in brain.",
            section="Tissue and cell expression",
            source_type=SourceType.expression_database,
            grade=EvidenceGrade.B,
        ),
        ev(
            "Allen Brain Atlas",
            AssertionType.expression,
            "1051154",
            "Allen HBA probe 1051154 detects SREBF2 in human brain.",
            section="Brain expression",
            source_type=SourceType.expression_database,
            grade=EvidenceGrade.B,
        ),
        ev(
            "DropViz",
            AssertionType.cell_type_expression,
            "neuron",
            "DropViz: SREBF2 is detected in neuronal clusters.",
            section="snRNA-seq cell type",
            source_type=SourceType.expression_database,
            grade=EvidenceGrade.B,
        ),
        ev(
            "GEO",
            AssertionType.perturbation,
            "GDS5046",
            "GEO Profiles: SREBF2 altered under neuronal perturbation (GDS5046).",
            section="GEO perturbation",
        ),
        ev(
            "Harmonizome",
            AssertionType.transcription_factor_association,
            "SP1",
            "Harmonizome: SP1 is associated with SREBF2 regulation.",
            section="Transcription factor",
        ),
        ev(
            "STRING",
            AssertionType.ppi,
            "SCAP",
            "STRING: SREBF2 interacts with SCAP.",
            section="Protein-protein interaction",
            source_type=SourceType.interaction_database,
        ),
        ev(
            "BioGRID",
            AssertionType.ppi,
            "SCAP",
            "BioGRID: SREBF2–SCAP physical interaction reported.",
            section="Protein-protein interaction",
            source_type=SourceType.interaction_database,
        ),
        ev(
            "CTD",
            AssertionType.chemical_interaction,
            "D000077185",
            "CTD: SREBF2 is annotated with chemical interactions including statins.",
            section="CTD perturbation",
            source_type=SourceType.chemical_database,
        ),
        ev(
            "ChEMBL",
            AssertionType.chemical_interaction,
            "CHEMBL2363040",
            "ChEMBL: chemical tools associated with SREBP2 biology.",
            section="Chemical tools",
            source_type=SourceType.chemical_database,
        ),
        ev(
            "GTEx",
            AssertionType.expression,
            "eqtl-brain",
            "GTEx brain eQTL signals alter SREBF2 expression.",
            section="eQTL brain tissue",
            fact_type="eqtl",
            source_type=SourceType.expression_database,
            grade=EvidenceGrade.B,
        ),
        ev(
            "ClinVar",
            AssertionType.variant_association,
            "VCV000012345",
            "ClinVar: SREBF2 variants with reported clinical significance.",
            section="ClinVar",
            source_type=SourceType.genetic_database,
            grade=EvidenceGrade.B,
        ),
        ev(
            "OMIM",
            AssertionType.disease_association,
            "600481",
            "OMIM MIM number for SREBF2 is 600481.",
            section="OMIM",
            source_type=SourceType.genetic_database,
            grade=EvidenceGrade.A,
        ),
        ev(
            "Reactome",
            AssertionType.pathway_membership,
            "R-HSA-1655829",
            "Reactome: SREBF2 participates in regulation of cholesterol biosynthesis.",
            section="Pathways",
            source_type=SourceType.pathway_database,
        ),
        ev(
            "MGI",
            AssertionType.knockout_phenotype,
            "MGI:107585",
            "MGI:107585 knockout phenotypes reported for Srebf2.",
            section="Knockout",
            source_type=SourceType.model_organism_database,
            grade=EvidenceGrade.D,
        ),
        ev(
            "PubMed",
            AssertionType.literature_summary,
            "PMID:12345678",
            "Major labs publish on SREBF2 lipid metabolism and neurodegeneration.",
            section="Major labs",
            source_type=SourceType.literature,
            grade=EvidenceGrade.E,
        ),
        ev(
            "Antibodies",
            AssertionType.literature_summary,
            "ab30682",
            "Commercial antibodies available against SREBP2 (e.g. Abcam ab30682).",
            section="Commercial antibodies",
            source_type=SourceType.commercial_source,
            grade=EvidenceGrade.F,
        ),
        ev(
            "Patents",
            AssertionType.patent_claim,
            "US2010123456",
            "Patent literature references SREBF2 / SREBP2 as a therapeutic target.",
            section="Patents",
            source_type=SourceType.patent_database,
            grade=EvidenceGrade.E,
        ),
        ev(
            "NIH RePORTER",
            AssertionType.grant_project,
            "R01HL123456",
            "NIH RePORTER lists funded projects mentioning SREBF2.",
            section="NIH grants",
            source_type=SourceType.grant_database,
            grade=EvidenceGrade.E,
        ),
        ev(
            "ERC",
            AssertionType.grant_project,
            "ERC-2020-AdG",
            "ERC projects include SREBF2-related lipid biology aims.",
            section="ERC grants",
            source_type=SourceType.grant_database,
            grade=EvidenceGrade.E,
        ),
        # Intentionally unmapped — should appear only in JSON sidecar provenance.
        ev(
            "Mystery Widget",
            AssertionType.visual_observation,
            "widget-1",
            "Unmapped visual note that should not invent a new major section.",
            section="misc",
            source_type=SourceType.manual_visual_reference,
            grade=EvidenceGrade.F,
        ),
    ]


# --------------------------------------------------------------------------------------
# Fake ReportDocument → HTML / files
# --------------------------------------------------------------------------------------
def test_render_fake_document_html_matches_rancho_chrome():
    doc = _fake_report_document()
    html = render_rancho_html(doc, toc_page_numbers={"1": 4, "1a": 4})

    assert "FAKEG (CHR1)" in html
    assert "Gene Report" in html
    assert "Prepared for the CHDI Foundation" in html
    assert "Table of Contents" in html
    assert "toc-leader" in html
    assert ">4<" in html
    assert "1. General Gene Information" in html
    assert "Fake gene aliases" in html
    assert "Identifier table" in html
    assert "bgcolor=" in html  # light-green table header for Story/PDF
    assert REPORT_STYLE.table_header_bg in html
    assert REPORT_STYLE.green_major in html
    assert REPORT_STYLE.orange_sub in html
    assert "References" in html
    assert "Compiled List of Relevant Databases" in html
    assert "GeneCards" in html
    # Defaults: no provenance endnotes, no cover logos
    assert "Provenance endnotes" not in html
    assert 'class="cover-logos"' not in html
    assert "var(--" not in html


def test_write_fake_document_html_json_and_optional_pdf(tmp_path: Path):
    doc = _fake_report_document()
    paths = write_rancho_report(doc, output_dir=tmp_path, write_pdf=True)

    assert paths["html"].is_file()
    assert paths["json"].is_file()
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert payload["gene_symbol"] == "FAKEG"
    assert payload["cover"]["gene_symbol"] == "FAKEG"
    assert len(payload["sections"]) == 15
    # source_ids preserved in sidecar even without endnotes
    assert payload["sections"][0]["subsections"][0]["source_ids"]

    if "pdf" in paths:
        assert paths["pdf"].is_file()
        assert paths["pdf"].stat().st_size > 1000
    else:
        pytest.importorskip("fitz")


# --------------------------------------------------------------------------------------
# Real SREBF2 evidence → schema → HTML/PDF
# --------------------------------------------------------------------------------------
def test_build_document_from_srebf2_evidence_fills_all_15_sections():
    records = _srebf2_evidence()
    doc = build_report_document(
        dossier_run_id="srebf2-e2e",
        gene_symbol="SREBF2",
        evidence_records=records,
        curator="Rancho BioSciences",
        report_date="2026-07-22",
        references=["Horton JD et al. J Clin Invest. 2002."],
    )

    assert doc.cover.gene_line == "SREBF2 (CHR22)"
    assert len(doc.sections) == 15
    assert doc.unmapped_source_ids  # Mystery Widget
    mystery = next(
        r for r in records if r.source_name == "Mystery Widget"
    )
    assert mystery.source_id in doc.unmapped_source_ids

    populated = [s for s in doc.sections if s.status in {"populated", "partial"}]
    assert len(populated) == 15, [s.status for s in doc.sections]

    # Spot-check key SREBF2 anchors landed in expected slots
    aliases = doc.sections[0].subsections[0]
    assert any("6721" in (b.text or "") for b in aliases.blocks)
    assert any("ENSG00000198911" in (b.text or "") for b in aliases.blocks)

    string_sub = doc.sections[4].subsections[0]
    assert any("SCAP" in (b.text or "") for b in string_sub.blocks)

    omim = doc.sections[8].subsections[1]
    assert any("600481" in (b.text or "") for b in omim.blocks)


def test_build_and_write_srebf2_rancho_report_end_to_end(tmp_path: Path):
    records = _srebf2_evidence()
    doc, paths = build_and_write_rancho_report(
        dossier_run_id="srebf2-e2e",
        gene_symbol="SREBF2",
        evidence_records=records,
        curator="Rancho BioSciences",
        report_date="2026-07-22",
        chromosome="22",
        references=["Horton JD et al. J Clin Invest. 2002."],
        output_dir=tmp_path,
        include_endnotes=False,
        write_pdf=True,
        toc_page_numbers={"1": 4, "1a": 4, "5": 12},
        stamp_cover=False,
        show_cover_logos=False,
    )

    assert doc.gene_symbol == "SREBF2"
    assert paths["html"].name == "srebf2-e2e_rancho_report.html"
    assert paths["json"].name == "srebf2-e2e_rancho_report.json"

    html = paths["html"].read_text(encoding="utf-8")
    assert "SREBF2 (CHR22)" in html
    assert "Entrez Gene ID is 6721" in html
    assert "ENSG00000198911" in html
    assert "Q12772" in html
    assert "SCAP" in html
    assert "600481" in html
    assert "Horton JD" in html
    assert "Compiled List of Relevant Databases" in html
    assert "Provenance endnotes" not in html
    assert 'class="cover-logos"' not in html
    assert "toc-leader" in html
    assert REPORT_STYLE.green_major in html

    sidecar = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert sidecar["dossier_run_id"] == "srebf2-e2e"
    assert len(sidecar["sections"]) == 15
    assert sidecar["unmapped_source_ids"]
    # Provenance retained in JSON even when endnotes are off
    sec1a = sidecar["sections"][0]["subsections"][0]
    assert sec1a["source_ids"]
    assert sec1a["blocks"][0]["source_ids"]

    if "pdf" in paths:
        pdf = paths["pdf"]
        assert pdf.is_file() and pdf.stat().st_size > 2000
        fitz = pytest.importorskip("fitz")
        with fitz.open(pdf) as pdf_doc:
            assert pdf_doc.page_count >= 2
            # Body pages should carry stamped chrome (URL + page number)
            body = pdf_doc[1]
            text = body.get_text()
            assert "RanchoBioSciences" in text or "www.RanchoBioSciences.com" in text
    else:
        pytest.importorskip("fitz")


def test_optional_endnotes_and_cover_logos_can_be_enabled(tmp_path: Path):
    records = _srebf2_evidence()[:3]
    _doc, paths = build_and_write_rancho_report(
        dossier_run_id="srebf2-opts",
        gene_symbol="SREBF2",
        evidence_records=records,
        output_dir=tmp_path,
        include_endnotes=True,
        show_cover_logos=True,
        write_pdf=False,
    )
    html = paths["html"].read_text(encoding="utf-8")
    assert "Provenance endnotes" in html
    assert 'class="cover-logos"' in html
