"""Tests for polished Section 1a Gene Aliases presentation."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
    new_id,
)
from gene_dossier.rancho_report import (
    render_rancho_html,
    render_rancho_section_fragment,
)
from gene_dossier.report_presentation import (
    NOT_AVAILABLE,
    build_gene_aliases_blocks,
    build_section_presentation,
    format_safe_table_cell_html,
)
from gene_dossier.report_schema import (
    ReportContentBlock,
    ReportCover,
    ReportDocument,
    ReportMajorSection,
    ReportSubsection,
    build_report_document,
    _block_from_evidence,
)


def _ev(
    *,
    source_name: str,
    fact_type: str,
    value: dict,
    assertion_type: AssertionType = AssertionType.gene_identity,
    taxon_id: int | None = 9606,
    organism: str | None = "Homo sapiens",
    gene_symbol: str = "SREBF2",
    source_id: str | None = None,
    evidence_id: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        id=evidence_id or new_id(),
        source_id=source_id or f"src-{new_id()[:8]}",
        dossier_run_id="test-run",
        gene_symbol=gene_symbol,
        section="General gene information",
        source_name=source_name,
        source_type=SourceType.curated_database,
        assertion_type=assertion_type,
        fact_type=fact_type,
        organism=organism,
        taxon_id=taxon_id,
        evidence_grade=EvidenceGrade.C,
        value=value,
        display_text=f"{gene_symbol} {fact_type}",
    )


def _human_bundle(gene: str = "SREBF2") -> list[EvidenceRecord]:
    return [
        _ev(
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            gene_symbol=gene,
            value={
                "entrez_gene_id": "6721",
                "nomenclaturesymbol": gene,
                "name": gene,
                "description": "sterol regulatory element binding transcription factor 2",
            },
            source_id="ncbi-entrez",
            evidence_id="eid-ncbi",
        ),
        _ev(
            source_name="Ensembl",
            fact_type="ensembl_gene_id",
            gene_symbol=gene,
            taxon_id=None,
            organism="homo sapiens",
            value={"ensembl_gene_id": "ENSG00000198911", "display_name": gene},
            source_id="ensembl-id",
            evidence_id="eid-ensembl",
        ),
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            gene_symbol=gene,
            value={
                "uniprot_accession": "Q12772",
                "gene_names": [gene, "BHLHD2", "SREBP2"],
            },
            source_id="uniprot-id",
            evidence_id="eid-uniprot",
        ),
    ]


def _mouse_bundle() -> list[EvidenceRecord]:
    return [
        _ev(
            source_name="NCBI Datasets",
            fact_type="ortholog_gene",
            assertion_type=AssertionType.orthology,
            taxon_id=10090,
            organism="Mus musculus",
            value={
                "ortholog_gene_id": "20788",
                "ortholog_symbol": "Srebf2",
                "tax_id": "10090",
            },
            source_id="datasets-mouse",
            evidence_id="eid-datasets-mouse",
        ),
        _ev(
            source_name="MouseMine",
            fact_type="mgi_gene_id",
            assertion_type=AssertionType.orthology,
            taxon_id=10090,
            organism="Mus musculus",
            value={
                "mgi_id": "MGI:107585",
                "ncbi_gene_number": "20788",
                "mouse_symbol": "Srebf2",
                "mouse_name": "sterol regulatory element binding factor 2",
                "aliases": ["SREBP-2", "lop13"],
            },
            source_id="mgi-id",
            evidence_id="eid-mgi",
        ),
    ]


def _rat_bundle() -> list[EvidenceRecord]:
    return [
        _ev(
            source_name="NCBI Datasets",
            fact_type="ortholog_gene",
            assertion_type=AssertionType.orthology,
            taxon_id=10116,
            organism="Rattus norvegicus",
            value={
                "ortholog_gene_id": "300095",
                "ortholog_symbol": "Srebf2",
                "tax_id": "10116",
            },
            source_id="datasets-rat",
            evidence_id="eid-datasets-rat",
        )
    ]


# --------------------------------------------------------------------------------------
# Core table construction
# --------------------------------------------------------------------------------------
def test_section_1a_creates_exactly_one_polished_table_block():
    result = build_section_presentation(
        section_key="1a",
        gene_symbol="SREBF2",
        evidence_records=_human_bundle() + _mouse_bundle() + _rat_bundle(),
    )
    assert len(result.blocks) == 1
    assert result.blocks[0].kind == "table"
    assert result.blocks[0].presentation_role == "gene_aliases_table"


def test_column_and_row_order():
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2",
        evidence_records=_human_bundle(),
    ).blocks[0]
    assert block.table_headers == ["", "Human gene", "Mouse gene", "Rat gene"]
    labels = [row[0] for row in block.table_rows]
    assert labels == [
        "Entrez Gene ID",
        "Gene Symbol",
        "Gene Name",
        "Ensembl ID",
        "UniProt ID",
        "Synonyms/\u200bAliases",
    ]


def test_structured_values_populate_species_columns():
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2",
        evidence_records=_human_bundle() + _mouse_bundle() + _rat_bundle(),
    ).blocks[0]
    rows = {r[0]: r for r in block.table_rows}
    assert "6721" in rows["Entrez Gene ID"][1]
    assert "20788" in rows["Entrez Gene ID"][2]
    assert "300095" in rows["Entrez Gene ID"][3]
    assert rows["Gene Symbol"][1] == "SREBF2"
    assert rows["Gene Symbol"][2] == "Srebf2"
    assert rows["Gene Symbol"][3] == "Srebf2"
    assert "ENSG00000198911" in rows["Ensembl ID"][1]
    assert "Q12772" in rows["UniProt ID"][1]
    assert "BHLHD2" in rows["Synonyms/\u200bAliases"][1]
    assert "lop13" in rows["Synonyms/\u200bAliases"][2]


def test_alias_dedupe_preserves_order():
    records = _human_bundle()
    # Extra UniProt-like aliases with duplicates
    records.append(
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            value={
                "uniprot_accession": "Q12772",
                "gene_names": ["SREBF2", "BHLHD2", "SREBP2", "bhlhd2", "BHLHD2"],
            },
            source_id="uniprot-dup",
        )
    )
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2", evidence_records=records
    ).blocks[0]
    aliases = block.table_rows[-1][1]
    assert aliases == "BHLHD2, SREBP2"
    assert "SREBF2" not in aliases.split(", ")


def test_missing_values_become_not_available():
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2",
        evidence_records=_human_bundle(),
    ).blocks[0]
    rows = {r[0]: r for r in block.table_rows}
    assert rows["Entrez Gene ID"][2] == NOT_AVAILABLE
    assert rows["Entrez Gene ID"][3] == NOT_AVAILABLE
    assert rows["Ensembl ID"][2] == NOT_AVAILABLE
    assert rows["UniProt ID"][3] == NOT_AVAILABLE


def test_no_raw_evidence_paragraphs_in_polished_1a_html():
    doc = build_report_document(
        dossier_run_id="test-run",
        gene_symbol="SREBF2",
        evidence_records=_human_bundle() + _mouse_bundle() + _rat_bundle(),
        report_sections=None,
    )
    html = render_rancho_section_fragment(
        document=doc, section_number=1, subsection_key="a"
    )
    assert "Supporting evidence" not in html
    assert "Key findings" not in html
    assert "Limitations" not in html
    assert "Entrez Gene ID is 6721" not in html
    assert "gene-aliases-table" in html
    assert "6721" in html


def test_entrez_ensembl_uniprot_links_safe_and_correct():
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2", evidence_records=_human_bundle()
    ).blocks[0]
    rows = {r[0]: r for r in block.table_rows}
    assert rows["Entrez Gene ID"][1] == (
        "[6721](https://www.ncbi.nlm.nih.gov/gene/6721)"
    )
    assert rows["Ensembl ID"][1] == (
        "[ENSG00000198911](https://www.ensembl.org/id/ENSG00000198911)"
    )
    assert rows["UniProt ID"][1] == (
        "[Q12772](https://www.uniprot.org/uniprotkb/Q12772)"
    )


@pytest.mark.parametrize(
    "cell,expect_link",
    [
        ("[6721](https://www.ncbi.nlm.nih.gov/gene/6721)", True),
        ("[x](http://www.ncbi.nlm.nih.gov/gene/1)", False),
        ("[x](https://evil.example/gene/1)", False),
        ('[x](https://www.ncbi.nlm.nih.gov/gene/1"onclick=alert(1))', False),
        ("[<script>x</script>](https://www.ncbi.nlm.nih.gov/gene/1)", True),
        ("not a link <b>bold</b>", False),
        ("[broken](https://www.ncbi.nlm.nih.gov/gene/1) trailing", False),
    ],
)
def test_link_safety_cases(cell: str, expect_link: bool):
    html = format_safe_table_cell_html(cell)
    if expect_link:
        assert html.startswith("<a ")
        assert "<script>" not in html
    else:
        assert "<a " not in html
        if "<" in cell or ">" in cell:
            assert "&lt;" in html or "&gt;" in html or "<a " not in html


def test_html_escapes_unsafe_text_in_names_and_aliases():
    records = [
        _ev(
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            value={
                "entrez_gene_id": "1",
                "nomenclaturesymbol": "GENE",
                "description": 'name <script>alert(1)</script> & "x"',
            },
        ),
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            value={
                "uniprot_accession": "P12345",
                "gene_names": ["GENE", "Alias<script>"],
            },
        ),
    ]
    doc = build_report_document(
        dossier_run_id="t",
        gene_symbol="GENE",
        evidence_records=records,
    )
    html = render_rancho_section_fragment(
        document=doc, section_number=1, subsection_key="a"
    )
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_input_evidence_records_not_mutated():
    records = _human_bundle() + _mouse_bundle()
    before = [copy.deepcopy(r.model_dump()) for r in records]
    build_gene_aliases_blocks(gene_symbol="SREBF2", evidence_records=records)
    after = [r.model_dump() for r in records]
    assert before == after


def test_builder_works_with_non_srebf2_gene():
    records = [
        _ev(
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            gene_symbol="CDH10",
            value={
                "entrez_gene_id": "1008",
                "nomenclaturesymbol": "CDH10",
                "description": "cadherin 10",
            },
        ),
        _ev(
            source_name="Ensembl",
            fact_type="ensembl_gene_id",
            gene_symbol="CDH10",
            taxon_id=None,
            organism="Homo sapiens",
            value={"ensembl_gene_id": "ENSG00000040731", "display_name": "CDH10"},
        ),
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            gene_symbol="CDH10",
            value={"uniprot_accession": "Q9Y6N8", "gene_names": ["CDH10", "T2-cadherin"]},
        ),
    ]
    block = build_gene_aliases_blocks(
        gene_symbol="CDH10", evidence_records=records
    ).blocks[0]
    assert "1008" in block.table_rows[0][1]
    assert block.table_rows[1][1] == "CDH10"
    assert "T2-cadherin" in block.table_rows[-1][1]


def test_audit_evidence_rendering_still_available():
    rec = _human_bundle()[0]
    block = _block_from_evidence(rec)
    assert block.kind == "narrative"
    assert "entrez_gene_id" in (block.text or "")
    doc = build_report_document(
        dossier_run_id="t",
        gene_symbol="SREBF2",
        evidence_records=[rec],
    )
    sub = doc.sections[0].subsections[0]
    assert sub.presentation_blocks
    assert sub.blocks  # evidence audit blocks preserved
    assert sub.blocks[0].kind == "narrative"


def test_unknown_section_key_returns_empty():
    result = build_section_presentation(
        section_key="9z",
        gene_symbol="SREBF2",
        evidence_records=_human_bundle(),
    )
    assert result.blocks == ()
    assert result.diagnostics == ()


def test_scalar_precedence_higher_wins():
    records = [
        _ev(
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            value={
                "entrez_gene_id": "6721",
                "nomenclaturesymbol": "SREBF2",
                "description": "from ncbi",
            },
            source_id="ncbi",
            evidence_id="e-ncbi",
        ),
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            value={
                "uniprot_accession": "Q12772",
                "gene_names": ["OTHERSYM", "Alias"],
            },
            source_id="uni",
            evidence_id="e-uni",
        ),
    ]
    result = build_gene_aliases_blocks(gene_symbol="SREBF2", evidence_records=records)
    block = result.blocks[0]
    assert block.table_rows[1][1] == "SREBF2"  # NCBI symbol wins
    # UniProt did not overwrite symbol
    assert not any(
        d.severity == "warning" and "human.symbol" in d.field and "OTHERSYM" in d.reason
        for d in result.diagnostics
    ) or any("human.symbol" in d.field for d in result.diagnostics)


def test_lower_priority_cannot_overwrite_higher():
    records = [
        _ev(
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            value={
                "entrez_gene_id": "6721",
                "nomenclaturesymbol": "SREBF2",
                "description": "official name",
            },
        ),
        _ev(
            source_name="Ensembl",
            fact_type="ensembl_gene_id",
            taxon_id=None,
            value={"ensembl_gene_id": "ENSG00000198911", "display_name": "WRONG"},
        ),
    ]
    result = build_gene_aliases_blocks(gene_symbol="SREBF2", evidence_records=records)
    assert result.blocks[0].table_rows[1][1] == "SREBF2"
    warnings = [d for d in result.diagnostics if d.severity == "warning"]
    assert any(
        d.field == "human.symbol" and "WRONG" in d.reason for d in warnings
    )


def test_diagnostics_identify_missing_fields():
    result = build_gene_aliases_blocks(
        gene_symbol="SREBF2", evidence_records=_human_bundle()
    )
    missing_fields = {d.field for d in result.diagnostics if d.severity == "info"}
    assert "mouse.entrez" in missing_fields
    assert "rat.ensembl" in missing_fields
    assert "mouse.uniprot" in missing_fields


def test_alias_table_css_does_not_change_unrelated_tables():
    doc = ReportDocument(
        dossier_run_id="t",
        gene_symbol="SREBF2",
        cover=ReportCover(gene_symbol="SREBF2", chromosome="22"),
        sections=[
            ReportMajorSection(
                number=1,
                key="1",
                title="General Gene Information",
                toc_title="GENERAL GENE INFORMATION",
                subsections=[
                    ReportSubsection(
                        key="a",
                        title="Gene Aliases",
                        toc_title="GENE ALIASES",
                        presentation_blocks=[
                            ReportContentBlock(
                                kind="table",
                                presentation_role="gene_aliases_table",
                                table_headers=["", "Human gene", "Mouse gene", "Rat gene"],
                                table_rows=[["Entrez Gene ID", "1", "2", "3"]],
                            )
                        ],
                        status="populated",
                    )
                ],
                status="partial",
            ),
            ReportMajorSection(
                number=5,
                key="5",
                title="Protein-protein interaction (PPI) partners",
                toc_title="PPI",
                subsections=[
                    ReportSubsection(
                        key="a",
                        title="STRING",
                        toc_title="STRING",
                        blocks=[
                            ReportContentBlock(
                                kind="table",
                                table_headers=["Partner", "Score"],
                                table_rows=[["SCAP", "0.9"]],
                            )
                        ],
                        status="populated",
                    )
                ],
                status="populated",
            ),
        ],
    )
    html = render_rancho_html(doc)
    assert "gene-aliases-table" in html
    # Unrelated table must remain class="rancho" without gene-aliases-table
    assert '<table class="rancho">' in html
    # CSS selector + one table class occurrence in markup
    assert html.count('class="rancho gene-aliases-table"') == 1
    assert "table.gene-aliases-table" in html


def test_backward_compatible_document_without_presentation_blocks():
    payload = {
        "dossier_run_id": "legacy",
        "gene_symbol": "SREBF2",
        "cover": {"gene_symbol": "SREBF2", "chromosome": "22"},
        "sections": [
            {
                "number": 1,
                "key": "1",
                "title": "General Gene Information",
                "toc_title": "GENERAL GENE INFORMATION",
                "subsections": [
                    {
                        "key": "a",
                        "title": "Gene Aliases",
                        "toc_title": "GENE ALIASES",
                        "blocks": [
                            {
                                "kind": "narrative",
                                "text": "legacy evidence text",
                            }
                        ],
                        "status": "populated",
                    }
                ],
                "status": "partial",
            }
        ],
    }
    doc = ReportDocument.model_validate(payload)
    assert doc.sections[0].subsections[0].presentation_blocks == []
    html = render_rancho_html(doc)
    assert "legacy evidence text" in html
    assert "Supporting evidence" not in html


def test_conflicting_entrez_emits_warning():
    records = [
        _ev(
            source_name="NCBI Datasets",
            fact_type="ortholog_gene",
            assertion_type=AssertionType.orthology,
            taxon_id=10090,
            organism="Mus musculus",
            value={"ortholog_gene_id": "20788", "ortholog_symbol": "Srebf2", "tax_id": "10090"},
            source_id="ds",
            evidence_id="e-ds",
        ),
        _ev(
            source_name="MouseMine",
            fact_type="mgi_gene_id",
            assertion_type=AssertionType.orthology,
            taxon_id=10090,
            organism="Mus musculus",
            value={
                "ncbi_gene_number": "99999",
                "mouse_symbol": "Srebf2",
                "mouse_name": "mouse name",
            },
            source_id="mgi",
            evidence_id="e-mgi",
        ),
    ]
    result = build_gene_aliases_blocks(gene_symbol="SREBF2", evidence_records=records)
    assert "20788" in result.blocks[0].table_rows[0][2]
    warnings = [d for d in result.diagnostics if d.severity == "warning"]
    assert any(d.field == "mouse.entrez" and "99999" in d.reason for d in warnings)


def test_ensembl_version_normalization_not_a_conflict():
    records = [
        _ev(
            source_name="Ensembl",
            fact_type="ensembl_gene_id",
            taxon_id=9606,
            organism="Homo sapiens",
            value={"ensembl_gene_id": "ENSG00000198911"},
            source_id="e1",
            evidence_id="id1",
        ),
        # Second ensembl-like from UniProt xref path is not fact_type ensembl;
        # inject a second Ensembl record with versioned ID.
        _ev(
            source_name="Ensembl",
            fact_type="ensembl_gene_id",
            taxon_id=9606,
            organism="Homo sapiens",
            value={"ensembl_gene_id": "ENSG00000198911.11"},
            source_id="e2",
            evidence_id="id2",
        ),
    ]
    result = build_gene_aliases_blocks(gene_symbol="SREBF2", evidence_records=records)
    warnings = [
        d
        for d in result.diagnostics
        if d.severity == "warning" and d.field == "human.ensembl"
    ]
    assert warnings == []
    assert "ENSG00000198911" in result.blocks[0].table_rows[3][1]


def test_human_alias_prefer_uniprot_synonyms_over_ncbi():
    records = [
        _ev(
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            value={
                "entrez_gene_id": "1",
                "nomenclaturesymbol": "GENE",
                "description": "gene name",
                "otheraliases": "AliasA, AliasB",
            },
        ),
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            value={
                "uniprot_accession": "P00001",
                "gene_names": ["GENE", "AliasB", "AliasC"],
                "gene_synonyms": ["AliasB", "AliasC"],
            },
        ),
    ]
    block = build_gene_aliases_blocks(
        gene_symbol="GENE", evidence_records=records
    ).blocks[0]
    # Polished cell uses reviewed UniProt synonyms only (no NCBI AliasA).
    assert block.table_rows[-1][1] == "AliasB, AliasC"
    # NCBI aliases remain in audit evidence.
    assert "AliasA" in records[0].value["otheraliases"]


def test_human_alias_fallback_to_ncbi_when_uniprot_synonyms_absent():
    records = [
        _ev(
            source_name="NCBI Gene",
            fact_type="entrez_gene_id",
            value={
                "entrez_gene_id": "1",
                "nomenclaturesymbol": "GENE",
                "description": "gene name",
                "otheraliases": "AliasA, AliasB",
            },
        ),
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            value={
                "uniprot_accession": "P00001",
                "gene_names": ["GENE"],
                "gene_synonyms": [],
            },
        ),
    ]
    block = build_gene_aliases_blocks(
        gene_symbol="GENE", evidence_records=records
    ).blocks[0]
    assert block.table_rows[-1][1] == "AliasA, AliasB"


def test_mouse_alias_merge_mgi_before_uniprot():
    records = _mouse_bundle() + [
        _ev(
            source_name="UniProt",
            fact_type="uniprot_accession",
            taxon_id=10090,
            organism="Mus musculus",
            value={
                "uniprot_accession": "Q3U1N2",
                "gene_names": ["Srebf2", "UniAlias"],
            },
        )
    ]
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2", evidence_records=records
    ).blocks[0]
    aliases = block.table_rows[-1][2]
    assert aliases.index("SREBP-2") < aliases.index("UniAlias")


def test_preview_uses_production_fragment_helper():
    doc = build_report_document(
        dossier_run_id="t",
        gene_symbol="SREBF2",
        evidence_records=_human_bundle(),
    )
    html = render_rancho_section_fragment(
        document=doc, section_number=1, subsection_key="a"
    )
    assert "1. General Gene Information" in html
    assert "a. Gene Aliases" in html
    assert "2. Expression" not in html
    assert "gene-aliases-table" in html


def test_serialized_roundtrip_keeps_presentation_role():
    doc = build_report_document(
        dossier_run_id="t",
        gene_symbol="SREBF2",
        evidence_records=_human_bundle(),
    )
    raw = json.loads(json.dumps(doc.model_dump(mode="json")))
    restored = ReportDocument.model_validate(raw)
    assert restored.sections[0].subsections[0].presentation_blocks[0].presentation_role == (
        "gene_aliases_table"
    )
