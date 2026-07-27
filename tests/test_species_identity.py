"""Offline fixtures and tests for multi-species gene identity normalization."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from gene_dossier.models import ToolResult
from gene_dossier.normalize.gene_identity import (
    ensembl_gene_id_from_xref,
    normalize_ensembl,
    normalize_ncbi_gene,
    normalize_uniprot,
    split_alias_field,
)
from gene_dossier.report_presentation import NOT_AVAILABLE, build_section_presentation
from gene_dossier.tools.ncbi_gene import select_safe_gene_match
from gene_dossier.tools.uniprot import summarize_entry

FIXTURES = Path(__file__).parent / "fixtures" / "species_identity"
RUN_ID = "species-identity-test-run"


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text())


def _identity_result(
    *,
    source_name: str,
    endpoint_name: str,
    payload: dict,
    request_url: str,
    query_symbol: str = "SREBF2",
    resolved_symbol: str | None = None,
    extra_params: dict | None = None,
) -> ToolResult:
    resolved = resolved_symbol or query_symbol
    params = {
        "query_symbol": query_symbol,
        "resolved_symbol": resolved,
        "species_identity": True,
    }
    if extra_params:
        params.update(extra_params)
    data = dict(payload)
    data.setdefault("query_gene_symbol", query_symbol)
    data.setdefault("species_gene_symbol", resolved)
    data.setdefault("query_symbol", query_symbol)
    data.setdefault("resolved_symbol", resolved)
    return ToolResult(
        source_name=source_name,
        endpoint_name=endpoint_name,
        success=True,
        gene_symbol=query_symbol,
        request_url=request_url,
        request_params=params,
        data=data,
    )


def _ncbi_result(
    payload: dict,
    *,
    query_symbol: str = "SREBF2",
    resolved_symbol: str | None = None,
) -> ToolResult:
    resolved = resolved_symbol or (
        (payload.get("selected_summary") or {}).get("nomenclaturesymbol")
        or query_symbol
    )
    return _identity_result(
        source_name="NCBI Gene",
        endpoint_name="lookup_gene",
        payload=payload,
        request_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene",
        query_symbol=query_symbol,
        resolved_symbol=resolved,
        extra_params={"organism": payload.get("organism")},
    )


def _ensembl_result(
    payload: dict,
    *,
    resolved_symbol: str,
    query_symbol: str = "SREBF2",
) -> ToolResult:
    return _identity_result(
        source_name="Ensembl",
        endpoint_name="lookup_symbol",
        payload=payload,
        request_url=(
            f"https://rest.ensembl.org/lookup/symbol/"
            f"{payload['species']}/{resolved_symbol}?content-type=application/json"
        ),
        query_symbol=query_symbol,
        resolved_symbol=resolved_symbol,
        extra_params={
            "content-type": "application/json",
            "species": payload["species"],
        },
    )


def _uniprot_result(
    payload: dict,
    *,
    resolved_symbol: str,
    query_symbol: str = "SREBF2",
) -> ToolResult:
    return _identity_result(
        source_name="UniProt",
        endpoint_name="search_reviewed",
        payload=payload,
        request_url="https://rest.uniprot.org/uniprotkb/search?query=reviewed:true",
        query_symbol=query_symbol,
        resolved_symbol=resolved_symbol,
        extra_params={"organism_id": payload.get("organism_id")},
    )


# ---------------------------------------------------------------------------
# NCBI
# ---------------------------------------------------------------------------
def test_ncbi_human_esummary_normalization():
    payload = _load("ncbi_human_lookup.json")
    original = copy.deepcopy(payload)
    records = normalize_ncbi_gene(_ncbi_result(payload), dossier_run_id=RUN_ID)
    assert len(records) == 1
    rec = records[0]
    assert rec.taxon_id == 9606
    assert rec.value["entrez_gene_id"] == "6721"
    assert rec.value["nomenclaturesymbol"] == "SREBF2"
    assert "sterol regulatory element binding transcription factor 2" in str(
        rec.value["gene_name"]
    ).lower()
    assert "BHLHD2" in rec.value["aliases"]
    assert "SREBP2" in rec.value["aliases"]
    assert payload == original


def test_ncbi_mouse_esummary_normalization():
    payload = _load("ncbi_mouse_lookup.json")
    records = normalize_ncbi_gene(
        _ncbi_result(payload, resolved_symbol="Srebf2"), dossier_run_id=RUN_ID
    )
    rec = records[0]
    assert rec.taxon_id == 10090
    assert rec.value["entrez_gene_id"] == "20788"
    assert rec.value["nomenclaturesymbol"] == "Srebf2"
    assert "SREBP2" in rec.value["aliases"]
    assert "bHLHd2" in rec.value["aliases"]


def test_ncbi_rat_candidate_selection_prefers_300095():
    summaries = _load("ncbi_rat_summaries.json")
    selected, warnings = select_safe_gene_match(
        summaries, "Srebf2", expected_taxid=10116
    )
    assert selected == "300095"
    assert warnings == []


def test_retired_rat_404651_redirects_and_loses_to_300095():
    summaries = _load("ncbi_rat_summaries.json")
    retired = summaries["404651"]
    assert str(retired.get("currentid") or "").strip() == "300095"
    selected, _ = select_safe_gene_match(summaries, "Srebf2", expected_taxid=10116)
    assert selected == "300095"
    assert selected != "404651"
    assert selected != "499505"


def test_ambiguous_active_records_produce_warning():
    summaries = _load("ncbi_ambiguous_summaries.json")
    selected, warnings = select_safe_gene_match(
        summaries, "FAKEG", expected_taxid=9606
    )
    assert selected is None
    assert any(w.startswith("ambiguous_safe_matches:") for w in warnings)


def test_mouse_aliases_split_and_normalize():
    raw = "SREBP-2, SREBP2, SREBP2gc, bHLHd2, lop13, nuc"
    aliases = split_alias_field(raw)
    assert aliases == ["SREBP-2", "SREBP2", "SREBP2gc", "bHLHd2", "lop13", "nuc"]


def test_rat_aliases_split_and_normalize():
    raw = "SREBP-2, SREBP2, Srebf2_retired"
    aliases = split_alias_field(raw)
    assert aliases == ["SREBP-2", "SREBP2", "Srebf2_retired"]


# ---------------------------------------------------------------------------
# Ensembl
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "fixture,symbol,tax,ens_id",
    [
        ("ensembl_human.json", "SREBF2", 9606, "ENSG00000198911"),
        ("ensembl_mouse.json", "Srebf2", 10090, "ENSMUSG00000022463"),
        ("ensembl_rat.json", "Srebf2", 10116, "ENSRNOG00000007400"),
    ],
)
def test_ensembl_species_lookup_normalization(fixture, symbol, tax, ens_id):
    payload = _load(fixture)
    original = copy.deepcopy(payload)
    records = normalize_ensembl(
        _ensembl_result(payload, resolved_symbol=symbol), dossier_run_id=RUN_ID
    )
    gene_recs = [r for r in records if r.fact_type == "ensembl_gene_id"]
    assert len(gene_recs) == 1
    rec = gene_recs[0]
    assert rec.taxon_id == tax
    assert rec.value["ensembl_gene_id"] == ens_id
    assert rec.value["tax_id"] == tax
    assert "api_key" not in (rec.value.get("source_url") or "")
    assert payload == original


# ---------------------------------------------------------------------------
# UniProt
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "fixture,symbol,tax,accession",
    [
        ("uniprot_human.json", "SREBF2", 9606, "Q12772"),
        ("uniprot_mouse.json", "Srebf2", 10090, "Q3U1N2"),
        ("uniprot_rat.json", "Srebf2", 10116, "Q3T1I5"),
    ],
)
def test_uniprot_species_reviewed_normalization(fixture, symbol, tax, accession):
    payload = _load(fixture)
    original = copy.deepcopy(payload)
    records = normalize_uniprot(
        _uniprot_result(payload, resolved_symbol=symbol), dossier_run_id=RUN_ID
    )
    acc = [r for r in records if r.fact_type == "uniprot_accession"]
    assert len(acc) == 1
    rec = acc[0]
    assert rec.taxon_id == tax
    assert rec.value["uniprot_accession"] == accession
    assert rec.value["reviewed"] is True
    assert rec.value["tax_id"] == tax
    assert payload == original


def test_species_taxon_ids_remain_attached_to_all_records():
    human = normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_human_lookup.json")), dossier_run_id=RUN_ID
    )
    mouse = normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_mouse_lookup.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    rat_ens = normalize_ensembl(
        _ensembl_result(_load("ensembl_rat.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    mouse_up = normalize_uniprot(
        _uniprot_result(_load("uniprot_mouse.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    assert all(r.taxon_id == 9606 for r in human)
    assert all(r.taxon_id == 10090 for r in mouse)
    assert all(r.taxon_id == 10116 for r in rat_ens)
    assert all(r.taxon_id == 10090 for r in mouse_up)


def test_direct_ensembl_outranks_uniprot_ensembl_xrefs():
    ens = normalize_ensembl(
        _ensembl_result(_load("ensembl_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    up = normalize_uniprot(
        _uniprot_result(_load("uniprot_human_conflicting_ensembl.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    records = ens + up
    result = build_section_presentation(
        section_key="1a",
        gene_symbol="SREBF2",
        evidence_records=records,
    )
    table = next(
        b for b in result.blocks if getattr(b, "presentation_role", None) == "gene_aliases_table"
    )
    ensembl_row = table.table_rows[3]
    human_cell = ensembl_row[1]
    assert "ENSG00000198911" in human_cell
    assert "ENSG99999999999" not in human_cell
    assert any("UniProt xref" in d.reason for d in result.diagnostics)


def test_different_uniprot_ensembl_xrefs_preserved_as_secondary():
    up = normalize_uniprot(
        _uniprot_result(_load("uniprot_human_conflicting_ensembl.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    xrefs = [r for r in up if r.fact_type == "ensembl_xref"]
    assert xrefs
    assert any(r.value.get("ensembl_gene_id") == "ENSG99999999999" for r in xrefs)
    assert all(r.value.get("primary") is False for r in xrefs)


def test_api_credentials_do_not_appear_in_source_urls():
    tr = ToolResult(
        source_name="NCBI Gene",
        endpoint_name="lookup_gene",
        success=True,
        gene_symbol="SREBF2",
        request_url=(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            "?db=gene&id=6721&api_key=SECRETKEY"
        ),
        data=_load("ncbi_human_lookup.json"),
    )
    original = copy.deepcopy(tr.data)
    records = normalize_ncbi_gene(tr, dossier_run_id=RUN_ID)
    url = records[0].value.get("source_url") or ""
    assert "SECRETKEY" not in url
    assert "api_key" not in url
    assert tr.data == original


def test_section_1a_fills_all_18_species_identity_cells():
    records = []
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_human_lookup.json")), dossier_run_id=RUN_ID
    )
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_mouse_lookup.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_rat_lookup.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_ensembl(
        _ensembl_result(_load("ensembl_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_ensembl(
        _ensembl_result(_load("ensembl_mouse.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_ensembl(
        _ensembl_result(_load("ensembl_rat.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_uniprot(
        _uniprot_result(_load("uniprot_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_uniprot(
        _uniprot_result(_load("uniprot_mouse.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_uniprot(
        _uniprot_result(_load("uniprot_rat.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    # Snapshot ids before presentation
    before = [(r.id, copy.deepcopy(r.value)) for r in records]
    result = build_section_presentation(
        section_key="1a", gene_symbol="SREBF2", evidence_records=records
    )
    table = next(
        b for b in result.blocks if getattr(b, "presentation_role", None) == "gene_aliases_table"
    )
    # 6 field rows × 3 species cells (skip label col)
    missing = []
    for row in table.table_rows:
        for cell in row[1:]:
            text = str(cell).strip()
            if not text or text == NOT_AVAILABLE:
                missing.append((row[0], text))
    assert not missing, f"missing cells: {missing}"
    after = [(r.id, r.value) for r in records]
    assert before == after


def test_input_payloads_and_evidence_not_mutated():
    payload = _load("ncbi_human_lookup.json")
    original_payload = copy.deepcopy(payload)
    tr = _ncbi_result(payload)
    records = normalize_ncbi_gene(tr, dossier_run_id=RUN_ID)
    original_records = [copy.deepcopy(r.model_dump(mode="json")) for r in records]
    build_section_presentation(
        section_key="1a", gene_symbol="SREBF2", evidence_records=records
    )
    assert payload == original_payload
    assert [r.model_dump(mode="json") for r in records] == original_records


def test_existing_human_identity_behavior_remains():
    records = normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_human_lookup.json")), dossier_run_id=RUN_ID
    )
    records += normalize_ensembl(
        _ensembl_result(_load("ensembl_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_uniprot(
        _uniprot_result(_load("uniprot_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    result = build_section_presentation(
        section_key="1a", gene_symbol="SREBF2", evidence_records=records
    )
    table = next(
        b for b in result.blocks if getattr(b, "presentation_role", None) == "gene_aliases_table"
    )
    human_cells = [row[1] for row in table.table_rows]
    joined = " | ".join(str(c) for c in human_cells)
    assert "6721" in joined
    assert "SREBF2" in joined
    assert "ENSG00000198911" in joined
    assert "Q12772" in joined


def test_uniprot_raw_geneid_property_extraction():
    raw = _load("uniprot_human_raw_entry.json")
    original = copy.deepcopy(raw)
    summary = summarize_entry(raw)
    assert summary["ensembl_gene_ids"] == ["ENSG00000198911"]
    assert summary["ensembl_xrefs"] == ["ENSG00000198911"]
    assert "ENST00000389809" in summary["ensembl_transcript_ids"]
    assert "ENST00000457199" in summary["ensembl_transcript_ids"]
    assert all(not t.startswith("ENSG") for t in summary["ensembl_transcript_ids"])
    assert raw == original


def test_transcript_ids_not_treated_as_gene_ids():
    assert ensembl_gene_id_from_xref("ENST00000389809.8") is None
    assert ensembl_gene_id_from_xref("ENSMUST00000023104.8") is None
    assert ensembl_gene_id_from_xref("ENSRNOT00000009912.6") is None
    assert ensembl_gene_id_from_xref("ENSG00000198911.14") == "ENSG00000198911"


def test_query_gene_symbol_retained_on_all_species_records():
    records = []
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_human_lookup.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_mouse_lookup.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_rat_lookup.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_ensembl(
        _ensembl_result(_load("ensembl_mouse.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_uniprot(
        _uniprot_result(_load("uniprot_rat.json"), resolved_symbol="Srebf2"),
        dossier_run_id=RUN_ID,
    )
    assert records
    assert all(r.gene_symbol == "SREBF2" for r in records)
    mouse = next(r for r in records if r.taxon_id == 10090 and r.source_name == "NCBI Gene")
    rat = next(r for r in records if r.taxon_id == 10116 and r.source_name == "NCBI Gene")
    human = next(r for r in records if r.taxon_id == 9606 and r.source_name == "NCBI Gene")
    assert human.official_symbol == "SREBF2"
    assert mouse.official_symbol == "Srebf2"
    assert rat.official_symbol == "Srebf2"
    assert mouse.value["species_gene_symbol"] == "Srebf2"
    assert mouse.value["query_gene_symbol"] == "SREBF2"
    assert rat.value["species_gene_symbol"] == "Srebf2"


def test_human_polished_aliases_use_reviewed_uniprot_synonyms():
    records = []
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_human_lookup.json")), dossier_run_id=RUN_ID
    )
    # NCBI has BHLHD2, SREBP2; also inject an extra NCBI-only alias into evidence
    records[0].value["aliases"] = list(records[0].value["aliases"]) + ["NCBI_ONLY_ALIAS"]
    records[0].value["otheraliases"] = list(records[0].value["aliases"])
    records += normalize_uniprot(
        _uniprot_result(_load("uniprot_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    result = build_section_presentation(
        section_key="1a", gene_symbol="SREBF2", evidence_records=records
    )
    table = next(b for b in result.blocks if b.presentation_role == "gene_aliases_table")
    human_aliases = table.table_rows[5][1]
    assert human_aliases == "BHLHD2, SREBP2"
    assert "NCBI_ONLY_ALIAS" not in human_aliases
    # Audit evidence still retains NCBI aliases
    ncbi = next(r for r in records if r.source_name == "NCBI Gene")
    assert "NCBI_ONLY_ALIAS" in ncbi.value["aliases"]


def test_section_preview_pdf_and_png_when_pymupdf_available(tmp_path: Path):
    fitz = pytest.importorskip("fitz")
    from gene_dossier.rancho_report import (
        rasterize_pdf_page_to_png,
        render_rancho_pdf,
        render_rancho_section_fragment,
    )
    from gene_dossier.report_schema import build_report_document

    records = []
    records += normalize_ncbi_gene(
        _ncbi_result(_load("ncbi_human_lookup.json")), dossier_run_id=RUN_ID
    )
    records += normalize_ensembl(
        _ensembl_result(_load("ensembl_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    records += normalize_uniprot(
        _uniprot_result(_load("uniprot_human.json"), resolved_symbol="SREBF2"),
        dossier_run_id=RUN_ID,
    )
    doc = build_report_document(
        dossier_run_id=RUN_ID,
        gene_symbol="SREBF2",
        evidence_records=records,
    )
    html = render_rancho_section_fragment(
        document=doc, section_number=1, subsection_key="a"
    )
    pdf_path = tmp_path / "SREBF2_1a_gene_aliases.pdf"
    png_path = tmp_path / "SREBF2_1a_gene_aliases.png"
    written_pdf = render_rancho_pdf(html, pdf_path, stamp_page_chrome=False)
    assert written_pdf is not None and written_pdf.is_file()
    written_png = rasterize_pdf_page_to_png(written_pdf, png_path, dpi=150)
    assert written_png is not None and written_png.is_file()
    with fitz.open(written_pdf) as pdf_doc:
        assert pdf_doc.page_count >= 1
