"""Offline tests for Section 1c official CDD/PDBe asset helpers."""

from __future__ import annotations

from typing import Any

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
    new_id,
)
from gene_dossier.section_1c import (
    cdd_domain_rows,
    cdd_master_lineage,
    choose_pdbe_official_image,
    parse_cdd_family_html,
    protein_seeds_by_species,
    rank_pdb_candidates,
    select_authoritative_protein_seed,
)


def _ev(**kwargs: Any) -> EvidenceRecord:
    defaults = dict(
        id=new_id(),
        source_id=f"src-{new_id()[:8]}",
        dossier_run_id="run-1c",
        gene_symbol="SREBF2",
        section="Known structure / domains",
        source_name="CDD",
        source_type=SourceType.structure_database,
        assertion_type=AssertionType.protein_structure,
        fact_type="conserved_domain_hit",
        evidence_grade=EvidenceGrade.E,
        taxon_id=9606,
        organism="Homo sapiens",
        value={},
        display_text="x",
    )
    defaults.update(kwargs)
    return EvidenceRecord(**defaults)


def test_protein_seed_uses_uniprot_length_not_observed_spans():
    seed = select_authoritative_protein_seed(
        gene_symbol="SREBF2",
        gene_ids={},
        evidence_records=[
            _ev(
                source_name="UniProt",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="uniprot_accession",
                evidence_grade=EvidenceGrade.C,
                value={
                    "uniprot_accession": "Q12772",
                    "reviewed": True,
                    "taxon_id": 9606,
                    "protein_length": 1141,
                    "refseq_protein_accessions": ["NP_004590.2"],
                },
            ),
            _ev(value={"from_residue": "326", "to_residue": "402"}),
        ],
    )
    assert seed.uniprot_accession == "Q12772"
    assert seed.refseq_protein == "NP_004590.2"
    assert seed.protein_length == 1141
    rows = cdd_domain_rows([_ev(value={"from_residue": "326", "to_residue": "402"})], protein_length=None)
    assert rows[0]["coverage"] is None


def test_cdd_master_lineage_tracks_master_and_extended_ids():
    lineage = cdd_master_lineage(
        {
            "master_cdsid": "QM3-qcdsearch-MASTER",
            "hits_request_cdsid": "QM3-qcdsearch-MASTER-HITS",
            "features_request_cdsid": "QM3-qcdsearch-MASTER-FEATS",
        }
    )
    assert lineage["master_cdsid"] == "QM3-qcdsearch-MASTER"
    assert lineage["hits_request_cdsid"] == "QM3-qcdsearch-MASTER-HITS"
    assert lineage["features_request_cdsid"] == "QM3-qcdsearch-MASTER-FEATS"
    assert lineage["graphical_result_master_cdsid"] == "QM3-qcdsearch-MASTER"
    assert lineage["same_master_job"] is True
    assert lineage["query_index"] == 0


def test_cdd_family_html_parser_is_versioned_and_preserves_thumbnails():
    raw = """
    <html><body>
      <h1>cd11304 Cadherin_repeat</h1>
      <p>PSSM-ID: 206637</p>
      <p>Summary: Cadherins are glycoproteins involved in Ca2+-mediated cell adhesion.</p>
      <img src="/Structure/cdd/cdThumbnail.cgi?uid=206637">
      <a href="/Structure/cdd/cdThumbnail.cgi?uid=206637&feature=1">feature</a>
    </body></html>
    """
    parsed = parse_cdd_family_html(
        raw,
        requested_uid="206637",
        source_page_url="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=206637",
    )
    assert parsed["parser_name"] == "ncbi_cdd_family_html"
    assert parsed["parser_version"] == "1"
    assert parsed["requested_uid"] == "206637"
    assert parsed["canonical_accession"].lower() == "cd11304"
    assert parsed["pssm_id"] == "206637"
    assert parsed["unparsed_sections_preserved"] is True
    assert len(parsed["thumbnail_urls"]) == 2


def test_cdd_family_parser_rejects_alignment_thumbnail_before_domain_context():
    raw = """
    <html><body>
      <h1>cd18922 bHLHzip_SREBP2</h1>
      <p>PSSM-ID: 381492</p>
      <div class="domain alignment">
        <p>Domain sequence alignment logo</p>
        <img src="/Structure/cdd/cdThumbnail.cgi?uid=381492&seqgraphic=1">
      </div>
    </body></html>
    """
    parsed = parse_cdd_family_html(
        raw,
        requested_uid="381492",
        source_page_url="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=381492",
    )
    candidate = parsed["thumbnail_candidates"][0]
    assert candidate["classified_role"] == "alignment_or_sequence_thumbnail"
    assert candidate["rejection_reason"] == "alignment or sequence thumbnail"


def test_pdb_ranking_requires_exact_selected_accession_mapping_and_null_coverage():
    payload = {
        "structure_summaries": [
            {"pdb_id": "2bbb", "chain_id": "A", "unp_start": 1, "unp_end": 50, "coverage": 0.5, "resolution": 2.0},
            {"pdb_id": "1aaa", "chain_id": "B", "unp_start": 10, "unp_end": 80, "coverage": 0.7, "resolution": 3.0},
            {"pdb_id": "3ccc", "chain_id": "C", "unp_start": 1, "unp_end": 200, "coverage": 0.9, "resolution": 1.5},
        ],
        "uniprot_mappings": {
            "2bbb": {"2bbb": {"UniProt": {"QX1111": {"mappings": [{"chain_id": "A", "unp_start": 1, "unp_end": 50}]}}}},
            "1aaa": {"1aaa": {"UniProt": {"QX1111": {"mappings": [{"chain_id": "B", "unp_start": 10, "unp_end": 80}]}}}},
            "3ccc": {"3ccc": {"UniProt": {"OTHER": {"mappings": [{"chain_id": "C", "unp_start": 1, "unp_end": 200}]}}}},
        },
    }
    ranked = rank_pdb_candidates(payload, selected_uniprot_accession="QX1111", protein_length=None)
    assert ranked[0].pdb_id == "1aaa"
    assert ranked[0].selected is True
    assert ranked[0].calculated_coverage is None
    rejected = {candidate.pdb_id: candidate for candidate in ranked if not candidate.selected}
    assert "missing exact mapping" in rejected["3ccc"].rejection_reasons[0]
    assert rejected["2bbb"].provider_coverage == 0.5


def test_species_seeds_are_tiered_human_mouse_rat():
    records = [
        _ev(source_name="UniProt", fact_type="uniprot_accession", value={"uniprot_accession": "P70408", "taxon_id": 10090, "reviewed": True}),
        _ev(source_name="UniProt", fact_type="uniprot_accession", value={"uniprot_accession": "Q9Y6N8", "taxon_id": 9606, "reviewed": True}),
        _ev(source_name="UniProt", fact_type="uniprot_accession", value={"uniprot_accession": "Q920Q8", "taxon_id": 10116, "reviewed": True}),
    ]
    seeds = protein_seeds_by_species(records)
    assert [seed.uniprot_accession for seed in seeds] == ["Q9Y6N8", "P70408", "Q920Q8"]


def test_pdbe_image_inventory_role_policy_prefers_assembly_for_full_structure():
    inventory = {
        "images": [
            {"basename": "1ukl_deposited_chain_front", "suffixes": ["png"]},
            {"basename": "1ukl_assembly-1_chain_front", "suffixes": [".png"]},
            {"basename": "1ukl_entity-2_front", "suffixes": ["png"]},
        ]
    }
    choice = choose_pdbe_official_image(inventory, pdb_id="1ukl", preferred_assembly_id="1")
    assert choice["selected_image_name"] == "1ukl_assembly-1_chain_front.png"
    assert choice["selection_reason"] == "preferred biological assembly image available"
    assert "1ukl_deposited_chain_front.png" in choice["available_image_candidates"]


def test_forbidden_rendering_helpers_are_not_public_api():
    import gene_dossier.section_1c as s1c

    assert not hasattr(s1c, "render_pymol_png")
    assert not hasattr(s1c, "render_mmcif_projection_png")
    assert not hasattr(s1c, "parse_sifts_residue_mappings")
    assert not hasattr(s1c, "build_domain_architecture_svg")
