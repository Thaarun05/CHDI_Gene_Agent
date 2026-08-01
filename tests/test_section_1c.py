"""Offline tests for Section 1c official CDD/PDBe asset helpers."""

from __future__ import annotations

from typing import Any

from gene_dossier.config import Settings
from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    RawArtifact,
    SourceType,
    new_id,
)
from gene_dossier.report_presentation import build_known_structure_blocks
from gene_dossier.section_1c import (
    cdd_domain_rows,
    cdd_master_lineage,
    choose_pdbe_official_image,
    node_generate_section_1c_derived_artifacts,
    parse_cdd_family_html,
    protein_seeds_by_species,
    rank_pdb_candidates,
    select_authoritative_protein_seed,
)
import gene_dossier.section_1c as s1c


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


def _settings(tmp_path) -> Settings:
    return Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")


def _base_state(*, gene: str = "SREBF2", length: int | None = 1141) -> dict[str, Any]:
    evidence = []
    if length is not None:
        evidence.append(
            _ev(
                gene_symbol=gene,
                source_name="UniProt",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="uniprot_accession",
                evidence_grade=EvidenceGrade.C,
                value={
                    "uniprot_accession": "Q12772" if gene == "SREBF2" else "Q9Y6N8",
                    "reviewed": True,
                    "taxon_id": 9606,
                    "protein_length": length,
                    "refseq_protein_accessions": ["NP_004590.2"],
                },
            )
        )
    return {
        "dossier_run_id": "run-arch-1c",
        "gene_symbol": gene,
        "run_type": "section_bundle",
        "selected_section_keys": ["1c"],
        "gene_ids": {},
        "evidence_records": evidence,
        "api_runs": [],
        "raw_artifacts": [],
        "errors": [],
        "tool_results": [],
    }


def _fake_official_success(*_args: Any, **_kwargs: Any):
    artifact = RawArtifact(
        dossier_run_id="run-arch-1c",
        source_name="CDD",
        artifact_type="image",
        file_path="/tmp/official.png",
        content_hash="a" * 64,
    )
    value = {
        "status": "success",
        "origin": "official",
        "relative_path": "tests/official-architecture.png",
        "media_type": "image/png",
        "width": 800,
        "height": 120,
        "sha256": artifact.content_hash,
        "artifact_class": "official",
        "artifact_origin": "ncbi_cdsearch_official_graphic",
        "figure_raw_artifact_id": artifact.id,
    }
    return [], [], artifact, {"relative_path": value["relative_path"]}, value


def _fake_official_failure(*_args: Any, **_kwargs: Any):
    return [], [], None, None, {"status": "unavailable", "reason": "official capture failed"}


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
    assert not hasattr(s1c, "render_pymol_png")
    assert not hasattr(s1c, "render_mmcif_projection_png")
    assert not hasattr(s1c, "parse_sifts_residue_mappings")
    assert not hasattr(s1c, "build_domain_architecture_svg")
    assert hasattr(s1c, "_render_cdd_architecture_fallback")
    assert "_render_cdd_architecture_fallback" not in getattr(s1c, "__all__", [])


def test_architecture_policy_official_succeeds(monkeypatch, tmp_path):
    monkeypatch.setattr(s1c, "_discover_cdd_architecture", _fake_official_success)
    state = _base_state()
    state["evidence_records"].append(
        _ev(
            value={
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "from_residue": 343,
                "to_residue": 403,
                "superfamily": "cl00081",
            }
        )
    )
    state["tool_results"] = [
        type(
            "TR",
            (),
            {
                "source_name": "CDD",
                "success": True,
                "data": {"master_cdsid": "QM3-master", "hit_summaries": [{}]},
            },
        )()
    ]
    out = node_generate_section_1c_derived_artifacts(state, settings=_settings(tmp_path), persist_db=False)
    arch = out["section_1c"]["cdd_architecture"]
    assert out["section_1c"]["rendering_status"]["cdd_architecture"] == "success"
    assert arch["origin"] == "official"
    assert arch["selected_fact_type"] == "cdd_official_architecture_figure"
    facts = [r.fact_type for r in out["evidence_records"]]
    assert "cdd_official_architecture_figure" in facts
    assert "cdd_architecture_figure" not in facts


def test_architecture_policy_derived_when_official_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(s1c, "_discover_cdd_architecture", _fake_official_failure)
    state = _base_state()
    state["evidence_records"].append(
        _ev(
            value={
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "from_residue": 343,
                "to_residue": 403,
                "superfamily": "cl00081",
            }
        )
    )
    state["tool_results"] = [
        type(
            "TR",
            (),
            {
                "source_name": "CDD",
                "success": True,
                "data": {"master_cdsid": "QM3-master", "hit_summaries": [{}]},
            },
        )()
    ]
    out = node_generate_section_1c_derived_artifacts(state, settings=_settings(tmp_path), persist_db=False)
    arch = out["section_1c"]["cdd_architecture"]
    assert out["section_1c"]["rendering_status"]["cdd_architecture"] == "success"
    assert arch["origin"] == "derived"
    assert arch["selected_fact_type"] == "cdd_architecture_figure"
    assert arch["derivation_type"] == "cdd_domain_architecture_render"
    assert arch["artifact_class"] == "derived"
    assert arch.get("artifact_origin") != "ncbi_cdsearch_official_graphic"
    derived = next(r for r in out["evidence_records"] if r.fact_type == "cdd_architecture_figure")
    assert derived.value["artifact_class"] == "derived"


def test_architecture_policy_protein_length_unavailable(monkeypatch, tmp_path):
    monkeypatch.setattr(s1c, "_discover_cdd_architecture", _fake_official_failure)
    state = _base_state(length=None)
    state["evidence_records"].append(
        _ev(value={"domain_accession": "cd18922", "from_residue": 343, "to_residue": 403})
    )
    state["tool_results"] = [
        type(
            "TR",
            (),
            {
                "source_name": "CDD",
                "success": True,
                "data": {"master_cdsid": "QM3-master", "hit_summaries": [{}]},
            },
        )()
    ]
    out = node_generate_section_1c_derived_artifacts(state, settings=_settings(tmp_path), persist_db=False)
    arch = out["section_1c"]["cdd_architecture"]
    assert out["section_1c"]["rendering_status"]["cdd_architecture"] == "unavailable"
    assert arch["reason"] == "protein_length_unavailable"
    assert "protein_length_unavailable" in arch["reasons"]
    assert "official capture failed" in arch["reasons"]
    assert not any(r.fact_type.endswith("architecture_figure") for r in out["evidence_records"])


def test_architecture_policy_both_fail_preserves_reasons(monkeypatch, tmp_path):
    monkeypatch.setattr(s1c, "_discover_cdd_architecture", _fake_official_failure)

    def _boom(**_kwargs: Any):
        raise RuntimeError("derived renderer exploded")

    monkeypatch.setattr(s1c, "_render_cdd_architecture_fallback", _boom)
    state = _base_state()
    state["evidence_records"].append(
        _ev(value={"domain_accession": "cd18922", "from_residue": 343, "to_residue": 403})
    )
    state["tool_results"] = [
        type(
            "TR",
            (),
            {
                "source_name": "CDD",
                "success": True,
                "data": {"master_cdsid": "QM3-master", "hit_summaries": [{}]},
            },
        )()
    ]
    out = node_generate_section_1c_derived_artifacts(state, settings=_settings(tmp_path), persist_db=False)
    arch = out["section_1c"]["cdd_architecture"]
    assert arch["status"] == "unavailable"
    assert arch["origin"] == "unavailable"
    assert any("official capture failed" in r for r in arch["reasons"])
    assert any("derived renderer exploded" in r for r in arch["reasons"])


def test_presentation_includes_architecture_after_cdd_link_and_prefers_official():
    from io import BytesIO

    from PIL import Image

    from gene_dossier.config import get_settings

    root = get_settings().raw_data_path
    path = root / "tests" / "section-1c-arch-policy.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = BytesIO()
    Image.new("RGB", (64, 32), (20, 40, 80)).save(buf, format="PNG")
    path.write_bytes(buf.getvalue())
    rel = str(path.relative_to(root))

    domains = [
        _ev(
            value={
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "from_residue": 343,
                "to_residue": 403,
            }
        )
    ]
    official = _ev(
        fact_type="cdd_official_architecture_figure",
        value={"relative_path": rel, "width": 64, "height": 32, "artifact_class": "official"},
    )
    derived = _ev(
        fact_type="cdd_architecture_figure",
        value={
            "relative_path": rel,
            "width": 64,
            "height": 32,
            "artifact_class": "derived",
            "derivation_type": "cdd_domain_architecture_render",
        },
    )
    both = build_known_structure_blocks(
        gene_symbol="SREBF2",
        evidence_records=[*domains, official, derived],
    )
    roles = [b.presentation_role for b in both.blocks]
    assert "section_1c_cdd_link" in roles
    assert "section_1c_domain_architecture_figure" in roles
    assert roles.index("section_1c_cdd_link") < roles.index("section_1c_domain_architecture_figure")
    arch_block = next(
        b for b in both.blocks if b.presentation_role == "section_1c_domain_architecture_figure"
    )
    assert arch_block.evidence_record_ids == [official.id]

    derived_only = build_known_structure_blocks(
        gene_symbol="SREBF2",
        evidence_records=[*domains, derived],
    )
    roles2 = [b.presentation_role for b in derived_only.blocks]
    assert "section_1c_domain_architecture_figure" in roles2
    assert roles2.index("section_1c_cdd_link") < roles2.index("section_1c_domain_architecture_figure")


def test_hits_exist_architecture_not_silently_omitted_in_presentation():
    from io import BytesIO

    from PIL import Image

    from gene_dossier.config import get_settings

    root = get_settings().raw_data_path
    path = root / "tests" / "section-1c-arch-present.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = BytesIO()
    Image.new("RGB", (80, 40), (10, 100, 40)).save(buf, format="PNG")
    path.write_bytes(buf.getvalue())
    rel = str(path.relative_to(root))
    result = build_known_structure_blocks(
        gene_symbol="CDH10",
        evidence_records=[
            _ev(
                gene_symbol="CDH10",
                value={
                    "domain_accession": "cd11304",
                    "domain_short_name": "Cadherin_repeat",
                    "from_residue": 170,
                    "to_residue": 280,
                },
            ),
            _ev(
                gene_symbol="CDH10",
                fact_type="cdd_architecture_figure",
                value={
                    "relative_path": rel,
                    "width": 80,
                    "height": 40,
                    "artifact_class": "derived",
                    "derivation_type": "cdd_domain_architecture_render",
                },
            ),
        ],
    )
    assert any(b.presentation_role == "section_1c_domain_architecture_figure" for b in result.blocks)


def test_srebf2_shaped_derived_architecture_has_specific_and_superfamily_lanes(tmp_path):
    artifact, _meta, value = s1c._render_cdd_architecture_fallback(
        dossier_run_id="run-srebf2-arch",
        gene_symbol="SREBF2",
        protein_length=1141,
        domains=[
            {
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "from_residue": 343,
                "to_residue": 403,
                "superfamily": "cl00081",
            }
        ],
        features=[],
        settings=_settings(tmp_path),
        persist_db=False,
    )
    assert artifact is not None
    assert value["status"] == "success"
    assert value["specific_hit_count"] == 1
    assert value["superfamily_hit_count"] == 1
    svg, _w, _h = s1c._build_cdd_architecture_svg_bytes(
        protein_length=1141,
        gene_symbol="SREBF2",
        specific_hits=[
            {
                "domain_accession": "cd18922",
                "domain_short_name": "bHLHzip_SREBP2",
                "from_residue": 343,
                "to_residue": 403,
            }
        ],
        superfamily_hits=[
            {
                "domain_accession": "cl00081",
                "domain_short_name": "cl00081",
                "from_residue": 343,
                "to_residue": 403,
            }
        ],
        feature_markers=[],
    )
    text = svg.decode("utf-8")
    assert "Specific hits" in text
    assert "Superfamilies" in text
    assert "bHLHzip_SREBP2" in text


def test_cdh10_shaped_derived_architecture_keeps_repeat_instances_and_features(tmp_path):
    domains = [
        {
            "domain_accession": "cd11304",
            "domain_short_name": "Cadherin_repeat",
            "from_residue": 80,
            "to_residue": 160,
        },
        {
            "domain_accession": "cd11304",
            "domain_short_name": "Cadherin_repeat",
            "from_residue": 170,
            "to_residue": 260,
        },
        {
            "domain_accession": "cd11304",
            "domain_short_name": "Cadherin_repeat",
            "from_residue": 270,
            "to_residue": 360,
        },
        {
            "domain_accession": "pfam01049",
            "domain_short_name": "Cadherin_C",
            "from_residue": 720,
            "to_residue": 780,
        },
    ]
    features = [
        {
            "domain_accession": "cd11304",
            "feature_label": "Ca2+ binding site [ion binding site]",
            "query_residues": "E171, M172",
        }
    ]
    artifact, _meta, value = s1c._render_cdd_architecture_fallback(
        dossier_run_id="run-cdh10-arch",
        gene_symbol="CDH10",
        protein_length=788,
        domains=domains,
        features=features,
        settings=_settings(tmp_path),
        persist_db=False,
    )
    assert artifact is not None
    assert value["specific_hit_count"] == 4
    assert value["feature_marker_count"] >= 2
    specific, _sf, markers = s1c._architecture_lane_rows(domains, features)
    assert len(specific) == 4
    assert {m["position"] for m in markers} >= {171, 172}
