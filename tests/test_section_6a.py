"""Offline tests for Section 6a CTD aggregation and safeguards."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest

from gene_dossier.models import AssertionType
from gene_dossier.report_schema import resolve_report_slot
from gene_dossier.section_6a import (
    PARSER_VERSION,
    aggregate_by_chemical_id,
    filter_target_rows,
    pubmed_occurrence_count,
    rank_top_chemicals,
    render_top_chemicals_png,
    validate_top_chemicals_png,
    write_ctd_workbook,
)
from gene_dossier.section_bundle import (
    DEFAULT_SECTION_BUNDLE_KEYS,
    SECTION_SOURCE_DEPENDENCIES,
    SUPPORTED_SECTION_BUNDLE_KEYS,
)
from gene_dossier.tools import ctd as ctd_client

FIXTURE = Path(__file__).parent / "fixtures" / "ctd" / "CTD_chem_gene_ixns_mini.tsv.gz"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_supported_not_default_and_owns_http():
    """Section 6a freeze contract: supported/opt-in only; defaults unchanged."""
    assert "6a" in SUPPORTED_SECTION_BUNDLE_KEYS
    assert SUPPORTED_SECTION_BUNDLE_KEYS[-1] == "6a"
    assert "6a" not in DEFAULT_SECTION_BUNDLE_KEYS
    assert SECTION_SOURCE_DEPENDENCIES["6a"] == set()
    # Frozen: do not promote 6a into the default bundle.
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


def test_pubmed_occurrence_count():
    assert pubmed_occurrence_count("") == 0
    assert pubmed_occurrence_count(None) == 0
    assert pubmed_occurrence_count("111") == 1
    assert pubmed_occurrence_count("111|222") == 2
    assert pubmed_occurrence_count("111||222|") == 2


def test_filter_malformed_source_rows_do_not_fail_source():
    content = FIXTURE.read_bytes()
    result = filter_target_rows(
        content, target_gene_id=6721, target_symbol="SREBF2"
    )
    assert result["ok"] is True
    assert result["malformed_source_row_count"] >= 1
    assert result["filtered_target_row_count"] >= 3
    assert not result["target_row_validation_errors"]


def test_target_symbol_mismatch_fails_closed_for_gene():
    content = FIXTURE.read_bytes()
    result = filter_target_rows(
        content, target_gene_id=1008, target_symbol="CDH10"
    )
    assert result["ok"] is False
    assert result["error_type"] == "target_identity_mismatch"
    assert result["target_row_validation_errors"]
    assert any(
        e.get("reason") == "gene_symbol_mismatch"
        for e in result["target_row_validation_errors"]
    )


def test_cdh10_valproic_acid_interaction_count_nine_without_mismatch_row():
    """Construct CDH10 rows matching the approved fixture arithmetic."""
    rows = [
        {
            "ChemicalName": "Valproic Acid",
            "ChemicalID": "C004307",
            "CasRN": "99-66-1",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "mRNA",
            "Organism": "Homo sapiens",
            "OrganismID": "9606",
            "Interaction": "i1",
            "InteractionActions": "increases",
            "PubMedIDs": "111",
            "source_order": 0,
            "_pubmed_occurrence_count": 1,
        },
        {
            "ChemicalName": "Valproic Acid",
            "ChemicalID": "C004307",
            "CasRN": "99-66-1",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "protein",
            "Organism": "Homo sapiens",
            "OrganismID": "9606",
            "Interaction": "i2",
            "InteractionActions": "decreases",
            "PubMedIDs": "222",
            "source_order": 1,
            "_pubmed_occurrence_count": 1,
        },
        {
            "ChemicalName": "Valproic Acid",
            "ChemicalID": "C004307",
            "CasRN": "99-66-1",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "mRNA",
            "Organism": "Homo sapiens",
            "OrganismID": "9606",
            "Interaction": "i3",
            "InteractionActions": "affects",
            "PubMedIDs": "333|444",
            "source_order": 2,
            "_pubmed_occurrence_count": 2,
        },
        {
            "ChemicalName": "Valproic Acid",
            "ChemicalID": "C004307",
            "CasRN": "99-66-1",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "protein",
            "Organism": "Mus musculus",
            "OrganismID": "10090",
            "Interaction": "i4",
            "InteractionActions": "affects",
            "PubMedIDs": "555",
            "source_order": 3,
            "_pubmed_occurrence_count": 1,
        },
        {
            "ChemicalName": "Valproic Acid",
            "ChemicalID": "C004307",
            "CasRN": "99-66-1",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "mRNA",
            "Organism": "Homo sapiens",
            "OrganismID": "9606",
            "Interaction": "i5",
            "InteractionActions": "affects",
            "PubMedIDs": "666|777|888|999",
            "source_order": 4,
            "_pubmed_occurrence_count": 4,
        },
        {
            "ChemicalName": "Benzo(a)pyrene",
            "ChemicalID": "C006380",
            "CasRN": "50-32-8",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "mRNA",
            "Organism": "Homo sapiens",
            "OrganismID": "9606",
            "Interaction": "b1",
            "InteractionActions": "affects",
            "PubMedIDs": "10|20",
            "source_order": 5,
            "_pubmed_occurrence_count": 2,
        },
        {
            "ChemicalName": "Benzo(a)pyrene",
            "ChemicalID": "C006380",
            "CasRN": "50-32-8",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "protein",
            "Organism": "Homo sapiens",
            "OrganismID": "9606",
            "Interaction": "b2",
            "InteractionActions": "affects",
            "PubMedIDs": "30|40",
            "source_order": 6,
            "_pubmed_occurrence_count": 2,
        },
        {
            "ChemicalName": (
                "4-(5-benzo(1,3)dioxol-5-yl-4-pyridin-2-yl-1H-imidazol-2-yl)benzamide"
            ),
            "ChemicalID": "C559067",
            "CasRN": "",
            "GeneSymbol": "CDH10",
            "GeneID": "1008",
            "GeneForms": "protein",
            "Organism": "Homo sapiens",
            "OrganismID": "9606",
            "Interaction": "c1",
            "InteractionActions": "affects",
            "PubMedIDs": "1|2|3|4|5",
            "source_order": 7,
            "_pubmed_occurrence_count": 5,
        },
    ]
    chemicals = aggregate_by_chemical_id(rows)
    by_id = {c["chemical_id"]: c for c in chemicals}
    assert by_id["C004307"]["interaction_count"] == 9
    assert by_id["C006380"]["interaction_count"] == 4
    assert by_id["C559067"]["interaction_count"] == 5
    assert "Mus musculus" in by_id["C004307"]["organisms"]
    top = rank_top_chemicals(chemicals, limit=10)
    assert top[0]["chemical_id"] == "C004307"
    assert top[0]["display_rank"] == 1
    assert top[1]["chemical_id"] == "C559067"
    assert top[2]["chemical_id"] == "C006380"


def test_srebf2_fixture_top_counts():
    result = filter_target_rows(
        FIXTURE.read_bytes(), target_gene_id=6721, target_symbol="SREBF2"
    )
    assert result["ok"]
    chemicals = aggregate_by_chemical_id(result["filtered_rows"])
    by_id = {c["chemical_id"]: c for c in chemicals}
    assert by_id["D004041"]["interaction_count"] == 25  # Dietary Fats
    assert by_id["C006780"]["interaction_count"] == 24  # BPA
    assert by_id["C532164"]["interaction_count"] == 22  # fatostatin
    assert by_id["C999999"]["interaction_count"] == 0  # empty PMID
    top = rank_top_chemicals(chemicals, limit=10)
    assert [c["chemical_id"] for c in top[:3]] == ["D004041", "C006780", "C532164"]


def test_chart_rank1_at_top_without_mutating_ranked(tmp_path: Path):
    chemicals = [
        {"chemical_id": "C2", "chemical_name": "Beta", "interaction_count": 5},
        {"chemical_id": "C1", "chemical_name": "Alpha", "interaction_count": 10},
    ]
    ranked = rank_top_chemicals(chemicals, limit=10)
    assert ranked[0]["chemical_id"] == "C1"
    ranked_before = [dict(r) for r in ranked]
    png, meta = render_top_chemicals_png(ranked)
    assert meta["rank_1_at_top"] is True
    assert meta.get("title") in (None, "")
    assert ranked == ranked_before
    check = validate_top_chemicals_png(png)
    assert check["ok"] is True
    assert meta.get("bar_color") == "#F66400"
    (tmp_path / "fig.png").write_bytes(png)


def test_workbook_and_slot(tmp_path: Path):
    result = filter_target_rows(
        FIXTURE.read_bytes(), target_gene_id=6721, target_symbol="SREBF2"
    )
    chemicals = aggregate_by_chemical_id(result["filtered_rows"])
    top = rank_top_chemicals(chemicals, limit=10)
    path = tmp_path / "SREBF2_CTD.xlsx"
    write_ctd_workbook(
        path,
        gene_symbol="SREBF2",
        gene_id=6721,
        filtered_rows=result["filtered_rows"],
        chemicals=chemicals,
        top_chemicals=top,
        provenance=[("Parser version", PARSER_VERSION)],
    )
    assert path.is_file()
    digest = _sha(path)
    assert len(digest) == 64

    from gene_dossier.models import EvidenceRecord, EvidenceGrade, SourceType
    from gene_dossier.source_ids import make_source_id

    rec = EvidenceRecord(
        source_id=make_source_id(
            "CTD", "SREBF2", AssertionType.chemical_interaction, "x"
        ),
        dossier_run_id="run",
        gene_symbol="SREBF2",
        section="CTD perturbations",
        subsection="Comparative Toxicogenomics Database",
        source_name="CTD",
        source_type=SourceType.chemical_database,
        assertion_type=AssertionType.chemical_interaction,
        fact_type="section_6a_summary",
        evidence_grade=EvidenceGrade.C,
        value={},
        display_text="x",
    )
    slot = resolve_report_slot(rec)
    assert slot is not None
    assert slot.major_key == "6"
    assert slot.subsection_key == "a"


def test_no_batch_query_or_playwright_in_section_6a_module():
    src = Path("src/gene_dossier/section_6a.py").read_text(encoding="utf-8")
    assert "batchQuery" not in src
    assert "basicQuery" not in src
    assert "playwright" not in src.lower()
    assert "download_chem_gene_ixns_bulk" in src or "ctd_client" in src


def test_png_rejects_uniform_blank_image():
    from PIL import Image

    buf = __import__("io").BytesIO()
    Image.new("RGB", (500, 300), color=(255, 255, 255)).save(buf, format="PNG")
    check = validate_top_chemicals_png(buf.getvalue())
    assert check["ok"] is False
    assert check["reason"] == "blank_or_uniform"


def test_rank_top_chemicals_tie_break_detects_wrong_name_order():
    chemicals = [
        {"chemical_id": "C2", "chemical_name": "Zeta", "interaction_count": 9},
        {"chemical_id": "C1", "chemical_name": "Alpha", "interaction_count": 9},
        {"chemical_id": "C3", "chemical_name": "Beta", "interaction_count": 5},
    ]
    ranked = rank_top_chemicals(chemicals, limit=10)
    assert [c["chemical_id"] for c in ranked] == ["C1", "C2", "C3"]
    assert ranked[0]["chemical_name"] == "Alpha"


def test_accepted_pointer_reuse_skips_http(tmp_path, monkeypatch):
    """Pre-seeded accepted CTD pointer is reused with zero HTTP and same provenance."""
    from gene_dossier.config import Settings
    from gene_dossier.section_6a import PARSER_VERSION, _resolve_ctd_bulk_source
    from gene_dossier.section_6a_sources import (
        OFFICIAL_URL,
        SOURCE_KEY,
        accept_source,
        load_accepted_source,
        paths_for,
    )

    content = FIXTURE.read_bytes()
    calls = {"n": 0}

    def fake_download(*, settings=None, url=None):
        calls["n"] += 1
        raise AssertionError("CTD HTTP must not run when accepted pointer exists")

    monkeypatch.setattr(ctd_client, "download_chem_gene_ixns_bulk", fake_download)
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")
    paths = paths_for(tmp_path / "out")

    attempt = paths.new_source_attempt(SOURCE_KEY)
    artifact = attempt / "CTD_chem_gene_ixns.tsv.gz"
    artifact.write_bytes(content)
    digest = hashlib.sha256(content).hexdigest()
    accept_source(
        paths,
        source_key=SOURCE_KEY,
        attempt_dir=attempt,
        artifact_path=artifact,
        official_url=OFFICIAL_URL,
        sha256=digest,
        byte_size=len(content),
        validation={"required_columns_ok": True},
        extra={
            "api_run_id": "api-seed-reuse-1",
            "raw_artifact_id": "raw-seed-reuse-1",
            "source_attempt_id": attempt.name,
            "ctd_report_created": "Thu Jul 30 14:03:48 EDT 2026",
            "retrieval_timestamp": "2026-08-09T13:01:15+00:00",
            "parser_version": PARSER_VERSION,
            "dossier_run_id": "seed-reuse-run",
        },
    )
    seeded = load_accepted_source(paths)
    assert seeded is not None
    assert seeded["api_run_id"] == "api-seed-reuse-1"
    assert seeded["raw_artifact_id"] == "raw-seed-reuse-1"

    api_runs: list = []
    raws: list = []
    tools: list = []
    reused = _resolve_ctd_bulk_source(
        paths=paths,
        force_refresh=False,
        promote_source=False,
        dossier_run_id="gene-run-2",
        gene_symbol="SREBF2",
        settings=settings,
        persist_db=False,
        api_runs=api_runs,
        raw_artifacts=raws,
        tool_results=tools,
    )
    assert reused.ok is True
    assert reused.origin == "accepted_pointer"
    assert reused.api_run_id == "api-seed-reuse-1"
    assert reused.raw_artifact_id == "raw-seed-reuse-1"
    assert reused.source_attempt_id == attempt.name
    assert calls["n"] == 0
    assert tools == []


def test_no_db_acquisition_does_not_create_or_replace_accepted_pointer(
    tmp_path, monkeypatch
):
    """--no-db / persist_db=False must never pin accepted/sources/ctd_chem_gene_ixns.json."""
    import importlib.util

    from gene_dossier.config import Settings
    from gene_dossier.models import ToolResult
    from gene_dossier.section_6a import PARSER_VERSION, download_persist_and_pin_ctd_bulk
    from gene_dossier.section_6a_sources import (
        OFFICIAL_URL,
        SOURCE_KEY,
        accept_source,
        load_accepted_source,
        paths_for,
    )

    content = FIXTURE.read_bytes()

    def fake_download(*, settings=None, url=None):
        return ToolResult(
            source_name="CTD",
            endpoint_name="download_chem_gene_ixns_bulk",
            success=True,
            gene_symbol="",
            request_url=ctd_client.CHEM_GENE_IXNS_BULK_URL,
            request_params={},
            status_code=200,
            data={
                "content": content,
                "sha256": hashlib.sha256(content).hexdigest(),
                "byte_size": len(content),
                "final_url": ctd_client.CHEM_GENE_IXNS_BULK_URL,
                "content_type": "application/gzip",
            },
        )

    monkeypatch.setattr(ctd_client, "download_chem_gene_ixns_bulk", fake_download)
    out_root = tmp_path / "out"
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=out_root)
    paths = paths_for(out_root)
    pointer_path = paths.accepted_source_pointer(SOURCE_KEY)
    assert not pointer_path.exists()

    # First-miss --no-db: local attempt ok, no accepted pointer.
    first = download_persist_and_pin_ctd_bulk(
        paths=paths,
        dossier_run_id="nodb-run-1",
        settings=settings,
        gene_symbol="CTD_BULK",
        persist_db=False,
        promote=False,
        origin="acquire_script",
    )
    assert first.ok is True
    assert first.attempt_dir
    assert Path(first.attempt_dir).exists()
    assert first.pointer is None
    assert not pointer_path.exists()
    assert load_accepted_source(paths) is None

    # Seed an independent pre-existing accepted fixture (not via no-DB promotion).
    seed_attempt = paths.new_source_attempt(SOURCE_KEY)
    seed_artifact = seed_attempt / "CTD_chem_gene_ixns.tsv.gz"
    seed_artifact.write_bytes(content)
    digest = hashlib.sha256(content).hexdigest()
    accept_source(
        paths,
        source_key=SOURCE_KEY,
        attempt_dir=seed_attempt,
        artifact_path=seed_artifact,
        official_url=OFFICIAL_URL,
        sha256=digest,
        byte_size=len(content),
        validation={"required_columns_ok": True},
        extra={
            "api_run_id": "api-seed-existing-1",
            "raw_artifact_id": "raw-seed-existing-1",
            "source_attempt_id": seed_attempt.name,
            "ctd_report_created": "Thu Jul 30 14:03:48 EDT 2026",
            "retrieval_timestamp": "2026-08-09T13:01:15+00:00",
            "parser_version": PARSER_VERSION,
            "dossier_run_id": "seed-existing-run",
        },
    )
    before = pointer_path.read_text(encoding="utf-8")
    before_meta = load_accepted_source(paths)
    assert before_meta is not None
    assert before_meta["api_run_id"] == "api-seed-existing-1"
    assert before_meta["raw_artifact_id"] == "raw-seed-existing-1"

    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "acquire_ctd_chem_gene_ixns.py"
    )
    spec = importlib.util.spec_from_file_location(
        "acquire_ctd_chem_gene_ixns", script_path
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "get_settings", lambda: settings)
    rc = mod.main(
        [
            "--output-root",
            str(out_root),
            "--force-refresh",
            "--no-db",
        ]
    )
    assert rc == 0
    assert pointer_path.read_text(encoding="utf-8") == before
    after_meta = load_accepted_source(paths)
    assert after_meta is not None
    assert after_meta["api_run_id"] == before_meta["api_run_id"]
    assert after_meta["raw_artifact_id"] == before_meta["raw_artifact_id"]
    assert after_meta["sha256"] == before_meta["sha256"]


def test_helper_promote_true_persist_db_false_does_not_pin(tmp_path, monkeypatch, caplog):
    """Helper-level invariant: promote=True cannot pin when persist_db=False."""
    import logging

    from gene_dossier.config import Settings
    from gene_dossier.models import ToolResult
    from gene_dossier.section_6a import download_persist_and_pin_ctd_bulk
    from gene_dossier.section_6a_sources import (
        SOURCE_KEY,
        load_accepted_source,
        paths_for,
    )

    content = FIXTURE.read_bytes()

    def fake_download(*, settings=None, url=None):
        return ToolResult(
            source_name="CTD",
            endpoint_name="download_chem_gene_ixns_bulk",
            success=True,
            gene_symbol="",
            request_url=ctd_client.CHEM_GENE_IXNS_BULK_URL,
            request_params={},
            status_code=200,
            data={
                "content": content,
                "sha256": hashlib.sha256(content).hexdigest(),
                "byte_size": len(content),
                "final_url": ctd_client.CHEM_GENE_IXNS_BULK_URL,
                "content_type": "application/gzip",
            },
        )

    monkeypatch.setattr(ctd_client, "download_chem_gene_ixns_bulk", fake_download)
    settings = Settings(raw_data_dir=tmp_path / "raw", output_dir=tmp_path / "out")
    paths = paths_for(tmp_path / "out")
    pointer_path = paths.accepted_source_pointer(SOURCE_KEY)

    with caplog.at_level(logging.WARNING, logger="gene_dossier.section_6a"):
        payload = download_persist_and_pin_ctd_bulk(
            paths=paths,
            dossier_run_id="blocked-promote-run",
            settings=settings,
            gene_symbol="CTD_BULK",
            persist_db=False,
            promote=True,
            origin="acquire_script",
        )

    assert payload.ok is True
    assert payload.api_run_id
    assert payload.raw_artifact_id
    assert payload.attempt_dir
    assert Path(payload.attempt_dir).exists()
    assert payload.pointer is None
    assert not pointer_path.exists()
    assert load_accepted_source(paths) is None
    assert any(
        "persist_db=False" in rec.message and "skipping accepted source pin" in rec.message
        for rec in caplog.records
    )


def test_accept_source_rejects_null_provenance_ids(tmp_path):
    from gene_dossier.section_6a_sources import accept_source, paths_for

    paths = paths_for(tmp_path / "out")
    attempt = paths.new_source_attempt()
    artifact = attempt / "CTD_chem_gene_ixns.tsv.gz"
    artifact.write_bytes(b"x")
    with pytest.raises(ValueError, match="api_run_id"):
        accept_source(
            paths,
            source_key="ctd_chem_gene_ixns",
            attempt_dir=attempt,
            artifact_path=artifact,
            official_url="https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz",
            sha256=hashlib.sha256(b"x").hexdigest(),
            byte_size=1,
            validation={},
            extra={"api_run_id": None, "raw_artifact_id": None},
        )


def test_evaluate_top10_full_ranking_matches_workbook(tmp_path):
    from gene_dossier.section_6a import evaluate_section_6a_complete

    chemicals = [
        {"chemical_id": "C2", "chemical_name": "Zeta", "interaction_count": 9},
        {"chemical_id": "C1", "chemical_name": "Alpha", "interaction_count": 9},
        {"chemical_id": "C3", "chemical_name": "Beta", "interaction_count": 5},
    ]
    # Fake filtered rows so Interaction_Records pubmed sums match
    rows = []
    for chem in chemicals:
        for _ in range(chem["interaction_count"]):
            rows.append(
                {
                    "ChemicalName": chem["chemical_name"],
                    "ChemicalID": chem["chemical_id"],
                    "CasRN": "",
                    "GeneSymbol": "GENE",
                    "GeneID": "1",
                    "GeneForms": "mRNA",
                    "Organism": "Homo sapiens",
                    "OrganismID": "9606",
                    "Interaction": "i",
                    "InteractionActions": "affects",
                    "PubMedIDs": "1",
                    "source_order": len(rows),
                }
            )
    aggregated = aggregate_by_chemical_id(
        [{**r, "_pubmed_occurrence_count": 1} for r in rows]
    )
    top = rank_top_chemicals(aggregated, limit=10)
    attempt = tmp_path / "attempt"
    (attempt / "supplementary").mkdir(parents=True)
    xlsx = attempt / "supplementary" / "GENE_CTD.xlsx"
    write_ctd_workbook(
        xlsx,
        gene_symbol="GENE",
        gene_id=1,
        filtered_rows=rows,
        chemicals=aggregated,
        top_chemicals=top,
        provenance=[],
    )
    png, _ = render_top_chemicals_png(top)
    fig = attempt / "figures" / "GENE_CTD_top_chemicals.png"
    fig.parent.mkdir(parents=True)
    fig.write_bytes(png)
    status = {
        "rendering_status": {
            "scientific_status": "success",
            "visual_status": "success",
            "presentation_status": "success",
        },
        "summary": {
            "scientific_status": "success",
            "visual_status": "success",
            "presentation_status": "success",
            "filtered_target_row_count": len(rows),
            "unique_chemical_count": len(aggregated),
            "supplementary_xlsx": "GENE_CTD.xlsx",
            "supplementary_xlsx_sha256": _sha(xlsx),
            "top_chemicals_figure_relative_path": "figures/GENE_CTD_top_chemicals.png",
            "top_chemicals_figure_sha256": hashlib.sha256(png).hexdigest(),
            "top_chemicals": [
                {
                    "rank": c["display_rank"],
                    "chemical_name": c["chemical_name"],
                    "chemical_id": c["chemical_id"],
                    "interaction_count": c["interaction_count"],
                }
                for c in top
            ],
        },
    }
    html = (
        "<h2>6. Information on which perturbations affect the gene</h2>"
        "<h3>a. Comparative Toxicogenomics Database</h3>"
        '<figure class="rancho-figure section-6a-top-chemicals-figure"></figure>'
        "Comparative Toxicogenomics Database (CTD)"
    )
    pdf = tmp_path / "out.pdf"
    pdf.write_bytes(b"%PDF-1.4")
    result = evaluate_section_6a_complete(
        status=status, html_text=html, pdf_path=pdf, attempt_dir=attempt
    )
    assert result["complete"] is True
    assert result["checks"]["top_chemicals_workbook_matches_expected_ranking"]["passed"]
    assert result["checks"]["top_chemicals_summary_matches_expected_ranking"]["passed"]
