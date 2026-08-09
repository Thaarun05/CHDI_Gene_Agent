"""Offline tests for CTD bulk helpers used by Section 6a."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from gene_dossier.tools import ctd as ctd_client

FIXTURE = Path(__file__).parent / "fixtures" / "ctd" / "CTD_chem_gene_ixns_mini.tsv.gz"


def test_legacy_batch_query_contract_preserved():
    assert ctd_client.BATCH_QUERY_URL.endswith("batchQuery.go")
    assert ctd_client.DEFAULT_INPUT_TYPE == "gene"
    assert ctd_client.DEFAULT_REPORT == "cgixns"
    assert "Input" in ctd_client.EXPECTED_COLUMNS
    assert callable(ctd_client.batch_query)
    assert callable(ctd_client.fetch_chemical_gene_interactions)


def test_bulk_header_required_columns_open_schema():
    ok = ctd_client.validate_bulk_header(list(ctd_client.BULK_REQUIRED_COLUMNS) + ["ExtraCol"])
    assert ok["ok"] is True
    assert ok["extra_columns"] == ["ExtraCol"]
    missing = ctd_client.validate_bulk_header(["ChemicalName", "GeneID"])
    assert missing["ok"] is False
    assert "ChemicalID" in missing["missing_required_columns"]


def test_iter_bulk_tsv_rows_fixture_metadata_and_extras():
    content = FIXTURE.read_bytes()
    meta, rows = ctd_client.iter_bulk_tsv_rows(content)
    assert meta["ok"] is True
    assert meta["ctd_report_created"]
    assert "ExtraCol" in meta["extra_columns"]
    assert "GeneForms" in meta["fieldnames"]
    all_rows = list(rows)
    assert len(all_rows) >= 10
    assert any(r.get("ChemicalName") == "Valproic Acid" for r in all_rows)


def test_iter_bulk_rejects_bad_gzip():
    with pytest.raises(ValueError, match="gzip"):
        ctd_client.iter_bulk_tsv_rows(b"not-a-gzip")


def test_iter_bulk_comment_fields_header():
    body = (
        "# Report created: Thu Jul 30 14:03:48 EDT 2026\n"
        "# Fields:\n"
        "# ChemicalName\tChemicalID\tCasRN\tGeneSymbol\tGeneID\tGeneForms\t"
        "Organism\tOrganismID\tInteraction\tInteractionActions\tPubMedIDs\n"
        "#\n"
        "ChemA\tC1\t\tGENE\t1\tmRNA\tHomo sapiens\t9606\ti\taffects\t11\n"
    )
    gz = gzip.compress(body.encode("utf-8"))
    meta, rows = ctd_client.iter_bulk_tsv_rows(gz)
    assert meta["header_origin"] == "comment_fields"
    assert meta["ctd_report_created"] == "Thu Jul 30 14:03:48 EDT 2026"
    data = list(rows)
    assert len(data) == 1
    assert data[0]["ChemicalID"] == "C1"
    assert data[0]["GeneID"] == "1"


def test_iter_bulk_rejects_missing_required_column():
    body = (
        "# Report created: test\n"
        "ChemicalName\tChemicalID\tGeneSymbol\tGeneID\n"
        "x\tC1\tGENE\t1\n"
    )
    gz = gzip.compress(body.encode("utf-8"))
    with pytest.raises(ValueError, match="missing_required_columns"):
        ctd_client.iter_bulk_tsv_rows(gz)
