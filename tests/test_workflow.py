"""Tests for LangGraph dossier workflow orchestration (offline / no network)."""

from __future__ import annotations

from pathlib import Path

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult
from gene_dossier.workflow import (
    _client_opentargets,
    _client_reactome,
    extract_gene_ids_from_tool_result,
    node_resolve_gene_identity,
    run_gene_dossier_full_api_pass,
)


def _ncbi_tool_result() -> ToolResult:
    return ToolResult(
        source_name="NCBI Gene",
        endpoint_name="lookup_gene",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/ncbi",
        data={
            "gene_symbol": "SREBF2",
            "selected_gene_id": "6721",
            "selection_method": "exact_symbol",
            "selected_summary": {
                "nomenclaturesymbol": "SREBF2",
                "chromosome": "22",
                "uid": "6721",
            },
            "candidate_ids": ["6721"],
        },
    )


def _ensembl_tool_result() -> ToolResult:
    return ToolResult(
        source_name="Ensembl",
        endpoint_name="lookup_symbol",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/ensembl",
        data={
            "ensembl_id": "ENSG00000198911",
            "gene_symbol": "SREBF2",
        },
    )


def _uniprot_tool_result() -> ToolResult:
    return ToolResult(
        source_name="UniProt",
        endpoint_name="search_reviewed",
        success=True,
        gene_symbol="SREBF2",
        request_url="https://example.test/uniprot",
        data={
            "selected_accession": "Q12772",
            "gene_symbol": "SREBF2",
        },
    )


def test_client_reactome_requires_uniprot_accession():
    settings = get_settings()
    missing = _client_reactome(
        gene_symbol="SREBF2", gene_ids={}, settings=settings
    )
    assert missing.success is False
    assert missing.error_type == "missing_identifier"

    # With accession present, do not claim missing_identifier (network not required here).
    # We only assert the guard path; live Reactome call is out of scope for unit tests.


def test_client_opentargets_requires_ensembl_id():
    settings = get_settings()
    missing = _client_opentargets(
        gene_symbol="SREBF2", gene_ids={}, settings=settings
    )
    assert missing.success is False
    assert missing.error_type == "missing_identifier"
    assert "ensembl_id" in (missing.error_message or "")


def test_extract_gene_ids_from_identity_tool_results():
    gene_ids: dict = {}
    for tr in (_ncbi_tool_result(), _ensembl_tool_result(), _uniprot_tool_result()):
        gene_ids = extract_gene_ids_from_tool_result(tr, gene_ids)
    assert gene_ids["entrez_gene_id"] == "6721"
    assert gene_ids["chromosome"] == "22"
    assert gene_ids["ensembl_id"] == "ENSG00000198911"
    assert gene_ids["uniprot_accession"] == "Q12772"
    assert gene_ids["official_symbol"] == "SREBF2"


def test_preloaded_identity_results_populate_gene_ids_offline():
    """call_network=False must still harvest IDs from preloaded identity ToolResults."""
    settings = Settings()
    state = {
        "gene_symbol": "SREBF2",
        "dossier_run_id": "offline-id-1",
        "gene_ids": {},
        "tool_results": [
            _ncbi_tool_result(),
            _ensembl_tool_result(),
            _uniprot_tool_result(),
        ],
        "errors": [],
    }
    updated = node_resolve_gene_identity(
        state, settings=settings, call_network=False  # type: ignore[arg-type]
    )
    ids = updated["gene_ids"]
    assert ids["entrez_gene_id"] == "6721"
    assert ids["ensembl_id"] == "ENSG00000198911"
    assert ids["uniprot_accession"] == "Q12772"
    assert ids["chromosome"] == "22"
    assert updated.get("official_symbol") == "SREBF2"


def test_offline_srebf2_smoke_persist_db_false(tmp_path: Path):
    result = run_gene_dossier_full_api_pass(
        "SREBF2",
        call_network=False,
        preloaded_tool_results=[
            _ncbi_tool_result(),
            _ensembl_tool_result(),
            _uniprot_tool_result(),
        ],
        output_dir=tmp_path,
        write_pdf=False,
        persist_db=False,
        force_deterministic=True,
        dossier_run_id="wf-offline-srebf2",
    )
    assert result.status == "completed"
    assert result.gene_ids.get("entrez_gene_id") == "6721"
    assert result.gene_ids.get("ensembl_id") == "ENSG00000198911"
    assert result.gene_ids.get("uniprot_accession") == "Q12772"
    assert result.output_paths.get("rancho_html")
    assert result.output_paths.get("coverage_markdown")
    assert result.output_paths.get("debug_markdown")
    assert result.evidence_records  # NCBI identity normalized
