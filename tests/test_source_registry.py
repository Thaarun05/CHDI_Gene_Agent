"""Tests for the full source registry (no network required)."""

from __future__ import annotations

from gene_dossier.models import SourceStatus, SourceType
from gene_dossier.source_registry import (
    SourceDefinition,
    SourcePriority,
    get_all_sources,
    get_source,
    get_sources_by_priority,
    iter_sources,
    list_source_names,
    registry_summary,
)


EXPECTED_A = {
    "NCBI Gene",
    "PubMed",
    "UniProt",
    "Ensembl",
    "GTEx",
    "STRING",
    "Reactome",
    "ClinVar",
    "Open Targets",
    "MouseMine",
    "CTD",
    "ChEMBL",
    "PubChem",
    "NIH RePORTER",
}

EXPECTED_B = {
    "GEO",
    "Harmonizome",
    "BioGRID",
    "WikiPathways",
    "AlphaFold",
    "PDBe",
    "CDD",
    "NCBI Datasets",
    "UCSC",
}

EXPECTED_C = {
    "Allen Brain Atlas",
    "BrainRNASeq",
    "Patents",
    "Antibodies",
    "OMIM",
    "DrugBank",
    "NCATS",
    "ERC Grants",
    "HDinHD",
}


def test_registry_counts_and_priorities():
    sources = get_all_sources()
    assert len(sources) == 32
    assert {s.name for s in get_sources_by_priority(SourcePriority.A)} == EXPECTED_A
    assert {s.name for s in get_sources_by_priority("B")} == EXPECTED_B
    assert {s.name for s in get_sources_by_priority("C")} == EXPECTED_C


def test_no_duplicate_source_names():
    names = list_source_names()
    assert len(names) == len(set(names))
    assert names == [s.name for s in get_all_sources()]


def test_get_source_case_insensitive():
    src = get_source("ncbi gene")
    assert src is not None
    assert src.name == "NCBI Gene"
    assert src.priority is SourcePriority.A
    assert src.source_type is SourceType.curated_database
    assert "NCBI_API_KEY" in src.optional_keys
    assert get_source("does-not-exist") is None


def test_required_key_sources():
    biogrid = get_source("BioGRID")
    assert biogrid is not None
    assert biogrid.required_keys == ["BIOGRID_ACCESSKEY"]
    assert biogrid.default_status is SourceStatus.requires_key

    omim = get_source("OMIM")
    assert omim is not None
    assert omim.required_keys == ["OMIM_API_KEY"]
    assert omim.default_status is SourceStatus.requires_key

    patents = get_source("Patents")
    assert patents is not None
    assert patents.required_keys == ["SERPAPI_API_KEY"]
    assert patents.default_status is SourceStatus.requires_key


def test_manual_and_deferred_sources():
    antibodies = get_source("Antibodies")
    assert antibodies is not None
    assert antibodies.default_status is SourceStatus.manual

    for name in (
        "Allen Brain Atlas",
        "BrainRNASeq",
        "DrugBank",
        "NCATS",
        "ERC Grants",
        "HDinHD",
    ):
        src = get_source(name)
        assert src is not None
        assert src.default_status is SourceStatus.deferred


def test_hdinhd_architecture_placeholder():
    src = get_source("HDinHD")
    assert src is not None
    assert src.source_type is SourceType.hd_specific_database
    assert src.priority is SourcePriority.C
    assert src.notes is not None
    assert "MCP" in src.notes


def test_clients_not_marked_implemented_yet():
    # Phase 4 only declares the map; tools/ clients come later.
    assert all(not s.client_implemented for s in get_all_sources())
    assert all(not s.normalizer_implemented for s in get_all_sources())
    assert list(iter_sources(implemented_clients_only=True)) == []


def test_priority_a_sources_have_client_and_normalizer_modules():
    for src in get_sources_by_priority(SourcePriority.A):
        assert src.client_module, f"{src.name} missing client_module"
        assert src.normalizer_module, f"{src.name} missing normalizer_module"
        assert src.report_sections, f"{src.name} missing report_sections"
        assert src.default_status is SourceStatus.not_implemented


def test_registry_summary_shape():
    summary = registry_summary()
    assert summary["total"] == 32
    assert summary["by_priority"] == {"A": 14, "B": 9, "C": 9}
    assert summary["by_default_status"]["requires_key"] == 3
    assert summary["by_default_status"]["manual"] == 1
    assert summary["by_default_status"]["deferred"] == 6
    assert summary["by_default_status"]["not_implemented"] == 22


def test_get_all_sources_returns_copies():
    a = get_all_sources()
    b = get_all_sources()
    assert a[0] is not b[0]
    assert isinstance(a[0], SourceDefinition)
    a[0].notes = "mutated"
    assert get_all_sources()[0].notes != "mutated"
