"""Tests for MouseMine client input validation and identifier helpers."""

from __future__ import annotations

from gene_dossier.models import ToolResult
from gene_dossier.tools import mousemine


def test_normalize_ncbi_gene_number_accepts_int_and_digit_string():
    assert mousemine._normalize_ncbi_gene_number(20788) == "20788"
    assert mousemine._normalize_ncbi_gene_number("20788") == "20788"
    assert mousemine._normalize_ncbi_gene_number(" NCBIGene:20788 ") == "20788"


def test_normalize_ncbi_gene_number_rejects_blank_and_non_numeric():
    assert mousemine._normalize_ncbi_gene_number(None) is None
    assert mousemine._normalize_ncbi_gene_number("") is None
    assert mousemine._normalize_ncbi_gene_number("   ") is None
    assert mousemine._normalize_ncbi_gene_number("MGI:107585") is None


def test_fetch_mouse_annotations_requires_identifier():
    result = mousemine.fetch_mouse_annotations(gene_symbol="")
    assert result.success is False
    assert result.error_type == "invalid_request"
    assert "mgi_id" in (result.error_message or "")
    assert "ncbi_gene_number" in (result.error_message or "")


def test_fetch_mouse_annotations_rejects_non_numeric_ncbi():
    result = mousemine.fetch_mouse_annotations(
        ncbi_gene_number="not-a-gene-id", gene_symbol="SREBF2"
    )
    # gene_symbol is present, so symbol lookup may proceed; force no-symbol path:
    result = mousemine.fetch_mouse_annotations(ncbi_gene_number="not-a-gene-id")
    assert result.success is False
    assert result.error_type == "invalid_request"
    assert "numeric mouse Entrez ID" in (result.error_message or "")
    assert "not-a-gene-id" in (result.error_message or "")


def test_fetch_mouse_annotations_prefers_ncbi_when_mgi_missing(monkeypatch):
    """When mgi_id is missing, ncbi_gene_number must drive gene lookup."""
    calls: list[str | int] = []

    def fake_gene_lookup(ncbi_gene_number, *, gene_symbol="", settings=None):
        calls.append(ncbi_gene_number)
        return ToolResult(
            source_name="MouseMine",
            endpoint_name="gene_lookup",
            success=True,
            gene_symbol=gene_symbol or str(ncbi_gene_number),
            request_url="https://example.test/mousemine",
            request_params={"ncbi_gene_number": str(ncbi_gene_number)},
            status_code=200,
            data={
                "results": [
                    [
                        "MGI:107585",
                        "Srebf2",
                        "sterol regulatory element binding factor 2",
                        "Mus musculus/domesticus",
                        "20788",
                    ]
                ],
                "views": list(mousemine.GENE_LOOKUP_VIEWS),
            },
        )

    def fake_alleles(mgi_id, *, gene_symbol="", settings=None):
        return ToolResult(
            source_name="MouseMine",
            endpoint_name="alleles",
            success=True,
            gene_symbol=gene_symbol or mgi_id,
            request_url="https://example.test/alleles",
            request_params={"mgi_id": mgi_id},
            status_code=200,
            data={"results": [], "views": list(mousemine.ALLELE_VIEWS)},
        )

    monkeypatch.setattr(mousemine, "gene_lookup", fake_gene_lookup)
    monkeypatch.setattr(mousemine, "alleles", fake_alleles)
    monkeypatch.setattr(mousemine, "allele_phenotypes", fake_alleles)
    monkeypatch.setattr(mousemine, "stocks_carried_by", fake_alleles)

    result = mousemine.fetch_mouse_annotations(
        ncbi_gene_number=20788, mgi_id=None, gene_symbol="SREBF2"
    )
    assert calls == [20788] or calls == ["20788"]
    assert result.success is True
    assert result.data["mgi_id"] == "MGI:107585"
    assert result.data["ncbi_gene_number"] == "20788"


def test_prefer_mouse_mgi_id_filters_mus_musculus():
    rows = [
        {
            "Gene.primaryIdentifier": "100037309",
            "Gene.symbol": "srebf2",
            "Gene.organism.name": "Danio rerio",
            "Gene.ncbiGeneNumber": "100037309",
        },
        {
            "Gene.primaryIdentifier": "MGI:107585",
            "Gene.symbol": "Srebf2",
            "Gene.organism.name": "Mus musculus/domesticus",
            "Gene.ncbiGeneNumber": "20788",
        },
    ]
    assert mousemine.prefer_mouse_mgi_id(rows) == "MGI:107585"
