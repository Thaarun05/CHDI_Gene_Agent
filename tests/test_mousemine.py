"""Tests for MouseMine client input validation and identifier helpers."""

from __future__ import annotations

from gene_dossier.config import Settings
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


def test_fetch_mouse_annotations_uses_validated_request_contract(monkeypatch):
    payloads = [
        {
            "results": [
                [
                    "MGI:107585",
                    "Srebf2",
                    "sterol regulatory element binding factor 2",
                    "Mus musculus/domesticus",
                    "20788",
                ]
            ]
        },
        {"results": []},
        {"results": []},
        {"results": []},
    ]
    calls: list[tuple[str, dict[str, str]]] = []

    class _Response:
        status_code = 200
        text = ""
        is_success = True

        def __init__(self, payload):
            self.payload = payload

        def json(self):
            return self.payload

    class _Client:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, url, params=None):
            calls.append((url, dict(params or {})))
            return _Response(payloads.pop(0))

    monkeypatch.setattr(mousemine.httpx, "Client", _Client)
    result = mousemine.fetch_mouse_annotations(
        ncbi_gene_number=20788,
        gene_symbol="SREBF2",
        settings=Settings(http_timeout_seconds=1),
    )

    assert result.success
    assert result.source_name == "MouseMine"
    assert result.endpoint_name == "fetch_mouse_annotations"
    assert result.request_params["ncbi_gene_number"] == "20788"
    assert result.request_params["mgi_id"] == "MGI:107585"
    assert len(calls) == 4
    assert all(url == mousemine.MOUSEMINE_RESULTS_URL for url, _ in calls)
    assert all(params["format"] == "json" for _, params in calls)
    assert 'constraint path="Gene.ncbiGeneNumber" op="=" value="20788"' in calls[0][1]["query"]
    assert all(
        'constraint path="Allele.feature.primaryIdentifier" op="=" value="MGI:107585"'
        in params["query"]
        for _, params in calls[1:]
    )
    assert result.request_url.startswith(mousemine.MOUSEMINE_RESULTS_URL)
