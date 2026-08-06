"""Unit tests for Harmonizome Section 4a parsing and selection."""

from __future__ import annotations

import json
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gene_dossier.models import ToolResult
from gene_dossier.tools import harmonizome as hz
from gene_dossier.tools.harmonizome_section4a import (
    CURATED_TF_DATASET_ORDER,
    CURATED_TF_DATASETS,
    PARSE_PARTIAL,
    PREDICTED_TF_DATASET_ORDER,
    PREDICTED_TF_DATASETS,
    SECTION_4A_TF_DATASETS,
    absolute_harmonizome_href,
    collect_section_4a_from_payload,
    collect_section_4a_harmonizome,
    display_dataset_label,
    parse_chea_binding,
    parse_encode_binding,
    round_robin_select,
    split_gene_set_name,
)
from gene_dossier.workflow import WorkflowTransientContext

FIXTURES = Path(__file__).parent / "fixtures" / "harmonizome"


def _payload() -> dict:
    return json.loads((FIXTURES / "gene_associations_genex.json").read_text(encoding="utf-8"))


def test_ordered_tuples_not_sets_for_iteration_contract() -> None:
    assert isinstance(CURATED_TF_DATASET_ORDER, tuple)
    assert isinstance(PREDICTED_TF_DATASET_ORDER, tuple)
    assert CURATED_TF_DATASETS == set(CURATED_TF_DATASET_ORDER)
    assert PREDICTED_TF_DATASETS == set(PREDICTED_TF_DATASET_ORDER)
    assert SECTION_4A_TF_DATASETS == CURATED_TF_DATASETS | PREDICTED_TF_DATASETS
    # Membership sets must not be used as ordered sequences.
    assert CURATED_TF_DATASET_ORDER[0].startswith("ENCODE")
    assert CURATED_TF_DATASET_ORDER[2].startswith("ChEA")
    assert PREDICTED_TF_DATASET_ORDER[0].startswith("JASPAR")


def test_split_gene_set_name_last_slash() -> None:
    attr, ds = split_gene_set_name(
        "TEAD4_HepG2_hg19_1/ENCODE Transcription Factor Binding Site Profiles"
    )
    assert attr == "TEAD4_HepG2_hg19_1"
    assert ds == "ENCODE Transcription Factor Binding Site Profiles"


def test_parse_encode_binding_and_blank_organism() -> None:
    parsed = parse_encode_binding("TEAD4_HepG2_hg19_1")
    assert parsed["association"] == "TEAD4"
    assert parsed["tissue_cells"] == "HepG2"
    assert parsed["genome_build"] == "hg19_1"
    assert parsed["organism"] is None
    assert parsed["organism_audit"] == "Human"
    assert parsed["parse_status"] == "parsed_complete"


def test_parse_chea_hyphenated_context() -> None:
    parsed = parse_chea_binding("FOXP2-12345-Frontal-Cortex-Human")
    assert parsed["association"] == "FOXP2"
    assert parsed["tissue_cells"] == "Frontal-Cortex"
    assert parsed["organism"] == "Human"
    assert parsed["genome_build"] is None
    assert parsed["pubmed_id"] == "12345"


def test_parse_chea_unknown_organism_token_not_titlecased() -> None:
    parsed = parse_chea_binding("TF-123456-T-cell-line")
    assert parsed["association"] == "TF"
    assert parsed["pubmed_id"] == "123456"
    assert parsed["tissue_cells"] == "T-cell"
    assert parsed["organism"] is None
    assert parsed["organism_audit"] is None
    assert parsed["organism_token_unparsed"] == "line"
    assert parsed["parse_status"] == PARSE_PARTIAL


def test_display_dataset_label_chea_to_chea() -> None:
    assert display_dataset_label("ChEA Transcription Factor Targets").startswith("CHEA ")
    assert display_dataset_label("ENCODE Transcription Factor Targets").startswith("ENCODE ")


def test_absolute_href() -> None:
    href = absolute_harmonizome_href("/api/1.0/gene_set/TF1/ENCODE")
    assert href == "https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/TF1/ENCODE"


def test_allowlist_excludes_out_of_scope_and_counts_by_order() -> None:
    collection = collect_section_4a_from_payload(_payload(), query_gene="GENEX")
    assert "JASPAR Predicted Transcription Factor Targets 2025" in {
        r["dataset_name"] for r in collection["out_of_scope_summary"]
    }
    assert list(collection["curated_counts"]) == list(CURATED_TF_DATASET_ORDER)
    assert list(collection["predicted_counts"]) == list(PREDICTED_TF_DATASET_ORDER)
    assert collection["curated_total"] == sum(collection["curated_counts"].values())
    assert collection["predicted_total"] == sum(collection["predicted_counts"].values())


def test_round_robin_uses_ordered_tuple_not_set() -> None:
    by_ds = {
        CURATED_TF_DATASET_ORDER[0]: [
            {"association": "A1", "dataset_name": CURATED_TF_DATASET_ORDER[0], "parse_status": "parsed_complete", "source_order": 0},
            {"association": "A2", "dataset_name": CURATED_TF_DATASET_ORDER[0], "parse_status": "parsed_complete", "source_order": 1},
        ],
        CURATED_TF_DATASET_ORDER[1]: [
            {"association": "B1", "dataset_name": CURATED_TF_DATASET_ORDER[1], "parse_status": "parsed_complete", "source_order": 2},
        ],
        CURATED_TF_DATASET_ORDER[2]: [
            {"association": "C1", "dataset_name": CURATED_TF_DATASET_ORDER[2], "parse_status": "parsed_complete", "source_order": 3},
        ],
        CURATED_TF_DATASET_ORDER[3]: [
            {"association": "D1", "dataset_name": CURATED_TF_DATASET_ORDER[3], "parse_status": "parsed_complete", "source_order": 4},
        ],
    }
    selected = round_robin_select(by_ds, CURATED_TF_DATASET_ORDER, max_rows=4)
    assert [r["association"] for r in selected] == ["A1", "B1", "C1", "D1"]


def test_curated_display_encode_organism_blank_and_chea_label() -> None:
    collection = collect_section_4a_from_payload(_payload(), query_gene="GENEX", max_displayed_curated=14)
    encode_rows = [r for r in collection["curated_display"] if "ENCODE" in r["dataset"] and "Binding" in r["dataset"]]
    assert encode_rows
    assert all(r["organism"] == "" for r in encode_rows)
    chea_rows = [r for r in collection["curated_display"] if r["dataset"].startswith("CHEA")]
    assert chea_rows


def test_collect_never_calls_gene_set() -> None:
    fake = ToolResult(
        source_name=hz.SOURCE_NAME,
        endpoint_name="gene_associations",
        gene_symbol="GENEX",
        request_url="https://maayanlab.cloud/Harmonizome/api/1.0/gene/GENEX?showAssociations=true",
        request_params={"gene_symbol": "GENEX", "showAssociations": "true"},
        success=True,
        status_code=200,
        data=_payload(),
    )
    mock_gene_set = MagicMock()
    with patch.object(hz, "gene_associations", return_value=fake):
        with patch.object(hz, "gene_set_associations", mock_gene_set):
            result = collect_section_4a_harmonizome("GENEX")
    mock_gene_set.assert_not_called()
    assert result["scientific_status"] == "success"
    assert result["collection"]["curated_total"] > 0


def test_collect_accepts_injected_gene_associations_fn() -> None:
    fake = ToolResult(
        source_name=hz.SOURCE_NAME,
        endpoint_name="gene_associations",
        gene_symbol="GENEX",
        request_url="https://example",
        request_params={},
        success=True,
        status_code=200,
        data=_payload(),
    )
    mock_fn = MagicMock(return_value=fake)
    mock_gene_set = MagicMock()
    with patch.object(hz, "gene_set_associations", mock_gene_set):
        result = collect_section_4a_harmonizome(
            "GENEX", gene_associations_fn=mock_fn
        )
    mock_fn.assert_called_once()
    mock_gene_set.assert_not_called()
    assert result["scientific_status"] == "success"


def test_concurrent_section_4a_collect_does_not_mutate_shared_client() -> None:
    fake = ToolResult(
        source_name=hz.SOURCE_NAME,
        endpoint_name="gene_associations",
        gene_symbol="GENEX",
        request_url="https://example",
        request_params={},
        success=True,
        status_code=200,
        data=_payload(),
    )
    original_gene_set = hz.gene_set_associations
    original_gene = hz.gene_associations
    barrier = threading.Barrier(4)

    def _worker() -> dict:
        barrier.wait(timeout=5)
        return collect_section_4a_harmonizome(
            "GENEX", gene_associations_fn=lambda *a, **k: fake
        )

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(_worker) for _ in range(4)]
        results = [f.result(timeout=30) for f in futures]
    assert hz.gene_set_associations is original_gene_set
    assert hz.gene_associations is original_gene
    assert all(r["scientific_status"] == "success" for r in results)


def test_request_identity_cache_shares_identical_gene_associations() -> None:
    payload = _payload()
    body = json.dumps(payload).encode("utf-8")
    transient = WorkflowTransientContext()
    calls: list[str] = []

    class _Resp:
        status_code = 200
        content = body
        url = "https://maayanlab.cloud/Harmonizome/api/1.0/gene/GENEX?showAssociations=true"
        history: list = []
        headers = {"content-type": "application/json"}
        is_success = True

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, params=None):
            calls.append(url)
            return _Resp()

    with patch("gene_dossier.tools.harmonizome.httpx.Client", _Client):
        first = hz.gene_associations("GENEX", transient=transient)
        second = hz.gene_associations("GENEX", transient=transient)
        # Different Harmonizome request identity must still execute.
        other = hz.gene_set_associations(
            "TF1",
            "ENCODE Transcription Factor Targets",
            gene_symbol="GENEX",
            transient=transient,
        )
    assert first is second
    assert len(calls) == 2
    assert first.success and other.success
    meta = first.data["_harmonizome_meta"]
    assert meta["requested_url"] != meta.get("final_url") or meta["requested_url"]
    assert meta["response_body_sha256"]
    assert meta["response_byte_length"] == len(body)
    cached = transient.get_cached_request(meta["request_identity"])
    assert cached["response_bytes"] == body
    assert cached["tool_result"] is first


def test_gene_mismatch_fail_closed() -> None:
    payload = _payload()
    payload["symbol"] = "OTHER"
    fake = ToolResult(
        source_name=hz.SOURCE_NAME,
        endpoint_name="gene_associations",
        gene_symbol="GENEX",
        request_url="https://example",
        request_params={},
        success=True,
        status_code=200,
        data=payload,
    )
    with patch.object(hz, "gene_associations", return_value=fake):
        result = collect_section_4a_harmonizome("GENEX")
    assert result["scientific_status"] == "gene_mismatch"
    assert result["collection"] is None
