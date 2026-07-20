"""Tests for content-addressed raw artifact storage (no network required)."""

from __future__ import annotations

import json
from pathlib import Path

from gene_dossier.raw_store import RawStore, compute_hash


def test_compute_hash_is_sha256_hex():
    digest = compute_hash(b"hello")
    assert len(digest) == 64
    assert digest == compute_hash(b"hello")
    assert digest != compute_hash(b"hello!")


def test_save_json_round_trip_and_layout(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    data = {"idlist": ["6721"], "count": 1}
    art = store.save_json(
        "run1",
        "NCBI Gene",
        data,
        original_url="https://eutils.ncbi.nlm.nih.gov/example",
        filename_hint="esearch",
    )

    rel = Path(art.file_path).relative_to(tmp_path)
    assert rel.parts[0] == "run1"
    assert rel.parts[1] == "ncbi-gene"
    assert rel.name.startswith("esearch-")
    assert rel.suffix == ".json"

    assert art.artifact_type == "json"
    assert art.source_name == "NCBI Gene"
    assert art.original_url is not None
    assert Path(art.file_path).exists()
    assert store.load_json(art) == data


def test_json_hash_is_canonical_sorted_keys(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    a = store.save_json("run1", "NCBI Gene", {"b": 2, "a": 1}, filename_hint="x")
    b = store.save_json("run1", "NCBI Gene", {"a": 1, "b": 2}, filename_hint="x")
    assert a.file_path == b.file_path
    assert a.content_hash == b.content_hash

    expected = compute_hash(
        json.dumps({"a": 1, "b": 2}, sort_keys=True, ensure_ascii=False, indent=2).encode(
            "utf-8"
        )
    )
    assert a.content_hash == expected


def test_idempotent_same_content_same_path(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    a = store.save_json("run1", "UniProt", {"accession": "Q12772"})
    b = store.save_json("run1", "UniProt", {"accession": "Q12772"})
    assert a.file_path == b.file_path
    assert a.content_hash == b.content_hash


def test_different_content_different_path(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    a = store.save_json("run1", "UniProt", {"accession": "Q12772"})
    b = store.save_json("run1", "UniProt", {"accession": "P42858"})
    assert a.file_path != b.file_path
    assert a.content_hash != b.content_hash


def test_save_text_tsv(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    tsv = "ChemicalName\tGeneSymbol\nstatin\tSREBF2\n"
    art = store.save_text(
        "run1",
        "CTD",
        tsv,
        extension="tsv",
        artifact_type="tsv",
        filename_hint="batch",
    )
    assert Path(art.file_path).suffix == ".tsv"
    assert art.artifact_type == "tsv"
    assert store.load_text(art) == tsv


def test_save_bytes_round_trip(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    payload = b"\x00\x01\xffPNG"
    art = store.save_bytes(
        "run1",
        "AlphaFold",
        payload,
        extension="png",
        artifact_type="image",
    )
    assert Path(art.file_path).suffix == ".png"
    assert art.artifact_type == "image"
    assert store.load_bytes(art) == payload


def test_verify_detects_missing_and_tampering(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    art = store.save_json("run1", "NCBI Gene", {"ok": True})
    assert store.verify(art) is True

    Path(art.file_path).write_bytes(b"tampered")
    assert store.verify(art) is False

    Path(art.file_path).unlink()
    assert store.verify(art) is False


def test_api_run_id_and_notes_are_recorded(tmp_path: Path):
    store = RawStore(base_dir=tmp_path)
    art = store.save_json(
        "run1",
        "PubMed",
        {"ids": ["1"]},
        api_run_id="api-abc",
        notes="esearch sample",
    )
    assert art.api_run_id == "api-abc"
    assert art.notes == "esearch sample"
    assert art.dossier_run_id == "run1"
