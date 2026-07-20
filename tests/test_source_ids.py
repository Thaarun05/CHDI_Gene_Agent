"""Tests for deterministic source_id generation (no network required)."""

from __future__ import annotations

import pytest

from gene_dossier.models import AssertionType
from gene_dossier.source_ids import (
    is_valid_source_id,
    make_source_id,
    parse_source_id,
    slugify,
)


# --------------------------------------------------------------------------------------
# slugify
# --------------------------------------------------------------------------------------
def test_slugify_basic():
    assert slugify("NCBI Gene") == "ncbi-gene"
    assert slugify("PMID:12345678") == "pmid-12345678"
    assert slugify("  Homo sapiens!!  ") == "homo-sapiens"
    assert slugify("") == ""


def test_slugify_allow_underscore_preserves_underscores():
    assert slugify("gene_identity", allow_underscore=True) == "gene_identity"
    # Spaces still become hyphens even when underscores are kept.
    assert slugify("NCBI Gene", allow_underscore=True) == "ncbi-gene"


def test_slugify_uses_enum_value():
    assert slugify(AssertionType.ppi, allow_underscore=True) == "ppi"


# --------------------------------------------------------------------------------------
# make_source_id: determinism / idempotency
# --------------------------------------------------------------------------------------
def test_make_source_id_shape():
    sid = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6721")
    assert sid == "ncbi-gene:srebf2:gene_identity:6721"


def test_make_source_id_is_deterministic():
    a = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6721")
    b = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6721")
    assert a == b


def test_make_source_id_case_and_enum_invariant():
    a = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6721")
    b = make_source_id("ncbi gene", "srebf2", AssertionType.gene_identity, "6721")
    assert a == b


def test_make_source_id_multiple_key_parts():
    sid = make_source_id("STRING", "SREBF2", "ppi", "SREBF2", "SCAP")
    assert sid == "string:srebf2:ppi:srebf2-scap"


# --------------------------------------------------------------------------------------
# make_source_id: hashing fallbacks
# --------------------------------------------------------------------------------------
def test_long_key_falls_back_to_stable_hash():
    long_key = "x" * 100
    a = make_source_id("STRING", "SREBF2", "ppi", long_key)
    b = make_source_id("STRING", "SREBF2", "ppi", long_key)
    assert a == b
    key = parse_source_id(a)["key"]
    assert key.startswith("h-")
    assert len(key) == len("h-") + 10


def test_symbol_only_key_falls_back_to_hash():
    sid = make_source_id("CTD", "SREBF2", "chemical_interaction", "%%%")
    assert parse_source_id(sid)["key"].startswith("h-")


def test_different_keys_produce_different_ids():
    a = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6721")
    b = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6722")
    assert a != b


# --------------------------------------------------------------------------------------
# validation / parsing
# --------------------------------------------------------------------------------------
def test_is_valid_source_id():
    sid = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6721")
    assert is_valid_source_id(sid)
    assert not is_valid_source_id("a:b")
    assert not is_valid_source_id("a:b:c:")  # blank trailing component
    assert not is_valid_source_id("")
    assert not is_valid_source_id(None)  # type: ignore[arg-type]


def test_parse_source_id_round_trip():
    sid = make_source_id("PubMed", "SREBF2", "literature_summary", "PMID:12345678")
    parsed = parse_source_id(sid)
    assert parsed == {
        "source": "pubmed",
        "gene": "srebf2",
        "assertion": "literature_summary",
        "key": "pmid-12345678",
    }


def test_parse_invalid_raises():
    with pytest.raises(ValueError):
        parse_source_id("not-a-valid-id")


# --------------------------------------------------------------------------------------
# errors on empty components
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "args",
    [
        ("", "SREBF2", "gene_identity", "6721"),
        ("NCBI Gene", "", "gene_identity", "6721"),
        ("NCBI Gene", "SREBF2", "", "6721"),
        ("NCBI Gene", "SREBF2", "gene_identity", ""),
    ],
)
def test_empty_component_raises(args):
    with pytest.raises(ValueError):
        make_source_id(*args)
