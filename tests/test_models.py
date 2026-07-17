"""Tests for the core domain models and enums (no network required)."""

from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

from gene_dossier import models as m


# --------------------------------------------------------------------------------------
# Enums
# --------------------------------------------------------------------------------------
def test_enum_member_counts_match_spec():
    assert [g.value for g in m.EvidenceGrade] == ["A", "B", "C", "D", "E", "F"]
    assert len(list(m.SourceType)) == 15
    assert len(list(m.AssertionType)) == 18
    assert len(list(m.SourceStatus)) == 8


def test_enums_are_string_valued():
    assert m.EvidenceGrade.A == "A"
    assert m.SourceStatus.requires_key.value == "requires_key"
    assert m.SourceType.curated_database == "curated_database"


# --------------------------------------------------------------------------------------
# Defaults / helpers
# --------------------------------------------------------------------------------------
def test_new_id_is_unique():
    assert m.new_id() != m.new_id()


def test_utcnow_is_timezone_aware():
    assert m.utcnow().tzinfo is not None


def test_dossier_run_defaults():
    run = m.DossierRun(gene_symbol="SREBF2")
    assert run.id  # auto-generated
    assert run.run_type == "full_dossier"
    assert run.status == "created"
    assert run.completed_at is None
    assert run.started_at.tzinfo is not None


def test_mutable_defaults_not_shared():
    a = m.ReportSection(dossier_run_id="r1", section_name="s")
    b = m.ReportSection(dossier_run_id="r1", section_name="s")
    a.source_ids.append("x")
    assert b.source_ids == []  # separate list instances


# --------------------------------------------------------------------------------------
# Provenance chain
# --------------------------------------------------------------------------------------
def test_evidence_record_full_chain():
    run = m.DossierRun(gene_symbol="SREBF2")
    ev = m.EvidenceRecord(
        source_id="ncbi_gene:6721:identity",
        dossier_run_id=run.id,
        gene_symbol="SREBF2",
        section="General gene information",
        source_name="NCBI Gene",
        source_type=m.SourceType.curated_database,
        assertion_type=m.AssertionType.gene_identity,
        fact_type="entrez_id",
        organism="Homo sapiens",
        taxon_id=9606,
        evidence_grade=m.EvidenceGrade.C,
        value={"entrez_id": "6721"},
        display_text="SREBF2 Entrez Gene ID is 6721.",
        api_run_id="api1",
        raw_artifact_id="art1",
    )
    assert ev.evidence_grade is m.EvidenceGrade.C
    assert ev.assertion_type is m.AssertionType.gene_identity
    assert ev.value["entrez_id"] == "6721"
    assert ev.manual_review_required is False


def test_tool_result_captures_failure_without_raising():
    tr = m.ToolResult(
        source_name="biogrid",
        endpoint_name="interactions",
        success=False,
        gene_symbol="SREBF2",
        request_url="https://webservice.thebiogrid.org/interactions",
        error_type="requires_key",
        error_message="BIOGRID_ACCESSKEY missing",
    )
    assert tr.success is False
    assert tr.data is None
    assert tr.error_type == "requires_key"


def test_source_coverage_result_defaults():
    cov = m.SourceCoverageResult(
        dossier_run_id="r1",
        source_name="BioGRID",
        status=m.SourceStatus.requires_key,
    )
    assert cov.evidence_record_count == 0
    assert cov.report_sections_supported == []


# --------------------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------------------
def test_verification_result_defaults_and_literals():
    vr = m.VerificationResult(
        claim_id="c1",
        source_id_presence_passed=True,
        source_exists_passed=True,
    )
    assert vr.verdict == "pass"
    assert vr.semantic_support == "pass"
    assert vr.needs_human_review is False


def test_verification_rejects_invalid_verdict():
    with pytest.raises(ValidationError):
        m.VerificationResult(
            claim_id="c1",
            source_id_presence_passed=True,
            source_exists_passed=True,
            verdict="definitely",  # not in Verdict4
        )


# --------------------------------------------------------------------------------------
# Serialization / validation
# --------------------------------------------------------------------------------------
def test_json_round_trip_serializes_enums_as_values():
    ev = m.EvidenceRecord(
        source_id="s1",
        dossier_run_id="r1",
        gene_symbol="SREBF2",
        section="Pathways",
        source_name="Reactome",
        source_type=m.SourceType.pathway_database,
        assertion_type=m.AssertionType.pathway_membership,
        fact_type="pathway",
        evidence_grade=m.EvidenceGrade.C,
        display_text="SREBF2 participates in the SREBP pathway.",
    )
    dumped = ev.model_dump(mode="json")
    assert dumped["evidence_grade"] == "C"
    assert dumped["source_type"] == "pathway_database"
    # Reconstructs from JSON string.
    restored = m.EvidenceRecord.model_validate(json.loads(json.dumps(dumped)))
    assert restored.source_id == ev.source_id
    assert restored.evidence_grade is m.EvidenceGrade.C


def test_missing_required_field_raises():
    with pytest.raises(ValidationError):
        m.EvidenceRecord(  # type: ignore[call-arg]
            source_id="s1",
            dossier_run_id="r1",
            gene_symbol="SREBF2",
            # section missing
            source_name="Reactome",
            source_type=m.SourceType.pathway_database,
            assertion_type=m.AssertionType.pathway_membership,
            fact_type="pathway",
            evidence_grade=m.EvidenceGrade.C,
            display_text="x",
        )
