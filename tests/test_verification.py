"""Tests for rule-based claim verification (no network / no LLM)."""

from __future__ import annotations

from gene_dossier.models import (
    AssertionType,
    Claim,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.source_ids import make_source_id
from gene_dossier.verification import (
    CAUSAL_PHRASES,
    find_causal_language,
    grade_rank,
    has_disease_evidence,
    index_evidence_by_source_id,
    resolve_claim_grade,
    strongest_grade,
    verify_claim,
    verify_claims,
)


def _evidence(
    *,
    source_name: str = "NCBI Gene",
    assertion_type: AssertionType = AssertionType.gene_identity,
    fact_type: str = "entrez_gene_id",
    grade: EvidenceGrade = EvidenceGrade.A,
    key: str = "6721",
    section: str = "General",
    source_type: SourceType = SourceType.curated_database,
    manual_review_required: bool = False,
) -> EvidenceRecord:
    sid = make_source_id(source_name, "SREBF2", assertion_type, key)
    return EvidenceRecord(
        source_id=sid,
        dossier_run_id="run1",
        gene_symbol="SREBF2",
        section=section,
        source_name=source_name,
        source_type=source_type,
        assertion_type=assertion_type,
        fact_type=fact_type,
        evidence_grade=grade,
        display_text="test evidence",
        value={},
        manual_review_required=manual_review_required,
    )


def _claim(
    text: str,
    source_ids: list[str] | None = None,
    *,
    grade: EvidenceGrade | None = None,
    claim_id: str | None = None,
) -> Claim:
    kwargs: dict = {
        "dossier_run_id": "run1",
        "claim_text": text,
        "source_ids": source_ids or [],
        "evidence_grade": grade,
    }
    if claim_id is not None:
        kwargs["id"] = claim_id
    return Claim(**kwargs)


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
def test_find_causal_language_detects_phrases():
    hits = find_causal_language("SREBF2 causes HD and is a therapeutic target.")
    assert "causes" in hits
    assert "therapeutic target" in hits
    # Empty / non-causal
    assert find_causal_language("") == []
    assert find_causal_language("SREBF2 is expressed in brain.") == []


def test_find_causal_language_covers_mvp_phrases():
    for phrase in CAUSAL_PHRASES:
        assert phrase in find_causal_language(f"Gene {phrase} disease.")


def test_grade_rank_and_strongest_grade():
    assert grade_rank(EvidenceGrade.A) > grade_rank(EvidenceGrade.B)
    assert grade_rank(None) == 0
    assert strongest_grade([EvidenceGrade.C, EvidenceGrade.A, EvidenceGrade.F]) is EvidenceGrade.A
    assert strongest_grade([]) is None
    assert strongest_grade([None, None]) is None


def test_index_evidence_keeps_first_source_id():
    a = _evidence(key="1")
    b = a.model_copy(update={"display_text": "duplicate sid"})
    indexed = index_evidence_by_source_id([a, b])
    assert indexed[a.source_id] is a


def test_resolve_claim_grade_prefers_claim_then_evidence():
    ev = _evidence(grade=EvidenceGrade.C)
    claim = _claim("ok", [ev.source_id], grade=EvidenceGrade.B)
    assert resolve_claim_grade(claim, [ev]) is EvidenceGrade.B
    claim_no_grade = _claim("ok", [ev.source_id])
    assert resolve_claim_grade(claim_no_grade, [ev]) is EvidenceGrade.C
    assert resolve_claim_grade(claim_no_grade, []) is None


def test_has_disease_evidence():
    gene = _evidence()
    disease = _evidence(
        source_name="OMIM",
        assertion_type=AssertionType.disease_association,
        fact_type="omim_phenotype_map",
        key="123",
        section="ClinVar",
        source_type=SourceType.genetic_database,
    )
    variant = _evidence(
        source_name="ClinVar",
        assertion_type=AssertionType.variant_association,
        fact_type="clinvar_variant",
        key="456",
        section="ClinVar",
        source_type=SourceType.genetic_database,
    )
    assert has_disease_evidence([gene]) is False
    assert has_disease_evidence([disease]) is True
    assert has_disease_evidence([variant]) is True


# --------------------------------------------------------------------------------------
# Source presence / existence
# --------------------------------------------------------------------------------------
def test_missing_source_ids_fails():
    result = verify_claim(_claim("SREBF2 is a transcription factor."), [])
    assert result.verdict == "fail"
    assert result.source_id_presence_passed is False
    assert result.source_exists_passed is False
    assert result.semantic_support == "fail"
    assert "no source_id" in (result.reason or "").lower()


def test_missing_cited_evidence_fails():
    sid = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "missing")
    result = verify_claim(_claim("SREBF2 is annotated.", [sid]), [])
    assert result.verdict == "fail"
    assert result.source_id_presence_passed is True
    assert result.source_exists_passed is False
    assert "missing from evidence" in (result.reason or "").lower()


def test_valid_citation_passes():
    ev = _evidence()
    result = verify_claim(_claim("SREBF2 encodes a transcription factor.", [ev.source_id]), [ev])
    assert result.verdict == "pass"
    assert result.source_id_presence_passed is True
    assert result.source_exists_passed is True
    assert result.needs_human_review is False
    assert result.causal_language_check == "pass"
    assert result.evidence_strength_check == "pass"


# --------------------------------------------------------------------------------------
# Causal overclaim
# --------------------------------------------------------------------------------------
def test_causal_no_disease_evidence_human_review_even_if_grade_a():
    """Grade A gene identity alone is not enough for causal language."""
    ev = _evidence(grade=EvidenceGrade.A)
    result = verify_claim(
        _claim("SREBF2 causes Huntington disease.", [ev.source_id], grade=EvidenceGrade.A),
        [ev],
    )
    assert result.verdict == "human_review"
    assert result.needs_human_review is True
    assert result.causal_language_check == "warning"
    assert result.evidence_strength_check == "warning"
    assert "no cited disease" in (result.reason or "").lower()


def test_causal_disease_evidence_below_a_human_review():
    ev = _evidence(
        source_name="Open Targets",
        assertion_type=AssertionType.disease_association,
        fact_type="ot_association",
        grade=EvidenceGrade.C,
        key="EFO_0000001",
        section="ClinVar",
        source_type=SourceType.genetic_database,
    )
    result = verify_claim(
        _claim("SREBF2 drives disease progression.", [ev.source_id]),
        [ev],
    )
    assert result.verdict == "human_review"
    assert result.needs_human_review is True
    assert result.evidence_strength_check == "warning"
    assert "below A" in (result.reason or "")


def test_causal_grade_a_disease_evidence_is_warning_only():
    ev = _evidence(
        source_name="OMIM",
        assertion_type=AssertionType.disease_association,
        fact_type="omim_phenotype_map",
        grade=EvidenceGrade.A,
        key="123",
        section="ClinVar",
        source_type=SourceType.genetic_database,
    )
    result = verify_claim(
        _claim("This variant is pathogenic.", [ev.source_id], grade=EvidenceGrade.A),
        [ev],
    )
    assert result.verdict == "warning"
    assert result.needs_human_review is False
    assert result.causal_language_check == "warning"
    assert result.evidence_strength_check == "pass"


def test_fail_takes_precedence_over_causal_human_review():
    result = verify_claim(_claim("Gene causes disease.", []), [])
    assert result.verdict == "fail"
    assert result.needs_human_review is True
    assert result.causal_language_check == "warning"
    assert result.evidence_strength_check == "warning"


def test_variant_association_counts_as_disease_evidence():
    ev = _evidence(
        source_name="ClinVar",
        assertion_type=AssertionType.variant_association,
        fact_type="clinvar_variant",
        grade=EvidenceGrade.A,
        key="rs1",
        section="ClinVar",
        source_type=SourceType.genetic_database,
    )
    result = verify_claim(
        _claim("Variant is pathogenic for HD.", [ev.source_id], grade=EvidenceGrade.A),
        [ev],
    )
    assert result.verdict == "warning"
    assert result.needs_human_review is False


# --------------------------------------------------------------------------------------
# Manual review flag on evidence
# --------------------------------------------------------------------------------------
def test_manual_review_evidence_produces_warning():
    ev = _evidence(manual_review_required=True)
    result = verify_claim(
        _claim("SREBF2 is linked to lipid metabolism.", [ev.source_id]),
        [ev],
    )
    assert result.verdict == "warning"
    assert result.needs_human_review is True
    assert "manual review" in (result.reason or "").lower()


def test_manual_review_does_not_override_fail():
    sid = make_source_id("NCBI Gene", "SREBF2", "gene_identity", "gone")
    # Presence ok but existence fails; causal also present.
    result = verify_claim(_claim("Gene causes disease.", [sid]), [])
    assert result.verdict == "fail"


# --------------------------------------------------------------------------------------
# Batch + dict index path
# --------------------------------------------------------------------------------------
def test_verify_claims_batch_and_dict_index():
    ok = _evidence(key="ok")
    disease = _evidence(
        source_name="OMIM",
        assertion_type=AssertionType.disease_association,
        fact_type="omim_phenotype_map",
        grade=EvidenceGrade.A,
        key="pheno",
        section="ClinVar",
        source_type=SourceType.genetic_database,
    )
    claims = [
        _claim("SREBF2 encodes SREBP2.", [ok.source_id], claim_id="c-pass"),
        _claim("SREBF2 causes HD.", [ok.source_id], claim_id="c-review"),
        _claim("Variant is pathogenic.", [disease.source_id], grade=EvidenceGrade.A, claim_id="c-warn"),
        _claim("Unsupported.", [], claim_id="c-fail"),
    ]
    results = verify_claims(claims, [ok, disease])
    by_claim = {r.claim_id: r for r in results}
    assert by_claim["c-pass"].verdict == "pass"
    assert by_claim["c-review"].verdict == "human_review"
    assert by_claim["c-warn"].verdict == "warning"
    assert by_claim["c-fail"].verdict == "fail"

    # Pre-built dict index is accepted by verify_claim.
    indexed = index_evidence_by_source_id([ok])
    r = verify_claim(_claim("SREBF2 is annotated.", [ok.source_id]), indexed)
    assert r.verdict == "pass"
