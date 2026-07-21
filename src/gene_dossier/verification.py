"""Rule-based claim verification (MVP).

Checks that every claim cites existing evidence ``source_id`` values and flags
causal overclaiming relative to evidence grades. Does **not** call the network
or an LLM.

Rules (from IMPLEMENTATION_PLAN §8):
- Every claim must cite >= 1 ``source_id``
- Every cited ``source_id`` must exist in the provided evidence records
- Flag causal language (causes, proves, drives, therapeutic target, …)
- Causal language + evidence grade below A, **or** no disease-association evidence
  -> warning / human_review
- Return :class:`~gene_dossier.models.VerificationResult` objects
"""

from __future__ import annotations

import re
from typing import Iterable

from gene_dossier.models import (
    AssertionType,
    Claim,
    EvidenceGrade,
    EvidenceRecord,
    Verdict3,
    Verdict4,
    VerificationResult,
)

# Phrases from the platform verification MVP. Matched case-insensitively as
# whole-word / phrase patterns where practical.
CAUSAL_PHRASES: tuple[str, ...] = (
    "causes",
    "cause",
    "caused",
    "proves",
    "prove",
    "proven",
    "drives",
    "drive",
    "driven",
    "therapeutic target",
    "clinically validated",
    "pathogenic",
    "disease-modifying",
    "disease modifying",
)

_GRADE_RANK: dict[EvidenceGrade, int] = {
    EvidenceGrade.A: 6,
    EvidenceGrade.B: 5,
    EvidenceGrade.C: 4,
    EvidenceGrade.D: 3,
    EvidenceGrade.E: 2,
    EvidenceGrade.F: 1,
}

_DISEASE_ASSERTIONS = {
    AssertionType.disease_association,
    AssertionType.variant_association,
}


def _compile_causal_patterns() -> list[re.Pattern[str]]:
    patterns: list[re.Pattern[str]] = []
    for phrase in CAUSAL_PHRASES:
        escaped = re.escape(phrase)
        # Multi-word phrases: allow flexible whitespace.
        escaped = escaped.replace(r"\ ", r"\s+")
        patterns.append(re.compile(rf"\b{escaped}\b", re.IGNORECASE))
    return patterns


_CAUSAL_PATTERNS = _compile_causal_patterns()


def find_causal_language(text: str) -> list[str]:
    """Return distinct causal phrases found in ``text`` (lowercase canonical)."""
    if not text:
        return []
    found: list[str] = []
    seen: set[str] = set()
    for phrase, pattern in zip(CAUSAL_PHRASES, _CAUSAL_PATTERNS, strict=True):
        if pattern.search(text) and phrase not in seen:
            seen.add(phrase)
            found.append(phrase)
    return found


def grade_rank(grade: EvidenceGrade | None) -> int:
    """Numeric strength for comparing evidence grades (higher = stronger)."""
    if grade is None:
        return 0
    return _GRADE_RANK.get(grade, 0)


def strongest_grade(grades: Iterable[EvidenceGrade | None]) -> EvidenceGrade | None:
    """Return the strongest (best) grade among ``grades``, or None if empty."""
    best: EvidenceGrade | None = None
    best_rank = 0
    for grade in grades:
        if grade is None:
            continue
        rank = grade_rank(grade)
        if rank > best_rank:
            best = grade
            best_rank = rank
    return best


def index_evidence_by_source_id(
    evidence_records: Iterable[EvidenceRecord],
) -> dict[str, EvidenceRecord]:
    """Map ``source_id`` -> first matching evidence record."""
    out: dict[str, EvidenceRecord] = {}
    for record in evidence_records:
        sid = record.source_id
        if sid and sid not in out:
            out[sid] = record
    return out


def resolve_claim_grade(
    claim: Claim,
    cited_records: list[EvidenceRecord],
) -> EvidenceGrade | None:
    """Prefer claim.evidence_grade; else strongest grade among cited evidence."""
    if claim.evidence_grade is not None:
        return claim.evidence_grade
    return strongest_grade(r.evidence_grade for r in cited_records)


def has_disease_evidence(cited_records: list[EvidenceRecord]) -> bool:
    """True if any cited record is a disease/variant association assertion."""
    return any(r.assertion_type in _DISEASE_ASSERTIONS for r in cited_records)


def verify_claim(
    claim: Claim,
    evidence_records: Iterable[EvidenceRecord] | dict[str, EvidenceRecord],
) -> VerificationResult:
    """Verify one claim against evidence records.

    ``evidence_records`` may be a list/iterable of :class:`EvidenceRecord` or a
    pre-built ``source_id`` -> record map.
    """
    if isinstance(evidence_records, dict):
        by_id = evidence_records
    else:
        by_id = index_evidence_by_source_id(evidence_records)

    source_ids = [str(s) for s in (claim.source_ids or []) if s]
    presence_ok = len(source_ids) >= 1
    missing = [sid for sid in source_ids if sid not in by_id]
    exists_ok = presence_ok and not missing

    cited = [by_id[sid] for sid in source_ids if sid in by_id]
    causal_hits = find_causal_language(claim.claim_text or "")
    grade = resolve_claim_grade(claim, cited)
    disease_ok = has_disease_evidence(cited)

    # Defaults.
    semantic: Verdict3 = "pass"
    causal_check: Verdict3 = "pass"
    strength_check: Verdict3 = "pass"
    verdict: Verdict4 = "pass"
    needs_review = False
    reasons: list[str] = []

    if not presence_ok:
        semantic = "fail"
        verdict = "fail"
        reasons.append("Claim cites no source_id values.")
    elif not exists_ok:
        semantic = "fail"
        verdict = "fail"
        reasons.append(
            "Cited source_id(s) missing from evidence records: "
            + ", ".join(missing)
        )

    if causal_hits:
        causal_check = "warning"
        reasons.append(
            "Causal language detected: " + ", ".join(sorted(set(causal_hits)))
        )
        strong_enough = grade_rank(grade) >= grade_rank(EvidenceGrade.A)

        if not disease_ok:
            # Causal language without disease/variant evidence always needs review,
            # even when the strongest cited grade is A (e.g. gene identity only).
            strength_check = "warning"
            needs_review = True
            reasons.append(
                "Causal language with no cited disease/variant association evidence."
            )
            if verdict != "fail":
                verdict = "human_review"
        elif not strong_enough:
            strength_check = "warning"
            needs_review = True
            reasons.append(
                f"Causal language with evidence grade "
                f"{grade.value if grade else 'unknown'} (below A); "
                "human review recommended."
            )
            if verdict != "fail":
                verdict = "human_review"
        else:
            # Disease/variant evidence present and grade A — soft warning only.
            if verdict == "pass":
                verdict = "warning"

    # Soft flag: any cited evidence already marked for manual review.
    if exists_ok and any(r.manual_review_required for r in cited):
        needs_review = True
        reasons.append(
            "One or more cited evidence records require manual review."
        )
        if verdict == "pass":
            verdict = "warning"

    if verdict == "pass" and not reasons:
        reasons.append("Claim cites existing evidence; no causal overclaim flags.")

    return VerificationResult(
        claim_id=claim.id,
        source_id_presence_passed=presence_ok,
        source_exists_passed=exists_ok,
        semantic_support=semantic,
        causal_language_check=causal_check,
        evidence_strength_check=strength_check,
        verdict=verdict,
        reason="; ".join(reasons),
        needs_human_review=needs_review,
    )


def verify_claims(
    claims: Iterable[Claim],
    evidence_records: Iterable[EvidenceRecord],
) -> list[VerificationResult]:
    """Verify many claims against a shared evidence set."""
    by_id = index_evidence_by_source_id(evidence_records)
    return [verify_claim(claim, by_id) for claim in claims]


__all__ = [
    "CAUSAL_PHRASES",
    "find_causal_language",
    "grade_rank",
    "strongest_grade",
    "index_evidence_by_source_id",
    "resolve_claim_grade",
    "has_disease_evidence",
    "verify_claim",
    "verify_claims",
]
