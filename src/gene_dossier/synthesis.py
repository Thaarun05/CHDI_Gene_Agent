"""Report-section synthesis from evidence records (MVP).

Uses LangChain for optional LLM section writing when an OpenAI or Anthropic
API key is configured. Without keys (or on LLM failure), falls back to
deterministic markdown built only from :class:`~gene_dossier.models.EvidenceRecord`
fields — the LLM is never treated as a source of truth.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any, Iterable, Literal

from pydantic import BaseModel, Field

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import (
    Claim,
    EvidenceGrade,
    EvidenceRecord,
    ReportSection,
)

logger = logging.getLogger(__name__)

# Canonical CHDI-style section order (IMPLEMENTATION_PLAN §7).
CHDI_REPORT_SECTIONS: tuple[str, ...] = (
    "General gene information",
    "Gene aliases and identifiers",
    "Conservation / orthologs",
    "Known structure / domains",
    "AlphaFold / PDBe / CDD",
    "Homologues",
    "Tissue and cell expression",
    "GEO perturbations",
    "Transcription factors",
    "Protein-protein interactions",
    "CTD perturbations",
    "Chemical tools",
    "eQTLs",
    "ClinVar / OMIM / Open Targets / SNPs",
    "Pathways",
    "Knockouts / model phenotypes",
    "Major labs / literature",
    "Antibodies",
    "Patents",
    "NIH/ERC grants",
    "Missing / deferred / manual sources",
    "Verification warnings",
)

# Filled later by coverage / verification layers, not by evidence synthesis.
_META_SECTIONS = frozenset(
    {
        "Missing / deferred / manual sources",
        "Verification warnings",
    }
)

SYNTHESIS_SYSTEM_PROMPT = """\
You are a cautious biomedical analyst writing one section of a Huntington's \
disease gene dossier STRICTLY from the evidence provided.

Rules:
1. Use ONLY the supplied evidence. Do not invent facts, identifiers, or source_ids.
2. Every factual sentence must cite one or more provided source_id values using \
the exact form [source_id].
3. Do not upgrade associations, GWAS hits, expression, or computational scores \
into causal disease claims.
4. Prefer cautious wording (associated, reported, annotated) over causes / proves \
/ therapeutic target unless the cited evidence grade is A and is disease-related.
5. If evidence is thin or conflicting, say so briefly.
6. Return structured output: markdown content plus discrete claims, each with \
claim_text and the source_ids that support it.
"""

SynthesisMode = Literal["deterministic", "llm", "empty", "deferred"]


class SectionClaimDraft(BaseModel):
    """LLM / draft claim before becoming a :class:`Claim`."""

    claim_text: str
    source_ids: list[str] = Field(default_factory=list)
    evidence_grade: EvidenceGrade | None = None


class SectionDraft(BaseModel):
    """Structured LLM output for a single report section."""

    content_markdown: str = ""
    claims: list[SectionClaimDraft] = Field(default_factory=list)


class SynthesisResult(BaseModel):
    """Bundle of report sections and claims for one dossier run."""

    dossier_run_id: str
    gene_symbol: str
    sections: list[ReportSection] = Field(default_factory=list)
    claims: list[Claim] = Field(default_factory=list)
    mode: SynthesisMode = "deterministic"
    notes: list[str] = Field(default_factory=list)


# --------------------------------------------------------------------------------------
# Evidence formatting / grouping
# --------------------------------------------------------------------------------------
def group_evidence_by_section(
    evidence_records: Iterable[EvidenceRecord],
) -> dict[str, list[EvidenceRecord]]:
    """Group records by ``section`` (insertion order preserved within groups)."""
    grouped: dict[str, list[EvidenceRecord]] = defaultdict(list)
    for record in evidence_records:
        section = (record.section or "").strip() or "Unsectioned"
        grouped[section].append(record)
    return dict(grouped)


def format_evidence_block(records: list[EvidenceRecord], *, max_records: int = 80) -> str:
    """Serialize evidence for prompts / deterministic rendering."""
    lines: list[str] = []
    for record in records[:max_records]:
        grade = record.evidence_grade.value if record.evidence_grade else "?"
        sub = f" / {record.subsection}" if record.subsection else ""
        review = " [manual_review]" if record.manual_review_required else ""
        text = (record.display_text or "").strip() or "(no display_text)"
        lines.append(
            f"- source_id={record.source_id} | grade={grade} | "
            f"{record.source_name} | {record.assertion_type.value}{sub}{review}\n"
            f"  {text}"
        )
    if len(records) > max_records:
        lines.append(f"- … {len(records) - max_records} additional evidence records omitted.")
    return "\n".join(lines)


def _allowed_source_ids(records: list[EvidenceRecord]) -> set[str]:
    return {r.source_id for r in records if r.source_id}


def _filter_claim_source_ids(
    source_ids: list[str],
    allowed: set[str],
) -> list[str]:
    """Keep only source_ids present in the section evidence (deduped, order preserved)."""
    out: list[str] = []
    seen: set[str] = set()
    for sid in source_ids:
        s = str(sid).strip()
        if not s or s not in allowed or s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def _strongest_grade_for_ids(
    source_ids: list[str],
    by_id: dict[str, EvidenceRecord],
) -> EvidenceGrade | None:
    best: EvidenceGrade | None = None
    best_rank = -1
    rank = {"A": 6, "B": 5, "C": 4, "D": 3, "E": 2, "F": 1}
    for sid in source_ids:
        rec = by_id.get(sid)
        if rec is None or rec.evidence_grade is None:
            continue
        r = rank.get(rec.evidence_grade.value, 0)
        if r > best_rank:
            best_rank = r
            best = rec.evidence_grade
    return best


# --------------------------------------------------------------------------------------
# Deterministic synthesis
# --------------------------------------------------------------------------------------
def synthesize_section_deterministic(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    section_name: str,
    records: list[EvidenceRecord],
) -> tuple[ReportSection, list[Claim]]:
    """Build a section + one claim per evidence row from display_text only."""
    if section_name in _META_SECTIONS:
        section = ReportSection(
            dossier_run_id=dossier_run_id,
            section_name=section_name,
            content_markdown=(
                f"_Deferred: `{section_name}` is populated by coverage / "
                "verification layers, not evidence synthesis._\n"
            ),
            source_ids=[],
            status="deferred",
        )
        return section, []

    if not records:
        section = ReportSection(
            dossier_run_id=dossier_run_id,
            section_name=section_name,
            content_markdown=(
                f"**{section_name}**\n\n"
                f"_No evidence records available for {gene_symbol} in this section._\n"
            ),
            source_ids=[],
            status="empty",
        )
        return section, []

    by_id = {r.source_id: r for r in records if r.source_id}
    bullets: list[str] = []
    claims: list[Claim] = []
    source_ids: list[str] = []
    seen_sid: set[str] = set()

    for record in records:
        sid = record.source_id
        if sid and sid not in seen_sid:
            seen_sid.add(sid)
            source_ids.append(sid)
        text = (record.display_text or "").strip() or "(no display_text)"
        grade = record.evidence_grade.value if record.evidence_grade else "?"
        cite = f" [{sid}]" if sid else ""
        bullets.append(f"- ({grade}) {text}{cite}")
        if sid:
            claims.append(
                Claim(
                    dossier_run_id=dossier_run_id,
                    claim_text=text,
                    source_ids=[sid],
                    evidence_grade=record.evidence_grade,
                    claim_type=record.assertion_type.value,
                )
            )

    body = "\n".join(bullets)
    markdown = (
        f"**{section_name}** ({gene_symbol})\n\n"
        f"{body}\n\n"
        f"_Synthesized deterministically from {len(records)} evidence record(s)._\n"
    )
    section = ReportSection(
        dossier_run_id=dossier_run_id,
        section_name=section_name,
        content_markdown=markdown,
        source_ids=source_ids,
        status="deterministic",
    )
    for claim in claims:
        claim.section_id = section.id
        if claim.evidence_grade is None:
            claim.evidence_grade = _strongest_grade_for_ids(claim.source_ids, by_id)
    return section, claims


# --------------------------------------------------------------------------------------
# LLM synthesis (LangChain)
# --------------------------------------------------------------------------------------
def build_chat_model(settings: Settings | None = None) -> Any | None:
    """Return a LangChain chat model if an LLM API key is configured, else None.

    Provider import/construction failures are logged and return ``None`` so
    callers fall back to deterministic synthesis instead of crashing.
    """
    cfg = settings or get_settings()
    if not cfg.has_llm():
        return None

    if cfg.has_key("ANTHROPIC_API_KEY"):
        try:
            from langchain_anthropic import ChatAnthropic

            model_name = (cfg.default_llm_model or "claude-sonnet-4-20250514").strip()
            return ChatAnthropic(
                model=model_name,
                api_key=cfg.anthropic_api_key,
                temperature=0,
            )
        except Exception as exc:  # noqa: BLE001 — soft-fail to deterministic
            logger.warning("Failed to build Anthropic chat model: %s", exc)
            return None

    if cfg.has_key("OPENAI_API_KEY"):
        try:
            from langchain_openai import ChatOpenAI

            model_name = (cfg.default_llm_model or "gpt-4o-mini").strip()
            return ChatOpenAI(
                model=model_name,
                api_key=cfg.openai_api_key,
                temperature=0,
            )
        except Exception as exc:  # noqa: BLE001 — soft-fail to deterministic
            logger.warning("Failed to build OpenAI chat model: %s", exc)
            return None

    return None


def _invoke_section_llm(
    *,
    model: Any,
    gene_symbol: str,
    section_name: str,
    records: list[EvidenceRecord],
) -> SectionDraft:
    """Call LangChain structured output for one section."""
    from langchain_core.prompts import ChatPromptTemplate

    evidence_block = format_evidence_block(records)
    prompt = ChatPromptTemplate.from_messages(
        [
            ("system", SYNTHESIS_SYSTEM_PROMPT),
            (
                "human",
                "Gene: {gene_symbol}\n"
                "Section: {section_name}\n\n"
                "Evidence (use only these source_ids):\n{evidence_block}\n\n"
                "Write the section markdown and extract discrete claims.",
            ),
        ]
    )
    chain = prompt | model.with_structured_output(SectionDraft)
    result = chain.invoke(
        {
            "gene_symbol": gene_symbol,
            "section_name": section_name,
            "evidence_block": evidence_block,
        }
    )
    if isinstance(result, SectionDraft):
        return result
    if isinstance(result, dict):
        return SectionDraft.model_validate(result)
    raise TypeError(f"Unexpected structured output type: {type(result)!r}")


def synthesize_section_llm(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    section_name: str,
    records: list[EvidenceRecord],
    model: Any,
) -> tuple[ReportSection, list[Claim]]:
    """LLM section write; raises on failure (caller may fall back)."""
    if section_name in _META_SECTIONS or not records:
        return synthesize_section_deterministic(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            section_name=section_name,
            records=records,
        )

    draft = _invoke_section_llm(
        model=model,
        gene_symbol=gene_symbol,
        section_name=section_name,
        records=records,
    )
    allowed = _allowed_source_ids(records)
    by_id = {r.source_id: r for r in records if r.source_id}

    content = (draft.content_markdown or "").strip()
    if not content:
        raise ValueError("LLM returned empty section markdown")

    section_source_ids: list[str] = []
    seen: set[str] = set()
    claims: list[Claim] = []

    for draft_claim in draft.claims:
        sids = _filter_claim_source_ids(draft_claim.source_ids, allowed)
        if not sids:
            continue
        text = (draft_claim.claim_text or "").strip()
        if not text:
            continue
        for sid in sids:
            if sid not in seen:
                seen.add(sid)
                section_source_ids.append(sid)
        grade = draft_claim.evidence_grade or _strongest_grade_for_ids(sids, by_id)
        claims.append(
            Claim(
                dossier_run_id=dossier_run_id,
                claim_text=text,
                source_ids=sids,
                evidence_grade=grade,
                claim_type="llm_section",
            )
        )

    # If the model wrote prose but no usable claims, fall back to evidence claims.
    if not claims:
        _, det_claims = synthesize_section_deterministic(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            section_name=section_name,
            records=records,
        )
        claims = det_claims
        section_source_ids = list(
            dict.fromkeys(sid for c in claims for sid in c.source_ids)
        )

    section = ReportSection(
        dossier_run_id=dossier_run_id,
        section_name=section_name,
        content_markdown=content + "\n",
        source_ids=section_source_ids,
        status="llm",
    )
    for claim in claims:
        claim.section_id = section.id
    return section, claims


# --------------------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------------------
def synthesize_section(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    section_name: str,
    records: list[EvidenceRecord],
    settings: Settings | None = None,
    force_deterministic: bool = False,
    model: Any | None = None,
) -> tuple[ReportSection, list[Claim], SynthesisMode]:
    """Synthesize one section; prefer LLM when available, else deterministic."""
    cfg = settings or get_settings()
    if (
        not force_deterministic
        and section_name not in _META_SECTIONS
        and records
    ):
        chat = model if model is not None else build_chat_model(cfg)
        if chat is not None:
            try:
                section, claims = synthesize_section_llm(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    section_name=section_name,
                    records=records,
                    model=chat,
                )
                return section, claims, "llm"
            except Exception as exc:  # noqa: BLE001 — soft-fail to deterministic
                logger.warning(
                    "LLM synthesis failed for section %r; using deterministic fallback: %s",
                    section_name,
                    exc,
                )

    section, claims = synthesize_section_deterministic(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        section_name=section_name,
        records=records,
    )
    mode: SynthesisMode
    if section.status == "empty":
        mode = "empty"
    elif section.status == "deferred":
        mode = "deferred"
    else:
        mode = "deterministic"
    return section, claims, mode


def synthesize_dossier(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    evidence_records: Iterable[EvidenceRecord],
    settings: Settings | None = None,
    force_deterministic: bool = False,
    section_names: Iterable[str] | None = None,
) -> SynthesisResult:
    """Synthesize all CHDI sections (or a provided subset) for one gene run.

    Without LLM keys, every content section is deterministic. Meta sections are
    deferred placeholders. Unknown evidence sections (not in the CHDI list) are
    appended after the canonical list so nothing is silently dropped.
    """
    cfg = settings or get_settings()
    records = list(evidence_records)
    grouped = group_evidence_by_section(records)

    canonical = list(section_names) if section_names is not None else list(CHDI_REPORT_SECTIONS)
    extras = [name for name in grouped if name not in canonical and name != "Unsectioned"]
    ordered = canonical + extras
    if "Unsectioned" in grouped and "Unsectioned" not in ordered:
        ordered.append("Unsectioned")

    chat = None if force_deterministic else build_chat_model(cfg)
    notes: list[str] = []
    if chat is None and not force_deterministic and not cfg.has_llm():
        notes.append("No LLM API key configured; using deterministic synthesis.")
    elif force_deterministic:
        notes.append("force_deterministic=True; skipped LLM synthesis.")

    sections: list[ReportSection] = []
    claims: list[Claim] = []
    used_llm = False
    used_det = False

    for name in ordered:
        section_records = grouped.get(name, [])
        section, section_claims, mode = synthesize_section(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            section_name=name,
            records=section_records,
            settings=cfg,
            force_deterministic=force_deterministic,
            model=chat,
        )
        sections.append(section)
        claims.extend(section_claims)
        if mode == "llm":
            used_llm = True
        elif mode == "deterministic":
            used_det = True

    if used_llm:
        overall_mode: SynthesisMode = "llm"
        if used_det:
            notes.append("Some sections used LLM; others used deterministic fallback.")
    else:
        overall_mode = "deterministic"

    return SynthesisResult(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        sections=sections,
        claims=claims,
        mode=overall_mode,
        notes=notes,
    )


__all__ = [
    "CHDI_REPORT_SECTIONS",
    "SYNTHESIS_SYSTEM_PROMPT",
    "SectionClaimDraft",
    "SectionDraft",
    "SynthesisResult",
    "group_evidence_by_section",
    "format_evidence_block",
    "build_chat_model",
    "synthesize_section_deterministic",
    "synthesize_section_llm",
    "synthesize_section",
    "synthesize_dossier",
]
