"""Report-section synthesis from evidence records (MVP).

Uses LangChain for optional LLM section writing when an OpenAI, NVIDIA NIM,
Anthropic, or Google Gemini API key is configured. Without keys (or on LLM failure), falls back to
deterministic markdown built only from :class:`~gene_dossier.models.EvidenceRecord`
fields — the LLM is never treated as a source of truth. Multiple providers are
tried in order at invocation time so one failing vendor does not block synthesis.
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
You are a Rancho BioSciences / CHDI-style gene-dossier writer preparing one \
section of a Huntington's disease research gene report.

Tone and style:
- Concise, professional biomedical prose suitable for a CHDI/Rancho dossier.
- Prefer measured language (associated, reported, annotated, curated) over \
causal claims.
- Do not write marketing copy, speculation, or therapeutic recommendations.

Evidence rules (non-negotiable):
1. Use ONLY the supplied evidence records. Do not invent biology, identifiers, \
assays, pathways, diseases, labs, patents, or source_ids.
2. Do not claim facts outside the supplied evidence.
3. Do not upgrade associations, GWAS hits, expression patterns, or computational \
scores into causal or clinically validated disease statements.
4. Prefer "associated / reported / annotated" over "causes / proves / drives / \
therapeutic target / pathogenic / disease-modifying" unless the cited evidence \
grade is A and explicitly disease-related.
5. If evidence for this section is thin, conflicting, or absent, say that \
evidence was unavailable or limited from this run. Do not invent filler.

Output (structured):
- section_id: short slug for the section (may echo the section name).
- subsection_id: optional subsection slug, else null/omit.
- summary_paragraphs: 1–3 short paragraphs of dossier prose; cite evidence with \
exact [source_id] tokens when stating facts.
- key_findings: short bullet-ready findings grounded in the evidence.
- claims: discrete factual claims. Prefer supporting_evidence_ids (list of \
exact source_id values from the evidence). You may also set source_ids; when \
supporting_evidence_ids is present it takes precedence. Every claim MUST cite \
at least one valid supplied source_id. Omit claims you cannot support.
- limitations: brief caveats (coverage gaps, conflicting evidence, grade limits).
- content_markdown: optional; leave empty unless you must supply preformatted \
markdown. The system will render structured fields into markdown when this is \
blank.
"""

SynthesisMode = Literal["deterministic", "llm", "empty", "deferred"]
LlmPurpose = Literal["section", "planner", "answer"]


class SectionClaimDraft(BaseModel):
    """LLM / draft claim before becoming a :class:`Claim`."""

    claim_text: str
    source_ids: list[str] = Field(default_factory=list)
    supporting_evidence_ids: list[str] = Field(default_factory=list)
    evidence_grade: EvidenceGrade | None = None


class SectionDraft(BaseModel):
    """Structured LLM output for a single Rancho/CHDI report section."""

    section_id: str = ""
    subsection_id: str | None = None
    summary_paragraphs: list[str] = Field(default_factory=list)
    key_findings: list[str] = Field(default_factory=list)
    claims: list[SectionClaimDraft] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    content_markdown: str = ""


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


def _claim_draft_source_ids(draft_claim: SectionClaimDraft) -> list[str]:
    """Prefer supporting_evidence_ids when present; else fall back to source_ids."""
    supporting = [
        str(sid).strip()
        for sid in (draft_claim.supporting_evidence_ids or [])
        if str(sid).strip()
    ]
    if supporting:
        return supporting
    return [
        str(sid).strip()
        for sid in (draft_claim.source_ids or [])
        if str(sid).strip()
    ]


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


def _render_section_draft_markdown(
    draft: SectionDraft,
    *,
    gene_symbol: str,
    section_name: str,
) -> str:
    """Render structured Rancho draft fields into ReportSection markdown."""
    preformatted = (draft.content_markdown or "").strip()
    if preformatted:
        return preformatted

    lines: list[str] = [f"**{section_name}** ({gene_symbol})", ""]
    paragraphs = [p.strip() for p in (draft.summary_paragraphs or []) if p and str(p).strip()]
    if paragraphs:
        lines.extend(paragraphs)
        lines.append("")

    findings = [f.strip() for f in (draft.key_findings or []) if f and str(f).strip()]
    if findings:
        lines.append("**Key findings**")
        lines.append("")
        for item in findings:
            lines.append(f"- {item}")
        lines.append("")

    limitations = [lim.strip() for lim in (draft.limitations or []) if lim and str(lim).strip()]
    if limitations:
        lines.append("**Limitations**")
        lines.append("")
        for item in limitations:
            lines.append(f"- {item}")
        lines.append("")

    body = "\n".join(lines).strip()
    return body


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
                f"_Evidence for {gene_symbol} in this section was unavailable "
                "from this run._\n"
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
DEFAULT_NIM_BASE_URL = "https://integrate.api.nvidia.com/v1"
DEFAULT_NIM_MODEL = "meta/llama-3.1-8b-instruct"
DEFAULT_OPENAI_MODEL = "gpt-4o-mini"
DEFAULT_ANTHROPIC_MODEL = "claude-sonnet-4-20250514"
DEFAULT_GOOGLE_GEMINI_MODEL = "gemini-3.5-flash"

PROVIDER_OPENAI = "openai"
PROVIDER_NVIDIA_NIM = "nvidia_nim"
PROVIDER_ANTHROPIC = "anthropic"
PROVIDER_GOOGLE_GEMINI = "google_gemini"


class LlmModelCandidate(BaseModel):
    """One constructed chat model plus its provider label (for fallback logging)."""

    provider: str
    model: Any


def _normalized_provider(settings: Settings) -> str:
    raw = (settings.default_llm_provider or "").strip().lower()
    if raw in {
        PROVIDER_OPENAI,
        PROVIDER_NVIDIA_NIM,
        PROVIDER_ANTHROPIC,
        PROVIDER_GOOGLE_GEMINI,
    }:
        return raw
    return ""


def _provider_order(settings: Settings) -> list[str]:
    """Ordered provider names for construction / invocation fallback."""
    preferred = _normalized_provider(settings)
    auto = [
        PROVIDER_OPENAI,
        PROVIDER_NVIDIA_NIM,
        PROVIDER_ANTHROPIC,
        PROVIDER_GOOGLE_GEMINI,
    ]
    if not preferred:
        return auto
    rest = [p for p in auto if p != preferred]
    return [preferred, *rest]


def _sanitize_llm_error(exc: BaseException, *, limit: int = 200) -> str:
    """Short error text without dumping secrets or full env."""
    text = f"{type(exc).__name__}: {exc}"
    collapsed = " ".join(text.split())
    if len(collapsed) <= limit:
        return collapsed
    return collapsed[: limit - 3] + "..."


def _build_openai_chat(settings: Settings, *, purpose: LlmPurpose = "section") -> Any | None:
    if not settings.has_key("openai_api_key"):
        return None
    try:
        from langchain_openai import ChatOpenAI

        if purpose in {"planner", "answer"}:
            model_name = (settings.openai_model or "").strip() or "gpt-5.6-terra"
            reasoning_effort = (
                settings.openai_planner_reasoning_effort
                if purpose == "planner"
                else settings.openai_answer_reasoning_effort
            )
            max_tokens = (
                settings.openai_planner_max_output_tokens
                if purpose == "planner"
                else settings.openai_answer_max_output_tokens
            )
            timeout = settings.openai_timeout_seconds
            use_responses_api = True
        else:
            model_name = (settings.default_llm_model or DEFAULT_OPENAI_MODEL).strip()
            reasoning_effort = None
            max_tokens = None
            timeout = settings.http_timeout_seconds
            use_responses_api = None
        api_key = settings.openai_api_key
        if hasattr(api_key, "get_secret_value"):
            api_key = api_key.get_secret_value()
        kwargs: dict[str, Any] = {
            "model": model_name,
            "api_key": api_key,
            "timeout": timeout,
            "max_retries": 0,
        }
        if purpose == "section":
            kwargs["temperature"] = 0
        else:
            kwargs.update(
                {
                    "reasoning_effort": reasoning_effort,
                    "max_completion_tokens": max_tokens,
                    "use_responses_api": use_responses_api,
                    "store": False,
                    "service_tier": "default",
                }
            )
        base_url = (settings.openai_base_url or "").strip()
        if base_url:
            kwargs["base_url"] = base_url
        return ChatOpenAI(**kwargs)
    except Exception as exc:  # noqa: BLE001 — soft-fail to next provider
        logger.warning(
            "Failed to build OpenAI chat model: %s", _sanitize_llm_error(exc)
        )
        return None


def _build_nvidia_nim_chat(settings: Settings) -> Any | None:
    if not settings.has_key("nvidia_nim_api_key"):
        return None
    try:
        from langchain_openai import ChatOpenAI

        base_url = (settings.nvidia_nim_base_url or "").strip() or DEFAULT_NIM_BASE_URL
        model_name = (
            (settings.nvidia_nim_model or "").strip()
            or (settings.default_llm_model or "").strip()
            or DEFAULT_NIM_MODEL
        )
        api_key = settings.nvidia_nim_api_key
        if hasattr(api_key, "get_secret_value"):
            api_key = api_key.get_secret_value()
        return ChatOpenAI(
            model=model_name,
            api_key=api_key,
            base_url=base_url,
            temperature=0,
            timeout=settings.nvidia_nim_timeout_seconds,
            max_retries=0,
        )
    except Exception as exc:  # noqa: BLE001 — soft-fail to next provider
        logger.warning(
            "Failed to build NVIDIA NIM chat model: %s", _sanitize_llm_error(exc)
        )
        return None


def _build_anthropic_chat(settings: Settings) -> Any | None:
    if not settings.has_key("anthropic_api_key"):
        return None
    try:
        from langchain_anthropic import ChatAnthropic

        model_name = (settings.default_llm_model or DEFAULT_ANTHROPIC_MODEL).strip()
        return ChatAnthropic(
            model=model_name,
            api_key=settings.anthropic_api_key,
            temperature=0,
            timeout=settings.http_timeout_seconds,
            max_retries=0,
        )
    except Exception as exc:  # noqa: BLE001 — soft-fail to next provider
        logger.warning(
            "Failed to build Anthropic chat model: %s", _sanitize_llm_error(exc)
        )
        return None


def _build_google_gemini_chat(settings: Settings) -> Any | None:
    if not settings.has_key("google_api_key"):
        return None
    try:
        from langchain_google_genai import ChatGoogleGenerativeAI

        model_name = (settings.google_gemini_model or "").strip() or DEFAULT_GOOGLE_GEMINI_MODEL
        return ChatGoogleGenerativeAI(
            model=model_name,
            google_api_key=settings.google_api_key,
            thinking_level="low",
            timeout=settings.http_timeout_seconds,
            max_retries=0,
        )
    except Exception as exc:  # noqa: BLE001 — soft-fail to next provider
        logger.warning(
            "Failed to build Google Gemini chat model: %s", _sanitize_llm_error(exc)
        )
        return None


def build_chat_model_candidates(
    settings: Settings | None = None,
    *,
    purpose: LlmPurpose = "section",
) -> list[LlmModelCandidate]:
    """Build ordered LLM candidates for invocation-time fallback.

    Skips providers with missing keys. Construction failures are logged (no
    secrets) and the next provider is tried. Never raises.
    """
    cfg = settings or get_settings()
    builders = {
        PROVIDER_OPENAI: _build_openai_chat,
        PROVIDER_NVIDIA_NIM: _build_nvidia_nim_chat,
        PROVIDER_ANTHROPIC: _build_anthropic_chat,
        PROVIDER_GOOGLE_GEMINI: _build_google_gemini_chat,
    }
    out: list[LlmModelCandidate] = []
    for provider in _provider_order(cfg):
        builder = builders.get(provider)
        if builder is None:
            continue
        if provider == PROVIDER_OPENAI:
            model = builder(cfg, purpose=purpose)
        else:
            model = builder(cfg)
        if model is None:
            continue
        out.append(LlmModelCandidate(provider=provider, model=model))
    return out


def build_chat_model(settings: Settings | None = None) -> Any | None:
    """Return the first available chat model, or None.

    Prefer :func:`build_chat_model_candidates` when invocation-time fallback
    across providers is needed. Supports OpenAI, NVIDIA NIM (OpenAI-compatible),
    Anthropic, and Google Gemini.
    """
    candidates = build_chat_model_candidates(settings)
    if not candidates:
        return None
    return candidates[0].model


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
                "Evidence (use only these source_ids; do not invent others):\n"
                "{evidence_block}\n\n"
                "Return structured Rancho section fields: summary_paragraphs, "
                "key_findings, claims (with supporting_evidence_ids), and "
                "limitations. Leave content_markdown empty unless necessary.",
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

    content = _render_section_draft_markdown(
        draft, gene_symbol=gene_symbol, section_name=section_name
    ).strip()
    if not content:
        raise ValueError("LLM returned empty section markdown")

    section_source_ids: list[str] = []
    seen: set[str] = set()
    claims: list[Claim] = []

    for draft_claim in draft.claims:
        text = (draft_claim.claim_text or "").strip()
        if not text:
            continue
        sids = _filter_claim_source_ids(_claim_draft_source_ids(draft_claim), allowed)
        if not sids:
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
        subsection_name=draft.subsection_id or None,
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
    model_candidates: list[LlmModelCandidate] | None = None,
) -> tuple[ReportSection, list[Claim], SynthesisMode]:
    """Synthesize one section; try LLM providers in order, else deterministic.

    When ``model`` is provided (tests / injected client), it is tried alone.
    Otherwise ``model_candidates`` or :func:`build_chat_model_candidates` are
    tried at **invocation** time so one failing provider does not block others.
    """
    cfg = settings or get_settings()
    if (
        not force_deterministic
        and section_name not in _META_SECTIONS
        and records
    ):
        if model is not None:
            candidates = [LlmModelCandidate(provider="injected", model=model)]
        elif model_candidates is not None:
            candidates = list(model_candidates)
        else:
            candidates = build_chat_model_candidates(cfg)

        for candidate in candidates:
            try:
                section, claims = synthesize_section_llm(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    section_name=section_name,
                    records=records,
                    model=candidate.model,
                )
                return section, claims, "llm"
            except Exception as exc:  # noqa: BLE001 — try next provider
                logger.warning(
                    "LLM provider %s failed for section %r; trying next: %s",
                    candidate.provider,
                    section_name,
                    _sanitize_llm_error(exc),
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

    Configured LLM providers are tried in order per section; invocation failures
    fall through to the next provider before deterministic fallback.
    """
    cfg = settings or get_settings()
    records = list(evidence_records)
    grouped = group_evidence_by_section(records)

    canonical = list(section_names) if section_names is not None else list(CHDI_REPORT_SECTIONS)
    extras = [name for name in grouped if name not in canonical and name != "Unsectioned"]
    ordered = canonical + extras
    if "Unsectioned" in grouped and "Unsectioned" not in ordered:
        ordered.append("Unsectioned")

    candidates: list[LlmModelCandidate] = (
        [] if force_deterministic else build_chat_model_candidates(cfg)
    )
    notes: list[str] = []
    if not candidates and not force_deterministic and not cfg.has_llm():
        notes.append("No LLM API key configured; using deterministic synthesis.")
    elif force_deterministic:
        notes.append("force_deterministic=True; skipped LLM synthesis.")
    elif not candidates and cfg.has_llm():
        notes.append(
            "LLM keys present but no chat model could be constructed; "
            "using deterministic synthesis."
        )

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
            model_candidates=candidates,
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
    "DEFAULT_NIM_BASE_URL",
    "DEFAULT_GOOGLE_GEMINI_MODEL",
    "PROVIDER_GOOGLE_GEMINI",
    "LlmModelCandidate",
    "SectionClaimDraft",
    "SectionDraft",
    "SynthesisResult",
    "group_evidence_by_section",
    "format_evidence_block",
    "build_chat_model",
    "build_chat_model_candidates",
    "synthesize_section_deterministic",
    "synthesize_section_llm",
    "synthesize_section",
    "synthesize_dossier",
]
