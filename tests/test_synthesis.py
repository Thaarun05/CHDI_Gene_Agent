"""Tests for deterministic and LLM-assisted Rancho section synthesis."""

from __future__ import annotations

from gene_dossier.config import Settings
from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
)
from gene_dossier.synthesis import (
    DEFAULT_NIM_BASE_URL,
    LlmModelCandidate,
    SYNTHESIS_SYSTEM_PROMPT,
    SectionClaimDraft,
    SectionDraft,
    build_chat_model_candidates,
    format_evidence_block,
    synthesize_section,
    synthesize_section_deterministic,
    synthesize_section_llm,
)


def _llm_settings(**overrides) -> Settings:
    """Settings isolated from ``.env`` for provider / candidate tests."""
    base = {
        "openai_api_key": None,
        "openai_base_url": None,
        "anthropic_api_key": None,
        "nvidia_nim_api_key": None,
        "nvidia_nim_base_url": None,
        "nvidia_nim_model": None,
        "default_llm_model": None,
        "default_llm_provider": None,
    }
    base.update(overrides)
    return Settings(_env_file=None, **base)


class _TaggedModel:
    """Marker object so invocation fakes can distinguish providers."""

    def __init__(self, tag: str) -> None:
        self.tag = tag


def _evidence(
    *,
    source_id: str,
    section: str,
    display_text: str,
    source_name: str = "NCBI Gene",
) -> EvidenceRecord:
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section=section,
        source_name=source_name,
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="test_fact",
        evidence_grade=EvidenceGrade.C,
        display_text=display_text,
        value={},
    )


def test_deterministic_fallback_still_works():
    records = [
        _evidence(
            source_id="sid-ncbi-1",
            section="General gene information",
            display_text="SREBF2 Entrez Gene ID is 6721.",
        )
    ]
    section, claims, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        force_deterministic=True,
        model=object(),  # must be ignored when force_deterministic=True
    )
    assert mode == "deterministic"
    assert section.status == "deterministic"
    assert "6721" in section.content_markdown
    assert "sid-ncbi-1" in section.content_markdown
    assert len(claims) == 1
    assert claims[0].source_ids == ["sid-ncbi-1"]
    assert "Rancho BioSciences" in SYNTHESIS_SYSTEM_PROMPT


def test_prompt_evidence_block_contains_only_relevant_records(monkeypatch):
    relevant = _evidence(
        source_id="sid-path-1",
        section="Pathways",
        display_text="SREBF2 participates in cholesterol biosynthesis.",
        source_name="Reactome",
    )
    irrelevant = _evidence(
        source_id="sid-ppi-99",
        section="Protein-protein interactions",
        display_text="SREBF2 interacts with FOXO1.",
        source_name="STRING",
    )
    captured: dict = {}

    def fake_invoke(*, model, gene_symbol, section_name, records):
        captured["records"] = list(records)
        captured["evidence_block"] = format_evidence_block(records)
        return SectionDraft(
            section_id="pathways",
            summary_paragraphs=["SREBF2 is annotated in cholesterol pathways [sid-path-1]."],
            key_findings=["Cholesterol biosynthesis pathway membership [sid-path-1]."],
            claims=[
                SectionClaimDraft(
                    claim_text="SREBF2 is in cholesterol biosynthesis.",
                    supporting_evidence_ids=["sid-path-1"],
                )
            ],
            limitations=["Limited to curated pathway membership in this run."],
        )

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    section, claims = synthesize_section_llm(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="Pathways",
        records=[relevant],
        model=object(),
    )
    assert [r.source_id for r in captured["records"]] == ["sid-path-1"]
    assert "sid-path-1" in captured["evidence_block"]
    assert "sid-ppi-99" not in captured["evidence_block"]
    assert irrelevant.source_id not in section.content_markdown
    assert section.status == "llm"
    assert claims[0].source_ids == ["sid-path-1"]


def test_empty_section_unavailable_language_no_claims():
    section, claims = synthesize_section_deterministic(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="Antibodies",
        records=[],
    )
    assert section.status == "empty"
    assert claims == []
    assert "unavailable from this run" in section.content_markdown.lower()

    section2, claims2, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="Antibodies",
        records=[],
        force_deterministic=False,
        model=object(),
    )
    assert mode == "empty"
    assert claims2 == []
    assert "unavailable from this run" in section2.content_markdown.lower()


def test_supporting_evidence_ids_become_claim_source_ids(monkeypatch):
    records = [
        _evidence(
            source_id="sid-a",
            section="General gene information",
            display_text="Official symbol SREBF2.",
        )
    ]

    def fake_invoke(*, model, gene_symbol, section_name, records):
        return SectionDraft(
            section_id="general",
            summary_paragraphs=["SREBF2 is the official symbol [sid-a]."],
            claims=[
                SectionClaimDraft(
                    claim_text="Official symbol is SREBF2.",
                    source_ids=["should-be-ignored"],
                    supporting_evidence_ids=["sid-a"],
                )
            ],
        )

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    section, claims = synthesize_section_llm(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        model=object(),
    )
    assert section.status == "llm"
    assert len(claims) == 1
    assert claims[0].source_ids == ["sid-a"]
    assert "should-be-ignored" not in claims[0].source_ids


def test_claim_without_ids_is_omitted(monkeypatch):
    records = [
        _evidence(
            source_id="sid-a",
            section="General gene information",
            display_text="Official symbol SREBF2.",
        )
    ]

    def fake_invoke(*, model, gene_symbol, section_name, records):
        return SectionDraft(
            section_id="general",
            summary_paragraphs=["SREBF2 summary grounded in evidence [sid-a]."],
            key_findings=["Symbol annotated [sid-a]."],
            claims=[
                SectionClaimDraft(
                    claim_text="Unsupported claim with no IDs.",
                    source_ids=[],
                    supporting_evidence_ids=[],
                ),
                SectionClaimDraft(
                    claim_text="Official symbol is SREBF2.",
                    supporting_evidence_ids=["sid-a"],
                ),
            ],
        )

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    _, claims = synthesize_section_llm(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        model=object(),
    )
    assert len(claims) == 1
    assert claims[0].claim_text == "Official symbol is SREBF2."
    assert claims[0].source_ids == ["sid-a"]


def test_claim_with_unknown_source_id_is_omitted(monkeypatch):
    records = [
        _evidence(
            source_id="sid-a",
            section="General gene information",
            display_text="Official symbol SREBF2.",
        )
    ]

    def fake_invoke(*, model, gene_symbol, section_name, records):
        return SectionDraft(
            section_id="general",
            summary_paragraphs=["SREBF2 summary [sid-a]."],
            claims=[
                SectionClaimDraft(
                    claim_text="Invented claim.",
                    supporting_evidence_ids=["sid-unknown"],
                ),
                SectionClaimDraft(
                    claim_text="",
                    supporting_evidence_ids=["sid-a"],
                ),
            ],
        )

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    section, claims = synthesize_section_llm(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        model=object(),
    )
    # LLM claims rejected → reuse deterministic evidence-backed claims.
    assert section.status == "llm"
    assert "SREBF2 summary" in section.content_markdown
    assert len(claims) == 1
    assert claims[0].source_ids == ["sid-a"]
    assert claims[0].claim_text == "Official symbol SREBF2."


def test_fake_structured_llm_draft_renders_markdown_status_llm(monkeypatch):
    records = [
        _evidence(
            source_id="sid-expr-1",
            section="Tissue and cell expression",
            display_text="SREBF2 is expressed in brain cortex.",
            source_name="GTEx",
        )
    ]

    def fake_invoke(*, model, gene_symbol, section_name, records):
        return SectionDraft(
            section_id="expression",
            subsection_id="tissue",
            summary_paragraphs=[
                "GTEx annotates SREBF2 expression in cortex [sid-expr-1].",
                "No causal disease inference is made from expression alone.",
            ],
            key_findings=["Cortex expression reported [sid-expr-1]."],
            claims=[
                SectionClaimDraft(
                    claim_text="SREBF2 expression is reported in cortex.",
                    supporting_evidence_ids=["sid-expr-1"],
                )
            ],
            limitations=["Expression is associative, not causal."],
            content_markdown="",
        )

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    section, claims, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="Tissue and cell expression",
        records=records,
        force_deterministic=False,
        model=object(),
    )
    assert mode == "llm"
    assert section.status == "llm"
    assert section.subsection_name == "tissue"
    assert "Key findings" in section.content_markdown
    assert "Limitations" in section.content_markdown
    assert "cortex" in section.content_markdown.lower()
    assert len(claims) == 1
    assert claims[0].source_ids == ["sid-expr-1"]


def test_meta_section_remains_deferred():
    section, claims = synthesize_section_deterministic(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="Missing / deferred / manual sources",
        records=[],
    )
    assert section.status == "deferred"
    assert claims == []
    assert "Deferred" in section.content_markdown


def test_nvidia_nim_preferred_constructs_chat_openai_with_nim_url(monkeypatch):
    captured: dict = {}

    class FakeChatOpenAI:
        def __init__(self, **kwargs):
            captured.clear()
            captured.update(kwargs)

    monkeypatch.setattr("langchain_openai.ChatOpenAI", FakeChatOpenAI)
    cfg = _llm_settings(
        default_llm_provider="nvidia_nim",
        nvidia_nim_api_key="nim-test-key",
        nvidia_nim_base_url="https://integrate.api.nvidia.com/v1",
        nvidia_nim_model="meta/llama-3.1-70b-instruct",
    )
    candidates = build_chat_model_candidates(cfg)
    assert len(candidates) == 1
    assert candidates[0].provider == "nvidia_nim"
    assert captured["api_key"] == "nim-test-key"
    assert captured["base_url"] == "https://integrate.api.nvidia.com/v1"
    assert captured["model"] == "meta/llama-3.1-70b-instruct"
    assert captured["temperature"] == 0


def test_nvidia_nim_base_url_defaults_when_omitted(monkeypatch):
    captured: dict = {}

    class FakeChatOpenAI:
        def __init__(self, **kwargs):
            captured.clear()
            captured.update(kwargs)

    monkeypatch.setattr("langchain_openai.ChatOpenAI", FakeChatOpenAI)
    cfg = _llm_settings(
        default_llm_provider="nvidia_nim",
        nvidia_nim_api_key="nim-test-key",
        nvidia_nim_base_url=None,
        nvidia_nim_model="meta/llama-3.1-8b-instruct",
    )
    candidates = build_chat_model_candidates(cfg)
    assert len(candidates) == 1
    assert candidates[0].provider == "nvidia_nim"
    assert captured["base_url"] == DEFAULT_NIM_BASE_URL
    assert captured["base_url"] == "https://integrate.api.nvidia.com/v1"


def test_openai_preferred_uses_openai_path_not_nim_url(monkeypatch):
    builds: list[dict] = []

    class FakeChatOpenAI:
        def __init__(self, **kwargs):
            builds.append(dict(kwargs))

    monkeypatch.setattr("langchain_openai.ChatOpenAI", FakeChatOpenAI)
    cfg = _llm_settings(
        default_llm_provider="openai",
        openai_api_key="openai-test-key",
        openai_base_url=None,
        nvidia_nim_api_key="nim-test-key",
        nvidia_nim_base_url="https://integrate.api.nvidia.com/v1",
        nvidia_nim_model="meta/llama-3.1-8b-instruct",
        default_llm_model="gpt-4o-mini",
    )
    candidates = build_chat_model_candidates(cfg)
    assert [c.provider for c in candidates] == ["openai", "nvidia_nim"]
    openai_kwargs = builds[0]
    assert openai_kwargs["api_key"] == "openai-test-key"
    assert openai_kwargs["model"] == "gpt-4o-mini"
    assert "base_url" not in openai_kwargs
    nim_kwargs = builds[1]
    assert nim_kwargs["base_url"] == "https://integrate.api.nvidia.com/v1"


def test_no_keys_returns_no_candidates_and_deterministic_fallback():
    cfg = _llm_settings()
    assert cfg.has_llm() is False
    assert build_chat_model_candidates(cfg) == []

    records = [
        _evidence(
            source_id="sid-ncbi-1",
            section="General gene information",
            display_text="SREBF2 Entrez Gene ID is 6721.",
        )
    ]
    section, claims, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        settings=cfg,
        force_deterministic=False,
    )
    assert mode == "deterministic"
    assert section.status == "deterministic"
    assert len(claims) == 1


def test_force_deterministic_skips_all_llm_providers(monkeypatch):
    called = {"n": 0}

    def boom(*, model, gene_symbol, section_name, records):
        called["n"] += 1
        raise AssertionError("LLM must not be invoked when force_deterministic=True")

    monkeypatch.setattr("gene_dossier.synthesis._invoke_section_llm", boom)
    records = [
        _evidence(
            source_id="sid-ncbi-1",
            section="General gene information",
            display_text="SREBF2 Entrez Gene ID is 6721.",
        )
    ]
    candidates = [
        LlmModelCandidate(provider="openai", model=_TaggedModel("openai")),
        LlmModelCandidate(provider="nvidia_nim", model=_TaggedModel("nim")),
    ]
    section, claims, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        force_deterministic=True,
        model_candidates=candidates,
    )
    assert called["n"] == 0
    assert mode == "deterministic"
    assert section.status == "deterministic"
    assert len(claims) == 1


def test_anthropic_invoke_failure_falls_back_to_openai(monkeypatch):
    records = [
        _evidence(
            source_id="sid-a",
            section="General gene information",
            display_text="Official symbol SREBF2.",
        )
    ]
    tried: list[str] = []

    def fake_invoke(*, model, gene_symbol, section_name, records):
        tried.append(model.tag)
        if model.tag == "anthropic":
            raise RuntimeError("credit balance too low")
        return SectionDraft(
            section_id="general",
            summary_paragraphs=["SREBF2 is the official symbol [sid-a]."],
            claims=[
                SectionClaimDraft(
                    claim_text="Official symbol is SREBF2.",
                    supporting_evidence_ids=["sid-a"],
                )
            ],
        )

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    _, _, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        model_candidates=[
            LlmModelCandidate(provider="anthropic", model=_TaggedModel("anthropic")),
            LlmModelCandidate(provider="openai", model=_TaggedModel("openai")),
        ],
    )
    assert tried == ["anthropic", "openai"]
    assert mode == "llm"


def test_openai_invoke_failure_falls_back_to_nvidia_nim(monkeypatch):
    records = [
        _evidence(
            source_id="sid-a",
            section="General gene information",
            display_text="Official symbol SREBF2.",
        )
    ]
    tried: list[str] = []

    def fake_invoke(*, model, gene_symbol, section_name, records):
        tried.append(model.tag)
        if model.tag == "openai":
            raise RuntimeError("openai unavailable")
        return SectionDraft(
            section_id="general",
            summary_paragraphs=["SREBF2 is the official symbol [sid-a]."],
            claims=[
                SectionClaimDraft(
                    claim_text="Official symbol is SREBF2.",
                    supporting_evidence_ids=["sid-a"],
                )
            ],
        )

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    _, _, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        model_candidates=[
            LlmModelCandidate(provider="openai", model=_TaggedModel("openai")),
            LlmModelCandidate(provider="nvidia_nim", model=_TaggedModel("nim")),
        ],
    )
    assert tried == ["openai", "nim"]
    assert mode == "llm"


def test_all_providers_raise_during_invocation_uses_deterministic(monkeypatch):
    records = [
        _evidence(
            source_id="sid-a",
            section="General gene information",
            display_text="Official symbol SREBF2.",
        )
    ]
    tried: list[str] = []

    def fake_invoke(*, model, gene_symbol, section_name, records):
        tried.append(model.tag)
        raise RuntimeError(f"{model.tag} failed")

    monkeypatch.setattr(
        "gene_dossier.synthesis._invoke_section_llm", fake_invoke
    )
    section, claims, mode = synthesize_section(
        dossier_run_id="synth-run",
        gene_symbol="SREBF2",
        section_name="General gene information",
        records=records,
        model_candidates=[
            LlmModelCandidate(provider="openai", model=_TaggedModel("openai")),
            LlmModelCandidate(provider="nvidia_nim", model=_TaggedModel("nim")),
            LlmModelCandidate(provider="anthropic", model=_TaggedModel("anthropic")),
        ],
    )
    assert tried == ["openai", "nim", "anthropic"]
    assert mode == "deterministic"
    assert section.status == "deterministic"
    assert len(claims) == 1
    assert claims[0].source_ids == ["sid-a"]
