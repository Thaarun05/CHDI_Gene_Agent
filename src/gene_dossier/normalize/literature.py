"""Normalize PubMed ToolResults into literature EvidenceRecords.

Consumes successful ``search_hd_literature`` (or compatible) payloads. Does **not**
call the network.

Rules:
- One EvidenceRecord per PMID with ESummary metadata
- Do **not** treat a search hit as proof of a gene–HD causal link
- Do **not** invent abstracts, claim strength, or lab affiliations
- Preserve the client caveat in confidence notes
"""

from __future__ import annotations

from typing import Any

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
    ToolResult,
)
from gene_dossier.source_ids import make_source_id

SECTION_LITERATURE = "Major labs / literature"

CAVEAT = "Do not assume every hit strongly supports the gene–HD link."


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    evidence_grade: EvidenceGrade,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    raw_response_pointer: str | None = None,
) -> EvidenceRecord:
    """Build one PubMed literature EvidenceRecord."""
    source_id = make_source_id("PubMed", gene_symbol, AssertionType.literature_summary, key)
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_LITERATURE,
        subsection=subsection,
        source_name="PubMed",
        source_type=SourceType.literature,
        assertion_type=AssertionType.literature_summary,
        fact_type=fact_type,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
        raw_response_pointer=raw_response_pointer,
    )


def _esummary_uid_map(esummary_payload: Any) -> dict[str, dict[str, Any]]:
    """Map PMID -> ESummary entry from a raw ESummary JSON body."""
    payload = _as_dict(esummary_payload)
    result = payload.get("result")
    if not isinstance(result, dict):
        return {}
    out: dict[str, dict[str, Any]] = {}
    for uid, entry in result.items():
        if uid == "uids" or not isinstance(entry, dict):
            continue
        out[str(uid)] = entry
    return out


def _author_names(entry: dict[str, Any]) -> list[str]:
    """Extract author display names from an ESummary entry."""
    authors = entry.get("authors") or []
    if not isinstance(authors, list):
        return []
    names: list[str] = []
    for author in authors:
        if isinstance(author, dict) and author.get("name"):
            names.append(str(author["name"]))
        elif isinstance(author, str) and author.strip():
            names.append(author.strip())
    return names


def _summarize_article(pmid: str, entry: dict[str, Any] | None) -> dict[str, Any]:
    """Pull stable bibliographic fields from one ESummary entry (or stubs)."""
    if not entry:
        return {
            "pmid": str(pmid),
            "title": None,
            "authors": [],
            "source": None,
            "fulljournalname": None,
            "pubdate": None,
            "epubdate": None,
            "doi": None,
            "elocationid": None,
            "pubtype": [],
        }
    authors = _author_names(entry)
    pubtypes = entry.get("pubtype") or []
    if not isinstance(pubtypes, list):
        pubtypes = [str(pubtypes)] if pubtypes else []
    elocation = entry.get("elocationid")
    doi = None
    if isinstance(elocation, str) and "doi:" in elocation.lower():
        doi = elocation.split(":", 1)[-1].strip()
    articleids = entry.get("articleids") or []
    if doi is None and isinstance(articleids, list):
        for item in articleids:
            if isinstance(item, dict) and str(item.get("idtype") or "").lower() == "doi":
                if item.get("value"):
                    doi = str(item["value"])
                    break
    return {
        "pmid": str(pmid),
        "title": entry.get("title"),
        "authors": authors,
        "source": entry.get("source"),
        "fulljournalname": entry.get("fulljournalname"),
        "pubdate": entry.get("pubdate"),
        "epubdate": entry.get("epubdate"),
        "doi": doi,
        "elocationid": elocation,
        "pubtype": [str(p) for p in pubtypes],
    }


def normalize_pubmed(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize PubMed ``search_hd_literature`` payloads into literature records.

    Each PMID becomes one :class:`EvidenceRecord` with grade F and
    ``manual_review_required=True``. Search membership is not treated as causal
    gene–HD evidence.
    """
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    if not gene_symbol:
        return []

    pmids = data.get("pmids") or []
    if not isinstance(pmids, list):
        pmids = []
    pmids = [str(p) for p in pmids if p is not None and str(p).strip()]
    if not pmids:
        return []

    uid_map = _esummary_uid_map(data.get("esummary"))
    search_term = data.get("search_term")
    caveat = str(data.get("caveat") or CAVEAT)
    count = data.get("count")

    records: list[EvidenceRecord] = []
    for pmid in pmids:
        entry = uid_map.get(pmid)
        summary = _summarize_article(pmid, entry)
        title = summary.get("title")
        pubdate = summary.get("pubdate")
        journal = summary.get("fulljournalname") or summary.get("source")
        display_bits = [f"PMID {pmid}"]
        if title:
            display_bits.append(str(title).rstrip("."))
        if journal or pubdate:
            meta = ", ".join(x for x in (journal, pubdate) if x)
            if meta:
                display_bits.append(f"({meta})")
        display_text = ". ".join(display_bits) + "."
        if not title:
            display_text = (
                f"PMID {pmid} retrieved by PubMed HD literature search"
                f" for {gene_symbol}."
            )

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                fact_type="pubmed_article",
                key=f"PMID:{pmid}",
                value={
                    **summary,
                    "search_term": search_term,
                    "search_count": count,
                    "hd_mesh_search": True,
                    "caveat": caveat,
                },
                display_text=display_text,
                evidence_grade=EvidenceGrade.F,
                subsection="PubMed HD search",
                confidence_notes=caveat,
                manual_review_required=True,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_literature(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch literature normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    if source == "PubMed":
        return normalize_pubmed(
            tool_result,
            dossier_run_id=dossier_run_id,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    return []


__all__ = [
    "normalize_pubmed",
    "normalize_literature",
]
