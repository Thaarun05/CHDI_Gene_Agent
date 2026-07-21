"""Normalize gene-identity and orthology ToolResults into EvidenceRecords.

Consumes successful client payloads from NCBI Gene, Ensembl, UniProt, UCSC, and
NCBI Datasets (orthologs). Does **not** call the network.

Skips ambiguous / unresolved identity selections (e.g. NCBI
``selection_method="ambiguous"``, UCSC ``selection_method="ambiguous"``).
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

SECTION_GENERAL = "General gene information"
SECTION_HOMOLOGUES = "Homologues"


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    source_type: SourceType,
    assertion_type: AssertionType,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    evidence_grade: EvidenceGrade,
    section: str,
    organism: str | None = None,
    species: str | None = None,
    taxon_id: int | None = None,
    official_symbol: str | None = None,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    raw_response_pointer: str | None = None,
) -> EvidenceRecord:
    """Build one EvidenceRecord with a deterministic source_id."""
    source_id = make_source_id(source_name, gene_symbol, assertion_type, key)
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=official_symbol or gene_symbol,
        section=section,
        subsection=subsection,
        source_name=source_name,
        source_type=source_type,
        assertion_type=assertion_type,
        fact_type=fact_type,
        organism=organism,
        species=species,
        taxon_id=taxon_id,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
        raw_response_pointer=raw_response_pointer,
    )


def normalize_ncbi_gene(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize NCBI Gene ``lookup_gene`` / identity payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    method = data.get("selection_method")
    gene_id = data.get("selected_gene_id")
    if method == "ambiguous" or not gene_id:
        return []

    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    organism = data.get("organism")
    taxon_id = data.get("expected_taxid")
    try:
        taxon_id_int = int(taxon_id) if taxon_id is not None else None
    except (TypeError, ValueError):
        taxon_id_int = None

    summary = data.get("selected_summary") if isinstance(data.get("selected_summary"), dict) else {}
    official = str(summary.get("nomenclaturesymbol") or summary.get("name") or gene_symbol)
    records = [
        _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            official_symbol=official,
            source_name="NCBI Gene",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="entrez_gene_id",
            key=str(gene_id),
            value={
                "entrez_gene_id": str(gene_id),
                "selection_method": method,
                "name": summary.get("name"),
                "nomenclaturesymbol": summary.get("nomenclaturesymbol"),
                "nomenclaturestatus": summary.get("nomenclaturestatus"),
                "description": summary.get("description") or summary.get("summary"),
                "chromosome": summary.get("chromosome"),
                "maplocation": summary.get("maplocation"),
            },
            display_text=f"{official} Entrez Gene ID is {gene_id}.",
            evidence_grade=EvidenceGrade.C,
            section=SECTION_GENERAL,
            organism=str(organism) if organism else None,
            taxon_id=taxon_id_int,
            confidence_notes="Selected via safe NCBI Gene symbol/taxid/status checks.",
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    ]
    return records


def normalize_ensembl(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize Ensembl ``lookup_symbol`` payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    summary = data.get("summary") if isinstance(data.get("summary"), dict) else {}
    if not summary and data.get("id"):
        summary = data
    ensembl_id = summary.get("ensembl_gene_id") or summary.get("id")
    if not ensembl_id:
        return []

    gene_symbol = tool_result.gene_symbol or str(
        summary.get("display_name") or data.get("gene_symbol") or ""
    )
    official = str(summary.get("display_name") or gene_symbol)
    chrom = summary.get("seq_region_name")
    start = summary.get("start")
    end = summary.get("end")
    transcript = summary.get("canonical_transcript")
    species = summary.get("species")

    records = [
        _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            official_symbol=official,
            source_name="Ensembl",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ensembl_gene_id",
            key=str(ensembl_id),
            value={
                "ensembl_gene_id": str(ensembl_id),
                "display_name": official,
                "biotype": summary.get("biotype"),
                "seq_region_name": chrom,
                "start": start,
                "end": end,
                "strand": summary.get("strand"),
                "assembly_name": summary.get("assembly_name"),
                "canonical_transcript": transcript,
                "description": summary.get("description"),
            },
            display_text=f"{official} Ensembl gene ID is {ensembl_id}.",
            evidence_grade=EvidenceGrade.C,
            section=SECTION_GENERAL,
            organism=str(species).replace("_", " ") if species else None,
            species=str(species) if species else None,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    ]
    if transcript:
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                official_symbol=official,
                source_name="Ensembl",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="canonical_transcript",
                key=str(transcript),
                value={"canonical_transcript": str(transcript), "ensembl_gene_id": str(ensembl_id)},
                display_text=(
                    f"{official} canonical transcript is {transcript}."
                ),
                evidence_grade=EvidenceGrade.C,
                section=SECTION_GENERAL,
                organism=str(species).replace("_", " ") if species else None,
                species=str(species) if species else None,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    if chrom is not None and start is not None and end is not None:
        loc_key = f"{chrom}:{start}-{end}"
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                official_symbol=official,
                source_name="Ensembl",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="genomic_location",
                key=loc_key,
                value={
                    "seq_region_name": chrom,
                    "start": start,
                    "end": end,
                    "strand": summary.get("strand"),
                    "assembly_name": summary.get("assembly_name"),
                },
                display_text=f"{official} is located at {loc_key}.",
                evidence_grade=EvidenceGrade.C,
                section=SECTION_GENERAL,
                organism=str(species).replace("_", " ") if species else None,
                species=str(species) if species else None,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_uniprot(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize UniProt search / summarized entry payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    entries = data.get("entries") or data.get("results") or []
    if isinstance(data.get("selected_entry"), dict):
        entries = [data["selected_entry"]]
    if not isinstance(entries, list):
        # Single summarized entry at top level
        if data.get("accession") or data.get("primaryAccession"):
            entries = [data]
        else:
            entries = []

    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    selected_accession = data.get("selected_accession")
    records: list[EvidenceRecord] = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        accession = (
            entry.get("accession")
            or entry.get("primaryAccession")
            or entry.get("uni_protkb_id")
        )
        if not accession:
            continue
        # Prefer the client's selected accession when present; skip other hits.
        if selected_accession and str(accession) != str(selected_accession):
            continue
        protein_name = entry.get("protein_name") or entry.get("uniProtkbId")
        organism = entry.get("organism_name")
        organism_id = entry.get("organism_id")
        try:
            taxon_id = int(organism_id) if organism_id is not None else None
        except (TypeError, ValueError):
            taxon_id = None
        gene_names = entry.get("gene_names") or []
        official = gene_symbol
        if isinstance(gene_names, list) and gene_names:
            first = gene_names[0]
            if isinstance(first, str):
                official = first
            elif isinstance(first, dict):
                official = str(first.get("geneName") or first.get("value") or gene_symbol)

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                official_symbol=str(official),
                source_name="UniProt",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="uniprot_accession",
                key=str(accession),
                value={
                    "uniprot_accession": str(accession),
                    "protein_name": protein_name,
                    "gene_names": gene_names,
                    "organism_name": organism,
                    "organism_id": organism_id,
                    "ensembl_xrefs": entry.get("ensembl_xrefs"),
                },
                display_text=(
                    f"{official} UniProt accession is {accession}"
                    + (f" ({protein_name})." if protein_name else ".")
                ),
                evidence_grade=EvidenceGrade.C,
                section=SECTION_GENERAL,
                organism=str(organism) if organism else None,
                taxon_id=taxon_id,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_ucsc(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize UCSC ``fetch_gene_region`` payloads (safe matches only)."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    if data.get("selection_method") != "matched":
        return []
    chrom = data.get("chrom")
    start = data.get("start")
    end = data.get("end")
    if not chrom or start is None or end is None:
        return []

    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    genome = data.get("genome") or "hg38"
    position = data.get("position") or f"{chrom}:{start}-{end}"
    browser_url = data.get("browser_url")
    return [
        _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            source_name="UCSC",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="genomic_region",
            key=str(position),
            value={
                "genome": genome,
                "chrom": chrom,
                "start": start,
                "end": end,
                "position": position,
                "browser_url": browser_url,
                "selection_method": "matched",
            },
            display_text=f"{gene_symbol} UCSC/{genome} region is {position}.",
            evidence_grade=EvidenceGrade.C,
            section=SECTION_GENERAL,
            subsection="UCSC genome browser",
            organism="Homo sapiens" if str(genome).startswith("hg") else None,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    ]


def normalize_ncbi_datasets_orthologs(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize NCBI Datasets ortholog summaries into orthology evidence."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    query_gene_id = data.get("gene_id")
    summaries = data.get("ortholog_summaries") or []
    if not isinstance(summaries, list):
        summaries = []

    records: list[EvidenceRecord] = []
    for row in summaries:
        if not isinstance(row, dict):
            continue
        orth_id = row.get("gene_id")
        orth_symbol = row.get("symbol")
        if not orth_id:
            continue
        tax_id = row.get("tax_id")
        try:
            taxon_id = int(tax_id) if tax_id is not None else None
        except (TypeError, ValueError):
            taxon_id = None
        scientific = row.get("scientific_name")
        key = f"{query_gene_id or gene_symbol}-{orth_id}"
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="NCBI Datasets",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.orthology,
                fact_type="ortholog_gene",
                key=str(key),
                value={
                    "query_gene_id": query_gene_id,
                    "ortholog_gene_id": orth_id,
                    "ortholog_symbol": orth_symbol,
                    "tax_id": tax_id,
                    "common_name": row.get("common_name"),
                    "scientific_name": scientific,
                    "chromosomes": row.get("chromosomes"),
                },
                display_text=(
                    f"{gene_symbol} has NCBI Datasets ortholog "
                    f"{orth_symbol or orth_id}"
                    + (f" ({scientific})" if scientific else "")
                    + f" [Gene ID {orth_id}]."
                ),
                evidence_grade=EvidenceGrade.C,
                section=SECTION_HOMOLOGUES,
                organism=str(scientific) if scientific else None,
                taxon_id=taxon_id,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_gene_identity(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    if source == "NCBI Gene":
        return normalize_ncbi_gene(tool_result, **kwargs)
    if source == "Ensembl":
        return normalize_ensembl(tool_result, **kwargs)
    if source == "UniProt":
        return normalize_uniprot(tool_result, **kwargs)
    if source == "UCSC":
        return normalize_ucsc(tool_result, **kwargs)
    if source == "NCBI Datasets":
        return normalize_ncbi_datasets_orthologs(tool_result, **kwargs)
    return []


__all__ = [
    "normalize_ncbi_gene",
    "normalize_ensembl",
    "normalize_uniprot",
    "normalize_ucsc",
    "normalize_ncbi_datasets_orthologs",
    "normalize_gene_identity",
]
