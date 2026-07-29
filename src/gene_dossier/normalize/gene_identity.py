"""Normalize gene-identity and orthology ToolResults into EvidenceRecords.

Consumes successful client payloads from NCBI Gene, Ensembl, UniProt, UCSC, and
NCBI Datasets (orthologs). Does **not** call the network.

Skips ambiguous / unresolved identity selections (e.g. NCBI
``selection_method="ambiguous"``, UCSC ``selection_method="ambiguous"``).
"""

from __future__ import annotations

import re
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

ENSEMBL_SPECIES_TAXON: dict[str, int] = {
    "homo_sapiens": 9606,
    "mus_musculus": 10090,
    "rattus_norvegicus": 10116,
}

TAXON_SCIENTIFIC: dict[int, str] = {
    9606: "Homo sapiens",
    10090: "Mus musculus",
    10116: "Rattus norvegicus",
}

TAXON_COMMON: dict[int, str] = {
    9606: "human",
    10090: "mouse",
    10116: "rat",
}


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def split_alias_field(raw: Any) -> list[str]:
    """Split NCBI/UniProt alias fields into trimmed unique tokens."""
    if raw is None:
        return []
    items: list[str] = []
    if isinstance(raw, list):
        for entry in raw:
            if entry is None:
                continue
            if isinstance(entry, dict):
                text = entry.get("value") or entry.get("synonym") or entry.get("name")
                if text:
                    items.append(str(text))
            else:
                items.extend(re.split(r"[,;|]", str(entry)))
    else:
        items.extend(re.split(r"[,;|]", str(raw)))
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        token = item.strip()
        if not token:
            continue
        key = token.casefold()
        if key in seen:
            continue
        seen.add(key)
        out.append(token)
    return out


def ensembl_gene_id_from_xref(xref: str) -> str | None:
    """Extract an Ensembl gene ID from a UniProt Ensembl cross-reference.

    Rejects transcript (ENST/ENSMUST/ENSRNOT) and protein IDs.
    """
    base = str(xref or "").strip().split(".", 1)[0]
    if not base:
        return None
    if base.startswith(("ENST", "ENSMUST", "ENSRNOT", "ENSP", "ENSMUSP", "ENSRNOP")):
        return None
    if base.startswith(("ENSG", "ENSMUSG", "ENSRNOG")):
        return base
    return None


def _dossier_identity_symbols(
    tool_result: ToolResult,
    data: dict[str, Any],
    *,
    species_official: str | None = None,
) -> tuple[str, str]:
    """Return ``(query_gene_symbol, species_gene_symbol)`` for identity records.

    The dossier query symbol (e.g. ``SREBF2``) is preserved on every species
    EvidenceRecord. The species-specific symbol (e.g. ``Srebf2``) is stored
    separately as ``official_symbol`` / ``value.species_gene_symbol``.
    """
    params = (
        tool_result.request_params
        if isinstance(tool_result.request_params, dict)
        else {}
    )
    query = (
        params.get("query_symbol")
        or data.get("query_gene_symbol")
        or data.get("query_symbol")
        or tool_result.gene_symbol
        or ""
    )
    species = (
        params.get("resolved_symbol")
        or data.get("species_gene_symbol")
        or data.get("resolved_symbol")
        or species_official
        or query
    )
    query_s = str(query).strip()
    species_s = str(species).strip() or query_s
    return query_s, species_s


def _with_identity_symbol_fields(
    value: dict[str, Any],
    *,
    query_gene_symbol: str,
    species_gene_symbol: str,
) -> dict[str, Any]:
    """Attach dossier query / species symbol fields without mutating callers."""
    out = dict(value)
    out["query_gene_symbol"] = query_gene_symbol
    out["species_gene_symbol"] = species_gene_symbol
    return out


def _safe_source_url(tool_result: ToolResult) -> str | None:
    """Return request URL with credential-bearing query params stripped."""
    url = (tool_result.request_url or "").strip()
    if not url:
        return None
    for marker in ("api_key=", "apiKey=", "accesskey=", "api-key="):
        if marker in url:
            return url.split("?", 1)[0]
    return url


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
    summary = (
        data.get("selected_summary")
        if isinstance(data.get("selected_summary"), dict)
        else {}
    )
    organism_block = (
        summary.get("organism") if isinstance(summary.get("organism"), dict) else {}
    )
    if taxon_id is None:
        taxon_id = organism_block.get("taxid")
    try:
        taxon_id_int = int(taxon_id) if taxon_id is not None else None
    except (TypeError, ValueError):
        taxon_id_int = None

    scientific_name = (
        organism_block.get("scientificname")
        or organism_block.get("scientificName")
        or (TAXON_SCIENTIFIC.get(taxon_id_int) if taxon_id_int else None)
        or (str(organism) if organism else None)
    )
    common_name = (
        organism_block.get("commonname")
        or organism_block.get("commonName")
        or (TAXON_COMMON.get(taxon_id_int) if taxon_id_int else None)
    )

    official = str(
        summary.get("nomenclaturesymbol") or summary.get("name") or gene_symbol
    )
    query_symbol, species_symbol = _dossier_identity_symbols(
        tool_result, data, species_official=official
    )
    gene_name = (
        summary.get("nomenclaturename")
        or summary.get("description")
        or summary.get("summary")
    )
    aliases = split_alias_field(summary.get("otheraliases"))
    other_designations = split_alias_field(summary.get("otherdesignations"))
    source_url = _safe_source_url(tool_result)

    value = _with_identity_symbol_fields(
        {
            "entrez_gene_id": str(gene_id),
            "gene_symbol": species_symbol,
            "gene_name": gene_name,
            "selection_method": method,
            "name": summary.get("name"),
            "nomenclaturesymbol": summary.get("nomenclaturesymbol"),
            "nomenclaturename": summary.get("nomenclaturename"),
            "nomenclaturestatus": summary.get("nomenclaturestatus"),
            "description": gene_name,
            "otheraliases": aliases,
            "aliases": aliases,
            "otherdesignations": other_designations,
            "tax_id": taxon_id_int,
            "taxon_id": taxon_id_int,
            "scientific_name": scientific_name,
            "common_name": common_name,
            "status": summary.get("status"),
            "currentid": summary.get("currentid"),
            "chromosome": summary.get("chromosome"),
            "maplocation": summary.get("maplocation"),
            "source_url": source_url,
            "raw_artifact_id": raw_artifact_id,
        },
        query_gene_symbol=query_symbol,
        species_gene_symbol=species_symbol,
    )

    return [
        _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=query_symbol,
            official_symbol=species_symbol,
            source_name="NCBI Gene",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="entrez_gene_id",
            key=str(gene_id),
            value=value,
            display_text=f"{species_symbol} Entrez Gene ID is {gene_id}.",
            evidence_grade=EvidenceGrade.C,
            section=SECTION_GENERAL,
            organism=str(scientific_name)
            if scientific_name
            else (str(organism) if organism else None),
            taxon_id=taxon_id_int,
            confidence_notes="Selected via safe NCBI Gene symbol/taxid/status checks.",
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    ]


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
    ensembl_id = (
        summary.get("ensembl_gene_id")
        or summary.get("id")
        or data.get("ensembl_id")
        or data.get("id")
    )
    if not ensembl_id:
        return []

    display = str(summary.get("display_name") or "").strip()
    query_symbol, species_symbol = _dossier_identity_symbols(
        tool_result, data, species_official=display or None
    )
    if display:
        species_symbol = display
    chrom = summary.get("seq_region_name")
    start = summary.get("start")
    end = summary.get("end")
    transcript = summary.get("canonical_transcript")
    species = summary.get("species") or data.get("species")
    taxon_id_int = None
    if isinstance(species, str) and species in ENSEMBL_SPECIES_TAXON:
        taxon_id_int = ENSEMBL_SPECIES_TAXON[species]
    scientific_name = (
        TAXON_SCIENTIFIC.get(taxon_id_int)
        if taxon_id_int
        else (str(species).replace("_", " ") if species else None)
    )
    common_name = TAXON_COMMON.get(taxon_id_int) if taxon_id_int else None
    source_url = _safe_source_url(tool_result)
    raw = data.get("raw") if isinstance(data.get("raw"), dict) else {}
    version = summary.get("version") or raw.get("version")
    source = summary.get("source") or raw.get("source") or "ensembl"

    value = _with_identity_symbol_fields(
        {
        "ensembl_gene_id": str(ensembl_id),
        "id": str(ensembl_id),
        "display_name": species_symbol,
        "gene_symbol": species_symbol,
        "biotype": summary.get("biotype"),
        "seq_region_name": chrom,
        "start": start,
        "end": end,
        "strand": summary.get("strand"),
        "assembly_name": summary.get("assembly_name"),
        "canonical_transcript": transcript,
        "description": summary.get("description"),
        "species": species,
        "tax_id": taxon_id_int,
        "taxon_id": taxon_id_int,
        "scientific_name": scientific_name,
        "common_name": common_name,
        "version": version,
        "source": source,
        "source_url": source_url,
        "raw_artifact_id": raw_artifact_id,
        },
        query_gene_symbol=query_symbol,
        species_gene_symbol=species_symbol,
    )

    records = [
        _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=query_symbol,
            official_symbol=species_symbol,
            source_name="Ensembl",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ensembl_gene_id",
            key=str(ensembl_id),
            value=value,
            display_text=f"{species_symbol} Ensembl gene ID is {ensembl_id}.",
            evidence_grade=EvidenceGrade.C,
            section=SECTION_GENERAL,
            organism=scientific_name,
            species=str(species) if species else None,
            taxon_id=taxon_id_int,
            confidence_notes="Direct Ensembl symbol lookup (primary Ensembl gene ID).",
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    ]
    if transcript:
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=query_symbol,
                official_symbol=species_symbol,
                source_name="Ensembl",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="canonical_transcript",
                key=str(transcript),
                value=_with_identity_symbol_fields(
                    {
                    "canonical_transcript": str(transcript),
                    "ensembl_gene_id": str(ensembl_id),
                    "tax_id": taxon_id_int,
                    "taxon_id": taxon_id_int,
                    "scientific_name": scientific_name,
                    "source_url": source_url,
                    },
                    query_gene_symbol=query_symbol,
                    species_gene_symbol=species_symbol,
                ),
                display_text=f"{species_symbol} canonical transcript is {transcript}.",
                evidence_grade=EvidenceGrade.C,
                section=SECTION_GENERAL,
                organism=scientific_name,
                species=str(species) if species else None,
                taxon_id=taxon_id_int,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    if chrom is not None and start is not None and end is not None:
        loc_key = f"{chrom}:{start}-{end}"
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=query_symbol,
                official_symbol=species_symbol,
                source_name="Ensembl",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="genomic_location",
                key=loc_key,
                value=_with_identity_symbol_fields(
                    {
                    "seq_region_name": chrom,
                    "start": start,
                    "end": end,
                    "strand": summary.get("strand"),
                    "assembly_name": summary.get("assembly_name"),
                    "tax_id": taxon_id_int,
                    "taxon_id": taxon_id_int,
                    "scientific_name": scientific_name,
                    "source_url": source_url,
                    },
                    query_gene_symbol=query_symbol,
                    species_gene_symbol=species_symbol,
                ),
                display_text=f"{species_symbol} is located at {loc_key}.",
                evidence_grade=EvidenceGrade.C,
                section=SECTION_GENERAL,
                organism=scientific_name,
                species=str(species) if species else None,
                taxon_id=taxon_id_int,
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
        if data.get("accession") or data.get("primaryAccession"):
            entries = [data]
        else:
            entries = []

    selected_accession = data.get("selected_accession")
    default_organism_id = data.get("organism_id")
    source_url = _safe_source_url(tool_result)
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
        if selected_accession and str(accession) != str(selected_accession):
            continue
        protein_name = entry.get("protein_name") or entry.get("uniProtkbId")
        organism = entry.get("organism_name")
        organism_id = entry.get("organism_id")
        if organism_id is None:
            organism_id = default_organism_id
        try:
            taxon_id = int(organism_id) if organism_id is not None else None
        except (TypeError, ValueError):
            taxon_id = None
        scientific_name = organism or (
            TAXON_SCIENTIFIC.get(taxon_id) if taxon_id else None
        )
        common_name = TAXON_COMMON.get(taxon_id) if taxon_id else None
        gene_names = entry.get("gene_names") or []
        if not isinstance(gene_names, list):
            gene_names = []
        str_names = [str(n).strip() for n in gene_names if str(n).strip()]
        official_from_entry = ""
        if str_names:
            official_from_entry = str_names[0]
        elif gene_names:
            first = gene_names[0]
            if isinstance(first, str):
                official_from_entry = first
            elif isinstance(first, dict):
                official_from_entry = str(
                    first.get("geneName") or first.get("value") or ""
                )
        query_symbol, species_symbol = _dossier_identity_symbols(
            tool_result, data, species_official=official_from_entry or None
        )
        if official_from_entry:
            species_symbol = official_from_entry
        if entry.get("gene_synonyms"):
            aliases = split_alias_field(entry.get("gene_synonyms"))
        else:
            aliases = [n for n in str_names[1:] if n]
        reviewed = entry.get("reviewed")
        if reviewed is None:
            reviewed = True
        ensembl_gene_ids = entry.get("ensembl_gene_ids") or entry.get("ensembl_xrefs") or []
        if not isinstance(ensembl_gene_ids, list):
            ensembl_gene_ids = []
        # Keep only true gene IDs (never transcripts).
        ensembl_xrefs = []
        for xref in ensembl_gene_ids:
            gene_ens = ensembl_gene_id_from_xref(str(xref))
            if gene_ens and gene_ens not in ensembl_xrefs:
                ensembl_xrefs.append(gene_ens)
        ensembl_transcript_ids = entry.get("ensembl_transcript_ids") or []
        if not isinstance(ensembl_transcript_ids, list):
            ensembl_transcript_ids = []
        refseq_protein_accessions = entry.get("refseq_protein_accessions") or []
        if not isinstance(refseq_protein_accessions, list):
            refseq_protein_accessions = []
        protein_length = (
            entry.get("protein_length")
            or entry.get("sequence_length")
            or entry.get("length")
        )
        try:
            protein_length_int = (
                int(protein_length) if protein_length is not None else None
            )
        except (TypeError, ValueError):
            protein_length_int = None

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=query_symbol,
                official_symbol=str(species_symbol),
                source_name="UniProt",
                source_type=SourceType.curated_database,
                assertion_type=AssertionType.gene_identity,
                fact_type="uniprot_accession",
                key=str(accession),
                value=_with_identity_symbol_fields(
                    {
                    "uniprot_accession": str(accession),
                    "primaryAccession": str(accession),
                    "reviewed": bool(reviewed),
                    "protein_name": protein_name,
                    "protein_length": protein_length_int,
                    "sequence_length": protein_length_int,
                    "gene_names": str_names or gene_names,
                    "gene_name": str(species_symbol),
                    "aliases": aliases,
                    "gene_synonyms": aliases,
                    "refseq_protein_accessions": [
                        str(x).strip()
                        for x in refseq_protein_accessions
                        if str(x).strip()
                    ],
                    "organism_name": organism,
                    "organism_id": organism_id,
                    "tax_id": taxon_id,
                    "taxon_id": taxon_id,
                    "scientific_name": scientific_name,
                    "common_name": common_name,
                    "ensembl_xrefs": ensembl_xrefs,
                    "ensembl_gene_ids": ensembl_xrefs,
                    "ensembl_transcript_ids": ensembl_transcript_ids,
                    "source_url": source_url,
                    "raw_artifact_id": raw_artifact_id,
                    },
                    query_gene_symbol=query_symbol,
                    species_gene_symbol=str(species_symbol),
                ),
                display_text=(
                    f"{species_symbol} UniProt accession is {accession}"
                    + (f" ({protein_name})." if protein_name else ".")
                ),
                evidence_grade=EvidenceGrade.C,
                section=SECTION_GENERAL,
                organism=str(scientific_name) if scientific_name else None,
                taxon_id=taxon_id,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
        for xref in ensembl_xrefs:
            gene_ens = ensembl_gene_id_from_xref(str(xref))
            if not gene_ens:
                continue
            records.append(
                _record(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=query_symbol,
                    official_symbol=str(species_symbol),
                    source_name="UniProt",
                    source_type=SourceType.curated_database,
                    assertion_type=AssertionType.gene_identity,
                    fact_type="ensembl_xref",
                    key=f"{accession}:{gene_ens}",
                    value=_with_identity_symbol_fields(
                        {
                        "ensembl_gene_id": gene_ens,
                        "ensembl_xref": str(xref),
                        "uniprot_accession": str(accession),
                        "primary": False,
                        "source": "uniprot_xref",
                        "tax_id": taxon_id,
                        "taxon_id": taxon_id,
                        "scientific_name": scientific_name,
                        "common_name": common_name,
                        "source_url": source_url,
                        "raw_artifact_id": raw_artifact_id,
                        },
                        query_gene_symbol=query_symbol,
                        species_gene_symbol=str(species_symbol),
                    ),
                    display_text=(
                        f"{species_symbol} UniProt {accession} cross-references "
                        f"Ensembl {gene_ens}."
                    ),
                    evidence_grade=EvidenceGrade.C,
                    section=SECTION_GENERAL,
                    organism=str(scientific_name) if scientific_name else None,
                    taxon_id=taxon_id,
                    confidence_notes=(
                        "Secondary Ensembl evidence from UniProt xref; "
                        "direct Ensembl lookup remains primary for Section 1a."
                    ),
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
    """Normalize UCSC payloads.

    Conservation ToolResults (search + track and/or figure) are delegated to
    ``normalize.ucsc_conservation``. Legacy matched region-only payloads still
    emit a single ``genomic_region`` record for backward compatibility.
    """
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    has_conservation_shape = isinstance(data.get("search"), dict) or isinstance(
        data.get("track_data"), dict
    ) or isinstance(data.get("figure"), dict)
    if has_conservation_shape and (
        data.get("search") is not None or data.get("track_data") is not None
    ):
        from gene_dossier.normalize.ucsc_conservation import normalize_ucsc_conservation

        return normalize_ucsc_conservation(
            tool_result,
            dossier_run_id=dossier_run_id,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )

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
    "ENSEMBL_SPECIES_TAXON",
    "TAXON_SCIENTIFIC",
    "TAXON_COMMON",
    "split_alias_field",
    "ensembl_gene_id_from_xref",
    "normalize_ncbi_gene",
    "normalize_ensembl",
    "normalize_uniprot",
    "normalize_ucsc",
    "normalize_ncbi_datasets_orthologs",
    "normalize_gene_identity",
]
