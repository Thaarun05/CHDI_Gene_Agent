"""Polished Rancho/CHDI section presentation builders.

Separates human-facing presentation blocks from audit evidence blocks.
Currently implements Section 1a (Gene Aliases) and Section 1b (UCSC conservation).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Sequence

from gene_dossier.models import AssertionType, EvidenceRecord
from gene_dossier.report_schema import ReportContentBlock

NOT_AVAILABLE = "Not available"

TAXON_HUMAN = 9606
TAXON_MOUSE = 10090
TAXON_RAT = 10116

_ENTREZ_RE = re.compile(r"^\d+$")
_ENSEMBL_RE = re.compile(r"^ENS(?:[A-Z]{0,4})?G\d+(?:\.\d+)?$", re.IGNORECASE)
_ENSEMBL_BASE_RE = re.compile(r"^(ENS(?:[A-Z]{0,4})?G\d+)", re.IGNORECASE)
# Canonical UniProt accession (not isoform); isoform has trailing -N.
_UNIPROT_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
_UNIPROT_ISOFORM_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-\d+$"
)

_CELL_LINK_RE = re.compile(r"^\[([^\]]+)\]\((https://[^)]+)\)$")

ALLOWED_LINK_HOSTS = frozenset(
    {
        "www.ncbi.nlm.nih.gov",
        "www.ensembl.org",
        "www.uniprot.org",
        "genome.ucsc.edu",
        "alphafold.ebi.ac.uk",
        "www.alphafold.ebi.ac.uk",
        "alphafold.com",
        "www.alphafold.com",
        "www.orthodb.org",
        "orthodb.org",
        "www.gtexportal.org",
        "gtexportal.org",
        "hbatlas.org",
        "www.hbatlas.org",
        "brain-map.org",
        "www.brain-map.org",
        "celltypes.brain-map.org",
        "transcriptomics.brain-map.org",
        "knowledge.brain-map.org",
        "portal.brain-map.org",
    }
)

UCSC_STABLE_INTRO = (
    "The UCSC Genome Browser hosts the human genome alongside multiple genome "
    "assemblies and a large collection of vertebrate and model-organism annotations. "
    "It provides a graphical interface for viewing, analyzing, and downloading "
    "genomic data."
)


@dataclass(frozen=True)
class PresentationDiagnostic:
    """Diagnostic outside the polished report body."""

    field: str
    reason: str
    severity: Literal["info", "warning"]


@dataclass(frozen=True)
class SectionPresentationResult:
    """Polished blocks plus non-rendered diagnostics."""

    blocks: tuple[ReportContentBlock, ...]
    diagnostics: tuple[PresentationDiagnostic, ...]


@dataclass
class _Candidate:
    display: str
    source_name: str
    source_id: str | None = None
    evidence_id: str | None = None


@dataclass
class _SpeciesFields:
    entrez: _Candidate | None = None
    symbol: _Candidate | None = None
    name: _Candidate | None = None
    ensembl: _Candidate | None = None
    uniprot: _Candidate | None = None
    aliases: list[str] = field(default_factory=list)
    source_ids: list[str] = field(default_factory=list)
    evidence_ids: list[str] = field(default_factory=list)


def build_section_presentation(
    *,
    section_key: str,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build polished presentation blocks for a known section key.

    Unknown keys return an empty result without raising.
    """
    key = (section_key or "").strip().lower()
    if key in {"1a", "1.a", "gene_aliases", "gene-aliases"}:
        return build_gene_aliases_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
        )
    if key in {
        "1b",
        "1.b",
        "conservation",
        "ucsc_conservation",
        "gene_conservation",
    }:
        return build_conservation_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
        )
    if key in {
        "1c",
        "1.c",
        "known_structure",
        "known-structure",
        "known_structure_domains",
    }:
        return build_known_structure_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
        )
    if key in {
        "1d",
        "1.d",
        "alphafold",
        "alphafold_prediction",
        "alphafold-protein-structure-prediction",
    }:
        return build_alphafold_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    if key in {
        "1e",
        "1.e",
        "homologues",
        "homologues_in_model_animals",
        "orthologs",
    }:
        return build_homologues_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    if key in {
        "2a",
        "2.a",
        "tissue_specific",
        "tissue-specific",
        "tissue_specific_information",
    }:
        return build_tissue_specific_information_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    if key in {
        "2b",
        "2.b",
        "barres",
        "barres_lab",
        "barres_lab_rna_seq",
        "barres_lab_rna-seq_brain_specific_expression_data",
    }:
        return build_barres_brain_expression_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    if key in {
        "2c",
        "2.c",
        "snrna_seq",
        "snrna-seq",
        "cell_type_database",
        "snrna_seq_gene_expression_in_cell_type_database",
    }:
        return build_section_2c_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    if key in {
        "3a",
        "3.a",
        "geo_profiles",
        "geo-profiles",
        "geo_profiles_search",
        "geo_profiles_search_focusing_on_brain_and_or_neurons",
    }:
        return build_section_3a_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    if key in {
        "4a",
        "4.a",
        "harmonizome",
        "transcription_factors",
        "harmonizome_integrated_knowledge_about_genes_and_proteins",
    }:
        return build_section_4a_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    if key in {
        "5a",
        "5.a",
        "string",
        "string_db",
        "protein_protein_interaction",
        "ppi_partners_string",
    }:
        return build_section_5a_blocks(
            gene_symbol=gene_symbol,
            evidence_records=evidence_records,
            section_status=section_status,
        )
    return SectionPresentationResult(blocks=(), diagnostics=())


def build_gene_aliases_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
) -> SectionPresentationResult:
    """Build the Section 1a Human/Mouse/Rat identity table."""
    records = list(evidence_records)
    diagnostics: list[PresentationDiagnostic] = []

    human = _build_human(records, diagnostics)
    mouse = _build_mouse(records, diagnostics)
    rat = _build_rat(records, diagnostics)

    for species_key, species in (("human", human), ("mouse", mouse), ("rat", rat)):
        _emit_missing(species_key, species, diagnostics)

    headers = ["", "Human gene", "Mouse gene", "Rat gene"]
    rows = [
        [
            "Entrez Gene ID",
            _linked_entrez(human.entrez),
            _linked_entrez(mouse.entrez),
            _linked_entrez(rat.entrez),
        ],
        [
            "Gene Symbol",
            _plain(human.symbol),
            _plain(mouse.symbol),
            _plain(rat.symbol),
        ],
        [
            "Gene Name",
            _plain(human.name),
            _plain(mouse.name),
            _plain(rat.name),
        ],
        [
            "Ensembl ID",
            _linked_ensembl(human.ensembl),
            _linked_ensembl(mouse.ensembl),
            _linked_ensembl(rat.ensembl),
        ],
        [
            "UniProt ID",
            _linked_uniprot(human.uniprot),
            _linked_uniprot(mouse.uniprot),
            _linked_uniprot(rat.uniprot),
        ],
        [
            # Zero-width space after slash: wrap at "/" never mid-"Synonyms".
            "Synonyms/\u200bAliases",
            _aliases_cell(human.aliases),
            _aliases_cell(mouse.aliases),
            _aliases_cell(rat.aliases),
        ],
    ]

    source_ids = list(
        dict.fromkeys(human.source_ids + mouse.source_ids + rat.source_ids)
    )
    evidence_ids = list(
        dict.fromkeys(human.evidence_ids + mouse.evidence_ids + rat.evidence_ids)
    )

    block = ReportContentBlock(
        kind="table",
        presentation_role="gene_aliases_table",
        table_headers=headers,
        table_rows=rows,
        source_ids=source_ids,
        evidence_record_ids=evidence_ids,
    )
    # gene_symbol kept for API symmetry / future captions; table is evidence-driven.
    _ = gene_symbol
    return SectionPresentationResult(
        blocks=(block,),
        diagnostics=tuple(diagnostics),
    )


# --------------------------------------------------------------------------------------
# Species builders
# --------------------------------------------------------------------------------------
def _build_human(
    records: list[EvidenceRecord],
    diagnostics: list[PresentationDiagnostic],
) -> _SpeciesFields:
    out = _SpeciesFields()
    ncbi = [
        r
        for r in records
        if r.source_name == "NCBI Gene"
        and r.assertion_type == AssertionType.gene_identity
        and _is_human(r)
    ]
    ensembl = [
        r
        for r in records
        if r.source_name == "Ensembl"
        and r.assertion_type == AssertionType.gene_identity
        and r.fact_type == "ensembl_gene_id"
        and _is_human(r)
    ]
    uniprot = [
        r
        for r in records
        if r.source_name == "UniProt"
        and r.assertion_type == AssertionType.gene_identity
        and r.fact_type == "uniprot_accession"
        and _is_human(r)
    ]

    ncbi_alias_batches: list[tuple[EvidenceRecord, list[str]]] = []

    # Entrez / symbol / name: NCBI Gene → Ensembl/UniProt cross-ref if missing
    for rec in ncbi:
        _select_scalar(
            out,
            "human.entrez",
            "entrez",
            _cand_from_value(rec, "entrez_gene_id"),
            diagnostics,
        )
        _select_scalar(
            out,
            "human.symbol",
            "symbol",
            _cand_from_value(rec, "nomenclaturesymbol")
            or _cand_from_value(rec, "name"),
            diagnostics,
        )
        _select_scalar(
            out,
            "human.name",
            "name",
            _cand_from_value(rec, "nomenclaturename")
            or _cand_from_value(rec, "gene_name")
            or _cand_from_value(rec, "description"),
            diagnostics,
        )
        # Retain NCBI aliases for audit; polished cell prefers UniProt synonyms.
        aliases = _extract_alias_list(
            rec.value or {}, keys=("otheraliases", "synonyms", "aliases")
        )
        if aliases:
            ncbi_alias_batches.append((rec, aliases))

    for rec in ensembl:
        _select_scalar(
            out,
            "human.ensembl",
            "ensembl",
            _cand_from_value(rec, "ensembl_gene_id")
            or _cand_from_value(rec, "id"),
            diagnostics,
        )
        _select_scalar(
            out,
            "human.symbol",
            "symbol",
            _cand_from_value(rec, "display_name"),
            diagnostics,
        )

    uniprot_synonyms_used = False
    for rec in uniprot:
        _select_scalar(
            out,
            "human.uniprot",
            "uniprot",
            _cand_from_value(rec, "uniprot_accession")
            or _cand_from_value(rec, "primaryAccession"),
            diagnostics,
        )
        names = rec.value.get("gene_names") if isinstance(rec.value, dict) else None
        if isinstance(names, list) and names:
            str_names = [str(n).strip() for n in names if str(n).strip()]
            if out.symbol is None and str_names:
                _select_scalar(
                    out,
                    "human.symbol",
                    "symbol",
                    _Candidate(
                        display=str_names[0],
                        source_name=rec.source_name,
                        source_id=rec.source_id,
                        evidence_id=rec.id,
                    ),
                    diagnostics,
                )
        # Polished Human aliases: reviewed UniProt synonyms only when present.
        synonyms = _extract_alias_list(
            rec.value or {}, keys=("gene_synonyms",)
        )
        if not synonyms and isinstance(names, list):
            str_names = [str(n).strip() for n in names if str(n).strip()]
            symbol = out.symbol.display if out.symbol else None
            synonyms = [n for n in str_names if n != symbol]
        if synonyms:
            uniprot_synonyms_used = True
            out.aliases = []
            _merge_aliases(out, synonyms, rec)

    if not uniprot_synonyms_used:
        for rec, aliases in ncbi_alias_batches:
            _merge_aliases(out, aliases, rec)

    _note_ensembl_xref_conflicts(
        out, records, taxon=TAXON_HUMAN, field_prefix="human", diagnostics=diagnostics
    )
    return out


def _build_mouse(
    records: list[EvidenceRecord],
    diagnostics: list[PresentationDiagnostic],
) -> _SpeciesFields:
    out = _SpeciesFields()
    ncbi = [
        r
        for r in records
        if r.source_name == "NCBI Gene"
        and r.assertion_type == AssertionType.gene_identity
        and _taxon_of(r) == TAXON_MOUSE
    ]
    ensembl = [
        r
        for r in records
        if r.source_name == "Ensembl"
        and r.assertion_type == AssertionType.gene_identity
        and r.fact_type == "ensembl_gene_id"
        and _taxon_of(r) == TAXON_MOUSE
    ]
    orthologs = [
        r
        for r in records
        if r.fact_type == "ortholog_gene" and _taxon_of(r) == TAXON_MOUSE
    ]
    mgi = [r for r in records if r.fact_type == "mgi_gene_id"]

    # Primary species identity: NCBI Gene → Datasets → MGI
    for rec in ncbi:
        _select_scalar(
            out,
            "mouse.entrez",
            "entrez",
            _cand_from_value(rec, "entrez_gene_id"),
            diagnostics,
        )
        _select_scalar(
            out,
            "mouse.symbol",
            "symbol",
            _cand_from_value(rec, "nomenclaturesymbol")
            or _cand_from_value(rec, "gene_symbol")
            or _cand_from_value(rec, "name"),
            diagnostics,
        )
        _select_scalar(
            out,
            "mouse.name",
            "name",
            _cand_from_value(rec, "nomenclaturename")
            or _cand_from_value(rec, "gene_name")
            or _cand_from_value(rec, "description"),
            diagnostics,
        )
        aliases = _extract_alias_list(
            rec.value or {}, keys=("otheraliases", "aliases", "synonyms")
        )
        if aliases:
            _merge_aliases(out, aliases, rec)

    # Entrez / symbol: Datasets → MGI cross-check
    for rec in orthologs:
        _select_scalar(
            out,
            "mouse.entrez",
            "entrez",
            _cand_from_value(rec, "ortholog_gene_id"),
            diagnostics,
        )
        _select_scalar(
            out,
            "mouse.symbol",
            "symbol",
            _cand_from_value(rec, "ortholog_symbol"),
            diagnostics,
        )
        # name from datasets only if present
        _select_scalar(
            out,
            "mouse.name",
            "name",
            _cand_from_value(rec, "ortholog_name")
            or _cand_from_value(rec, "name"),
            diagnostics,
        )
        for key in ("ensembl_gene_id", "ensembl_id"):
            _select_scalar(
                out,
                "mouse.ensembl",
                "ensembl",
                _cand_from_value(rec, key),
                diagnostics,
            )
        for key in ("uniprot_accession", "uniprot_id"):
            _select_scalar(
                out,
                "mouse.uniprot",
                "uniprot",
                _cand_from_value(rec, key),
                diagnostics,
            )
    # Collect Datasets aliases to merge after MGI (species-specific order).
    datasets_alias_batches: list[tuple[EvidenceRecord, list[str]]] = []
    for rec in orthologs:
        aliases = _extract_alias_list(
            rec.value or {}, keys=("aliases", "synonyms", "otheraliases")
        )
        if aliases:
            datasets_alias_batches.append((rec, aliases))

    for rec in mgi:
        # Cross-check Entrez/symbol (lower priority than Datasets when already set)
        _select_scalar(
            out,
            "mouse.entrez",
            "entrez",
            _cand_from_value(rec, "ncbi_gene_number"),
            diagnostics,
        )
        _select_scalar(
            out,
            "mouse.symbol",
            "symbol",
            _cand_from_value(rec, "mouse_symbol"),
            diagnostics,
        )
        # Name: MGI preferred over datasets when NCBI Gene name is absent.
        mgi_name = _cand_from_value(rec, "mouse_name")
        if mgi_name is not None:
            if out.name is None:
                _select_scalar(out, "mouse.name", "name", mgi_name, diagnostics)
            elif _normalize_scalar("name", out.name.display) != _normalize_scalar(
                "name", mgi_name.display
            ):
                # Do not override NCBI Gene when already selected.
                if out.name.source_name == "NCBI Gene":
                    _track(out, mgi_name)
                else:
                    rejected = out.name
                    out.name = mgi_name
                    _track(out, mgi_name)
                    diagnostics.append(
                        PresentationDiagnostic(
                            field="mouse.name",
                            reason=(
                                f"selected {mgi_name.source_name}={mgi_name.display!r} "
                                f"(evidence_id={mgi_name.evidence_id}); "
                                f"rejected {rejected.source_name}={rejected.display!r} "
                                f"(evidence_id={rejected.evidence_id})"
                            ),
                            severity="warning",
                        )
                    )
            else:
                _track(out, mgi_name)

        aliases = _extract_alias_list(
            rec.value or {},
            keys=("aliases", "synonyms", "otheraliases", "mouse_aliases"),
        )
        _merge_aliases(out, aliases, rec)

    # Ortholog NCBI-style aliases after MGI
    for rec, aliases in datasets_alias_batches:
        _merge_aliases(out, aliases, rec)

    for rec in ensembl:
        _select_scalar(
            out,
            "mouse.ensembl",
            "ensembl",
            _cand_from_value(rec, "ensembl_gene_id")
            or _cand_from_value(rec, "id"),
            diagnostics,
        )

    # UniProt mouse ortholog evidence if present
    for rec in records:
        if (
            rec.source_name == "UniProt"
            and _taxon_of(rec) == TAXON_MOUSE
            and rec.fact_type == "uniprot_accession"
        ):
            _select_scalar(
                out,
                "mouse.uniprot",
                "uniprot",
                _cand_from_value(rec, "uniprot_accession"),
                diagnostics,
            )
            names = (rec.value or {}).get("gene_names")
            if isinstance(names, list):
                str_names = [str(n).strip() for n in names if str(n).strip()]
                symbol = out.symbol.display if out.symbol else None
                alias_part = [n for n in str_names if n != symbol]
                _merge_aliases(out, alias_part, rec)

    _note_ensembl_xref_conflicts(
        out, records, taxon=TAXON_MOUSE, field_prefix="mouse", diagnostics=diagnostics
    )
    return out


def _build_rat(
    records: list[EvidenceRecord],
    diagnostics: list[PresentationDiagnostic],
) -> _SpeciesFields:
    out = _SpeciesFields()
    ncbi = [
        r
        for r in records
        if r.source_name == "NCBI Gene"
        and r.assertion_type == AssertionType.gene_identity
        and _taxon_of(r) == TAXON_RAT
    ]
    ensembl = [
        r
        for r in records
        if r.source_name == "Ensembl"
        and r.assertion_type == AssertionType.gene_identity
        and r.fact_type == "ensembl_gene_id"
        and _taxon_of(r) == TAXON_RAT
    ]
    orthologs = [
        r
        for r in records
        if r.fact_type == "ortholog_gene" and _taxon_of(r) == TAXON_RAT
    ]

    for rec in ncbi:
        _select_scalar(
            out,
            "rat.entrez",
            "entrez",
            _cand_from_value(rec, "entrez_gene_id"),
            diagnostics,
        )
        _select_scalar(
            out,
            "rat.symbol",
            "symbol",
            _cand_from_value(rec, "nomenclaturesymbol")
            or _cand_from_value(rec, "gene_symbol")
            or _cand_from_value(rec, "name"),
            diagnostics,
        )
        _select_scalar(
            out,
            "rat.name",
            "name",
            _cand_from_value(rec, "nomenclaturename")
            or _cand_from_value(rec, "gene_name")
            or _cand_from_value(rec, "description"),
            diagnostics,
        )
        aliases = _extract_alias_list(
            rec.value or {}, keys=("otheraliases", "aliases", "synonyms")
        )
        if aliases:
            _merge_aliases(out, aliases, rec)

    for rec in orthologs:
        _select_scalar(
            out,
            "rat.entrez",
            "entrez",
            _cand_from_value(rec, "ortholog_gene_id"),
            diagnostics,
        )
        _select_scalar(
            out,
            "rat.symbol",
            "symbol",
            _cand_from_value(rec, "ortholog_symbol"),
            diagnostics,
        )
        _select_scalar(
            out,
            "rat.name",
            "name",
            _cand_from_value(rec, "ortholog_name")
            or _cand_from_value(rec, "name"),
            diagnostics,
        )
        for key in ("ensembl_gene_id", "ensembl_id"):
            _select_scalar(
                out,
                "rat.ensembl",
                "ensembl",
                _cand_from_value(rec, key),
                diagnostics,
            )
        for key in ("uniprot_accession", "uniprot_id"):
            _select_scalar(
                out,
                "rat.uniprot",
                "uniprot",
                _cand_from_value(rec, key),
                diagnostics,
            )
        aliases = _extract_alias_list(
            rec.value or {}, keys=("aliases", "synonyms", "otheraliases")
        )
        _merge_aliases(out, aliases, rec)

    for rec in ensembl:
        _select_scalar(
            out,
            "rat.ensembl",
            "ensembl",
            _cand_from_value(rec, "ensembl_gene_id")
            or _cand_from_value(rec, "id"),
            diagnostics,
        )

    for rec in records:
        if (
            rec.source_name == "UniProt"
            and _taxon_of(rec) == TAXON_RAT
            and rec.fact_type == "uniprot_accession"
        ):
            _select_scalar(
                out,
                "rat.uniprot",
                "uniprot",
                _cand_from_value(rec, "uniprot_accession"),
                diagnostics,
            )
            names = (rec.value or {}).get("gene_names")
            if isinstance(names, list):
                str_names = [str(n).strip() for n in names if str(n).strip()]
                symbol = out.symbol.display if out.symbol else None
                _merge_aliases(out, [n for n in str_names if n != symbol], rec)

    _note_ensembl_xref_conflicts(
        out, records, taxon=TAXON_RAT, field_prefix="rat", diagnostics=diagnostics
    )
    return out


def _note_ensembl_xref_conflicts(
    species: _SpeciesFields,
    records: list[EvidenceRecord],
    *,
    taxon: int,
    field_prefix: str,
    diagnostics: list[PresentationDiagnostic],
) -> None:
    """Warn when UniProt Ensembl xrefs disagree with direct Ensembl lookup."""
    direct = species.ensembl
    if direct is None:
        return
    direct_norm = _normalize_scalar("ensembl", direct.display)
    for rec in records:
        if (
            rec.source_name != "UniProt"
            or rec.fact_type != "ensembl_xref"
            or _taxon_of(rec) != taxon
        ):
            continue
        xref = _cand_from_value(rec, "ensembl_gene_id")
        if xref is None:
            continue
        xref_norm = _normalize_scalar("ensembl", xref.display)
        if xref_norm and direct_norm and xref_norm != direct_norm:
            diagnostics.append(
                PresentationDiagnostic(
                    field=f"{field_prefix}.ensembl",
                    reason=(
                        f"kept direct Ensembl={direct.display!r}; "
                        f"preserved secondary UniProt xref={xref.display!r} "
                        f"(evidence_id={xref.evidence_id})"
                    ),
                    severity="warning",
                )
            )


# --------------------------------------------------------------------------------------
# Selection / merge helpers
# --------------------------------------------------------------------------------------
def _select_scalar(
    species: _SpeciesFields,
    field_path: str,
    attr: str,
    candidate: _Candidate | None,
    diagnostics: list[PresentationDiagnostic],
) -> None:
    if candidate is None or not str(candidate.display).strip():
        return
    current: _Candidate | None = getattr(species, attr)
    if current is None:
        setattr(species, attr, candidate)
        _track(species, candidate)
        return

    cur_norm = _normalize_scalar(attr, current.display)
    new_norm = _normalize_scalar(attr, candidate.display)
    if cur_norm == new_norm:
        _track(species, candidate)
        return

    # Higher-priority already set; reject lower-priority with warning
    diagnostics.append(
        PresentationDiagnostic(
            field=field_path,
            reason=(
                f"selected {current.source_name}={current.display!r} "
                f"(evidence_id={current.evidence_id}); "
                f"rejected {candidate.source_name}={candidate.display!r} "
                f"(evidence_id={candidate.evidence_id})"
            ),
            severity="warning",
        )
    )


def _merge_aliases(
    species: _SpeciesFields,
    aliases: Sequence[str],
    rec: EvidenceRecord,
) -> None:
    if not aliases:
        return
    symbol = species.symbol.display if species.symbol else None
    merged = _dedupe_aliases(species.aliases + list(aliases), official_symbol=symbol)
    if merged != species.aliases:
        species.aliases = merged
        _track(
            species,
            _Candidate(
                display=",".join(aliases),
                source_name=rec.source_name,
                source_id=rec.source_id,
                evidence_id=rec.id,
            ),
        )


def _dedupe_aliases(
    aliases: Sequence[str],
    *,
    official_symbol: str | None,
) -> list[str]:
    out: list[str] = []
    seen_exact: set[str] = set()
    seen_ci: set[str] = set()
    official_ci = (official_symbol or "").strip().lower()
    for raw in aliases:
        text = str(raw).strip()
        if not text:
            continue
        if official_ci and text.lower() == official_ci:
            continue
        if text in seen_exact:
            continue
        if text.lower() in seen_ci:
            continue
        seen_exact.add(text)
        seen_ci.add(text.lower())
        out.append(text)
    return out


def _emit_missing(
    species_key: str,
    species: _SpeciesFields,
    diagnostics: list[PresentationDiagnostic],
) -> None:
    mapping = {
        "entrez": species.entrez,
        "symbol": species.symbol,
        "name": species.name,
        "ensembl": species.ensembl,
        "uniprot": species.uniprot,
    }
    for field_name, value in mapping.items():
        if value is None:
            diagnostics.append(
                PresentationDiagnostic(
                    field=f"{species_key}.{field_name}",
                    reason=(
                        f"no structured {field_name} evidence for {species_key} "
                        "in the supplied EvidenceRecords"
                    ),
                    severity="info",
                )
            )
    if not species.aliases:
        diagnostics.append(
            PresentationDiagnostic(
                field=f"{species_key}.aliases",
                reason=(
                    f"no structured aliases for {species_key} "
                    "in the supplied EvidenceRecords"
                ),
                severity="info",
            )
        )


# --------------------------------------------------------------------------------------
# Normalization / URL helpers
# --------------------------------------------------------------------------------------
def _normalize_scalar(kind: str, value: str) -> str:
    text = (value or "").strip()
    if kind == "entrez":
        digits = re.sub(r"\D", "", text)
        return digits.lstrip("0") or digits
    if kind == "ensembl":
        match = _ENSEMBL_BASE_RE.match(text)
        return (match.group(1) if match else text).upper()
    if kind == "uniprot":
        # Isoform accessions stay distinct from canonical (do not strip -N)
        return text.upper()
    if kind == "name":
        return re.sub(r"\s+", " ", text).casefold()
    if kind == "symbol":
        return text  # preserve case; whitespace already trimmed
    return text


def _cand_from_value(rec: EvidenceRecord, key: str) -> _Candidate | None:
    value = rec.value or {}
    if not isinstance(value, dict):
        return None
    raw = value.get(key)
    if raw is None:
        return None
    text = str(raw).strip()
    if not text:
        return None
    return _Candidate(
        display=text,
        source_name=rec.source_name,
        source_id=rec.source_id,
        evidence_id=rec.id,
    )


def _extract_alias_list(value: dict[str, Any], *, keys: Sequence[str]) -> list[str]:
    out: list[str] = []
    for key in keys:
        raw = value.get(key)
        if raw is None:
            continue
        if isinstance(raw, list):
            out.extend(str(x).strip() for x in raw if str(x).strip())
        elif isinstance(raw, str):
            # split common separators
            parts = re.split(r"[,;|]", raw)
            out.extend(p.strip() for p in parts if p.strip())
    return out


def _track(species: _SpeciesFields, cand: _Candidate) -> None:
    if cand.source_id and cand.source_id not in species.source_ids:
        species.source_ids.append(cand.source_id)
    if cand.evidence_id and cand.evidence_id not in species.evidence_ids:
        species.evidence_ids.append(cand.evidence_id)


def _taxon_of(rec: EvidenceRecord) -> int | None:
    if rec.taxon_id is not None:
        try:
            return int(rec.taxon_id)
        except (TypeError, ValueError):
            pass
    value = rec.value or {}
    if isinstance(value, dict):
        for key in ("tax_id", "taxon_id", "organism_id"):
            raw = value.get(key)
            if raw is None:
                continue
            try:
                return int(raw)
            except (TypeError, ValueError):
                continue
    organism = (rec.organism or rec.species or "").lower().replace("_", " ")
    if "homo sapiens" in organism or organism == "human":
        return TAXON_HUMAN
    if "mus musculus" in organism or organism == "mouse":
        return TAXON_MOUSE
    if "rattus norvegicus" in organism or "norway rat" in organism or organism == "rat":
        return TAXON_RAT
    return None


def _is_human(rec: EvidenceRecord) -> bool:
    tax = _taxon_of(rec)
    return tax == TAXON_HUMAN


def _plain(cand: _Candidate | None) -> str:
    return cand.display if cand else NOT_AVAILABLE


def _aliases_cell(aliases: list[str]) -> str:
    if not aliases:
        return NOT_AVAILABLE
    return ", ".join(aliases)


def _linked_entrez(cand: _Candidate | None) -> str:
    if cand is None:
        return NOT_AVAILABLE
    text = cand.display.strip()
    if not _ENTREZ_RE.match(text):
        return text
    return f"[{text}](https://www.ncbi.nlm.nih.gov/gene/{text})"


def _linked_ensembl(cand: _Candidate | None) -> str:
    if cand is None:
        return NOT_AVAILABLE
    text = cand.display.strip()
    if not _ENSEMBL_RE.match(text):
        return text
    # Prefer unversioned ID in Gene Aliases table display/link target
    base = _ENSEMBL_BASE_RE.match(text)
    display = base.group(1) if base else text
    return f"[{display}](https://www.ensembl.org/id/{display})"


def _linked_uniprot(cand: _Candidate | None) -> str:
    if cand is None:
        return NOT_AVAILABLE
    text = cand.display.strip()
    upper = text.upper()
    if _UNIPROT_ISOFORM_RE.match(upper):
        # Link isoform pages with isoform accession
        return f"[{text}](https://www.uniprot.org/uniprotkb/{upper})"
    if not _UNIPROT_RE.match(upper):
        return text
    return f"[{text}](https://www.uniprot.org/uniprotkb/{upper})"


def format_safe_table_cell_html(cell: str) -> str:
    """Render a table cell: whole-cell canonical https link or escaped text.

    Exported for renderer/tests. Not a general Markdown parser.
    """
    from html import escape

    text = cell if cell is not None else ""
    match = _CELL_LINK_RE.fullmatch(text.strip())
    if not match:
        return escape(text)

    label, url = match.group(1), match.group(2)
    if not _is_allowed_link(url):
        return escape(text)

    # Classify link kind for CSS (entrez/ensembl vs uniprot)
    host = _link_host(url)
    css = "id-link"
    if host == "www.uniprot.org":
        css = "id-link id-link-uniprot"
    elif host == "www.ncbi.nlm.nih.gov":
        css = "id-link id-link-entrez"
    elif host == "www.ensembl.org":
        css = "id-link id-link-ensembl"

    return (
        f'<a class="{css}" href="{escape(url, quote=True)}">'
        f"{escape(label)}</a>"
    )


def _link_host(url: str) -> str:
    from urllib.parse import urlparse

    try:
        return (urlparse(url).hostname or "").lower()
    except Exception:  # noqa: BLE001
        return ""


def _is_allowed_link(url: str) -> bool:
    from urllib.parse import urlparse

    if any(ch in url for ch in ('"', "'", "<", ">", " ")):
        return False
    try:
        parsed = urlparse(url)
    except Exception:  # noqa: BLE001
        return False
    if parsed.scheme != "https":
        return False
    host = (parsed.hostname or "").lower()
    if host not in ALLOWED_LINK_HOSTS:
        return False
    # Reject unexpected userinfo / credentials
    if parsed.username or parsed.password:
        return False
    return True


def _ucsc_records(records: Sequence[EvidenceRecord]) -> list[EvidenceRecord]:
    return [r for r in records if (r.source_name or "") == "UCSC"]


def _first_fact(records: Sequence[EvidenceRecord], fact_type: str) -> EvidenceRecord | None:
    matches = [r for r in records if r.fact_type == fact_type]
    if not matches:
        return None
    if fact_type != "ucsc_conservation_figure":
        return matches[0]
    # Prefer a figure whose managed bytes still validate; then newest created_at.
    from gene_dossier.ucsc_figure import resolve_artifact_path, sha256_hex

    validated: list[EvidenceRecord] = []
    for rec in matches:
        value = rec.value if isinstance(rec.value, dict) else {}
        rel = value.get("relative_path") or value.get("local_artifact_path")
        expected = value.get("sha256") or value.get("content_hash")
        if not rel or not expected:
            continue
        try:
            path = resolve_artifact_path(str(rel))
        except ValueError:
            continue
        if path.is_file() and sha256_hex(path.read_bytes()) == expected:
            validated.append(rec)
    pool = validated or matches
    return sorted(
        pool,
        key=lambda r: getattr(r, "created_at", None) or 0,
        reverse=True,
    )[0]


def _human_gene_name(records: Sequence[EvidenceRecord]) -> str | None:
    """Prefer structured human identity gene name for the transcript line."""
    for preferred_sources in (("NCBI Gene",), ("UniProt",), ("Ensembl",)):
        for rec in records:
            if rec.source_name not in preferred_sources:
                continue
            if rec.taxon_id not in {None, 9606}:
                org = (rec.organism or "").lower()
                if org and org not in {"homo sapiens", "human"}:
                    continue
            value = rec.value if isinstance(rec.value, dict) else {}
            name = (
                value.get("nomenclaturename")
                or value.get("gene_name")
                or value.get("official_full_name")
                or value.get("description")
            )
            if isinstance(name, str) and name.strip() and name.strip().upper() != (
                rec.gene_symbol or ""
            ).upper():
                # Prefer long names over symbols.
                if " " in name.strip() or len(name.strip()) > 12:
                    return name.strip()
            if isinstance(name, str) and name.strip():
                candidate = name.strip()
                if candidate.upper() != (rec.official_symbol or rec.gene_symbol or "").upper():
                    return candidate
    return None


def _resolve_figure_path(value: dict[str, Any]) -> tuple[str | None, list[PresentationDiagnostic]]:
    from gene_dossier.ucsc_figure import resolve_artifact_path, sha256_hex

    diags: list[PresentationDiagnostic] = []
    rel = value.get("relative_path") or value.get("local_artifact_path") or value.get("figure_path")
    if not rel:
        diags.append(
            PresentationDiagnostic("figure", "missing relative figure path", "warning")
        )
        return None, diags
    try:
        path = resolve_artifact_path(str(rel))
    except ValueError as exc:
        diags.append(PresentationDiagnostic("figure", str(exc), "warning"))
        return None, diags
    if not path.is_file():
        diags.append(
            PresentationDiagnostic("figure", f"managed figure missing: {rel}", "warning")
        )
        return None, diags
    expected = value.get("sha256") or value.get("content_hash")
    if expected:
        actual = sha256_hex(path.read_bytes())
        if actual != expected:
            diags.append(
                PresentationDiagnostic(
                    "figure",
                    "figure checksum mismatch; omitting image",
                    "warning",
                )
            )
            return None, diags
    # Persist portable relative path into the block (never home absolute).
    portable = value.get("relative_path") or str(rel)
    if portable.startswith("/Users/") or portable.startswith("C:\\") or portable.startswith("C:/"):
        diags.append(
            PresentationDiagnostic(
                "figure",
                "refusing machine-specific absolute figure path",
                "warning",
            )
        )
        return None, diags
    return str(portable), diags


def _transcript_selection_sentence(tx_val: dict[str, Any]) -> str:
    """Return flag-accurate transcript selection wording for Section 1b."""
    is_mane = bool(tx_val.get("is_mane_select"))
    is_ensembl = bool(tx_val.get("is_ensembl_canonical"))
    is_gencode_primary = bool(tx_val.get("is_gencode_primary"))
    is_canonical_tier = bool(tx_val.get("is_canonical_tier"))
    if is_mane and is_ensembl:
        return (
            "The transcript selected for display is both MANE Select and "
            "Ensembl canonical."
        )
    if is_mane:
        return "The MANE Select transcript was selected for display."
    if is_ensembl:
        return "The Ensembl canonical transcript was selected for display."
    if is_gencode_primary or is_canonical_tier:
        return (
            "The current canonical-tier transcript model was selected for display."
        )
    return "The highest-ranked current transcript model was selected for display."


def build_conservation_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
) -> SectionPresentationResult:
    """Build Section 1b polished UCSC conservation blocks."""
    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    ucsc = _ucsc_records(records)

    locus = _first_fact(ucsc, "ucsc_gene_locus")
    inventory = _first_fact(ucsc, "ucsc_transcript_inventory")
    transcript = _first_fact(ucsc, "ucsc_canonical_transcript")
    figure = _first_fact(ucsc, "ucsc_conservation_figure")

    if inventory is None:
        diagnostics.append(
            PresentationDiagnostic("search_summary", "missing transcript inventory", "warning")
        )
    if transcript is None:
        diagnostics.append(
            PresentationDiagnostic("canonical_transcript", "missing canonical transcript", "warning")
        )
    if figure is None:
        diagnostics.append(
            PresentationDiagnostic("figure", "missing conservation figure", "warning")
        )

    inv_val = inventory.value if inventory and isinstance(inventory.value, dict) else {}
    loc_val = locus.value if locus and isinstance(locus.value, dict) else {}
    tx_val = transcript.value if transcript and isinstance(transcript.value, dict) else {}
    fig_val = figure.value if figure and isinstance(figure.value, dict) else {}

    assembly = (
        loc_val.get("assembly_display_name")
        or tx_val.get("assembly_display_name")
        or "GRCh38/hg38"
    )
    exact_count = inv_val.get("exact_gene_transcript_count")
    release = inv_val.get("current_gencode_release")
    gene = gene_symbol.strip() or str(tx_val.get("requested_gene_symbol") or "")

    if exact_count is None:
        diagnostics.append(
            PresentationDiagnostic("exact_gene_count", "missing exact-gene transcript count", "warning")
        )
    if not release:
        diagnostics.append(
            PresentationDiagnostic("gencode_release", "missing current GENCODE release", "warning")
        )

    count_phrase = (
        str(exact_count)
        if exact_count is not None
        else "an unknown number of"
    )
    if release and str(release).upper().startswith("GENCODE"):
        release_phrase = f"current {release} annotation"
    elif release and str(release).upper().startswith("V"):
        release_phrase = f"current GENCODE {release} annotation"
    elif release:
        release_phrase = f"current GENCODE {release} annotation"
    else:
        release_phrase = "current GENCODE annotation"

    selection = _transcript_selection_sentence(tx_val)
    dynamic = (
        f"A query of the {assembly} assembly identified {count_phrase} "
        f"{gene} transcript models in the {release_phrase} within the selected "
        f"locus. {selection}"
    )
    combined_narrative = f"{UCSC_STABLE_INTRO} {dynamic}"

    source_ids = list(
        dict.fromkeys(
            sid
            for rec in (locus, inventory, transcript, figure)
            if rec is not None and rec.source_id
            for sid in [rec.source_id]
        )
    )
    evidence_ids = list(
        dict.fromkeys(
            eid
            for rec in (locus, inventory, transcript, figure)
            if rec is not None and rec.id
            for eid in [rec.id]
        )
    )

    blocks: list[ReportContentBlock] = [
        ReportContentBlock(
            kind="narrative",
            text=combined_narrative,
            presentation_role=None,
            source_ids=source_ids,
            evidence_record_ids=evidence_ids,
        ),
    ]

    tx_id = str(tx_val.get("transcript_id") or "").strip()
    display_pos = str(
        tx_val.get("display_position")
        or loc_val.get("display_position")
        or ""
    ).strip()
    gene_name = _human_gene_name(records) or str(
        tx_val.get("gene_name")
        or tx_val.get("source_gene_symbol")
        or gene
    )
    # Prefer full name from identity when available
    full_name = _human_gene_name(records)
    if full_name and full_name.upper() != gene.upper():
        gene_name = full_name

    label = f"{gene} ({tx_id}) - {display_pos} - Homo sapiens {gene_name}"
    from gene_dossier.ucsc_figure import build_safe_hgtracks_url, is_safe_ucsc_browser_url

    url = tx_val.get("browser_url")
    if not url or not is_safe_ucsc_browser_url(str(url)):
        url = build_safe_hgtracks_url(
            genome=str(tx_val.get("genome") or loc_val.get("genome") or "hg38"),
            display_position=display_pos,
            transcript_id=tx_id or None,
        )
    if url and is_safe_ucsc_browser_url(str(url)) and tx_id and display_pos:
        blocks.append(
            ReportContentBlock(
                kind="link",
                text=label,
                links=[{"label": label, "url": str(url)}],
                source_ids=source_ids,
                evidence_record_ids=evidence_ids,
            )
        )
    else:
        # Invalid transcript id / URL → escaped plain narrative (no unsafe link).
        diagnostics.append(
            PresentationDiagnostic(
                "transcript_link",
                "UCSC transcript link omitted or rendered as plain text",
                "info",
            )
        )
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=label,
                source_ids=source_ids,
                evidence_record_ids=evidence_ids,
            )
        )

    if figure is not None and fig_val:
        resolved, fig_diags = _resolve_figure_path(fig_val)
        diagnostics.extend(fig_diags)
        if resolved:
            alt = (
                f"UCSC Genome Browser conservation tracks for {gene} at "
                f"{display_pos or fig_val.get('display_position') or 'selected locus'}"
            )
            blocks.append(
                ReportContentBlock(
                    kind="figure",
                    figure_path=resolved,
                    figure_caption=alt,
                    text=alt,
                    presentation_role="ucsc_conservation_figure",
                    source_ids=source_ids,
                    evidence_record_ids=evidence_ids,
                )
            )
            # source_note / caption are provenance — keep in diagnostics, not
            # as a second visible narrative block.
            note = fig_val.get("caption") or fig_val.get("source_note")
            if note:
                diagnostics.append(
                    PresentationDiagnostic(
                        "figure_note",
                        str(note),
                        "info",
                    )
                )

    # Need at least intro + dynamic to consider presentation nonempty for unknown genes
    if locus is None and inventory is None and transcript is None:
        return SectionPresentationResult(blocks=(), diagnostics=tuple(diagnostics))

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


def _fmt_percent(value: Any) -> str:
    try:
        if value is None:
            return NOT_AVAILABLE
        return f"{float(value) * 100:.1f}%"
    except (TypeError, ValueError):
        return NOT_AVAILABLE


def _fmt_span(start: Any, end: Any) -> str:
    if start is None or end is None:
        return NOT_AVAILABLE
    return f"{start}-{end}"


def _fmt_resolution(value: Any) -> str:
    try:
        if value is None:
            return NOT_AVAILABLE
        return f"{float(value):.2g} Angstrom"
    except (TypeError, ValueError):
        return NOT_AVAILABLE


def _section_1c_records(records: Sequence[EvidenceRecord]) -> list[EvidenceRecord]:
    return [
        rec
        for rec in records
        if rec.source_name in {"CDD", "PDBe", "UniProt"}
        and (
            rec.fact_type
            in {
                "conserved_domain_hit",
                "cdd_architecture_figure",
                "cdd_official_architecture_figure",
                "cdd_family_summary",
                "cdd_family_thumbnail",
                "cdd_conserved_feature",
                "cdd_feature_thumbnail",
                "pdb_structure",
                "pdb_candidate_selection",
                "pdb_coordinate_mapping",
                "pdb_assembly_figure",
                "pdb_domain_focus_figure",
                "pdb_official_structure_image",
                "sequence_feature",
                "uniprot_accession",
            }
        )
    ]


def _protein_length_from_records(records: Sequence[EvidenceRecord]) -> int | None:
    for rec in records:
        if rec.source_name != "UniProt" or rec.fact_type != "uniprot_accession":
            continue
        value = rec.value if isinstance(rec.value, dict) else {}
        for key in ("protein_length", "sequence_length", "length"):
            try:
                raw = value.get(key)
                if raw is not None:
                    return int(raw)
            except (TypeError, ValueError):
                continue
    return None


def _fact_records(records: Sequence[EvidenceRecord], fact_type: str) -> list[EvidenceRecord]:
    return [rec for rec in records if rec.fact_type == fact_type]


def _first_section_1c_fact(
    records: Sequence[EvidenceRecord],
    fact_type: str,
) -> EvidenceRecord | None:
    matches = _fact_records(records, fact_type)
    if not matches:
        return None
    return sorted(
        matches,
        key=lambda r: getattr(r, "created_at", None) or 0,
        reverse=True,
    )[0]


def _candidate_selection_value(records: Sequence[EvidenceRecord]) -> dict[str, Any]:
    rec = _first_section_1c_fact(records, "pdb_candidate_selection")
    return rec.value if rec is not None and isinstance(rec.value, dict) else {}


def _selected_candidate(selection: dict[str, Any]) -> dict[str, Any] | None:
    candidates = selection.get("candidates")
    if not isinstance(candidates, list):
        return None
    for candidate in candidates:
        if isinstance(candidate, dict) and candidate.get("selected"):
            return candidate
    for candidate in candidates:
        if isinstance(candidate, dict):
            return candidate
    return None


def _candidate_rows(selection: dict[str, Any]) -> list[list[str]]:
    candidates = selection.get("candidates")
    if not isinstance(candidates, list):
        return []
    rows: list[list[str]] = []
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        chains = ", ".join(str(x) for x in candidate.get("chain_ids") or [])
        spans = ", ".join(
            f"{span[0]}-{span[1]}"
            for span in candidate.get("mapped_spans") or []
            if isinstance(span, list) and len(span) == 2
        )
        status = "Selected" if candidate.get("selected") else "Rejected"
        reasons = "; ".join(str(x) for x in candidate.get("rejection_reasons") or [])
        if not reasons and not candidate.get("selected"):
            reasons = "lower ranked than selected candidate"
        rows.append(
            [
                str(candidate.get("pdb_id") or NOT_AVAILABLE).upper(),
                chains or NOT_AVAILABLE,
                str(candidate.get("experimental_method") or NOT_AVAILABLE),
                _fmt_resolution(candidate.get("resolution")),
                spans or NOT_AVAILABLE,
                _fmt_percent(candidate.get("calculated_coverage")),
                status if not reasons else f"{status}: {reasons}",
            ]
        )
    return rows


def _cdd_accession_url(accession: Any) -> str | None:
    acc = str(accession or "").strip()
    if not acc:
        return None
    return f"https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={acc}"


def _pdbe_entry_url(pdb_id: Any) -> str | None:
    pdb = str(pdb_id or "").strip().lower()
    if not pdb:
        return None
    return f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb}"


def _best_number(values: Sequence[Any], *, lowest: bool) -> str:
    parsed: list[float] = []
    for value in values:
        try:
            if value is not None and str(value).strip():
                parsed.append(float(str(value)))
        except (TypeError, ValueError):
            continue
    if not parsed:
        return NOT_AVAILABLE
    number = min(parsed) if lowest else max(parsed)
    return f"{number:.4g}"


def _domain_groups(domains: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], dict[str, Any]] = {}
    for row in domains:
        name = str(row.get("domain_short_name") or row.get("domain_accession") or "CDD domain")
        accession = str(row.get("domain_accession") or "")
        key = (name, accession)
        group = grouped.setdefault(
            key,
            {
                "domain_short_name": name,
                "domain_accession": accession,
                "domain_description": row.get("domain_description"),
                "rows": [],
            },
        )
        group["rows"].append(row)
    return sorted(
        grouped.values(),
        key=lambda item: min(
            (
                row.get("from_residue")
                for row in item["rows"]
                if row.get("from_residue") is not None
            ),
            default=10**9,
        ),
    )


def _domain_summary_text(
    *,
    gene_symbol: str,
    group: dict[str, Any],
) -> str:
    rows = [row for row in group.get("rows") or [] if isinstance(row, dict)]
    spans = ", ".join(
        _fmt_span(row.get("from_residue"), row.get("to_residue")) for row in rows
    )
    coverages = [
        _fmt_percent(row.get("coverage"))
        for row in rows
        if _fmt_percent(row.get("coverage")) != NOT_AVAILABLE
    ]
    coverage_text = ", ".join(coverages)
    name = str(group.get("domain_short_name") or "CDD domain")
    accession = str(group.get("domain_accession") or NOT_AVAILABLE)
    occurrence = "match" if len(rows) == 1 else f"{len(rows)} matches"
    text = (
        f"{name} ({accession}) has {occurrence} in {gene_symbol}"
        + (f" at residues {spans}" if spans else "")
        + "."
    )
    if coverage_text:
        text += f" Calculated coverage: {coverage_text}."
    text += (
        f" Best E-value: {_best_number([row.get('evalue') for row in rows], lowest=True)}; "
        f"best bit score: {_best_number([row.get('bitscore') for row in rows], lowest=False)}."
    )
    description = group.get("domain_description")
    if description:
        text += f" {description}"
    return text


def _source_and_evidence_ids(records: Sequence[EvidenceRecord]) -> tuple[list[str], list[str]]:
    return (
        list(dict.fromkeys(rec.source_id for rec in records if rec.source_id)),
        list(dict.fromkeys(rec.id for rec in records if rec.id)),
    )


def _safe_item_token(value: Any, *, fallback: str) -> str:
    token = re.sub(r"[^a-z0-9]+", "-", str(value or fallback).lower()).strip("-")
    return token or fallback


def _value_item_key(value: dict[str, Any], *, prefix: str, fallback: str) -> str:
    if value.get("presentation_item_key"):
        return str(value["presentation_item_key"])
    raw = (
        value.get("canonical_accession")
        or value.get("domain_accession")
        or value.get("pdb_id")
        or value.get("feature_label")
        or fallback
    )
    return f"{prefix}-{_safe_item_token(raw, fallback=fallback)}"


def _figure_block_from_record(
    rec: EvidenceRecord,
    *,
    role: str,
    caption: str,
    diagnostics: list[PresentationDiagnostic],
) -> ReportContentBlock | None:
    value = rec.value if isinstance(rec.value, dict) else {}
    resolved, fig_diags = _resolve_figure_path(value)
    diagnostics.extend(fig_diags)
    if not resolved:
        return None
    return ReportContentBlock(
        kind="figure",
        figure_path=resolved,
        figure_caption=caption,
        text=caption,
        presentation_role=role,  # type: ignore[arg-type]
        presentation_item_key=(
            str(value.get("presentation_item_key"))
            if isinstance(value, dict) and value.get("presentation_item_key")
            else None
        ),
        source_ids=[rec.source_id] if rec.source_id else [],
        evidence_record_ids=[rec.id] if rec.id else [],
    )


# Trailing CDD link wording used by the historical Rancho domain paragraphs.
SECTION_1C_CDD_LINK_LABEL = "(NCBI CDD Link)"

# Rancho body style sets the specific-domain lead phrase in bold and continues
# the sentence in normal weight. Renderers bold these phrases when a Section 1c
# domain summary opens with one; the stored block text stays plain.
SECTION_1C_BOLD_LEAD_PHRASES: tuple[str, ...] = (
    "basic Helix-Loop-Helix-zipper (bHLHzip) domain",
)

# Domain families whose polished block opens a new Section 1c page, matching the
# page-per-topic body flow of the original reports (C-terminal cytoplasmic
# region follows the extracellular repeat content on its own page).
_SECTION_1C_PAGE_LEADING_ACCESSIONS = frozenset({"pfam01049"})

_SECTION_1C_PDB_PAGE_ROLES = frozenset(
    {"section_1c_pdb_link", "section_1c_pdb_official_image"}
)


def _strip_visible_pssm(text: str) -> str:
    cleaned = re.sub(r"\s*\(PSSM\s*ID:\s*\d+\)", "", text, flags=re.IGNORECASE)
    cleaned = re.sub(r"\bPSSM[-\s]*ID:\s*\d+\b[:;,\s]*", "", cleaned, flags=re.IGNORECASE)
    return re.sub(r"\s+", " ", cleaned).strip()


def _record_accession(rec: EvidenceRecord) -> str:
    value = rec.value if isinstance(rec.value, dict) else {}
    return str(value.get("canonical_accession") or value.get("domain_accession") or "").lower()


def _family_record_for_accession(
    family_summaries: Sequence[EvidenceRecord],
    accession: str,
) -> EvidenceRecord | None:
    target = accession.lower()
    candidates: list[EvidenceRecord] = []
    for rec in family_summaries:
        value = rec.value if isinstance(rec.value, dict) else {}
        aliases = {
            _record_accession(rec),
            str(value.get("matched_query_domain_accession") or "").lower(),
            str(value.get("domain_accession") or "").lower(),
            str(value.get("canonical_accession") or "").lower(),
        }
        if target in aliases:
            candidates.append(rec)
    if not candidates:
        return None
    if target == "pfam01049":
        for rec in candidates:
            value = rec.value if isinstance(rec.value, dict) else {}
            if str(value.get("pssm_id") or "") == "426014":
                return rec
            haystack = " ".join(
                str(value.get(key) or "")
                for key in ("domain_short_name", "canonical_accession", "synopsis")
            ).lower()
            if "cadherin_c" in haystack or "cadherin cytoplasmic" in haystack:
                return rec
    return candidates[0]
    return None


def _thumbnail_records_for_key(
    thumbnails: Sequence[EvidenceRecord],
    item_key: str,
) -> list[EvidenceRecord]:
    return [
        rec
        for rec in thumbnails
        if isinstance(rec.value, dict)
        and str(rec.value.get("presentation_item_key") or "") == item_key
    ]


# NCBI CDD serves thin alignment/sequence strips (e.g. 100x24) from the same
# thumbnail endpoints as real Cn3D structure images, sometimes under a
# ``*_structure_thumbnail`` role. True structure thumbnails are square-ish
# (100x100, 300x300), so height and aspect ratio separate the two reliably.
_MIN_STRUCTURE_THUMBNAIL_HEIGHT = 40
_MAX_STRUCTURE_THUMBNAIL_ASPECT = 2.5
_NON_STRUCTURE_THUMBNAIL_ROLE_TOKENS = ("alignment", "sequence", "logo", "msa")
_CDD_UID_RE = re.compile(r"uid=(\d+)")


def _cdd_uid(value: dict[str, Any]) -> str:
    """CDD page uid for a family summary or thumbnail record, or ``""``."""
    for key in ("requested_uid", "pssm_id"):
        raw = str(value.get(key) or "").strip()
        if raw:
            return raw
    match = _CDD_UID_RE.search(str(value.get("source_url") or ""))
    return match.group(1) if match else ""


def _is_renderable_structure_thumbnail(
    rec: EvidenceRecord,
    *,
    family_uid: str = "",
) -> bool:
    """True when a CDD thumbnail is a real structure image for ``family_uid``.

    Rejects alignment/sequence/logo/MSA roles, thin strips, and thumbnails whose
    CDD uid differs from the family page whose synopsis is being rendered (pfam
    accessions can resolve to more than one CDD family page).
    """
    value = rec.value if isinstance(rec.value, dict) else {}
    role = str(value.get("classified_role") or "").strip().lower()
    if any(token in role for token in _NON_STRUCTURE_THUMBNAIL_ROLE_TOKENS):
        return False
    thumbnail_uid = _cdd_uid(value)
    if family_uid and thumbnail_uid and thumbnail_uid != family_uid:
        return False
    try:
        width = int(value.get("width"))  # type: ignore[arg-type]
        height = int(value.get("height"))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        # Dimensions unavailable: the role and uid checks are all we can apply.
        return True
    if width <= 0 or height <= 0:
        return True
    if height <= _MIN_STRUCTURE_THUMBNAIL_HEIGHT:
        return False
    return (width / height) < _MAX_STRUCTURE_THUMBNAIL_ASPECT


def _renderable_thumbnail_records(
    thumbnails: Sequence[EvidenceRecord],
    *,
    item_key: str,
    family_uid: str,
    diagnostics: list[PresentationDiagnostic],
) -> list[EvidenceRecord]:
    """Structure thumbnails for ``item_key``; rejected ones become diagnostics."""
    kept: list[EvidenceRecord] = []
    for rec in _thumbnail_records_for_key(thumbnails, item_key):
        if _is_renderable_structure_thumbnail(rec, family_uid=family_uid):
            kept.append(rec)
            continue
        diagnostics.append(
            PresentationDiagnostic(
                "section_1c_domain_thumbnail_omitted",
                (
                    f"Omitted non-structure CDD thumbnail for {item_key}; "
                    f"domain text and link retained."
                ),
                "warning",
            )
        )
    return kept


def _domain_visible_heading(
    *,
    accession: str,
    name: str,
    synopsis: str | None,
) -> str:
    acc = accession.lower()
    if acc == "cd18922":
        return (
            f"{SECTION_1C_BOLD_LEAD_PHRASES[0]} found in "
            "sterol regulatory element-binding protein 2 (SREBP2) and similar proteins:"
        )
    if acc == "cl00081":
        prefix = "Conserved Protein Domain Family bHLH_SF:"
        return f"{prefix}\n\n{synopsis}" if synopsis else prefix
    if acc == "cd11304":
        prefix = "CD11304: Cadherin_repeat:"
        return f"{prefix} {synopsis}" if synopsis else prefix
    if acc == "pfam01049":
        # Internal CDD aliases (CADH_Y-type_LIR) stay in audit only.
        prefix = "Cadherin_C:"
        return f"{prefix} {synopsis}" if synopsis else prefix
    prefix = f"{name}:"
    return f"{prefix} {synopsis}" if synopsis else prefix


def _should_suppress_visible_domain(accession: str, name: str, visible_accessions: set[str]) -> bool:
    """True for broad repeat-family prose already covered by a specific block."""
    acc = accession.lower()
    lower_name = name.lower()
    if "cd11304" in visible_accessions and (acc == "smart00112" or lower_name == "ca"):
        return True
    return False


def _feature_visible_label(value: dict[str, Any]) -> str | None:
    label = str(value.get("feature_label") or value.get("feature_name") or "").strip()
    if not label:
        return None
    return label.rstrip(".:")


def build_known_structure_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
) -> SectionPresentationResult:
    """Build Section 1c polished official CDD/PDBe structure and domain blocks."""
    diagnostics: list[PresentationDiagnostic] = []
    records = _section_1c_records(list(evidence_records))
    protein_length = _protein_length_from_records(records)

    from gene_dossier.section_1c import cdd_domain_rows

    domains = cdd_domain_rows(records, protein_length=protein_length)
    selection = _candidate_selection_value(records)
    selected = _selected_candidate(selection)
    architecture = _first_section_1c_fact(records, "cdd_official_architecture_figure")
    if architecture is None:
        architecture = _first_section_1c_fact(records, "cdd_architecture_figure")
    pdbe_figure = _first_section_1c_fact(records, "pdb_official_structure_image")
    family_summaries = _fact_records(records, "cdd_family_summary")
    family_thumbnails = _fact_records(records, "cdd_family_thumbnail")
    cdd_features = _fact_records(records, "cdd_conserved_feature")
    feature_thumbnails = _fact_records(records, "cdd_feature_thumbnail")

    if not domains:
        diagnostics.append(
            PresentationDiagnostic("cdd", "no CDD conserved-domain evidence", "warning")
        )
    if not selection:
        diagnostics.append(
            PresentationDiagnostic("pdbe", "no PDBe candidate-selection evidence", "warning")
        )

    if not domains and not selection and architecture is None and pdbe_figure is None:
        return SectionPresentationResult(blocks=(), diagnostics=tuple(diagnostics))

    contributing = [
        rec
        for rec in records
        if rec.fact_type
        in {
            "conserved_domain_hit",
            "cdd_official_architecture_figure",
            "cdd_architecture_figure",
            "cdd_family_summary",
            "cdd_family_thumbnail",
            "cdd_conserved_feature",
            "cdd_feature_thumbnail",
            "pdb_candidate_selection",
            "pdb_official_structure_image",
        }
    ]
    source_ids, evidence_ids = _source_and_evidence_ids(contributing)

    intro = (
        "The NCBI Structure Group and Conserved Domain Database (CDD) are "
        "protein annotation resources utilizing multiple sequence alignment "
        "models for identification of conserved domains in protein sequences. "
        "The Protein Data Bank, accessed here through PDBe, provides "
        "experimental 3D structure data for biological macromolecules. Below "
        "is a subset of the most relevant domains and structures, with links "
        "to additional information."
    )
    if protein_length is None and domains:
        intro += (
            " Canonical protein length was unavailable in this bundle, so "
            "calculated coverage is not reported."
        )

    blocks: list[ReportContentBlock] = [
        ReportContentBlock(
            kind="narrative",
            text=intro,
            source_ids=source_ids,
            evidence_record_ids=evidence_ids,
        )
    ]

    if domains:
        domain_records = _fact_records(records, "conserved_domain_hit")
        domain_source_ids, domain_evidence_ids = _source_and_evidence_ids(domain_records)
        first_domain_url = _cdd_accession_url(domains[0].get("domain_accession"))
        blocks.append(
            ReportContentBlock(
                kind="link",
                text="Conserved domains (NCBI CDD link):",
                presentation_role="section_1c_cdd_link",
                links=[
                    {
                        "label": "Conserved domains (NCBI CDD link):",
                        "url": first_domain_url
                        or "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml",
                    }
                ],
                source_ids=domain_source_ids,
                evidence_record_ids=domain_evidence_ids,
            )
        )

    if architecture is not None:
        block = _figure_block_from_record(
            architecture,
            role="section_1c_domain_architecture_figure",
            caption="Source: NCBI Conserved Domain Database",
            diagnostics=diagnostics,
        )
        if block is not None:
            blocks.append(block)

    if domains:
        for group in _domain_groups(domains):
            rows = [row for row in group.get("rows") or [] if isinstance(row, dict)]
            group_evidence_ids = list(
                dict.fromkeys(
                    str(row.get("evidence_record_id"))
                    for row in rows
                    if row.get("evidence_record_id")
                )
            )
            accession = str(group.get("domain_accession") or "")
            visible_accessions = {
                str(item.get("domain_accession") or "").lower()
                for item in _domain_groups(domains)
                if item.get("domain_accession")
            }
            domain_label = str(group.get("domain_short_name") or accession or "CDD domain")
            if _should_suppress_visible_domain(accession, domain_label, visible_accessions):
                diagnostics.append(
                    PresentationDiagnostic(
                        "section_1c_visible_domain_suppressed",
                        f"Suppressed broad or non-Rancho-match CDD block {domain_label} ({accession}) from polished Section 1c.",
                        "info",
                    )
                )
                continue
            accession_url = _cdd_accession_url(accession)
            # Rancho body style puts the CDD link at the end of the paragraph.
            links = (
                [{"label": SECTION_1C_CDD_LINK_LABEL, "url": accession_url}]
                if accession_url
                else []
            )
            if accession.lower() == "cd18922":
                links = []
            family = _family_record_for_accession(family_summaries, accession)
            summary_text = _domain_summary_text(gene_symbol=gene_symbol, group=group)
            item_key = f"domain-{_safe_item_token(accession, fallback='unknown')}"
            if family is not None and isinstance(family.value, dict) and family.value.get("synopsis"):
                name = str(
                    family.value.get("domain_short_name")
                    or group.get("domain_short_name")
                    or group.get("domain_accession")
                    or "CDD domain"
                )
                synopsis = _strip_visible_pssm(str(family.value["synopsis"]).strip())
                summary_text = _domain_visible_heading(
                    accession=accession,
                    name=name,
                    synopsis=synopsis,
                )
                item_key = _value_item_key(family.value, prefix="domain", fallback=accession or "unknown")
            elif accession.lower() == "cd18922":
                summary_text = _domain_visible_heading(accession=accession, name=domain_label, synopsis=None)
            matching_features = [
                feature
                for feature in cdd_features
                if isinstance(feature.value, dict)
                and str(feature.value.get("domain_accession") or "").lower() == accession.lower()
            ]

            def append_feature_blocks() -> None:
                if accession.lower() == "cd18922":
                    for feature in matching_features:
                        diagnostics.append(
                            PresentationDiagnostic(
                                "section_1c_feature_moved_to_audit",
                                f"Suppressed residue-level CDD feature {feature.id} from polished bHLHzip block.",
                                "info",
                            )
                        )
                    return
                for feature in matching_features:
                    value = feature.value if isinstance(feature.value, dict) else {}
                    feature_label = _feature_visible_label(value)
                    feature_type = str(value.get("feature_type") or "").strip()
                    description = str(value.get("description") or "").strip()
                    if not feature_label or not (feature_type or description):
                        continue
                    feature_item_key = _value_item_key(
                        value,
                        prefix="feature",
                        fallback=f"{accession}-{feature_label}",
                    )
                    feature_url = accession_url or _cdd_accession_url(value.get("domain_accession"))
                    feature_link_text = f"{feature_label} (NCBI CDD Link):"
                    blocks.append(
                        ReportContentBlock(
                            kind="link",
                            text=feature_link_text,
                            presentation_role="section_1c_feature_summary",
                            presentation_item_key=feature_item_key,
                            links=(
                                [{"label": feature_link_text, "url": feature_url}]
                                if feature_url
                                else []
                            ),
                            source_ids=[feature.source_id] if feature.source_id else [],
                            evidence_record_ids=[feature.id] if feature.id else [],
                        )
                    )
                    feature_thumb = next(
                        (
                            rec
                            for rec in feature_thumbnails
                            if isinstance(rec.value, dict)
                            and str(rec.value.get("presentation_item_key") or "") == feature_item_key
                        ),
                        None,
                    )
                    if feature_thumb is not None:
                        block = _figure_block_from_record(
                            feature_thumb,
                            role="section_1c_feature_thumbnail",
                            caption="Source: NCBI Conserved Domain Database",
                            diagnostics=diagnostics,
                        )
                        if block is not None:
                            if not block.presentation_item_key:
                                block = block.model_copy(update={"presentation_item_key": feature_item_key})
                            blocks.append(block)

            family_uid = (
                _cdd_uid(family.value)
                if family is not None and isinstance(family.value, dict)
                else ""
            )
            structure_thumbs = (
                []
                if accession.lower() == "cd18922"
                else _renderable_thumbnail_records(
                    family_thumbnails,
                    item_key=item_key,
                    family_uid=family_uid,
                    diagnostics=diagnostics,
                )
            )
            # A lone thumbnail sits below its paragraph; when a family page
            # supplies several, the first leads the block and the rest follow.
            leading_thumbs, trailing_thumbs = (
                ([structure_thumbs[0]], structure_thumbs[1:])
                if len(structure_thumbs) > 1
                else ([], structure_thumbs)
            )
            pending_page_break = accession.lower() in _SECTION_1C_PAGE_LEADING_ACCESSIONS

            def append_group_block(block: ReportContentBlock) -> None:
                nonlocal pending_page_break
                if pending_page_break:
                    block = block.model_copy(
                        update={"presentation_page_break_before": True}
                    )
                    pending_page_break = False
                blocks.append(block)

            def append_thumbnail_blocks(records: Sequence[EvidenceRecord]) -> None:
                for record in records:
                    block = _figure_block_from_record(
                        record,
                        role="section_1c_domain_thumbnail",
                        caption="Source: NCBI Conserved Domain Database",
                        diagnostics=diagnostics,
                    )
                    if block is None:
                        continue
                    if not block.presentation_item_key:
                        block = block.model_copy(
                            update={"presentation_item_key": item_key}
                        )
                    append_group_block(block)

            if accession.lower() == "cd11304":
                append_feature_blocks()
            append_thumbnail_blocks(leading_thumbs)
            append_group_block(
                ReportContentBlock(
                    kind="narrative",
                    text=summary_text,
                    presentation_role="section_1c_domain_summary",
                    presentation_item_key=item_key,
                    links=links,
                    source_ids=domain_source_ids,
                    evidence_record_ids=(
                        ([family.id] if family is not None else []) or group_evidence_ids or domain_evidence_ids
                    ),
                )
            )
            append_thumbnail_blocks(trailing_thumbs)
            if accession.lower() != "cd11304":
                append_feature_blocks()
            if accession.lower() == "cd18922":
                companion = _family_record_for_accession(family_summaries, "cl00081")
                if companion is not None and isinstance(companion.value, dict):
                    companion_value = companion.value
                    companion_key = _value_item_key(
                        companion_value,
                        prefix="domain",
                        fallback="cl00081",
                    )
                    synopsis = _strip_visible_pssm(str(companion_value.get("synopsis") or "").strip())
                    blocks.append(
                        ReportContentBlock(
                            kind="narrative",
                            text=_domain_visible_heading(
                                accession="cl00081",
                                name=str(companion_value.get("domain_short_name") or "bHLH_SF"),
                                synopsis=synopsis or None,
                            ),
                            presentation_role="section_1c_domain_summary",
                            presentation_item_key=companion_key,
                            links=[
                                {
                                    "label": "Conserved Protein Domain Family bHLH_SF:",
                                    "url": _cdd_accession_url("cl00081") or "#",
                                }
                            ],
                            source_ids=[companion.source_id] if companion.source_id else [],
                            evidence_record_ids=[companion.id] if companion.id else [],
                        )
                    )
                    for companion_thumb in _renderable_thumbnail_records(
                        family_thumbnails,
                        item_key=companion_key,
                        family_uid=_cdd_uid(companion_value),
                        diagnostics=diagnostics,
                    ):
                        block = _figure_block_from_record(
                            companion_thumb,
                            role="section_1c_domain_thumbnail",
                            caption="Source: NCBI Conserved Domain Database",
                            diagnostics=diagnostics,
                        )
                        if block is not None:
                            if not block.presentation_item_key:
                                block = block.model_copy(update={"presentation_item_key": companion_key})
                            blocks.append(block)
                else:
                    diagnostics.append(
                        PresentationDiagnostic(
                            "section_1c_missing_companion_superfamily",
                            "Supported bHLHzip hit lacked parsed bHLH_SF companion-family evidence.",
                            "warning",
                        )
                    )
    else:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=f"No CDD conserved-domain evidence was available for {gene_symbol}.",
                source_ids=source_ids,
                evidence_record_ids=evidence_ids,
            )
        )

    if selected:
        candidate_rec = _first_section_1c_fact(records, "pdb_candidate_selection")
        c_source, c_evidence = _source_and_evidence_ids(
            [candidate_rec] if candidate_rec is not None else []
        )
        pdb_id = str(selected.get("pdb_id") or "")
        common = str(selected.get("species_common_name") or "human").strip().lower()
        species_label = common.capitalize() if common else "Human"
        label = f"3D structures from PDB: {species_label} {gene_symbol} protein"
        expression_host = selected.get("expression_host")
        if expression_host and common and common != "human":
            label += f" expressed in {expression_host}"
        blocks.append(
            ReportContentBlock(
                kind="link",
                text=f"{label} (PDB link)",
                presentation_role="section_1c_pdb_link",
                presentation_item_key=f"pdb-{_safe_item_token(pdb_id, fallback='unknown')}",
                links=[
                    {
                        "label": f"{label} (PDB link)",
                        "url": _pdbe_entry_url(pdb_id) or "#",
                    }
                ],
                source_ids=c_source,
                evidence_record_ids=c_evidence,
            )
        )
    elif selection:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=(
                    f"No selected PDB experimental structure was available for "
                    f"the selected human UniProt accession in this {gene_symbol} bundle."
                ),
                source_ids=source_ids,
                evidence_record_ids=evidence_ids,
            )
        )

    if pdbe_figure is not None:
        value = pdbe_figure.value if isinstance(pdbe_figure.value, dict) else {}
        block = _figure_block_from_record(
            pdbe_figure,
            role="section_1c_pdb_official_image",
            caption=str(value.get("attribution") or f"Image source: PDBe, PDB {str(value.get('pdb_id') or '').upper()}"),
            diagnostics=diagnostics,
        )
        if block is not None:
            pdb_key = f"pdb-{_safe_item_token(value.get('pdb_id'), fallback='unknown')}"
            if not block.presentation_item_key:
                block = block.model_copy(update={"presentation_item_key": pdb_key})
            blocks.append(block)
    else:
        diagnostics.append(
            PresentationDiagnostic(
                "pdbe_official_image",
                "official PDBe static image unavailable",
                "info",
            )
        )

    # The PDB heading, image, and attribution always occupy their own page.
    for index, block in enumerate(blocks):
        if block.presentation_role in _SECTION_1C_PDB_PAGE_ROLES:
            blocks[index] = block.model_copy(
                update={"presentation_page_break_before": True}
            )
            break

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


SECTION_1D_CONFIDENCE_BLURB = (
    "AlphaFold produces a per-residue model confidence score (pLDDT) between 0 and 100. "
    "Some regions below 50 pLDDT may be unstructured in isolation."
)

SECTION_1D_CONFIDENCE_LEGEND_TEXT = (
    "Model Confidence\n"
    "Very high (pLDDT > 90)\n"
    "High (90 > pLDDT > 70)\n"
    "Low (70 > pLDDT > 50)\n"
    "Very low (pLDDT < 50)\n\n"
    f"{SECTION_1D_CONFIDENCE_BLURB}"
)


def build_alphafold_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 1d polished AlphaFold blocks (human image + species links)."""
    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    slots = list(status.get("species_slots") or [])

    predictions_by_species: dict[str, EvidenceRecord] = {}
    for rec in records:
        if rec.source_name != "AlphaFold" or rec.fact_type != "alphafold_species_prediction":
            continue
        value = rec.value if isinstance(rec.value, dict) else {}
        species_key = str(value.get("species_key") or rec.species or "").lower()
        if species_key:
            predictions_by_species[species_key] = rec

    capture = next(
        (
            rec
            for rec in records
            if rec.source_name == "AlphaFold"
            and rec.fact_type == "alphafold_official_viewer_capture"
        ),
        None,
    )

    # Default slot order when status is absent: human → rat → mouse from predictions.
    if not slots:
        for species_key, label in (("human", "Human"), ("rat", "Rat"), ("mouse", "Mouse")):
            rec = predictions_by_species.get(species_key)
            if rec is None:
                continue
            value = rec.value if isinstance(rec.value, dict) else {}
            slots.append(
                {
                    "species_key": species_key,
                    "species_label": label,
                    "display_symbol": value.get("display_symbol") or gene_symbol,
                    "accession": value.get("uniprot_accession"),
                    "status": "selected",
                    "model_entity_id": value.get("model_entity_id"),
                    "entry_url": value.get("entry_url"),
                    "presentation_item_key": value.get("presentation_item_key"),
                }
            )

    blocks: list[ReportContentBlock] = []
    human_rendered = False
    for slot in slots:
        species_key = str(slot.get("species_key") or "").lower()
        label = str(slot.get("species_label") or species_key.title())
        symbol = str(slot.get("display_symbol") or gene_symbol)
        status_code = str(slot.get("status") or "")
        item_key = str(
            slot.get("presentation_item_key")
            or f"alphafold-{species_key}-unavailable"
        )

        if status_code == "selected":
            rec = predictions_by_species.get(species_key)
            value = rec.value if rec and isinstance(rec.value, dict) else {}
            url = str(slot.get("entry_url") or value.get("entry_url") or "").strip()
            if not url:
                diagnostics.append(
                    PresentationDiagnostic(
                        f"{species_key}_link",
                        "selected prediction missing entry URL",
                        "warning",
                    )
                )
                blocks.append(
                    ReportContentBlock(
                        kind="narrative",
                        text=f"{label} {symbol}: AlphaFold prediction not available",
                        presentation_role="section_1d_species_status",
                        presentation_item_key=item_key,
                    )
                )
                continue
            link_block = ReportContentBlock(
                kind="link",
                text=f"{label} {symbol}: ",
                links=[
                    {
                        "label": "AlphaFold Protein Structure Link",
                        "url": url,
                    }
                ],
                presentation_role="section_1d_species_link",
                presentation_item_key=item_key,
                source_ids=[rec.source_id] if rec and rec.source_id else [],
                evidence_record_ids=[rec.id] if rec and rec.id else [],
            )
            blocks.append(link_block)

            if species_key == "human" and not human_rendered:
                human_rendered = True
                fig = None
                if capture is not None:
                    fig = _figure_block_from_record(
                        capture,
                        role="section_1d_human_structure_capture",
                        caption="",
                        diagnostics=diagnostics,
                    )
                    if fig is not None:
                        fig = fig.model_copy(
                            update={
                                "presentation_item_key": item_key,
                            }
                        )
                        blocks.append(fig)
                        blocks.append(
                            ReportContentBlock(
                                kind="narrative",
                                text=SECTION_1D_CONFIDENCE_LEGEND_TEXT,
                                presentation_role="section_1d_confidence_legend",
                                presentation_item_key=f"{item_key}-legend",
                            )
                        )
                if fig is None:
                    diagnostics.append(
                        PresentationDiagnostic(
                            "viewer_capture",
                            "human AlphaFold viewer capture unavailable",
                            "warning",
                        )
                    )
                    blocks.append(
                        ReportContentBlock(
                            kind="narrative",
                            text=(
                                "AlphaFold structure visualization temporarily "
                                "unavailable"
                            ),
                            presentation_role="section_1d_species_status",
                            presentation_item_key=f"{item_key}-visualization-unavailable",
                        )
                    )
            continue

        # Unavailable / failed — visible status line, no evidence refs.
        text = slot.get("message") or f"{label} {symbol}: AlphaFold prediction not available"
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=str(text),
                presentation_role="section_1d_species_status",
                presentation_item_key=item_key,
            )
        )
        diagnostics.append(
            PresentationDiagnostic(
                f"{species_key}_status",
                str(status_code or "unavailable"),
                "info",
            )
        )

    if not blocks:
        diagnostics.append(
            PresentationDiagnostic(
                "section_1d",
                "no AlphaFold species slots or predictions available",
                "warning",
            )
        )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


def build_homologues_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 1e polished NCBI Orthologs / OrthoDB blocks."""
    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    summary_rec = next(
        (
            rec
            for rec in records
            if rec.fact_type == "ortholog_collection_summary"
            and (
                rec.source_name == "NCBI Datasets"
                or rec.subsection == "Homologues in model animals"
            )
        ),
        None,
    )
    capture_rec = next(
        (rec for rec in records if rec.fact_type == "ortholog_table_capture"),
        None,
    )
    orthodb_rec = next(
        (
            rec
            for rec in records
            if rec.source_name == "OrthoDB" and rec.fact_type == "orthodb_gene_search"
        ),
        None,
    )

    summary_value: dict[str, Any] = {}
    if summary_rec is not None and isinstance(summary_rec.value, dict):
        summary_value = dict(summary_rec.value)
    elif isinstance(status.get("summary"), dict):
        summary_value = dict(status["summary"])

    if not summary_value:
        diagnostics.append(
            PresentationDiagnostic(
                "section_1e",
                "ortholog collection summary unavailable",
                "warning",
            )
        )
        return SectionPresentationResult(blocks=(), diagnostics=tuple(diagnostics))

    item_key = str(
        summary_value.get("presentation_item_key")
        or f"orthologs-{(gene_symbol or '').lower()}"
    )
    narrative = str(
        summary_value.get("narrative")
        or (summary_rec.display_text if summary_rec else "")
        or ""
    ).strip()
    table_status = str(
        summary_value.get("table_status")
        or (status.get("rendering_status") or {}).get("table_status")
        or "unavailable"
    )
    ncbi_url = str(summary_value.get("ncbi_url") or "").strip()
    orthodb_url = str(summary_value.get("orthodb_url") or "").strip()

    blocks: list[ReportContentBlock] = []
    if narrative:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=narrative,
                presentation_role="section_1e_narrative",
                presentation_item_key=item_key,
                source_ids=[summary_rec.source_id]
                if summary_rec and summary_rec.source_id
                else [],
                evidence_record_ids=[summary_rec.id]
                if summary_rec and summary_rec.id
                else [],
            )
        )

    capture_fig = None
    if table_status == "official_capture" and capture_rec is not None:
        capture_fig = _figure_block_from_record(
            capture_rec,
            role="section_1e_ortholog_capture",
            caption="",
            diagnostics=diagnostics,
        )
        if capture_fig is not None:
            capture_fig = capture_fig.model_copy(
                update={"presentation_item_key": item_key}
            )
            blocks.append(capture_fig)
        else:
            diagnostics.append(
                PresentationDiagnostic(
                    "viewer_capture",
                    "ortholog table capture path unresolved",
                    "warning",
                )
            )

    if capture_fig is None and table_status in {
        "complete_api_fallback",
        "partial_api_fallback",
        "official_capture",
    }:
        rows_raw = summary_value.get("fallback_rows") or []
        headers = ["Species", "Gene", "Description", "Gene ID"]
        rows: list[list[str]] = []
        for row in rows_raw:
            if not isinstance(row, dict):
                continue
            rows.append(
                [
                    str(row.get("species") or ""),
                    str(row.get("gene") or ""),
                    str(row.get("description") or ""),
                    str(row.get("gene_id") or ""),
                ]
            )
        if rows:
            blocks.append(
                ReportContentBlock(
                    kind="table",
                    table_headers=headers,
                    table_rows=rows,
                    presentation_role="section_1e_fallback_table",
                    presentation_item_key=f"{item_key}-fallback",
                    source_ids=[summary_rec.source_id]
                    if summary_rec and summary_rec.source_id
                    else [],
                    evidence_record_ids=[summary_rec.id]
                    if summary_rec and summary_rec.id
                    else [],
                )
            )
            if table_status == "official_capture":
                diagnostics.append(
                    PresentationDiagnostic(
                        "viewer_capture",
                        "falling back to structured ortholog table",
                        "warning",
                    )
                )
        elif table_status != "unavailable":
            diagnostics.append(
                PresentationDiagnostic(
                    "fallback_table",
                    "fallback table requested but no scoped rows available",
                    "warning",
                )
            )

    attr_links: list[dict[str, str]] = []
    if ncbi_url:
        attr_links.append({"label": "NCBI Orthologs", "url": ncbi_url})
    if orthodb_url:
        attr_links.append({"label": "OrthoDB", "url": orthodb_url})
    if attr_links:
        source_ids = []
        evidence_ids = []
        if summary_rec and summary_rec.source_id:
            source_ids.append(summary_rec.source_id)
        if summary_rec and summary_rec.id:
            evidence_ids.append(summary_rec.id)
        if orthodb_rec and orthodb_rec.source_id:
            source_ids.append(orthodb_rec.source_id)
        if orthodb_rec and orthodb_rec.id:
            evidence_ids.append(orthodb_rec.id)
        blocks.append(
            ReportContentBlock(
                kind="link",
                text="Sources: ",
                links=attr_links,
                presentation_role="section_1e_attribution",
                presentation_item_key=f"{item_key}-attribution",
                source_ids=list(dict.fromkeys(source_ids)),
                evidence_record_ids=list(dict.fromkeys(evidence_ids)),
            )
        )

    if not blocks:
        diagnostics.append(
            PresentationDiagnostic(
                "section_1e",
                "no homologues presentation blocks available",
                "warning",
            )
        )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


def build_tissue_specific_information_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 2a GTEx + HBT presentation blocks in golden order."""
    from gene_dossier.section_2a import (
        gtex_gene_url,
        gtex_intro_text,
        hbt_intro_text,
        hbt_link_text,
        hbt_pdf_url,
    )

    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    summary = dict(status.get("summary") or {})
    rendering = dict(status.get("rendering_status") or {})
    gene = (gene_symbol or "").strip()
    item_key = str(summary.get("presentation_item_key") or f"tissue-{(gene or '').lower()}")

    collection = next(
        (r for r in records if r.fact_type == "gtex_expression_collection_summary"),
        None,
    )
    all_fig_rec = next(
        (r for r in records if r.fact_type == "gtex_all_tissues_figure"),
        None,
    )
    brain_fig_rec = next(
        (r for r in records if r.fact_type == "gtex_brain_tissues_figure"),
        None,
    )
    hbt_fig_rec = next(
        (r for r in records if r.fact_type == "hbt_whole_brain_figure"),
        None,
    )
    hbt_pdf_rec = next(
        (r for r in records if r.fact_type == "hbt_whole_brain_pdf"),
        None,
    )
    gtex_ref = next(
        (r for r in records if r.fact_type == "gtex_gene_reference"),
        None,
    )

    gtex_url = str(summary.get("gtex_portal_url") or gtex_gene_url(gene))
    hbt_url = str(summary.get("hbt_pdf_url") or hbt_pdf_url(gene))
    blocks: list[ReportContentBlock] = []

    # 1. GTEx intro
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=gtex_intro_text(gene),
            presentation_role="section_2a_gtex_intro",
            presentation_item_key=item_key,
            source_ids=[gtex_ref.source_id]
            if gtex_ref and gtex_ref.source_id
            else ([collection.source_id] if collection and collection.source_id else []),
            evidence_record_ids=[gtex_ref.id]
            if gtex_ref and gtex_ref.id
            else ([collection.id] if collection and collection.id else []),
        )
    )

    # 2. All-tissues link
    blocks.append(
        ReportContentBlock(
            kind="link",
            text="",
            links=[{"label": f"GTEx: {gene} expression in all tissues", "url": gtex_url}],
            presentation_role="section_2a_gtex_all_tissues_link",
            presentation_item_key=item_key,
            source_ids=[collection.source_id]
            if collection and collection.source_id
            else [],
            evidence_record_ids=[collection.id] if collection and collection.id else [],
        )
    )

    # 3. All-tissues figure or status
    all_fig = None
    if all_fig_rec is not None:
        all_fig = _figure_block_from_record(
            all_fig_rec,
            role="section_2a_gtex_all_tissues_figure",
            caption="",
            diagnostics=diagnostics,
        )
        if all_fig is not None:
            all_fig = all_fig.model_copy(update={"presentation_item_key": item_key})
            blocks.append(all_fig)
    if all_fig is None:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="GTEx all-tissue expression plot temporarily unavailable.",
                presentation_role="section_2a_source_status",
                presentation_item_key=f"{item_key}-gtex-all-status",
            )
        )
        diagnostics.append(
            PresentationDiagnostic(
                "section_2a",
                "gtex all-tissue figure unavailable",
                "warning",
            )
        )

    # 4. Brain link (new page)
    blocks.append(
        ReportContentBlock(
            kind="link",
            text="",
            links=[{"label": f"GTEx: {gene} expression in brain", "url": gtex_url}],
            presentation_role="section_2a_gtex_brain_link",
            presentation_item_key=item_key,
            presentation_page_break_before=True,
            source_ids=[collection.source_id]
            if collection and collection.source_id
            else [],
            evidence_record_ids=[collection.id] if collection and collection.id else [],
        )
    )

    # 5. Brain figure or status
    brain_fig = None
    if brain_fig_rec is not None:
        brain_fig = _figure_block_from_record(
            brain_fig_rec,
            role="section_2a_gtex_brain_figure",
            caption="",
            diagnostics=diagnostics,
        )
        if brain_fig is not None:
            brain_fig = brain_fig.model_copy(update={"presentation_item_key": item_key})
            blocks.append(brain_fig)
    if brain_fig is None:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="GTEx brain-tissue expression plot temporarily unavailable.",
                presentation_role="section_2a_source_status",
                presentation_item_key=f"{item_key}-gtex-brain-status",
            )
        )

    # 6. HBT intro
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=hbt_intro_text(),
            presentation_role="section_2a_hbt_intro",
            presentation_item_key=item_key,
            source_ids=[hbt_pdf_rec.source_id]
            if hbt_pdf_rec and hbt_pdf_rec.source_id
            else [],
            evidence_record_ids=[hbt_pdf_rec.id]
            if hbt_pdf_rec and hbt_pdf_rec.id
            else [],
            links=[{"label": "HBT", "url": "https://hbatlas.org/"}],
        )
    )

    # 7. HBT gene-expression link
    blocks.append(
        ReportContentBlock(
            kind="link",
            text=hbt_link_text(gene),
            links=[{"label": "Link", "url": hbt_url}],
            presentation_role="section_2a_hbt_link",
            presentation_item_key=item_key,
            source_ids=[hbt_pdf_rec.source_id]
            if hbt_pdf_rec and hbt_pdf_rec.source_id
            else [],
            evidence_record_ids=[hbt_pdf_rec.id]
            if hbt_pdf_rec and hbt_pdf_rec.id
            else [],
        )
    )

    # 8. HBT figure (new page) or status
    hbt_fig = None
    if hbt_fig_rec is not None:
        hbt_fig = _figure_block_from_record(
            hbt_fig_rec,
            role="section_2a_hbt_figure",
            caption="",
            diagnostics=diagnostics,
        )
        if hbt_fig is not None:
            hbt_fig = hbt_fig.model_copy(
                update={
                    "presentation_item_key": item_key,
                    "presentation_page_break_before": True,
                }
            )
            blocks.append(hbt_fig)
    if hbt_fig is None:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="HBT developmental-expression figure temporarily unavailable.",
                presentation_role="section_2a_source_status",
                presentation_item_key=f"{item_key}-hbt-status",
                presentation_page_break_before=True,
            )
        )
        diagnostics.append(
            PresentationDiagnostic(
                "section_2a",
                "hbt figure unavailable",
                "warning",
            )
        )

    if rendering.get("overall") in {None, "empty"} and not any(
        r.fact_type.startswith("gtex_") or r.fact_type.startswith("hbt_")
        for r in records
    ):
        diagnostics.append(
            PresentationDiagnostic(
                "section_2a",
                "no tissue-specific evidence available",
                "warning",
            )
        )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


def build_barres_brain_expression_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 2b Allen/Barres presentation blocks in golden order."""
    from gene_dossier.section_2b import (
        BRAINRNASEQ_HOME_URL,
        BRAINRNASEQ_SOURCE_LINK_LABEL,
        CATEGORY_NOTE,
        NOT_DETERMINED,
        format_numeric_cell,
        section_2b_celltype_intro_text,
        section_2b_intro_text,
    )

    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    summary = dict(status.get("summary") or {})
    rendering = dict(status.get("rendering_status") or {})
    gene = (gene_symbol or "").strip()
    item_key = str(summary.get("presentation_item_key") or f"barres-{(gene or '').lower()}")

    summary_rec = next(
        (r for r in records if r.fact_type == "section_2b_summary_table"),
        None,
    )
    category_rec = next(
        (r for r in records if r.fact_type == "section_2b_category_status"),
        None,
    )
    figure_rec = next(
        (r for r in records if r.fact_type == "brainrnaseq_celltype_figure"),
        None,
    )
    selection_rec = next(
        (r for r in records if r.fact_type == "allen_human_microarray_selection_summary"),
        None,
    )
    pooled_rec = next(
        (r for r in records if r.fact_type == "allen_human_rnaseq_pooled_summary"),
        None,
    )
    mouse_tpm_rec = next(
        (r for r in records if r.fact_type == "brainrnaseq_mouse_average_tpm"),
        None,
    )

    summary_val = (
        summary_rec.value
        if summary_rec is not None and isinstance(summary_rec.value, dict)
        else {}
    )
    agilent = summary.get("average_human_brain_agilent_expression")
    if agilent is None:
        agilent = summary_val.get("average_human_brain_agilent_expression")
    human_tpm = summary.get("average_human_brain_rnaseq_expression_tpm")
    if human_tpm is None:
        human_tpm = summary_val.get("average_human_brain_rnaseq_expression_tpm")
        if human_tpm is None and pooled_rec is not None and isinstance(pooled_rec.value, dict):
            human_tpm = pooled_rec.value.get("mean_tpm")
    mouse_tpm = summary.get("average_mouse_brain_rnaseq_expression_tpm")
    if mouse_tpm is None:
        mouse_tpm = summary_val.get("average_mouse_brain_rnaseq_expression_tpm")
        if (
            mouse_tpm is None
            and mouse_tpm_rec is not None
            and isinstance(mouse_tpm_rec.value, dict)
        ):
            mouse_tpm = mouse_tpm_rec.value.get("mean_tpm")

    brs_url = str(summary.get("brainrnaseq_url") or BRAINRNASEQ_HOME_URL)
    blocks: list[ReportContentBlock] = []

    # 1. Intro
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=section_2b_intro_text(gene),
            presentation_role="section_2b_intro",
            presentation_item_key=item_key,
            source_ids=[summary_rec.source_id]
            if summary_rec and summary_rec.source_id
            else ([selection_rec.source_id] if selection_rec and selection_rec.source_id else []),
            evidence_record_ids=[summary_rec.id]
            if summary_rec and summary_rec.id
            else ([selection_rec.id] if selection_rec and selection_rec.id else []),
        )
    )

    # 2. Summary table (5 columns, Rancho green header via CSS)
    headers = [
        "Average Human Brain Agilent Expression",
        "Average Human Brain RNA-Seq Expression (TPM)",
        "Average Mouse Brain RNA-Seq Expression (TPM)",
        "Human Brain Expression Category",
        "Mouse Brain Expression Category",
    ]
    row = [
        format_numeric_cell(agilent if isinstance(agilent, (int, float)) else None),
        format_numeric_cell(human_tpm if isinstance(human_tpm, (int, float)) else None),
        format_numeric_cell(mouse_tpm if isinstance(mouse_tpm, (int, float)) else None),
        NOT_DETERMINED,
        NOT_DETERMINED,
    ]
    blocks.append(
        ReportContentBlock(
            kind="table",
            text="",
            table_headers=headers,
            table_rows=[row],
            presentation_role="section_2b_summary_table",
            presentation_item_key=item_key,
            source_ids=[summary_rec.source_id]
            if summary_rec and summary_rec.source_id
            else [],
            evidence_record_ids=[summary_rec.id] if summary_rec and summary_rec.id else [],
        )
    )

    # 3. Category note (only while unresolved)
    category_policy = "threshold_policy_unresolved"
    if summary.get("category_policy"):
        category_policy = str(summary.get("category_policy"))
    elif category_rec is not None and isinstance(category_rec.value, dict):
        category_policy = str(
            category_rec.value.get("policy") or "threshold_policy_unresolved"
        )
    if category_policy == "threshold_policy_unresolved":
        note = CATEGORY_NOTE
        if category_rec is not None and isinstance(category_rec.value, dict):
            note = str(category_rec.value.get("note") or CATEGORY_NOTE)
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=note,
                presentation_role="section_2b_category_status",
                presentation_item_key=item_key,
                source_ids=[category_rec.source_id]
                if category_rec and category_rec.source_id
                else [],
                evidence_record_ids=[category_rec.id]
                if category_rec and category_rec.id
                else [],
            )
        )

    # 4. Cell-type intro (new page)
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=section_2b_celltype_intro_text(gene),
            presentation_role="section_2b_celltype_intro",
            presentation_item_key=item_key,
            presentation_page_break_before=True,
            links=[{"label": "brainrnaseq.org", "url": brs_url}],
        )
    )

    # 5. Source link
    blocks.append(
        ReportContentBlock(
            kind="link",
            text=BRAINRNASEQ_SOURCE_LINK_LABEL,
            links=[{"label": BRAINRNASEQ_SOURCE_LINK_LABEL, "url": brs_url}],
            presentation_role="section_2b_source_link",
            presentation_item_key=item_key,
            source_ids=[figure_rec.source_id]
            if figure_rec and figure_rec.source_id
            else [],
            evidence_record_ids=[figure_rec.id] if figure_rec and figure_rec.id else [],
        )
    )

    # 6. Combined cell-type figure or status
    cell_fig = None
    if figure_rec is not None:
        cell_fig = _figure_block_from_record(
            figure_rec,
            role="section_2b_celltype_figure",
            caption="",
            diagnostics=diagnostics,
        )
        if cell_fig is not None:
            cell_fig = cell_fig.model_copy(update={"presentation_item_key": item_key})
            blocks.append(cell_fig)
    if cell_fig is None:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="Barres Lab cell-type expression figure temporarily unavailable.",
                presentation_role="section_2b_source_status",
                presentation_item_key=f"{item_key}-figure-status",
            )
        )
        diagnostics.append(
            PresentationDiagnostic(
                "section_2b",
                "cell-type figure unavailable",
                "warning",
            )
        )

    if rendering.get("overall") in {None, "empty"} and not any(
        r.fact_type.startswith(
            (
                "allen_",
                "brainrnaseq_",
                "section_2b_",
            )
        )
        for r in records
    ):
        diagnostics.append(
            PresentationDiagnostic(
                "section_2b",
                "no Barres/Allen brain-expression evidence available",
                "warning",
            )
        )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


@dataclass(frozen=True)
class _Section2cFigureSlot:
    """One optional Section 2c figure slot and its evidence record."""

    role: str
    item_key: str
    record: EvidenceRecord | None
    fallback_note: str


def build_section_2c_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 2c snRNA-Seq cell-type presentation blocks in golden order."""
    from gene_dossier.tools.allen_celltypes import (
        DATASET_HUMAN_M1,
        DATASET_MOUSE_CTX_HPF,
        VISUALIZATION_HEATMAP,
        VISUALIZATION_SCATTER,
    )

    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    summary = dict(status.get("summary") or {})
    rendering = dict(status.get("rendering_status") or {})
    item_key = str(
        summary.get("presentation_item_key")
        or f"celltype-{(gene_symbol or '').strip().lower()}"
    )

    human = dict(summary.get("human") or {})
    mouse = dict(summary.get("mouse") or {})
    dropviz = dict(summary.get("dropviz") or {})

    def _record(fact_type: str) -> EvidenceRecord | None:
        return next((r for r in records if r.fact_type == fact_type), None)

    def _ids(rec: EvidenceRecord | None) -> tuple[list[str], list[str]]:
        if rec is None:
            return [], []
        return (
            [rec.source_id] if rec.source_id else [],
            [rec.id] if rec.id else [],
        )

    def _explorer_figure_record(
        dataset: str,
        visualization: str,
    ) -> EvidenceRecord | None:
        for rec in records:
            if rec.fact_type != "allen_celltype_explorer_figure":
                continue
            value = rec.value if isinstance(rec.value, dict) else {}
            if (
                str(value.get("dataset") or "") == dataset
                and str(value.get("visualization") or "") == visualization
            ):
                return rec
        return None

    def _value_cell(entry: dict[str, Any], value_key: str) -> str:
        display = entry.get("value_display")
        if display is not None and str(display) != "":
            return str(display)
        value = entry.get(value_key)
        return "" if value is None else str(value)

    blocks: list[ReportContentBlock] = []

    def _append_link(url: Any, link_label: Any, *, item_suffix: str, source_ids, evidence_ids) -> None:
        if not url:
            return
        label = str(link_label or url)
        blocks.append(
            ReportContentBlock(
                kind="link",
                text=label,
                links=[{"label": label, "url": str(url)}],
                presentation_role="section_2c_source_link",
                presentation_item_key=f"{item_key}-{item_suffix}",
                source_ids=source_ids,
                evidence_record_ids=evidence_ids,
            )
        )

    def _append_narrative(
        text: str,
        *,
        role: str,
        source_ids,
        evidence_ids,
        page_break_before: bool = False,
        item_suffix: str | None = None,
    ) -> None:
        if not str(text or "").strip():
            return
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=str(text),
                presentation_role=role,  # type: ignore[arg-type]
                presentation_item_key=(
                    f"{item_key}-{item_suffix}" if item_suffix else item_key
                ),
                presentation_page_break_before=page_break_before,
                source_ids=source_ids,
                evidence_record_ids=evidence_ids,
            )
        )

    def _append_status_notes(branch: dict[str, Any], *, item_suffix: str) -> None:
        notes = [
            str(note)
            for note in (branch.get("figure_status_notes") or [])
            if str(note).strip()
        ]
        status_texts: list[str] = []
        status_note = branch.get("status_note")
        if status_note and str(status_note).strip():
            status_texts.append(str(status_note))
        status_texts.extend(notes)
        for text_value in status_texts:
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text=text_value,
                    presentation_role="section_2c_source_status",
                    presentation_item_key=f"{item_key}-{item_suffix}-status",
                )
            )

    def _append_figure_slots(
        figures: Sequence[_Section2cFigureSlot],
        *,
        page_break_first: bool = False,
    ) -> tuple[bool, list[_Section2cFigureSlot]]:
        missing: list[_Section2cFigureSlot] = []
        rendered = False
        for slot in figures:
            figure_block = None
            if slot.record is not None:
                figure_block = _figure_block_from_record(
                    slot.record,
                    role=slot.role,
                    caption="",
                    diagnostics=diagnostics,
                )
            if figure_block is None:
                missing.append(slot)
                continue
            update: dict[str, Any] = {"presentation_item_key": slot.item_key}
            if page_break_first and not rendered:
                update["presentation_page_break_before"] = True
            blocks.append(figure_block.model_copy(update=update))
            rendered = True
        return rendered, missing

    def _append_allen_branch(
        branch: dict[str, Any],
        *,
        record: EvidenceRecord | None,
        narrative_role: str,
        scatter_narrative_role: str,
        heatmap_narrative_role: str,
        table_role: str,
        item_suffix: str,
        figures: Sequence[_Section2cFigureSlot],
        page_break_before: bool,
    ) -> None:
        source_ids, evidence_ids = _ids(record)
        figures_complete = bool(branch.get("figures_complete"))
        _append_narrative(
            str(branch.get("narrative") or ""),
            role=narrative_role,
            source_ids=source_ids,
            evidence_ids=evidence_ids,
            page_break_before=page_break_before,
        )
        # Primary: dataset-specific source page. Secondary: gene explorer deep link.
        _append_link(
            branch.get("source_page_url") or branch.get("database_url"),
            branch.get("source_link_label") or branch.get("database_link_label"),
            item_suffix=f"{item_suffix}-source",
            source_ids=source_ids,
            evidence_ids=evidence_ids,
        )
        _append_link(
            branch.get("explorer_url"),
            branch.get("explorer_link_label"),
            item_suffix=f"{item_suffix}-explorer",
            source_ids=source_ids,
            evidence_ids=evidence_ids,
        )

        if figures_complete:
            scatter_slots = [slot for slot in figures if "scatter" in slot.role]
            heatmap_slots = [slot for slot in figures if "heatmap" in slot.role]
            _append_narrative(
                str(branch.get("scatter_interpretation") or ""),
                role=scatter_narrative_role,
                source_ids=source_ids,
                evidence_ids=evidence_ids,
                item_suffix=f"{item_suffix}-scatter-note",
            )
            _append_figure_slots(scatter_slots)
            _append_narrative(
                str(branch.get("heatmap_interpretation") or ""),
                role=heatmap_narrative_role,
                source_ids=source_ids,
                evidence_ids=evidence_ids,
                item_suffix=f"{item_suffix}-heatmap-note",
            )
            _append_figure_slots(heatmap_slots)
            _append_status_notes(branch, item_suffix=item_suffix)
            return

        # Figure failure path: structured narrative already above; top-five table;
        # then unavailable notes for missing figures.
        entries = [e for e in (branch.get("top") or []) if isinstance(e, dict)]
        if str(branch.get("analysis_status") or "") == "success" and entries:
            blocks.append(
                ReportContentBlock(
                    kind="table",
                    text="",
                    table_headers=[
                        str(branch.get("label_column_header") or ""),
                        str(branch.get("value_label") or ""),
                    ],
                    table_rows=[
                        [str(e.get("label") or ""), _value_cell(e, "value")]
                        for e in entries
                    ],
                    presentation_role=table_role,  # type: ignore[arg-type]
                    presentation_item_key=item_key,
                    source_ids=source_ids,
                    evidence_record_ids=evidence_ids,
                )
            )
        _, missing = _append_figure_slots(figures)
        _append_status_notes(branch, item_suffix=item_suffix)
        notes = [
            str(note)
            for note in (branch.get("figure_status_notes") or [])
            if str(note).strip()
        ]
        for index, slot in enumerate(missing):
            reason = notes[index] if index < len(notes) else slot.fallback_note
            diagnostics.append(PresentationDiagnostic("figure_note", reason, "info"))

    def _append_dropviz_branch(
        branch: dict[str, Any],
        *,
        record: EvidenceRecord | None,
        figures: Sequence[_Section2cFigureSlot],
        page_break_before: bool,
        figure_page_break_before: bool,
    ) -> bool:
        source_ids, evidence_ids = _ids(record)
        _append_narrative(
            str(branch.get("narrative") or ""),
            role="section_2c_dropviz_narrative",
            source_ids=source_ids,
            evidence_ids=evidence_ids,
            page_break_before=page_break_before,
        )
        population_note = str(branch.get("population_identifier_note") or "").strip()
        if population_note:
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text=population_note,
                    presentation_role="section_2c_source_status",
                    presentation_item_key=f"{item_key}-dropviz-population-note",
                    source_ids=source_ids,
                    evidence_record_ids=evidence_ids,
                )
            )
        entries = [
            e for e in (branch.get("top_populations") or []) if isinstance(e, dict)
        ]
        # Compact top-ten optional when ranking succeeded; keep for evidence fidelity.
        if str(branch.get("analysis_status") or "") == "success" and entries:
            blocks.append(
                ReportContentBlock(
                    kind="table",
                    text="",
                    table_headers=[
                        str(branch.get("label_column_header") or ""),
                        str(branch.get("value_label") or ""),
                    ],
                    table_rows=[
                        [
                            str(e.get("population_label") or ""),
                            _value_cell(e, "ranking_value"),
                        ]
                        for e in entries[:10]
                    ],
                    presentation_role="section_2c_dropviz_table",
                    presentation_item_key=item_key,
                    source_ids=source_ids,
                    evidence_record_ids=evidence_ids,
                )
            )
        rendered, missing = _append_figure_slots(
            figures, page_break_first=figure_page_break_before
        )
        attribution = str(branch.get("geo_attribution") or "").strip()
        if attribution:
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text=attribution,
                    presentation_role="section_2c_geo_attribution",
                    presentation_item_key=f"{item_key}-dropviz-geo",
                    source_ids=source_ids,
                    evidence_record_ids=evidence_ids,
                )
            )
        # Do not render DropViz homepage or saved-state links. GEO series is the
        # only polished quantitative source attribution (via geo_attribution).
        _append_status_notes(branch, item_suffix="dropviz")
        notes = [
            str(note)
            for note in (branch.get("figure_status_notes") or [])
            if str(note).strip()
        ]
        for index, slot in enumerate(missing):
            reason = notes[index] if index < len(notes) else slot.fallback_note
            diagnostics.append(PresentationDiagnostic("figure_note", reason, "info"))
        tsne_note = str(branch.get("regional_tsne_limitation_note") or "").strip()
        if tsne_note:
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text=tsne_note,
                    presentation_role="section_2c_source_status",
                    presentation_item_key=f"{item_key}-dropviz-tsne-limitation",
                    source_ids=source_ids,
                    evidence_record_ids=evidence_ids,
                )
            )
        return rendered

    def _figure_note(branch: dict[str, Any], visualization: str) -> str:
        label = str(branch.get("dataset_label") or branch.get("dataset") or "").strip()
        prefix = f"{label} " if label else ""
        figure_status = str(
            (branch.get("figure_statuses") or {}).get(visualization) or ""
        )
        if figure_status == "not_attempted_optional":
            return f"{prefix}{visualization} figure not attempted"
        return f"{prefix}{visualization} figure unavailable"

    summary_rec = _record("section_2c_summary")
    summary_sources, summary_evidence = _ids(summary_rec)
    human_rec = _record("allen_human_celltype_summary")
    mouse_rec = _record("allen_mouse_celltype_summary")
    dropviz_rec = _record("dropviz_top_populations")

    # 1. Section intro
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=str(summary.get("intro_text") or ""),
            presentation_role="section_2c_intro",
            presentation_item_key=item_key,
            source_ids=summary_sources,
            evidence_record_ids=summary_evidence,
        )
    )

    _append_allen_branch(
        human,
        record=human_rec,
        narrative_role="section_2c_human_narrative",
        scatter_narrative_role="section_2c_human_scatter_narrative",
        heatmap_narrative_role="section_2c_human_heatmap_narrative",
        table_role="section_2c_human_table",
        item_suffix="human",
        figures=(
            _Section2cFigureSlot(
                role="section_2c_human_scatter_figure",
                item_key=f"{item_key}-human-scatter",
                record=_explorer_figure_record(DATASET_HUMAN_M1, VISUALIZATION_SCATTER),
                fallback_note=_figure_note(human, VISUALIZATION_SCATTER),
            ),
            _Section2cFigureSlot(
                role="section_2c_human_heatmap_figure",
                item_key=f"{item_key}-human-heatmap",
                record=_explorer_figure_record(DATASET_HUMAN_M1, VISUALIZATION_HEATMAP),
                fallback_note=_figure_note(human, VISUALIZATION_HEATMAP),
            ),
        ),
        page_break_before=False,
    )

    _append_allen_branch(
        mouse,
        record=mouse_rec,
        narrative_role="section_2c_mouse_narrative",
        scatter_narrative_role="section_2c_mouse_scatter_narrative",
        heatmap_narrative_role="section_2c_mouse_heatmap_narrative",
        table_role="section_2c_mouse_table",
        item_suffix="mouse",
        figures=(
            _Section2cFigureSlot(
                role="section_2c_mouse_scatter_figure",
                item_key=f"{item_key}-mouse-scatter",
                record=_explorer_figure_record(
                    DATASET_MOUSE_CTX_HPF, VISUALIZATION_SCATTER
                ),
                fallback_note=_figure_note(mouse, VISUALIZATION_SCATTER),
            ),
            _Section2cFigureSlot(
                role="section_2c_mouse_heatmap_figure",
                item_key=f"{item_key}-mouse-heatmap",
                record=_explorer_figure_record(
                    DATASET_MOUSE_CTX_HPF, VISUALIZATION_HEATMAP
                ),
                fallback_note=_figure_note(mouse, VISUALIZATION_HEATMAP),
            ),
        ),
        page_break_before=True,
    )

    dropviz_figure_rendered = _append_dropviz_branch(
        dropviz,
        record=dropviz_rec,
        figures=(
            _Section2cFigureSlot(
                role="section_2c_dropviz_rank_figure",
                item_key=f"{item_key}-dropviz-rank",
                record=_record("dropviz_top_populations_figure"),
                fallback_note="DropViz population ranking figure unavailable",
            ),
        ),
        page_break_before=True,
        figure_page_break_before=True,
    )

    # Therapeutic implications; shares the ranking figure's page when there
    # is one, otherwise starts its own.
    therapeutic = str(summary.get("therapeutic_narrative") or "").strip()
    if therapeutic:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=therapeutic,
                presentation_role="section_2c_therapeutic_narrative",
                presentation_item_key=item_key,
                presentation_page_break_before=not dropviz_figure_rendered,
                source_ids=summary_sources,
                evidence_record_ids=summary_evidence,
            )
        )

    for branch, reason in (
        (human, "Allen Human M1 structured analysis unavailable"),
        (mouse, "Allen mouse cortex and hippocampus structured analysis unavailable"),
        (dropviz, "DropViz GSE116470 population ranking unavailable"),
    ):
        if str(branch.get("analysis_status") or "") != "success":
            diagnostics.append(PresentationDiagnostic("section_2c", reason, "warning"))

    for issue in summary.get("unresolved_issues") or []:
        if not isinstance(issue, dict):
            continue
        diagnostics.append(
            PresentationDiagnostic(
                "section_2c",
                f"{issue.get('field', '')}: {issue.get('reason', '')}",
                "warning",
            )
        )

    if rendering.get("overall") in {None, "empty"} and not any(
        r.fact_type.startswith(
            (
                "allen_celltype_",
                "allen_human_celltype",
                "allen_mouse_celltype",
                "dropviz_",
                "section_2c_",
            )
        )
        for r in records
    ):
        diagnostics.append(
            PresentationDiagnostic(
                "section_2c",
                "no single-cell / single-nucleus cell-type evidence available",
                "warning",
            )
        )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


def build_section_3a_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 3a GEO Profiles presentation blocks."""
    from gene_dossier.section_3a import (
        COMPARABILITY_NOTE,
        SCREENING_CAVEAT,
        SELECTION_POLICY,
        STATUS_NOT_ATTEMPTED,
        STATUS_SOURCE_UNAVAILABLE,
        STATUS_SUCCESS,
        build_intro_text,
    )

    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    summary = dict(status.get("summary") or {})
    rendering = dict(status.get("rendering_status") or {})
    item_key = str(
        summary.get("presentation_item_key")
        or f"geo-profiles-{(gene_symbol or '').strip().lower()}"
    )

    summary_rec = next((r for r in records if r.fact_type == "section_3a_summary"), None)
    profile_recs = [r for r in records if r.fact_type == "section_3a_profile"]
    profile_by_uid = {
        str((r.value or {}).get("profile_uid") or ""): r for r in profile_recs
    }

    def _ids(rec: EvidenceRecord | None) -> tuple[list[str], list[str]]:
        if rec is None:
            return [], []
        return (
            [rec.source_id] if rec.source_id else [],
            [rec.id] if rec.id else [],
        )

    blocks: list[ReportContentBlock] = []
    scientific = str(rendering.get("scientific_status") or summary.get("scientific_status") or "")
    visual = str(rendering.get("visual_status") or summary.get("visual_status") or "")

    intro = str(
        summary.get("intro_text")
        or build_intro_text(
            gene_symbol,
            exact_count=summary.get("exact_profile_count"),
            neural_count=summary.get("neural_profile_count"),
            subset_count=summary.get("subset_effect_profile_count"),
        )
    )
    src_ids, ev_ids = _ids(summary_rec)
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=intro,
            presentation_role="section_3a_intro",
            presentation_item_key=f"{item_key}-intro",
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
        )
    )

    if scientific == STATUS_SOURCE_UNAVAILABLE:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="GEO Profiles discovery or metadata enrichment was unavailable for this run.",
                presentation_role="section_3a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
        diagnostics.append(
            PresentationDiagnostic("section_3a", "source unavailable", "warning")
        )
        return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))

    selected = list(summary.get("selected_profiles") or [])
    if not selected and scientific != STATUS_SUCCESS:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="No potentially relevant GEO Profiles met the neural/perturbation screening criteria.",
                presentation_role="section_3a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
    elif selected:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="The following GEO Profiles were selected as potentially relevant:",
                presentation_role="section_3a_intro",
                presentation_item_key=f"{item_key}-selected-lead",
                source_ids=src_ids,
                evidence_record_ids=ev_ids,
            )
        )

    for index, profile in enumerate(selected):
        uid = str(profile.get("profile_uid") or "")
        rec = profile_by_uid.get(uid)
        p_src, p_ev = _ids(rec or summary_rec)
        title = str(profile.get("title") or f"GEO Profile {uid}")
        url = profile.get("profile_url")
        page_break = index > 0
        if url:
            blocks.append(
                ReportContentBlock(
                    kind="link",
                    text=title,
                    links=[{"label": title, "url": str(url)}],
                    presentation_role="section_3a_profile_title",
                    presentation_item_key=f"{item_key}-{uid}-title",
                    presentation_page_break_before=page_break,
                    source_ids=p_src,
                    evidence_record_ids=p_ev,
                )
            )
        else:
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text=title,
                    presentation_role="section_3a_profile_title",
                    presentation_item_key=f"{item_key}-{uid}-title",
                    presentation_page_break_before=page_break,
                    source_ids=p_src,
                    evidence_record_ids=p_ev,
                )
            )

        meta_lines: list[str] = []
        organism = profile.get("organism") or profile.get("taxon")
        if organism:
            meta_lines.append(f"Organism: {organism}")
        reporter = profile.get("reporter_line")
        if reporter:
            meta_lines.append(f"Reporter: {reporter}")
        samples = profile.get("sample_count")
        if samples is not None:
            meta_lines.append(f"Samples: {samples}")
        gse = profile.get("gse")
        if gse:
            meta_lines.append(f"Series: {gse}")
        if profile.get("subset_effect_flag"):
            meta_lines.append("Value subset effect: flagged (ranking signal only)")
        if meta_lines:
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text="\n".join(meta_lines),
                    presentation_role="section_3a_profile_metadata",
                    presentation_item_key=f"{item_key}-{uid}-metadata",
                    source_ids=p_src,
                    evidence_record_ids=p_ev,
                )
            )

        fig_path = None
        fig_val = rec.value if rec and isinstance(rec.value, dict) else {}
        # Prefer an existing on-disk absolute path so Rancho can embed the PNG.
        # ``_resolve_figure_path`` refuses machine-absolute paths for portability,
        # which is correct for audit JSON but not for local HTML/PDF rendering.
        for key in ("local_artifact_path", "relative_path"):
            candidate = fig_val.get(key) or profile.get(key)
            if candidate and Path(str(candidate)).is_file():
                fig_path = str(candidate)
                break
        if not fig_path:
            rel = profile.get("figure_relative_path")
            attempt_dir = str((status.get("audit") or {}).get("gene_attempt_dir") or "")
            if rel and attempt_dir:
                candidate = Path(attempt_dir) / str(rel)
                if candidate.is_file():
                    fig_path = str(candidate)
        if not fig_path and fig_val:
            resolved, fig_diags = _resolve_figure_path(fig_val)
            diagnostics.extend(fig_diags)
            fig_path = resolved
        graph_status = str(profile.get("graph_status") or "")
        if fig_path and profile.get("graph_ok"):
            blocks.append(
                ReportContentBlock(
                    kind="figure",
                    text=f"{gene_symbol} GEO Profile chart ({uid})",
                    figure_path=str(fig_path),
                    presentation_role="section_3a_profile_figure",
                    presentation_item_key=f"{item_key}-{uid}-figure",
                    source_ids=p_src,
                    evidence_record_ids=p_ev,
                )
            )
        elif visual == STATUS_NOT_ATTEMPTED:
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text="Figure acquisition was not attempted for this run.",
                    presentation_role="section_3a_profile_figure_status",
                    presentation_item_key=f"{item_key}-{uid}-figure-status",
                )
            )
        else:
            note = "Profile chart unavailable."
            if graph_status and graph_status not in {"success", "not_attempted_outside_shortlist"}:
                note = f"Profile chart unavailable ({graph_status})."
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text=note,
                    presentation_role="section_3a_profile_figure_status",
                    presentation_item_key=f"{item_key}-{uid}-figure-status",
                )
            )

    caveat = str(summary.get("screening_caveat") or SCREENING_CAVEAT)
    comparability = str(summary.get("comparability_note") or COMPARABILITY_NOTE)
    policy = str(summary.get("selection_policy") or SELECTION_POLICY)
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=f"{caveat} {comparability} {policy}",
            presentation_role="section_3a_caveat",
            presentation_item_key=f"{item_key}-caveat",
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
        )
    )

    if not selected and scientific == STATUS_SUCCESS:
        diagnostics.append(
            PresentationDiagnostic(
                "section_3a",
                "scientific success with no polished selected profiles",
                "warning",
            )
        )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


def build_section_4a_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 4a Harmonizome TF association presentation blocks."""
    from gene_dossier.section_4a import (
        CHEA_BLURB,
        ENCODE_BLURB,
        JASPAR_BLURB,
        MOTIFMAP_BLURB,
        SCIENTIFIC_CAVEAT,
        STATUS_SOURCE_UNAVAILABLE,
        STATUS_SUCCESS,
    )
    from gene_dossier.tools.harmonizome_section4a import gene_page_url

    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    summary = dict(status.get("summary") or {})
    rendering = dict(status.get("rendering_status") or {})
    gene = (
        str(summary.get("official_symbol") or gene_symbol or "").strip() or "this gene"
    )
    item_key = str(
        summary.get("presentation_item_key")
        or f"harmonizome-{gene.lower()}"
    )
    summary_rec = next((r for r in records if r.fact_type == "section_4a_summary"), None)
    supp_rec = next(
        (r for r in records if r.fact_type == "section_4a_supplementary_workbook"),
        None,
    )

    def _ids(rec: EvidenceRecord | None) -> tuple[list[str], list[str]]:
        if rec is None:
            return [], []
        return (
            [rec.source_id] if rec.source_id else [],
            [rec.id] if rec.id else [],
        )

    src_ids, ev_ids = _ids(summary_rec)
    blocks: list[ReportContentBlock] = []
    scientific = str(
        rendering.get("scientific_status") or summary.get("scientific_status") or ""
    )
    gene_url = str(summary.get("gene_page_url") or gene_page_url(gene))

    blurbs = dict(summary.get("source_blurbs") or {})
    source_text = " ".join(
        [
            str(blurbs.get("encode") or ENCODE_BLURB),
            str(blurbs.get("chea") or CHEA_BLURB),
            str(blurbs.get("motifmap") or MOTIFMAP_BLURB),
            str(blurbs.get("jaspar") or JASPAR_BLURB),
        ]
    )

    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=(
                "Transcription factors (TFs) modulate transcription of target genes. "
                "Associations curated in Harmonizome for this gene are summarized below."
            ),
            presentation_role="section_4a_intro",
            presentation_item_key=f"{item_key}-intro",
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
            links=[{"label": "Harmonizome", "url": gene_url}],
        )
    )
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=source_text,
            presentation_role="section_4a_source_description",
            presentation_item_key=f"{item_key}-sources",
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
        )
    )
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=str(summary.get("scientific_caveat") or SCIENTIFIC_CAVEAT),
            presentation_role="section_4a_scientific_caveat",
            presentation_item_key=f"{item_key}-caveat",
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
        )
    )

    if scientific == STATUS_SOURCE_UNAVAILABLE:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="Harmonizome gene associations were unavailable for this run.",
                presentation_role="section_4a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
        diagnostics.append(
            PresentationDiagnostic("section_4a", "source unavailable", "warning")
        )
        return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))

    if scientific in {"gene_mismatch", "no_associations"}:
        msg = (
            "Harmonizome returned a gene symbol that did not match the query."
            if scientific == "gene_mismatch"
            else "No allowlisted transcription-factor associations were found for this gene."
        )
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=msg,
                presentation_role="section_4a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
        return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))

    curated_total = int(summary.get("curated_total") or 0)
    predicted_total = int(summary.get("predicted_total") or 0)
    curated_display = list(summary.get("curated_display") or [])
    predicted_display = list(summary.get("predicted_display") or [])

    curated_count_text = (
        f"{gene} has {curated_total} transcription factor associations."
    )
    if curated_display:
        curated_count_text += " Selected associations are listed below."
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=curated_count_text,
            presentation_role="section_4a_curated_count",
            presentation_item_key=f"{item_key}-curated-count",
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
        )
    )

    curated_headers = [
        "Association",
        "Dataset",
        "Tissue/Cells",
        "Organism",
        "Genome Build",
    ]
    # Safe mid-table page splits with repeated headers (never split supplementary).
    chunk_size = 8
    for chunk_index, start in enumerate(range(0, max(len(curated_display), 0) or 0, chunk_size)):
        chunk = curated_display[start : start + chunk_size]
        if not chunk:
            break
        rows = [
            [
                str(r.get("association") or ""),
                str(r.get("dataset") or ""),
                str(r.get("tissue_cells") or ""),
                str(r.get("organism") or ""),
                str(r.get("genome_build") or ""),
            ]
            for r in chunk
        ]
        blocks.append(
            ReportContentBlock(
                kind="table",
                table_headers=curated_headers,
                table_rows=rows,
                presentation_role="section_4a_curated_table",
                presentation_item_key=f"{item_key}-curated-table-{chunk_index}",
                presentation_page_break_before=chunk_index > 0,
                source_ids=src_ids,
                evidence_record_ids=ev_ids,
            )
        )

    predicted_count_text = (
        f"{gene} has {predicted_total} predicted transcription factor associations."
    )
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=predicted_count_text,
            presentation_role="section_4a_predicted_count",
            presentation_item_key=f"{item_key}-predicted-count",
            presentation_page_break_before=bool(curated_display),
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
        )
    )

    predicted_headers = ["Predicted Association", "Dataset"]
    for chunk_index, start in enumerate(
        range(0, max(len(predicted_display), 0) or 0, chunk_size)
    ):
        chunk = predicted_display[start : start + chunk_size]
        if not chunk:
            break
        rows = [
            [
                str(r.get("predicted_association") or ""),
                str(r.get("dataset") or ""),
            ]
            for r in chunk
        ]
        blocks.append(
            ReportContentBlock(
                kind="table",
                table_headers=predicted_headers,
                table_rows=rows,
                presentation_role="section_4a_predicted_table",
                presentation_item_key=f"{item_key}-predicted-table-{chunk_index}",
                presentation_page_break_before=chunk_index > 0,
                source_ids=src_ids,
                evidence_record_ids=ev_ids,
            )
        )

    xlsx_name = summary.get("supplementary_xlsx") or f"{gene}_Harmonizome.xlsx"
    supp_ids, supp_ev = _ids(supp_rec or summary_rec)
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=(
                "The download of all transcription factor associations is available "
                f"as Supplementary Material ({xlsx_name})."
            ),
            presentation_role="section_4a_supplementary_note",
            presentation_item_key=f"{item_key}-supplementary",
            source_ids=supp_ids,
            evidence_record_ids=supp_ev,
        )
    )

    if scientific == STATUS_SUCCESS and curated_total == 0 and predicted_total == 0:
        diagnostics.append(
            PresentationDiagnostic(
                "section_4a",
                "scientific success with zero allowlisted associations",
                "warning",
            )
        )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


def build_section_5a_blocks(
    *,
    gene_symbol: str,
    evidence_records: Sequence[EvidenceRecord],
    section_status: dict[str, Any] | None = None,
) -> SectionPresentationResult:
    """Build Section 5a STRING network presentation blocks."""
    from gene_dossier.section_5a import (
        SCIENTIFIC_INTRO,
        STATUS_IDENTIFIER_AMBIGUOUS,
        STATUS_IDENTIFIER_NOT_FOUND,
        STATUS_NO_INTERACTIONS,
        STATUS_SOURCE_UNAVAILABLE,
        STATUS_SUCCESS,
        STRING_HOMEPAGE,
    )

    diagnostics: list[PresentationDiagnostic] = []
    records = list(evidence_records)
    status = section_status or {}
    summary = dict(status.get("summary") or {})
    rendering = dict(status.get("rendering_status") or {})
    gene = (
        str(summary.get("preferred_name") or summary.get("official_symbol") or gene_symbol or "")
        .strip()
        or "this gene"
    )
    item_key = str(summary.get("presentation_item_key") or f"string-{gene.lower()}")
    summary_rec = next((r for r in records if r.fact_type == "section_5a_summary"), None)
    fig_rec = next(
        (r for r in records if r.fact_type == "section_5a_network_figure"), None
    )
    supp_rec = next(
        (r for r in records if r.fact_type == "section_5a_supplementary_workbook"),
        None,
    )

    def _ids(rec: EvidenceRecord | None) -> tuple[list[str], list[str]]:
        if rec is None:
            return [], []
        return (
            [rec.source_id] if rec.source_id else [],
            [rec.id] if rec.id else [],
        )

    src_ids, ev_ids = _ids(summary_rec)
    blocks: list[ReportContentBlock] = []
    scientific = str(
        rendering.get("scientific_status") or summary.get("scientific_status") or ""
    )
    visual = str(rendering.get("visual_status") or summary.get("visual_status") or "")
    network_url = str(summary.get("string_network_url") or STRING_HOMEPAGE)

    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=str(summary.get("scientific_intro") or SCIENTIFIC_INTRO),
            presentation_role="section_5a_intro",
            presentation_item_key=f"{item_key}-intro",
            source_ids=src_ids,
            evidence_record_ids=ev_ids,
            links=[{"label": "STRING", "url": network_url}],
        )
    )

    if scientific == STATUS_SOURCE_UNAVAILABLE:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text="STRING network data were unavailable for this run.",
                presentation_role="section_5a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
        diagnostics.append(
            PresentationDiagnostic("section_5a", "source unavailable", "warning")
        )
        return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))

    if scientific == STATUS_IDENTIFIER_NOT_FOUND:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=f"STRING could not resolve an unambiguous protein identifier for {gene}.",
                presentation_role="section_5a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
        return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))

    if scientific == STATUS_IDENTIFIER_AMBIGUOUS:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=f"STRING returned an ambiguous identifier mapping for {gene}.",
                presentation_role="section_5a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
        return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))

    if scientific == STATUS_NO_INTERACTIONS:
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=(
                    f"No STRING functional associations were returned for {gene} "
                    "in the configured bounded network."
                ),
                presentation_role="section_5a_source_status",
                presentation_item_key=f"{item_key}-status",
            )
        )
        return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))

    if scientific == STATUS_SUCCESS:
        n = int(summary.get("unique_node_count") or 0)
        e = int(summary.get("unique_edge_count") or 0)
        d = int(summary.get("direct_query_edge_count") or 0)
        required = int(summary.get("required_score") or 400)
        threshold = required / 1000.0
        summary_text = (
            f"The companion bounded structured STRING network contains {n} proteins "
            f"and {e} STRING associations at a required combined-score threshold of "
            f"{threshold:.3f}. {d} associations directly involve {gene}; "
            f"the remaining {max(e - d, 0)} connect proteins within the displayed "
            "neighborhood."
        )
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=summary_text,
                presentation_role="section_5a_network_summary",
                presentation_item_key=f"{item_key}-network-summary",
                source_ids=src_ids,
                evidence_record_ids=ev_ids,
            )
        )
        blocks.append(
            ReportContentBlock(
                kind="narrative",
                text=(
                    "The official STRING evidence-style visualization below uses a "
                    "10-node first-shell and 20-node second-shell display."
                ),
                presentation_role="section_5a_network_summary",
                presentation_item_key=f"{item_key}-figure-scope",
                source_ids=src_ids,
                evidence_record_ids=ev_ids,
            )
        )

    xlsx_name = summary.get("supplementary_xlsx") or f"{gene}_STRING.xlsx"
    supp_ids, supp_ev = _ids(supp_rec or summary_rec)
    blocks.append(
        ReportContentBlock(
            kind="narrative",
            text=(
                "The complete association inventory for the displayed bounded network "
                f"is available as Supplementary Material ({xlsx_name})."
            ),
            presentation_role="section_5a_supplementary_note",
            presentation_item_key=f"{item_key}-supplementary",
            source_ids=supp_ids,
            evidence_record_ids=supp_ev,
        )
    )

    fig_val = dict(fig_rec.value) if fig_rec and isinstance(fig_rec.value, dict) else {}
    if not fig_val and summary.get("network_figure_local_path"):
        fig_val = {
            "local_artifact_path": summary.get("network_figure_local_path"),
            "relative_path": summary.get("network_figure_relative_path"),
        }
    if visual == STATUS_SUCCESS and fig_val:
        fig_ids, fig_ev = _ids(fig_rec or summary_rec)
        local = fig_val.get("local_artifact_path") or fig_val.get("relative_path")
        fig_path = None
        if local and Path(str(local)).is_file():
            fig_path = str(local)
        else:
            resolved, fig_diags = _resolve_figure_path(fig_val)
            diagnostics.extend(fig_diags)
            fig_path = resolved
        if fig_path:
            blocks.append(
                ReportContentBlock(
                    kind="figure",
                    figure_path=fig_path,
                    figure_caption=f"STRING evidence-style network for {gene}",
                    text=f"STRING evidence-style network for {gene}",
                    presentation_role="section_5a_network_figure",
                    presentation_item_key=f"{item_key}-network-figure",
                    source_ids=fig_ids,
                    evidence_record_ids=fig_ev,
                )
            )
            blocks.append(
                ReportContentBlock(
                    kind="narrative",
                    text=(
                        "Legend: nodes represent proteins (colored first-shell vs white "
                        "second-shell neighbors; empty vs filled markers are "
                        "STRING-provided structure indicators). Edges represent "
                        "functional associations with evidence-style categories: "
                        "curated databases, experimentally determined, gene neighborhood, "
                        "gene fusion, gene co-occurrence, text mining, co-expression, "
                        "and protein homology."
                    ),
                    presentation_role="section_5a_network_legend",
                    presentation_item_key=f"{item_key}-network-legend",
                    source_ids=fig_ids,
                    evidence_record_ids=fig_ev,
                )
            )

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


__all__ = [
    "ALLOWED_LINK_HOSTS",
    "NOT_AVAILABLE",
    "SECTION_1C_BOLD_LEAD_PHRASES",
    "SECTION_1C_CDD_LINK_LABEL",
    "SECTION_1D_CONFIDENCE_BLURB",
    "SECTION_1D_CONFIDENCE_LEGEND_TEXT",
    "UCSC_STABLE_INTRO",
    "PresentationDiagnostic",
    "SectionPresentationResult",
    "build_alphafold_blocks",
    "build_barres_brain_expression_blocks",
    "build_conservation_blocks",
    "build_gene_aliases_blocks",
    "build_homologues_blocks",
    "build_known_structure_blocks",
    "build_section_2c_blocks",
    "build_section_3a_blocks",
    "build_section_4a_blocks",
    "build_section_5a_blocks",
    "build_section_presentation",
    "build_tissue_specific_information_blocks",
    "format_safe_table_cell_html",
    "transcript_selection_sentence",
]


def transcript_selection_sentence(tx_val: dict[str, Any]) -> str:
    """Public alias for Section 1b transcript-selection wording."""
    return _transcript_selection_sentence(tx_val)
