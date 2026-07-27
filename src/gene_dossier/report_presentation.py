"""Polished Rancho/CHDI section presentation builders.

Separates human-facing presentation blocks from audit evidence blocks.
Currently implements Section 1a (Gene Aliases) and Section 1b (UCSC conservation).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
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
    return str(path), diags


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

    dynamic = (
        f"A query of the {assembly} assembly identified "
        f"{exact_count if exact_count is not None else 'an unknown number of'} "
        f"{gene} transcript models in the current "
        f"{release or 'GENCODE'} GENCODE annotation within the selected locus. "
        f"The MANE Select and Ensembl canonical transcript was selected for display."
    )
    # Avoid double "GENCODE GENCODE" when release already includes word
    if release and str(release).upper().startswith("GENCODE"):
        dynamic = (
            f"A query of the {assembly} assembly identified "
            f"{exact_count if exact_count is not None else 'an unknown number of'} "
            f"{gene} transcript models in the current {release} annotation "
            f"within the selected locus. The MANE Select and Ensembl canonical "
            f"transcript was selected for display."
        )
    elif release and str(release).upper().startswith("V"):
        dynamic = (
            f"A query of the {assembly} assembly identified "
            f"{exact_count if exact_count is not None else 'an unknown number of'} "
            f"{gene} transcript models in the current GENCODE {release} annotation "
            f"within the selected locus. The MANE Select and Ensembl canonical "
            f"transcript was selected for display."
        )

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
            text=UCSC_STABLE_INTRO,
            source_ids=source_ids,
            evidence_record_ids=evidence_ids,
        ),
        ReportContentBlock(
            kind="narrative",
            text=dynamic,
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
            caption = fig_val.get("caption") or fig_val.get("source_note")
            if caption:
                blocks.append(
                    ReportContentBlock(
                        kind="narrative",
                        text=str(caption),
                        source_ids=source_ids,
                        evidence_record_ids=evidence_ids,
                    )
                )

    # Need at least intro + dynamic to consider presentation nonempty for unknown genes
    if locus is None and inventory is None and transcript is None:
        return SectionPresentationResult(blocks=(), diagnostics=tuple(diagnostics))

    return SectionPresentationResult(blocks=tuple(blocks), diagnostics=tuple(diagnostics))


__all__ = [
    "ALLOWED_LINK_HOSTS",
    "NOT_AVAILABLE",
    "UCSC_STABLE_INTRO",
    "PresentationDiagnostic",
    "SectionPresentationResult",
    "build_conservation_blocks",
    "build_gene_aliases_blocks",
    "build_section_presentation",
    "format_safe_table_cell_html",
]
