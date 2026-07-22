"""Rancho BioSciences / CHDI gene-report schema and layout mapping.

The final dossier must match ``SREBF2_report.pdf`` (cover, TOC, 15 major
sections with lettered subsections, References, Compiled databases).

This module is the **layout contract**:

- defines the exact section tree and visual tokens from the reference PDF
- maps :class:`~gene_dossier.models.EvidenceRecord` objects into those slots
- builds a :class:`ReportDocument` for a polished renderer (HTML/PDF)

The markdown assembler in ``rendering.py`` remains a provenance/debug view only.
It is **not** the final CHDI/Rancho report format.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Iterable, Literal

from pydantic import BaseModel, Field

from gene_dossier.models import AssertionType, EvidenceRecord

# --------------------------------------------------------------------------------------
# Visual tokens sampled from SREBF2_report.pdf (Rancho / CHDI)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ReportStyle:
    """Brand colors and type sizes matching the reference PDF."""

    # Cover / major headers: rgb(122, 193, 67)
    green_major: str = "#7AC143"
    # Subsection headers / date: rgb(245, 130, 32)
    orange_sub: str = "#F58220"
    # Body / prepared-for text: rgb(89, 31, 0)
    brown_body: str = "#591F00"
    # Hyperlinks / IDs: rgb(246, 100, 0)
    orange_link: str = "#F66400"
    # Table header fill approx rgb(201, 230, 179)
    table_header_bg: str = "#C9E6B3"
    # Thin header rule under Rancho logo
    rule_green: str = "#B7DFA0"

    cover_title_pt: float = 48.0
    cover_subtitle_pt: float = 36.0
    major_header_pt: float = 20.0
    subsection_header_pt: float = 16.0
    body_pt: float = 12.0

    footer_url: str = "www.RanchoBioSciences.com"
    prepared_for: str = "Prepared for the CHDI Foundation"
    prepared_by: str = "by Rancho BioSciences"


REPORT_STYLE = ReportStyle()

BlockKind = Literal[
    "narrative",
    "table",
    "figure",
    "list",
    "link",
    "empty",
]


# --------------------------------------------------------------------------------------
# Canonical 15-section tree (exact TOC wording from SREBF2_report.pdf)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SubsectionSpec:
    """One lettered subsection under a major section."""

    key: str  # "a", "b", ...
    title: str  # display title as on body pages
    toc_title: str  # ALL-CAPS form used in the table of contents


@dataclass(frozen=True)
class MajorSectionSpec:
    """One of the 15 top-level Rancho/CHDI report sections."""

    number: int
    key: str  # "1" .. "15"
    title: str
    toc_title: str
    subsections: tuple[SubsectionSpec, ...] = ()


def _sub(key: str, title: str, toc_title: str | None = None) -> SubsectionSpec:
    return SubsectionSpec(key=key, title=title, toc_title=(toc_title or title).upper())


REPORT_SECTIONS: tuple[MajorSectionSpec, ...] = (
    MajorSectionSpec(
        number=1,
        key="1",
        title="General Gene Information",
        toc_title="GENERAL GENE INFORMATION",
        subsections=(
            _sub("a", "Gene Aliases", "GENE ALIASES"),
            _sub(
                "b",
                "Conservation of gene across species and conservation of region "
                "under peak in other genomes",
                "CONSERVATION OF GENE ACROSS SPECIES AND CONSERVATION OF REGION "
                "UNDER PEAK IN OTHER GENOMES",
            ),
            _sub("c", "Known structure", "KNOWN STRUCTURE"),
            _sub(
                "d",
                "AlphaFold protein structure prediction",
                "ALPHAFOLD PROTEIN STRUCTURE PREDICTION",
            ),
            _sub(
                "e",
                "Homologues in model animals",
                "HOMOLOGUES IN MODEL ANIMALS",
            ),
        ),
    ),
    MajorSectionSpec(
        number=2,
        key="2",
        title="Expression pattern by cell and tissue",
        toc_title="EXPRESSION PATTERN BY CELL AND TISSUE",
        subsections=(
            _sub("a", "Tissue-specific information", "TISSUE-SPECIFIC INFORMATION"),
            _sub(
                "b",
                "Barres Lab RNA-Seq brain specific expression data",
                "BARRES LAB RNA-SEQ BRAIN SPECIFIC EXPRESSION DATA",
            ),
            _sub(
                "c",
                "snRNA-Seq gene expression in cell type database",
                "SNRNA-SEQ GENE EXPRESSION IN CELL TYPE DATABASE",
            ),
        ),
    ),
    MajorSectionSpec(
        number=3,
        key="3",
        title="Perturbations in GEO that alter the gene",
        toc_title="PERTURBATIONS IN GEO THAT ALTER THE GENE",
        subsections=(
            _sub(
                "a",
                "GEO Profiles search focusing on brain and/or neurons",
                "GEO PROFILES SEARCH FOCUSING ON BRAIN AND/OR NEURONS",
            ),
        ),
    ),
    MajorSectionSpec(
        number=4,
        key="4",
        title="Transcription factors that drive the gene’s expression",
        toc_title="TRANSCRIPTION FACTORS THAT DRIVE THE GENE’S EXPRESSION",
        subsections=(
            _sub(
                "a",
                "Harmonizome integrated knowledge about genes & proteins",
                "HARMONIZOME INTEGRATED KNOWLEDGE ABOUT GENES & PROTEINS",
            ),
        ),
    ),
    MajorSectionSpec(
        number=5,
        key="5",
        title="Protein-protein interaction (PPI) partners",
        toc_title="PROTEIN-PROTEIN INTERACTION (PPI) PARTNERS",
        subsections=(
            _sub("a", "STRING", "STRING"),
            _sub("b", "BioGRID", "BIOGRID"),
        ),
    ),
    MajorSectionSpec(
        number=6,
        key="6",
        title="Information on which perturbations affect the gene",
        toc_title="INFORMATION ON WHICH PERTURBATIONS AFFECT THE GENE",
        subsections=(
            _sub(
                "a",
                "Comparative Toxicogenomics Database",
                "COMPARATIVE TOXICOGENOMICS DATABASE",
            ),
        ),
    ),
    MajorSectionSpec(
        number=7,
        key="7",
        title="Available chemical tools and their effects on the gene",
        toc_title="AVAILABLE CHEMICAL TOOLS AND THEIR EFFECTS ON THE GENE",
        subsections=(
            _sub(
                "a",
                "Queried in the following databases to identify small molecule inhibitors",
                "QUERIED IN THE FOLLOWING DATABASES TO IDENTIFY SMALL MOLECULE "
                "INHIBITORS",
            ),
            _sub(
                "b",
                "Tractability assessments for small molecule",
                "TRACTABILITY ASSESSMENTS FOR SMALL MOLECULE",
            ),
        ),
    ),
    MajorSectionSpec(
        number=8,
        key="8",
        title="eQTLs that alter the gene expression in the brain tissue",
        toc_title="EQTLS THAT ALTER THE GENE EXPRESSION IN THE BRAIN TISSUE",
        subsections=(
            _sub("a", "GTEx eQTL Browser", "GTEX EQTL BROWSER"),
        ),
    ),
    MajorSectionSpec(
        number=9,
        key="9",
        title="Known SNPs linked to function/disease",
        toc_title="KNOWN SNPS LINKED TO FUNCTION/DISEASE",
        subsections=(
            _sub("a", "ClinVar", "CLINVAR"),
            _sub("b", "OMIM", "OMIM"),
            _sub("c", "Open Targets", "OPEN TARGETS"),
            # Reference TOC lists SNPs3D after Open Targets.
            _sub("d", "SNPs3D", "SNPS3D"),
        ),
    ),
    MajorSectionSpec(
        number=10,
        key="10",
        title="Major pathways",
        toc_title="MAJOR PATHWAYS",
        subsections=(
            _sub(
                "a",
                "Pathways associated with the gene",
                "PATHWAYS ASSOCIATED WITH THE GENE",
            ),
        ),
    ),
    MajorSectionSpec(
        number=11,
        key="11",
        title="Available knockouts of the gene and their phenotype",
        toc_title="AVAILABLE KNOCKOUTS OF THE GENE AND THEIR PHENOTYPE",
        subsections=(
            _sub("a", "PubMed", "PUBMED"),
            _sub(
                "b",
                "Mouse Genome Informatics (MGI)",
                "MOUSE GENOME INFORMATICS (MGI)",
            ),
        ),
    ),
    MajorSectionSpec(
        number=12,
        key="12",
        title="Major labs working on the gene based on publications",
        toc_title="MAJOR LABS WORKING ON THE GENE BASED ON PUBLICATIONS",
        subsections=(),
    ),
    MajorSectionSpec(
        number=13,
        key="13",
        title="Commercial antibodies",
        toc_title="COMMERCIAL ANTIBODIES",
        subsections=(),
    ),
    MajorSectionSpec(
        number=14,
        key="14",
        title="Information on Patents",
        toc_title="INFORMATION ON PATENTS",
        subsections=(),
    ),
    MajorSectionSpec(
        number=15,
        key="15",
        title="NIH and ERC grants",
        toc_title="NIH AND ERC GRANTS",
        subsections=(
            _sub("a", "NIH RePORTER", "NIH REPORTER"),
            _sub("b", "ERC", "ERC"),
        ),
    ),
)

BACK_MATTER_TITLES: tuple[str, ...] = (
    "References",
    "Compiled List of Relevant Databases",
)

# Exact compiled database list from SREBF2_report.pdf pp. 38-39.
COMPILED_RELEVANT_DATABASES: tuple[tuple[str, str], ...] = (
    ("GeneCards", "https://www.genecards.org/"),
    ("USCS Genome Browser", "https://genome.ucsc.edu/"),
    ("NCBI Structure Group", "https://www.ncbi.nlm.nih.gov/Structure/index.shtml"),
    ("The RCSB Protein Data Bank", "https://www.rcsb.org/"),
    ("AlphaFold Protein Structure Database", "https://www.alphafold.ebi.ac.uk/"),
    (
        "NCBI Annotation Pipeline",
        "https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/",
    ),
    ("GTEx Portal", "https://gtexportal.org/home/"),
    ("The Human Brain Transcriptome (HBT) Project", "http://hbatlas.org/"),
    ("Allen Brain Map", "http://portal.brain-map.org/"),
    ("Barres Lab Brain RNA-Seq Portal", "http://www.brainrnaseq.org/"),
    ("DropViz", "http://dropviz.org/"),
    ("NCBI GEO Profiles", "https://www.ncbi.nlm.nih.gov/geoprofiles/"),
    ("Harmonizome", "https://maayanlab.cloud/Harmonizome/"),
    ("STRING", "https://string-db.org/"),
    ("BioGRID", "https://thebiogrid.org/"),
    ("Comparative Toxicogenomics Database", "http://ctdbase.org/"),
    ("Chembl", "https://www.ebi.ac.uk/chembl/"),
    ("DrugBank", "https://www.drugbank.ca/"),
    ("PubMed", "https://www.ncbi.nlm.nih.gov/pubmed/"),
    ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/"),
    ("NCATS Inxight", "https://drugs.ncats.io/about"),
    ("NCBI ClinVar", "https://www.ncbi.nlm.nih.gov/clinvar/"),
    ("Online Mendelian Inheritance in Man (OMIM)", "https://www.omim.org/"),
    ("Open Targets", "https://www.targetvalidation.org/"),
    ("Pathway Commons", "https://apps.pathwaycommons.org/"),
    ("Reactome", "https://reactome.org/"),
    ("Mouse Genomic Informatics (MGI)", "http://www.informatics.jax.org/"),
    (
        "NIH Research Portfolio Online Reporting Tools (RePORT)",
        "https://report.nih.gov/",
    ),
    ("European Research Council (ERC)", "https://erc.europa.eu/homepage"),
)


def get_major_section(number_or_key: int | str) -> MajorSectionSpec:
    """Look up a major section by number (1-15) or key string."""
    key = str(number_or_key)
    for spec in REPORT_SECTIONS:
        if spec.key == key or str(spec.number) == key:
            return spec
    raise KeyError(f"Unknown major section: {number_or_key!r}")


def iter_toc_entries() -> list[dict[str, Any]]:
    """Flattened TOC entries for rendering (majors + lettered subs + back matter)."""
    entries: list[dict[str, Any]] = []
    for major in REPORT_SECTIONS:
        entries.append(
            {
                "kind": "major",
                "number": major.number,
                "key": major.key,
                "title": major.toc_title,
                "display_title": major.title,
            }
        )
        for sub in major.subsections:
            entries.append(
                {
                    "kind": "subsection",
                    "number": major.number,
                    "major_key": major.key,
                    "key": sub.key,
                    "title": sub.toc_title,
                    "display_title": sub.title,
                    "slot": f"{major.key}{sub.key}",
                }
            )
    for title in BACK_MATTER_TITLES:
        entries.append(
            {
                "kind": "back_matter",
                "title": title.upper(),
                "display_title": title,
            }
        )
    return entries


# --------------------------------------------------------------------------------------
# Evidence -> slot mapping
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ReportSlot:
    """Address of a content cell in the Rancho layout."""

    major_key: str
    subsection_key: str | None = None

    @property
    def slot_id(self) -> str:
        if self.subsection_key:
            return f"{self.major_key}{self.subsection_key}"
        return self.major_key


def _norm_source(name: str) -> str:
    return " ".join((name or "").strip().lower().split())


def resolve_report_slot(record: EvidenceRecord) -> ReportSlot | None:
    """Map one evidence record into a Rancho major/subsection slot.

    Returns ``None`` for material that does not belong in the polished 15-section
    body (those stay in the provenance/debug view).
    """
    source = _norm_source(record.source_name)
    assertion = record.assertion_type
    section = (record.section or "").strip().lower()
    subsection = (record.subsection or "").strip().lower()
    fact = (record.fact_type or "").strip().lower()

    # --- Section 1: General Gene Information ---
    if source in {"ncbi gene", "ensembl"} and assertion is AssertionType.gene_identity:
        return ReportSlot("1", "a")
    if source == "uniprot" and assertion is AssertionType.gene_identity:
        return ReportSlot("1", "a")
    if source == "ucsc":
        return ReportSlot("1", "b")
    if source in {"cdd", "pdbe"}:
        return ReportSlot("1", "c")
    if source == "uniprot" and assertion in {
        AssertionType.protein_function,
        AssertionType.protein_structure,
        AssertionType.disease_association,
    }:
        return ReportSlot("1", "c")
    if source == "alphafold" or "alphafold" in section:
        return ReportSlot("1", "d")
    if source in {"ncbi datasets", "mousemine", "mouse mine"} and (
        assertion is AssertionType.orthology
    ):
        return ReportSlot("1", "e")
    if assertion is AssertionType.orthology and "homolog" in section:
        return ReportSlot("1", "e")

    # --- Section 2: Expression (GTEx eQTLs are section 8 — check first) ---
    if source == "gtex" and (
        "eqtl" in section or "eqtl" in fact or "eqtl" in subsection
    ):
        return ReportSlot("8", "a")
    if source == "gtex" and assertion is AssertionType.expression:
        return ReportSlot("2", "a")
    if assertion is AssertionType.cell_type_expression:
        return ReportSlot("2", "c")
    if source in {
        "brainrnaseq",
        "brain rnaseq",
        "brain rna-seq",
        "barres lab brain rna-seq",
        "allen",
        "allen brain",
        "allen hba",
        "allen brain atlas",
        "human brain transcriptome",
        "hbt",
        "dropviz",
    }:
        if (
            source == "dropviz"
            or "cell" in section
            or "snrna" in subsection
            or "single" in subsection
        ):
            return ReportSlot("2", "c")
        return ReportSlot("2", "b")
    if source == "harmonizome" and assertion is AssertionType.expression:
        return ReportSlot("2", "a")

    # --- Section 3: GEO perturbations ---
    if source == "geo" or (
        assertion is AssertionType.perturbation and "geo" in section
    ):
        return ReportSlot("3", "a")
    if assertion is AssertionType.perturbation and source == "geo":
        return ReportSlot("3", "a")

    # --- Section 4: Transcription factors ---
    if assertion is AssertionType.transcription_factor_association or (
        source == "harmonizome" and "transcription" in section
    ):
        return ReportSlot("4", "a")

    # --- Section 5: PPI ---
    if assertion is AssertionType.ppi or "protein-protein" in section:
        if source == "biogrid":
            return ReportSlot("5", "b")
        return ReportSlot("5", "a")

    # --- Section 6: CTD ---
    if source == "ctd" or (
        assertion is AssertionType.chemical_interaction and "ctd" in section
    ):
        return ReportSlot("6", "a")

    # --- Section 7: Chemical tools ---
    if source in {"chembl", "pubchem"}:
        return ReportSlot("7", "a")
    if assertion is AssertionType.chemical_interaction and source != "ctd":
        if "tractab" in fact or "tractability" in subsection:
            return ReportSlot("7", "b")
        return ReportSlot("7", "a")
    if source == "open targets":
        if "tractab" in fact or "tractability" in subsection or "tractab" in subsection:
            return ReportSlot("7", "b")
        if "probe" in fact or "chemical" in section:
            return ReportSlot("7", "a")

    # --- Section 8: eQTLs (non-GTEx fallback; GTEx handled above) ---
    if "eqtl" in section:
        return ReportSlot("8", "a")

    # --- Section 9: SNPs / disease ---
    if source == "clinvar" or (
        assertion is AssertionType.variant_association and source != "omim"
    ):
        return ReportSlot("9", "a")
    if source == "omim":
        return ReportSlot("9", "b")
    if source == "open targets" and assertion is AssertionType.disease_association:
        return ReportSlot("9", "c")
    if "snp" in source or source == "snps3d":
        return ReportSlot("9", "d")

    # --- Section 10: Pathways ---
    if assertion is AssertionType.pathway_membership or source in {
        "reactome",
        "wikipathways",
    }:
        return ReportSlot("10", "a")

    # --- Section 11: Knockouts ---
    if source in {"mousemine", "mouse mine", "mgi"} or (
        assertion is AssertionType.knockout_phenotype
    ):
        return ReportSlot("11", "b")
    if source == "pubmed" and (
        "knockout" in section or "phenotype" in section
    ):
        return ReportSlot("11", "a")

    # --- Section 13/14 before the broad literature catch-all ---
    if source == "antibodies" or "antibody" in section:
        return ReportSlot("13", None)
    if source == "patents" or assertion is AssertionType.patent_claim or "patent" in section:
        return ReportSlot("14", None)

    # --- Section 12: Major labs / literature ---
    if source == "pubmed" or assertion is AssertionType.literature_summary:
        return ReportSlot("12", None)
    if "literature" in section or "labs" in section:
        return ReportSlot("12", None)

    # --- Section 15: Grants ---
    if source == "erc" or (
        assertion is AssertionType.grant_project and "erc" in source
    ):
        return ReportSlot("15", "b")
    if source in {"nih reporter", "nih_reporter"} or (
        assertion is AssertionType.grant_project
    ):
        return ReportSlot("15", "a")

    # Coarse section-name fallback from normalizer section labels.
    section_fallbacks: list[tuple[str, ReportSlot]] = [
        ("general gene information", ReportSlot("1", "a")),
        ("gene aliases", ReportSlot("1", "a")),
        ("conservation", ReportSlot("1", "b")),
        ("known structure", ReportSlot("1", "c")),
        ("alphafold", ReportSlot("1", "d")),
        ("homologues", ReportSlot("1", "e")),
        ("tissue and cell expression", ReportSlot("2", "a")),
        ("geo perturbation", ReportSlot("3", "a")),
        ("transcription factor", ReportSlot("4", "a")),
        ("protein-protein", ReportSlot("5", "a")),
        ("ctd perturbation", ReportSlot("6", "a")),
        ("chemical tool", ReportSlot("7", "a")),
        ("eqtl", ReportSlot("8", "a")),
        ("clinvar", ReportSlot("9", "a")),
        ("pathway", ReportSlot("10", "a")),
        ("knockout", ReportSlot("11", "b")),
        ("major labs", ReportSlot("12", None)),
        ("literature", ReportSlot("12", None)),
        ("antibody", ReportSlot("13", None)),
        ("patent", ReportSlot("14", None)),
        ("grant", ReportSlot("15", "a")),
        ("nih", ReportSlot("15", "a")),
    ]
    for needle, slot in section_fallbacks:
        if needle in section:
            return slot

    return None


# --------------------------------------------------------------------------------------
# Document models
# --------------------------------------------------------------------------------------
class ReportCover(BaseModel):
    """Cover page fields matching the reference PDF."""

    gene_symbol: str
    chromosome: str | None = None
    title_line: str = "Gene Report"
    prepared_for: str = REPORT_STYLE.prepared_for
    prepared_by: str = REPORT_STYLE.prepared_by
    curator: str | None = None
    report_date: str | None = None

    @property
    def gene_line(self) -> str:
        if self.chromosome:
            chrom = self.chromosome.upper().removeprefix("CHR")
            return f"{self.gene_symbol} (CHR{chrom})"
        return self.gene_symbol


class ReportContentBlock(BaseModel):
    """One visual block inside a subsection (narrative, table, figure, ...)."""

    kind: BlockKind = "narrative"
    title: str | None = None
    text: str | None = None
    table_headers: list[str] = Field(default_factory=list)
    table_rows: list[list[str]] = Field(default_factory=list)
    figure_path: str | None = None
    figure_caption: str | None = None
    links: list[dict[str, str]] = Field(default_factory=list)
    source_ids: list[str] = Field(default_factory=list)
    evidence_record_ids: list[str] = Field(default_factory=list)


class ReportSubsection(BaseModel):
    """Lettered subsection with ordered content blocks + provenance ids."""

    key: str
    title: str
    toc_title: str
    blocks: list[ReportContentBlock] = Field(default_factory=list)
    source_ids: list[str] = Field(default_factory=list)
    evidence_record_ids: list[str] = Field(default_factory=list)
    status: Literal["populated", "empty"] = "empty"


class ReportMajorSection(BaseModel):
    """Numbered major section (1-15)."""

    number: int
    key: str
    title: str
    toc_title: str
    subsections: list[ReportSubsection] = Field(default_factory=list)
    blocks: list[ReportContentBlock] = Field(default_factory=list)
    source_ids: list[str] = Field(default_factory=list)
    evidence_record_ids: list[str] = Field(default_factory=list)
    status: Literal["populated", "empty", "partial"] = "empty"


class ReportDocument(BaseModel):
    """Full Rancho/CHDI dossier layout instance (not the debug markdown view)."""

    dossier_run_id: str
    gene_symbol: str
    cover: ReportCover
    sections: list[ReportMajorSection] = Field(default_factory=list)
    references: list[str] = Field(default_factory=list)
    compiled_databases: list[dict[str, str]] = Field(default_factory=list)
    unmapped_source_ids: list[str] = Field(default_factory=list)
    style: dict[str, Any] = Field(
        default_factory=lambda: REPORT_STYLE.__dict__.copy()
    )


def _empty_document_skeleton(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    cover: ReportCover,
) -> ReportDocument:
    sections: list[ReportMajorSection] = []
    for spec in REPORT_SECTIONS:
        subs = [
            ReportSubsection(
                key=sub.key,
                title=sub.title,
                toc_title=sub.toc_title,
                status="empty",
            )
            for sub in spec.subsections
        ]
        sections.append(
            ReportMajorSection(
                number=spec.number,
                key=spec.key,
                title=spec.title,
                toc_title=spec.toc_title,
                subsections=subs,
                status="empty",
            )
        )
    return ReportDocument(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        cover=cover,
        sections=sections,
        compiled_databases=[
            {"name": name, "url": url} for name, url in COMPILED_RELEVANT_DATABASES
        ],
    )


def _block_from_evidence(record: EvidenceRecord) -> ReportContentBlock:
    return ReportContentBlock(
        kind="narrative",
        text=(record.display_text or "").strip() or None,
        source_ids=[record.source_id] if record.source_id else [],
        evidence_record_ids=[record.id] if record.id else [],
    )


def infer_chromosome(evidence_records: Iterable[EvidenceRecord]) -> str | None:
    """Best-effort chromosome from evidence values / display text (e.g. chr22)."""
    pattern = re.compile(r"\bchr(?:omosome)?[\s_:-]*([0-9XYM]+)\b", re.IGNORECASE)
    for record in evidence_records:
        for blob in (record.display_text or "", str(record.value or "")):
            match = pattern.search(blob)
            if match:
                return match.group(1).upper()
        if isinstance(record.value, dict):
            for key in ("chromosome", "chrom", "chr"):
                if record.value.get(key):
                    return str(record.value[key]).upper().removeprefix("CHR")
    return None


def build_report_document(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    evidence_records: Iterable[EvidenceRecord],
    curator: str | None = None,
    report_date: str | None = None,
    chromosome: str | None = None,
    references: Iterable[str] | None = None,
) -> ReportDocument:
    """Bucket evidence into the Rancho 15-section layout.

    Does not invent narrative beyond ``display_text``. A later renderer turns this
    document into the polished PDF/HTML; ``source_id``s stay on every block for the
    JSON/provenance sidecar and optional endnotes.
    """
    records = list(evidence_records)
    cover = ReportCover(
        gene_symbol=gene_symbol,
        chromosome=chromosome or infer_chromosome(records),
        curator=curator,
        report_date=report_date,
    )
    doc = _empty_document_skeleton(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        cover=cover,
    )
    unmapped: list[str] = []
    slotted: dict[str, list[EvidenceRecord]] = defaultdict(list)

    for record in records:
        slot = resolve_report_slot(record)
        if slot is None:
            if record.source_id:
                unmapped.append(record.source_id)
            continue
        slotted[slot.slot_id].append(record)

    for major in doc.sections:
        if major.subsections:
            any_pop = False
            all_pop = True
            for sub in major.subsections:
                recs = slotted.get(f"{major.key}{sub.key}", [])
                if not recs:
                    all_pop = False
                    continue
                any_pop = True
                for rec in recs:
                    block = _block_from_evidence(rec)
                    sub.blocks.append(block)
                    if rec.source_id and rec.source_id not in sub.source_ids:
                        sub.source_ids.append(rec.source_id)
                    if rec.id and rec.id not in sub.evidence_record_ids:
                        sub.evidence_record_ids.append(rec.id)
                    if rec.source_id and rec.source_id not in major.source_ids:
                        major.source_ids.append(rec.source_id)
                    if rec.id and rec.id not in major.evidence_record_ids:
                        major.evidence_record_ids.append(rec.id)
                sub.status = "populated"
            if any_pop and all_pop:
                major.status = "populated"
            elif any_pop:
                major.status = "partial"
            else:
                major.status = "empty"
        else:
            recs = slotted.get(major.key, [])
            if not recs:
                major.status = "empty"
                continue
            for rec in recs:
                block = _block_from_evidence(rec)
                major.blocks.append(block)
                if rec.source_id and rec.source_id not in major.source_ids:
                    major.source_ids.append(rec.source_id)
                if rec.id and rec.id not in major.evidence_record_ids:
                    major.evidence_record_ids.append(rec.id)
            major.status = "populated"

    doc.unmapped_source_ids = list(dict.fromkeys(unmapped))
    if references is not None:
        doc.references = [str(r) for r in references]
    return doc


def cover_lines(cover: ReportCover) -> list[str]:
    """Ordered cover text lines for renderers (logos added by the visual renderer)."""
    lines = [
        cover.gene_line,
        cover.title_line,
        "",
        cover.prepared_for,
        cover.prepared_by,
    ]
    if cover.curator:
        lines.append(f"Curator: {cover.curator}")
    if cover.report_date:
        lines.append(cover.report_date)
    return lines


__all__ = [
    "REPORT_STYLE",
    "ReportStyle",
    "REPORT_SECTIONS",
    "BACK_MATTER_TITLES",
    "COMPILED_RELEVANT_DATABASES",
    "SubsectionSpec",
    "MajorSectionSpec",
    "ReportSlot",
    "ReportCover",
    "ReportContentBlock",
    "ReportSubsection",
    "ReportMajorSection",
    "ReportDocument",
    "get_major_section",
    "iter_toc_entries",
    "resolve_report_slot",
    "infer_chromosome",
    "build_report_document",
    "cover_lines",
]
