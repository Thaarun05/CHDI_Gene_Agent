"""Full source registry for the Gene Dossier Platform.

Declares every dossier source up front (Priority A / B / C) so the platform never
silently omits a source. The registry is the single map of:

- what sources exist
- their priority and type
- which report sections they support
- which API keys they require
- whether the client / normalizer is implemented yet
- default coverage status (``not_implemented``, ``requires_key``, ``manual``, ``deferred``)

Clients live in ``tools/``; normalizers in ``normalize/``. This module only declares
metadata — it does not call APIs.
"""

from __future__ import annotations

from enum import Enum
from typing import Iterable

from pydantic import BaseModel, Field

from .models import SourceStatus, SourceType


class SourcePriority(str, Enum):
    """Implementation priority for a registered source."""

    A = "A"  # full client + deep normalizer target
    B = "B"  # client + raw storage first; basic normalizer OK
    C = "C"  # scaffold; manual / semi-structured / requires_key / deferred


class SourceDefinition(BaseModel):
    """Metadata for one registered biomedical / dossier source."""

    name: str
    priority: SourcePriority
    source_type: SourceType
    description: str = ""
    # Module path under gene_dossier.tools (e.g. "ncbi_gene") when a client exists.
    client_module: str | None = None
    # Module path under gene_dossier.normalize when a normalizer exists.
    normalizer_module: str | None = None
    # Env var names required before the client can succeed (e.g. "BIOGRID_ACCESSKEY").
    required_keys: list[str] = Field(default_factory=list)
    # Optional keys that improve rate limits / quality but are not required.
    optional_keys: list[str] = Field(default_factory=list)
    # CHDI-style report sections this source can feed.
    report_sections: list[str] = Field(default_factory=list)
    # Default coverage status before a run attempts the source.
    default_status: SourceStatus = SourceStatus.not_implemented
    # True once a tools/ client file exists and is wired.
    client_implemented: bool = False
    # True once a normalize/ module exists and is wired.
    normalizer_implemented: bool = False
    notes: str | None = None


# --------------------------------------------------------------------------------------
# Canonical source map
# --------------------------------------------------------------------------------------
# Status policy at registry time (before clients exist):
# - Priority A/B public APIs: not_implemented (client pending)
# - Sources needing a key with no public fallback: requires_key or not_implemented + notes
# - Semi-structured / no clean API: manual or deferred

_SOURCES: list[SourceDefinition] = [
    # ----- Priority A -----
    SourceDefinition(
        name="NCBI Gene",
        priority=SourcePriority.A,
        source_type=SourceType.curated_database,
        description="Entrez Gene identity, aliases, summary (ESearch + ESummary).",
        client_module="ncbi_gene",
        normalizer_module="gene_identity",
        optional_keys=["NCBI_API_KEY"],
        report_sections=["General gene information", "Gene aliases and identifiers"],
        notes="Prefer exact official-symbol match; SREBF2 Entrez expected 6721.",
    ),
    SourceDefinition(
        name="PubMed",
        priority=SourcePriority.A,
        source_type=SourceType.literature,
        description="Literature search for gene + Huntington Disease MeSH.",
        client_module="pubmed",
        normalizer_module="literature",
        optional_keys=["NCBI_API_KEY"],
        report_sections=["Major labs / literature"],
        notes="Do not assume every hit strongly supports the gene–HD link.",
    ),
    SourceDefinition(
        name="Ensembl",
        priority=SourcePriority.A,
        source_type=SourceType.curated_database,
        description="Gene lookup by symbol; location, biotype, canonical transcript.",
        client_module="ensembl",
        normalizer_module="gene_identity",
        report_sections=["General gene information", "Gene aliases and identifiers"],
        notes="SREBF2 Ensembl expected ENSG00000198911.",
    ),
    SourceDefinition(
        name="UniProt",
        priority=SourcePriority.A,
        source_type=SourceType.curated_database,
        description="Reviewed human protein accession and function annotations.",
        client_module="uniprot",
        normalizer_module="protein",
        report_sections=[
            "General gene information",
            "Known structure / domains",
            "Gene aliases and identifiers",
        ],
        notes="Query gene_exact + organism_id:9606 + reviewed:true; SREBF2 = Q12772.",
    ),
    SourceDefinition(
        name="GTEx",
        priority=SourcePriority.A,
        source_type=SourceType.expression_database,
        description="Tissue median expression and eQTLs (GTEx API v2).",
        client_module="gtex",
        normalizer_module="expression",
        report_sections=["Tissue and cell expression", "eQTLs"],
        notes="Resolve gene first; SREBF2 GTEx GENCODE expected ENSG00000198911.11.",
    ),
    SourceDefinition(
        name="STRING",
        priority=SourcePriority.A,
        source_type=SourceType.interaction_database,
        description="Protein–protein interaction partners.",
        client_module="string_db",
        normalizer_module="ppi",
        report_sections=["Protein-protein interactions"],
        notes="Resolve STRING ID first; caller_identity=gene_dossier_platform.",
    ),
    SourceDefinition(
        name="Reactome",
        priority=SourcePriority.A,
        source_type=SourceType.pathway_database,
        description="Pathway membership via UniProt accession.",
        client_module="reactome",
        normalizer_module="pathways",
        report_sections=["Pathways"],
    ),
    SourceDefinition(
        name="ClinVar",
        priority=SourcePriority.A,
        source_type=SourceType.genetic_database,
        description="Variant / clinical significance summaries for the gene.",
        client_module="clinvar",
        normalizer_module="variants",
        optional_keys=["NCBI_API_KEY"],
        report_sections=["ClinVar / OMIM / Open Targets / SNPs"],
        notes="ESearch gene[gene] AND single_gene[prop]; then ESummary.",
    ),
    SourceDefinition(
        name="Open Targets",
        priority=SourcePriority.A,
        source_type=SourceType.genetic_database,
        description="Disease associations via GraphQL for Ensembl gene ID.",
        client_module="opentargets",
        normalizer_module="variants",
        report_sections=["ClinVar / OMIM / Open Targets / SNPs"],
    ),
    SourceDefinition(
        name="MouseMine",
        priority=SourcePriority.A,
        source_type=SourceType.model_organism_database,
        description="Mouse ortholog, MGI ID, and knockout / phenotype data.",
        client_module="mousemine",
        normalizer_module="model_organisms",
        report_sections=["Homologues", "Knockouts / model phenotypes", "Conservation / orthologs"],
        notes="SREBF2 mouse Entrez 20788; MGI:107585.",
    ),
    SourceDefinition(
        name="CTD",
        priority=SourcePriority.A,
        source_type=SourceType.chemical_database,
        description="Chemical–gene interactions (batch TSV).",
        client_module="ctd",
        normalizer_module="chemicals",
        report_sections=["CTD perturbations", "Chemical tools"],
        notes="Save raw TSV as artifact before normalizing.",
    ),
    SourceDefinition(
        name="ChEMBL",
        priority=SourcePriority.A,
        source_type=SourceType.chemical_database,
        description="Bioactivity / chemical tool annotations for the target.",
        client_module="chembl",
        normalizer_module="chemicals",
        report_sections=["Chemical tools"],
    ),
    SourceDefinition(
        name="PubChem",
        priority=SourcePriority.A,
        source_type=SourceType.chemical_database,
        description="Compound / substance lookups related to chemical tools.",
        client_module="pubchem",
        normalizer_module="chemicals",
        report_sections=["Chemical tools"],
    ),
    SourceDefinition(
        name="NIH RePORTER",
        priority=SourcePriority.A,
        source_type=SourceType.grant_database,
        description="NIH-funded projects mentioning the gene / disease.",
        client_module="nih_reporter",
        normalizer_module="grants",
        report_sections=["NIH/ERC grants"],
    ),
    # ----- Priority B -----
    SourceDefinition(
        name="GEO",
        priority=SourcePriority.B,
        source_type=SourceType.expression_database,
        description="Perturbation / expression series that alter the gene.",
        client_module="geo",
        normalizer_module="perturbation",
        optional_keys=["NCBI_API_KEY"],
        report_sections=["GEO perturbations"],
    ),
    SourceDefinition(
        name="Harmonizome",
        priority=SourcePriority.B,
        source_type=SourceType.expression_database,
        description="Aggregated functional / expression associations.",
        client_module="harmonizome",
        normalizer_module="expression",
        report_sections=["Tissue and cell expression", "Transcription factors"],
    ),
    SourceDefinition(
        name="BioGRID",
        priority=SourcePriority.B,
        source_type=SourceType.interaction_database,
        description="Curated protein and genetic interactions.",
        client_module="biogrid",
        normalizer_module="ppi",
        required_keys=["BIOGRID_ACCESSKEY"],
        report_sections=["Protein-protein interactions"],
        default_status=SourceStatus.requires_key,
        notes="Mark requires_key when BIOGRID_ACCESSKEY is missing.",
    ),
    SourceDefinition(
        name="WikiPathways",
        priority=SourcePriority.B,
        source_type=SourceType.pathway_database,
        description="Community pathway membership.",
        client_module="wikipathways",
        normalizer_module="pathways",
        report_sections=["Pathways"],
    ),
    SourceDefinition(
        name="AlphaFold",
        priority=SourcePriority.B,
        source_type=SourceType.structure_database,
        description="Predicted protein structure metadata / artifacts.",
        client_module="alphafold",
        normalizer_module="protein",
        report_sections=["AlphaFold / PDBe / CDD", "Known structure / domains"],
    ),
    SourceDefinition(
        name="PDBe",
        priority=SourcePriority.B,
        source_type=SourceType.structure_database,
        description="Experimental PDB structures for the protein.",
        client_module="pdbe",
        normalizer_module="protein",
        report_sections=["AlphaFold / PDBe / CDD", "Known structure / domains"],
    ),
    SourceDefinition(
        name="CDD",
        priority=SourcePriority.B,
        source_type=SourceType.structure_database,
        description="Conserved domain annotations (NCBI CDD).",
        client_module="cdd",
        normalizer_module="protein",
        optional_keys=["NCBI_API_KEY"],
        report_sections=["Known structure / domains", "AlphaFold / PDBe / CDD"],
    ),
    SourceDefinition(
        name="NCBI Datasets",
        priority=SourcePriority.B,
        source_type=SourceType.curated_database,
        description="Gene / ortholog dataset packages from NCBI Datasets API.",
        client_module="ncbi_datasets",
        normalizer_module="gene_identity",
        optional_keys=["NCBI_API_KEY"],
        report_sections=["General gene information", "Conservation / orthologs", "Homologues"],
    ),
    SourceDefinition(
        name="UCSC",
        priority=SourcePriority.B,
        source_type=SourceType.curated_database,
        description="Genome browser / conservation context.",
        client_module="ucsc",
        normalizer_module="gene_identity",
        report_sections=["Conservation / orthologs"],
    ),
    # ----- Priority C -----
    SourceDefinition(
        name="Allen Brain Atlas",
        priority=SourcePriority.C,
        source_type=SourceType.expression_database,
        description="Human / mouse brain expression (API or semi-structured).",
        client_module="allen_brain",
        normalizer_module="expression",
        report_sections=["Tissue and cell expression"],
        default_status=SourceStatus.deferred,
        notes="Treat as semi-structured/deferred if no clean API path is available.",
    ),
    SourceDefinition(
        name="BrainRNASeq",
        priority=SourcePriority.C,
        source_type=SourceType.expression_database,
        description="Brain cell-type expression reference.",
        client_module="brainrnaseq",
        normalizer_module="expression",
        report_sections=["Tissue and cell expression"],
        default_status=SourceStatus.deferred,
        notes="Semi-structured / deferred until a stable retrieval path exists.",
    ),
    SourceDefinition(
        name="Patents",
        priority=SourcePriority.C,
        source_type=SourceType.patent_database,
        description="Patent search related to the gene / target.",
        client_module="patents",
        required_keys=["SERPAPI_API_KEY"],
        report_sections=["Patents"],
        default_status=SourceStatus.requires_key,
        notes="Without SERPAPI_API_KEY mark requires_key or manual.",
    ),
    SourceDefinition(
        name="Antibodies",
        priority=SourcePriority.C,
        source_type=SourceType.commercial_source,
        description="Commercial antibody listings (semi-structured / manual).",
        client_module="antibodies",
        report_sections=["Antibodies"],
        default_status=SourceStatus.manual,
        notes="Manual / semi-structured source; do not invent product claims.",
    ),
    SourceDefinition(
        name="OMIM",
        priority=SourcePriority.C,
        source_type=SourceType.genetic_database,
        description="Mendelian disease associations (requires OMIM API key).",
        client_module="omim",
        normalizer_module="variants",
        required_keys=["OMIM_API_KEY"],
        report_sections=["ClinVar / OMIM / Open Targets / SNPs"],
        default_status=SourceStatus.requires_key,
        notes="Mark requires_key when OMIM_API_KEY is missing.",
    ),
    SourceDefinition(
        name="DrugBank",
        priority=SourcePriority.C,
        source_type=SourceType.chemical_database,
        description="Drug / target annotations (license / access constrained).",
        report_sections=["Chemical tools"],
        default_status=SourceStatus.deferred,
        notes="Deferred; license and access path TBD.",
    ),
    SourceDefinition(
        name="NCATS",
        priority=SourcePriority.C,
        source_type=SourceType.chemical_database,
        description="NCATS / translational chemical resources.",
        report_sections=["Chemical tools"],
        default_status=SourceStatus.deferred,
        notes="Deferred pending endpoint selection.",
    ),
    SourceDefinition(
        name="ERC Grants",
        priority=SourcePriority.C,
        source_type=SourceType.grant_database,
        description="European Research Council grant mentions (manual / deferred).",
        report_sections=["NIH/ERC grants"],
        default_status=SourceStatus.deferred,
        notes="No clean public API path yet; record as deferred/manual.",
    ),
    # Future HD-specific integration placeholder (architecture note only).
    SourceDefinition(
        name="HDinHD",
        priority=SourcePriority.C,
        source_type=SourceType.hd_specific_database,
        description="Huntington's disease-specific resources (future MCP integration).",
        report_sections=["Missing / deferred / manual sources"],
        default_status=SourceStatus.deferred,
        notes="TODO: HDinHD MCP integration — not implemented in MVP.",
    ),
]


def get_all_sources() -> list[SourceDefinition]:
    """Return a copy of every registered source definition."""
    return [s.model_copy() for s in _SOURCES]


def get_source(name: str) -> SourceDefinition | None:
    """Return the source definition matching ``name`` (case-insensitive), or None."""
    needle = name.strip().lower()
    for src in _SOURCES:
        if src.name.lower() == needle:
            return src.model_copy()
    return None


def get_sources_by_priority(priority: SourcePriority | str) -> list[SourceDefinition]:
    """Return all sources with the given priority (``A`` / ``B`` / ``C``)."""
    if isinstance(priority, str):
        priority = SourcePriority(priority.upper())
    return [s.model_copy() for s in _SOURCES if s.priority == priority]


def list_source_names() -> list[str]:
    """Return ordered source names as declared in the registry."""
    return [s.name for s in _SOURCES]


def iter_sources(
    *,
    priority: SourcePriority | None = None,
    implemented_clients_only: bool = False,
) -> Iterable[SourceDefinition]:
    """Iterate sources, optionally filtered by priority or client implementation."""
    for src in _SOURCES:
        if priority is not None and src.priority != priority:
            continue
        if implemented_clients_only and not src.client_implemented:
            continue
        yield src.model_copy()


def registry_summary() -> dict[str, int]:
    """Return counts by priority and default status (useful for smoke checks)."""
    by_priority = {p.value: 0 for p in SourcePriority}
    by_status: dict[str, int] = {}
    for src in _SOURCES:
        by_priority[src.priority.value] += 1
        key = src.default_status.value
        by_status[key] = by_status.get(key, 0) + 1
    return {"total": len(_SOURCES), "by_priority": by_priority, "by_default_status": by_status}


__all__ = [
    "SourcePriority",
    "SourceDefinition",
    "get_all_sources",
    "get_source",
    "get_sources_by_priority",
    "list_source_names",
    "iter_sources",
    "registry_summary",
]
