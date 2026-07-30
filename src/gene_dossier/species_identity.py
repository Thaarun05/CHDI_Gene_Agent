"""Species-specific gene identity resolution (Human / Mouse / Rat).

Resolves NCBI Gene, Ensembl, and UniProt identity for each target taxon.
Does not hardcode gene biology; callers supply a gene symbol.

Workflow:
1. NCBI Gene search + ESummary selection for the taxon.
2. Resolve the species-specific canonical symbol from NCBI (or fall back).
3. Ensembl lookup/symbol with that resolved symbol.
4. UniProt reviewed search with that resolved symbol + taxon ID.
   If Swiss-Prot returns no accession, fall back to an unreviewed
   (TrEMBL) UniProtKB search for the same symbol + taxon.

ToolResult.gene_symbol always retains the dossier query symbol. The
species-specific symbol is carried in request_params / payload metadata as
``resolved_symbol`` for normalizers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult
from gene_dossier.tools import ensembl, ncbi_gene, uniprot


@dataclass(frozen=True)
class SpeciesIdentitySpec:
    """One target species for identity resolution."""

    common_name: str
    taxon_id: int
    scientific_name: str
    ncbi_organism: str
    ensembl_species: str


SPECIES_IDENTITY_SPECS: tuple[SpeciesIdentitySpec, ...] = (
    SpeciesIdentitySpec(
        common_name="human",
        taxon_id=9606,
        scientific_name="Homo sapiens",
        ncbi_organism=ncbi_gene.ORGANISM_HUMAN,
        ensembl_species=ensembl.SPECIES_HUMAN,
    ),
    SpeciesIdentitySpec(
        common_name="mouse",
        taxon_id=10090,
        scientific_name="Mus musculus",
        ncbi_organism=ncbi_gene.ORGANISM_MOUSE,
        ensembl_species=ensembl.SPECIES_MOUSE,
    ),
    SpeciesIdentitySpec(
        common_name="rat",
        taxon_id=10116,
        scientific_name="Rattus norvegicus",
        ncbi_organism=ncbi_gene.ORGANISM_RAT,
        ensembl_species=ensembl.SPECIES_RAT,
    ),
)


def _tag_tool_result(
    result: ToolResult,
    *,
    endpoint_name: str,
    query_symbol: str,
    resolved_symbol: str,
    extra_params: dict[str, Any] | None = None,
) -> ToolResult:
    """Tag a ToolResult; always keep ``gene_symbol`` as the dossier query."""
    params = dict(result.request_params or {})
    params["species_identity"] = True
    params["query_symbol"] = query_symbol
    params["resolved_symbol"] = resolved_symbol
    if extra_params:
        params.update(extra_params)
    return result.model_copy(
        update={
            "endpoint_name": endpoint_name,
            "gene_symbol": query_symbol,
            "request_params": params,
        }
    )


def resolved_symbol_from_ncbi(ncbi_result: ToolResult, fallback: str) -> str:
    """Prefer NCBI official nomenclature symbol when a safe selection exists."""
    if not ncbi_result.success or not isinstance(ncbi_result.data, dict):
        return fallback
    summary = ncbi_result.data.get("selected_summary")
    if isinstance(summary, dict):
        official = summary.get("nomenclaturesymbol") or summary.get("name")
        if official and str(official).strip():
            return str(official).strip()
    return fallback


def covered_ncbi_taxons(tool_results: Iterable[ToolResult]) -> set[int]:
    """Return taxon IDs already covered by successful NCBI Gene identity results."""
    covered: set[int] = set()
    for tr in tool_results:
        if tr.source_name != "NCBI Gene" or not tr.success:
            continue
        data = tr.data if isinstance(tr.data, dict) else {}
        tax = data.get("expected_taxid")
        if tax is None and isinstance(data.get("selected_summary"), dict):
            organism = data["selected_summary"].get("organism") or {}
            if isinstance(organism, dict):
                tax = organism.get("taxid")
        try:
            if tax is not None:
                covered.add(int(tax))
        except (TypeError, ValueError):
            continue
    return covered


def fetch_species_identity_results(
    gene_symbol: str,
    *,
    settings: Settings | None = None,
    skip_taxons: Iterable[int] | None = None,
    specs: Iterable[SpeciesIdentitySpec] | None = None,
) -> list[ToolResult]:
    """Fetch NCBI Gene + Ensembl + UniProt identity ToolResults per species.

    For each species, NCBI runs first so Ensembl/UniProt use the resolved
    species-specific symbol (not an assumed human-case symbol).
    ``ToolResult.gene_symbol`` remains the dossier query symbol throughout.
    """
    cfg = settings or get_settings()
    query_symbol = gene_symbol.strip()
    skip = {int(t) for t in (skip_taxons or [])}
    results: list[ToolResult] = []
    for spec in specs or SPECIES_IDENTITY_SPECS:
        if spec.taxon_id in skip:
            continue

        ncbi = ncbi_gene.lookup_gene(
            query_symbol, organism=spec.ncbi_organism, settings=cfg
        )
        resolved = resolved_symbol_from_ncbi(ncbi, query_symbol)
        ncbi = _tag_tool_result(
            ncbi,
            endpoint_name=f"lookup_gene_{spec.common_name}",
            query_symbol=query_symbol,
            resolved_symbol=resolved,
            extra_params={
                "common_name": spec.common_name,
                "taxon_id": spec.taxon_id,
                "scientific_name": spec.scientific_name,
            },
        )
        if isinstance(ncbi.data, dict):
            data = dict(ncbi.data)
            data.setdefault("organism", spec.ncbi_organism)
            data.setdefault("expected_taxid", spec.taxon_id)
            data.setdefault("scientific_name", spec.scientific_name)
            data.setdefault("common_name", spec.common_name)
            data["query_gene_symbol"] = query_symbol
            data["species_gene_symbol"] = resolved
            data["query_symbol"] = query_symbol
            data["resolved_symbol"] = resolved
            ncbi = ncbi.model_copy(update={"data": data})
        results.append(ncbi)

        ens = ensembl.lookup_symbol(
            resolved, species=spec.ensembl_species, settings=cfg
        )
        ens = _tag_tool_result(
            ens,
            endpoint_name=f"lookup_symbol_{spec.common_name}",
            query_symbol=query_symbol,
            resolved_symbol=resolved,
            extra_params={
                "common_name": spec.common_name,
                "taxon_id": spec.taxon_id,
                "scientific_name": spec.scientific_name,
            },
        )
        if isinstance(ens.data, dict):
            data = dict(ens.data)
            data.setdefault("species", spec.ensembl_species)
            data.setdefault("taxon_id", spec.taxon_id)
            data.setdefault("scientific_name", spec.scientific_name)
            data["query_gene_symbol"] = query_symbol
            data["species_gene_symbol"] = resolved
            data["query_symbol"] = query_symbol
            data["resolved_symbol"] = resolved
            ens = ens.model_copy(update={"data": data})
        results.append(ens)

        up = uniprot.search_reviewed(
            resolved, organism_id=spec.taxon_id, settings=cfg
        )
        up = _tag_tool_result(
            up,
            endpoint_name=f"search_reviewed_{spec.common_name}",
            query_symbol=query_symbol,
            resolved_symbol=resolved,
            extra_params={
                "common_name": spec.common_name,
                "taxon_id": spec.taxon_id,
                "scientific_name": spec.scientific_name,
            },
        )
        if isinstance(up.data, dict):
            data = dict(up.data)
            data.setdefault("organism_id", spec.taxon_id)
            data.setdefault("scientific_name", spec.scientific_name)
            data.setdefault("common_name", spec.common_name)
            data["query_gene_symbol"] = query_symbol
            data["species_gene_symbol"] = resolved
            data["query_symbol"] = query_symbol
            data["resolved_symbol"] = resolved
            up = up.model_copy(update={"data": data})
        # When no Swiss-Prot hit exists for this taxon (e.g. rat CDH10 → TrEMBL
        # F1LR98), fall back to an unreviewed UniProtKB search so Section 1a/1d
        # share the same structured identity accession evidence.
        selected = None
        if isinstance(up.data, dict):
            selected = up.data.get("selected_accession")
        if up.success and not selected:
            fallback = uniprot.search_gene_symbol(
                resolved,
                organism_id=spec.taxon_id,
                reviewed=False,
                settings=cfg,
                endpoint_name=f"search_unreviewed_{spec.common_name}",
            )
            fallback = _tag_tool_result(
                fallback,
                endpoint_name=f"search_unreviewed_{spec.common_name}",
                query_symbol=query_symbol,
                resolved_symbol=resolved,
                extra_params={
                    "common_name": spec.common_name,
                    "taxon_id": spec.taxon_id,
                    "scientific_name": spec.scientific_name,
                    "reviewed_fallback": True,
                },
            )
            if isinstance(fallback.data, dict):
                data = dict(fallback.data)
                data.setdefault("organism_id", spec.taxon_id)
                data.setdefault("scientific_name", spec.scientific_name)
                data.setdefault("common_name", spec.common_name)
                data["query_gene_symbol"] = query_symbol
                data["species_gene_symbol"] = resolved
                data["query_symbol"] = query_symbol
                data["resolved_symbol"] = resolved
                data["reviewed_fallback"] = True
                fallback = fallback.model_copy(update={"data": data})
            if fallback.success and isinstance(fallback.data, dict) and fallback.data.get(
                "selected_accession"
            ):
                up = fallback
        results.append(up)

    return results


__all__ = [
    "SpeciesIdentitySpec",
    "SPECIES_IDENTITY_SPECS",
    "resolved_symbol_from_ncbi",
    "covered_ncbi_taxons",
    "fetch_species_identity_results",
]
