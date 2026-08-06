"""LangGraph orchestration for a full gene-dossier API pass.

Pipeline (each node is a small, testable function)::

    create_dossier_run
      -> resolve_gene_identity
      -> call_source_clients
      -> save_raw_artifacts
      -> normalize_evidence
      -> (optional) index_evidence_in_chroma   # soft-fail via retrieval.py
      -> build_report_sections
      -> verify_claims
      -> render_outputs

Rules:
- Failures in one source never abort the graph.
- LLM keys are optional; synthesis falls back to deterministic markdown.
- Final polished output is the Rancho/CHDI report; markdown remains a
  provenance/debug view.
- No invented facts: evidence comes only from ToolResult -> normalizers.
"""

from __future__ import annotations

import contextvars
import logging
import threading
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, TypedDict
from uuid import uuid4

from gene_dossier.config import Settings, get_settings
from gene_dossier.coverage import build_and_write_coverage, persist_coverage
from gene_dossier.db import (
    init_db,
    save_api_run,
    save_dossier_run,
    save_evidence_record,
    save_raw_artifact,
    session_scope,
)
from gene_dossier.models import (
    ApiRun,
    Claim,
    DossierRun,
    EvidenceRecord,
    RawArtifact,
    ReportSection,
    SourceCoverageResult,
    SourceStatus,
    ToolResult,
    VerificationResult,
    utcnow,
)
from gene_dossier.raw_store import RawStore
from gene_dossier.rancho_report import build_and_write_rancho_report
from gene_dossier.rendering import write_dossier_report
from gene_dossier.source_registry import get_all_sources, get_source
from gene_dossier.synthesis import SynthesisResult, synthesize_dossier
from gene_dossier.ucsc_figure import install_ucsc_api_key_log_redaction
from gene_dossier.verification import verify_claims

logger = logging.getLogger(__name__)
install_ucsc_api_key_log_redaction()

ClientFn = Callable[..., ToolResult]
NormalizerFn = Callable[..., list[EvidenceRecord]]

IDENTITY_SOURCES = ("NCBI Gene", "Ensembl", "UniProt")
_UCSC_TRANSIENT_TOKEN_KEY = "_transient_ucsc_figure_token"


@dataclass
class WorkflowTransientContext:
    """Per-compiled-graph transient payloads; never enters LangGraph state.

    Access is guarded by ``lock`` so concurrent invokes of the same compiled
    graph cannot overwrite or consume one another's live figure / binary
    artifact payloads.
    """

    live_figures: dict[str, Any] = field(default_factory=dict)
    request_cache: dict[str, Any] = field(default_factory=dict)
    lock: threading.Lock = field(default_factory=threading.Lock)

    def put_figure(self, token: str, payload: Any) -> None:
        with self.lock:
            self.live_figures[token] = payload

    def pop_figure(self, token: str) -> Any | None:
        with self.lock:
            return self.live_figures.pop(token, None)

    def put(self, dossier_run_id: str, token: str, payload: Any) -> str:
        """Store an arbitrary transient payload and return its run-scoped key."""
        key = token if token.startswith(f"{dossier_run_id}:") else f"{dossier_run_id}:{token}"
        with self.lock:
            self.live_figures[key] = payload
        return key

    def pop(self, dossier_run_id: str, token: str) -> Any | None:
        """Consume an arbitrary transient payload by run-scoped or local token."""
        key = token if token.startswith(f"{dossier_run_id}:") else f"{dossier_run_id}:{token}"
        with self.lock:
            return self.live_figures.pop(key, None)

    def get_cached_request(self, identity_key: str) -> Any | None:
        with self.lock:
            return self.request_cache.get(identity_key)

    def put_cached_request(self, identity_key: str, payload: Any) -> None:
        with self.lock:
            self.request_cache[identity_key] = payload

    def clear_run(self, dossier_run_id: str) -> None:
        """Drop any remaining payloads keyed for ``dossier_run_id``."""
        prefix = f"{dossier_run_id}:"
        with self.lock:
            stale = [k for k in self.live_figures if k.startswith(prefix)]
            for key in stale:
                self.live_figures.pop(key, None)


_ACTIVE_TRANSIENT: contextvars.ContextVar[WorkflowTransientContext | None] = (
    contextvars.ContextVar("gene_dossier_active_transient", default=None)
)


def bind_workflow_transient(
    transient: WorkflowTransientContext | None,
) -> contextvars.Token[WorkflowTransientContext | None]:
    """Bind ``transient`` for nested client calls in this context."""
    return _ACTIVE_TRANSIENT.set(transient)


def reset_workflow_transient(
    token: contextvars.Token[WorkflowTransientContext | None],
) -> None:
    _ACTIVE_TRANSIENT.reset(token)


def get_active_workflow_transient() -> WorkflowTransientContext | None:
    return _ACTIVE_TRANSIENT.get()


class DossierState(TypedDict, total=False):
    """Mutable LangGraph state for one dossier pass."""

    gene_symbol: str
    dossier_run_id: str
    official_symbol: str | None
    gene_ids: dict[str, Any]
    tool_results: list[ToolResult]
    api_runs: list[ApiRun]
    raw_artifacts: list[dict[str, Any]]
    evidence_records: list[EvidenceRecord]
    coverage: list[SourceCoverageResult]
    sections: list[ReportSection]
    claims: list[Claim]
    verification_results: list[VerificationResult]
    synthesis_mode: str | None
    synthesis_notes: list[str]
    output_paths: dict[str, str]
    errors: list[str]
    status: str


@dataclass
class DossierPassResult:
    """Convenience return object for :func:`run_gene_dossier_full_api_pass`."""

    gene_symbol: str
    dossier_run_id: str
    status: str
    gene_ids: dict[str, Any] = field(default_factory=dict)
    evidence_records: list[EvidenceRecord] = field(default_factory=list)
    coverage: list[SourceCoverageResult] = field(default_factory=list)
    sections: list[ReportSection] = field(default_factory=list)
    claims: list[Claim] = field(default_factory=list)
    verification_results: list[VerificationResult] = field(default_factory=list)
    output_paths: dict[str, Path] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    synthesis_mode: str | None = None
    synthesis_notes: list[str] = field(default_factory=list)

    @classmethod
    def from_state(cls, state: DossierState) -> DossierPassResult:
        paths = {
            key: Path(value)
            for key, value in (state.get("output_paths") or {}).items()
            if value
        }
        return cls(
            gene_symbol=state.get("gene_symbol") or "",
            dossier_run_id=state.get("dossier_run_id") or "",
            status=state.get("status") or "unknown",
            gene_ids=dict(state.get("gene_ids") or {}),
            evidence_records=list(state.get("evidence_records") or []),
            coverage=list(state.get("coverage") or []),
            sections=list(state.get("sections") or []),
            claims=list(state.get("claims") or []),
            verification_results=list(state.get("verification_results") or []),
            output_paths=paths,
            errors=list(state.get("errors") or []),
            synthesis_mode=state.get("synthesis_mode"),
            synthesis_notes=list(state.get("synthesis_notes") or []),
        )


def _safe_call_client(
    name: str,
    fn: ClientFn,
    *,
    gene_symbol: str,
    gene_ids: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """Invoke one client; never raise — convert unexpected errors into ToolResult."""
    try:
        return fn(gene_symbol=gene_symbol, gene_ids=gene_ids, settings=settings)
    except Exception as exc:  # noqa: BLE001
        logger.exception("Client %s raised unexpectedly", name)
        return ToolResult(
            source_name=name,
            endpoint_name="workflow_dispatch",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def _client_ncbi_gene(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import ncbi_gene

    return ncbi_gene.lookup_gene(gene_symbol, settings=settings)


def _client_ensembl(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import ensembl

    return ensembl.lookup_symbol(gene_symbol, settings=settings)


def _client_uniprot(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import uniprot

    return uniprot.search_reviewed(gene_symbol, settings=settings)


def _client_pubmed(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import pubmed

    return pubmed.search_hd_literature(gene_symbol, settings=settings)


def _client_gtex(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import gtex

    kwargs: dict[str, Any] = {"settings": settings}
    if gene_ids.get("gtex_gencode_id"):
        kwargs["gencode_id"] = gene_ids["gtex_gencode_id"]
    return gtex.fetch_expression_and_eqtl(gene_symbol, **kwargs)


def _client_string(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import string_db

    return string_db.fetch_interaction_partners(gene_symbol, settings=settings)


def _client_reactome(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import reactome

    accession = gene_ids.get("uniprot_accession")
    if not accession:
        return ToolResult(
            source_name="Reactome",
            endpoint_name="fetch_pathways",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="missing_identifier",
            error_message="Reactome requires uniprot_accession from identity resolution",
        )
    return reactome.fetch_pathways(
        str(accession), gene_symbol=gene_symbol, settings=settings
    )


def _client_clinvar(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import clinvar

    return clinvar.fetch_clinvar_variants(gene_symbol, settings=settings)


def _client_opentargets(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import opentargets

    ensembl_id = gene_ids.get("ensembl_id")
    if not ensembl_id:
        return ToolResult(
            source_name="Open Targets",
            endpoint_name="fetch_disease_associations",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="missing_identifier",
            error_message="Open Targets requires ensembl_id from identity resolution",
        )
    return opentargets.fetch_disease_associations(
        str(ensembl_id), gene_symbol=gene_symbol, settings=settings
    )


def _client_mousemine(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import mousemine

    mouse_id = gene_ids.get("mouse_entrez_id")
    return mousemine.fetch_mouse_annotations(
        gene_symbol=gene_symbol,
        ncbi_gene_number=int(mouse_id) if mouse_id else None,
        mgi_id=gene_ids.get("mgi_id"),
        settings=settings,
    )


def _client_ctd(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import ctd

    return ctd.fetch_chemical_gene_interactions(gene_symbol, settings=settings)


def _client_chembl(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import chembl

    return chembl.fetch_chemical_tools(gene_symbol, settings=settings)


def _client_pubchem(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import pubchem

    gene_id = gene_ids.get("entrez_gene_id")
    if not gene_id:
        return ToolResult(
            source_name="PubChem",
            endpoint_name="fetch_bioassays",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="missing_identifier",
            error_message="PubChem requires entrez_gene_id from identity resolution",
        )
    return pubchem.fetch_bioassays(
        str(gene_id), gene_symbol=gene_symbol, settings=settings
    )


def _client_nih_reporter(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import nih_reporter

    return nih_reporter.fetch_grants(gene_symbol, settings=settings)


def _client_geo(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import geo

    return geo.fetch_perturbations(gene_symbol, settings=settings)


def _client_harmonizome(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import harmonizome

    return harmonizome.fetch_tf_associations(
        gene_symbol,
        settings=settings,
        transient=get_active_workflow_transient(),
    )


def _client_biogrid(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import biogrid

    return biogrid.fetch_interactions(gene_symbol, settings=settings)


def _client_wikipathways(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import wikipathways

    entrez = gene_ids.get("entrez_gene_id")
    ensembl = gene_ids.get("ensembl_id")
    uniprot = gene_ids.get("uniprot_accession")
    return wikipathways.fetch_pathways(
        gene_symbol,
        entrez_ids=[entrez] if entrez else None,
        ensembl_ids=[ensembl] if ensembl else None,
        uniprot_ids=[uniprot] if uniprot else None,
        settings=settings,
    )


def _client_alphafold(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import alphafold

    accession = gene_ids.get("uniprot_accession")
    if not accession:
        return ToolResult(
            source_name="AlphaFold",
            endpoint_name="fetch_prediction",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="missing_identifier",
            error_message="AlphaFold requires uniprot_accession from identity resolution",
        )
    return alphafold.fetch_prediction(
        accession, gene_symbol=gene_symbol, settings=settings
    )


def _client_pdbe(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import pdbe

    accession = gene_ids.get("uniprot_accession")
    if not accession:
        return ToolResult(
            source_name="PDBe",
            endpoint_name="fetch_structures",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="missing_identifier",
            error_message="PDBe requires uniprot_accession from identity resolution",
        )
    return pdbe.fetch_structures(
        accession, gene_symbol=gene_symbol, settings=settings
    )


def _client_cdd(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import cdd

    query = gene_ids.get("refseq_protein") or gene_ids.get("uniprot_accession")
    if not query:
        return ToolResult(
            source_name="CDD",
            endpoint_name="fetch_domains",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="missing_identifier",
            error_message="CDD requires refseq_protein or uniprot_accession",
        )
    return cdd.fetch_domains(str(query), gene_symbol=gene_symbol, settings=settings)


def _client_ncbi_datasets(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import ncbi_datasets

    gene_id = gene_ids.get("entrez_gene_id")
    if not gene_id:
        return ToolResult(
            source_name="NCBI Datasets",
            endpoint_name="fetch_orthologs",
            success=False,
            gene_symbol=gene_symbol,
            request_url="",
            error_type="missing_identifier",
            error_message="NCBI Datasets requires entrez_gene_id",
        )
    return ncbi_datasets.fetch_orthologs(
        str(gene_id), gene_symbol=gene_symbol, settings=settings
    )


def _client_ucsc(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import ucsc

    return ucsc.fetch_gene_region(gene_symbol, settings=settings)


def _client_allen(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import allen_brain

    return allen_brain.fetch_hba_expression(gene_symbol, settings=settings)


def _client_brainrnaseq(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import brainrnaseq

    return brainrnaseq.fetch_gene_expression(gene_symbol, settings=settings)


def _client_patents(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import patents

    return patents.fetch_patents(gene_symbol, settings=settings)


def _client_antibodies(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import antibodies

    return antibodies.fetch_antibodies(gene_symbol, settings=settings)


def _client_omim(
    *, gene_symbol: str, gene_ids: dict[str, Any], settings: Settings
) -> ToolResult:
    from gene_dossier.tools import omim

    return omim.fetch_gene_entry(gene_symbol, settings=settings)


CLIENT_DISPATCH: dict[str, ClientFn] = {
    "NCBI Gene": _client_ncbi_gene,
    "Ensembl": _client_ensembl,
    "UniProt": _client_uniprot,
    "PubMed": _client_pubmed,
    "GTEx": _client_gtex,
    "STRING": _client_string,
    "Reactome": _client_reactome,
    "ClinVar": _client_clinvar,
    "Open Targets": _client_opentargets,
    "MouseMine": _client_mousemine,
    "CTD": _client_ctd,
    "ChEMBL": _client_chembl,
    "PubChem": _client_pubchem,
    "NIH RePORTER": _client_nih_reporter,
    "GEO": _client_geo,
    "Harmonizome": _client_harmonizome,
    "BioGRID": _client_biogrid,
    "WikiPathways": _client_wikipathways,
    "AlphaFold": _client_alphafold,
    "PDBe": _client_pdbe,
    "CDD": _client_cdd,
    "NCBI Datasets": _client_ncbi_datasets,
    "UCSC": _client_ucsc,
    "Allen Brain Atlas": _client_allen,
    "BrainRNASeq": _client_brainrnaseq,
    "Patents": _client_patents,
    "Antibodies": _client_antibodies,
    "OMIM": _client_omim,
}


def normalize_tool_result(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Run all applicable normalizers for one ToolResult; dedupe by source_id."""
    from gene_dossier.normalize.chemicals import normalize_chemicals
    from gene_dossier.normalize.expression import normalize_expression
    from gene_dossier.normalize.gene_identity import normalize_gene_identity
    from gene_dossier.normalize.grants import normalize_grants
    from gene_dossier.normalize.literature import normalize_literature
    from gene_dossier.normalize.model_organisms import normalize_model_organisms
    from gene_dossier.normalize.pathways import normalize_pathways
    from gene_dossier.normalize.perturbation import normalize_perturbation
    from gene_dossier.normalize.ppi import normalize_ppi
    from gene_dossier.normalize.protein import normalize_protein
    from gene_dossier.normalize.transcription_factors import (
        normalize_transcription_factors,
    )
    from gene_dossier.normalize.variants import normalize_variants

    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    dispatchers: list[NormalizerFn] = [
        normalize_gene_identity,
        normalize_protein,
        normalize_expression,
        normalize_literature,
        normalize_ppi,
        normalize_pathways,
        normalize_chemicals,
        normalize_variants,
        normalize_model_organisms,
        normalize_perturbation,
        normalize_transcription_factors,
        normalize_grants,
    ]
    records: list[EvidenceRecord] = []
    seen: set[str] = set()
    for fn in dispatchers:
        try:
            batch = fn(tool_result, **kwargs)
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "Normalizer %s failed for %s: %s",
                getattr(fn, "__name__", fn),
                tool_result.source_name,
                exc,
            )
            continue
        for rec in batch:
            if rec.source_id and rec.source_id in seen:
                continue
            if rec.source_id:
                seen.add(rec.source_id)
            records.append(rec)
    return records


def extract_gene_ids_from_tool_result(
    tool_result: ToolResult, gene_ids: dict[str, Any]
) -> dict[str, Any]:
    """Merge identity anchors discovered in a successful ToolResult into ``gene_ids``."""
    if not tool_result.success or not isinstance(tool_result.data, dict):
        return gene_ids
    data = tool_result.data
    updated = dict(gene_ids)
    source = tool_result.source_name

    if source == "NCBI Gene":
        gid = data.get("selected_gene_id")
        tax = data.get("expected_taxid")
        summary = data.get("selected_summary") or {}
        if not isinstance(summary, dict):
            summary = {}
        if tax is None:
            organism = summary.get("organism") or {}
            if isinstance(organism, dict):
                tax = organism.get("taxid")
        try:
            tax_int = int(tax) if tax is not None else None
        except (TypeError, ValueError):
            tax_int = None
        # Legacy human-only payloads omit taxid; treat as human.
        if tax_int is None:
            tax_int = 9606
        if gid:
            if tax_int == 9606:
                updated["entrez_gene_id"] = str(gid)
            elif tax_int == 10090:
                updated["mouse_entrez_id"] = str(gid)
            elif tax_int == 10116:
                updated["rat_entrez_id"] = str(gid)
        official = summary.get("nomenclaturesymbol") or summary.get("name")
        if official:
            if tax_int == 9606:
                updated["official_symbol"] = str(official)
            elif tax_int == 10090:
                updated["mouse_symbol"] = str(official)
            elif tax_int == 10116:
                updated["rat_symbol"] = str(official)
        chrom = summary.get("chromosome")
        if chrom and tax_int == 9606:
            updated["chromosome"] = str(chrom)
    elif source == "Ensembl":
        eid = data.get("ensembl_id") or data.get("id")
        if not eid and isinstance(data.get("summary"), dict):
            eid = data["summary"].get("ensembl_gene_id") or data["summary"].get("id")
        if not eid and isinstance(data.get("gene"), dict):
            eid = data["gene"].get("id")
        species = data.get("species")
        if not species and isinstance(data.get("summary"), dict):
            species = data["summary"].get("species")
        if eid:
            if species in (None, "", "homo_sapiens"):
                updated["ensembl_id"] = str(eid)
            elif species == "mus_musculus":
                updated["mouse_ensembl_id"] = str(eid)
            elif species == "rattus_norvegicus":
                updated["rat_ensembl_id"] = str(eid)
    elif source == "UniProt":
        acc = data.get("selected_accession") or data.get("primaryAccession")
        selected_entry: dict[str, Any] = {}
        entries = data.get("entries") or data.get("selected") or data.get("results") or []
        if not acc:
            if isinstance(entries, list) and entries and isinstance(entries[0], dict):
                selected_entry = entries[0]
                acc = selected_entry.get("primaryAccession") or selected_entry.get("accession")
        elif isinstance(entries, list):
            for entry in entries:
                if not isinstance(entry, dict):
                    continue
                entry_acc = entry.get("primaryAccession") or entry.get("accession")
                if str(entry_acc) == str(acc):
                    selected_entry = entry
                    break
            if not selected_entry and entries and isinstance(entries[0], dict):
                selected_entry = entries[0]
        organism_id = data.get("organism_id")
        if organism_id is None and selected_entry:
            organism_id = selected_entry.get("organism_id")
        try:
            tax_int = int(organism_id) if organism_id is not None else 9606
        except (TypeError, ValueError):
            tax_int = 9606
        if acc:
            if tax_int == 9606:
                updated["uniprot_accession"] = str(acc)
                protein_length = (
                    selected_entry.get("protein_length")
                    or selected_entry.get("sequence_length")
                    or selected_entry.get("length")
                )
                try:
                    if protein_length is not None:
                        updated["protein_length"] = int(protein_length)
                        updated["protein_length_source"] = "reviewed canonical UniProt sequence length"
                except (TypeError, ValueError):
                    pass
                refseq_accessions = selected_entry.get("refseq_protein_accessions") or []
                if isinstance(refseq_accessions, list):
                    selected_refseq = next(
                        (
                            str(x).strip()
                            for x in refseq_accessions
                            if str(x).strip().startswith(("NP_", "XP_", "YP_"))
                        ),
                        None,
                    )
                    if selected_refseq:
                        updated["refseq_protein"] = selected_refseq
                        updated["refseq_protein_accessions"] = [
                            str(x).strip()
                            for x in refseq_accessions
                            if str(x).strip()
                        ]
            elif tax_int == 10090:
                updated["mouse_uniprot_accession"] = str(acc)
            elif tax_int == 10116:
                updated["rat_uniprot_accession"] = str(acc)
    elif source == "NCBI Datasets":
        for row in data.get("ortholog_summaries") or []:
            if not isinstance(row, dict):
                continue
            tax = str(row.get("tax_id") or row.get("taxon_id") or "")
            organism = str(row.get("organism") or "").lower()
            if tax in {"10090", "Mus musculus"} or "mouse" in organism:
                mid = row.get("gene_id") or row.get("ncbi_gene_id")
                if mid:
                    updated["mouse_entrez_id"] = str(mid)
                mgi = row.get("mgi_id") or row.get("mgi")
                if mgi:
                    updated["mgi_id"] = str(mgi)
                break
    elif source == "MouseMine":
        mgi = data.get("mgi_id")
        if mgi:
            updated["mgi_id"] = str(mgi)
    elif source == "GTEx":
        gencode = data.get("gencode_id")
        if gencode:
            updated["gtex_gencode_id"] = str(gencode)
            if "ensembl_id" not in updated and "." in str(gencode):
                updated["ensembl_id"] = str(gencode).split(".", 1)[0]
            elif "ensembl_id" not in updated and str(gencode).startswith("ENSG"):
                updated["ensembl_id"] = str(gencode)
    elif source == "UCSC":
        region = data.get("region") or data.get("selected_region")
        if region:
            updated["ucsc_region"] = str(region)

    return updated


def node_create_dossier_run(
    state: DossierState, *, settings: Settings, persist_db: bool = True
) -> DossierState:
    """Create the DossierRun and optionally persist it."""
    gene_symbol = (state.get("gene_symbol") or "").strip()
    if not gene_symbol:
        return {
            **state,
            "status": "failed",
            "errors": list(state.get("errors") or []) + ["gene_symbol is required"],
        }

    run = DossierRun(
        gene_symbol=gene_symbol,
        run_type="full_api_pass",
        status="running",
    )
    if state.get("dossier_run_id"):
        run.id = state["dossier_run_id"]

    if persist_db:
        init_db()
        with session_scope() as session:
            save_dossier_run(session, run)

    return {
        **state,
        "dossier_run_id": run.id,
        "gene_symbol": gene_symbol,
        "gene_ids": dict(state.get("gene_ids") or {}),
        "tool_results": list(state.get("tool_results") or []),
        "api_runs": [],
        "raw_artifacts": [],
        "evidence_records": [],
        "coverage": [],
        "sections": [],
        "claims": [],
        "verification_results": [],
        "output_paths": {},
        "errors": list(state.get("errors") or []),
        "status": "running",
    }


def node_resolve_gene_identity(
    state: DossierState, *, settings: Settings, call_network: bool = True
) -> DossierState:
    """Call identity sources first and collect chained identifiers.

    Fetches NCBI Gene / Ensembl / UniProt for Human, Mouse, and Rat when
    networking is enabled. Species-specific symbols are resolved via NCBI
    before Ensembl/UniProt requests.
    """
    from gene_dossier.species_identity import (
        covered_ncbi_taxons,
        fetch_species_identity_results,
    )

    gene_symbol = state["gene_symbol"]
    gene_ids = dict(state.get("gene_ids") or {})
    tool_results = list(state.get("tool_results") or [])
    errors = list(state.get("errors") or [])

    # Always harvest IDs from preloaded / already-fetched identity ToolResults
    # before deciding whether to call the network.
    for existing in tool_results:
        gene_ids = extract_gene_ids_from_tool_result(existing, gene_ids)

    if not call_network:
        official = gene_ids.get("official_symbol")
        return {
            **state,
            "gene_ids": gene_ids,
            "tool_results": tool_results,
            "official_symbol": official or state.get("official_symbol"),
            "errors": errors,
        }

    skip_taxons = covered_ncbi_taxons(tool_results)
    # If only legacy human-only identity is present (no species_identity tag),
    # still fetch mouse/rat; skip human to avoid duplicate human calls.
    has_tagged_species = any(
        isinstance(tr.request_params, dict) and tr.request_params.get("species_identity")
        for tr in tool_results
    )
    if not has_tagged_species and any(
        tr.source_name == "NCBI Gene" and tr.success for tr in tool_results
    ):
        skip_taxons = set(skip_taxons) | {9606}

    try:
        species_results = fetch_species_identity_results(
            gene_symbol, settings=settings, skip_taxons=skip_taxons
        )
    except Exception as exc:  # noqa: BLE001 — identity must never abort the graph
        errors.append(f"species_identity: {type(exc).__name__}: {exc}")
        species_results = []

    for result in species_results:
        tool_results.append(result)
        gene_ids = extract_gene_ids_from_tool_result(result, gene_ids)
        if not result.success:
            errors.append(
                f"{result.source_name}/{result.endpoint_name}: "
                f"{result.error_message or result.error_type}"
            )
        elif (
            result.source_name == "NCBI Gene"
            and isinstance(result.data, dict)
            and result.data.get("selection_method") == "ambiguous"
        ):
            warnings = result.data.get("selection_warnings") or []
            errors.append(
                f"NCBI Gene/{result.endpoint_name}: ambiguous selection "
                f"({', '.join(str(w) for w in warnings) or 'no unique safe match'})"
            )

    official = gene_ids.get("official_symbol")
    return {
        **state,
        "gene_ids": gene_ids,
        "tool_results": tool_results,
        "official_symbol": official or state.get("official_symbol"),
        "errors": errors,
    }


def node_call_source_clients(
    state: DossierState,
    *,
    settings: Settings,
    call_network: bool = True,
    sources: Iterable[str] | None = None,
    transient: WorkflowTransientContext | None = None,
) -> DossierState:
    """Call remaining registered clients (skip identity sources already fetched)."""
    if not call_network:
        return state

    gene_symbol = state["gene_symbol"]
    gene_ids = dict(state.get("gene_ids") or {})
    tool_results = list(state.get("tool_results") or [])
    errors = list(state.get("errors") or [])
    already = {tr.source_name for tr in tool_results}

    wanted: set[str] | None = None
    if sources is not None:
        wanted = {s.strip().lower() for s in sources}

    for src in get_all_sources():
        name = src.name
        if name in IDENTITY_SOURCES:
            continue
        if wanted is not None and name.lower() not in wanted:
            continue
        if name in already:
            continue
        fn = CLIENT_DISPATCH.get(name)
        if fn is None:
            continue
        missing = [k for k in src.required_keys if not settings.has_key(k)]
        if missing:
            tool_results.append(
                ToolResult(
                    source_name=name,
                    endpoint_name="workflow_skip",
                    success=False,
                    gene_symbol=gene_symbol,
                    request_url="",
                    error_type="requires_key",
                    error_message=f"missing required key(s): {', '.join(missing)}",
                )
            )
            continue

        bind_token = bind_workflow_transient(transient)
        try:
            if name == "UCSC":
                from gene_dossier.tools import ucsc as ucsc_client

                execution = ucsc_client.fetch_gene_region_execution(
                    gene_symbol,
                    settings=settings,
                )
                result = execution.tool_result
                if transient is not None and execution.live_figure is not None:
                    data = dict(result.data) if isinstance(result.data, dict) else {}
                    token = f"{state['dossier_run_id']}:{uuid4().hex}"
                    data[_UCSC_TRANSIENT_TOKEN_KEY] = token
                    result.data = data
                    transient.put_figure(token, execution.live_figure)
            else:
                result = _safe_call_client(
                    name, fn, gene_symbol=gene_symbol, gene_ids=gene_ids, settings=settings
                )
        finally:
            reset_workflow_transient(bind_token)
        tool_results.append(result)
        gene_ids = extract_gene_ids_from_tool_result(result, gene_ids)
        if not result.success:
            errors.append(f"{name}: {result.error_message or result.error_type}")

    return {
        **state,
        "gene_ids": gene_ids,
        "tool_results": tool_results,
        "errors": errors,
    }


def _persist_ucsc_result_with_live_figure(
    *,
    result: ToolResult,
    dossier_run_id: str,
    gene_symbol: str,
    settings: Settings,
    store: RawStore,
    transient: WorkflowTransientContext | None,
    persist_db: bool,
) -> tuple[ToolResult, list[ApiRun], list[dict[str, Any]]]:
    """Persist combined UCSC JSON plus any live figure payload.

    Persistence contract (two-phase):
    - Phase 1 (this function): ApiRun + RawArtifact rows and the final image
      file are written together. On failure, newly created files for this
      attempt are cleaned up.
    - Phase 2 (``node_normalize_evidence``): ``ucsc_conservation_figure``
      EvidenceRecords are created later. A normalization failure may leave an
      orphan live figure artifact without a matching EvidenceRecord; this is
      accepted and consistent with the existing pipeline architecture.

    Live figure bytes are consumed from workflow-local transient context before
    control returns to ordinary serializable state.
    """
    data = result.data if isinstance(result.data, dict) else {}
    token = str(data.get(_UCSC_TRANSIENT_TOKEN_KEY) or "").strip()
    live_figure = transient.pop_figure(token) if transient and token else None
    if token and isinstance(data, dict):
        data.pop(_UCSC_TRANSIENT_TOKEN_KEY, None)
        result.data = data

    api_runs: list[ApiRun] = []
    raw_meta: list[dict[str, Any]] = []

    combined_api_run = ApiRun(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        source_name=result.source_name,
        endpoint_name=result.endpoint_name,
        request_url=result.request_url,
        request_params=dict(result.request_params or {}),
        status_code=result.status_code,
        success=result.success,
        error_type=result.error_type,
        error_message=result.error_message,
    )
    json_artifact = None

    figure_api_runs: list[ApiRun] = []
    figure_artifact: RawArtifact | None = None
    staged = None
    created_final = False
    api_run_ids_by_index: dict[int, str] = {}
    json_cleanup_path: Path | None = None
    json_created_by_attempt = False

    try:
        if live_figure is not None:
            from gene_dossier.ucsc_figure import stage_figure_tempfile, sha256_hex

            if not live_figure.request_chain:
                raise RuntimeError(
                    "Validated live UCSC figure is missing HTTP request provenance"
                )

            ext = (
                live_figure.media_type.split("/")[-1]
                if "/" in live_figure.media_type
                else "png"
            )
            staged = stage_figure_tempfile(
                dossier_run_id=dossier_run_id,
                content=live_figure.content,
                extension=ext,
                settings=settings,
            )

            for idx, attempt in enumerate(live_figure.request_chain):
                parent_idx = attempt.get("parent_request_index")
                parent_id = (
                    api_run_ids_by_index.get(int(parent_idx))
                    if parent_idx is not None
                    else None
                )
                request_params = dict(attempt.get("request_params") or {})
                if parent_id is not None:
                    request_params["parent_api_run_id"] = parent_id
                api = ApiRun(
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    source_name="UCSC",
                    endpoint_name=str(
                        attempt.get("endpoint_name") or "hgRenderTracks"
                    ),
                    method="GET",
                    request_url=str(attempt.get("request_url") or ""),
                    request_params=request_params,
                    status_code=attempt.get("status_code"),
                    success=bool(attempt.get("success")),
                    error_type=attempt.get("error_type"),
                    error_message=attempt.get("error_message"),
                )
                figure_api_runs.append(api)
                api_run_ids_by_index[idx] = api.id

            image_api_run = figure_api_runs[live_figure.image_request_index]
            figure_artifact = RawArtifact(
                dossier_run_id=dossier_run_id,
                api_run_id=image_api_run.id,
                source_name="UCSC",
                artifact_type="image",
                file_path=staged.relative_path,
                original_url=None,
                content_hash=staged.sha256,
                notes="programmatic_browser_render",
            )

            figure_value = dict(data.get("figure") or {})
            figure_value.update(
                {
                    "relative_path": staged.relative_path,
                    "local_artifact_path": staged.relative_path,
                    "media_type": staged.media_type,
                    "width": staged.width,
                    "height": staged.height,
                    "byte_size": staged.byte_size,
                    "sha256": staged.sha256,
                    "genome": live_figure.genome,
                    "display_position": live_figure.display_position,
                    "selected_transcript": live_figure.selected_transcript,
                    "retrieval_method": "programmatic_browser_render",
                    "origin_endpoint": "hgRenderTracks",
                    "api_key_used": True,
                    "api_key_persisted": False,
                    "track_preset_id": live_figure.track_preset_id,
                    "track_preset_version": live_figure.track_preset_version,
                    "track_params": dict(live_figure.track_params),
                    "pixel_width": live_figure.track_params.get("pix")
                    or figure_value.get("pixel_width"),
                    "figure_api_run_id": image_api_run.id,
                    "figure_raw_artifact_id": figure_artifact.id,
                    "wrapper_request_index": live_figure.wrapper_request_index,
                    "image_request_index": live_figure.image_request_index,
                }
            )
            data["figure"] = figure_value
            result.data = data

        if result.data is not None:
            # Detect newly created JSON files so a later failure can clean them up
            # without deleting a pre-existing reused content-addressed artifact.
            from gene_dossier.raw_store import compute_hash
            from gene_dossier.source_ids import slugify
            import json as _json

            payload_bytes = _json.dumps(
                result.data, sort_keys=True, ensure_ascii=False, indent=2
            ).encode("utf-8")
            digest = compute_hash(payload_bytes)
            hint = slugify(result.endpoint_name) if result.endpoint_name else ""
            stem = f"{hint}-" if hint else ""
            expected_json = (
                store._dir_for(dossier_run_id, result.source_name)
                / f"{stem}{digest[:12]}.json"
            )
            json_created_by_attempt = not expected_json.exists()
            json_artifact = store.save_json(
                dossier_run_id,
                result.source_name,
                result.data,
                api_run_id=combined_api_run.id,
                original_url=result.request_url or None,
                filename_hint=result.endpoint_name,
            )
            json_cleanup_path = Path(json_artifact.file_path)
            combined_api_run.raw_artifact_id = json_artifact.id
            result.raw_artifact_id = json_artifact.id
            raw_meta.append(json_artifact.model_dump(mode="json"))

        # Finalize the staged image even when persist_db=False so HTML can
        # resolve and embed the figure path.
        if live_figure is not None and staged is not None and figure_artifact is not None:
            from gene_dossier.ucsc_figure import sha256_hex

            if staged.temp_path is not None:
                staged.temp_path.replace(staged.final_absolute_path)
                created_final = True
            if not staged.final_absolute_path.is_file():
                raise RuntimeError("Live UCSC figure final file missing after move")
            if sha256_hex(staged.final_absolute_path.read_bytes()) != staged.sha256:
                raise RuntimeError("Live UCSC figure checksum mismatch after move")
            raw_meta.append(figure_artifact.model_dump(mode="json"))

        if persist_db:
            with session_scope() as session:
                save_api_run(session, combined_api_run)
                if json_artifact is not None:
                    save_raw_artifact(session, json_artifact)
                if live_figure is not None and figure_artifact is not None:
                    for api in figure_api_runs:
                        save_api_run(session, api)
                    save_raw_artifact(session, figure_artifact)

        api_runs.append(combined_api_run)
        api_runs.extend(figure_api_runs)
        return result, api_runs, raw_meta
    except Exception:
        if transient is not None:
            transient.clear_run(dossier_run_id)
        if staged is not None and staged.temp_path is not None and staged.temp_path.exists():
            staged.temp_path.unlink(missing_ok=True)
        if (
            staged is not None
            and created_final
            and not staged.existed_already
            and staged.final_absolute_path.exists()
        ):
            staged.final_absolute_path.unlink(missing_ok=True)
        if (
            json_created_by_attempt
            and json_cleanup_path is not None
            and json_cleanup_path.exists()
        ):
            json_cleanup_path.unlink(missing_ok=True)
        raise


def node_save_raw_artifacts(
    state: DossierState,
    *,
    settings: Settings,
    persist_db: bool = True,
    transient: WorkflowTransientContext | None = None,
) -> DossierState:
    """Persist ToolResult payloads to the raw store + ApiRun / RawArtifact rows."""
    store = RawStore(base_dir=settings.raw_data_path)
    dossier_run_id = state["dossier_run_id"]
    gene_symbol = state["gene_symbol"]
    api_runs: list[ApiRun] = list(state.get("api_runs") or [])
    raw_meta: list[dict[str, Any]] = list(state.get("raw_artifacts") or [])
    updated_results: list[ToolResult] = []

    for result in state.get("tool_results") or []:
        if result.source_name == "UCSC":
            try:
                persisted_result, runs, meta = _persist_ucsc_result_with_live_figure(
                    result=result,
                    dossier_run_id=dossier_run_id,
                    gene_symbol=gene_symbol,
                    settings=settings,
                    store=store,
                    transient=transient,
                    persist_db=persist_db,
                )
                api_runs.extend(runs)
                raw_meta.extend(meta)
                updated_results.append(persisted_result)
            except Exception as exc:  # noqa: BLE001
                logger.warning("UCSC raw/figure persist failed: %s", exc)
                if transient is not None:
                    transient.clear_run(dossier_run_id)
                updated_results.append(result)
            continue

        api_run = ApiRun(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            source_name=result.source_name,
            endpoint_name=result.endpoint_name,
            request_url=result.request_url,
            request_params=dict(result.request_params or {}),
            status_code=result.status_code,
            success=result.success,
            error_type=result.error_type,
            error_message=result.error_message,
        )
        artifact = None
        if result.data is not None:
            try:
                if isinstance(result.data, str):
                    artifact = store.save_text(
                        dossier_run_id,
                        result.source_name,
                        result.data,
                        api_run_id=api_run.id,
                        original_url=result.request_url or None,
                        filename_hint=result.endpoint_name,
                    )
                else:
                    artifact = store.save_json(
                        dossier_run_id,
                        result.source_name,
                        result.data,
                        api_run_id=api_run.id,
                        original_url=result.request_url or None,
                        filename_hint=result.endpoint_name,
                    )
                api_run.raw_artifact_id = artifact.id
                result.raw_artifact_id = artifact.id
                raw_meta.append(artifact.model_dump(mode="json"))
            except Exception as exc:  # noqa: BLE001
                logger.warning("Raw store failed for %s: %s", result.source_name, exc)

        api_runs.append(api_run)
        updated_results.append(result)

        if persist_db:
            try:
                with session_scope() as session:
                    save_api_run(session, api_run)
                    if artifact is not None:
                        save_raw_artifact(session, artifact)
            except Exception as exc:  # noqa: BLE001
                logger.warning("DB persist failed for %s: %s", result.source_name, exc)

    return {
        **state,
        "tool_results": updated_results,
        "api_runs": api_runs,
        "raw_artifacts": raw_meta,
    }


def node_normalize_evidence(
    state: DossierState, *, persist_db: bool = True
) -> DossierState:
    """Normalize successful ToolResults into EvidenceRecords."""
    dossier_run_id = state["dossier_run_id"]
    evidence: list[EvidenceRecord] = list(state.get("evidence_records") or [])
    seen = {e.source_id for e in evidence if e.source_id}
    api_by_artifact = {
        a.raw_artifact_id: a
        for a in (state.get("api_runs") or [])
        if a.raw_artifact_id
    }
    # Fallback for payloads without artifacts: match by (source, endpoint).
    api_by_source_endpoint: dict[tuple[str, str], Any] = {}
    for a in state.get("api_runs") or []:
        key = (a.source_name or "", a.endpoint_name or "")
        api_by_source_endpoint[key] = a

    for result in state.get("tool_results") or []:
        if not result.success:
            continue
        api_run = None
        if result.raw_artifact_id:
            api_run = api_by_artifact.get(result.raw_artifact_id)
        if api_run is None:
            api_run = api_by_source_endpoint.get(
                (result.source_name or "", result.endpoint_name or "")
            )
        batch = normalize_tool_result(
            result,
            dossier_run_id=dossier_run_id,
            api_run_id=api_run.id if api_run else None,
            raw_artifact_id=result.raw_artifact_id,
        )
        for rec in batch:
            if rec.source_id and rec.source_id in seen:
                continue
            if rec.source_id:
                seen.add(rec.source_id)
            evidence.append(rec)
            if persist_db:
                try:
                    with session_scope() as session:
                        save_evidence_record(session, rec)
                except Exception as exc:  # noqa: BLE001
                    logger.warning(
                        "DB evidence persist failed for %s: %s",
                        rec.source_id,
                        exc,
                    )

    return {
        **state,
        "evidence_records": evidence,
    }


def node_index_evidence_in_chroma(
    state: DossierState, *, settings: Settings
) -> DossierState:
    """Optionally index evidence in Chroma (soft-fail; never aborts the graph)."""
    notes = list(state.get("synthesis_notes") or [])
    records = list(state.get("evidence_records") or [])
    if not records:
        notes.append("No evidence records available for Chroma indexing.")
        return {**state, "synthesis_notes": notes}

    try:
        from gene_dossier.retrieval import index_evidence_in_chroma

        status = index_evidence_in_chroma(records, settings=settings)
        if status.available:
            notes.append(
                f"Chroma indexing complete: {status.indexed_count} evidence "
                f"records indexed ({status.backend})."
            )
        else:
            err = status.error or "unknown error"
            notes.append(f"Chroma indexing unavailable: {err}")
    except Exception as exc:  # noqa: BLE001
        logger.warning("Chroma indexing failed: %s", exc)
        notes.append(f"Chroma indexing failed: {exc}")

    return {**state, "synthesis_notes": notes}


def node_build_report_sections(
    state: DossierState, *, settings: Settings, force_deterministic: bool = True
) -> DossierState:
    """Synthesize CHDI-style sections + claims from evidence."""
    synthesis: SynthesisResult = synthesize_dossier(
        dossier_run_id=state["dossier_run_id"],
        gene_symbol=state["gene_symbol"],
        evidence_records=state.get("evidence_records") or [],
        settings=settings,
        force_deterministic=force_deterministic,
    )
    return {
        **state,
        "sections": list(synthesis.sections),
        "claims": list(synthesis.claims),
        "synthesis_mode": synthesis.mode,
        "synthesis_notes": list(state.get("synthesis_notes") or [])
        + list(synthesis.notes),
    }


def node_verify_claims(state: DossierState) -> DossierState:
    """Rule-based claim verification against evidence records."""
    results = verify_claims(
        state.get("claims") or [],
        state.get("evidence_records") or [],
    )
    return {**state, "verification_results": list(results)}


def _coverage_updates_from_state(state: DossierState) -> list[SourceCoverageResult]:
    """Build per-source coverage updates from tool results + evidence counts."""
    dossier_run_id = state["dossier_run_id"]
    evidence = state.get("evidence_records") or []
    counts: dict[str, int] = {}
    for rec in evidence:
        counts[rec.source_name] = counts.get(rec.source_name, 0) + 1

    artifact_by_source: dict[str, str] = {}
    for meta in state.get("raw_artifacts") or []:
        name = meta.get("source_name")
        path = meta.get("file_path")
        if name and path and name not in artifact_by_source:
            artifact_by_source[name] = path

    updates: list[SourceCoverageResult] = []
    for result in state.get("tool_results") or []:
        src_def = get_source(result.source_name)
        sections = list(src_def.report_sections) if src_def else []
        if result.error_type == "requires_key":
            status = SourceStatus.requires_key
        elif result.error_type in {
            "access_forbidden",
            "source_unavailable",
            "endpoint_unavailable",
        }:
            status = SourceStatus.deferred
        elif result.success and counts.get(result.source_name, 0) > 0:
            status = SourceStatus.success
        elif result.success:
            status = SourceStatus.partial
        else:
            status = SourceStatus.failed
        updates.append(
            SourceCoverageResult(
                dossier_run_id=dossier_run_id,
                source_name=result.source_name,
                status=status,
                raw_artifact_path=artifact_by_source.get(result.source_name),
                evidence_record_count=counts.get(result.source_name, 0),
                error_message=result.error_message,
                report_sections_supported=sections,
                notes=src_def.notes if src_def else None,
            )
        )
    return updates


def node_render_outputs(
    state: DossierState,
    *,
    settings: Settings,
    output_dir: Path | None = None,
    write_rancho: bool = True,
    write_pdf: bool = True,
    persist_db: bool = True,
) -> DossierState:
    """Write coverage, debug markdown, and polished Rancho report outputs."""
    out = Path(output_dir) if output_dir is not None else settings.output_path
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = dict(state.get("output_paths") or {})
    errors = list(state.get("errors") or [])
    coverage = list(state.get("coverage") or [])

    updates = _coverage_updates_from_state(state)
    try:
        coverage, cov_paths = build_and_write_coverage(
            state["dossier_run_id"],
            gene_symbol=state["gene_symbol"],
            updates=updates,
            settings=settings,
            output_dir=out,
            session=None,
        )
        if persist_db:
            with session_scope() as session:
                persist_coverage(session, coverage)
        paths["coverage_markdown"] = str(cov_paths["markdown"])
        paths["coverage_json"] = str(cov_paths["json"])
    except Exception as exc:  # noqa: BLE001
        logger.warning("Coverage write failed: %s", exc)
        errors.append(f"coverage: {exc}")

    try:
        debug_paths = write_dossier_report(
            gene_symbol=state["gene_symbol"],
            dossier_run_id=state["dossier_run_id"],
            sections=list(state.get("sections") or []),
            claims=state.get("claims") or [],
            coverage=coverage,
            verification_results=state.get("verification_results") or [],
            synthesis_mode=state.get("synthesis_mode"),
            synthesis_notes=state.get("synthesis_notes") or [],
            output_dir=out,
            settings=settings,
        )
        paths["debug_markdown"] = str(debug_paths["markdown"])
        paths["debug_json"] = str(debug_paths["json"])
    except Exception as exc:  # noqa: BLE001
        logger.warning("Debug report write failed: %s", exc)
        errors.append(f"debug_report: {exc}")

    if write_rancho:
        try:
            chromosome = (state.get("gene_ids") or {}).get("chromosome")
            _doc, rancho_paths = build_and_write_rancho_report(
                dossier_run_id=state["dossier_run_id"],
                gene_symbol=state["gene_symbol"],
                evidence_records=state.get("evidence_records") or [],
                report_sections=state.get("sections") or None,
                curator="Gene Dossier Platform",
                report_date=datetime.now(timezone.utc).strftime("%Y-%m-%d"),
                chromosome=str(chromosome) if chromosome else None,
                output_dir=out,
                settings=settings,
                include_endnotes=False,
                write_pdf=write_pdf,
                show_cover_logos=False,
            )
            for key, path in rancho_paths.items():
                paths[f"rancho_{key}"] = str(path)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Rancho report write failed: %s", exc)
            errors.append(f"rancho_report: {exc}")

    if persist_db:
        try:
            with session_scope() as session:
                run = DossierRun(
                    id=state["dossier_run_id"],
                    gene_symbol=state["gene_symbol"],
                    official_symbol=state.get("official_symbol"),
                    run_type="full_api_pass",
                    status="completed",
                    completed_at=utcnow(),
                )
                save_dossier_run(session, run)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to mark dossier run completed: %s", exc)

    return {
        **state,
        "coverage": coverage,
        "output_paths": paths,
        "errors": errors,
        "status": "completed",
    }


def build_dossier_graph(
    *,
    settings: Settings | None = None,
    output_dir: str | Path | None = None,
    sources: Iterable[str] | None = None,
    call_network: bool = True,
    force_deterministic: bool = True,
    write_rancho: bool = True,
    write_pdf: bool = True,
    persist_db: bool = True,
):
    """Compile the LangGraph StateGraph for a full API pass."""
    from langgraph.graph import END, START, StateGraph

    cfg = settings or get_settings()
    out = Path(output_dir) if output_dir is not None else None
    source_list = list(sources) if sources is not None else None
    transient = WorkflowTransientContext()

    graph = StateGraph(DossierState)
    graph.add_node(
        "create_dossier_run",
        lambda s: node_create_dossier_run(s, settings=cfg, persist_db=persist_db),
    )
    graph.add_node(
        "resolve_gene_identity",
        lambda s: node_resolve_gene_identity(
            s, settings=cfg, call_network=call_network
        ),
    )
    graph.add_node(
        "call_source_clients",
        lambda s: node_call_source_clients(
            s,
            settings=cfg,
            call_network=call_network,
            sources=source_list,
            transient=transient,
        ),
    )
    graph.add_node(
        "save_raw_artifacts",
        lambda s: node_save_raw_artifacts(
            s,
            settings=cfg,
            persist_db=persist_db,
            transient=transient,
        ),
    )
    graph.add_node(
        "normalize_evidence",
        lambda s: node_normalize_evidence(s, persist_db=persist_db),
    )
    graph.add_node(
        "index_evidence_in_chroma",
        lambda s: node_index_evidence_in_chroma(s, settings=cfg),
    )
    graph.add_node(
        "build_report_sections",
        lambda s: node_build_report_sections(
            s, settings=cfg, force_deterministic=force_deterministic
        ),
    )
    graph.add_node("verify_claims", node_verify_claims)
    graph.add_node(
        "render_outputs",
        lambda s: node_render_outputs(
            s,
            settings=cfg,
            output_dir=out,
            write_rancho=write_rancho,
            write_pdf=write_pdf,
            persist_db=persist_db,
        ),
    )

    graph.add_edge(START, "create_dossier_run")
    graph.add_edge("create_dossier_run", "resolve_gene_identity")
    graph.add_edge("resolve_gene_identity", "call_source_clients")
    graph.add_edge("call_source_clients", "save_raw_artifacts")
    graph.add_edge("save_raw_artifacts", "normalize_evidence")
    graph.add_edge("normalize_evidence", "index_evidence_in_chroma")
    graph.add_edge("index_evidence_in_chroma", "build_report_sections")
    graph.add_edge("build_report_sections", "verify_claims")
    graph.add_edge("verify_claims", "render_outputs")
    graph.add_edge("render_outputs", END)
    return graph.compile()


def run_gene_dossier_full_api_pass(
    gene_symbol: str,
    *,
    settings: Settings | None = None,
    output_dir: str | Path | None = None,
    dossier_run_id: str | None = None,
    sources: Iterable[str] | None = None,
    call_network: bool = True,
    force_deterministic: bool = True,
    write_rancho: bool = True,
    write_pdf: bool = True,
    persist_db: bool = True,
    preloaded_tool_results: Iterable[ToolResult] | None = None,
    gene_ids: dict[str, Any] | None = None,
) -> DossierPassResult:
    """Run the full LangGraph dossier pass for ``gene_symbol``.

    Set ``call_network=False`` and provide ``preloaded_tool_results`` for offline
    / unit tests. Without LLM keys, synthesis is deterministic regardless of
    ``force_deterministic``.
    """
    cfg = settings or get_settings()
    cfg.raw_data_path.mkdir(parents=True, exist_ok=True)
    cfg.output_path.mkdir(parents=True, exist_ok=True)

    initial: DossierState = {
        "gene_symbol": gene_symbol.strip(),
        "gene_ids": dict(gene_ids or {}),
        "tool_results": list(preloaded_tool_results or []),
        "errors": [],
        "status": "created",
    }
    if dossier_run_id:
        initial["dossier_run_id"] = dossier_run_id

    app = build_dossier_graph(
        settings=cfg,
        output_dir=output_dir,
        sources=sources,
        call_network=call_network,
        force_deterministic=force_deterministic,
        write_rancho=write_rancho,
        write_pdf=write_pdf,
        persist_db=persist_db,
    )
    final_state: DossierState = app.invoke(initial)
    return DossierPassResult.from_state(final_state)


__all__ = [
    "DossierState",
    "DossierPassResult",
    "CLIENT_DISPATCH",
    "normalize_tool_result",
    "extract_gene_ids_from_tool_result",
    "build_dossier_graph",
    "run_gene_dossier_full_api_pass",
    "node_create_dossier_run",
    "node_resolve_gene_identity",
    "node_call_source_clients",
    "node_save_raw_artifacts",
    "node_normalize_evidence",
    "node_index_evidence_in_chroma",
    "node_build_report_sections",
    "node_verify_claims",
    "node_render_outputs",
    "WorkflowTransientContext",
    "bind_workflow_transient",
    "get_active_workflow_transient",
    "reset_workflow_transient",
    "coverage_updates_from_state",
]


def coverage_updates_from_state(state: DossierState) -> list[SourceCoverageResult]:
    """Public alias for read-only coverage derivation (section-bundle audit)."""
    return _coverage_updates_from_state(state)
