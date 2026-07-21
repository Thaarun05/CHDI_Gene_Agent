"""Normalize model-organism ToolResults into EvidenceRecords.

Consumes successful MouseMine ``fetch_mouse_annotations`` payloads. Does **not**
call the network.

Rules:
- Do not invent alleles, phenotypes, or stock relationships
- Mouse / knockout evidence is grade D
- Phenotype terms come only from ontologyAnnotations fields in the payload
- MGI gene identity routes to Homologues; alleles/phenotypes/stocks to
  Knockouts / model phenotypes
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

SECTION_HOMOLOGUES = "Homologues"
SECTION_KNOCKOUTS = "Knockouts / model phenotypes"

DEFAULT_MAX_ALLELES = 100
DEFAULT_MAX_PHENOTYPES = 200
DEFAULT_MAX_STOCKS = 100


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _row_get(row: dict[str, Any], *keys: str) -> Any:
    """Return the first present value among InterMine dotted or short keys."""
    for key in keys:
        if key in row and row[key] not in (None, ""):
            return row[key]
    return None


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    assertion_type: AssertionType,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    section: str,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one MouseMine EvidenceRecord."""
    source_id = make_source_id("MouseMine", gene_symbol, assertion_type, key)
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=section,
        subsection=subsection,
        source_name="MouseMine",
        source_type=SourceType.model_organism_database,
        assertion_type=assertion_type,
        fact_type=fact_type,
        organism="Mus musculus",
        species="Mus musculus",
        taxon_id=10090,
        evidence_grade=EvidenceGrade.D,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def normalize_mousemine(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
    max_alleles: int = DEFAULT_MAX_ALLELES,
    max_phenotypes: int = DEFAULT_MAX_PHENOTYPES,
    max_stocks: int = DEFAULT_MAX_STOCKS,
) -> list[EvidenceRecord]:
    """Normalize MouseMine ``fetch_mouse_annotations`` payloads."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    mgi_id = data.get("mgi_id")
    if not gene_symbol:
        # Fall back to MGI feature symbol from gene summaries when present.
        for row in data.get("gene_summaries") or []:
            if isinstance(row, dict) and row.get("symbol"):
                gene_symbol = str(row["symbol"])
                break
    if not gene_symbol or not mgi_id:
        return []

    records: list[EvidenceRecord] = []
    common = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }

    # MGI ortholog / mouse gene identifier.
    mouse_symbol = None
    mouse_name = None
    for row in data.get("gene_summaries") or []:
        if isinstance(row, dict) and (
            str(row.get("mgi_id")) == str(mgi_id) or row.get("mgi_id")
        ):
            mouse_symbol = row.get("symbol")
            mouse_name = row.get("name")
            if str(row.get("mgi_id")) == str(mgi_id):
                break

    records.append(
        _record(
            **common,
            assertion_type=AssertionType.orthology,
            fact_type="mgi_gene_id",
            key=str(mgi_id),
            value={
                "mgi_id": mgi_id,
                "ncbi_gene_number": data.get("ncbi_gene_number"),
                "mouse_symbol": mouse_symbol,
                "mouse_name": mouse_name,
            },
            display_text=(
                f"{gene_symbol} mouse MGI ID is {mgi_id}"
                + (f" ({mouse_symbol})." if mouse_symbol else ".")
            ),
            section=SECTION_HOMOLOGUES,
            subsection="MouseMine MGI gene",
            confidence_notes=(
                "MouseMine MGI identifier for the mouse gene/ortholog; mouse "
                "model evidence is grade D."
            ),
        )
    )

    allele_rows = data.get("allele_rows") or []
    if not isinstance(allele_rows, list):
        allele_rows = []
    for idx, row in enumerate(allele_rows[: max(0, max_alleles)], start=1):
        if not isinstance(row, dict):
            continue
        allele_id = _row_get(
            row, "Allele.primaryIdentifier", "primaryIdentifier"
        )
        allele_symbol = _row_get(row, "Allele.symbol", "symbol")
        if not allele_id and not allele_symbol:
            continue
        allele_type = _row_get(row, "Allele.alleleType", "alleleType")
        allele_name = _row_get(row, "Allele.name", "name")
        key = str(allele_id or f"{allele_symbol}-{idx}")
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.knockout_phenotype,
                fact_type="mouse_allele",
                key=key,
                value={
                    "mgi_id": mgi_id,
                    "allele_id": allele_id,
                    "allele_symbol": allele_symbol,
                    "allele_name": allele_name,
                    "allele_type": allele_type,
                    "feature_mgi_id": _row_get(
                        row,
                        "Allele.feature.primaryIdentifier",
                        "feature.primaryIdentifier",
                    ),
                    "feature_symbol": _row_get(
                        row, "Allele.feature.symbol", "feature.symbol"
                    ),
                },
                display_text=(
                    f"{gene_symbol} mouse allele "
                    f"{allele_symbol or allele_id}"
                    + (f" ({allele_type})." if allele_type else ".")
                ),
                section=SECTION_KNOCKOUTS,
                subsection="MouseMine alleles",
            )
        )

    pheno_rows = data.get("phenotype_rows") or []
    if not isinstance(pheno_rows, list):
        pheno_rows = []
    for idx, row in enumerate(pheno_rows[: max(0, max_phenotypes)], start=1):
        if not isinstance(row, dict):
            continue
        allele_id = _row_get(
            row, "Allele.primaryIdentifier", "primaryIdentifier"
        )
        allele_symbol = _row_get(row, "Allele.symbol", "symbol")
        term_id = _row_get(
            row,
            "Allele.ontologyAnnotations.ontologyTerm.identifier",
            "ontologyAnnotations.ontologyTerm.identifier",
            "ontologyTerm.identifier",
        )
        term_name = _row_get(
            row,
            "Allele.ontologyAnnotations.ontologyTerm.name",
            "ontologyAnnotations.ontologyTerm.name",
            "ontologyTerm.name",
        )
        if not term_id and not term_name:
            continue
        key = "-".join(
            str(x)
            for x in (allele_id or allele_symbol or "allele", term_id or term_name, idx)
            if x is not None and str(x) != ""
        )
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.knockout_phenotype,
                fact_type="mouse_allele_phenotype",
                key=key,
                value={
                    "mgi_id": mgi_id,
                    "allele_id": allele_id,
                    "allele_symbol": allele_symbol,
                    "allele_type": _row_get(row, "Allele.alleleType", "alleleType"),
                    "ontology_term_id": term_id,
                    "ontology_term_name": term_name,
                    "caveat": (
                        "MouseMine ontology phenotype annotation; mouse model "
                        "evidence, not human disease proof."
                    ),
                },
                display_text=(
                    f"{gene_symbol} mouse allele "
                    f"{allele_symbol or allele_id} phenotype "
                    f"{term_name or term_id}."
                ),
                section=SECTION_KNOCKOUTS,
                subsection="MouseMine allele phenotypes",
                confidence_notes=(
                    "MouseMine ontology phenotype annotation; mouse model "
                    "evidence, not human disease proof."
                ),
                manual_review_required=True,
            )
        )

    stock_rows = data.get("stock_rows") or []
    if not isinstance(stock_rows, list):
        stock_rows = []
    for idx, row in enumerate(stock_rows[: max(0, max_stocks)], start=1):
        if not isinstance(row, dict):
            continue
        allele_id = _row_get(
            row, "Allele.primaryIdentifier", "primaryIdentifier"
        )
        allele_symbol = _row_get(row, "Allele.symbol", "symbol")
        stock_id = _row_get(
            row,
            "Allele.carriedBy.primaryIdentifier",
            "carriedBy.primaryIdentifier",
        )
        stock_symbol = _row_get(
            row, "Allele.carriedBy.symbol", "carriedBy.symbol"
        )
        stock_name = _row_get(row, "Allele.carriedBy.name", "carriedBy.name")
        if not stock_id and not stock_symbol and not stock_name:
            continue
        key = "-".join(
            str(x)
            for x in (
                allele_id or allele_symbol or "allele",
                stock_id or stock_symbol or stock_name,
                idx,
            )
            if x is not None and str(x) != ""
        )
        records.append(
            _record(
                **common,
                assertion_type=AssertionType.knockout_phenotype,
                fact_type="mouse_allele_stock",
                key=key,
                value={
                    "mgi_id": mgi_id,
                    "allele_id": allele_id,
                    "allele_symbol": allele_symbol,
                    "stock_id": stock_id,
                    "stock_symbol": stock_symbol,
                    "stock_name": stock_name,
                },
                display_text=(
                    f"{gene_symbol} mouse allele "
                    f"{allele_symbol or allele_id} carried by "
                    f"{stock_symbol or stock_name or stock_id}."
                ),
                section=SECTION_KNOCKOUTS,
                subsection="MouseMine stocks",
            )
        )

    return records


def normalize_model_organisms(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch model-organism normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    if source == "MouseMine":
        return normalize_mousemine(
            tool_result,
            dossier_run_id=dossier_run_id,
            api_run_id=api_run_id,
            raw_artifact_id=raw_artifact_id,
        )
    return []


__all__ = [
    "DEFAULT_MAX_ALLELES",
    "DEFAULT_MAX_PHENOTYPES",
    "DEFAULT_MAX_STOCKS",
    "normalize_mousemine",
    "normalize_model_organisms",
]
