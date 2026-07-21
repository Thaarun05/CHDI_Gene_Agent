"""MouseMine (InterMine) client for mouse gene / allele / phenotype data.

Resolves mouse NCBI Gene → MGI ID, then fetches alleles, ontology phenotype
annotations, and stocks. Does **not** normalize into evidence records — that
belongs in ``normalize/model_organisms.py``.

Key endpoint (validated)::

    GET https://www.mousemine.org/mousemine/service/query/results
        ?format=json&query=<XML>

Queries:
- Gene lookup by ``Gene.ncbiGeneNumber`` (SREBF2 mouse Entrez ``20788`` → ``MGI:107585``)
- Alleles by ``Allele.feature.primaryIdentifier`` (MGI ID)
- Allele phenotypes via ``Allele.ontologyAnnotations.ontologyTerm.*``
  (not ``phenotypeAnnotations`` — that path is not in the MouseMine model)
- Stocks / carriedBy by MGI ID

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode
from xml.sax.saxutils import escape as xml_escape

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "MouseMine"
MOUSEMINE_RESULTS_URL = (
    "https://www.mousemine.org/mousemine/service/query/results"
)

# Validated SREBF2 mouse anchors.
DEFAULT_MOUSE_NCBI_GENE = "20788"
DEFAULT_MGI_ID = "MGI:107585"

GENE_LOOKUP_VIEWS = (
    "Gene.primaryIdentifier",
    "Gene.symbol",
    "Gene.name",
    "Gene.organism.name",
    "Gene.ncbiGeneNumber",
)

ALLELE_VIEWS = (
    "Allele.primaryIdentifier",
    "Allele.symbol",
    "Allele.name",
    "Allele.alleleType",
    "Allele.feature.primaryIdentifier",
    "Allele.feature.symbol",
)

PHENOTYPE_VIEWS = (
    "Allele.primaryIdentifier",
    "Allele.symbol",
    "Allele.name",
    "Allele.alleleType",
    "Allele.ontologyAnnotations.ontologyTerm.identifier",
    "Allele.ontologyAnnotations.ontologyTerm.name",
)

STOCKS_VIEWS = (
    "Allele.primaryIdentifier",
    "Allele.symbol",
    "Allele.name",
    "Allele.alleleType",
    "Allele.carriedBy.primaryIdentifier",
    "Allele.carriedBy.symbol",
    "Allele.carriedBy.name",
)


def _tool_result(
    *,
    endpoint_name: str,
    gene_symbol: str,
    request_url: str,
    request_params: dict[str, Any],
    success: bool,
    status_code: int | None = None,
    data: Any | None = None,
    error_type: str | None = None,
    error_message: str | None = None,
) -> ToolResult:
    """Build a uniform :class:`ToolResult` for this source."""
    return ToolResult(
        source_name=SOURCE_NAME,
        endpoint_name=endpoint_name,
        success=success,
        gene_symbol=gene_symbol,
        request_url=request_url,
        request_params=request_params,
        status_code=status_code,
        data=data,
        error_type=error_type,
        error_message=error_message,
    )


def _build_query_xml(
    *,
    view: tuple[str, ...] | list[str],
    constraint_path: str,
    constraint_value: str,
    sort_order: str,
) -> str:
    """Build a MouseMine PathQuery XML string (validated shape)."""
    view_str = " ".join(view)
    value = xml_escape(str(constraint_value), {"\"": "&quot;", "'": "&apos;"})
    path = xml_escape(constraint_path, {"\"": "&quot;", "'": "&apos;"})
    sort = xml_escape(sort_order, {"\"": "&quot;", "'": "&apos;"})
    return (
        f'<query model="genomic" view="{view_str}" sortOrder="{sort}">'
        f'<constraint path="{path}" op="=" value="{value}"/>'
        f"</query>"
    )


def build_gene_lookup_query(ncbi_gene_number: str | int) -> str:
    """XML query: gene by NCBI Gene number."""
    return _build_query_xml(
        view=GENE_LOOKUP_VIEWS,
        constraint_path="Gene.ncbiGeneNumber",
        constraint_value=str(ncbi_gene_number),
        sort_order="Gene.primaryIdentifier asc",
    )


def build_alleles_query(mgi_id: str) -> str:
    """XML query: alleles for an MGI gene ID."""
    return _build_query_xml(
        view=ALLELE_VIEWS,
        constraint_path="Allele.feature.primaryIdentifier",
        constraint_value=mgi_id,
        sort_order="Allele.symbol asc",
    )


def build_allele_phenotypes_query(mgi_id: str) -> str:
    """XML query: allele ontology phenotype annotations for an MGI gene ID."""
    return _build_query_xml(
        view=PHENOTYPE_VIEWS,
        constraint_path="Allele.feature.primaryIdentifier",
        constraint_value=mgi_id,
        sort_order="Allele.symbol asc",
    )


def build_stocks_query(mgi_id: str) -> str:
    """XML query: allele stocks / carriedBy for an MGI gene ID."""
    return _build_query_xml(
        view=STOCKS_VIEWS,
        constraint_path="Allele.feature.primaryIdentifier",
        constraint_value=mgi_id,
        sort_order="Allele.symbol asc",
    )


def rows_as_dicts(payload: Any, fallback_views: tuple[str, ...] | list[str]) -> list[dict[str, Any]]:
    """Convert InterMine JSON ``results`` rows into dicts keyed by view path."""
    if not isinstance(payload, dict):
        return []
    results = payload.get("results")
    if not isinstance(results, list):
        return []

    views = payload.get("views") or payload.get("columnHeaders") or list(fallback_views)
    if not isinstance(views, list):
        views = list(fallback_views)

    out: list[dict[str, Any]] = []
    for row in results:
        if isinstance(row, dict):
            out.append(dict(row))
            continue
        if not isinstance(row, (list, tuple)):
            continue
        mapped: dict[str, Any] = {}
        for idx, value in enumerate(row):
            key = str(views[idx]) if idx < len(views) else f"col_{idx}"
            mapped[key] = value
        out.append(mapped)
    return out


def prefer_mgi_id(rows: list[dict[str, Any]]) -> str | None:
    """Pick an MGI primaryIdentifier from gene-lookup rows."""
    for row in rows:
        mgi = (
            row.get("Gene.primaryIdentifier")
            or row.get("primaryIdentifier")
            or row.get("mgi_id")
        )
        if mgi and str(mgi).startswith("MGI:"):
            return str(mgi)
    return None


def summarize_gene_row(row: dict[str, Any]) -> dict[str, Any]:
    """Light gene-lookup summary (not evidence)."""
    return {
        "mgi_id": row.get("Gene.primaryIdentifier") or row.get("primaryIdentifier"),
        "symbol": row.get("Gene.symbol") or row.get("symbol"),
        "name": row.get("Gene.name") or row.get("name"),
        "organism": row.get("Gene.organism.name") or row.get("organism"),
        "ncbi_gene_number": row.get("Gene.ncbiGeneNumber") or row.get("ncbiGeneNumber"),
    }


def _query_results(
    *,
    endpoint_name: str,
    gene_symbol: str,
    query_xml: str,
    settings: Settings,
    extra_params: dict[str, Any] | None = None,
) -> ToolResult:
    """GET MouseMine query/results and return :class:`ToolResult`."""
    params = {"format": "json", "query": query_xml}
    request_params = {
        "format": "json",
        "query": query_xml,
        **(extra_params or {}),
    }
    request_url = f"{MOUSEMINE_RESULTS_URL}?{urlencode(params)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(MOUSEMINE_RESULTS_URL, params=params)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            # InterMine may return 200 with an error object for bad PathQueries.
            if isinstance(payload, dict) and payload.get("error"):
                return _tool_result(
                    endpoint_name=endpoint_name,
                    gene_symbol=gene_symbol,
                    request_url=request_url,
                    request_params=request_params,
                    success=False,
                    status_code=response.status_code,
                    data=payload,
                    error_type="query_error",
                    error_message=str(payload.get("error")),
                )
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=request_params,
            success=False,
            status_code=response.status_code,
            data=payload,
            error_type="http_error",
            error_message=f"HTTP {response.status_code}",
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def gene_lookup(
    ncbi_gene_number: str | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Look up MouseMine gene rows by NCBI Gene number."""
    cfg = settings or get_settings()
    ncbi = str(ncbi_gene_number).strip()
    query_xml = build_gene_lookup_query(ncbi)
    return _query_results(
        endpoint_name="gene_lookup",
        gene_symbol=gene_symbol or ncbi,
        query_xml=query_xml,
        settings=cfg,
        extra_params={"ncbi_gene_number": ncbi},
    )


def alleles(
    mgi_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch alleles for an MGI gene ID."""
    cfg = settings or get_settings()
    mgi = mgi_id.strip()
    return _query_results(
        endpoint_name="alleles",
        gene_symbol=gene_symbol or mgi,
        query_xml=build_alleles_query(mgi),
        settings=cfg,
        extra_params={"mgi_id": mgi},
    )


def allele_phenotypes(
    mgi_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch allele ontology phenotype annotations for an MGI gene ID."""
    cfg = settings or get_settings()
    mgi = mgi_id.strip()
    return _query_results(
        endpoint_name="allele_phenotypes",
        gene_symbol=gene_symbol or mgi,
        query_xml=build_allele_phenotypes_query(mgi),
        settings=cfg,
        extra_params={"mgi_id": mgi},
    )


def stocks_carried_by(
    mgi_id: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch allele stocks / carriedBy for an MGI gene ID."""
    cfg = settings or get_settings()
    mgi = mgi_id.strip()
    return _query_results(
        endpoint_name="stocks_carried_by",
        gene_symbol=gene_symbol or mgi,
        query_xml=build_stocks_query(mgi),
        settings=cfg,
        extra_params={"mgi_id": mgi},
    )


def fetch_mouse_annotations(
    *,
    ncbi_gene_number: str | int | None = None,
    mgi_id: str | None = None,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Resolve MGI (if needed) and fetch alleles, phenotypes, and stocks.

    Provide ``mgi_id`` and/or ``ncbi_gene_number``. When only NCBI is given,
    runs gene lookup first. On success, ``data`` includes raw payloads plus
    row dicts for each query.

    Never raises.
    """
    cfg = settings or get_settings()
    resolved_mgi = mgi_id.strip() if mgi_id else None
    gene_payload: Any = None
    gene_rows: list[dict[str, Any]] = []
    gene_summaries: list[dict[str, Any]] = []
    lookup_url = MOUSEMINE_RESULTS_URL
    lookup_params: dict[str, Any] = {}

    if not resolved_mgi:
        if ncbi_gene_number is None:
            return _tool_result(
                endpoint_name="fetch_mouse_annotations",
                gene_symbol=gene_symbol,
                request_url=MOUSEMINE_RESULTS_URL,
                request_params={},
                success=False,
                error_type="invalid_request",
                error_message="Provide mgi_id and/or ncbi_gene_number",
            )
        lookup = gene_lookup(
            ncbi_gene_number, gene_symbol=gene_symbol, settings=cfg
        )
        lookup_url = lookup.request_url
        lookup_params = lookup.request_params
        if not lookup.success:
            return _tool_result(
                endpoint_name="fetch_mouse_annotations",
                gene_symbol=gene_symbol or str(ncbi_gene_number),
                request_url=lookup.request_url,
                request_params=lookup.request_params,
                success=False,
                status_code=lookup.status_code,
                data={"gene_lookup": lookup.data},
                error_type=lookup.error_type or "gene_lookup_failed",
                error_message=lookup.error_message or "MouseMine gene lookup failed",
            )
        gene_payload = lookup.data
        gene_rows = rows_as_dicts(gene_payload, GENE_LOOKUP_VIEWS)
        gene_summaries = [summarize_gene_row(r) for r in gene_rows]
        resolved_mgi = prefer_mgi_id(gene_rows)
        if not resolved_mgi:
            return _tool_result(
                endpoint_name="fetch_mouse_annotations",
                gene_symbol=gene_symbol or str(ncbi_gene_number),
                request_url=lookup.request_url,
                request_params=lookup.request_params,
                success=False,
                status_code=lookup.status_code,
                data={
                    "ncbi_gene_number": str(ncbi_gene_number),
                    "gene_lookup": gene_payload,
                    "gene_rows": gene_rows,
                    "gene_summaries": gene_summaries,
                },
                error_type="no_results",
                error_message=f"No MGI ID for NCBI Gene {ncbi_gene_number}",
            )

    allele_res = alleles(resolved_mgi, gene_symbol=gene_symbol, settings=cfg)
    pheno_res = allele_phenotypes(
        resolved_mgi, gene_symbol=gene_symbol, settings=cfg
    )
    stocks_res = stocks_carried_by(
        resolved_mgi, gene_symbol=gene_symbol, settings=cfg
    )

    # Prefer overall success only when all three MGI-scoped queries succeed.
    # Surface the first failure; still attach any successful sibling payloads.
    failures = [
        ("alleles", allele_res),
        ("allele_phenotypes", pheno_res),
        ("stocks_carried_by", stocks_res),
    ]
    for name, res in failures:
        if not res.success:
            return _tool_result(
                endpoint_name="fetch_mouse_annotations",
                gene_symbol=gene_symbol or resolved_mgi,
                request_url=res.request_url,
                request_params=res.request_params,
                success=False,
                status_code=res.status_code,
                data={
                    "ncbi_gene_number": (
                        str(ncbi_gene_number) if ncbi_gene_number is not None else None
                    ),
                    "mgi_id": resolved_mgi,
                    "gene_lookup": gene_payload,
                    "gene_rows": gene_rows,
                    "gene_summaries": gene_summaries,
                    "alleles": allele_res.data,
                    "allele_phenotypes": pheno_res.data,
                    "stocks_carried_by": stocks_res.data,
                    "failed_step": name,
                },
                error_type=res.error_type or f"{name}_failed",
                error_message=res.error_message or f"MouseMine {name} failed",
            )

    allele_rows = rows_as_dicts(allele_res.data, ALLELE_VIEWS)
    pheno_rows = rows_as_dicts(pheno_res.data, PHENOTYPE_VIEWS)
    stock_rows = rows_as_dicts(stocks_res.data, STOCKS_VIEWS)

    return _tool_result(
        endpoint_name="fetch_mouse_annotations",
        gene_symbol=gene_symbol or resolved_mgi,
        request_url=allele_res.request_url or lookup_url,
        request_params={
            "ncbi_gene_number": (
                str(ncbi_gene_number) if ncbi_gene_number is not None else None
            ),
            "mgi_id": resolved_mgi,
            **lookup_params,
        },
        success=True,
        status_code=allele_res.status_code,
        data={
            "ncbi_gene_number": (
                str(ncbi_gene_number) if ncbi_gene_number is not None else None
            ),
            "mgi_id": resolved_mgi,
            "gene_lookup": gene_payload,
            "gene_rows": gene_rows,
            "gene_summaries": gene_summaries,
            "alleles": allele_res.data,
            "allele_rows": allele_rows,
            "allele_count": len(allele_rows),
            "allele_phenotypes": pheno_res.data,
            "phenotype_rows": pheno_rows,
            "phenotype_count": len(pheno_rows),
            "stocks_carried_by": stocks_res.data,
            "stock_rows": stock_rows,
            "stock_count": len(stock_rows),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "MOUSEMINE_RESULTS_URL",
    "DEFAULT_MOUSE_NCBI_GENE",
    "DEFAULT_MGI_ID",
    "GENE_LOOKUP_VIEWS",
    "ALLELE_VIEWS",
    "PHENOTYPE_VIEWS",
    "STOCKS_VIEWS",
    "build_gene_lookup_query",
    "build_alleles_query",
    "build_allele_phenotypes_query",
    "build_stocks_query",
    "rows_as_dicts",
    "prefer_mgi_id",
    "summarize_gene_row",
    "gene_lookup",
    "alleles",
    "allele_phenotypes",
    "stocks_carried_by",
    "fetch_mouse_annotations",
]
