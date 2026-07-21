"""NIH RePORTER client (projects search).

Searches NIH-funded projects mentioning a gene (exact text search first;
optional broader pathway search only when terms are explicit or gene-mapped).
Does **not** normalize into evidence records — that belongs in
``normalize/grants.py``.

Key endpoint (validated)::

    POST https://api.reporter.nih.gov/v2/projects/search

Use exact search first. Broader pathway search only when title/abstract/terms
clearly connect to the gene biology (do not apply SREBF2 pathway terms globally).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "NIH RePORTER"
PROJECTS_SEARCH_URL = "https://api.reporter.nih.gov/v2/projects/search"

DEFAULT_LIMIT = 100
DEFAULT_OFFSET = 0
DEFAULT_SORT_FIELD = "project_start_date"
DEFAULT_SORT_ORDER = "desc"
DEFAULT_SEARCH_FIELD = "projecttitle,abstracttext,terms"
DEFAULT_OPERATOR = "or"

# Validated include_fields from CHDI_Data_APIs_Gene_Report_SREBF2.md §15.1 / §15.2.
DEFAULT_INCLUDE_FIELDS = [
    "ProjectTitle",
    "ProjectNum",
    "CoreProjectNum",
    "OrganizationName",
    "BudgetStartDate",
    "BudgetEndDate",
    "AgencyICAdmin",
    "AgencyICFundings",
    "FiscalYear",
    "AwardAmount",
    "PrincipalInvestigators",
    "ProjectDetailUrl",
    "AbstractText",
    "Terms",
]

# Gene-specific exact search tokens (quoted phrases kept intact).
GENE_SPECIFIC_EXACT_TERMS: dict[str, list[str]] = {
    "SREBF2": [
        "SREBF2",
        "SREBP2",
        "SREBP-2",
        "Srebf2",
        '"sterol regulatory element binding protein 2"',
        '"sterol regulatory element-binding protein 2"',
    ],
}

# Gene-specific broader pathway terms — only used when mapped or explicitly passed.
GENE_SPECIFIC_BROADER_TERMS: dict[str, list[str]] = {
    "SREBF2": [
        "SREBP",
        "cholesterol",
        "mevalonate",
        "sterol",
        "lipid metabolism",
        "cholesterol biosynthesis",
    ],
}


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


def _join_search_terms(terms: list[str]) -> str:
    """Join search tokens with spaces (quoted phrases preserved as-is)."""
    out: list[str] = []
    seen: set[str] = set()
    for term in terms:
        cleaned = str(term).strip()
        if not cleaned:
            continue
        key = cleaned.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(cleaned)
    return " ".join(out)


def build_exact_search_text(
    gene_symbol: str,
    aliases: list[str] | None = None,
) -> str:
    """Build exact advanced_text_search ``search_text`` for ``gene_symbol``.

    Includes the symbol, optional workflow ``aliases``, and gene-specific mapped
    terms only when present in :data:`GENE_SPECIFIC_EXACT_TERMS`.
    """
    symbol = gene_symbol.strip()
    terms: list[str] = [symbol] if symbol else []
    if aliases:
        terms.extend(str(a) for a in aliases)
    mapped = GENE_SPECIFIC_EXACT_TERMS.get(symbol.upper(), [])
    # Avoid duplicating the symbol when it is already first in the map.
    for token in mapped:
        if token.strip().upper() == symbol.upper() and symbol:
            continue
        terms.append(token)
    return _join_search_terms(terms)


def build_broader_search_text(
    gene_symbol: str,
    terms: list[str] | None = None,
) -> str | None:
    """Build broader pathway ``search_text``, or ``None`` if not applicable.

    Uses explicit ``terms`` when provided; otherwise gene-specific mapped terms
    only. Returns ``None`` when neither is available (e.g. ``HTT`` with no map).
    """
    symbol = gene_symbol.strip()
    if terms is not None:
        joined = _join_search_terms([str(t) for t in terms])
        return joined or None
    mapped = GENE_SPECIFIC_BROADER_TERMS.get(symbol.upper(), [])
    joined = _join_search_terms(list(mapped))
    return joined or None


def build_search_body(
    search_text: str,
    *,
    offset: int = DEFAULT_OFFSET,
    limit: int = DEFAULT_LIMIT,
    include_fields: list[str] | None = None,
    include_active_projects: bool = True,
    exclude_subprojects: bool = True,
    search_field: str = DEFAULT_SEARCH_FIELD,
    operator: str = DEFAULT_OPERATOR,
    sort_field: str = DEFAULT_SORT_FIELD,
    sort_order: str = DEFAULT_SORT_ORDER,
) -> dict[str, Any]:
    """Build a validated NIH RePORTER projects/search request body."""
    return {
        "criteria": {
            "advanced_text_search": {
                "operator": operator,
                "search_field": search_field,
                "search_text": search_text,
            },
            "include_active_projects": include_active_projects,
            "exclude_subprojects": exclude_subprojects,
        },
        "include_fields": list(include_fields or DEFAULT_INCLUDE_FIELDS),
        "offset": offset,
        "limit": limit,
        "sort_field": sort_field,
        "sort_order": sort_order,
    }


def summarize_project(project: dict[str, Any]) -> dict[str, Any]:
    """Extract key NIH RePORTER project fields (not evidence)."""
    org = project.get("organization") or {}
    if not isinstance(org, dict):
        org = {}
    agency = project.get("agency_ic_admin") or {}
    if not isinstance(agency, dict):
        agency = {}

    fundings = project.get("agency_ic_fundings") or []
    funding_names: list[str] = []
    if isinstance(fundings, list):
        for item in fundings:
            if isinstance(item, dict) and item.get("name"):
                funding_names.append(str(item["name"]))

    pis = project.get("principal_investigators") or []
    pi_names: list[str] = []
    if isinstance(pis, list):
        for pi in pis:
            if isinstance(pi, dict) and pi.get("full_name"):
                pi_names.append(str(pi["full_name"]))

    return {
        "project_title": project.get("project_title") or project.get("ProjectTitle"),
        "project_num": project.get("project_num") or project.get("ProjectNum"),
        "core_project_num": project.get("core_project_num")
        or project.get("CoreProjectNum"),
        "organization_name": org.get("org_name")
        or project.get("OrganizationName")
        or project.get("organization_name"),
        "budget_start": project.get("budget_start") or project.get("BudgetStartDate"),
        "budget_end": project.get("budget_end") or project.get("BudgetEndDate"),
        "agency_ic_admin": agency.get("name") or project.get("AgencyICAdmin"),
        "agency_ic_fundings": funding_names,
        "fiscal_year": project.get("fiscal_year") or project.get("FiscalYear"),
        "award_amount": project.get("award_amount") or project.get("AwardAmount"),
        "principal_investigators": pi_names,
        "project_detail_url": project.get("project_detail_url")
        or project.get("ProjectDetailUrl"),
        "abstract_text": project.get("abstract_text") or project.get("AbstractText"),
        "terms": project.get("terms") or project.get("Terms"),
    }


def _post_search(
    *,
    endpoint_name: str,
    gene_symbol: str,
    body: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """POST projects/search and return :class:`ToolResult`."""
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    # Provenance: store searchable criteria without dumping huge abstracts twice.
    request_params = {
        "endpoint": PROJECTS_SEARCH_URL,
        "operation": endpoint_name,
        "search_text": (
            ((body.get("criteria") or {}).get("advanced_text_search") or {}).get(
                "search_text"
            )
        ),
        "offset": body.get("offset"),
        "limit": body.get("limit"),
        "sort_field": body.get("sort_field"),
        "sort_order": body.get("sort_order"),
    }
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.post(
                PROJECTS_SEARCH_URL, json=body, headers=headers
            )
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=PROJECTS_SEARCH_URL,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=PROJECTS_SEARCH_URL,
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
            request_url=PROJECTS_SEARCH_URL,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=PROJECTS_SEARCH_URL,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=PROJECTS_SEARCH_URL,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def projects_search(
    search_text: str,
    *,
    gene_symbol: str = "",
    offset: int = DEFAULT_OFFSET,
    limit: int = DEFAULT_LIMIT,
    settings: Settings | None = None,
) -> ToolResult:
    """Run a NIH RePORTER projects search with the given ``search_text``."""
    cfg = settings or get_settings()
    body = build_search_body(search_text, offset=offset, limit=limit)
    return _post_search(
        endpoint_name="projects_search",
        gene_symbol=gene_symbol or search_text[:40],
        body=body,
        settings=cfg,
    )


def exact_search(
    gene_symbol: str,
    *,
    aliases: list[str] | None = None,
    search_text: str | None = None,
    offset: int = DEFAULT_OFFSET,
    limit: int = DEFAULT_LIMIT,
    settings: Settings | None = None,
) -> ToolResult:
    """Exact gene-oriented NIH RePORTER search (validated shape)."""
    cfg = settings or get_settings()
    text = search_text if search_text is not None else build_exact_search_text(
        gene_symbol, aliases=aliases
    )
    body = build_search_body(text, offset=offset, limit=limit)
    return _post_search(
        endpoint_name="exact_search",
        gene_symbol=gene_symbol,
        body=body,
        settings=cfg,
    )


def broader_search(
    gene_symbol: str,
    *,
    terms: list[str] | None = None,
    search_text: str | None = None,
    offset: int = DEFAULT_OFFSET,
    limit: int = DEFAULT_LIMIT,
    settings: Settings | None = None,
) -> ToolResult:
    """Broader pathway NIH RePORTER search when terms are available.

    Fails with ``invalid_request`` when no broader text can be built (no gene
    map and no explicit terms), so SREBF2 pathway terms are never applied to
    unrelated genes by default.
    """
    cfg = settings or get_settings()
    text = search_text
    if text is None:
        text = build_broader_search_text(gene_symbol, terms=terms)
    if not text:
        return _tool_result(
            endpoint_name="broader_search",
            gene_symbol=gene_symbol,
            request_url=PROJECTS_SEARCH_URL,
            request_params={"gene_symbol": gene_symbol},
            success=False,
            error_type="invalid_request",
            error_message=(
                f"No broader pathway search terms for {gene_symbol!r}; "
                "pass terms=... or add a gene-specific map entry"
            ),
        )
    body = build_search_body(text, offset=offset, limit=limit)
    return _post_search(
        endpoint_name="broader_search",
        gene_symbol=gene_symbol,
        body=body,
        settings=cfg,
    )


def fetch_grants(
    gene_symbol: str,
    *,
    aliases: list[str] | None = None,
    include_broader: bool = False,
    broader_terms: list[str] | None = None,
    limit: int = DEFAULT_LIMIT,
    settings: Settings | None = None,
) -> ToolResult:
    """Run exact search (required) and optional broader pathway search.

    On success, ``data`` includes exact (and optional broader) raw payloads plus
    project summaries. Broader search is off by default per API-map guidance.

    Never raises.
    """
    cfg = settings or get_settings()
    exact = exact_search(
        gene_symbol, aliases=aliases, limit=limit, settings=cfg
    )
    if not exact.success:
        return _tool_result(
            endpoint_name="fetch_grants",
            gene_symbol=gene_symbol,
            request_url=exact.request_url,
            request_params=exact.request_params,
            success=False,
            status_code=exact.status_code,
            data={"exact_search": exact.data},
            error_type=exact.error_type or "exact_search_failed",
            error_message=exact.error_message or "NIH RePORTER exact search failed",
        )

    exact_results = []
    if isinstance(exact.data, dict):
        raw_results = exact.data.get("results") or []
        if isinstance(raw_results, list):
            exact_results = [r for r in raw_results if isinstance(r, dict)]
    exact_summaries = [summarize_project(p) for p in exact_results]

    broader_payload: Any = None
    broader_summaries: list[dict[str, Any]] = []
    broader_search_text: str | None = None
    last_url = exact.request_url
    last_params = exact.request_params
    last_status = exact.status_code

    if include_broader:
        broader = broader_search(
            gene_symbol, terms=broader_terms, limit=limit, settings=cfg
        )
        last_url = broader.request_url
        last_params = broader.request_params
        last_status = broader.status_code
        if not broader.success:
            return _tool_result(
                endpoint_name="fetch_grants",
                gene_symbol=gene_symbol,
                request_url=broader.request_url,
                request_params=broader.request_params,
                success=False,
                status_code=broader.status_code,
                data={
                    "exact_search": exact.data,
                    "exact_summaries": exact_summaries,
                    "exact_count": len(exact_summaries),
                    "broader_search": broader.data,
                },
                error_type=broader.error_type or "broader_search_failed",
                error_message=broader.error_message
                or "NIH RePORTER broader search failed",
            )
        broader_payload = broader.data
        broader_search_text = (broader.request_params or {}).get("search_text")
        if isinstance(broader.data, dict):
            brows = broader.data.get("results") or []
            if isinstance(brows, list):
                broader_summaries = [
                    summarize_project(p) for p in brows if isinstance(p, dict)
                ]

    exact_text = build_exact_search_text(gene_symbol, aliases=aliases)
    return _tool_result(
        endpoint_name="fetch_grants",
        gene_symbol=gene_symbol,
        request_url=last_url,
        request_params={
            "gene_symbol": gene_symbol,
            "exact_search_text": exact_text,
            "include_broader": include_broader,
            "broader_search_text": broader_search_text,
            **last_params,
        },
        success=True,
        status_code=last_status,
        data={
            "gene_symbol": gene_symbol,
            "exact_search_text": exact_text,
            "exact_search": exact.data,
            "exact_summaries": exact_summaries,
            "exact_count": len(exact_summaries),
            "include_broader": include_broader,
            "broader_search_text": broader_search_text,
            "broader_search": broader_payload,
            "broader_summaries": broader_summaries,
            "broader_count": len(broader_summaries),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "PROJECTS_SEARCH_URL",
    "DEFAULT_LIMIT",
    "DEFAULT_INCLUDE_FIELDS",
    "GENE_SPECIFIC_EXACT_TERMS",
    "GENE_SPECIFIC_BROADER_TERMS",
    "build_exact_search_text",
    "build_broader_search_text",
    "build_search_body",
    "summarize_project",
    "projects_search",
    "exact_search",
    "broader_search",
    "fetch_grants",
]
