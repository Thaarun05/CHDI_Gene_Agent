"""GEO client (E-utilities: GEO Profiles → GDS).

Finds GEO Profiles for a gene, links them to GEO DataSets (GDS), and fetches
dataset summaries. Does **not** normalize into evidence records — that belongs
in ``normalize/perturbation.py``.

Validated chain (important)::

    geoprofiles ESearch → ELink to gds → GDS ESummary

``esummary.fcgi?db=geoprofiles`` failed in validation — do not use it.

Key endpoints::

    GET .../esearch.fcgi?db=geoprofiles&term=...
    GET .../elink.fcgi?dbfrom=geoprofiles&db=gds&id=...
    GET .../esummary.fcgi?db=gds&id=...

Uses ``NCBI_API_KEY`` when present (redacted in provenance fields).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any, Literal
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "GEO"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

DEFAULT_RETMAX = 50
ORGANISM_MOUSE = "Mus musculus"
ORGANISM_HUMAN = "Homo sapiens"

# CHDI/HD-oriented context filters from the validated GEO search plan.
BRAIN_CONTEXT = (
    "(brain OR neuron OR hippocampus OR cortex OR cerebellum OR striatum)"
)
PERTURBATION_CONTEXT = (
    "(stress OR fluoxetine OR antidepressant OR paraquat OR mutant OR "
    "knockout OR treatment OR exposed OR disease)"
)

SearchContext = Literal["broad", "brain", "perturbation"]


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


def _with_api_key(params: dict[str, Any], settings: Settings) -> dict[str, Any]:
    """Return a copy of ``params`` with ``api_key`` added when configured."""
    out = dict(params)
    if settings.has_key("NCBI_API_KEY"):
        out["api_key"] = settings.ncbi_api_key
    return out


def _safe_params(params: dict[str, Any]) -> dict[str, Any]:
    """Redact API key for provenance logging."""
    out = {k: v for k, v in params.items() if k != "api_key"}
    if "api_key" in params:
        out["api_key"] = "***"
    return out


def build_geoprofiles_term(
    gene_symbol: str,
    *,
    organism: str = ORGANISM_MOUSE,
    context: SearchContext = "broad",
) -> str:
    """Build a GEO Profiles ESearch term.

    Contexts match the validated plan:
    - ``broad`` — gene + organism
    - ``brain`` — plus brain/neuron tissue terms
    - ``perturbation`` — brain terms plus stress/treatment/knockout terms
    """
    symbol = gene_symbol.strip()
    base = f'{symbol} AND "{organism}"[Organism]'
    if context == "brain":
        return f"{base} AND {BRAIN_CONTEXT}"
    if context == "perturbation":
        return f"{base} AND {BRAIN_CONTEXT} AND {PERTURBATION_CONTEXT}"
    return base


def _request(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET an E-utilities endpoint; return :class:`ToolResult` (never raises)."""
    query = _with_api_key(params, settings)
    url = f"{EUTILS_BASE}/{path}"
    safe = _safe_params(query)
    request_url = f"{url}?{urlencode(safe)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
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
            request_params=safe,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def esearch_geoprofiles(
    gene_symbol: str,
    *,
    organism: str = ORGANISM_MOUSE,
    context: SearchContext = "broad",
    term: str | None = None,
    retmax: int = DEFAULT_RETMAX,
    settings: Settings | None = None,
) -> ToolResult:
    """ESearch GEO Profiles for ``gene_symbol``."""
    cfg = settings or get_settings()
    search_term = term if term is not None else build_geoprofiles_term(
        gene_symbol, organism=organism, context=context
    )
    params = {
        "db": "geoprofiles",
        "term": search_term,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": "relevance",
    }
    return _request(
        endpoint_name="esearch_geoprofiles",
        gene_symbol=gene_symbol,
        path="esearch.fcgi",
        params=params,
        settings=cfg,
    )


def elink_profiles_to_gds(
    profile_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """ELink GEO Profiles IDs to GEO DataSets (``db=gds``)."""
    cfg = settings or get_settings()
    if isinstance(profile_ids, (list, tuple)):
        id_str = ",".join(str(i).strip() for i in profile_ids if str(i).strip())
    else:
        id_str = str(profile_ids).strip()
    if not id_str:
        return _tool_result(
            endpoint_name="elink_profiles_to_gds",
            gene_symbol=gene_symbol,
            request_url=f"{EUTILS_BASE}/elink.fcgi",
            request_params={"dbfrom": "geoprofiles", "db": "gds", "id": ""},
            success=False,
            error_type="invalid_request",
            error_message="ELink requires at least one GEO Profiles ID",
        )
    params = {
        "dbfrom": "geoprofiles",
        "db": "gds",
        "id": id_str,
        "retmode": "json",
    }
    return _request(
        endpoint_name="elink_profiles_to_gds",
        gene_symbol=gene_symbol or id_str,
        path="elink.fcgi",
        params=params,
        settings=cfg,
    )


def esummary_gds(
    gds_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """ESummary for GEO DataSets (``db=gds``). Do not use ``db=geoprofiles``."""
    cfg = settings or get_settings()
    if isinstance(gds_ids, (list, tuple)):
        id_str = ",".join(str(i).strip() for i in gds_ids if str(i).strip())
    else:
        id_str = str(gds_ids).strip()
    if not id_str:
        return _tool_result(
            endpoint_name="esummary_gds",
            gene_symbol=gene_symbol,
            request_url=f"{EUTILS_BASE}/esummary.fcgi",
            request_params={"db": "gds", "id": ""},
            success=False,
            error_type="invalid_request",
            error_message="GDS ESummary requires at least one ID",
        )
    params = {
        "db": "gds",
        "id": id_str,
        "retmode": "json",
    }
    return _request(
        endpoint_name="esummary_gds",
        gene_symbol=gene_symbol or id_str,
        path="esummary.fcgi",
        params=params,
        settings=cfg,
    )


def extract_id_list(esearch_result: ToolResult) -> list[str]:
    """Return GEO Profiles IDs from a successful ESearch :class:`ToolResult`."""
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return []
    result = esearch_result.data.get("esearchresult") or {}
    ids = result.get("idlist") or []
    return [str(i) for i in ids]


def extract_gds_uids(elink_payload: Any) -> list[str]:
    """Extract unique GDS UIDs from an ELink JSON payload.

    Only reads ``linksetdbs[*].links``. Root-level ``linkset["ids"]`` are source
    GEO Profiles IDs and must not be treated as GDS UIDs.
    """
    if not isinstance(elink_payload, dict):
        return []
    linksets = elink_payload.get("linksets") or []
    if not isinstance(linksets, list):
        return []
    out: list[str] = []
    seen: set[str] = set()

    def _add(uid: Any) -> None:
        s = str(uid).strip()
        if s and s not in seen:
            seen.add(s)
            out.append(s)

    for linkset in linksets:
        if not isinstance(linkset, dict):
            continue
        for db in linkset.get("linksetdbs") or []:
            if not isinstance(db, dict):
                continue
            links = db.get("links") or []
            if not isinstance(links, list) or not links:
                continue
            dbto = str(db.get("dbto") or "").lower()
            # Prefer dbto == "gds"; allow missing dbto only when there is no
            # clear mismatch to another database.
            if dbto and dbto != "gds":
                continue
            for uid in links:
                _add(uid)
    return out


def summarize_gds(uid: str, entry: dict[str, Any]) -> dict[str, Any]:
    """Extract key GDS ESummary fields for later normalization (not evidence)."""
    samples = entry.get("samples") or []
    sample_rows: list[dict[str, Any]] = []
    if isinstance(samples, list):
        for sample in samples:
            if isinstance(sample, dict):
                sample_rows.append(
                    {
                        "accession": sample.get("accession"),
                        "title": sample.get("title"),
                    }
                )
    return {
        "gds_uid": str(uid),
        "accession": entry.get("accession"),
        "title": entry.get("title"),
        "summary": entry.get("summary"),
        "gpl": entry.get("gpl"),
        "gse": entry.get("gse"),
        "taxon": entry.get("taxon"),
        "gdstype": entry.get("gdstype"),
        "valtype": entry.get("valtype"),
        "ssinfo": entry.get("ssinfo"),
        "subsetinfo": entry.get("subsetinfo"),
        "n_samples": entry.get("n_samples"),
        "samples": sample_rows,
        "pubmedids": list(entry.get("pubmedids") or []),
        "ftplink": entry.get("ftplink"),
    }


def _gds_uid_map(esummary_payload: Any) -> dict[str, dict[str, Any]]:
    """Map GDS UID → summary dict from an ESummary payload."""
    if not isinstance(esummary_payload, dict):
        return {}
    result = esummary_payload.get("result") or {}
    out: dict[str, dict[str, Any]] = {}
    for uid in result.get("uids") or []:
        entry = result.get(str(uid))
        if isinstance(entry, dict):
            out[str(uid)] = entry
    return out


def fetch_perturbations(
    gene_symbol: str,
    *,
    organism: str = ORGANISM_MOUSE,
    context: SearchContext = "perturbation",
    retmax: int = DEFAULT_RETMAX,
    settings: Settings | None = None,
) -> ToolResult:
    """GEO Profiles ESearch → ELink to GDS → GDS ESummary.

    Default ``context="perturbation"`` matches the validated HD/brain
    perturbation search. Use ``broad`` or ``brain`` to widen results.

    On success, ``data`` includes profile IDs, GDS UIDs, raw payloads, and
    GDS summaries. Never raises.
    """
    cfg = settings or get_settings()
    term = build_geoprofiles_term(
        gene_symbol, organism=organism, context=context
    )
    search = esearch_geoprofiles(
        gene_symbol,
        organism=organism,
        context=context,
        term=term,
        retmax=retmax,
        settings=cfg,
    )
    if not search.success:
        return _tool_result(
            endpoint_name="fetch_perturbations",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=False,
            status_code=search.status_code,
            data={"esearch": search.data, "search_term": term},
            error_type=search.error_type or "esearch_failed",
            error_message=search.error_message or "GEO Profiles ESearch failed",
        )

    profile_ids = extract_id_list(search)
    if not profile_ids:
        return _tool_result(
            endpoint_name="fetch_perturbations",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=True,
            status_code=search.status_code,
            data={
                "gene_symbol": gene_symbol,
                "organism": organism,
                "context": context,
                "search_term": term,
                "profile_ids": [],
                "gds_uids": [],
                "esearch": search.data,
                "elink": None,
                "esummary": None,
                "gds_summaries": [],
            },
        )

    link = elink_profiles_to_gds(
        profile_ids, gene_symbol=gene_symbol, settings=cfg
    )
    if not link.success:
        return _tool_result(
            endpoint_name="fetch_perturbations",
            gene_symbol=gene_symbol,
            request_url=link.request_url,
            request_params=link.request_params,
            success=False,
            status_code=link.status_code,
            data={
                "gene_symbol": gene_symbol,
                "search_term": term,
                "profile_ids": profile_ids,
                "esearch": search.data,
                "elink": link.data,
            },
            error_type=link.error_type or "elink_failed",
            error_message=link.error_message or "GEO Profiles→GDS ELink failed",
        )

    gds_uids = extract_gds_uids(link.data)
    if not gds_uids:
        return _tool_result(
            endpoint_name="fetch_perturbations",
            gene_symbol=gene_symbol,
            request_url=link.request_url,
            request_params=link.request_params,
            success=True,
            status_code=link.status_code,
            data={
                "gene_symbol": gene_symbol,
                "organism": organism,
                "context": context,
                "search_term": term,
                "profile_ids": profile_ids,
                "gds_uids": [],
                "esearch": search.data,
                "elink": link.data,
                "esummary": None,
                "gds_summaries": [],
            },
        )

    summary = esummary_gds(gds_uids, gene_symbol=gene_symbol, settings=cfg)
    if not summary.success:
        return _tool_result(
            endpoint_name="fetch_perturbations",
            gene_symbol=gene_symbol,
            request_url=summary.request_url,
            request_params=summary.request_params,
            success=False,
            status_code=summary.status_code,
            data={
                "gene_symbol": gene_symbol,
                "search_term": term,
                "profile_ids": profile_ids,
                "gds_uids": gds_uids,
                "esearch": search.data,
                "elink": link.data,
                "esummary": summary.data,
            },
            error_type=summary.error_type or "esummary_failed",
            error_message=summary.error_message or "GEO GDS ESummary failed",
        )

    uid_map = _gds_uid_map(summary.data)
    gds_summaries = [
        summarize_gds(uid, entry) for uid, entry in uid_map.items()
    ]
    return _tool_result(
        endpoint_name="fetch_perturbations",
        gene_symbol=gene_symbol,
        request_url=summary.request_url,
        request_params={
            "organism": organism,
            "context": context,
            "search_term": term,
            "profile_ids": profile_ids,
            "gds_uids": gds_uids,
        },
        success=True,
        status_code=summary.status_code,
        data={
            "gene_symbol": gene_symbol,
            "organism": organism,
            "context": context,
            "search_term": term,
            "profile_ids": profile_ids,
            "gds_uids": gds_uids,
            "esearch": search.data,
            "elink": link.data,
            "esummary": summary.data,
            "gds_summaries": gds_summaries,
            "gds_count": len(gds_summaries),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "EUTILS_BASE",
    "DEFAULT_RETMAX",
    "ORGANISM_MOUSE",
    "ORGANISM_HUMAN",
    "BRAIN_CONTEXT",
    "PERTURBATION_CONTEXT",
    "build_geoprofiles_term",
    "esearch_geoprofiles",
    "elink_profiles_to_gds",
    "esummary_gds",
    "extract_id_list",
    "extract_gds_uids",
    "summarize_gds",
    "fetch_perturbations",
]
