"""NCBI Gene client (E-utilities ESearch + ESummary).

Retrieves Entrez Gene candidates and summaries. Does **not** normalize into evidence
records — that belongs in ``normalize/gene_identity.py``.

Rules:
- Prefer exact official-symbol match; never blindly trust the first hit.
- Use ``NCBI_API_KEY`` when present (higher rate limits).
- Never raise: all failures return :class:`~gene_dossier.models.ToolResult`.

For SREBF2 / Homo sapiens, the expected Entrez Gene ID is ``6721``.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "NCBI Gene"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ORGANISM_HUMAN = "Homo sapiens"
ORGANISM_MOUSE = "Mus musculus"
ORGANISM_RAT = "Rattus norvegicus"


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


def _get_json(
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET an E-utilities endpoint and return JSON as :class:`ToolResult`."""
    query = _with_api_key(params, settings)
    url = f"{EUTILS_BASE}/{path}"
    # Full URL with query string for provenance (api_key redacted in params copy for logs).
    safe_params = {k: v for k, v in query.items() if k != "api_key"}
    if "api_key" in query:
        safe_params["api_key"] = "***"
    request_url = f"{url}?{urlencode(query)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query)
        try:
            payload = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:2000]}
        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=safe_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe_params,
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
            request_params=safe_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def esearch(
    gene_symbol: str,
    *,
    organism: str = ORGANISM_HUMAN,
    settings: Settings | None = None,
) -> ToolResult:
    """Search NCBI Gene for ``gene_symbol`` in ``organism``.

    Term format: ``{symbol}[Gene Name] AND {organism}[Organism]``.
    """
    cfg = settings or get_settings()
    term = f"{gene_symbol}[Gene Name] AND {organism}[Organism]"
    params = {
        "db": "gene",
        "term": term,
        "retmode": "json",
        "sort": "relevance",
    }
    return _get_json("esearch", gene_symbol, "esearch.fcgi", params, cfg)


def esummary(
    gene_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch NCBI Gene ESummary for one or more Entrez Gene IDs."""
    cfg = settings or get_settings()
    if isinstance(gene_ids, (list, tuple)):
        id_str = ",".join(str(g) for g in gene_ids)
    else:
        id_str = str(gene_ids)
    params = {
        "db": "gene",
        "id": id_str,
        "retmode": "json",
    }
    return _get_json("esummary", gene_symbol or id_str, "esummary.fcgi", params, cfg)


def extract_id_list(esearch_result: ToolResult) -> list[str]:
    """Return Entrez IDs from a successful ESearch :class:`ToolResult`."""
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return []
    result = esearch_result.data.get("esearchresult") or {}
    ids = result.get("idlist") or []
    return [str(i) for i in ids]


def _summary_uid_map(esummary_result: ToolResult) -> dict[str, dict[str, Any]]:
    """Map Entrez ID → summary dict from an ESummary payload."""
    if not esummary_result.success or not isinstance(esummary_result.data, dict):
        return {}
    result = esummary_result.data.get("result") or {}
    out: dict[str, dict[str, Any]] = {}
    for uid in result.get("uids") or []:
        entry = result.get(str(uid))
        if isinstance(entry, dict):
            out[str(uid)] = entry
    return out


def prefer_exact_symbol_match(
    summaries: dict[str, dict[str, Any]],
    gene_symbol: str,
) -> str | None:
    """Return the Entrez ID whose ``name`` equals ``gene_symbol`` (case-insensitive).

    If multiple exact matches exist, prefer the first. Returns ``None`` if none match.
    """
    target = gene_symbol.strip().upper()
    for uid, entry in summaries.items():
        name = str(entry.get("name") or "").strip().upper()
        if name == target:
            return uid
    return None


def lookup_gene(
    gene_symbol: str,
    *,
    organism: str = ORGANISM_HUMAN,
    settings: Settings | None = None,
) -> ToolResult:
    """ESearch + ESummary with exact official-symbol preference.

    On success, ``data`` is a dict::

        {
          "gene_symbol": ...,
          "organism": ...,
          "selected_gene_id": "6721" | null,
          "selection_method": "exact_symbol" | "first_candidate" | "none",
          "candidate_ids": [...],
          "esearch": <raw esearch json>,
          "esummary": <raw esummary json> | null,
          "selected_summary": <summary dict> | null,
        }

    Never raises.
    """
    cfg = settings or get_settings()
    search = esearch(gene_symbol, organism=organism, settings=cfg)
    if not search.success:
        return _tool_result(
            endpoint_name="lookup_gene",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=False,
            status_code=search.status_code,
            data={"esearch": search.data, "organism": organism},
            error_type=search.error_type or "esearch_failed",
            error_message=search.error_message or "ESearch failed",
        )

    candidate_ids = extract_id_list(search)
    if not candidate_ids:
        return _tool_result(
            endpoint_name="lookup_gene",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=False,
            status_code=search.status_code,
            data={
                "gene_symbol": gene_symbol,
                "organism": organism,
                "selected_gene_id": None,
                "selection_method": "none",
                "candidate_ids": [],
                "esearch": search.data,
                "esummary": None,
                "selected_summary": None,
            },
            error_type="no_results",
            error_message=f"No NCBI Gene hits for {gene_symbol} / {organism}",
        )

    # Cap summaries to avoid huge payloads; exact match is usually in the top hits.
    ids_for_summary = candidate_ids[:20]
    summary = esummary(ids_for_summary, gene_symbol=gene_symbol, settings=cfg)
    if not summary.success:
        return _tool_result(
            endpoint_name="lookup_gene",
            gene_symbol=gene_symbol,
            request_url=summary.request_url,
            request_params=summary.request_params,
            success=False,
            status_code=summary.status_code,
            data={
                "gene_symbol": gene_symbol,
                "organism": organism,
                "selected_gene_id": None,
                "selection_method": "none",
                "candidate_ids": candidate_ids,
                "esearch": search.data,
                "esummary": summary.data,
                "selected_summary": None,
            },
            error_type=summary.error_type or "esummary_failed",
            error_message=summary.error_message or "ESummary failed",
        )

    uid_map = _summary_uid_map(summary)
    exact_id = prefer_exact_symbol_match(uid_map, gene_symbol)
    if exact_id is not None:
        selected_id = exact_id
        method = "exact_symbol"
    else:
        selected_id = candidate_ids[0]
        method = "first_candidate"

    selected_summary = uid_map.get(selected_id)
    return _tool_result(
        endpoint_name="lookup_gene",
        gene_symbol=gene_symbol,
        request_url=summary.request_url,
        request_params={
            "organism": organism,
            "candidate_ids": candidate_ids,
            "selected_gene_id": selected_id,
            "selection_method": method,
        },
        success=True,
        status_code=summary.status_code,
        data={
            "gene_symbol": gene_symbol,
            "organism": organism,
            "selected_gene_id": selected_id,
            "selection_method": method,
            "candidate_ids": candidate_ids,
            "esearch": search.data,
            "esummary": summary.data,
            "selected_summary": selected_summary,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "EUTILS_BASE",
    "ORGANISM_HUMAN",
    "ORGANISM_MOUSE",
    "ORGANISM_RAT",
    "esearch",
    "esummary",
    "extract_id_list",
    "prefer_exact_symbol_match",
    "lookup_gene",
]
