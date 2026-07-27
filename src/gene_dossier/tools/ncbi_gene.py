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

# Expected NCBI taxonomy IDs for selection safety checks.
ORGANISM_TAXID = {
    ORGANISM_HUMAN: 9606,
    ORGANISM_MOUSE: 10090,
    ORGANISM_RAT: 10116,
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
    # Provenance logging: redact api_key in both request_params and request_url.
    # The live HTTP call still uses the full ``query`` (with real key when present).
    safe_params = {k: v for k, v in query.items() if k != "api_key"}
    if "api_key" in query:
        safe_params["api_key"] = "***"
    request_url = f"{url}?{urlencode(safe_params)}"
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


def expected_taxid_for_organism(organism: str) -> int | None:
    """Return the NCBI taxonomy ID for a known organism name, if mapped."""
    return ORGANISM_TAXID.get(organism) or ORGANISM_TAXID.get(organism.strip())


def _entry_symbol_matches(entry: dict[str, Any], gene_symbol: str) -> bool:
    """True if ``name`` or ``nomenclaturesymbol`` matches ``gene_symbol`` (case-insensitive)."""
    target = gene_symbol.strip().upper()
    if not target:
        return False
    name = str(entry.get("name") or "").strip().upper()
    nomen = str(entry.get("nomenclaturesymbol") or "").strip().upper()
    return target in {name, nomen} and bool(target)


def _entry_taxid(entry: dict[str, Any]) -> int | None:
    """Extract organism taxid from an ESummary entry when present."""
    organism = entry.get("organism") or {}
    if not isinstance(organism, dict):
        return None
    raw = organism.get("taxid")
    if raw is None:
        return None
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _is_safe_gene_entry(
    entry: dict[str, Any],
    gene_symbol: str,
    *,
    expected_taxid: int | None = None,
) -> bool:
    """Return True if ``entry`` is a safe match for ``gene_symbol``.

    Prefers records where:
    - ``name`` or ``nomenclaturesymbol`` matches the requested symbol
    - ``organism.taxid`` matches ``expected_taxid`` when both are available
    - ``nomenclaturestatus`` is Official when available
    - ``currentid`` is empty when available (avoids retired/replaced records)
    - ``status`` is empty / not retired when available
    """
    if not _entry_symbol_matches(entry, gene_symbol):
        return False

    if expected_taxid is not None:
        taxid = _entry_taxid(entry)
        if taxid is not None and taxid != expected_taxid:
            return False

    nomen_status = entry.get("nomenclaturestatus")
    if nomen_status is not None and str(nomen_status).strip():
        if str(nomen_status).strip().lower() != "official":
            return False

    current_id = entry.get("currentid")
    if current_id is not None and str(current_id).strip():
        return False

    status = entry.get("status")
    if status is not None and str(status).strip():
        if str(status).strip().lower() in {"retired", "replaced", "secondary"}:
            return False

    return True


def prefer_safe_gene_match(
    summaries: dict[str, dict[str, Any]],
    gene_symbol: str,
    *,
    expected_taxid: int | None = None,
) -> str | None:
    """Return the best safe Entrez ID for ``gene_symbol``, or ``None`` if ambiguous.

    Does **not** fall back to the first candidate. See :func:`_is_safe_gene_entry`.
    When multiple equally ranked safe matches exist, returns ``None`` (ambiguous)
    instead of silently taking the first hit.
    """
    selected, _warnings = select_safe_gene_match(
        summaries, gene_symbol, expected_taxid=expected_taxid
    )
    return selected


def select_safe_gene_match(
    summaries: dict[str, dict[str, Any]],
    gene_symbol: str,
    *,
    expected_taxid: int | None = None,
) -> tuple[str | None, list[str]]:
    """Select a safe Entrez ID and return structured selection warnings.

    Returns ``(selected_uid_or_none, warnings)``. Warnings are emitted when no
    safe match exists or when multiple equally ranked safe matches remain.
    """
    safe: list[tuple[int, str]] = []
    for uid, entry in summaries.items():
        if not _is_safe_gene_entry(entry, gene_symbol, expected_taxid=expected_taxid):
            continue
        nomen_status = str(entry.get("nomenclaturestatus") or "").strip().lower()
        # Lower rank = better; Official first.
        rank = 0 if nomen_status == "official" else 1
        safe.append((rank, str(uid)))
    if not safe:
        return None, ["no_safe_ncbi_gene_match"]
    safe.sort(key=lambda item: (item[0], item[1]))
    best_rank = safe[0][0]
    top = [uid for rank, uid in safe if rank == best_rank]
    if len(top) > 1:
        return None, [f"ambiguous_safe_matches:{','.join(top)}"]
    return top[0], []


def prefer_exact_symbol_match(
    summaries: dict[str, dict[str, Any]],
    gene_symbol: str,
    *,
    expected_taxid: int | None = None,
) -> str | None:
    """Return a safe Entrez ID match (symbol + organism/status checks).

    Kept for callers; delegates to :func:`prefer_safe_gene_match`.
    """
    return prefer_safe_gene_match(
        summaries, gene_symbol, expected_taxid=expected_taxid
    )


def lookup_gene(
    gene_symbol: str,
    *,
    organism: str = ORGANISM_HUMAN,
    settings: Settings | None = None,
) -> ToolResult:
    """ESearch + ESummary with safe official-symbol preference.

    On success, ``data`` is a dict::

        {
          "gene_symbol": ...,
          "organism": ...,
          "selected_gene_id": "6721" | null,
          "selection_method": "exact_symbol" | "ambiguous" | "none",
          "candidate_ids": [...],
          "esearch": <raw esearch json>,
          "esummary": <raw esummary json> | null,
          "selected_summary": <summary dict> | null,
        }

    Never blindly trusts the first ESearch hit. If no safe match exists,
    ``selection_method`` is ``"ambiguous"`` and ``selected_gene_id`` is ``null``.

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
    expected_taxid = expected_taxid_for_organism(organism)
    selected_id, selection_warnings = select_safe_gene_match(
        uid_map, gene_symbol, expected_taxid=expected_taxid
    )
    if selected_id is not None:
        method = "exact_symbol"
        selected_summary = uid_map.get(selected_id)
    else:
        # Do not fall back to first_candidate — leave unresolved for human/downstream review.
        selected_id = None
        method = "ambiguous"
        selected_summary = None

    return _tool_result(
        endpoint_name="lookup_gene",
        gene_symbol=gene_symbol,
        request_url=summary.request_url,
        request_params={
            "organism": organism,
            "expected_taxid": expected_taxid,
            "candidate_ids": candidate_ids,
            "selected_gene_id": selected_id,
            "selection_method": method,
        },
        success=True,
        status_code=summary.status_code,
        data={
            "gene_symbol": gene_symbol,
            "organism": organism,
            "expected_taxid": expected_taxid,
            "selected_gene_id": selected_id,
            "selection_method": method,
            "selection_warnings": selection_warnings,
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
    "ORGANISM_TAXID",
    "esearch",
    "esummary",
    "extract_id_list",
    "expected_taxid_for_organism",
    "prefer_safe_gene_match",
    "select_safe_gene_match",
    "prefer_exact_symbol_match",
    "lookup_gene",
]
