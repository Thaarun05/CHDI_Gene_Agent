"""ClinVar client (E-utilities ESearch + ESummary; optional EFetch).

Retrieves ClinVar variant IDs and summaries for a gene. Does **not** normalize into
evidence records — that belongs in ``normalize/variants.py``.

Key endpoints (validated)::

    GET .../esearch.fcgi?db=clinvar&term={gene}[gene] AND single_gene[prop]
        &retmode=json&retmax=500
    GET .../esummary.fcgi?db=clinvar&id={clinvar_ids}&retmode=json
    GET .../efetch.fcgi?db=clinvar&id={ids}&rettype=vcv&retmode=xml   (optional)

Rules:
- Prefer ESummary for dossier tables; EFetch VCV XML is optional.
- Chunk long comma-separated ID lists for ESummary/EFetch.
- Use ``NCBI_API_KEY`` when present (redacted in provenance fields).
- Never raise: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "ClinVar"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

DEFAULT_RETMAX = 500
# Keep request URLs under practical limits; API map notes long ID lists must be chunked.
ESUMMARY_CHUNK_SIZE = 100


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


def build_search_term(gene_symbol: str) -> str:
    """Build the validated ClinVar gene search term."""
    return f"{gene_symbol}[gene] AND single_gene[prop]"


def _normalize_ids(ids: str | list[str] | int) -> list[str]:
    """Normalize ClinVar ID input to a flat list of ID strings."""
    if isinstance(ids, (list, tuple)):
        return [str(i).strip() for i in ids if str(i).strip()]
    if isinstance(ids, int):
        return [str(ids)]
    text = str(ids).strip()
    if not text:
        return []
    return [part.strip() for part in text.split(",") if part.strip()]


def _chunked(items: list[str], size: int) -> list[list[str]]:
    """Split ``items`` into consecutive chunks of at most ``size``."""
    if size <= 0:
        return [items] if items else []
    return [items[i : i + size] for i in range(0, len(items), size)]


def _request(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
    expect_json: bool = True,
) -> ToolResult:
    """GET an E-utilities endpoint; return :class:`ToolResult` (never raises)."""
    query = _with_api_key(params, settings)
    url = f"{EUTILS_BASE}/{path}"
    # Provenance logging: redact api_key in both request_params and request_url.
    # The live HTTP call still uses the full ``query`` (with real key when present).
    safe = _safe_params(query)
    request_url = f"{url}?{urlencode(safe)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query)
        if expect_json:
            try:
                payload: Any = response.json()
            except ValueError:
                payload = {"raw_text": response.text[:4000]}
        else:
            payload = {
                "raw_text": response.text,
                "content_type": response.headers.get("content-type"),
            }

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


def esearch(
    gene_symbol: str,
    *,
    term: str | None = None,
    retmax: int = DEFAULT_RETMAX,
    settings: Settings | None = None,
) -> ToolResult:
    """Search ClinVar for variants linked to ``gene_symbol`` (single-gene prop)."""
    cfg = settings or get_settings()
    search_term = term if term is not None else build_search_term(gene_symbol)
    params = {
        "db": "clinvar",
        "term": search_term,
        "retmode": "json",
        "retmax": str(retmax),
    }
    return _request(
        endpoint_name="esearch",
        gene_symbol=gene_symbol,
        path="esearch.fcgi",
        params=params,
        settings=cfg,
        expect_json=True,
    )


def esummary(
    clinvar_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    chunk_size: int = ESUMMARY_CHUNK_SIZE,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch ClinVar ESummary for one or more IDs (chunked when needed).

    On success with multiple chunks, ``data`` merges ``result`` maps and records
    ``uids`` / ``chunk_count`` for provenance. Single-chunk responses keep the
    raw ESummary JSON shape.
    """
    cfg = settings or get_settings()
    ids = _normalize_ids(clinvar_ids)
    if not ids:
        return _tool_result(
            endpoint_name="esummary",
            gene_symbol=gene_symbol or "",
            request_url=f"{EUTILS_BASE}/esummary.fcgi",
            request_params={"db": "clinvar", "id": ""},
            success=False,
            error_type="invalid_request",
            error_message="ClinVar ESummary requires at least one ID",
        )

    chunks = _chunked(ids, chunk_size)
    if len(chunks) == 1:
        params = {
            "db": "clinvar",
            "id": ",".join(chunks[0]),
            "retmode": "json",
        }
        return _request(
            endpoint_name="esummary",
            gene_symbol=gene_symbol or ",".join(chunks[0]),
            path="esummary.fcgi",
            params=params,
            settings=cfg,
            expect_json=True,
        )

    merged_result: dict[str, Any] = {"uids": []}
    chunk_payloads: list[Any] = []
    last_url = f"{EUTILS_BASE}/esummary.fcgi"
    last_params: dict[str, Any] = {"db": "clinvar", "chunk_count": len(chunks)}
    last_status: int | None = None

    for chunk in chunks:
        params = {
            "db": "clinvar",
            "id": ",".join(chunk),
            "retmode": "json",
        }
        part = _request(
            endpoint_name="esummary",
            gene_symbol=gene_symbol or ",".join(chunk),
            path="esummary.fcgi",
            params=params,
            settings=cfg,
            expect_json=True,
        )
        last_url = part.request_url
        last_params = part.request_params
        last_status = part.status_code
        if not part.success:
            return _tool_result(
                endpoint_name="esummary",
                gene_symbol=gene_symbol or "",
                request_url=part.request_url,
                request_params={
                    **part.request_params,
                    "requested_ids": ids,
                    "chunk_count": len(chunks),
                },
                success=False,
                status_code=part.status_code,
                data={"partial_chunks": chunk_payloads, "failed_chunk": part.data},
                error_type=part.error_type or "esummary_failed",
                error_message=part.error_message or "ClinVar ESummary chunk failed",
            )
        chunk_payloads.append(part.data)
        if isinstance(part.data, dict):
            result = part.data.get("result") or {}
            for uid in result.get("uids") or []:
                uid_s = str(uid)
                merged_result["uids"].append(uid_s)
                entry = result.get(uid_s) or result.get(uid)
                if isinstance(entry, dict):
                    merged_result[uid_s] = entry

    return _tool_result(
        endpoint_name="esummary",
        gene_symbol=gene_symbol or "",
        request_url=last_url,
        request_params={
            **last_params,
            "requested_ids": ids,
            "chunk_count": len(chunks),
        },
        success=True,
        status_code=last_status,
        data={"result": merged_result, "chunk_count": len(chunks)},
    )


def efetch(
    clinvar_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    rettype: str = "vcv",
    retmode: str = "xml",
    chunk_size: int = ESUMMARY_CHUNK_SIZE,
    settings: Settings | None = None,
) -> ToolResult:
    """Optional ClinVar EFetch (default VCV XML). ESummary is preferred for tables."""
    cfg = settings or get_settings()
    ids = _normalize_ids(clinvar_ids)
    if not ids:
        return _tool_result(
            endpoint_name="efetch",
            gene_symbol=gene_symbol or "",
            request_url=f"{EUTILS_BASE}/efetch.fcgi",
            request_params={"db": "clinvar", "id": ""},
            success=False,
            error_type="invalid_request",
            error_message="ClinVar EFetch requires at least one ID",
        )

    # Single request for small lists; concatenate XML chunks when needed.
    chunks = _chunked(ids, chunk_size)
    texts: list[str] = []
    last_url = f"{EUTILS_BASE}/efetch.fcgi"
    last_params: dict[str, Any] = {}
    last_status: int | None = None

    for chunk in chunks:
        params = {
            "db": "clinvar",
            "id": ",".join(chunk),
            "rettype": rettype,
            "retmode": retmode,
        }
        part = _request(
            endpoint_name="efetch",
            gene_symbol=gene_symbol or ",".join(chunk),
            path="efetch.fcgi",
            params=params,
            settings=cfg,
            expect_json=False,
        )
        last_url = part.request_url
        last_params = part.request_params
        last_status = part.status_code
        if not part.success:
            return _tool_result(
                endpoint_name="efetch",
                gene_symbol=gene_symbol or "",
                request_url=part.request_url,
                request_params={
                    **part.request_params,
                    "requested_ids": ids,
                    "chunk_count": len(chunks),
                },
                success=False,
                status_code=part.status_code,
                data=part.data,
                error_type=part.error_type or "efetch_failed",
                error_message=part.error_message or "ClinVar EFetch chunk failed",
            )
        if isinstance(part.data, dict):
            texts.append(str(part.data.get("raw_text") or ""))

    return _tool_result(
        endpoint_name="efetch",
        gene_symbol=gene_symbol or "",
        request_url=last_url,
        request_params={
            **last_params,
            "requested_ids": ids,
            "chunk_count": len(chunks),
        },
        success=True,
        status_code=last_status,
        data={
            "raw_text": "\n".join(texts),
            "content_type": "application/xml",
            "chunk_count": len(chunks),
        },
    )


def extract_id_list(esearch_result: ToolResult) -> list[str]:
    """Return ClinVar IDs from a successful ESearch :class:`ToolResult`."""
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return []
    result = esearch_result.data.get("esearchresult") or {}
    ids = result.get("idlist") or []
    return [str(i) for i in ids]


def extract_count(esearch_result: ToolResult) -> int | None:
    """Return ESearch hit count when present."""
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return None
    result = esearch_result.data.get("esearchresult") or {}
    raw = result.get("count")
    if raw is None:
        return None
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _summary_uid_map(esummary_payload: Any) -> dict[str, dict[str, Any]]:
    """Map ClinVar UID → summary dict from an ESummary payload."""
    if not isinstance(esummary_payload, dict):
        return {}
    result = esummary_payload.get("result") or {}
    out: dict[str, dict[str, Any]] = {}
    for uid in result.get("uids") or []:
        entry = result.get(str(uid))
        if isinstance(entry, dict):
            out[str(uid)] = entry
    return out


def summarize_variant(uid: str, entry: dict[str, Any]) -> dict[str, Any]:
    """Extract key ClinVar ESummary fields for later normalization (not evidence)."""
    variation_set = entry.get("variation_set") or []
    first_var: dict[str, Any] = {}
    if isinstance(variation_set, list) and variation_set and isinstance(variation_set[0], dict):
        first_var = variation_set[0]

    genes = entry.get("genes") or []
    first_gene: dict[str, Any] = {}
    if isinstance(genes, list) and genes and isinstance(genes[0], dict):
        first_gene = genes[0]

    germline = entry.get("germline_classification") or {}
    if not isinstance(germline, dict):
        germline = {}

    trait_names: list[str] = []
    for trait in entry.get("trait_set") or []:
        if isinstance(trait, dict) and trait.get("trait_name"):
            trait_names.append(str(trait["trait_name"]))

    locs: list[dict[str, Any]] = []
    for loc in first_var.get("variation_loc") or []:
        if isinstance(loc, dict):
            locs.append(
                {
                    "assembly_name": loc.get("assembly_name"),
                    "chr": loc.get("chr"),
                    "start": loc.get("start"),
                    "stop": loc.get("stop"),
                }
            )

    return {
        "uid": str(uid),
        "accession": entry.get("accession"),
        "title": entry.get("title"),
        "obj_type": entry.get("obj_type"),
        "variation_name": first_var.get("variation_name"),
        "measure_id": first_var.get("measure_id"),
        "cdna_change": first_var.get("cdna_change"),
        "canonical_spdi": first_var.get("canonical_spdi"),
        "variation_locs": locs,
        "gene_symbol": first_gene.get("symbol"),
        "gene_id": first_gene.get("geneid"),
        "germline_classification": germline.get("description"),
        "review_status": germline.get("review_status"),
        "last_evaluated": germline.get("last_evaluated"),
        "trait_names": trait_names,
        "molecular_consequence_list": list(entry.get("molecular_consequence_list") or []),
        "protein_change": entry.get("protein_change"),
    }


def fetch_clinvar_variants(
    gene_symbol: str,
    *,
    retmax: int = DEFAULT_RETMAX,
    include_efetch: bool = False,
    settings: Settings | None = None,
) -> ToolResult:
    """ESearch + ESummary for ClinVar variants of ``gene_symbol``.

    On success, ``data`` includes::

        {
          "gene_symbol": ...,
          "search_term": ...,
          "clinvar_ids": [...],
          "count": int | None,
          "esearch": <raw>,
          "esummary": <raw> | null,
          "variant_summaries": [...],
          "efetch": <optional xml wrapper> | null,
        }

    Never raises. EFetch is off by default (ESummary is enough for dossier tables).
    """
    cfg = settings or get_settings()
    term = build_search_term(gene_symbol)
    search = esearch(gene_symbol, term=term, retmax=retmax, settings=cfg)
    if not search.success:
        return _tool_result(
            endpoint_name="fetch_clinvar_variants",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=False,
            status_code=search.status_code,
            data={"esearch": search.data, "search_term": term},
            error_type=search.error_type or "esearch_failed",
            error_message=search.error_message or "ClinVar ESearch failed",
        )

    clinvar_ids = extract_id_list(search)
    count = extract_count(search)
    if not clinvar_ids:
        return _tool_result(
            endpoint_name="fetch_clinvar_variants",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=search.request_params,
            success=True,
            status_code=search.status_code,
            data={
                "gene_symbol": gene_symbol,
                "search_term": term,
                "clinvar_ids": [],
                "count": count if count is not None else 0,
                "esearch": search.data,
                "esummary": None,
                "variant_summaries": [],
                "efetch": None,
            },
        )

    summary = esummary(clinvar_ids, gene_symbol=gene_symbol, settings=cfg)
    if not summary.success:
        return _tool_result(
            endpoint_name="fetch_clinvar_variants",
            gene_symbol=gene_symbol,
            request_url=summary.request_url,
            request_params=summary.request_params,
            success=False,
            status_code=summary.status_code,
            data={
                "gene_symbol": gene_symbol,
                "search_term": term,
                "clinvar_ids": clinvar_ids,
                "count": count,
                "esearch": search.data,
                "esummary": summary.data,
                "variant_summaries": [],
                "efetch": None,
            },
            error_type=summary.error_type or "esummary_failed",
            error_message=summary.error_message or "ClinVar ESummary failed",
        )

    uid_map = _summary_uid_map(summary.data)
    variant_summaries = [
        summarize_variant(uid, entry) for uid, entry in uid_map.items()
    ]

    efetch_payload: Any = None
    if include_efetch:
        fetched = efetch(clinvar_ids, gene_symbol=gene_symbol, settings=cfg)
        if not fetched.success:
            return _tool_result(
                endpoint_name="fetch_clinvar_variants",
                gene_symbol=gene_symbol,
                request_url=fetched.request_url,
                request_params=fetched.request_params,
                success=False,
                status_code=fetched.status_code,
                data={
                    "gene_symbol": gene_symbol,
                    "search_term": term,
                    "clinvar_ids": clinvar_ids,
                    "count": count,
                    "esearch": search.data,
                    "esummary": summary.data,
                    "variant_summaries": variant_summaries,
                    "efetch": fetched.data,
                },
                error_type=fetched.error_type or "efetch_failed",
                error_message=fetched.error_message or "ClinVar EFetch failed",
            )
        efetch_payload = fetched.data

    return _tool_result(
        endpoint_name="fetch_clinvar_variants",
        gene_symbol=gene_symbol,
        request_url=summary.request_url,
        request_params={
            "search_term": term,
            "clinvar_ids": clinvar_ids,
            "count": count,
            "include_efetch": include_efetch,
        },
        success=True,
        status_code=summary.status_code,
        data={
            "gene_symbol": gene_symbol,
            "search_term": term,
            "clinvar_ids": clinvar_ids,
            "count": count,
            "esearch": search.data,
            "esummary": summary.data,
            "variant_summaries": variant_summaries,
            "efetch": efetch_payload,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "EUTILS_BASE",
    "DEFAULT_RETMAX",
    "ESUMMARY_CHUNK_SIZE",
    "build_search_term",
    "esearch",
    "esummary",
    "efetch",
    "extract_id_list",
    "extract_count",
    "summarize_variant",
    "fetch_clinvar_variants",
]
