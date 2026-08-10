"""PubTator3 client (entity autocomplete, relations, search, BioC-JSON export).

Does **not** normalize into evidence records — that belongs in Section 7a
orchestration / ``normalize/``.

Key endpoints::

    GET .../entity/autocomplete/?query={q}&concept={concept}&limit={n}
    GET .../relations?e1={entityId}&type={type}&e2={entity_type}
    GET .../search/?text=relations:{type}|{e1}|{e2}&page={page}
    GET .../publications/export/biocjson?pmids={ids}

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "PubTator3"
BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"

DEFAULT_AUTOCOMPLETE_LIMIT = 10
DEFAULT_SEARCH_PAGE = 1


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


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    url = f"{BASE}/{path.lstrip('/')}"
    query = {k: str(v) for k, v in params.items() if v is not None and str(v) != ""}
    request_url = f"{url}?{urlencode(query)}" if query else url
    headers = {"Accept": "application/json"}
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url, params=query, headers=headers)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params=query,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
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
            request_params=query,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def entity_autocomplete(
    query: str,
    *,
    concept: str = "gene",
    limit: int = DEFAULT_AUTOCOMPLETE_LIMIT,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Find PubTator entity IDs for free-text ``query``."""
    cfg = settings or get_settings()
    q = str(query).strip()
    params = {
        "query": q,
        "concept": concept,
        "limit": limit,
    }
    return _request_json(
        endpoint_name="entity_autocomplete",
        gene_symbol=gene_symbol or q,
        path="entity/autocomplete/",
        params=params,
        settings=cfg,
    )


def relations(
    e1: str,
    type: str,  # noqa: A002 — matches PubTator query param name
    e2: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Query related entities (counts / entity pairs; typically no PMIDs)."""
    cfg = settings or get_settings()
    params = {
        "e1": str(e1).strip(),
        "type": str(type).strip(),
        "e2": str(e2).strip(),
    }
    return _request_json(
        endpoint_name="relations",
        gene_symbol=gene_symbol or str(e1).strip(),
        path="relations",
        params=params,
        settings=cfg,
    )


def search_relations(
    relation_type: str,
    entity1: str,
    entity2: str,
    page: int = DEFAULT_SEARCH_PAGE,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Search literature for a typed relation between two entities.

    Calls ``GET /search/?text=relations:{type}|{e1}|{e2}&page={page}``.
    Caller controls entity order (commonly chemical then gene).
    """
    cfg = settings or get_settings()
    text = (
        f"relations:{str(relation_type).strip()}"
        f"|{str(entity1).strip()}"
        f"|{str(entity2).strip()}"
    )
    params = {
        "text": text,
        "page": int(page),
    }
    return _request_json(
        endpoint_name="search_relations",
        gene_symbol=gene_symbol or str(entity2).strip() or str(entity1).strip(),
        path="search/",
        params=params,
        settings=cfg,
    )


def fetch_publication_annotations(
    pmids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    full: bool = False,
    settings: Settings | None = None,
) -> ToolResult:
    """Export PubTator annotations as BioC-JSON for one or more PMIDs.

    Endpoint called::

        GET {BASE}/publications/export/biocjson?pmids={comma_ids}[&full=true]
    """
    cfg = settings or get_settings()
    if isinstance(pmids, (list, tuple)):
        id_str = ",".join(str(p).strip() for p in pmids if str(p).strip())
    else:
        id_str = str(pmids).strip()
    if not id_str:
        return _tool_result(
            endpoint_name="fetch_publication_annotations",
            gene_symbol=gene_symbol,
            request_url=f"{BASE}/publications/export/biocjson",
            request_params={"pmids": ""},
            success=False,
            error_type="invalid_request",
            error_message="PubTator3 biocjson export requires at least one PMID",
        )
    params: dict[str, Any] = {"pmids": id_str}
    if full:
        params["full"] = "true"
    return _request_json(
        endpoint_name="fetch_publication_annotations",
        gene_symbol=gene_symbol or id_str,
        path="publications/export/biocjson",
        params=params,
        settings=cfg,
    )


def biocjson_documents(payload: Any) -> list[dict[str, Any]]:
    """Return document dicts from a PubTator3 BioC-JSON payload."""
    if isinstance(payload, list):
        return [d for d in payload if isinstance(d, dict)]
    if not isinstance(payload, dict):
        return []
    for key in ("PubTator3", "documents", "document"):
        value = payload.get(key)
        if isinstance(value, list):
            return [d for d in value if isinstance(d, dict)]
        if isinstance(value, dict):
            return [value]
    if "passages" in payload or "annotations" in payload:
        return [payload]
    return []


def iter_biocjson_annotations(document: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten passage- and document-level annotations from one BioC document."""
    out: list[dict[str, Any]] = []
    top = document.get("annotations")
    if isinstance(top, list):
        out.extend(a for a in top if isinstance(a, dict))
    passages = document.get("passages") or []
    if isinstance(passages, list):
        for passage in passages:
            if not isinstance(passage, dict):
                continue
            anns = passage.get("annotations") or []
            if isinstance(anns, list):
                out.extend(a for a in anns if isinstance(a, dict))
    return out


def extract_chemical_entities(payload: Any) -> list[dict[str, Any]]:
    """Extract Chemical entity mentions from a BioC-JSON export payload.

    Each item includes pmid, text, identifier/accession when present, and
    location metadata. Annotation of a chemical is entity identification only.
    """
    chemicals: list[dict[str, Any]] = []
    for doc in biocjson_documents(payload):
        pmid = doc.get("pmid") or doc.get("id") or doc.get("_id")
        for ann in iter_biocjson_annotations(doc):
            infons = ann.get("infons") if isinstance(ann.get("infons"), dict) else {}
            ann_type = str(infons.get("type") or "").strip().lower()
            biotype = str(infons.get("biotype") or "").strip().lower()
            if ann_type != "chemical" and biotype != "chemical":
                continue
            chemicals.append(
                {
                    "pmid": str(pmid) if pmid is not None else None,
                    "text": ann.get("text"),
                    "identifier": infons.get("identifier"),
                    "normalized_id": infons.get("normalized_id"),
                    "accession": infons.get("accession"),
                    "name": infons.get("name"),
                    "database": infons.get("database"),
                    "valid": infons.get("valid"),
                    "locations": ann.get("locations"),
                    "infons": infons,
                }
            )
    return chemicals


def extract_search_pmids(search_payload: Any) -> list[str]:
    """Return PMID strings from a PubTator3 ``/search/`` response."""
    if not isinstance(search_payload, dict):
        return []
    results = search_payload.get("results") or []
    if not isinstance(results, list):
        return []
    out: list[str] = []
    seen: set[str] = set()
    for row in results:
        if not isinstance(row, dict):
            continue
        pmid = row.get("pmid")
        if pmid is None:
            pmid = row.get("_id")
        if pmid is None:
            continue
        key = str(pmid).strip()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


__all__ = [
    "SOURCE_NAME",
    "BASE",
    "DEFAULT_AUTOCOMPLETE_LIMIT",
    "DEFAULT_SEARCH_PAGE",
    "entity_autocomplete",
    "relations",
    "search_relations",
    "fetch_publication_annotations",
    "biocjson_documents",
    "iter_biocjson_annotations",
    "extract_chemical_entities",
    "extract_search_pmids",
]
