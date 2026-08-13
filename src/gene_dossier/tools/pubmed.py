"""PubMed client (E-utilities ESearch + ESummary + EFetch).

Searches literature for a gene in the context of Huntington's disease. Does **not**
normalize into evidence records — that belongs in ``normalize/literature.py``.

Default controlled HD-literature search term (per platform spec)::

    ("{gene_symbol}"[Title/Abstract] OR <validated aliases...>) AND
    ("Huntington Disease"[MeSH Terms] OR <HD title/abstract context...>)

Rules:
- Do not assume every hit strongly supports a gene–HD link.
- Use ``NCBI_API_KEY`` when present.
- Never raise: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "PubMed"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

DEFAULT_RETMAX = 50
HD_MESH = '"Huntington Disease"[MeSH Terms]'
MAX_GENE_ALIASES = 8
_UNSAFE_ALIAS_STOPLIST = {
    "AND",
    "CAG",
    "DNA",
    "GENE",
    "HD",
    "HUMAN",
    "MOUSE",
    "NOT",
    "NUC",
    "OR",
    "PROTEIN",
    "RAT",
    "RNA",
    "SNP",
}


@dataclass(frozen=True)
class PubMedGeneTerms:
    canonical_symbol: str
    aliases_used: tuple[str, ...] = ()
    full_name_used: str | None = None


def _clean_pubmed_term(term: Any, *, allow_spaces: bool = False) -> str | None:
    text = str(term or "").strip()
    if not text:
        return None
    text = re.sub(r"\s+", " ", text)
    if len(text) > 90:
        return None
    pattern = r"^[A-Za-z0-9][A-Za-z0-9_.\-/ ]*[A-Za-z0-9]$" if allow_spaces else r"^[A-Za-z0-9][A-Za-z0-9_.\-]*[A-Za-z0-9]$"
    if not re.match(pattern, text):
        return None
    if not allow_spaces:
        if len(text) < 3 or text.upper() in _UNSAFE_ALIAS_STOPLIST or text.isdigit():
            return None
    elif len(text) < 6:
        return None
    return text


def _quote_title_abstract_term(term: str) -> str:
    escaped = term.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"[Title/Abstract]'


def _dedupe_terms(terms: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for term in terms:
        key = term.casefold()
        if key in seen:
            continue
        seen.add(key)
        out.append(term)
    return out


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


def build_gene_title_abstract_clause(
    canonical_symbol: str,
    aliases: list[str] | tuple[str, ...] | None = None,
    *,
    full_name: str | None = None,
    max_aliases: int = MAX_GENE_ALIASES,
) -> tuple[str, PubMedGeneTerms]:
    """Build a bounded PubMed Title/Abstract clause from validated identity terms."""
    canonical = _clean_pubmed_term(canonical_symbol) or str(canonical_symbol).strip()
    alias_terms = [
        clean
        for alias in list(aliases or [])
        if (clean := _clean_pubmed_term(alias)) is not None
    ]
    alias_terms = _dedupe_terms([term for term in alias_terms if term.casefold() != canonical.casefold()])
    alias_terms = alias_terms[:max_aliases]
    full = _clean_pubmed_term(full_name, allow_spaces=True) if full_name else None
    all_terms = _dedupe_terms([canonical, *alias_terms, *([full] if full else [])])
    rendered = [_quote_title_abstract_term(term) for term in all_terms]
    clause = rendered[0] if len(rendered) == 1 else f"({' OR '.join(rendered)})"
    return clause, PubMedGeneTerms(
        canonical_symbol=canonical,
        aliases_used=tuple(alias_terms),
        full_name_used=full,
    )


def build_hd_context_clause() -> str:
    """Return the deterministic HD context clause for controlled HD-literature search."""
    return (
        f"({HD_MESH} OR "
        '"Huntington disease"[Title/Abstract] OR '
        '"Huntington\'s disease"[Title/Abstract] OR '
        '"huntingtin"[Title/Abstract])'
    )


def build_hd_search_term(
    gene_symbol: str,
    *,
    aliases: list[str] | tuple[str, ...] | None = None,
    full_name: str | None = None,
) -> tuple[str, PubMedGeneTerms]:
    """Build the controlled gene-alias × Huntington disease PubMed search term."""
    gene_clause, gene_terms = build_gene_title_abstract_clause(
        gene_symbol,
        aliases,
        full_name=full_name,
    )
    return f"{gene_clause} AND {build_hd_context_clause()}", gene_terms


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
            payload = {"raw_text": response.text, "content_type": response.headers.get("content-type")}

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
    """Search PubMed. Default term is gene Title/Abstract + HD MeSH."""
    cfg = settings or get_settings()
    search_term = term if term is not None else build_hd_search_term(gene_symbol)[0]
    params = {
        "db": "pubmed",
        "term": search_term,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": "relevance",
    }
    return _request(
        endpoint_name="esearch",
        gene_symbol=gene_symbol,
        path="esearch.fcgi",
        params=params,
        settings=cfg,
        expect_json=True,
    )


def esearch_custom(
    term: str,
    retmax: int = DEFAULT_RETMAX,
    *,
    gene_symbol: str = "",
    sort: str = "relevance",
    retstart: int = 0,
    settings: Settings | None = None,
) -> ToolResult:
    """ESearch with an arbitrary term (no HD-MeSH default)."""
    cfg = settings or get_settings()
    search_term = str(term).strip()
    params = {
        "db": "pubmed",
        "term": search_term,
        "retmode": "json",
        "retmax": str(retmax),
        "retstart": str(retstart),
        "sort": sort,
    }
    return _request(
        endpoint_name="esearch_custom",
        gene_symbol=gene_symbol or search_term[:80],
        path="esearch.fcgi",
        params=params,
        settings=cfg,
        expect_json=True,
    )


def esummary(
    pubmed_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PubMed ESummary metadata for one or more PMIDs."""
    cfg = settings or get_settings()
    if isinstance(pubmed_ids, (list, tuple)):
        id_str = ",".join(str(p) for p in pubmed_ids)
    else:
        id_str = str(pubmed_ids)
    params = {
        "db": "pubmed",
        "id": id_str,
        "retmode": "json",
    }
    return _request(
        endpoint_name="esummary",
        gene_symbol=gene_symbol or id_str,
        path="esummary.fcgi",
        params=params,
        settings=cfg,
        expect_json=True,
    )


def efetch(
    pubmed_ids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    rettype: str = "abstract",
    retmode: str = "xml",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PubMed records (default: abstract XML) for one or more PMIDs."""
    cfg = settings or get_settings()
    if isinstance(pubmed_ids, (list, tuple)):
        id_str = ",".join(str(p) for p in pubmed_ids)
    else:
        id_str = str(pubmed_ids)
    params: dict[str, Any] = {
        "db": "pubmed",
        "id": id_str,
        "retmode": retmode,
    }
    if rettype:
        params["rettype"] = rettype
    return _request(
        endpoint_name="efetch",
        gene_symbol=gene_symbol or id_str,
        path="efetch.fcgi",
        params=params,
        settings=cfg,
        expect_json=False,
    )


def efetch_abstracts(
    pmids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PubMed abstract XML for one or more PMIDs (Section 7a helper)."""
    cfg = settings or get_settings()
    if isinstance(pmids, (list, tuple)):
        id_list = [str(p).strip() for p in pmids if str(p).strip()]
        id_str = ",".join(id_list)
    else:
        id_str = str(pmids).strip()
        id_list = [id_str] if id_str else []
    if not id_str:
        return _tool_result(
            endpoint_name="efetch_abstracts",
            gene_symbol=gene_symbol,
            request_url=f"{EUTILS_BASE}/efetch.fcgi",
            request_params={"db": "pubmed", "id": "", "rettype": "abstract", "retmode": "xml"},
            success=False,
            error_type="invalid_request",
            error_message="efetch_abstracts requires at least one PMID",
        )
    result = efetch(
        id_list,
        gene_symbol=gene_symbol,
        rettype="abstract",
        retmode="xml",
        settings=cfg,
    )
    if result.endpoint_name == "efetch":
        # Preserve transport fields but label the endpoint for Section 7a provenance.
        return _tool_result(
            endpoint_name="efetch_abstracts",
            gene_symbol=result.gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=result.success,
            status_code=result.status_code,
            data=result.data,
            error_type=result.error_type,
            error_message=result.error_message,
        )
    return result


def extract_id_list(esearch_result: ToolResult) -> list[str]:
    """Return PMIDs from a successful ESearch :class:`ToolResult`."""
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return []
    result = esearch_result.data.get("esearchresult") or {}
    ids = result.get("idlist") or []
    return [str(i) for i in ids]


def search_hd_literature(
    gene_symbol: str,
    *,
    aliases: list[str] | tuple[str, ...] | None = None,
    full_name: str | None = None,
    retmax: int = DEFAULT_RETMAX,
    fetch_abstracts: bool = True,
    settings: Settings | None = None,
) -> ToolResult:
    """Run HD-focused ESearch, then ESummary (+ optional EFetch abstracts).

    On success, ``data`` is::

        {
          "gene_symbol": ...,
          "search_term": ...,
          "pmids": [...],
          "count": int | None,
          "esearch": <raw>,
          "esummary": <raw> | null,
          "efetch": <raw xml wrapper> | null,
          "caveat": "Do not assume every hit strongly supports the gene–HD link.",
        }

    Never raises.
    """
    cfg = settings or get_settings()
    term, gene_terms = build_hd_search_term(
        gene_symbol,
        aliases=aliases,
        full_name=full_name,
    )
    search = esearch(gene_symbol, term=term, retmax=retmax, settings=cfg)
    provenance_params = {
        "db": "pubmed",
        "search_term": term,
        "final_search_term": term,
        "canonical_symbol": gene_terms.canonical_symbol,
        "aliases_used": list(gene_terms.aliases_used),
        "full_name_used": gene_terms.full_name_used,
        "retmax": str(retmax),
        "sort": "relevance",
        "fetch_abstracts": fetch_abstracts,
    }
    if not search.success:
        return _tool_result(
            endpoint_name="search_hd_literature",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=provenance_params,
            success=False,
            status_code=search.status_code,
            data={
                "esearch": search.data,
                "search_term": term,
                "final_search_term": term,
                "canonical_symbol": gene_terms.canonical_symbol,
                "aliases_used": list(gene_terms.aliases_used),
                "full_name_used": gene_terms.full_name_used,
            },
            error_type=search.error_type or "esearch_failed",
            error_message=search.error_message or "PubMed ESearch failed",
        )

    pmids = extract_id_list(search)
    count_raw = None
    if isinstance(search.data, dict):
        count_raw = (search.data.get("esearchresult") or {}).get("count")

    if not pmids:
        return _tool_result(
            endpoint_name="search_hd_literature",
            gene_symbol=gene_symbol,
            request_url=search.request_url,
            request_params=provenance_params,
            success=True,
            status_code=search.status_code,
            data={
                "gene_symbol": gene_symbol,
                "search_term": term,
                "final_search_term": term,
                "canonical_symbol": gene_terms.canonical_symbol,
                "aliases_used": list(gene_terms.aliases_used),
                "full_name_used": gene_terms.full_name_used,
                "pmids": [],
                "count": int(count_raw) if count_raw is not None else 0,
                "esearch": search.data,
                "esummary": None,
                "efetch": None,
                "caveat": "Do not assume every hit strongly supports the gene–HD link.",
            },
        )

    ids_for_fetch = pmids[:retmax]
    summary = esummary(ids_for_fetch, gene_symbol=gene_symbol, settings=cfg)
    fetch: ToolResult | None = None
    if fetch_abstracts:
        fetch = efetch(ids_for_fetch, gene_symbol=gene_symbol, settings=cfg)

    # Success if we got summaries; efetch failure is recorded but does not void the search.
    ok = summary.success
    error_type = None if ok else (summary.error_type or "esummary_failed")
    error_message = None if ok else (summary.error_message or "PubMed ESummary failed")

    return _tool_result(
        endpoint_name="search_hd_literature",
        gene_symbol=gene_symbol,
        request_url=search.request_url,
        request_params={
            **provenance_params,
            "pmids": ids_for_fetch,
        },
        success=ok,
        status_code=summary.status_code or search.status_code,
        data={
            "gene_symbol": gene_symbol,
            "search_term": term,
            "final_search_term": term,
            "canonical_symbol": gene_terms.canonical_symbol,
            "aliases_used": list(gene_terms.aliases_used),
            "full_name_used": gene_terms.full_name_used,
            "pmids": ids_for_fetch,
            "count": int(count_raw) if count_raw is not None else len(pmids),
            "esearch": search.data,
            "esummary": summary.data,
            "efetch": fetch.data if fetch is not None else None,
            "efetch_success": fetch.success if fetch is not None else None,
            "caveat": "Do not assume every hit strongly supports the gene–HD link.",
        },
        error_type=error_type,
        error_message=error_message,
    )



def efetch_pubmed_xml(
    pmids: str | list[str] | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch full PubMed XML (includes ChemicalList when present)."""
    cfg = settings or get_settings()
    if isinstance(pmids, (list, tuple)):
        id_list = [str(p).strip() for p in pmids if str(p).strip()]
    else:
        id_list = [str(pmids).strip()] if str(pmids).strip() else []
    if not id_list:
        return _tool_result(
            endpoint_name="efetch_pubmed_xml",
            gene_symbol=gene_symbol,
            request_url=f"{EUTILS_BASE}/efetch.fcgi",
            request_params={"db": "pubmed", "id": "", "retmode": "xml"},
            success=False,
            error_type="invalid_request",
            error_message="efetch_pubmed_xml requires at least one PMID",
        )
    result = efetch(
        id_list,
        gene_symbol=gene_symbol,
        rettype="",
        retmode="xml",
        settings=cfg,
    )
    return _tool_result(
        endpoint_name="efetch_pubmed_xml",
        gene_symbol=result.gene_symbol,
        request_url=result.request_url,
        request_params=result.request_params,
        success=result.success,
        status_code=result.status_code,
        data=result.data,
        error_type=result.error_type,
        error_message=result.error_message,
    )



def extract_medline_chemicals(raw_xml: str) -> dict[str, list[dict[str, str]]]:
    """Parse PubMed XML ChemicalList entries keyed by PMID."""
    import re as _re
    out: dict[str, list[dict[str, str]]] = {}
    if not raw_xml:
        return out
    # Split roughly by PubmedArticle
    articles = _re.split(r"</PubmedArticle>", raw_xml)
    for art in articles:
        m_pmid = _re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
        if not m_pmid:
            continue
        pmid = m_pmid.group(1)
        chems: list[dict[str, str]] = []
        for cm in _re.finditer(
            r"<Chemical>.*?<NameOfSubstance[^>]*>(.*?)</NameOfSubstance>.*?</Chemical>",
            art,
            _re.I | _re.S,
        ):
            name = _re.sub(r"<[^>]+>", "", cm.group(1)).strip()
            if name:
                chems.append({"name": name, "source": "pubmed_chemical_list"})
        if chems:
            out[pmid] = chems
    return out


def extract_pmid_title_abstract(raw_xml: str) -> dict[str, dict[str, str]]:
    """Return {pmid: {title, abstract}} from PubMed XML."""
    import re as _re
    out: dict[str, dict[str, str]] = {}
    if not raw_xml:
        return out
    articles = _re.split(r"</PubmedArticle>", raw_xml)
    for art in articles:
        m_pmid = _re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
        if not m_pmid:
            continue
        pmid = m_pmid.group(1)
        m_title = _re.search(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", art, _re.I | _re.S)
        title = _re.sub(r"<[^>]+>", "", m_title.group(1)).strip() if m_title else ""
        abstracts = _re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, _re.I | _re.S)
        abstract = " ".join(_re.sub(r"<[^>]+>", "", a) for a in abstracts)
        out[pmid] = {"title": title, "abstract": abstract}
    return out


__all__ = [
    "SOURCE_NAME",
    "EUTILS_BASE",
    "DEFAULT_RETMAX",
    "HD_MESH",
    "PubMedGeneTerms",
    "build_gene_title_abstract_clause",
    "build_hd_context_clause",
    "build_hd_search_term",
    "esearch",
    "esearch_custom",
    "esummary",
    "efetch",
    "efetch_abstracts",
    "extract_id_list",
    "efetch_pubmed_xml",
    "extract_medline_chemicals",
    "extract_pmid_title_abstract",
    "search_hd_literature",
]
