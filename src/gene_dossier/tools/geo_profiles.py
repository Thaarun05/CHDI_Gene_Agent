"""GEO Profiles client for Section 3a (brain/neuron perturbation screening).

Independent of ``tools/geo.py``. Never calls geoprofiles EFetch XML.
Never uses profileIcon.cgi for polished charts.
"""

from __future__ import annotations

import io
import re
from datetime import datetime, timezone
from typing import Any, Sequence
from urllib.parse import parse_qs, urlencode, urljoin, urlparse

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "GEO Profiles"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_HOST = "www.ncbi.nlm.nih.gov"
PROFILE_PAGE_TMPL = "https://www.ncbi.nlm.nih.gov/geoprofiles/{uid}"
GRAPH_CGI_PATH = "/geo/tools/profileGraph.cgi"
ICON_CGI_PATH = "/geo/tools/profileIcon.cgi"

DEFAULT_MAX_DISCOVERY = 500
DEFAULT_MAX_SELECTED = 6
DEFAULT_ESUMMARY_BATCH = 100
DEFAULT_RETMAX_PAGE = 100

GRAPH_STATUS_SUCCESS = "success"
GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE = "not_attempted_outside_shortlist"
GRAPH_STATUS_NOT_ATTEMPTED_OPTIONAL = "not_attempted_optional"
GRAPH_STATUS_FAILED = "failed"

NCBI_HOME_URL = "https://www.ncbi.nlm.nih.gov/"
NCBI_APPROVED_HOSTS = frozenset({NCBI_HOST, "ncbi.nlm.nih.gov"})
USER_AGENT = "GeneDossier/0.1.0"
MIN_CHART_WIDTH = 400
MIN_CHART_HEIGHT = 200

NCBI_BLOCK_TOKENS = (
    "unusual browser activity",
    "possible abnormal browsing activity",
    "misuse.ncbi.nlm.nih.gov",
    "blocking.shtml",
    "your request was blocked",
    "not enough information for us to verify",
    "hhs vulnerability disclosure",
)
# Present on normal NCBI chrome; only counts with another block signal.
NCBI_BLOCK_CORROBORATION_ONLY = frozenset(
    {
        "hhs vulnerability disclosure",
    }
)

NEURAL_TERMS = (
    "brain[All Fields]",
    "neuron[All Fields]",
    "neuronal[All Fields]",
    "hippocampus[All Fields]",
    "cortex[All Fields]",
    "cerebellum[All Fields]",
    "striatum[All Fields]",
    '"motor neuron"[All Fields]',
    "motoneuron[All Fields]",
)

PERTURBATION_HINTS = (
    "stress", "treatment", "treated", "knockout", "knockdown", "mutant",
    "deficiency", "effect", "response", "exposed", "infection", "therapy",
    "drug", "hormone", "toxin", "paraquat", "fluoxetine", "lps",
    "lipopolysaccharide", "alcohol", "hiv", "antiretroviral", "thyroid",
)
COMPARATOR_HINTS = (
    "control", "untreated", "wild type", "wild-type", "vehicle", "saline",
    "sham", "placebo",
)


def normalize_gds_uid(value: str | int | None) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    match = re.match(r"^(?:GDS)?(\d+)$", text, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def format_gds_accession(value: str | int | None) -> str:
    uid = normalize_gds_uid(value)
    return f"GDS{uid}" if uid else ""


def build_exact_gene_symbol_query(gene_symbol: str) -> str:
    return f"{gene_symbol.strip()}[Gene Symbol]"


def build_neural_context_query(gene_symbol: str) -> str:
    joined = " OR ".join(NEURAL_TERMS)
    return f"{gene_symbol.strip()}[Gene Symbol] AND ({joined})"


def build_subset_effect_query(gene_symbol: str) -> str:
    return f'{build_neural_context_query(gene_symbol)} AND "value subset effect"[Flag Type]'


def _tool_result(**kwargs: Any) -> ToolResult:
    return ToolResult(source_name=SOURCE_NAME, **kwargs)


def _with_api_key(params: dict[str, Any], settings: Settings) -> dict[str, Any]:
    out = dict(params)
    if settings.has_key("NCBI_API_KEY"):
        out["api_key"] = settings.ncbi_api_key
    return out


def _safe_params(params: dict[str, Any]) -> dict[str, Any]:
    out = {k: v for k, v in params.items() if k != "api_key"}
    if "api_key" in params:
        out["api_key"] = "***"
    return out


def _request(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
    accept: str | None = None,
    headers: dict[str, str] | None = None,
    referer: str | None = None,
) -> ToolResult:
    is_eutils = path.endswith(".fcgi")
    query = _with_api_key(params, settings) if is_eutils else dict(params)
    request_base = f"{EUTILS_BASE}/{path}" if is_eutils else path
    safe = _safe_params(query) if is_eutils else {}
    request_url = f"{request_base}?{urlencode(safe)}" if safe else request_base
    req_headers: dict[str, str] = {}
    if not is_eutils:
        req_headers["User-Agent"] = USER_AGENT
        req_headers["Accept-Language"] = "en-US,en;q=0.9"
    if accept:
        req_headers["Accept"] = accept
    if referer:
        req_headers["Referer"] = referer
    if headers:
        req_headers.update(headers)
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds, follow_redirects=True) as client:
            response = client.get(
                request_base,
                params=query if is_eutils else None,
                headers=req_headers,
            )
        content_type = (response.headers.get("content-type") or "").lower()
        if is_eutils or "json" in content_type:
            try:
                payload: Any = response.json()
            except ValueError:
                payload = {
                    "raw_text": response.text[:8000],
                    "content_type": content_type,
                    "content_bytes": response.content,
                }
        else:
            payload = {
                "content_type": content_type,
                "content_bytes": response.content,
                "raw_text": response.text[:8000] if "html" in content_type or "text/" in content_type else "",
            }
        if not is_eutils and isinstance(payload, dict):
            payload = {
                **payload,
                "requested_url": request_url,
                "final_url": str(response.url),
                "redirect_history": [str(item.url) for item in response.history],
                "content_type": content_type,
            }
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=response.is_success,
            status_code=response.status_code,
            data=payload,
            error_type=None if response.is_success else "http_error",
            error_message=None if response.is_success else f"HTTP {response.status_code}",
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
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def extract_id_list(esearch_result: ToolResult) -> list[str]:
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return []
    ids = (esearch_result.data.get("esearchresult") or {}).get("idlist") or []
    return [str(i) for i in ids if str(i).strip()]


def extract_count(esearch_result: ToolResult) -> int | None:
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return None
    try:
        return int((esearch_result.data.get("esearchresult") or {}).get("count"))
    except (TypeError, ValueError):
        return None


def extract_querytranslation(esearch_result: ToolResult) -> str | None:
    if not esearch_result.success or not isinstance(esearch_result.data, dict):
        return None
    value = (esearch_result.data.get("esearchresult") or {}).get("querytranslation")
    return str(value) if value else None


def esearch_geoprofiles(
    *,
    gene_symbol: str,
    term: str,
    retmax: int = 0,
    retstart: int = 0,
    sort: str | None = "relevance",
    settings: Settings | None = None,
) -> ToolResult:
    cfg = settings or get_settings()
    params: dict[str, Any] = {
        "db": "geoprofiles",
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
        "retstart": str(retstart),
    }
    if sort:
        params["sort"] = sort
    return _request(
        endpoint_name="esearch_geoprofiles",
        gene_symbol=gene_symbol,
        path="esearch.fcgi",
        params=params,
        settings=cfg,
    )


def page_geoprofile_ids(
    *,
    gene_symbol: str,
    term: str,
    max_ids: int,
    settings: Settings | None = None,
) -> dict[str, Any]:
    cfg = settings or get_settings()
    collected: list[str] = []
    total_count: int | None = None
    querytranslation: str | None = None
    runs: list[ToolResult] = []
    truncated = False
    retstart = 0
    while len(collected) < max_ids:
        page_size = min(DEFAULT_RETMAX_PAGE, max_ids - len(collected))
        result = esearch_geoprofiles(
            gene_symbol=gene_symbol,
            term=term,
            retmax=page_size,
            retstart=retstart,
            settings=cfg,
        )
        runs.append(result)
        if not result.success:
            break
        if total_count is None:
            total_count = extract_count(result)
        if querytranslation is None:
            querytranslation = extract_querytranslation(result)
        ids = extract_id_list(result)
        if not ids:
            break
        for uid in ids:
            if uid not in collected:
                collected.append(uid)
            if len(collected) >= max_ids:
                break
        retstart += len(ids)
        if total_count is not None and retstart >= total_count:
            break
        if len(ids) < page_size:
            break
    if total_count is not None and total_count > len(collected):
        truncated = True
    return {
        "ids": collected,
        "count": total_count if total_count is not None else len(collected),
        "querytranslation": querytranslation,
        "truncated": truncated,
        "tool_results": runs,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }


def esummary_geoprofiles(
    profile_ids: Sequence[str],
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    cfg = settings or get_settings()
    id_str = ",".join(str(i).strip() for i in profile_ids if str(i).strip())
    if not id_str:
        return _tool_result(
            endpoint_name="esummary_geoprofiles",
            gene_symbol=gene_symbol,
            request_url=f"{EUTILS_BASE}/esummary.fcgi",
            request_params={"db": "geoprofiles", "id": ""},
            success=False,
            error_type="invalid_request",
            error_message="Profiles ESummary requires at least one ID",
        )
    return _request(
        endpoint_name="esummary_geoprofiles",
        gene_symbol=gene_symbol or id_str.split(",")[0],
        path="esummary.fcgi",
        params={"db": "geoprofiles", "id": id_str, "retmode": "json"},
        settings=cfg,
    )


def esummary_gds(
    gds_uids: Sequence[str],
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    cfg = settings or get_settings()
    normalized = [normalize_gds_uid(uid) for uid in gds_uids]
    id_str = ",".join(dict.fromkeys(uid for uid in normalized if uid))
    if not id_str:
        return _tool_result(
            endpoint_name="esummary_gds",
            gene_symbol=gene_symbol,
            request_url=f"{EUTILS_BASE}/esummary.fcgi",
            request_params={"db": "gds", "id": ""},
            success=False,
            error_type="invalid_request",
            error_message="GDS ESummary requires at least one numeric GDS UID",
        )
    return _request(
        endpoint_name="esummary_gds",
        gene_symbol=gene_symbol or id_str.split(",")[0],
        path="esummary.fcgi",
        params={"db": "gds", "id": id_str, "retmode": "json"},
        settings=cfg,
    )


def elink_profile_to_gds(
    profile_uid: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    cfg = settings or get_settings()
    uid = str(profile_uid).strip()
    if not uid:
        return _tool_result(
            endpoint_name="elink_profile_to_gds",
            gene_symbol=gene_symbol,
            request_url=f"{EUTILS_BASE}/elink.fcgi",
            request_params={"dbfrom": "geoprofiles", "db": "gds", "id": ""},
            success=False,
            error_type="invalid_request",
            error_message="ELink requires a profile UID",
        )
    return _request(
        endpoint_name="elink_profile_to_gds",
        gene_symbol=gene_symbol or uid,
        path="elink.fcgi",
        params={"dbfrom": "geoprofiles", "db": "gds", "id": uid, "retmode": "json"},
        settings=cfg,
    )


def extract_elink_gds_uids(elink_result: ToolResult) -> list[str]:
    if not elink_result.success or not isinstance(elink_result.data, dict):
        return []
    out: list[str] = []
    for linkset in elink_result.data.get("linksets") or []:
        for db in linkset.get("linksetdbs") or []:
            if str(db.get("linkto") or db.get("dbto") or "").lower() != "gds":
                continue
            for link in db.get("links") or []:
                uid = normalize_gds_uid(link)
                if uid:
                    out.append(uid)
    return out


def parse_profile_esummary_document(doc: dict[str, Any]) -> dict[str, Any]:
    fields = (
        "uid", "gds", "gpl", "title", "taxon", "gdstype", "valtype", "idref",
        "genename", "genedesc", "geneid", "gbacc", "ptacc", "cloneid", "spotid",
        "vmin", "vmax", "groups", "erank", "evalue", "abscall", "aflag", "rstd", "rmean",
    )
    out = {k: doc.get(k) for k in fields if k in doc}
    for alt in ("geneid", "GeneID", "entrezgene"):
        if doc.get(alt) is not None:
            out["entrez_gene_id"] = doc.get(alt)
            break
    out["gds_uid"] = normalize_gds_uid(out.get("gds"))
    out["gds_accession"] = format_gds_accession(out.get("gds"))
    out["profile_uid"] = str(out.get("uid") or "").strip()
    if out["profile_uid"]:
        out["profile_page_url"] = PROFILE_PAGE_TMPL.format(uid=out["profile_uid"])
    return out


def parse_profile_esummary_payload(payload: dict[str, Any]) -> list[dict[str, Any]]:
    result = payload.get("result") or {}
    docs: list[dict[str, Any]] = []
    for uid in [str(u) for u in (result.get("uids") or [])]:
        doc = result.get(uid)
        if isinstance(doc, dict):
            parsed = parse_profile_esummary_document(doc)
            if not parsed.get("profile_uid"):
                parsed["profile_uid"] = uid
                parsed["uid"] = uid
                parsed["profile_page_url"] = PROFILE_PAGE_TMPL.format(uid=uid)
            docs.append(parsed)
    return docs


def parse_gds_esummary_document(doc: dict[str, Any]) -> dict[str, Any]:
    samples = doc.get("samples") or doc.get("samplelist") or []
    if isinstance(samples, dict):
        sample_records = list(samples.values())
    elif isinstance(samples, list):
        sample_records = samples
    else:
        sample_records = []
    explicit = doc.get("n_samples") or doc.get("nsamples") or doc.get("samplecount")
    try:
        sample_count = int(explicit) if explicit is not None else len(sample_records)
    except (TypeError, ValueError):
        sample_count = len(sample_records)
    accession = str(doc.get("accession") or doc.get("gds") or "").strip()
    gds_uid = normalize_gds_uid(accession) or normalize_gds_uid(doc.get("uid"))
    return {
        "uid": str(doc.get("uid") or gds_uid or ""),
        "gds_uid": gds_uid,
        "gds_accession": format_gds_accession(gds_uid or accession),
        "gse_accession": doc.get("gse") or doc.get("gseaccession"),
        "gpl_accession": doc.get("gpl") or doc.get("platform"),
        "title": doc.get("title"),
        "summary": doc.get("summary"),
        "organism": doc.get("taxon") or doc.get("organism"),
        "dataset_type": doc.get("gdstype") or doc.get("entrytype"),
        "value_type": doc.get("valtype") or doc.get("valuetype"),
        "sample_count": sample_count,
        "samples": sample_records,
        "pubmed_id": doc.get("pubmedids") or doc.get("pubmed_id") or doc.get("pubmed"),
        "platform_title": doc.get("platformtitle") or doc.get("platform_title"),
        "pdat": doc.get("pdat") or doc.get("pdate") or doc.get("pubdate"),
        "bioproject": doc.get("bioproject"),
        "subsetinfo": doc.get("subsetinfo") or doc.get("subsets"),
        "raw": doc,
    }


def parse_gds_esummary_payload(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result = payload.get("result") or {}
    out: dict[str, dict[str, Any]] = {}
    for uid in result.get("uids") or []:
        doc = result.get(str(uid))
        if isinstance(doc, dict):
            parsed = parse_gds_esummary_document(doc)
            out[str(parsed.get("gds_uid") or uid)] = parsed
    return out


def genes_match(requested: str, returned: str | None) -> bool:
    return (requested or "").strip().lower() == (str(returned or "").strip().lower())


def _text_blob(*parts: Any) -> str:
    return " ".join(str(p).lower() for p in parts if p)


def _has_any(text: str, hints: Sequence[str]) -> bool:
    return any(h in text for h in hints)


def assess_neural_context(profile: dict[str, Any], gds: dict[str, Any] | None) -> bool:
    blob = _text_blob(
        profile.get("title"), profile.get("taxon"),
        (gds or {}).get("title"), (gds or {}).get("summary"), (gds or {}).get("organism"),
    )
    return _has_any(blob, (
        "brain", "neuron", "neuronal", "hippocamp", "cortex", "cerebell",
        "striatum", "motoneuron", "motor neuron", "dentate", "ganglia",
        "white matter", "frontal",
    ))


def assess_perturbation_design(profile: dict[str, Any], gds: dict[str, Any] | None) -> bool:
    blob = _text_blob(
        profile.get("title"), (gds or {}).get("title"),
        (gds or {}).get("summary"), (gds or {}).get("subsetinfo"),
    )
    return _has_any(blob, PERTURBATION_HINTS) and _has_any(blob, COMPARATOR_HINTS)


def score_profile(
    profile: dict[str, Any],
    *,
    gds: dict[str, Any] | None,
    subset_effect_flag: bool,
    graph_ok: bool | None = None,
) -> dict[str, Any]:
    components: dict[str, int] = {}
    title = str(profile.get("title") or "")
    blob = _text_blob(title, (gds or {}).get("title"), (gds or {}).get("summary"))
    if assess_neural_context(profile, gds):
        components["neural_context"] = 25
        if _has_any(title.lower(), ("brain", "neuron", "hippocamp", "cortex", "cerebell", "striatum")):
            components["neural_in_title"] = 10
    if assess_perturbation_design(profile, gds):
        components["perturbation_comparator"] = 25
        if _has_any(title.lower(), PERTURBATION_HINTS):
            components["perturbation_in_title"] = 10
    if profile.get("gds_uid") and profile.get("gpl") and profile.get("idref"):
        components["reporter_metadata"] = 15
    if gds and gds.get("sample_count"):
        components["sample_metadata"] = 8
    if subset_effect_flag:
        components["subset_effect_flag"] = 12
    if gds and (gds.get("pubmed_id") or gds.get("pubmedids")):
        components["pubmed_linked"] = 6
    if _has_any(blob, ("disease", "therapy", "toxic", "knockout", "knockdown", "hormone", "infection", "drug", "stress")):
        components["mechanistic_relevance"] = 8
    if graph_ok is True:
        components["validated_chart"] = 15
    score = int(sum(components.values()))
    return {"base_score": score, "score_components": components, "final_score": score}


def eligibility_status(
    profile: dict[str, Any],
    *,
    requested_gene: str,
    gds: dict[str, Any] | None,
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    if not genes_match(requested_gene, profile.get("genename")):
        return "gene_mismatch", ["returned genename does not match requested symbol"]
    if not profile.get("gds_uid"):
        reasons.append("missing_gds")
    if not profile.get("gpl"):
        reasons.append("missing_gpl")
    if not profile.get("idref"):
        reasons.append("missing_idref")
    if not assess_neural_context(profile, gds):
        reasons.append("neural_context_unresolved")
    if not assess_perturbation_design(profile, gds):
        reasons.append("perturbation_comparator_unresolved")
    if not (profile.get("title") and profile.get("taxon") and profile.get("genename")):
        reasons.append("incomplete_presentation_metadata")
    if reasons:
        return "ineligible", reasons
    return "eligible", []


def _diversity_keys(profile: dict[str, Any], gds: dict[str, Any] | None) -> dict[str, str]:
    organism = str(profile.get("taxon") or (gds or {}).get("organism") or "unknown").lower()
    title = str(profile.get("title") or "").lower()
    if "hiv" in title or "hand" in title:
        category = "disease_hiv"
    elif "alcohol" in title:
        category = "toxicant_alcohol"
    elif "stress" in title:
        category = "stress"
    elif "fluoxetine" in title or "antidepressant" in title:
        category = "drug_treatment"
    elif "knock" in title or "mutant" in title or "deficien" in title:
        category = "genetic_perturbation"
    elif "hormone" in title or "thyroid" in title:
        category = "hormone"
    elif "lps" in title or "immune" in title:
        category = "immune_challenge"
    else:
        category = "other_perturbation"
    region = "general_neural"
    for token, label in (
        ("hippocamp", "hippocampus"), ("cortex", "cortex"), ("cerebell", "cerebellum"),
        ("dentate", "dentate"), ("motoneuron", "motoneuron"), ("motor neuron", "motoneuron"),
        ("striatum", "striatum"), ("white matter", "white_matter"), ("basal gangli", "basal_ganglia"),
    ):
        if token in title:
            region = label
            break
    return {"organism": organism, "category": category, "region": region}


def _sort_key(item: dict[str, Any]) -> tuple:
    return (
        -int(item.get("final_score") or 0),
        0 if item.get("graph_ok") else 1,
        str(item.get("pdat") or "0000"),
        str(item.get("profile_uid") or ""),
    )


def build_diversity_shortlist(
    ranked: Sequence[dict[str, Any]],
    *,
    max_items: int,
) -> list[dict[str, Any]]:
    if max_items <= 0:
        return []
    selected: list[dict[str, Any]] = []
    seen_gds: set[str] = set()
    used_org: set[str] = set()
    used_cat: set[str] = set()
    used_region: set[str] = set()

    def _try_add(item: dict[str, Any], *, require_novel: bool) -> bool:
        gds = str(item.get("gds_uid") or "")
        if gds and gds in seen_gds:
            existing = [s for s in selected if s.get("gds_uid") == gds]
            if existing and existing[0].get("idref") == item.get("idref"):
                return False
        keys = item.get("diversity_keys") or {}
        if require_novel and selected:
            novel = (
                keys.get("organism") not in used_org
                or keys.get("category") not in used_cat
                or keys.get("region") not in used_region
            )
            if not novel:
                return False
        selected.append(item)
        if gds:
            seen_gds.add(gds)
        used_org.add(str(keys.get("organism") or ""))
        used_cat.add(str(keys.get("category") or ""))
        used_region.add(str(keys.get("region") or ""))
        return True

    ordered = sorted(ranked, key=_sort_key)
    for item in ordered:
        if len(selected) >= max_items:
            break
        _try_add(item, require_novel=True)
    for item in ordered:
        if len(selected) >= max_items:
            break
        if item in selected:
            continue
        _try_add(item, require_novel=False)
    return selected


def expected_graph_id(gds_uid: str, idref: str) -> str:
    return f"{format_gds_accession(gds_uid)}:{idref}"


def detect_ncbi_block_page(
    *,
    final_url: str | None = None,
    page_title: str | None = None,
    body_text: str | None = None,
    html: str | None = None,
) -> tuple[bool, str | None]:
    """Return (blocked, reason) for NCBI misuse / unusual-activity pages."""
    url = str(final_url or "")
    host = (urlparse(url).hostname or "").lower()
    if host == "misuse.ncbi.nlm.nih.gov":
        return True, "misuse_host"
    if "blocking.shtml" in url.lower():
        return True, "blocking_shtml"
    title = str(page_title or "").strip()
    if title and title.casefold() in {"error", "ncbi error", "access denied", "blocked"}:
        return True, "error_title"
    haystacks = (
        title,
        str(body_text or ""),
        str(html or ""),
        url,
    )
    joined = "\n".join(haystacks).casefold()
    matched: list[str] = []
    for token in NCBI_BLOCK_TOKENS:
        if token.casefold() in joined:
            matched.append(token)
    primary = [t for t in matched if t.casefold() not in NCBI_BLOCK_CORROBORATION_ONLY]
    if primary:
        return True, f"block_token:{primary[0]}"
    # Corroboration-only tokens (e.g. HHS footer) require another signal.
    if matched and (
        host == "misuse.ncbi.nlm.nih.gov"
        or "blocking.shtml" in url.lower()
        or title.casefold() in {"error", "ncbi error", "access denied", "blocked"}
    ):
        return True, f"block_token:{matched[0]}"
    return False, None


def validate_graph_url(url: str, *, gds_uid: str, idref: str) -> tuple[bool, str | None]:
    parsed = urlparse(url)
    host = (parsed.hostname or "").lower()
    if host not in NCBI_APPROVED_HOSTS:
        return False, "non_ncbi_host"
    if parsed.username or parsed.password:
        return False, "userinfo_forbidden"
    path = parsed.path or ""
    if ICON_CGI_PATH in path:
        return False, "graph_thumbnail_rejected"
    if GRAPH_CGI_PATH not in path:
        return False, "graph_link_missing"
    qs = parse_qs(parsed.query)
    raw_id = (qs.get("ID") or qs.get("id") or [""])[0]
    if not raw_id:
        return False, "missing_id_parameter"
    if raw_id != expected_graph_id(gds_uid, idref):
        return False, "graph_identifier_mismatch"
    return True, None


def discover_graph_url_from_html(html: str, *, gds_uid: str, idref: str) -> dict[str, Any]:
    pattern = re.compile(
        r'href=["\']([^"\']*?/geo/tools/profileGraph\.cgi\?ID=[^"\']+)["\']',
        flags=re.IGNORECASE,
    )
    matches = pattern.findall(html or "")
    last_err = None
    for href in matches:
        absolute = urljoin(f"https://{NCBI_HOST}/", href)
        ok, err = validate_graph_url(absolute, gds_uid=gds_uid, idref=idref)
        if ok:
            return {"ok": True, "url": absolute, "origin": "discovered_from_profile_html", "error_type": None}
        last_err = err
    if matches:
        return {"ok": False, "url": None, "origin": "discovered_from_profile_html", "error_type": last_err}
    return {"ok": False, "url": None, "origin": None, "error_type": "graph_link_missing"}


def construct_graph_url(gds_uid: str, idref: str) -> str:
    return f"https://{NCBI_HOST}{GRAPH_CGI_PATH}?ID={expected_graph_id(gds_uid, idref)}"


def _image_metrics(content: bytes) -> dict[str, Any]:
    try:
        from PIL import Image
        img = Image.open(io.BytesIO(content))
        width, height = img.size
        rgb = img.convert("RGB")
        pixels = list(rgb.getdata())
        sample = pixels[:: max(1, len(pixels) // 5000)] or pixels
        nonwhite = sum(1 for p in sample if p[0] < 250 or p[1] < 250 or p[2] < 250)
        # NCBI misuse pages are typically large red/black banners; chart pages have
        # mixed scientific content. Track dominant dark-red fraction as a soft signal.
        dark_red = sum(
            1
            for p in sample
            if p[0] > 140 and p[1] < 80 and p[2] < 80
        )
        near_black = sum(1 for p in sample if p[0] < 40 and p[1] < 40 and p[2] < 40)
        return {
            "ok": True,
            "width": width,
            "height": height,
            "nonwhite_fraction": round(nonwhite / max(1, len(sample)), 4),
            "dark_red_fraction": round(dark_red / max(1, len(sample)), 4),
            "near_black_fraction": round(near_black / max(1, len(sample)), 4),
        }
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": str(exc)}


def validate_chart_bytes(content: bytes) -> tuple[bool, str | None, dict[str, Any]]:
    """Dimension/nonblank checks only. Not sufficient for graph_status=success."""
    metrics = _image_metrics(content)
    if not metrics.get("ok"):
        return False, "graph_invalid_content_type", metrics
    if int(metrics.get("width") or 0) < MIN_CHART_WIDTH or int(metrics.get("height") or 0) < MIN_CHART_HEIGHT:
        return False, "graph_too_small", metrics
    if float(metrics.get("nonwhite_fraction") or 0) < 0.01:
        return False, "graph_blank", metrics
    return True, None, metrics


def _crop_outer_blank_margins(content: bytes) -> bytes:
    """Crop only outer near-white margins; preserve legend/condition bars."""
    try:
        from PIL import Image
        img = Image.open(io.BytesIO(content)).convert("RGB")
        width, height = img.size
        pixels = img.load()
        def _blank(x: int, y: int) -> bool:
            r, g, b = pixels[x, y]
            return r >= 250 and g >= 250 and b >= 250
        top = 0
        while top < height - 1 and all(_blank(x, top) for x in range(width)):
            top += 1
        bottom = height - 1
        while bottom > top and all(_blank(x, bottom) for x in range(width)):
            bottom -= 1
        left = 0
        while left < width - 1 and all(_blank(left, y) for y in range(top, bottom + 1)):
            left += 1
        right = width - 1
        while right > left and all(_blank(right, y) for y in range(top, bottom + 1)):
            right -= 1
        # Keep a small pad so legend/axes are not clipped.
        pad = 2
        box = (
            max(0, left - pad),
            max(0, top - pad),
            min(width, right + 1 + pad),
            min(height, bottom + 1 + pad),
        )
        if box[2] - box[0] < MIN_CHART_WIDTH or box[3] - box[1] < MIN_CHART_HEIGHT:
            return content
        cropped = img.crop(box)
        buf = io.BytesIO()
        cropped.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:  # noqa: BLE001
        return content


def validate_graph_page_identity(
    *,
    gds_uid: str,
    idref: str,
    final_url: str | None,
    page_title: str | None = None,
    body_text: str | None = None,
    html: str | None = None,
    require_graph_path: bool = True,
) -> tuple[bool, str | None, dict[str, Any]]:
    """URL/title/block identity checks without requiring image bytes."""
    checks: dict[str, Any] = {
        "final_host_ok": False,
        "graph_id_ok": False,
        "block_page_absent": False,
        "page_title_ok": False,
        "thumbnail_absent": True,
        "identity_ok": False,
    }
    blocked, block_reason = detect_ncbi_block_page(
        final_url=final_url,
        page_title=page_title,
        body_text=body_text,
        html=html,
    )
    checks["block_page_absent"] = not blocked
    if blocked:
        return False, "graph_http_blocked", {**checks, "block_reason": block_reason}

    title = str(page_title or "").strip()
    title_cf = title.casefold()
    checks["page_title_ok"] = title_cf not in {"error", "ncbi error", "access denied", "blocked"}
    if not checks["page_title_ok"]:
        return False, "graph_http_blocked", checks

    parsed = urlparse(str(final_url or ""))
    host = (parsed.hostname or "").lower()
    checks["final_host_ok"] = host in NCBI_APPROVED_HOSTS
    if not checks["final_host_ok"]:
        return False, "non_ncbi_host", checks

    path = parsed.path or ""
    if ICON_CGI_PATH in path:
        checks["thumbnail_absent"] = False
        return False, "graph_thumbnail_rejected", checks
    if require_graph_path and GRAPH_CGI_PATH not in path:
        return False, "graph_link_missing", checks
    if require_graph_path or GRAPH_CGI_PATH in path:
        qs = parse_qs(parsed.query)
        raw_id = (qs.get("ID") or qs.get("id") or [""])[0]
        checks["graph_id_ok"] = raw_id == expected_graph_id(gds_uid, idref)
        if not checks["graph_id_ok"]:
            return False, "graph_identifier_mismatch", checks
    else:
        checks["graph_id_ok"] = False
        return False, "graph_link_missing", checks

    checks["identity_ok"] = True
    return True, None, checks


def validate_graph_capture(
    content: bytes | None,
    *,
    gds_uid: str,
    idref: str,
    final_url: str | None,
    page_title: str | None = None,
    body_text: str | None = None,
    html: str | None = None,
    capture_method: str | None = None,
    require_graph_path: bool = True,
) -> tuple[bool, str | None, dict[str, Any]]:
    """Positive identity + image checks required for graph_status=success."""
    ok_id, err_id, checks = validate_graph_page_identity(
        gds_uid=gds_uid,
        idref=idref,
        final_url=final_url,
        page_title=page_title,
        body_text=body_text,
        html=html,
        require_graph_path=require_graph_path,
    )
    checks = {
        **checks,
        "dimensions_ok": False,
        "nonblank_ok": False,
        "capture_method": capture_method,
    }
    if not ok_id:
        return False, err_id, checks

    if capture_method == "generic_body" and not checks["identity_ok"]:
        return False, "graph_capture_identity_missing", checks

    if not isinstance(content, (bytes, bytearray)) or not content:
        return False, "graph_invalid_content_type", checks

    ok, err, metrics = validate_chart_bytes(bytes(content))
    checks["dimensions_ok"] = bool(ok) or err != "graph_too_small"
    checks["nonblank_ok"] = bool(ok) or err != "graph_blank"
    checks.update({f"image_{k}": v for k, v in metrics.items() if k != "ok"})
    if not ok:
        return False, err or "graph_invalid_content_type", checks

    # Soft reject banner-like captures even if URL somehow matched (defense in depth).
    if (
        float(metrics.get("dark_red_fraction") or 0) > 0.20
        and float(metrics.get("near_black_fraction") or 0) > 0.40
    ):
        return False, "graph_http_blocked", {**checks, "banner_like": True}

    return True, None, checks


def _is_rejected_img_url(url: str) -> bool:
    lower = (url or "").lower()
    if ICON_CGI_PATH.lower() in lower:
        return True
    if any(tok in lower for tok in ("logo", "spacer", "pixel.gif", "header", "footer", "icon")):
        return True
    return False


def fetch_profile_html(profile_uid: str, *, gene_symbol: str, settings: Settings | None = None) -> ToolResult:
    cfg = settings or get_settings()
    return _request(
        endpoint_name="profile_html",
        gene_symbol=gene_symbol,
        path=PROFILE_PAGE_TMPL.format(uid=profile_uid),
        params={},
        settings=cfg,
        accept="text/html",
    )


def _launch_playwright_chromium(pw: Any) -> tuple[Any, str]:
    try:
        browser = pw.chromium.launch(headless=True, channel="chrome")
        return browser, "chrome"
    except Exception:  # noqa: BLE001
        browser = pw.chromium.launch(headless=True)
        return browser, "chromium"


class GeoProfilesBrowserSession:
    """Reusable Chromium context for shortlist chart acquisition."""

    def __init__(self, *, viewport: tuple[int, int] = (1280, 900), headless: bool = True) -> None:
        self.viewport = {"width": int(viewport[0]), "height": int(viewport[1])}
        self.headless = headless
        self._pw: Any = None
        self._browser: Any = None
        self._context: Any = None
        self.browser_channel: str | None = None
        self.browser_version: str | None = None
        self.home_ok = False
        self._capture_count = 0

    def __enter__(self) -> "GeoProfilesBrowserSession":
        from playwright.sync_api import sync_playwright

        self._pw = sync_playwright().start()
        self._browser, self.browser_channel = _launch_playwright_chromium(self._pw)
        self.browser_version = str(self._browser.version)
        self._context = self._browser.new_context(
            viewport=self.viewport,
            user_agent=USER_AGENT,
            locale="en-US",
        )
        page = self._context.new_page()
        try:
            page.goto(NCBI_HOME_URL, wait_until="domcontentloaded", timeout=60000)
            blocked, reason = detect_ncbi_block_page(
                final_url=page.url,
                page_title=page.title(),
                body_text=page.locator("body").inner_text(timeout=5000) if page.locator("body").count() else "",
                html=page.content(),
            )
            if blocked:
                raise RuntimeError(f"ncbi_home_blocked:{reason}")
            self.home_ok = True
        finally:
            page.close()
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        try:
            if self._context is not None:
                self._context.close()
        except Exception:  # noqa: BLE001
            pass
        try:
            if self._browser is not None:
                self._browser.close()
        except Exception:  # noqa: BLE001
            pass
        try:
            if self._pw is not None:
                self._pw.stop()
        except Exception:  # noqa: BLE001
            pass
        self._context = None
        self._browser = None
        self._pw = None

    @property
    def capture_count(self) -> int:
        return self._capture_count

    def capture_profile_graph(
        self,
        *,
        profile_uid: str,
        gds_uid: str,
        idref: str,
        graph_url: str,
    ) -> dict[str, Any]:
        if self._context is None:
            return {
                "acquisition_method": "playwright_profile_navigation",
                "graph_status": GRAPH_STATUS_FAILED,
                "error_type": "browser_unavailable",
                "image_bytes": None,
                "tool_results": [],
            }
        out: dict[str, Any] = {
            "acquisition_method": "playwright_profile_navigation",
            "browser_channel": self.browser_channel,
            "browser_version": self.browser_version,
            "viewport": dict(self.viewport),
            "profile_page_requested_url": PROFILE_PAGE_TMPL.format(uid=profile_uid),
            "profile_page_final_url": None,
            "graph_requested_url": graph_url,
            "graph_final_url": None,
            "page_title": None,
            "graph_url_origin": None,
            "capture_selector": None,
            "capture_method": None,
            "validation_checks": {},
            "graph_status": GRAPH_STATUS_FAILED,
            "image_bytes": None,
            "content_type": None,
            "image_width": None,
            "image_height": None,
            "error_type": None,
            "tool_results": [],
        }
        page = self._context.new_page()
        popup_page = None
        try:
            profile_url = PROFILE_PAGE_TMPL.format(uid=profile_uid)
            page.goto(profile_url, wait_until="domcontentloaded", timeout=60000)
            out["profile_page_final_url"] = page.url
            profile_html = page.content()
            profile_title = page.title()
            try:
                profile_body = page.locator("body").inner_text(timeout=5000)
            except Exception:  # noqa: BLE001
                profile_body = ""
            blocked, reason = detect_ncbi_block_page(
                final_url=page.url,
                page_title=profile_title,
                body_text=profile_body,
                html=profile_html,
            )
            if blocked:
                out["error_type"] = "graph_http_blocked"
                out["page_title"] = profile_title
                out["validation_checks"] = {"block_reason": reason, "block_page_absent": False}
                return out
            if str(profile_uid) not in str(page.url) and str(profile_uid) not in profile_html:
                out["error_type"] = "profile_page_mismatch"
                return out

            expected_id = expected_graph_id(gds_uid, idref)
            ok_url, url_err = validate_graph_url(graph_url, gds_uid=gds_uid, idref=idref)
            if not ok_url:
                out["error_type"] = url_err or "graph_link_missing"
                return out

            # Prefer the exact validated full-graph anchor on the profile page.
            href_selector = f'a[href*="profileGraph.cgi"][href*="{expected_id}"]'
            link = page.locator(href_selector)
            if link.count() == 0:
                # Fall back to any validated profileGraph href discovered in HTML.
                discovered = discover_graph_url_from_html(
                    profile_html, gds_uid=gds_uid, idref=idref
                )
                if not discovered.get("ok"):
                    out["error_type"] = discovered.get("error_type") or "graph_link_missing"
                    return out
                graph_url = str(discovered["url"])
                out["graph_requested_url"] = graph_url
                out["graph_url_origin"] = "discovered_from_profile_html"
                link = page.locator(href_selector)
                if link.count() == 0:
                    out["error_type"] = "graph_link_missing"
                    return out
            else:
                out["graph_url_origin"] = "discovered_from_profile_html"

            target = link.first
            href = target.get_attribute("href") or ""
            absolute = urljoin(page.url, href)
            ok_href, href_err = validate_graph_url(absolute, gds_uid=gds_uid, idref=idref)
            if not ok_href:
                out["error_type"] = href_err or "graph_identifier_mismatch"
                return out
            out["graph_requested_url"] = absolute

            graph_page = page
            popup_page = None
            clicked = False
            try:
                with self._context.expect_page(timeout=3000) as popup_info:
                    target.click(timeout=15000)
                    clicked = True
                popup_page = popup_info.value
                graph_page = popup_page
            except Exception:  # noqa: BLE001
                # Same-tab navigation. Re-click if the popup waiter aborted before click.
                if not clicked:
                    try:
                        target.click(timeout=15000)
                        clicked = True
                    except Exception:  # noqa: BLE001
                        pass
                try:
                    page.wait_for_load_state("domcontentloaded", timeout=30000)
                except Exception:  # noqa: BLE001
                    pass
                graph_page = page

            try:
                graph_page.wait_for_load_state("domcontentloaded", timeout=30000)
            except Exception:  # noqa: BLE001
                pass

            out["graph_final_url"] = graph_page.url
            out["page_title"] = graph_page.title()
            try:
                body_text = graph_page.locator("body").inner_text(timeout=5000)
            except Exception:  # noqa: BLE001
                body_text = ""
            html = graph_page.content()
            blocked, reason = detect_ncbi_block_page(
                final_url=graph_page.url,
                page_title=out["page_title"],
                body_text=body_text,
                html=html,
            )
            if blocked:
                out["error_type"] = "graph_http_blocked"
                out["validation_checks"] = {"block_reason": reason, "block_page_absent": False}
                out["image_bytes"] = None
                return out

            identity_ok, identity_err, identity_checks = validate_graph_page_identity(
                gds_uid=gds_uid,
                idref=idref,
                final_url=graph_page.url,
                page_title=out["page_title"],
                body_text=body_text,
                html=html,
            )
            if not identity_ok:
                out["error_type"] = identity_err or "graph_http_blocked"
                out["validation_checks"] = identity_checks
                out["image_bytes"] = None
                return out

            png, selector, method = self._screenshot_graph_element(graph_page)
            if png is None:
                out["error_type"] = "capture_failed"
                out["validation_checks"] = identity_checks
                return out
            out["capture_selector"] = selector
            out["capture_method"] = method
            png = _crop_outer_blank_margins(png)
            ok, err, checks = validate_graph_capture(
                png,
                gds_uid=gds_uid,
                idref=idref,
                final_url=graph_page.url,
                page_title=out["page_title"],
                body_text=body_text,
                html=html,
                capture_method=method,
            )
            out["validation_checks"] = checks
            out["content_type"] = "image/png"
            out["image_width"] = checks.get("image_width")
            out["image_height"] = checks.get("image_height")
            if ok:
                out["graph_status"] = GRAPH_STATUS_SUCCESS
                out["image_bytes"] = png
                out["error_type"] = None
                self._capture_count += 1
            else:
                out["graph_status"] = GRAPH_STATUS_FAILED
                out["error_type"] = err or "capture_failed"
                out["image_bytes"] = None
            return out
        except Exception as exc:  # noqa: BLE001
            out["graph_status"] = GRAPH_STATUS_FAILED
            out["error_type"] = "capture_failed"
            out["error_message"] = str(exc)
            out["image_bytes"] = None
            return out
        finally:
            try:
                if popup_page is not None:
                    popup_page.close()
            except Exception:  # noqa: BLE001
                pass
            try:
                page.close()
            except Exception:  # noqa: BLE001
                pass

    def _screenshot_graph_element(self, page: Any) -> tuple[bytes | None, str | None, str | None]:
        """Capture the largest validated scientific chart element; never body-first union."""
        candidates: list[tuple[int, str, Any, str]] = []
        for selector, method in (
            ("img", "img"),
            ("canvas", "canvas"),
            ("svg", "svg"),
            ("#graphic, .graphic, #profileGraph, .profile-graph", "graph_container"),
        ):
            loc = page.locator(selector)
            try:
                count = loc.count()
            except Exception:  # noqa: BLE001
                count = 0
            for i in range(count):
                el = loc.nth(i)
                try:
                    if not el.is_visible():
                        continue
                    box = el.bounding_box()
                    if not box:
                        continue
                    w, h = float(box["width"]), float(box["height"])
                    if w < MIN_CHART_WIDTH or h < MIN_CHART_HEIGHT:
                        continue
                    src = ""
                    try:
                        src = el.get_attribute("src") or el.get_attribute("href") or ""
                    except Exception:  # noqa: BLE001
                        src = ""
                    if src and _is_rejected_img_url(src):
                        continue
                    area = int(w * h)
                    candidates.append((area, f"{selector} >> nth={i}", el, method))
                except Exception:  # noqa: BLE001
                    continue
        candidates.sort(key=lambda item: item[0], reverse=True)
        for _area, sel, el, method in candidates:
            try:
                png = el.screenshot(type="png")
                ok, _, _ = validate_chart_bytes(png)
                if ok:
                    return png, sel, method
            except Exception:  # noqa: BLE001
                continue

        # Body capture only after the caller already validated graph URL identity.
        try:
            body = page.locator("body")
            if body.count():
                png = body.screenshot(type="png")
                ok, _, _ = validate_chart_bytes(png)
                if ok:
                    return png, "body", "validated_body"
        except Exception:  # noqa: BLE001
            pass
        return None, None, None


def acquire_profile_chart(
    *,
    gene_symbol: str,
    profile_uid: str,
    gds_uid: str,
    idref: str,
    profile_html: str | None = None,
    settings: Settings | None = None,
    allow_playwright: bool = True,
    browser_session: GeoProfilesBrowserSession | None = None,
) -> dict[str, Any]:
    cfg = settings or get_settings()
    out: dict[str, Any] = {
        "profile_uid": profile_uid,
        "requested_gene_symbol": gene_symbol,
        "gds_accession": format_gds_accession(gds_uid),
        "idref": idref,
        "graph_status": GRAPH_STATUS_FAILED,
        "graph_url_origin": None,
        "graph_requested_url": None,
        "graph_final_url": None,
        "acquisition_method": None,
        "content_type": None,
        "image_bytes": None,
        "image_width": None,
        "image_height": None,
        "validation_checks": {},
        "error_type": None,
        "tool_results": [],
    }
    # Reject blocked profile HTML before any chart work.
    if profile_html:
        blocked, reason = detect_ncbi_block_page(html=profile_html, body_text=profile_html)
        if blocked:
            out["error_type"] = "graph_http_blocked"
            out["validation_checks"] = {"block_reason": reason, "block_page_absent": False}
            return out

    discovery = discover_graph_url_from_html(profile_html or "", gds_uid=gds_uid, idref=idref)
    if discovery.get("ok"):
        graph_url = str(discovery["url"])
        origin = "discovered_from_profile_html"
    else:
        graph_url = construct_graph_url(gds_uid, idref)
        origin = "constructed_from_validated_metadata"
        ok, err = validate_graph_url(graph_url, gds_uid=gds_uid, idref=idref)
        if not ok:
            out["error_type"] = err or discovery.get("error_type") or "graph_link_missing"
            return out
    out["graph_url_origin"] = origin
    out["graph_requested_url"] = graph_url

    profile_page_url = PROFILE_PAGE_TMPL.format(uid=profile_uid)
    result = _request(
        endpoint_name="profile_graph",
        gene_symbol=gene_symbol,
        path=graph_url,
        params={},
        settings=cfg,
        accept="image/*,*/*",
        referer=profile_page_url,
    )
    out["tool_results"].append(result)
    direct_blocked = False
    if isinstance(result.data, dict):
        final_url = str(result.data.get("final_url") or "")
        content_type = str(result.data.get("content_type") or "")
        content = result.data.get("content_bytes")
        raw_text = str(result.data.get("raw_text") or "")
        redirect_history = list(result.data.get("redirect_history") or [])
        out["graph_final_url"] = final_url or None
        blocked, reason = detect_ncbi_block_page(
            final_url=final_url,
            html=raw_text,
            body_text=raw_text,
        )
        if not blocked:
            for hist in redirect_history:
                blocked, reason = detect_ncbi_block_page(final_url=str(hist))
                if blocked:
                    break
        if blocked:
            direct_blocked = True
            out["error_type"] = "graph_http_blocked"
            out["validation_checks"] = {"block_reason": reason, "block_page_absent": False}
            out["image_bytes"] = None
        elif result.success and isinstance(content, (bytes, bytearray)) and content_type.startswith("image/"):
            ok, err, checks = validate_graph_capture(
                bytes(content),
                gds_uid=gds_uid,
                idref=idref,
                final_url=final_url or graph_url,
                capture_method="direct_image",
            )
            out["validation_checks"] = checks
            out["content_type"] = content_type
            out["image_width"] = checks.get("image_width")
            out["image_height"] = checks.get("image_height")
            if ok:
                out["graph_status"] = GRAPH_STATUS_SUCCESS
                out["acquisition_method"] = "direct_image"
                out["graph_final_url"] = final_url or graph_url
                out["image_bytes"] = _crop_outer_blank_margins(bytes(content))
                return out
            out["error_type"] = err
            out["image_bytes"] = None
        elif result.success and isinstance(content, (bytes, bytearray)) and "html" in content_type:
            html = raw_text or bytes(content).decode("utf-8", errors="replace")
            blocked, reason = detect_ncbi_block_page(
                final_url=final_url,
                html=html,
                body_text=html,
            )
            if blocked:
                direct_blocked = True
                out["error_type"] = "graph_http_blocked"
                out["validation_checks"] = {"block_reason": reason, "block_page_absent": False}
                out["image_bytes"] = None
            else:
                # Never accept the first arbitrary <img>; enumerate and validate.
                img_pattern = re.compile(
                    r'<img[^>]+src=["\']([^"\']+)["\']',
                    flags=re.IGNORECASE,
                )
                accepted = False
                for src in img_pattern.findall(html):
                    if _is_rejected_img_url(src):
                        out["error_type"] = "graph_thumbnail_rejected"
                        continue
                    img_url = urljoin(final_url or graph_url, src)
                    img_host = (urlparse(img_url).hostname or "").lower()
                    if img_host not in NCBI_APPROVED_HOSTS:
                        continue
                    img_result = _request(
                        endpoint_name="profile_graph_image",
                        gene_symbol=gene_symbol,
                        path=img_url,
                        params={},
                        settings=cfg,
                        referer=final_url or graph_url,
                    )
                    out["tool_results"].append(img_result)
                    img_data = img_result.data if isinstance(img_result.data, dict) else {}
                    img_bytes = img_data.get("content_bytes")
                    img_final = str(img_data.get("final_url") or img_url)
                    img_blocked, img_reason = detect_ncbi_block_page(final_url=img_final)
                    if img_blocked:
                        out["error_type"] = "graph_http_blocked"
                        out["validation_checks"] = {
                            "block_reason": img_reason,
                            "block_page_absent": False,
                        }
                        continue
                    if img_result.success and isinstance(img_bytes, (bytes, bytearray)):
                        # Embedded chart images may not keep profileGraph.cgi as final URL;
                        # require the parent graph page identity instead.
                        parent_ok, parent_err, parent_checks = validate_graph_page_identity(
                            gds_uid=gds_uid,
                            idref=idref,
                            final_url=final_url or graph_url,
                            html=html,
                        )
                        if not parent_ok:
                            out["error_type"] = parent_err
                            out["validation_checks"] = parent_checks
                            continue
                        ok, err, metrics = validate_chart_bytes(bytes(img_bytes))
                        checks = {
                            **parent_checks,
                            "dimensions_ok": bool(ok) or err != "graph_too_small",
                            "nonblank_ok": bool(ok) or err != "graph_blank",
                            **{f"image_{k}": v for k, v in metrics.items() if k != "ok"},
                        }
                        out["validation_checks"] = checks
                        out["content_type"] = img_data.get("content_type")
                        out["image_width"] = metrics.get("width")
                        out["image_height"] = metrics.get("height")
                        if ok:
                            out["graph_status"] = GRAPH_STATUS_SUCCESS
                            out["acquisition_method"] = "html_embedded_image"
                            out["graph_final_url"] = final_url or graph_url
                            out["image_bytes"] = _crop_outer_blank_margins(bytes(img_bytes))
                            accepted = True
                            return out
                        out["error_type"] = err
                if not accepted and not out.get("error_type"):
                    out["error_type"] = "graph_invalid_content_type"
                out["image_bytes"] = None
        elif not result.success:
            out["error_type"] = out.get("error_type") or "graph_http_error"
            out["image_bytes"] = None
        else:
            out["error_type"] = out.get("error_type") or "graph_invalid_content_type"
            out["image_bytes"] = None
    else:
        out["error_type"] = "graph_http_error"
        out["image_bytes"] = None

    if out.get("graph_status") == GRAPH_STATUS_SUCCESS:
        return out

    if allow_playwright:
        if browser_session is not None:
            browser = browser_session.capture_profile_graph(
                profile_uid=profile_uid,
                gds_uid=gds_uid,
                idref=idref,
                graph_url=graph_url,
            )
        else:
            try:
                with GeoProfilesBrowserSession() as session:
                    browser = session.capture_profile_graph(
                        profile_uid=profile_uid,
                        gds_uid=gds_uid,
                        idref=idref,
                        graph_url=graph_url,
                    )
            except Exception as exc:  # noqa: BLE001
                browser = {
                    "acquisition_method": "playwright_profile_navigation",
                    "graph_status": GRAPH_STATUS_FAILED,
                    "error_type": "browser_unavailable",
                    "error_message": str(exc),
                    "image_bytes": None,
                    "tool_results": [],
                }
        for k, v in browser.items():
            if k == "tool_results":
                out["tool_results"].extend(v or [])
            elif k == "error_type" or v is not None:
                out[k] = v
        if out.get("graph_status") != GRAPH_STATUS_SUCCESS:
            out["image_bytes"] = None
            if direct_blocked and out.get("error_type") in {None, "capture_failed", "graph_http_error"}:
                out["error_type"] = "graph_http_blocked"
        return out

    if direct_blocked:
        out["error_type"] = "graph_http_blocked"
        out["image_bytes"] = None
    return out


def _html_text(result: ToolResult) -> str:
    if not result.success or not isinstance(result.data, dict):
        return ""
    text = str(result.data.get("raw_text") or "")
    if text:
        return text
    raw = result.data.get("content_bytes")
    if isinstance(raw, (bytes, bytearray)):
        return bytes(raw).decode("utf-8", errors="replace")
    return ""


def _apply_chart(cand: dict[str, Any], chart: dict[str, Any]) -> None:
    cand["graph_status"] = chart.get("graph_status") or GRAPH_STATUS_FAILED
    cand["graph_ok"] = cand["graph_status"] == GRAPH_STATUS_SUCCESS
    cand["graph_artifact"] = {k: v for k, v in chart.items() if k not in {"image_bytes", "tool_results"}}
    cand["acquisition_method"] = chart.get("acquisition_method")
    cand["graph_requested_url"] = chart.get("graph_requested_url")
    cand["graph_final_url"] = chart.get("graph_final_url")
    cand["graph_url_origin"] = chart.get("graph_url_origin")
    cand["validation_checks"] = dict(chart.get("validation_checks") or {})
    cand["image_width"] = chart.get("image_width")
    cand["image_height"] = chart.get("image_height")
    cand["graph_error_type"] = chart.get("error_type")
    if chart.get("image_bytes") and cand["graph_ok"]:
        cand["graph_image_bytes"] = chart["image_bytes"]
    else:
        cand.pop("graph_image_bytes", None)
    rescored = score_profile(
        cand,
        gds=cand.get("gds_metadata"),
        subset_effect_flag=bool(cand.get("subset_effect_flag")),
        graph_ok=cand["graph_ok"],
    )
    cand.update(rescored)


def collect_section_3a_profiles(
    gene_symbol: str,
    *,
    max_discovery_profiles: int = DEFAULT_MAX_DISCOVERY,
    max_selected_profiles: int = DEFAULT_MAX_SELECTED,
    max_chart_candidates: int | None = None,
    attempt_figures: bool = True,
    settings: Settings | None = None,
) -> dict[str, Any]:
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    chart_cap = max_chart_candidates if max_chart_candidates is not None else max(max_selected_profiles * 3, 15)
    exact_query = build_exact_gene_symbol_query(symbol)
    neural_query = build_neural_context_query(symbol)
    subset_query = build_subset_effect_query(symbol)

    exact = esearch_geoprofiles(gene_symbol=symbol, term=exact_query, retmax=0, sort=None, settings=cfg)
    neural_page = page_geoprofile_ids(gene_symbol=symbol, term=neural_query, max_ids=max_discovery_profiles, settings=cfg)
    subset_page = page_geoprofile_ids(gene_symbol=symbol, term=subset_query, max_ids=max_discovery_profiles, settings=cfg)
    neural_ids = list(neural_page["ids"])
    subset_ids = list(subset_page["ids"])
    subset_set = set(subset_ids)
    union_ids = list(dict.fromkeys([*neural_ids, *subset_ids]))

    search_status = "success"
    if not exact.success and not neural_ids and not any(r.success for r in neural_page["tool_results"]):
        search_status = "source_unavailable"

    profiles_by_uid: dict[str, dict[str, Any]] = {}
    tool_results: list[ToolResult] = [exact, *neural_page["tool_results"], *subset_page["tool_results"]]
    for start in range(0, len(union_ids), DEFAULT_ESUMMARY_BATCH):
        batch = union_ids[start:start + DEFAULT_ESUMMARY_BATCH]
        summary = esummary_geoprofiles(batch, gene_symbol=symbol, settings=cfg)
        tool_results.append(summary)
        if summary.success and isinstance(summary.data, dict):
            for doc in parse_profile_esummary_payload(summary.data):
                uid = str(doc.get("profile_uid") or "")
                if uid:
                    profiles_by_uid[uid] = doc

    gds_uids = sorted({str(d.get("gds_uid")) for d in profiles_by_uid.values() if d.get("gds_uid")})
    gds_by_uid: dict[str, dict[str, Any]] = {}
    for start in range(0, len(gds_uids), DEFAULT_ESUMMARY_BATCH):
        batch = gds_uids[start:start + DEFAULT_ESUMMARY_BATCH]
        gsum = esummary_gds(batch, gene_symbol=symbol, settings=cfg)
        tool_results.append(gsum)
        if gsum.success and isinstance(gsum.data, dict):
            gds_by_uid.update(parse_gds_esummary_payload(gsum.data))

    candidates: list[dict[str, Any]] = []
    for uid in union_ids:
        profile = dict(profiles_by_uid.get(uid) or {"profile_uid": uid, "uid": uid})
        gds_uid = str(profile.get("gds_uid") or "")
        gds = gds_by_uid.get(gds_uid)
        link_status = "not_attempted"
        if gds_uid:
            link = elink_profile_to_gds(uid, gene_symbol=symbol, settings=cfg)
            tool_results.append(link)
            if link.success:
                linked = extract_elink_gds_uids(link)
                link_status = "profile_gds_mismatch" if linked and normalize_gds_uid(linked[0]) != gds_uid else "validated"
            else:
                link_status = "unavailable"
        status, reasons = eligibility_status(profile, requested_gene=symbol, gds=gds)
        if link_status == "profile_gds_mismatch":
            status = "profile_gds_mismatch"
            reasons = [*reasons, "profile_gds_mismatch"]
        subset_flag = uid in subset_set
        scored = score_profile(profile, gds=gds, subset_effect_flag=subset_flag, graph_ok=None)
        candidates.append({
            **profile,
            "gds_metadata": gds,
            "subset_effect_flag": subset_flag,
            "eligibility_status": status,
            "rejection_reasons": reasons,
            "link_validation_status": link_status,
            "diversity_keys": _diversity_keys(profile, gds),
            "graph_status": GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE,
            "graph_ok": False,
            "pdat": (gds or {}).get("pdat"),
            **scored,
            "selected": False,
            "selection_rank": None,
            "in_chart_shortlist": False,
        })

    eligible = [c for c in candidates if c.get("eligibility_status") == "eligible"]
    eligible_sorted = sorted(eligible, key=_sort_key)
    by_uid = {str(c["profile_uid"]): c for c in candidates}

    def _fetch_and_apply_chart(
        live: dict[str, Any],
        *,
        browser_session: GeoProfilesBrowserSession | None,
    ) -> None:
        uid = str(live["profile_uid"])
        html_result = fetch_profile_html(uid, gene_symbol=symbol, settings=cfg)
        tool_results.append(html_result)
        chart = acquire_profile_chart(
            gene_symbol=symbol,
            profile_uid=uid,
            gds_uid=str(live.get("gds_uid") or ""),
            idref=str(live.get("idref") or ""),
            profile_html=_html_text(html_result),
            settings=cfg,
            browser_session=browser_session,
        )
        tool_results.extend(chart.get("tool_results") or [])
        live["in_chart_shortlist"] = True
        _apply_chart(live, chart)

    if attempt_figures:
        shortlist = build_diversity_shortlist(eligible_sorted, max_items=chart_cap)
        shortlist_uids = {str(c["profile_uid"]) for c in shortlist}
        for cand in candidates:
            if cand["profile_uid"] in shortlist_uids:
                cand["in_chart_shortlist"] = True
                cand["graph_status"] = "pending_shortlist"
            else:
                cand["graph_status"] = GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE

        def _organism_count(pool: Sequence[dict[str, Any]], organism: str | None) -> int:
            return sum(
                1
                for s in pool
                if (s.get("diversity_keys") or {}).get("organism") == organism
            )

        def _run_figure_acquisition(
            browser_session: GeoProfilesBrowserSession | None,
        ) -> list[dict[str, Any]]:
            for item in shortlist:
                _fetch_and_apply_chart(
                    by_uid[str(item["profile_uid"])],
                    browser_session=browser_session,
                )
            shortlist_pool_local = [by_uid[str(s["profile_uid"])] for s in shortlist]
            selected_local = build_diversity_shortlist(
                sorted(shortlist_pool_local, key=_sort_key),
                max_items=max_selected_profiles,
            )
            selected_uids = {str(s["profile_uid"]) for s in selected_local}
            selected_orgs = {
                (s.get("diversity_keys") or {}).get("organism") for s in selected_local
            }
            shortlist_uid_set = {str(s["profile_uid"]) for s in shortlist_pool_local}

            # Late diversity: acquire chart before selecting outside-shortlist picks.
            for cand in eligible_sorted:
                org = (cand.get("diversity_keys") or {}).get("organism")
                uid = str(cand["profile_uid"])
                if not org or org in selected_orgs or uid in selected_uids:
                    continue
                live = by_uid[uid]
                if live.get("graph_status") == GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE:
                    _fetch_and_apply_chart(live, browser_session=browser_session)
                if live.get("graph_status") == GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE:
                    continue
                if len(selected_local) < max_selected_profiles:
                    selected_local.append(live)
                    selected_uids.add(uid)
                    selected_orgs.add(org)
                    continue
                replaceable = sorted(
                    (
                        s
                        for s in selected_local
                        if _organism_count(
                            selected_local, (s.get("diversity_keys") or {}).get("organism")
                        )
                        > 1
                    ),
                    key=_sort_key,
                )
                if not replaceable:
                    continue
                victim = replaceable[-1]
                selected_local = [s for s in selected_local if s is not victim]
                selected_local.append(live)
                selected_uids = {str(s["profile_uid"]) for s in selected_local}
                selected_orgs = {
                    (s.get("diversity_keys") or {}).get("organism") for s in selected_local
                }
                _ = shortlist_uid_set
            return selected_local

        try:
            with GeoProfilesBrowserSession() as browser_session:
                selected = _run_figure_acquisition(browser_session)
        except Exception:  # noqa: BLE001
            selected = _run_figure_acquisition(None)
    else:
        for cand in candidates:
            cand["graph_status"] = GRAPH_STATUS_NOT_ATTEMPTED_OPTIONAL
            cand["graph_ok"] = False
        shortlist_pool = eligible_sorted
        selected = build_diversity_shortlist(
            sorted(shortlist_pool, key=_sort_key),
            max_items=max_selected_profiles,
        )

    selected_final: list[dict[str, Any]] = []
    for item in selected:
        if attempt_figures and item.get("graph_status") == GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE:
            # Hard gate: polished selection never keeps this audit-only status.
            continue
        selected_final.append(item)
    selected_final = selected_final[:max_selected_profiles]
    for rank, item in enumerate(selected_final, start=1):
        item["selected"] = True
        item["selection_rank"] = rank
    selected_uids = {str(s["profile_uid"]) for s in selected_final}
    for cand in candidates:
        if cand["profile_uid"] not in selected_uids:
            cand["selected"] = False
            cand["selection_rank"] = None

    if search_status != "success":
        scientific_status = "source_unavailable"
    elif not eligible:
        scientific_status = "no_relevant_profiles"
    else:
        scientific_status = "success"

    if not attempt_figures:
        visual_status = "not_attempted_optional"
    elif not selected_final:
        visual_status = "unavailable"
    else:
        ok_count = sum(1 for s in selected_final if s.get("graph_status") == GRAPH_STATUS_SUCCESS)
        if ok_count == len(selected_final):
            visual_status = "success"
        elif ok_count == 0:
            visual_status = "unavailable"
        else:
            visual_status = "partial"

    return {
        "gene_symbol": symbol,
        "exact_query": exact_query,
        "exact_profile_count": extract_count(exact),
        "exact_querytranslation": extract_querytranslation(exact),
        "neural_query": neural_query,
        "neural_profile_count": neural_page["count"],
        "neural_profile_ids": neural_ids,
        "neural_querytranslation": neural_page.get("querytranslation"),
        "subset_effect_query": subset_query,
        "subset_effect_profile_count": subset_page["count"],
        "subset_effect_profile_ids": subset_ids,
        "subset_effect_querytranslation": subset_page.get("querytranslation"),
        "candidate_union_ids": union_ids,
        "candidate_union_count": len(union_ids),
        "candidate_retrieval_truncated": bool(neural_page.get("truncated") or subset_page.get("truncated")),
        "max_discovery_profiles": max_discovery_profiles,
        "max_chart_candidates": chart_cap,
        "max_selected_profiles": max_selected_profiles,
        "attempt_figures": attempt_figures,
        "candidates": candidates,
        "selected_profiles": selected_final,
        "selected_profile_count": len(selected_final),
        "rejected_candidate_count": len(candidates) - len(selected_final),
        "search_status": search_status,
        "scientific_status": scientific_status,
        "visual_status": visual_status,
        "tool_results": tool_results,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
    }


__all__ = [
    "DEFAULT_MAX_DISCOVERY", "DEFAULT_MAX_SELECTED",
    "GRAPH_STATUS_FAILED", "GRAPH_STATUS_NOT_ATTEMPTED_OPTIONAL",
    "GRAPH_STATUS_NOT_ATTEMPTED_OUTSIDE", "GRAPH_STATUS_SUCCESS",
    "NCBI_BLOCK_TOKENS", "SOURCE_NAME",
    "GeoProfilesBrowserSession",
    "acquire_profile_chart", "build_diversity_shortlist",
    "build_exact_gene_symbol_query", "build_neural_context_query", "build_subset_effect_query",
    "collect_section_3a_profiles", "construct_graph_url", "detect_ncbi_block_page",
    "discover_graph_url_from_html",
    "eligibility_status", "elink_profile_to_gds", "esearch_geoprofiles", "esummary_gds",
    "esummary_geoprofiles", "expected_graph_id", "extract_count", "extract_elink_gds_uids",
    "extract_id_list", "fetch_profile_html", "format_gds_accession", "genes_match",
    "normalize_gds_uid", "page_geoprofile_ids", "parse_gds_esummary_payload",
    "parse_profile_esummary_payload", "score_profile", "validate_chart_bytes",
    "validate_graph_capture", "validate_graph_page_identity", "validate_graph_url",
]
