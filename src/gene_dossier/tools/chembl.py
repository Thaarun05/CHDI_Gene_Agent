"""ChEMBL client for chemical-tool / bioactivity annotations.

Searches targets and assays, then fetches activities for selected assay IDs.
Does **not** normalize into evidence records — that belongs in
``normalize/chemicals.py``.

Key endpoints (validated)::

    GET https://www.ebi.ac.uk/chembl/api/data/target/search.json?q={symbol}
    GET https://www.ebi.ac.uk/chembl/api/data/assay.json
        ?description__icontains={term}&limit=100
    GET https://www.ebi.ac.uk/chembl/api/data/activity.json
        ?assay_chembl_id__in={ids}&limit=1000

NOTE: Do not append a trailing ``?`` after ``limit=`` (invalid).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "ChEMBL"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

DEFAULT_ASSAY_LIMIT = 100
DEFAULT_ACTIVITY_LIMIT = 1000
ORGANISM_HUMAN = "Homo sapiens"

# Gene-specific assay description aliases (do not apply globally to every gene).
GENE_SPECIFIC_ASSAY_TERMS: dict[str, list[str]] = {
    "SREBF2": ["SREBP2", "sterol regulatory element-binding protein"],
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


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    path: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a ChEMBL JSON path and return :class:`ToolResult`."""
    url = f"{CHEMBL_BASE}/{path.lstrip('/')}"
    # Ensure params are strings; never leave a dangling '?' in the URL.
    query = {k: str(v) for k, v in params.items() if v is not None}
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


def summarize_target(target: dict[str, Any]) -> dict[str, Any]:
    """Extract key ChEMBL target fields (not evidence)."""
    components = target.get("target_components") or []
    accessions: list[str] = []
    synonyms: list[Any] = []
    if isinstance(components, list):
        for comp in components:
            if not isinstance(comp, dict):
                continue
            acc = comp.get("accession")
            if acc:
                accessions.append(str(acc))
            syns = comp.get("target_component_synonyms")
            if isinstance(syns, list):
                synonyms.extend(syns)
    return {
        "target_chembl_id": target.get("target_chembl_id"),
        "pref_name": target.get("pref_name"),
        "organism": target.get("organism"),
        "target_type": target.get("target_type"),
        "accessions": accessions,
        "synonyms": synonyms,
    }


def summarize_assay(assay: dict[str, Any]) -> dict[str, Any]:
    """Extract key ChEMBL assay fields (not evidence)."""
    return {
        "assay_chembl_id": assay.get("assay_chembl_id"),
        "description": assay.get("description"),
        "assay_type": assay.get("assay_type"),
        "assay_organism": assay.get("assay_organism"),
        "assay_cell_type": assay.get("assay_cell_type"),
        "target_chembl_id": assay.get("target_chembl_id"),
        "relationship_type": assay.get("relationship_type"),
        "confidence_score": assay.get("confidence_score"),
        "confidence_description": assay.get("confidence_description"),
        "document_chembl_id": assay.get("document_chembl_id"),
    }


def summarize_activity(activity: dict[str, Any]) -> dict[str, Any]:
    """Extract key ChEMBL activity fields (not evidence)."""
    return {
        "molecule_chembl_id": activity.get("molecule_chembl_id"),
        "canonical_smiles": activity.get("canonical_smiles"),
        "standard_type": activity.get("standard_type"),
        "standard_relation": activity.get("standard_relation"),
        "standard_value": activity.get("standard_value"),
        "standard_units": activity.get("standard_units"),
        "pchembl_value": activity.get("pchembl_value"),
        "assay_chembl_id": activity.get("assay_chembl_id"),
        "document_chembl_id": activity.get("document_chembl_id"),
    }


def select_target_id(
    targets: list[Any],
    gene_symbol: str,
    *,
    organism: str = ORGANISM_HUMAN,
) -> tuple[str | None, str]:
    """Select a ChEMBL target ID by name/synonym match.

    Returns ``(target_chembl_id | None, selection_method)`` where method is:
    - ``"matched"`` — exactly one best name/synonym match
    - ``"ambiguous"`` — multiple equally best name/synonym matches
    - ``"not_found"`` — no name/synonym match

    Does **not** fall back to the first human/any target.
    """
    target = gene_symbol.strip().upper()
    if not target:
        return None, "not_found"

    scored: list[tuple[int, str]] = []
    for row in targets:
        if not isinstance(row, dict):
            continue
        tid = row.get("target_chembl_id")
        if not tid:
            continue
        pref = str(row.get("pref_name") or "").upper()
        org = str(row.get("organism") or "")
        syn_hit = False
        for comp in row.get("target_components") or []:
            if not isinstance(comp, dict):
                continue
            for syn in comp.get("target_component_synonyms") or []:
                # Synonym entries may be dicts with ``component_synonym`` or plain strings.
                if isinstance(syn, dict):
                    syn_text = str(
                        syn.get("component_synonym") or syn.get("synonym") or ""
                    ).upper()
                else:
                    syn_text = str(syn).upper()
                if syn_text == target or target in syn_text.split():
                    syn_hit = True
                    break
            if syn_hit:
                break

        name_hit = target == pref or target in pref.split() or target in pref
        if not (name_hit or syn_hit):
            continue
        # Lower rank = better: exact human match first.
        rank = 0
        if org != organism:
            rank += 2
        if not name_hit:
            rank += 1
        scored.append((rank, str(tid)))

    if not scored:
        return None, "not_found"

    scored.sort(key=lambda item: item[0])
    best_rank = scored[0][0]
    best_ids = [tid for rank, tid in scored if rank == best_rank]
    # Deduplicate while preserving order.
    unique_best: list[str] = []
    for tid in best_ids:
        if tid not in unique_best:
            unique_best.append(tid)
    if len(unique_best) > 1:
        return None, "ambiguous"
    return unique_best[0], "matched"


def prefer_target_id(
    targets: list[Any],
    gene_symbol: str,
    *,
    organism: str = ORGANISM_HUMAN,
) -> str | None:
    """Return a matched target ID, or ``None`` when not found / ambiguous.

    See :func:`select_target_id`. Never falls back to the first candidate.
    """
    selected, _method = select_target_id(
        targets, gene_symbol, organism=organism
    )
    return selected


def resolve_authoritative_target(
    targets: list[Any],
    *,
    uniprot_accession: str,
    gene_symbol: str,
    organism: str = ORGANISM_HUMAN,
) -> tuple[str | None, str, str]:
    """Resolve a human ChEMBL target by authoritative UniProt accession.

    Prefers organism-matching ``SINGLE PROTEIN`` rows whose
    ``target_components[].accession`` equals ``uniprot_accession``.

    Returns ``(target_chembl_id | None, method, detail)`` where method is one of:
    - ``uniprot_single_protein`` — unique human SINGLE PROTEIN UniProt match
    - ``uniprot_match`` — unique UniProt match that is not SINGLE PROTEIN
    - ``ambiguous`` — multiple equally preferred matches
    - ``not_found`` — no UniProt accession match

    Never falls back to name/synonym scoring or the first search hit.
    """
    accession = str(uniprot_accession or "").strip().upper()
    symbol = str(gene_symbol or "").strip()
    if not accession:
        return None, "not_found", "missing uniprot_accession"

    single_protein: list[str] = []
    other_matches: list[str] = []
    for row in targets:
        if not isinstance(row, dict):
            continue
        tid = row.get("target_chembl_id")
        if not tid:
            continue
        org = str(row.get("organism") or "")
        if organism and org and org != organism:
            continue
        components = row.get("target_components") or []
        if not isinstance(components, list):
            continue
        has_accession = False
        for comp in components:
            if not isinstance(comp, dict):
                continue
            acc = str(comp.get("accession") or "").strip().upper()
            if acc == accession:
                has_accession = True
                break
        if not has_accession:
            continue
        target_type = str(row.get("target_type") or "").strip().upper()
        tid_s = str(tid)
        if target_type == "SINGLE PROTEIN":
            if tid_s not in single_protein:
                single_protein.append(tid_s)
        else:
            if tid_s not in other_matches:
                other_matches.append(tid_s)

    if len(single_protein) == 1:
        return (
            single_protein[0],
            "uniprot_single_protein",
            f"{accession} SINGLE PROTEIN ({symbol or 'gene'})",
        )
    if len(single_protein) > 1:
        return (
            None,
            "ambiguous",
            f"multiple SINGLE PROTEIN matches for {accession}: {single_protein}",
        )
    if len(other_matches) == 1:
        return (
            other_matches[0],
            "uniprot_match",
            f"{accession} non-SINGLE-PROTEIN match ({symbol or 'gene'})",
        )
    if len(other_matches) > 1:
        return (
            None,
            "ambiguous",
            f"multiple UniProt matches for {accession}: {other_matches}",
        )
    return None, "not_found", f"no human UniProt match for {accession}"


def default_assay_search_terms(
    gene_symbol: str,
    aliases: list[str] | None = None,
) -> list[str]:
    """Return assay ``description__icontains`` terms for ``gene_symbol``.

    Includes:
    - ``gene_symbol``
    - ``aliases`` passed by the workflow/normalizer (when provided)
    - gene-specific mapped terms only when ``gene_symbol`` is in
      :data:`GENE_SPECIFIC_ASSAY_TERMS`

    Example: ``default_assay_search_terms("HTT")`` → ``["HTT"]``.
    """
    symbol = gene_symbol.strip()
    terms: list[str] = []

    def _add(term: str) -> None:
        cleaned = term.strip()
        if not cleaned:
            return
        for existing in terms:
            if existing.upper() == cleaned.upper():
                return
        terms.append(cleaned)

    _add(symbol)
    if aliases:
        for alias in aliases:
            _add(str(alias))
    for mapped in GENE_SPECIFIC_ASSAY_TERMS.get(symbol.upper(), []):
        _add(mapped)
    # Also allow exact-case map keys.
    if symbol.upper() != symbol:
        for mapped in GENE_SPECIFIC_ASSAY_TERMS.get(symbol, []):
            _add(mapped)
    return terms


def target_search(
    gene_symbol: str,
    *,
    settings: Settings | None = None,
) -> ToolResult:
    """Search ChEMBL targets for ``gene_symbol``."""
    cfg = settings or get_settings()
    params = {"q": gene_symbol.strip()}
    return _request_json(
        endpoint_name="target_search",
        gene_symbol=gene_symbol,
        path="target/search.json",
        params=params,
        settings=cfg,
    )


def assay_search(
    term: str,
    *,
    gene_symbol: str = "",
    limit: int = DEFAULT_ASSAY_LIMIT,
    settings: Settings | None = None,
) -> ToolResult:
    """Search ChEMBL assays by description substring."""
    cfg = settings or get_settings()
    params = {
        "description__icontains": term,
        "limit": limit,
    }
    return _request_json(
        endpoint_name="assay_search",
        gene_symbol=gene_symbol or term,
        path="assay.json",
        params=params,
        settings=cfg,
    )


def activities(
    assay_chembl_ids: str | list[str],
    *,
    gene_symbol: str = "",
    limit: int = DEFAULT_ACTIVITY_LIMIT,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch ChEMBL activities for one or more assay IDs (``assay_chembl_id__in``)."""
    cfg = settings or get_settings()
    if isinstance(assay_chembl_ids, (list, tuple)):
        id_str = ",".join(str(a).strip() for a in assay_chembl_ids if str(a).strip())
    else:
        id_str = str(assay_chembl_ids).strip()
    if not id_str:
        return _tool_result(
            endpoint_name="activities",
            gene_symbol=gene_symbol,
            request_url=f"{CHEMBL_BASE}/activity.json",
            request_params={"assay_chembl_id__in": "", "limit": str(limit)},
            success=False,
            error_type="invalid_request",
            error_message="ChEMBL activities require at least one assay_chembl_id",
        )
    params = {
        "assay_chembl_id__in": id_str,
        "limit": limit,
    }
    return _request_json(
        endpoint_name="activities",
        gene_symbol=gene_symbol or id_str,
        path="activity.json",
        params=params,
        settings=cfg,
    )


def activities_by_target(
    target_chembl_id: str,
    gene_symbol: str = "",
    limit: int = DEFAULT_ACTIVITY_LIMIT,
    offset: int = 0,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch activities filtered by exact ``target_chembl_id``."""
    cfg = settings or get_settings()
    tid = str(target_chembl_id).strip()
    if not tid:
        return _tool_result(
            endpoint_name="activities_by_target",
            gene_symbol=gene_symbol,
            request_url=f"{CHEMBL_BASE}/activity.json",
            request_params={"target_chembl_id": "", "limit": str(limit), "offset": str(offset)},
            success=False,
            error_type="invalid_request",
            error_message="ChEMBL activities_by_target requires target_chembl_id",
        )
    params = {
        "target_chembl_id": tid,
        "limit": limit,
        "offset": offset,
    }
    return _request_json(
        endpoint_name="activities_by_target",
        gene_symbol=gene_symbol or tid,
        path="activity.json",
        params=params,
        settings=cfg,
    )


def assays_by_target(
    target_chembl_id: str,
    gene_symbol: str = "",
    limit: int = DEFAULT_ASSAY_LIMIT,
    offset: int = 0,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch assays filtered by exact ``target_chembl_id`` (one HTTP page)."""
    cfg = settings or get_settings()
    tid = str(target_chembl_id).strip()
    if not tid:
        return _tool_result(
            endpoint_name="assays_by_target",
            gene_symbol=gene_symbol,
            request_url=f"{CHEMBL_BASE}/assay.json",
            request_params={"target_chembl_id": "", "limit": str(limit), "offset": str(offset)},
            success=False,
            error_type="invalid_request",
            error_message="ChEMBL assays_by_target requires target_chembl_id",
        )
    params = {
        "target_chembl_id": tid,
        "limit": limit,
        "offset": offset,
    }
    return _request_json(
        endpoint_name="assays_by_target",
        gene_symbol=gene_symbol or tid,
        path="assay.json",
        params=params,
        settings=cfg,
    )


def is_direct_assay_relationship(assay: dict[str, Any] | None) -> bool:
    """True only when ChEMBL assay metadata explicitly marks a direct target link."""
    if not isinstance(assay, dict):
        return False
    rel = str(assay.get("relationship_type") or "").strip().upper()
    return rel in {"D", "DIRECT"}


def fetch_exact_target_activities(
    target_chembl_id: str,
    *,
    gene_symbol: str = "",
    limit: int = DEFAULT_ACTIVITY_LIMIT,
    max_pages: int = 50,
    settings: Settings | None = None,
) -> ToolResult:
    """Paginate exact-target ``activity.json?target_chembl_id=`` until exhausted.

    On success, ``data`` includes merged ``activities``, ``activity_summaries``,
    page count, and the last page payload metadata.
    """
    cfg = settings or get_settings()
    tid = str(target_chembl_id).strip()
    if not tid:
        return _tool_result(
            endpoint_name="fetch_exact_target_activities",
            gene_symbol=gene_symbol,
            request_url=f"{CHEMBL_BASE}/activity.json",
            request_params={"target_chembl_id": ""},
            success=False,
            error_type="invalid_request",
            error_message="fetch_exact_target_activities requires target_chembl_id",
        )

    all_activities: list[dict[str, Any]] = []
    pages: list[Any] = []
    offset = 0
    last_url = f"{CHEMBL_BASE}/activity.json"
    last_params: dict[str, Any] = {}
    last_status: int | None = None
    total_count: int | None = None
    page_idx = 0

    while page_idx < max(1, int(max_pages)):
        page = activities_by_target(
            tid,
            gene_symbol=gene_symbol,
            limit=limit,
            offset=offset,
            settings=cfg,
        )
        last_url = page.request_url
        last_params = page.request_params
        last_status = page.status_code
        if not page.success:
            return _tool_result(
                endpoint_name="fetch_exact_target_activities",
                gene_symbol=gene_symbol or tid,
                request_url=page.request_url,
                request_params=page.request_params,
                success=False,
                status_code=page.status_code,
                data={
                    "target_chembl_id": tid,
                    "activities": all_activities,
                    "pages": pages,
                    "failed_offset": offset,
                    "page_error": page.data,
                },
                error_type=page.error_type or "activities_by_target_failed",
                error_message=page.error_message or "ChEMBL exact-target activities failed",
            )
        pages.append(page.data)
        payload = page.data if isinstance(page.data, dict) else {}
        page_meta = payload.get("page_meta") if isinstance(payload.get("page_meta"), dict) else {}
        if total_count is None and page_meta.get("total_count") is not None:
            try:
                total_count = int(page_meta["total_count"])
            except (TypeError, ValueError):
                total_count = None
        batch = payload.get("activities") or []
        if not isinstance(batch, list):
            batch = []
        for row in batch:
            if isinstance(row, dict):
                all_activities.append(row)
        if not batch:
            break
        offset += len(batch)
        if total_count is not None and offset >= total_count:
            break
        if page_meta.get("next") in (None, "", False):
            break
        page_idx += 1

    summaries = [summarize_activity(a) for a in all_activities]
    return _tool_result(
        endpoint_name="fetch_exact_target_activities",
        gene_symbol=gene_symbol or tid,
        request_url=last_url,
        request_params={
            "target_chembl_id": tid,
            "limit": limit,
            "pages_fetched": len(pages),
            **last_params,
        },
        success=True,
        status_code=last_status,
        data={
            "target_chembl_id": tid,
            "activities": all_activities,
            "activity_summaries": summaries,
            "activity_count": len(all_activities),
            "total_count": total_count if total_count is not None else len(all_activities),
            "pages_fetched": len(pages),
            "page_payloads": pages,
        },
    )


def fetch_chemical_tools(
    gene_symbol: str,
    *,
    assay_terms: list[str] | None = None,
    aliases: list[str] | None = None,
    assay_limit: int = DEFAULT_ASSAY_LIMIT,
    activity_limit: int = DEFAULT_ACTIVITY_LIMIT,
    settings: Settings | None = None,
) -> ToolResult:
    """Target search → assay searches → activities for discovered assay IDs.

    On success, ``data`` includes preferred ``target_chembl_id``,
    ``target_selection_method`` (``matched`` / ``not_found`` / ``ambiguous``),
    target/assay/activity summaries, and raw payloads.

    When ``assay_terms`` is omitted, terms come from
    :func:`default_assay_search_terms` (symbol + optional ``aliases`` +
    gene-specific map only).

    Never raises.
    """
    cfg = settings or get_settings()
    targets_res = target_search(gene_symbol, settings=cfg)
    if not targets_res.success:
        return _tool_result(
            endpoint_name="fetch_chemical_tools",
            gene_symbol=gene_symbol,
            request_url=targets_res.request_url,
            request_params=targets_res.request_params,
            success=False,
            status_code=targets_res.status_code,
            data={"target_search": targets_res.data},
            error_type=targets_res.error_type or "target_search_failed",
            error_message=targets_res.error_message or "ChEMBL target search failed",
        )

    target_payload = targets_res.data if isinstance(targets_res.data, dict) else {}
    targets = target_payload.get("targets") or []
    if not isinstance(targets, list):
        targets = []
    target_summaries = [
        summarize_target(t) for t in targets if isinstance(t, dict)
    ]
    preferred_target, target_selection_method = select_target_id(
        targets, gene_symbol
    )

    terms = (
        assay_terms
        if assay_terms is not None
        else default_assay_search_terms(gene_symbol, aliases=aliases)
    )
    assay_raw_by_term: dict[str, Any] = {}
    assay_by_id: dict[str, dict[str, Any]] = {}
    last_assay_url = targets_res.request_url
    last_assay_params: dict[str, Any] = {}
    last_status = targets_res.status_code

    for term in terms:
        term = str(term).strip()
        if not term:
            continue
        ares = assay_search(
            term, gene_symbol=gene_symbol, limit=assay_limit, settings=cfg
        )
        last_assay_url = ares.request_url
        last_assay_params = ares.request_params
        last_status = ares.status_code
        assay_raw_by_term[term] = ares.data
        if not ares.success:
            return _tool_result(
                endpoint_name="fetch_chemical_tools",
                gene_symbol=gene_symbol,
                request_url=ares.request_url,
                request_params=ares.request_params,
                success=False,
                status_code=ares.status_code,
                data={
                    "target_search": targets_res.data,
                    "target_chembl_id": preferred_target,
                    "target_selection_method": target_selection_method,
                    "target_summaries": target_summaries,
                    "assay_terms": terms,
                    "assay_searches": assay_raw_by_term,
                    "failed_assay_term": term,
                },
                error_type=ares.error_type or "assay_search_failed",
                error_message=ares.error_message
                or f"ChEMBL assay search failed for term {term!r}",
            )
        assays = (ares.data or {}).get("assays") if isinstance(ares.data, dict) else None
        if isinstance(assays, list):
            for assay in assays:
                if not isinstance(assay, dict):
                    continue
                aid = assay.get("assay_chembl_id")
                if aid:
                    assay_by_id[str(aid)] = assay

    assay_ids = sorted(assay_by_id.keys())
    assay_summaries = [summarize_assay(assay_by_id[aid]) for aid in assay_ids]

    activity_payload: Any = None
    activity_summaries: list[dict[str, Any]] = []
    if assay_ids:
        act = activities(
            assay_ids,
            gene_symbol=gene_symbol,
            limit=activity_limit,
            settings=cfg,
        )
        last_assay_url = act.request_url
        last_assay_params = act.request_params
        last_status = act.status_code
        if not act.success:
            return _tool_result(
                endpoint_name="fetch_chemical_tools",
                gene_symbol=gene_symbol,
                request_url=act.request_url,
                request_params=act.request_params,
                success=False,
                status_code=act.status_code,
                data={
                    "target_search": targets_res.data,
                    "target_chembl_id": preferred_target,
                    "target_selection_method": target_selection_method,
                    "target_summaries": target_summaries,
                    "assay_terms": terms,
                    "assay_searches": assay_raw_by_term,
                    "assay_chembl_ids": assay_ids,
                    "assay_summaries": assay_summaries,
                    "activities": act.data,
                },
                error_type=act.error_type or "activities_failed",
                error_message=act.error_message or "ChEMBL activities failed",
            )
        activity_payload = act.data
        acts = (
            activity_payload.get("activities")
            if isinstance(activity_payload, dict)
            else None
        )
        if isinstance(acts, list):
            activity_summaries = [
                summarize_activity(a) for a in acts if isinstance(a, dict)
            ]

    return _tool_result(
        endpoint_name="fetch_chemical_tools",
        gene_symbol=gene_symbol,
        request_url=last_assay_url,
        request_params={
            "gene_symbol": gene_symbol,
            "target_chembl_id": preferred_target,
            "target_selection_method": target_selection_method,
            "assay_terms": terms,
            "assay_chembl_ids": assay_ids,
            **last_assay_params,
        },
        success=True,
        status_code=last_status,
        data={
            "gene_symbol": gene_symbol,
            "target_chembl_id": preferred_target,
            "target_selection_method": target_selection_method,
            "target_search": targets_res.data,
            "target_summaries": target_summaries,
            "assay_terms": terms,
            "assay_searches": assay_raw_by_term,
            "assay_chembl_ids": assay_ids,
            "assay_summaries": assay_summaries,
            "assay_count": len(assay_ids),
            "activities": activity_payload,
            "activity_summaries": activity_summaries,
            "activity_count": len(activity_summaries),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "CHEMBL_BASE",
    "DEFAULT_ASSAY_LIMIT",
    "DEFAULT_ACTIVITY_LIMIT",
    "ORGANISM_HUMAN",
    "GENE_SPECIFIC_ASSAY_TERMS",
    "summarize_target",
    "summarize_assay",
    "summarize_activity",
    "select_target_id",
    "prefer_target_id",
    "resolve_authoritative_target",
    "default_assay_search_terms",
    "target_search",
    "assay_search",
    "activities",
    "activities_by_target",
    "assays_by_target",
    "is_direct_assay_relationship",
    "fetch_exact_target_activities",
    "fetch_chemical_tools",
]
