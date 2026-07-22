"""NCBI Conserved Domain Database (CDD) client via Batch CD-Search (bwrpsb).

Submits a protein query, polls until complete, then retrieves domain hits.
Does **not** normalize into evidence records — that belongs in
``normalize/protein.py``.

Key endpoints (validated)::

    GET https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
        ?queries={refseq}&db=cdd&smode=auto&useid1=true&maxhit=250
        &filter=true&evalue=0.01&tdata=hits
    GET .../bwrpsb.cgi?cdsid={cdsid}&tdata=hits&dmode=full&qdefl=true&cddefl=true
    GET .../bwrpsb.cgi?cdsid={cdsid}&tdata=aligns&alnfmt=json   (optional)

NOTE: Batch CDD returns text/status for submit/poll when ``tdata=hits`` is set;
without ``tdata`` NCBI often returns HTML UI. Poll until ``status=0`` (completed).
``alnfmt=json`` only applies when ``tdata=aligns``.

SREBF2 RefSeq protein example: ``NP_004590.2`` (domains include ``cd18922`` /
``bHLHzip_SREBP2``).

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import csv
import io
import re
import time
from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "CDD"
BWRPSB_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"

DEFAULT_DB = "cdd"
DEFAULT_MAXHIT = 250
DEFAULT_EVALUE = 0.01
DEFAULT_QUERY_SREBF2 = "NP_004590.2"

# Polling defaults for Batch CD-Search job completion.
DEFAULT_POLL_INTERVAL_SECONDS = 2.0
DEFAULT_MAX_POLLS = 30

# status=0 means search completed successfully (NCBI Batch CD-Search).
STATUS_COMPLETED = 0

# Programmatic Batch CD-Search returns QM3-qcdsearch-* Search-IDs.
_CDSID_RE = re.compile(r"(QM3-qcdsearch-[A-Za-z0-9-]+)")
_HTML_ID_CDSID_RE = re.compile(
    r'id=["\']id_cdsid["\'][^>]*>([^<]+)<',
    re.IGNORECASE,
)
_HTML_CTRL_HANDLE_RE = re.compile(
    r'ctrlHandle\s*=\s*["\']([^"\']+)["\']',
    re.IGNORECASE,
)
_HTML_NAME_CDSID_RE = re.compile(
    r'name=["\']cdsid["\'][^>]*value=["\']([^"\']+)["\']',
    re.IGNORECASE,
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


def _raw_preview(raw_text: str, *, limit: int = 200) -> str:
    """Whitespace-collapsed short preview for soft-fail error messages."""
    collapsed = re.sub(r"\s+", " ", (raw_text or "")).strip()
    if len(collapsed) <= limit:
        return collapsed
    return collapsed[: limit - 3] + "..."


def _extract_cdsid_candidate(raw: str) -> str | None:
    """Return a QM3-qcdsearch-* token from ``raw`` when present."""
    match = _CDSID_RE.search(raw or "")
    if not match:
        return None
    return match.group(1)


def _extract_cdsid_from_html(raw_text: str) -> str | None:
    """Extract cdsid from NCBI Batch CD-Search HTML UI fallbacks."""
    text = raw_text or ""
    for pattern in (_HTML_ID_CDSID_RE, _HTML_CTRL_HANDLE_RE, _HTML_NAME_CDSID_RE):
        match = pattern.search(text)
        if not match:
            continue
        candidate = _extract_cdsid_candidate(match.group(1).strip())
        if candidate:
            return candidate
    # Last-resort scan of the whole document for a QM3 token (never invents).
    return _extract_cdsid_candidate(text)


def parse_status_text(raw_text: str) -> dict[str, Any]:
    """Parse Batch CD-Search status/text (or HTML UI) for ``cdsid`` and ``status``.

    Text formats::

        #cdsid\\tQM3-qcdsearch-...
        #status\\t3\\tmsg\\tJob is still running
        #status\\t0

    When a completed job also emits a warning status line, keep the numeric
    status and record every status value in ``status_lines``::

        #status\\t0
        #status\\tWarning: Too many queries.\\tmsg\\t...

    HTML fallbacks (when NCBI returns the UI instead of text/plain)::

        <td id="id_cdsid" ...>QM3-qcdsearch-...</td>
        var ctrlHandle = "QM3-qcdsearch-...";
        <input name="cdsid" value="QM3-qcdsearch-...">
    """
    cdsid: str | None = None
    status: int | None = None
    extras: dict[str, str] = {}
    status_lines: list[str] = []
    for line in (raw_text or "").splitlines():
        line = line.strip()
        if not line.startswith("#"):
            continue
        # Formats seen: "#cdsid\tVALUE" or "#cdsid = VALUE"
        body = line.lstrip("#").strip()
        key = ""
        value = ""
        if "\t" in body:
            key, value = body.split("\t", 1)
        elif "=" in body:
            key, value = body.split("=", 1)
        else:
            parts = body.split(None, 1)
            if len(parts) == 2:
                key, value = parts
            else:
                continue
        key = key.strip().lower()
        value = value.strip()
        extras[key] = value
        if key == "cdsid":
            cdsid = _extract_cdsid_candidate(value) or (value or None)
        elif key == "status":
            # Preserve every status line (numeric completion + later warnings).
            status_lines.append(value)
            # Multi-field: "#status\t3\tmsg\tJob is still running"
            token = value.split("\t", 1)[0].strip().split(None, 1)[0] if value else ""
            if token.isdigit() or (token.startswith("-") and token[1:].isdigit()):
                status = int(token)
            # Non-numeric warning lines must not overwrite an existing numeric status.

    if not cdsid:
        cdsid = _extract_cdsid_from_html(raw_text or "")

    return {
        "cdsid": cdsid,
        "status": status,
        "fields": extras,
        "status_lines": status_lines,
    }


def parse_hits_text(raw_text: str) -> list[dict[str, str]]:
    """Parse completed CDD hits tabular text into row dicts.

    Header detection is best-effort across NCBI format variants. Raw text should
    still be preserved for artifacts.
    """
    lines = [
        line
        for line in (raw_text or "").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    if not lines:
        return []

    # Prefer a tab-delimited table with a recognizable header.
    header_idx = None
    for i, line in enumerate(lines):
        lower = line.lower()
        if "\t" in line and (
            "accession" in lower
            or "short name" in lower
            or "e-value" in lower
            or "evalue" in lower
        ):
            header_idx = i
            break

    if header_idx is None:
        # Fallback: treat first non-comment line as header if tabbed.
        if "\t" in lines[0]:
            header_idx = 0
        else:
            return []

    table_text = "\n".join(lines[header_idx:])
    reader = csv.DictReader(io.StringIO(table_text), delimiter="\t")
    rows: list[dict[str, str]] = []
    for row in reader:
        if not isinstance(row, dict):
            continue
        cleaned = {
            str(k).strip(): ("" if v is None else str(v).strip())
            for k, v in row.items()
            if k is not None
        }
        if any(v for v in cleaned.values()):
            rows.append(cleaned)
    return rows


def summarize_hit(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key CDD hit fields across common column name variants."""

    def _get(*names: str) -> Any:
        lower_map = {str(k).strip().lower(): v for k, v in row.items()}
        for name in names:
            if name.lower() in lower_map:
                return lower_map[name.lower()]
        return None

    return {
        "query_accession": _get(
            "query", "query accession", "query id", "queryacc"
        ),
        "domain_accession": _get(
            "accession", "domain accession", "hit accession", "pssm-id", "pssm_id"
        ),
        "domain_short_name": _get("short name", "short_name", "domain short name"),
        "domain_description": _get(
            "description", "domain description", "definition", "cddefl"
        ),
        "from_residue": _get("from", "from residue", "start", "query start"),
        "to_residue": _get("to", "to residue", "end", "stop", "query stop"),
        "evalue": _get("e-value", "evalue", "e_value"),
        "bitscore": _get("bitscore", "bit score", "bit_score"),
        "superfamily": _get(
            "superfamily", "superfamily accession", "hit type"
        ),
        "raw": dict(row),
    }


def _get_text(
    *,
    endpoint_name: str,
    gene_symbol: str,
    params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET bwrpsb.cgi and return text payload in :class:`ToolResult`."""
    query = {k: str(v) for k, v in params.items() if v is not None}
    request_url = f"{BWRPSB_URL}?{urlencode(query)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(BWRPSB_URL, params=query)
        text = response.text
        payload = {"raw_text": text, "content_type": response.headers.get("content-type")}
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


def submit_search(
    queries: str,
    *,
    gene_symbol: str = "",
    db: str = DEFAULT_DB,
    maxhit: int = DEFAULT_MAXHIT,
    evalue: float = DEFAULT_EVALUE,
    settings: Settings | None = None,
) -> ToolResult:
    """Submit a Batch CD-Search job; return ``cdsid`` / ``status`` when parseable."""
    cfg = settings or get_settings()
    query_id = queries.strip()
    if not query_id:
        return _tool_result(
            endpoint_name="submit_search",
            gene_symbol=gene_symbol,
            request_url=BWRPSB_URL,
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="CDD submit requires a non-empty queries value",
        )
    params = {
        "queries": query_id,
        "db": db,
        "smode": "auto",
        "useid1": "true",
        "maxhit": maxhit,
        "filter": "true",
        "evalue": evalue,
        # Without tdata, NCBI often returns the HTML UI instead of #cdsid text.
        "tdata": "hits",
    }
    result = _get_text(
        endpoint_name="submit_search",
        gene_symbol=gene_symbol or query_id,
        params=params,
        settings=cfg,
    )
    if not result.success:
        return result

    raw_text = ""
    if isinstance(result.data, dict):
        raw_text = str(result.data.get("raw_text") or "")
    parsed = parse_status_text(raw_text)
    if not parsed.get("cdsid"):
        preview = _raw_preview(raw_text)
        return _tool_result(
            endpoint_name="submit_search",
            gene_symbol=gene_symbol or query_id,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"raw_text": raw_text, "parsed": parsed},
            error_type="parse_error",
            error_message=(
                "Could not parse cdsid from CDD submit response"
                + (f" (preview: {preview})" if preview else "")
            ),
        )
    return _tool_result(
        endpoint_name="submit_search",
        gene_symbol=gene_symbol or query_id,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "queries": query_id,
            "cdsid": parsed["cdsid"],
            "status": parsed["status"],
            "parsed": parsed,
            "raw_text": raw_text,
        },
    )


def poll_status(
    cdsid: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Poll Batch CD-Search status for an existing ``cdsid``."""
    cfg = settings or get_settings()
    sid = cdsid.strip()
    params = {"cdsid": sid, "tdata": "hits"}
    result = _get_text(
        endpoint_name="poll_status",
        gene_symbol=gene_symbol or sid,
        params=params,
        settings=cfg,
    )
    if not result.success:
        return result
    raw_text = ""
    if isinstance(result.data, dict):
        raw_text = str(result.data.get("raw_text") or "")
    parsed = parse_status_text(raw_text)
    return _tool_result(
        endpoint_name="poll_status",
        gene_symbol=gene_symbol or sid,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "cdsid": parsed.get("cdsid") or sid,
            "status": parsed.get("status"),
            "parsed": parsed,
            "raw_text": raw_text,
        },
    )


def retrieve_hits(
    cdsid: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Retrieve conserved domain hits for a completed ``cdsid``."""
    cfg = settings or get_settings()
    sid = cdsid.strip()
    params = {
        "cdsid": sid,
        "tdata": "hits",
        "dmode": "full",
        "qdefl": "true",
        "cddefl": "true",
    }
    result = _get_text(
        endpoint_name="retrieve_hits",
        gene_symbol=gene_symbol or sid,
        params=params,
        settings=cfg,
    )
    if not result.success:
        return result
    raw_text = ""
    if isinstance(result.data, dict):
        raw_text = str(result.data.get("raw_text") or "")
    rows = parse_hits_text(raw_text)
    return _tool_result(
        endpoint_name="retrieve_hits",
        gene_symbol=gene_symbol or sid,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "cdsid": sid,
            "raw_text": raw_text,
            "hit_rows": rows,
            "hit_summaries": [summarize_hit(r) for r in rows],
            "hit_count": len(rows),
        },
    )


def retrieve_aligns(
    cdsid: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Optional residue-level alignments (``tdata=aligns&alnfmt=json``)."""
    cfg = settings or get_settings()
    sid = cdsid.strip()
    params = {
        "cdsid": sid,
        "tdata": "aligns",
        "alnfmt": "json",
    }
    query = {k: str(v) for k, v in params.items()}
    request_url = f"{BWRPSB_URL}?{urlencode(query)}"
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(BWRPSB_URL, params=query)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {
                "raw_text": response.text[:4000],
                "content_type": response.headers.get("content-type"),
            }
        if response.is_success:
            return _tool_result(
                endpoint_name="retrieve_aligns",
                gene_symbol=gene_symbol or sid,
                request_url=request_url,
                request_params=query,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="retrieve_aligns",
            gene_symbol=gene_symbol or sid,
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
            endpoint_name="retrieve_aligns",
            gene_symbol=gene_symbol or sid,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="retrieve_aligns",
            gene_symbol=gene_symbol or sid,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="retrieve_aligns",
            gene_symbol=gene_symbol or sid,
            request_url=request_url,
            request_params=query,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_domains(
    queries: str,
    *,
    gene_symbol: str = "",
    include_aligns: bool = False,
    poll_interval_seconds: float = DEFAULT_POLL_INTERVAL_SECONDS,
    max_polls: int = DEFAULT_MAX_POLLS,
    settings: Settings | None = None,
) -> ToolResult:
    """Submit → poll until ``status=0`` → retrieve hits (+ optional aligns).

    Never raises. Preserves raw submit/poll/hits text for artifacts.
    """
    cfg = settings or get_settings()
    submitted = submit_search(queries, gene_symbol=gene_symbol, settings=cfg)
    if not submitted.success:
        return _tool_result(
            endpoint_name="fetch_domains",
            gene_symbol=gene_symbol or queries,
            request_url=submitted.request_url,
            request_params=submitted.request_params,
            success=False,
            status_code=submitted.status_code,
            data={"submit": submitted.data},
            error_type=submitted.error_type or "submit_failed",
            error_message=submitted.error_message or "CDD submit failed",
        )

    cdsid = (submitted.data or {}).get("cdsid") if isinstance(submitted.data, dict) else None
    status = (submitted.data or {}).get("status") if isinstance(submitted.data, dict) else None
    poll_history: list[Any] = [submitted.data]
    last_url = submitted.request_url
    last_params = submitted.request_params
    last_status_code = submitted.status_code

    polls = 0
    while status != STATUS_COMPLETED and polls < max_polls:
        time.sleep(max(0.0, poll_interval_seconds))
        polls += 1
        polled = poll_status(str(cdsid), gene_symbol=gene_symbol, settings=cfg)
        last_url = polled.request_url
        last_params = polled.request_params
        last_status_code = polled.status_code
        poll_history.append(polled.data)
        if not polled.success:
            return _tool_result(
                endpoint_name="fetch_domains",
                gene_symbol=gene_symbol or queries,
                request_url=polled.request_url,
                request_params=polled.request_params,
                success=False,
                status_code=polled.status_code,
                data={
                    "queries": queries,
                    "cdsid": cdsid,
                    "submit": submitted.data,
                    "poll_history": poll_history,
                },
                error_type=polled.error_type or "poll_failed",
                error_message=polled.error_message or "CDD status poll failed",
            )
        status = (polled.data or {}).get("status") if isinstance(polled.data, dict) else None
        cdsid = (polled.data or {}).get("cdsid") or cdsid

    if status != STATUS_COMPLETED:
        return _tool_result(
            endpoint_name="fetch_domains",
            gene_symbol=gene_symbol or queries,
            request_url=last_url,
            request_params=last_params,
            success=False,
            status_code=last_status_code,
            data={
                "queries": queries,
                "cdsid": cdsid,
                "status": status,
                "submit": submitted.data,
                "poll_history": poll_history,
            },
            error_type="timeout",
            error_message=(
                f"CDD search did not complete after {max_polls} polls "
                f"(last status={status!r})"
            ),
        )

    hits = retrieve_hits(str(cdsid), gene_symbol=gene_symbol, settings=cfg)
    if not hits.success:
        return _tool_result(
            endpoint_name="fetch_domains",
            gene_symbol=gene_symbol or queries,
            request_url=hits.request_url,
            request_params=hits.request_params,
            success=False,
            status_code=hits.status_code,
            data={
                "queries": queries,
                "cdsid": cdsid,
                "submit": submitted.data,
                "poll_history": poll_history,
                "hits": hits.data,
            },
            error_type=hits.error_type or "hits_failed",
            error_message=hits.error_message or "CDD hits retrieve failed",
        )

    aligns_payload: Any = None
    if include_aligns:
        aligns = retrieve_aligns(str(cdsid), gene_symbol=gene_symbol, settings=cfg)
        if not aligns.success:
            return _tool_result(
                endpoint_name="fetch_domains",
                gene_symbol=gene_symbol or queries,
                request_url=aligns.request_url,
                request_params=aligns.request_params,
                success=False,
                status_code=aligns.status_code,
                data={
                    "queries": queries,
                    "cdsid": cdsid,
                    "submit": submitted.data,
                    "poll_history": poll_history,
                    "hits": hits.data,
                    "aligns": aligns.data,
                },
                error_type=aligns.error_type or "aligns_failed",
                error_message=aligns.error_message or "CDD aligns retrieve failed",
            )
        aligns_payload = aligns.data
        last_url = aligns.request_url
        last_params = aligns.request_params
        last_status_code = aligns.status_code
    else:
        last_url = hits.request_url
        last_params = hits.request_params
        last_status_code = hits.status_code

    hit_data = hits.data if isinstance(hits.data, dict) else {}
    return _tool_result(
        endpoint_name="fetch_domains",
        gene_symbol=gene_symbol or queries,
        request_url=last_url,
        request_params={
            "queries": queries,
            "cdsid": cdsid,
            "include_aligns": include_aligns,
            **(last_params or {}),
        },
        success=True,
        status_code=last_status_code,
        data={
            "queries": queries,
            "gene_symbol": gene_symbol or None,
            "cdsid": cdsid,
            "status": status,
            "submit": submitted.data,
            "poll_history": poll_history,
            "hits": hits.data,
            "hit_rows": hit_data.get("hit_rows") or [],
            "hit_summaries": hit_data.get("hit_summaries") or [],
            "hit_count": hit_data.get("hit_count") or 0,
            "aligns": aligns_payload,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "BWRPSB_URL",
    "DEFAULT_DB",
    "DEFAULT_MAXHIT",
    "DEFAULT_EVALUE",
    "DEFAULT_QUERY_SREBF2",
    "STATUS_COMPLETED",
    "parse_status_text",
    "parse_hits_text",
    "summarize_hit",
    "submit_search",
    "poll_status",
    "retrieve_hits",
    "retrieve_aligns",
    "fetch_domains",
]
