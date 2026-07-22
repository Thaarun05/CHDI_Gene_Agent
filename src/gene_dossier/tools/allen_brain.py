"""Allen Brain Atlas client (Human Brain Atlas microarray).

Looks up HBA microarray probes for a gene and fetches expression for one probe
at a time. Priority C scaffold: the expression service is encoding-sensitive.
Does **not** normalize into evidence records — that belongs in
``normalize/expression.py``.

Key endpoints (validated)::

    GET https://api.brain-map.org/api/v2/data/query.json?criteria=...
        model::Probe,rma::criteria,[probe_type$eq'DNA'],
        products[abbreviation$eq'HumanMA'],gene[acronym$eq'{symbol}'],
        rma::options[only$eq'probes.id,probes.name,...']

    GET https://api.brain-map.org/api/v2/data/query.json?criteria=...
        service::human_microarray_expression[probes$eq'{probe_id}']

NOTE: Validate / fetch expression one probe at a time. Criteria encoding is
strict — probe IDs must be single-quoted in the service clause, and the
criteria string is passed as a single query parameter for the HTTP client to
encode. ``fetch_hba_expression`` retries transient Allen expression soft-fails
(``success=false`` / Informatics service errors) with a short backoff.

SREBF2 validated probe IDs: ``1051154``, ``1067243``, ``1051153``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import re
import time
from typing import Any
from urllib.parse import urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Allen Brain Atlas"
ALLEN_API_BASE = "https://api.brain-map.org/api/v2"
QUERY_URL = f"{ALLEN_API_BASE}/data/query.json"

# Validated SREBF2 HBA probe IDs (do not apply globally to other genes).
DEFAULT_PROBES_SREBF2 = (1051154, 1067243, 1051153)

# Conservative defaults for transient Allen expression-service flakiness.
DEFAULT_MAX_EXPRESSION_ATTEMPTS = 3
DEFAULT_RETRY_SLEEP_SECONDS = 1.0


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


def _payload_preview(payload: Any, *, limit: int = 200) -> str:
    """Whitespace-collapsed short preview for soft-fail error messages."""
    if isinstance(payload, dict):
        msg = payload.get("msg") or payload.get("message") or payload
        text = str(msg)
    else:
        text = str(payload or "")
    collapsed = re.sub(r"\s+", " ", text).strip()
    if len(collapsed) <= limit:
        return collapsed
    return collapsed[: limit - 3] + "..."


def _allen_error_message(payload: Any, *, fallback: str) -> str:
    """Build a clear error_message from an Allen JSON body."""
    msg = None
    if isinstance(payload, dict):
        msg = payload.get("msg") or payload.get("message")
    if msg is not None and str(msg).strip():
        base = str(msg).strip()
        # Avoid duplicating the same text as a preview when msg is already clear.
        return base
    preview = _payload_preview(payload)
    if preview:
        return f"{fallback} (preview: {preview})"
    return fallback


def build_probe_lookup_criteria(gene_symbol: str) -> str:
    """Build RMA criteria for HBA DNA probes matching ``gene_symbol``."""
    symbol = gene_symbol.strip()
    return (
        "model::Probe,rma::criteria,"
        "[probe_type$eq'DNA'],"
        "products[abbreviation$eq'HumanMA'],"
        f"gene[acronym$eq'{symbol}'],"
        "rma::options[only$eq'probes.id,probes.name,genes.acronym,genes.name,genes.entrez_id']"
    )


def build_expression_criteria(probe_id: str | int) -> str:
    """Build service criteria for one HBA microarray probe.

    Probe IDs must be single-quoted in the RMA service clause; unquoted IDs
    return HTTP 200 with ``success=false`` ("Informatics service request failed").
    """
    pid = str(probe_id).strip()
    return f"service::human_microarray_expression[probes$eq'{pid}']"


def summarize_probe(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key probe fields from an Allen probe lookup row (not evidence)."""
    gene = row.get("gene") if isinstance(row.get("gene"), dict) else {}
    return {
        "id": row.get("id"),
        "name": row.get("name"),
        "gene_acronym": gene.get("acronym") or row.get("acronym"),
        "gene_name": gene.get("name"),
        "entrez_id": gene.get("entrez_id") or row.get("entrez_id"),
    }


def _request_query(
    *,
    endpoint_name: str,
    gene_symbol: str,
    criteria: str,
    settings: Settings,
    extra_params: dict[str, Any] | None = None,
) -> ToolResult:
    """GET Allen ``query.json`` with a criteria string (never raises)."""
    params = {"criteria": criteria, **(extra_params or {})}
    # httpx encodes the criteria value; request_url uses the same encoding.
    request_url = f"{QUERY_URL}?{urlencode(params)}"
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(QUERY_URL, params=params)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            # Allen often returns HTTP 200 with success=false in the JSON body.
            if isinstance(payload, dict) and payload.get("success") is False:
                return _tool_result(
                    endpoint_name=endpoint_name,
                    gene_symbol=gene_symbol,
                    request_url=request_url,
                    request_params={"criteria": criteria, **(extra_params or {})},
                    success=False,
                    status_code=response.status_code,
                    data=payload,
                    error_type="api_error",
                    error_message=_allen_error_message(
                        payload, fallback="Allen API returned success=false"
                    ),
                )
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=request_url,
                request_params={"criteria": criteria, **(extra_params or {})},
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params={"criteria": criteria, **(extra_params or {})},
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
            request_params={"criteria": criteria, **(extra_params or {})},
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params={"criteria": criteria, **(extra_params or {})},
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=request_url,
            request_params={"criteria": criteria, **(extra_params or {})},
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def probe_lookup(
    gene_symbol: str,
    *,
    settings: Settings | None = None,
) -> ToolResult:
    """Look up Allen HBA DNA microarray probes for ``gene_symbol``."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol:
        return _tool_result(
            endpoint_name="probe_lookup",
            gene_symbol=gene_symbol,
            request_url=QUERY_URL,
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="gene_symbol is required for Allen probe lookup",
        )
    criteria = build_probe_lookup_criteria(symbol)
    return _request_query(
        endpoint_name="probe_lookup",
        gene_symbol=symbol,
        criteria=criteria,
        settings=cfg,
    )


def microarray_expression(
    probe_id: str | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch HBA microarray expression for a single probe ID."""
    cfg = settings or get_settings()
    pid = str(probe_id).strip()
    if not pid:
        return _tool_result(
            endpoint_name="microarray_expression",
            gene_symbol=gene_symbol,
            request_url=QUERY_URL,
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="probe_id is required",
        )
    criteria = build_expression_criteria(pid)
    return _request_query(
        endpoint_name="microarray_expression",
        gene_symbol=gene_symbol or pid,
        criteria=criteria,
        settings=cfg,
        extra_params={"probe_id": pid},
    )


def extract_probe_ids(probe_payload: Any) -> list[int]:
    """Extract probe IDs from a probe_lookup JSON payload."""
    if not isinstance(probe_payload, dict):
        return []
    msg = probe_payload.get("msg")
    if not isinstance(msg, list):
        return []
    out: list[int] = []
    for row in msg:
        if not isinstance(row, dict):
            continue
        pid = row.get("id")
        if pid is None:
            continue
        try:
            out.append(int(pid))
        except (TypeError, ValueError):
            continue
    return out


def _is_transient_expression_failure(result: ToolResult) -> bool:
    """True when an expression soft-fail looks like Allen upstream flakiness."""
    if result.success or result.error_type != "api_error":
        return False
    msg = result.error_message or ""
    if "Informatics service request failed" in msg:
        return True
    if "success=false" in msg.lower():
        return True
    if isinstance(result.data, dict) and result.data.get("success") is False:
        return True
    return False


def _expression_attempt_record(
    *,
    attempt: int,
    result: ToolResult,
) -> dict[str, Any]:
    """Serialize one microarray_expression attempt for debugging."""
    return {
        "attempt": attempt,
        "success": result.success,
        "error_type": result.error_type,
        "error_message": result.error_message,
        "status_code": result.status_code,
        "request_url": result.request_url,
        "data": result.data,
    }


def fetch_hba_expression(
    gene_symbol: str,
    *,
    probe_ids: list[str | int] | None = None,
    max_probes: int = 3,
    max_expression_attempts: int = DEFAULT_MAX_EXPRESSION_ATTEMPTS,
    retry_sleep_seconds: float = DEFAULT_RETRY_SLEEP_SECONDS,
    settings: Settings | None = None,
) -> ToolResult:
    """Probe lookup, then expression for up to ``max_probes`` probes (one at a time).

    If ``probe_ids`` is provided, skip lookup and use those IDs. On success,
    ``data`` includes probe summaries and per-probe expression payloads.

    Transient Allen ``success=false`` / Informatics failures on
    ``microarray_expression`` are retried with a short backoff. Prior failed
    attempts are preserved under ``expression_attempts`` (and on ultimately
    failed probes under ``failed_probes``).

    Never raises.
    """
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    probe_payload: Any = None
    probe_summaries: list[dict[str, Any]] = []
    resolved_ids: list[int] = []
    attempts_allowed = max(1, int(max_expression_attempts))
    sleep_seconds = max(0.0, float(retry_sleep_seconds))

    if probe_ids is not None:
        for pid in probe_ids:
            try:
                resolved_ids.append(int(pid))
            except (TypeError, ValueError):
                continue
    else:
        probes = probe_lookup(symbol, settings=cfg)
        if not probes.success:
            return _tool_result(
                endpoint_name="fetch_hba_expression",
                gene_symbol=symbol,
                request_url=probes.request_url,
                request_params=probes.request_params,
                success=False,
                status_code=probes.status_code,
                data={"probe_lookup": probes.data},
                error_type=probes.error_type or "probe_lookup_failed",
                error_message=probes.error_message or "Allen probe lookup failed",
            )
        probe_payload = probes.data
        resolved_ids = extract_probe_ids(probes.data)
        msg = (probes.data or {}).get("msg") if isinstance(probes.data, dict) else []
        if isinstance(msg, list):
            probe_summaries = [
                summarize_probe(row) for row in msg if isinstance(row, dict)
            ]

    selected = resolved_ids[: max(0, max_probes)]
    expressions: dict[str, Any] = {}
    expression_attempts: dict[str, list[dict[str, Any]]] = {}
    failed_probes: list[dict[str, Any]] = []
    last_url = QUERY_URL
    last_params: dict[str, Any] = {
        "gene_symbol": symbol,
        "probe_ids": selected,
        "max_probes": max_probes,
        "max_expression_attempts": attempts_allowed,
        "retry_sleep_seconds": sleep_seconds,
    }
    last_status: int | None = None

    for pid in selected:
        attempts: list[dict[str, Any]] = []
        expr: ToolResult | None = None
        for attempt in range(1, attempts_allowed + 1):
            expr = microarray_expression(pid, gene_symbol=symbol, settings=cfg)
            last_url = expr.request_url
            last_params = expr.request_params
            last_status = expr.status_code
            attempts.append(_expression_attempt_record(attempt=attempt, result=expr))
            if expr.success:
                expressions[str(pid)] = expr.data
                expression_attempts[str(pid)] = attempts
                break
            if not _is_transient_expression_failure(expr) or attempt >= attempts_allowed:
                failed_probes.append(
                    {
                        "probe_id": pid,
                        "error_type": expr.error_type,
                        "error_message": expr.error_message,
                        "status_code": expr.status_code,
                        "request_url": expr.request_url,
                        "data": expr.data,
                        "expression_attempts": attempts,
                    }
                )
                expression_attempts[str(pid)] = attempts
                break
            if sleep_seconds > 0:
                time.sleep(sleep_seconds)

    if not expressions:
        first = failed_probes[0] if failed_probes else None
        detail = (
            (first or {}).get("error_message")
            or "Allen microarray expression failed for all selected probes"
        )
        return _tool_result(
            endpoint_name="fetch_hba_expression",
            gene_symbol=symbol,
            request_url=(first or {}).get("request_url") or last_url,
            request_params=last_params,
            success=False,
            status_code=(first or {}).get("status_code") or last_status,
            data={
                "gene_symbol": symbol,
                "probe_lookup": probe_payload,
                "probe_summaries": probe_summaries,
                "probe_ids": selected,
                "expressions": expressions,
                "failed_probes": failed_probes,
                "expression_attempts": expression_attempts,
            },
            error_type=(first or {}).get("error_type") or "expression_failed",
            error_message=detail,
        )

    return _tool_result(
        endpoint_name="fetch_hba_expression",
        gene_symbol=symbol,
        request_url=last_url,
        request_params={
            "gene_symbol": symbol,
            "probe_ids": selected,
            "max_probes": max_probes,
            "max_expression_attempts": attempts_allowed,
            "retry_sleep_seconds": sleep_seconds,
        },
        success=True,
        status_code=last_status,
        data={
            "gene_symbol": symbol,
            "probe_lookup": probe_payload,
            "probe_summaries": probe_summaries,
            "probe_ids": selected,
            "probe_count": len(selected),
            "expressions": expressions,
            "failed_probes": failed_probes,
            "expression_attempts": expression_attempts,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "ALLEN_API_BASE",
    "QUERY_URL",
    "DEFAULT_PROBES_SREBF2",
    "DEFAULT_MAX_EXPRESSION_ATTEMPTS",
    "DEFAULT_RETRY_SLEEP_SECONDS",
    "build_probe_lookup_criteria",
    "build_expression_criteria",
    "summarize_probe",
    "probe_lookup",
    "microarray_expression",
    "extract_probe_ids",
    "fetch_hba_expression",
]
