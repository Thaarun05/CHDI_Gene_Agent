"""NCBI Datasets taxonomy lineage client (Section 1e scope membership).

Endpoint (batched path ids)::

    GET https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{id1,id2,...}

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Iterable, Sequence
from urllib.parse import quote

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult
from gene_dossier.tools.ncbi_datasets import DATASETS_BASE, _headers, _safe_params

SOURCE_NAME = "NCBI Datasets"
DEFAULT_TAXONOMY_BATCH_SIZE = 20
DEFAULT_TAXONOMY_MAX_CONCURRENCY = 5
DEFAULT_TAXONOMY_MAX_ATTEMPTS = 4
DEFAULT_TAXONOMY_RETRY_SLEEP_SECONDS = 0.75
DEFAULT_TAXONOMY_RESOLVE_ROUNDS = 3
_TRANSIENT_HTTP_STATUS_CODES = frozenset({408, 425, 429, 500, 502, 503, 504})


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


def extract_taxonomy_nodes(payload: Any) -> list[dict[str, Any]]:
    if not isinstance(payload, dict):
        return []
    nodes = payload.get("taxonomy_nodes")
    if not isinstance(nodes, list):
        return []
    return [node for node in nodes if isinstance(node, dict)]


def lineage_tax_ids(taxonomy: dict[str, Any]) -> list[int]:
    raw = taxonomy.get("lineage")
    out: list[int] = []
    if isinstance(raw, list):
        for item in raw:
            try:
                out.append(int(item))
            except (TypeError, ValueError):
                continue
    tax_id = taxonomy.get("tax_id")
    try:
        tid = int(tax_id)
    except (TypeError, ValueError):
        return out
    if tid not in out:
        out.append(tid)
    return out


def membership_from_lineage(
    lineage: Sequence[int],
    *,
    tax_id: int,
    scope_tax_id: int,
) -> bool | None:
    """Return True/False membership, or None when lineage is empty/unknown."""
    if not lineage:
        return None
    ids = set(int(x) for x in lineage)
    ids.add(int(tax_id))
    return int(scope_tax_id) in ids


def _is_transient_taxonomy_failure(result: ToolResult) -> bool:
    if result.success:
        return False
    if result.error_type == "timeout":
        return True
    if result.error_type == "http_error":
        if result.status_code is None:
            return True
        return int(result.status_code) in _TRANSIENT_HTTP_STATUS_CODES
    # Unexpected transport / client exceptions are usually worth one more try.
    return result.error_type not in {"invalid_request"}


def _fetch_taxonomy_taxon_once(
    tax_ids: Sequence[str | int],
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch taxonomy metadata for one or more taxon IDs (single attempt)."""
    cfg = settings or get_settings()
    cleaned = [str(tid).strip() for tid in tax_ids if str(tid).strip()]
    if not cleaned:
        return _tool_result(
            endpoint_name="taxonomy_taxon",
            gene_symbol=gene_symbol,
            request_url=f"{DATASETS_BASE}/taxonomy/taxon/",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="tax_ids are required",
        )
    joined = ",".join(cleaned)
    path = f"taxonomy/taxon/{quote(joined, safe=',')}"
    url = f"{DATASETS_BASE}/{path}"
    headers = _headers(cfg)
    api_key_used = "api-key" in headers
    safe = _safe_params({"taxons": joined}, api_key_used=api_key_used)
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(url, headers=headers)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}
        if response.is_success:
            return _tool_result(
                endpoint_name="taxonomy_taxon",
                gene_symbol=gene_symbol,
                request_url=url,
                request_params=safe,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="taxonomy_taxon",
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=safe,
            success=False,
            status_code=response.status_code,
            data=payload,
            error_type="http_error",
            error_message=f"HTTP {response.status_code}",
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name="taxonomy_taxon",
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=safe,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="taxonomy_taxon",
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=safe,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name="taxonomy_taxon",
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=safe,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_taxonomy_taxon(
    tax_ids: Sequence[str | int],
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
    max_attempts: int = DEFAULT_TAXONOMY_MAX_ATTEMPTS,
    retry_sleep_seconds: float = DEFAULT_TAXONOMY_RETRY_SLEEP_SECONDS,
) -> ToolResult:
    """Fetch taxonomy metadata for one or more taxon IDs (comma-batched path).

    Transient HTTP / timeout failures are retried with exponential backoff.
    """
    attempts_allowed = max(1, int(max_attempts))
    sleep_seconds = max(0.0, float(retry_sleep_seconds))
    attempts: list[dict[str, Any]] = []
    result: ToolResult | None = None

    for attempt in range(1, attempts_allowed + 1):
        result = _fetch_taxonomy_taxon_once(
            tax_ids, gene_symbol=gene_symbol, settings=settings
        )
        attempts.append(
            {
                "attempt": attempt,
                "success": bool(result.success),
                "status_code": result.status_code,
                "error_type": result.error_type,
                "error_message": result.error_message,
            }
        )
        if result.success:
            if len(attempts) > 1:
                return result.model_copy(
                    update={
                        "request_params": {
                            **(result.request_params or {}),
                            "request_attempts": attempts,
                            "retry_count": len(attempts) - 1,
                        }
                    }
                )
            return result
        if not _is_transient_taxonomy_failure(result) or attempt >= attempts_allowed:
            break
        if sleep_seconds > 0:
            time.sleep(sleep_seconds * (2 ** (attempt - 1)))

    assert result is not None
    if len(attempts) > 1:
        return result.model_copy(
            update={
                "request_params": {
                    **(result.request_params or {}),
                    "request_attempts": attempts,
                    "retry_count": len(attempts) - 1,
                }
            }
        )
    return result


def _apply_taxonomy_result(
    tr: ToolResult,
    *,
    membership: dict[str, bool | None],
    seen: set[str],
    scope_tax_id: int,
) -> dict[str, int]:
    """Merge one batch ToolResult into membership. Returns counters for this result."""
    resolved = 0
    unknown = 0
    failed = 0
    if not tr.success or not isinstance(tr.data, dict):
        failed += 1
        # Do not permanently mark failed taxa here; caller retries unresolved IDs.
        return {
            "resolved": resolved,
            "unknown": unknown,
            "failed": failed,
        }

    for node in extract_taxonomy_nodes(tr.data):
        taxonomy = node.get("taxonomy")
        if not isinstance(taxonomy, dict):
            continue
        try:
            tid = str(int(taxonomy.get("tax_id")))
        except (TypeError, ValueError):
            continue
        lineage = lineage_tax_ids(taxonomy)
        try:
            tax_int = int(tid)
        except (TypeError, ValueError):
            membership[tid] = None
            seen.add(tid)
            unknown += 1
            continue
        member = membership_from_lineage(
            lineage, tax_id=tax_int, scope_tax_id=scope_tax_id
        )
        membership[tid] = member
        seen.add(tid)
        if member is None:
            unknown += 1
        else:
            resolved += 1
    return {
        "resolved": resolved,
        "unknown": unknown,
        "failed": failed,
    }


def _pending_tax_ids(
    requested: Sequence[str],
    *,
    membership: dict[str, bool | None],
    seen: set[str],
) -> list[str]:
    """Taxa still needing a successful network response (or cache hit)."""
    out: list[str] = []
    for tid in requested:
        if tid in seen:
            continue
        if membership.get(tid) is not None:
            continue
        out.append(tid)
    return out


def taxonomy_retrieval_complete(
    *,
    requested_count: int,
    resolved_count: int,
    unresolved_count: int,
    failed_request_count: int,
) -> bool:
    """True when every requested taxon has a True/False membership and no request failures remain."""
    if int(failed_request_count) > 0:
        return False
    if int(unresolved_count) > 0:
        return False
    if int(requested_count) == 0:
        return True
    return int(resolved_count) >= int(requested_count)


def resolve_taxonomy_memberships(
    tax_ids: Iterable[str | int],
    *,
    scope_tax_id: int,
    gene_symbol: str = "",
    cache: dict[str, bool | None] | None = None,
    batch_size: int = DEFAULT_TAXONOMY_BATCH_SIZE,
    max_concurrency: int = DEFAULT_TAXONOMY_MAX_CONCURRENCY,
    max_attempts: int = DEFAULT_TAXONOMY_MAX_ATTEMPTS,
    retry_sleep_seconds: float = DEFAULT_TAXONOMY_RETRY_SLEEP_SECONDS,
    resolve_rounds: int = DEFAULT_TAXONOMY_RESOLVE_ROUNDS,
    settings: Settings | None = None,
) -> tuple[dict[str, bool | None], list[ToolResult], dict[str, Any]]:
    """Resolve scope membership for distinct tax IDs with batching + cache.

    Failed / incomplete batches are retried with backoff across resolve rounds
    before unresolved taxa are marked unknown. Returns
    ``(membership_by_tax_id, tool_results, audit)``.
    """
    cfg = settings or get_settings()
    membership: dict[str, bool | None] = dict(cache or {})
    results: list[ToolResult] = []
    requested = sorted({str(tid).strip() for tid in tax_ids if str(tid).strip()})
    seen: set[str] = {
        tid for tid in requested if membership.get(tid) is not None
    }
    initial_cache_hits = len(seen)
    rounds_run = 0
    batch_attempts = 0
    retry_batch_count = 0
    total_request_retries = 0
    failed_after_retries = 0
    resolved_delta = 0
    unknown_from_empty_lineage = 0

    def _run_batch(batch: list[str]) -> ToolResult:
        return fetch_taxonomy_taxon(
            batch,
            gene_symbol=gene_symbol,
            settings=cfg,
            max_attempts=max_attempts,
            retry_sleep_seconds=retry_sleep_seconds,
        )

    rounds_allowed = max(1, int(resolve_rounds))
    for round_idx in range(1, rounds_allowed + 1):
        pending = _pending_tax_ids(requested, membership=membership, seen=seen)
        if not pending:
            break
        rounds_run = round_idx
        # Shrink batch size on later rounds to reduce partial/timeout failures.
        effective_batch = max(1, batch_size if round_idx == 1 else max(1, batch_size // 2))
        if round_idx >= 3:
            effective_batch = 1
        batches = [
            pending[i : i + effective_batch]
            for i in range(0, len(pending), effective_batch)
        ]
        if round_idx > 1:
            retry_batch_count += len(batches)
            if retry_sleep_seconds > 0:
                time.sleep(float(retry_sleep_seconds) * (2 ** (round_idx - 2)))

        round_results: list[ToolResult] = []
        workers = max(1, min(max_concurrency, len(batches)))
        if workers == 1:
            for batch in batches:
                round_results.append(_run_batch(batch))
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_run_batch, batch): batch for batch in batches}
                for future in as_completed(futures):
                    round_results.append(future.result())

        for tr in round_results:
            batch_attempts += 1
            params = tr.request_params or {}
            total_request_retries += int(params.get("retry_count") or 0)
            counters = _apply_taxonomy_result(
                tr, membership=membership, seen=seen, scope_tax_id=scope_tax_id
            )
            resolved_delta += counters["resolved"]
            unknown_from_empty_lineage += counters["unknown"]
            if counters["failed"]:
                failed_after_retries += 1
            results.append(tr)

    # Exhausted retries: mark still-missing taxa as unknown so callers can audit.
    unresolved = _pending_tax_ids(requested, membership=membership, seen=seen)
    for tid in unresolved:
        membership[tid] = None

    resolved_final = sum(
        1 for tid in requested if membership.get(tid) is not None
    )
    unknown_final = sum(1 for tid in requested if membership.get(tid) is None)
    # If everything resolved, prior failed-round counters should not block completeness.
    if unknown_final == 0:
        failed_after_retries = 0
    complete = taxonomy_retrieval_complete(
        requested_count=len(requested),
        resolved_count=resolved_final,
        unresolved_count=unknown_final,
        failed_request_count=failed_after_retries,
    )

    audit = {
        "requested_tax_ids": requested,
        "requested_count": len(requested),
        "cache_hits": initial_cache_hits,
        "resolved_count": resolved_final,
        "unknown_count": unknown_final,
        "unresolved_count": unknown_final,
        "failed_request_count": failed_after_retries,
        "batch_count": batch_attempts,
        "batch_size": batch_size,
        "max_concurrency": max_concurrency,
        "max_attempts": max_attempts,
        "retry_sleep_seconds": retry_sleep_seconds,
        "resolve_rounds_run": rounds_run,
        "retry_batch_count": retry_batch_count,
        "request_retry_count": total_request_retries,
        "resolved_from_network": resolved_delta,
        "empty_lineage_unknown_count": unknown_from_empty_lineage,
        "taxonomy_complete": complete,
        "scope_tax_id": scope_tax_id,
    }
    return membership, results, audit

__all__ = [
    "SOURCE_NAME",
    "DEFAULT_TAXONOMY_BATCH_SIZE",
    "DEFAULT_TAXONOMY_MAX_CONCURRENCY",
    "DEFAULT_TAXONOMY_MAX_ATTEMPTS",
    "DEFAULT_TAXONOMY_RETRY_SLEEP_SECONDS",
    "DEFAULT_TAXONOMY_RESOLVE_ROUNDS",
    "extract_taxonomy_nodes",
    "lineage_tax_ids",
    "membership_from_lineage",
    "taxonomy_retrieval_complete",
    "fetch_taxonomy_taxon",
    "resolve_taxonomy_memberships",
]
