"""Harmonizome client (Ma'ayan Lab).

Fetches gene associations and optional gene-set detail. Does **not** normalize
into evidence records — that belongs in ``normalize/expression.py``.

Key endpoints (validated)::

    GET https://maayanlab.cloud/Harmonizome/api/1.0/gene/{symbol}?showAssociations=true
    GET https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/{attribute}/{dataset}
        ?showAssociations=true

TF / regulator datasets to prefer for the dossier table::

    ENCODE Transcription Factor Binding Site Profiles
    ENCODE Transcription Factor Targets
    ChEA Transcription Factor Binding Site Profiles
    ChEA Transcription Factor Targets
    JASPAR Predicted Transcription Factor Targets
    MotifMap Predicted Transcription Factor Targets

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import json
from typing import Any
from urllib.parse import quote, urlencode

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Harmonizome"
HARMONIZOME_BASE = "https://maayanlab.cloud/Harmonizome/api/1.0"

# Validated TF / regulator dataset names from the API map §4.1.
TF_DATASET_NAMES = (
    "ENCODE Transcription Factor Binding Site Profiles",
    "ENCODE Transcription Factor Targets",
    "ChEA Transcription Factor Binding Site Profiles",
    "ChEA Transcription Factor Targets",
    "JASPAR Predicted Transcription Factor Targets",
    "MotifMap Predicted Transcription Factor Targets",
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


def normalize_harmonizome_request_identity(
    *,
    method: str,
    url_path: str,
    gene_symbol: str,
    query_params: dict[str, Any] | None = None,
) -> str:
    """Stable identity for request dedupe: method + path + gene + sorted params."""
    params = dict(query_params or {})
    sorted_params = tuple(sorted((str(k), str(v)) for k, v in params.items()))
    path = str(url_path or "").rstrip("/")
    return "|".join(
        [
            str(method or "GET").upper(),
            path,
            str(gene_symbol or "").strip().upper(),
            repr(sorted_params),
        ]
    )


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    url: str,
    request_params: dict[str, Any],
    settings: Settings,
    query_params: dict[str, Any] | None = None,
    transient: Any | None = None,
) -> ToolResult:
    """GET a Harmonizome JSON URL and return :class:`ToolResult`.

    Computes SHA-256 of ``response.content`` before parsing. Exact response
    bytes are retained on the transient request cache when available.
    """
    import hashlib
    from datetime import datetime, timezone

    params = dict(query_params or {})
    requested_url = url if not params else f"{url}?{urlencode(params)}"
    identity = normalize_harmonizome_request_identity(
        method="GET",
        url_path=url,
        gene_symbol=gene_symbol,
        query_params=params,
    )
    if transient is not None:
        cached = transient.get_cached_request(identity)
        if isinstance(cached, dict) and cached.get("tool_result") is not None:
            return cached["tool_result"]

    try:
        retrieved_at = datetime.now(timezone.utc).isoformat()
        with httpx.Client(timeout=settings.http_timeout_seconds, follow_redirects=True) as client:
            response = client.get(url, params=params or None)
        final_url = str(response.url)
        content = bytes(response.content)
        raw_sha = hashlib.sha256(content).hexdigest()
        redirect_history = [str(r.url) for r in response.history]
        content_type = response.headers.get("content-type")
        decoding_method = "utf-8"
        replacement_count = 0
        decode_warning = None
        payload: Any
        try:
            text = content.decode("utf-8")
            payload = json.loads(text)
            decoding_method = "utf-8"
        except UnicodeDecodeError:
            text = content.decode("utf-8", errors="replace")
            replacement_count = text.count("\ufffd")
            decoding_method = "utf-8_replace"
            if replacement_count:
                decode_warning = (
                    f"UTF-8 decoding replaced {replacement_count} invalid byte sequence(s) "
                    "with U+FFFD; exact response bytes are preserved separately."
                )
            try:
                payload = json.loads(text)
            except Exception:  # noqa: BLE001
                payload = {"raw_text": text[:4000]}
                decoding_method = "utf-8_replace_unparsed"
        except ValueError:
            # Valid UTF-8 that is not JSON
            try:
                payload = {"raw_text": content.decode("utf-8", errors="replace")[:4000]}
            except Exception:  # noqa: BLE001
                payload = {"raw_text": str(content[:4000])}
            decoding_method = "utf-8_non_json"

        meta = {
            "requested_url": requested_url,
            "final_url": final_url,
            "redirect_history": redirect_history,
            "status_code": response.status_code,
            "content_type": content_type,
            "retrieved_at": retrieved_at,
            "response_body_sha256": raw_sha,
            "response_byte_length": len(content),
            "decoding_method": decoding_method,
            "utf8_replacement_char_count": replacement_count,
            "decode_warning": decode_warning,
            "request_identity": identity,
            "method": "GET",
        }
        if isinstance(payload, dict):
            # Strip prior meta if re-parsed; keep exact response bytes out of JSON payload.
            payload = {k: v for k, v in payload.items() if k != "_harmonizome_meta"}
            payload = {**payload, "_harmonizome_meta": meta}
        else:
            payload = {"value": payload, "_harmonizome_meta": meta}

        result = _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            # Requested URL only; final_url lives in _harmonizome_meta.
            request_url=requested_url,
            request_params=request_params,
            success=bool(response.is_success),
            status_code=response.status_code,
            data=payload,
            error_type=None if response.is_success else "http_error",
            error_message=None if response.is_success else f"HTTP {response.status_code}",
        )
        if transient is not None:
            transient.put_cached_request(
                identity,
                {
                    "tool_result": result,
                    "response_bytes": content,
                    "meta": meta,
                },
            )
        return result
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=requested_url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def _dataset_name_from_row(row: dict[str, Any]) -> str | None:
    """Resolve dataset name from live ``geneSet.name`` or legacy ``dataset``."""
    gene_set = row.get("geneSet")
    if isinstance(gene_set, dict):
        name = str(gene_set.get("name") or "").strip()
        if "/" in name:
            dataset = name.rsplit("/", 1)[-1].strip()
            if dataset:
                return dataset
    dataset = row.get("dataset") or {}
    if isinstance(dataset, dict):
        value = dataset.get("name")
        return str(value) if value else None
    return None


def summarize_association(row: dict[str, Any]) -> dict[str, Any]:
    """Extract key fields from one Harmonizome association (not evidence)."""
    gene = row.get("gene") or {}
    if not isinstance(gene, dict):
        gene = {}
    dataset = row.get("dataset") or {}
    if not isinstance(dataset, dict):
        dataset = {}
    attribute = row.get("attribute") or {}
    if not isinstance(attribute, dict):
        attribute = {}
    gene_set = row.get("geneSet") if isinstance(row.get("geneSet"), dict) else {}
    gene_set_name = str(gene_set.get("name") or "").strip()
    attribute_name = attribute.get("name")
    dataset_name = dataset.get("name")
    attribute_href = attribute.get("href") or gene_set.get("href")
    if gene_set_name and "/" in gene_set_name:
        attr_part, ds_part = gene_set_name.rsplit("/", 1)
        attribute_name = attribute_name or attr_part
        dataset_name = dataset_name or ds_part
    return {
        "associated_gene_symbol": gene.get("symbol"),
        "associated_gene_name": gene.get("name"),
        "dataset_name": dataset_name,
        "attribute_name": attribute_name,
        "attribute_href": attribute_href,
        "threshold_value": row.get("thresholdValue"),
        "standardized_value": row.get("standardizedValue"),
    }


def is_tf_dataset(dataset_name: str | None) -> bool:
    """True if ``dataset_name`` is one of the validated TF/regulator datasets."""
    if not dataset_name:
        return False
    return dataset_name in TF_DATASET_NAMES


def filter_tf_associations(
    associations: list[Any],
) -> list[dict[str, Any]]:
    """Return association dicts whose dataset is in :data:`TF_DATASET_NAMES`."""
    out: list[dict[str, Any]] = []
    for row in associations:
        if not isinstance(row, dict):
            continue
        name = _dataset_name_from_row(row)
        if is_tf_dataset(name):
            out.append(row)
    return out


def gene_associations(
    gene_symbol: str,
    *,
    show_associations: bool = True,
    settings: Settings | None = None,
    transient: Any | None = None,
) -> ToolResult:
    """Fetch Harmonizome gene record with associations."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    params = {"showAssociations": "true" if show_associations else "false"}
    path = f"{HARMONIZOME_BASE}/gene/{quote(symbol, safe='')}"
    return _request_json(
        endpoint_name="gene_associations",
        gene_symbol=symbol,
        url=path,
        request_params={"gene_symbol": symbol, **params},
        settings=cfg,
        query_params=params,
        transient=transient,
    )


def gene_set_associations(
    attribute_name: str,
    dataset_name: str,
    *,
    gene_symbol: str = "",
    show_associations: bool = True,
    settings: Settings | None = None,
    transient: Any | None = None,
) -> ToolResult:
    """Fetch Harmonizome gene-set association detail.

    Binding-site profiles use names like ``TEAD4_HepG2_hg19_1``, not just
    ``TEAD4``. Diagnostic helper only — Section 4a must not call this.
    """
    cfg = settings or get_settings()
    attribute = attribute_name.strip()
    dataset = dataset_name.strip()
    params = {"showAssociations": "true" if show_associations else "false"}
    # Dataset path segments use '+' for spaces in the validated Postman URL.
    attr_seg = quote(attribute, safe="")
    dataset_seg = quote(dataset, safe="").replace("%20", "+")
    path = f"{HARMONIZOME_BASE}/gene_set/{attr_seg}/{dataset_seg}"
    return _request_json(
        endpoint_name="gene_set_associations",
        gene_symbol=gene_symbol or attribute,
        url=path,
        request_params={
            "attribute_name": attribute,
            "dataset_name": dataset,
            **params,
        },
        settings=cfg,
        query_params=params,
        transient=transient,
    )


def fetch_gene_associations(
    gene_symbol: str,
    *,
    tf_only: bool = False,
    settings: Settings | None = None,
    transient: Any | None = None,
) -> ToolResult:
    """Fetch gene associations with light summaries.

    When ``tf_only=True``, keep only associations from
    :data:`TF_DATASET_NAMES` (for the TF/regulator table).

    On success, ``data`` includes::

        {
          "gene_symbol": ...,
          "name": ...,
          "synonyms": ...,
          "associations": <raw list>,
          "association_summaries": [...],
          "tf_associations": [...],      # filtered raw (always computed)
          "tf_summaries": [...],
          "association_count": N,
          "tf_count": M,
          "raw": <full payload>,
        }
    """
    cfg = settings or get_settings()
    result = gene_associations(gene_symbol, settings=cfg, transient=transient)
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_gene_associations",
            gene_symbol=gene_symbol,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={"raw": result.data},
            error_type=result.error_type or "gene_associations_failed",
            error_message=result.error_message
            or "Harmonizome gene associations failed",
        )

    payload = result.data if isinstance(result.data, dict) else {}
    associations = payload.get("associations") or []
    if not isinstance(associations, list):
        associations = []

    tf_raw = filter_tf_associations(associations)
    all_summaries = [
        summarize_association(row) for row in associations if isinstance(row, dict)
    ]
    tf_summaries = [summarize_association(row) for row in tf_raw]

    if tf_only:
        selected = tf_raw
        selected_summaries = tf_summaries
    else:
        selected = [row for row in associations if isinstance(row, dict)]
        selected_summaries = all_summaries

    return _tool_result(
        endpoint_name="fetch_gene_associations",
        gene_symbol=gene_symbol,
        request_url=result.request_url,
        request_params={**result.request_params, "tf_only": tf_only},
        success=True,
        status_code=result.status_code,
        data={
            "gene_symbol": payload.get("symbol") or gene_symbol,
            "name": payload.get("name"),
            "synonyms": payload.get("synonyms"),
            "associations": selected,
            "association_summaries": selected_summaries,
            "tf_associations": tf_raw,
            "tf_summaries": tf_summaries,
            "association_count": len(selected),
            "tf_count": len(tf_raw),
            "tf_only": tf_only,
            "raw": result.data,
        },
    )


def fetch_tf_associations(
    gene_symbol: str,
    *,
    settings: Settings | None = None,
    transient: Any | None = None,
) -> ToolResult:
    """Convenience wrapper: gene associations filtered to TF/regulator datasets."""
    return fetch_gene_associations(
        gene_symbol, tf_only=True, settings=settings, transient=transient
    )


# Section 4a helpers (ordered tuples, parsers, collectors) live alongside this
# client so production collection shares one module family without breaking the
# pre-4a public I/O contract above.
from gene_dossier.tools.harmonizome_section4a import (  # noqa: E402
    CURATED_TF_DATASET_ORDER,
    CURATED_TF_DATASETS,
    PARSER_VERSION,
    PREDICTED_TF_DATASET_ORDER,
    PREDICTED_TF_DATASETS,
    SECTION_4A_TF_DATASETS,
    SUPPLEMENTARY_SCOPE,
    collect_section_4a_from_payload,
    collect_section_4a_harmonizome,
    display_dataset_label,
    gene_page_url,
)


__all__ = [
    "SOURCE_NAME",
    "HARMONIZOME_BASE",
    "TF_DATASET_NAMES",
    "normalize_harmonizome_request_identity",
    "CURATED_TF_DATASET_ORDER",
    "CURATED_TF_DATASETS",
    "PREDICTED_TF_DATASET_ORDER",
    "PREDICTED_TF_DATASETS",
    "SECTION_4A_TF_DATASETS",
    "PARSER_VERSION",
    "SUPPLEMENTARY_SCOPE",
    "summarize_association",
    "is_tf_dataset",
    "filter_tf_associations",
    "gene_associations",
    "gene_set_associations",
    "fetch_gene_associations",
    "fetch_tf_associations",
    "collect_section_4a_from_payload",
    "collect_section_4a_harmonizome",
    "display_dataset_label",
    "gene_page_url",
]
