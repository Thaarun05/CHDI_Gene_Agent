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


def _request_json(
    *,
    endpoint_name: str,
    gene_symbol: str,
    url: str,
    request_params: dict[str, Any],
    settings: Settings,
) -> ToolResult:
    """GET a Harmonizome JSON URL and return :class:`ToolResult`."""
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name=endpoint_name,
                gene_symbol=gene_symbol,
                request_url=url,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=request_params,
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
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name=endpoint_name,
            gene_symbol=gene_symbol,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


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
    return {
        "associated_gene_symbol": gene.get("symbol"),
        "associated_gene_name": gene.get("name"),
        "dataset_name": dataset.get("name"),
        "attribute_name": attribute.get("name"),
        "attribute_href": attribute.get("href"),
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
        dataset = row.get("dataset") or {}
        name = dataset.get("name") if isinstance(dataset, dict) else None
        if is_tf_dataset(name):
            out.append(row)
    return out


def gene_associations(
    gene_symbol: str,
    *,
    show_associations: bool = True,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch Harmonizome gene record with associations."""
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    params = {"showAssociations": "true" if show_associations else "false"}
    path = f"{HARMONIZOME_BASE}/gene/{quote(symbol, safe='')}"
    url = f"{path}?{urlencode(params)}"
    return _request_json(
        endpoint_name="gene_associations",
        gene_symbol=symbol,
        url=url,
        request_params={"gene_symbol": symbol, **params},
        settings=cfg,
    )


def gene_set_associations(
    attribute_name: str,
    dataset_name: str,
    *,
    gene_symbol: str = "",
    show_associations: bool = True,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch Harmonizome gene-set association detail.

    Binding-site profiles use names like ``TEAD4_HepG2_hg19_1``, not just
    ``TEAD4``.
    """
    cfg = settings or get_settings()
    attribute = attribute_name.strip()
    dataset = dataset_name.strip()
    params = {"showAssociations": "true" if show_associations else "false"}
    # Dataset path segments use '+' for spaces in the validated Postman URL.
    attr_seg = quote(attribute, safe="")
    dataset_seg = quote(dataset, safe="").replace("%20", "+")
    path = f"{HARMONIZOME_BASE}/gene_set/{attr_seg}/{dataset_seg}"
    url = f"{path}?{urlencode(params)}"
    return _request_json(
        endpoint_name="gene_set_associations",
        gene_symbol=gene_symbol or attribute,
        url=url,
        request_params={
            "attribute_name": attribute,
            "dataset_name": dataset,
            **params,
        },
        settings=cfg,
    )


def fetch_gene_associations(
    gene_symbol: str,
    *,
    tf_only: bool = False,
    settings: Settings | None = None,
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
    result = gene_associations(gene_symbol, settings=cfg)
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
) -> ToolResult:
    """Convenience wrapper: gene associations filtered to TF/regulator datasets."""
    return fetch_gene_associations(gene_symbol, tf_only=True, settings=settings)


__all__ = [
    "SOURCE_NAME",
    "HARMONIZOME_BASE",
    "TF_DATASET_NAMES",
    "summarize_association",
    "is_tf_dataset",
    "filter_tf_associations",
    "gene_associations",
    "gene_set_associations",
    "fetch_gene_associations",
    "fetch_tf_associations",
]
