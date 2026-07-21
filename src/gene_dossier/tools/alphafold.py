"""AlphaFold DB client (EBI).

Fetches predicted structure metadata for a UniProt accession. Does **not**
normalize into evidence records — that belongs in ``normalize/protein.py``.

Key endpoint (validated)::

    GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}

Build viewer link (not returned in JSON)::

    https://alphafold.ebi.ac.uk/entry/{uniprot_accession}

SREBF2 anchors: human ``Q12772``, mouse ``Q3U1N2``, rat ``Q3T1I5``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "AlphaFold"
ALPHAFOLD_API_BASE = "https://alphafold.ebi.ac.uk/api"
ALPHAFOLD_ENTRY_BASE = "https://alphafold.ebi.ac.uk/entry"

# Validated UniProt accessions for SREBF2 across species.
DEFAULT_ACCESSION_HUMAN = "Q12772"
DEFAULT_ACCESSION_MOUSE = "Q3U1N2"
DEFAULT_ACCESSION_RAT = "Q3T1I5"


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


def entry_url(uniprot_accession: str) -> str:
    """Build the AlphaFold browser entry URL for a UniProt accession."""
    return f"{ALPHAFOLD_ENTRY_BASE}/{quote(uniprot_accession.strip(), safe='')}"


def summarize_prediction(entry: dict[str, Any]) -> dict[str, Any]:
    """Extract key AlphaFold prediction fields (not evidence)."""
    accession = entry.get("uniprotAccession") or entry.get("uniprot_accession")
    return {
        "entry_id": entry.get("entryId") or entry.get("entry_id"),
        "model_entity_id": entry.get("modelEntityId") or entry.get("model_entity_id"),
        "uniprot_accession": accession,
        "uniprot_id": entry.get("uniprotId") or entry.get("uniprot_id"),
        "uniprot_description": entry.get("uniprotDescription")
        or entry.get("uniprot_description"),
        "organism_scientific_name": entry.get("organismScientificName")
        or entry.get("organism_scientific_name"),
        "tax_id": entry.get("taxId") or entry.get("tax_id"),
        "global_metric_value": entry.get("globalMetricValue")
        or entry.get("global_metric_value"),
        "fraction_plddt_very_low": entry.get("fractionPlddtVeryLow")
        or entry.get("fraction_plddt_very_low"),
        "fraction_plddt_low": entry.get("fractionPlddtLow")
        or entry.get("fraction_plddt_low"),
        "fraction_plddt_confident": entry.get("fractionPlddtConfident")
        or entry.get("fraction_plddt_confident"),
        "fraction_plddt_very_high": entry.get("fractionPlddtVeryHigh")
        or entry.get("fraction_plddt_very_high"),
        "latest_version": entry.get("latestVersion") or entry.get("latest_version"),
        "model_created_date": entry.get("modelCreatedDate")
        or entry.get("model_created_date"),
        "pdb_url": entry.get("pdbUrl") or entry.get("pdb_url"),
        "cif_url": entry.get("cifUrl") or entry.get("cif_url"),
        "bcif_url": entry.get("bcifUrl") or entry.get("bcif_url"),
        "pae_image_url": entry.get("paeImageUrl") or entry.get("pae_image_url"),
        "plddt_doc_url": entry.get("plddtDocUrl") or entry.get("plddt_doc_url"),
        "pae_doc_url": entry.get("paeDocUrl") or entry.get("pae_doc_url"),
        "entry_url": entry_url(str(accession)) if accession else None,
    }


def prediction(
    uniprot_accession: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch AlphaFold prediction metadata for a UniProt accession.

    On success, ``data`` is the raw JSON list from `/api/prediction/{accession}`.
    """
    cfg = settings or get_settings()
    accession = uniprot_accession.strip()
    if not accession:
        return _tool_result(
            endpoint_name="prediction",
            gene_symbol=gene_symbol,
            request_url=f"{ALPHAFOLD_API_BASE}/prediction/",
            request_params={},
            success=False,
            error_type="invalid_request",
            error_message="UniProt accession is required",
        )

    url = f"{ALPHAFOLD_API_BASE}/prediction/{quote(accession, safe='')}"
    request_params = {"uniprot_accession": accession}
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(url)
        try:
            payload: Any = response.json()
        except ValueError:
            payload = {"raw_text": response.text[:4000]}

        if response.is_success:
            return _tool_result(
                endpoint_name="prediction",
                gene_symbol=gene_symbol or accession,
                request_url=url,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data=payload,
            )
        return _tool_result(
            endpoint_name="prediction",
            gene_symbol=gene_symbol or accession,
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
            endpoint_name="prediction",
            gene_symbol=gene_symbol or accession,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="prediction",
            gene_symbol=gene_symbol or accession,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="prediction",
            gene_symbol=gene_symbol or accession,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def fetch_prediction(
    uniprot_accession: str,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch AlphaFold prediction and attach light summaries + entry URL.

    On success, ``data`` includes::

        {
          "uniprot_accession": ...,
          "gene_symbol": ...,
          "entry_url": ...,
          "predictions": <raw list>,
          "prediction_summaries": [...],
          "prediction_count": N,
        }
    """
    cfg = settings or get_settings()
    accession = uniprot_accession.strip()
    result = prediction(accession, gene_symbol=gene_symbol, settings=cfg)
    if not result.success:
        return _tool_result(
            endpoint_name="fetch_prediction",
            gene_symbol=gene_symbol or accession,
            request_url=result.request_url,
            request_params=result.request_params,
            success=False,
            status_code=result.status_code,
            data={
                "uniprot_accession": accession,
                "predictions": result.data,
                "entry_url": entry_url(accession) if accession else None,
            },
            error_type=result.error_type or "prediction_failed",
            error_message=result.error_message or "AlphaFold prediction failed",
        )

    raw = result.data
    entries = raw if isinstance(raw, list) else []
    summaries = [
        summarize_prediction(entry) for entry in entries if isinstance(entry, dict)
    ]
    return _tool_result(
        endpoint_name="fetch_prediction",
        gene_symbol=gene_symbol or accession,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "uniprot_accession": accession,
            "gene_symbol": gene_symbol or None,
            "entry_url": entry_url(accession),
            "predictions": raw,
            "prediction_summaries": summaries,
            "prediction_count": len(summaries),
        },
    )


__all__ = [
    "SOURCE_NAME",
    "ALPHAFOLD_API_BASE",
    "ALPHAFOLD_ENTRY_BASE",
    "DEFAULT_ACCESSION_HUMAN",
    "DEFAULT_ACCESSION_MOUSE",
    "DEFAULT_ACCESSION_RAT",
    "entry_url",
    "summarize_prediction",
    "prediction",
    "fetch_prediction",
]
