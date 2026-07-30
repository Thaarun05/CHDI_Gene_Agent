"""AlphaFold DB client (EBI).

Fetches predicted structure metadata for a UniProt accession. Does **not**
normalize into evidence records — that belongs in ``normalize/protein.py`` /
``section_1d.py``.

Key endpoint (validated)::

    GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}

Polished entry link (selected model)::

    https://alphafold.ebi.ac.uk/entry/{modelEntityId}

Prefer post-deprecation field names (``modelEntityId``, ``isUniProtReviewed``,
``sequenceStart``, …) with legacy aliases as compatibility fallbacks only.
``paeImageUrl`` is optional and must not be required.

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

# Acceptance anchors (not production gene branches).
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
    """Build the AlphaFold browser entry URL for a UniProt accession (fallback)."""
    return f"{ALPHAFOLD_ENTRY_BASE}/{quote(uniprot_accession.strip(), safe='')}"


def entry_url_for_model(model_entity_id: str) -> str:
    """Build the official AlphaFold entry URL for a selected model entity id."""
    return f"{ALPHAFOLD_ENTRY_BASE}/{quote(str(model_entity_id).strip(), safe='')}"


def model_entity_id(prediction: dict[str, Any]) -> str | None:
    """Prefer ``modelEntityId``; fall back to deprecated ``entryId``."""
    value = prediction.get("modelEntityId") or prediction.get("model_entity_id")
    if value:
        return str(value)
    legacy = prediction.get("entryId") or prediction.get("entry_id")
    return str(legacy) if legacy else None


def sequence_start(prediction: dict[str, Any]) -> Any:
    return prediction.get("sequenceStart") or prediction.get("uniprotStart")


def sequence_end(prediction: dict[str, Any]) -> Any:
    return prediction.get("sequenceEnd") or prediction.get("uniprotEnd")


def sequence_value(prediction: dict[str, Any]) -> Any:
    return prediction.get("sequence") or prediction.get("uniprotSequence")


def is_uniprot_reviewed(prediction: dict[str, Any]) -> bool | None:
    reviewed = prediction.get("isUniProtReviewed")
    if reviewed is None:
        reviewed = prediction.get("isReviewed")
    if reviewed is None:
        return None
    return bool(reviewed)


def is_uniprot_reference_proteome(prediction: dict[str, Any]) -> bool | None:
    reference = prediction.get("isUniProtReferenceProteome")
    if reference is None:
        reference = prediction.get("isReferenceProteome")
    if reference is None:
        return None
    return bool(reference)


def summarize_prediction(entry: dict[str, Any]) -> dict[str, Any]:
    """Extract key AlphaFold prediction fields (not evidence)."""
    accession = entry.get("uniprotAccession") or entry.get("uniprot_accession")
    mid = model_entity_id(entry)
    return {
        "model_entity_id": mid,
        "entry_id": entry.get("entryId") or entry.get("entry_id"),
        "uniprot_accession": accession,
        "uniprot_id": entry.get("uniprotId") or entry.get("uniprot_id"),
        "uniprot_description": entry.get("uniprotDescription")
        or entry.get("uniprot_description"),
        "organism_scientific_name": entry.get("organismScientificName")
        or entry.get("organism_scientific_name"),
        "tax_id": entry.get("taxId") or entry.get("tax_id"),
        "entity_type": entry.get("entityType") or entry.get("entity_type"),
        "is_complex": entry.get("isComplex")
        if "isComplex" in entry
        else entry.get("is_complex"),
        "is_uniprot": entry.get("isUniProt")
        if "isUniProt" in entry
        else entry.get("is_uniprot"),
        "is_uniprot_reviewed": is_uniprot_reviewed(entry),
        "is_uniprot_reference_proteome": is_uniprot_reference_proteome(entry),
        "provider_id": entry.get("providerId") or entry.get("provider_id"),
        "tool_used": entry.get("toolUsed") or entry.get("tool_used"),
        "sequence_start": sequence_start(entry),
        "sequence_end": sequence_end(entry),
        "sequence": sequence_value(entry),
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
        # Optional / may be absent after AFDB API deprecations — never required.
        "pae_image_url": entry.get("paeImageUrl") or entry.get("pae_image_url"),
        "plddt_doc_url": entry.get("plddtDocUrl") or entry.get("plddt_doc_url"),
        "pae_doc_url": entry.get("paeDocUrl") or entry.get("pae_doc_url"),
        "entry_url": entry_url_for_model(mid) if mid else (
            entry_url(str(accession)) if accession else None
        ),
    }


def _as_bool(value: Any) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return None


def _as_int(value: Any) -> int | None:
    try:
        if value is None or value == "":
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def select_canonical_monomer_prediction(
    predictions: Any,
    requested_accession: str,
    expected_taxon_id: int | None = None,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Select the exact canonical UniProt monomer for Rancho Section 1d.

    Hard requirements must all pass; otherwise return ``(None, diagnostics)``.
    Preference ranking only runs among hard-qualified candidates and never lets
    a higher ``latestVersion`` on an isoform beat the canonical F1 identity.
    """
    diagnostics: list[dict[str, Any]] = []
    requested = (requested_accession or "").strip()
    if not requested:
        diagnostics.append(
            {
                "code": "missing_requested_accession",
                "message": "No UniProt accession provided for selection",
                "severity": "error",
            }
        )
        return None, diagnostics

    if not isinstance(predictions, list):
        diagnostics.append(
            {
                "code": "invalid_predictions_payload",
                "message": "AlphaFold predictions payload is not a list",
                "severity": "error",
            }
        )
        return None, diagnostics

    expected_f1 = f"AF-{requested}-F1"
    qualified: list[tuple[tuple[Any, ...], dict[str, Any]]] = []

    for index, raw in enumerate(predictions):
        if not isinstance(raw, dict):
            diagnostics.append(
                {
                    "code": "rejected_non_object",
                    "index": index,
                    "message": "Prediction entry is not an object",
                    "severity": "info",
                }
            )
            continue

        accession = str(
            raw.get("uniprotAccession") or raw.get("uniprot_accession") or ""
        ).strip()
        mid = model_entity_id(raw)
        entity_type = raw.get("entityType") or raw.get("entity_type")
        is_complex = _as_bool(
            raw.get("isComplex") if "isComplex" in raw else raw.get("is_complex")
        )
        is_uniprot = _as_bool(
            raw.get("isUniProt") if "isUniProt" in raw else raw.get("is_uniprot")
        )
        tax_id = _as_int(raw.get("taxId") or raw.get("tax_id"))
        reviewed = is_uniprot_reviewed(raw)
        reference = is_uniprot_reference_proteome(raw)
        provider = str(raw.get("providerId") or raw.get("provider_id") or "")
        tool_used = str(raw.get("toolUsed") or raw.get("tool_used") or "")
        latest_version = _as_int(
            raw.get("latestVersion") or raw.get("latest_version")
        )

        reject_reason: str | None = None
        if accession != requested:
            reject_reason = "accession_mismatch"
        elif mid is None:
            reject_reason = "missing_model_entity_id"
        elif entity_type is not None and str(entity_type).strip().lower() != "protein":
            reject_reason = "entity_type_not_protein"
        elif is_uniprot is not True:
            reject_reason = "not_uniprot"
        elif is_complex is True:
            reject_reason = "is_complex"
        elif (
            expected_taxon_id is not None
            and tax_id is not None
            and tax_id != int(expected_taxon_id)
        ):
            reject_reason = "taxon_mismatch"

        if reject_reason:
            diagnostics.append(
                {
                    "code": f"rejected_{reject_reason}",
                    "index": index,
                    "model_entity_id": mid,
                    "uniprot_accession": accession or None,
                    "tax_id": tax_id,
                    "severity": "info",
                }
            )
            continue

        # Identity first, then soft preferences. latestVersion is last and only
        # among already identity-qualified candidates (never beats non-F1).
        sort_key = (
            0 if mid == expected_f1 else 1,
            0 if reviewed is True else 1,
            0 if reference is True else 1,
            0 if provider == "GDM" else 1,
            0 if "alphafold monomer" in tool_used.lower() else 1,
            -(latest_version if latest_version is not None else -1),
            mid or "",
        )
        qualified.append((sort_key, raw))

    if not qualified:
        diagnostics.append(
            {
                "code": "no_hard_qualified_candidate",
                "message": (
                    f"No canonical UniProt monomer matched hard requirements "
                    f"for {requested}"
                ),
                "severity": "warning",
            }
        )
        return None, diagnostics

    qualified.sort(key=lambda item: item[0])
    selected = qualified[0][1]
    diagnostics.append(
        {
            "code": "selected_canonical_monomer",
            "model_entity_id": model_entity_id(selected),
            "uniprot_accession": requested,
            "latest_version": selected.get("latestVersion")
            or selected.get("latest_version"),
            "severity": "info",
        }
    )
    return selected, diagnostics


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

    On success, ``data`` includes the full raw ``predictions`` list plus
    ``prediction_summaries``. Selection of a canonical monomer is performed by
    callers (Section 1d), not here.
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
    "entry_url_for_model",
    "model_entity_id",
    "sequence_start",
    "sequence_end",
    "sequence_value",
    "is_uniprot_reviewed",
    "is_uniprot_reference_proteome",
    "summarize_prediction",
    "select_canonical_monomer_prediction",
    "prediction",
    "fetch_prediction",
]
