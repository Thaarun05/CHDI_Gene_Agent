"""PubChem PUG REST client (BioAssay).

Fetches assay AIDs for an NCBI Gene ID, assay descriptions, and optional assay
CSV rows. Does **not** normalize into evidence records — that belongs in
``normalize/chemicals.py``.

Key endpoints (validated)::

    GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/geneid/{geneid}/aids/JSON
    GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/description/JSON
    GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/CSV

For SREBF2, the expected Entrez Gene ID is ``6721``.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import csv
import io
from typing import Any

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "PubChem"
PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

DEFAULT_MAX_DESCRIPTIONS = 25
DEFAULT_GENE_ID_SREBF2 = "6721"


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


def _request(
    *,
    endpoint_name: str,
    gene_symbol: str,
    url: str,
    request_params: dict[str, Any],
    settings: Settings,
    expect_json: bool = True,
) -> ToolResult:
    """GET a PubChem PUG URL and return :class:`ToolResult`."""
    try:
        with httpx.Client(timeout=settings.http_timeout_seconds) as client:
            response = client.get(url)
        if expect_json:
            try:
                payload: Any = response.json()
            except ValueError:
                payload = {"raw_text": response.text[:4000]}
        else:
            payload = {
                "raw_text": response.text,
                "content_type": response.headers.get("content-type"),
            }

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


def extract_aid_list(aids_payload: Any) -> list[int]:
    """Return AID integers from an ``.../aids/JSON`` payload."""
    if not isinstance(aids_payload, dict):
        return []
    ident = aids_payload.get("IdentifierList") or {}
    if not isinstance(ident, dict):
        return []
    aids = ident.get("AID") or []
    out: list[int] = []
    if not isinstance(aids, list):
        return out
    for aid in aids:
        try:
            out.append(int(aid))
        except (TypeError, ValueError):
            continue
    return out


def summarize_assay_description(payload: Any) -> dict[str, Any]:
    """Extract key fields from an assay description JSON payload (not evidence)."""
    if not isinstance(payload, dict):
        return {}
    containers = payload.get("PC_AssayContainer") or []
    if not isinstance(containers, list) or not containers:
        return {"raw_keys": sorted(payload.keys())}
    first = containers[0] if isinstance(containers[0], dict) else {}
    assay = first.get("assay") if isinstance(first, dict) else {}
    if not isinstance(assay, dict):
        assay = {}
    descr = assay.get("descr") if isinstance(assay.get("descr"), dict) else {}
    aid_obj = descr.get("aid") if isinstance(descr.get("aid"), dict) else {}
    return {
        "aid": aid_obj.get("id"),
        "name": descr.get("name"),
        "comment": descr.get("comment"),
    }


def parse_assay_csv(raw_csv: str) -> list[dict[str, str]]:
    """Parse PubChem assay CSV into row dicts (still not evidence records)."""
    text = (raw_csv or "").lstrip("\ufeff").strip()
    if not text:
        return []
    reader = csv.DictReader(io.StringIO(text))
    rows: list[dict[str, str]] = []
    for row in reader:
        if not isinstance(row, dict):
            continue
        cleaned = {
            str(k): ("" if v is None else str(v))
            for k, v in row.items()
            if k is not None
        }
        if any(v.strip() for v in cleaned.values()):
            rows.append(cleaned)
    return rows


def summarize_activity_row(row: dict[str, Any]) -> dict[str, Any]:
    """Extract common PubChem assay CSV activity fields (not evidence)."""
    return {
        "pubchem_cid": row.get("PUBCHEM_CID"),
        "activity_outcome": row.get("PUBCHEM_ACTIVITY_OUTCOME"),
        "standard_type": row.get("Standard Type"),
        "standard_value": row.get("Standard Value"),
        "standard_units": row.get("Standard Units"),
    }


def aids_by_geneid(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """List PubChem BioAssay AIDs for an NCBI Gene ID."""
    cfg = settings or get_settings()
    gid = str(gene_id).strip()
    url = f"{PUG_BASE}/assay/target/geneid/{gid}/aids/JSON"
    return _request(
        endpoint_name="aids_by_geneid",
        gene_symbol=gene_symbol or gid,
        url=url,
        request_params={"gene_id": gid},
        settings=cfg,
        expect_json=True,
    )


def assay_description(
    aid: str | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PubChem assay description JSON for one AID."""
    cfg = settings or get_settings()
    aid_s = str(aid).strip()
    url = f"{PUG_BASE}/assay/aid/{aid_s}/description/JSON"
    return _request(
        endpoint_name="assay_description",
        gene_symbol=gene_symbol or aid_s,
        url=url,
        request_params={"aid": aid_s},
        settings=cfg,
        expect_json=True,
    )


def assay_csv(
    aid: str | int,
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch PubChem assay CSV for one AID (raw text preserved)."""
    cfg = settings or get_settings()
    aid_s = str(aid).strip()
    url = f"{PUG_BASE}/assay/aid/{aid_s}/CSV"
    result = _request(
        endpoint_name="assay_csv",
        gene_symbol=gene_symbol or aid_s,
        url=url,
        request_params={"aid": aid_s},
        settings=cfg,
        expect_json=False,
    )
    if not result.success:
        return result
    raw = ""
    content_type = None
    if isinstance(result.data, dict):
        raw = str(result.data.get("raw_text") or "")
        content_type = result.data.get("content_type")
    return _tool_result(
        endpoint_name="assay_csv",
        gene_symbol=gene_symbol or aid_s,
        request_url=result.request_url,
        request_params=result.request_params,
        success=True,
        status_code=result.status_code,
        data={
            "aid": aid_s,
            "raw_csv": raw,
            "content_type": content_type,
        },
    )


def fetch_bioassays(
    gene_id: str | int,
    *,
    gene_symbol: str = "",
    max_descriptions: int = DEFAULT_MAX_DESCRIPTIONS,
    csv_aids: list[str | int] | None = None,
    settings: Settings | None = None,
) -> ToolResult:
    """Fetch AIDs for a gene, optional descriptions, and optional assay CSV.

    On success, ``data`` includes::

        {
          "gene_id": ...,
          "gene_symbol": ...,
          "aids": [...],
          "aid_count": N,
          "aids_raw": <JSON>,
          "descriptions": {aid: raw JSON},
          "description_summaries": [...],
          "assay_csv": {aid: {raw_csv, rows, activity_summaries}},
        }

    Descriptions are capped by ``max_descriptions`` (0 skips them).
    CSV is fetched only for ``csv_aids`` when provided (preserve raw for artifacts).

    Never raises.
    """
    cfg = settings or get_settings()
    gid = str(gene_id).strip()
    aids_res = aids_by_geneid(gid, gene_symbol=gene_symbol, settings=cfg)
    if not aids_res.success:
        return _tool_result(
            endpoint_name="fetch_bioassays",
            gene_symbol=gene_symbol or gid,
            request_url=aids_res.request_url,
            request_params=aids_res.request_params,
            success=False,
            status_code=aids_res.status_code,
            data={"gene_id": gid, "aids_raw": aids_res.data},
            error_type=aids_res.error_type or "aids_failed",
            error_message=aids_res.error_message or "PubChem AID lookup failed",
        )

    aids = extract_aid_list(aids_res.data)
    descriptions: dict[str, Any] = {}
    description_summaries: list[dict[str, Any]] = []
    last_url = aids_res.request_url
    last_params = aids_res.request_params
    last_status = aids_res.status_code

    for aid in aids[: max(0, max_descriptions)]:
        desc = assay_description(aid, gene_symbol=gene_symbol, settings=cfg)
        last_url = desc.request_url
        last_params = desc.request_params
        last_status = desc.status_code
        if not desc.success:
            return _tool_result(
                endpoint_name="fetch_bioassays",
                gene_symbol=gene_symbol or gid,
                request_url=desc.request_url,
                request_params=desc.request_params,
                success=False,
                status_code=desc.status_code,
                data={
                    "gene_id": gid,
                    "aids": aids,
                    "aids_raw": aids_res.data,
                    "descriptions": descriptions,
                    "description_summaries": description_summaries,
                    "failed_aid": aid,
                },
                error_type=desc.error_type or "description_failed",
                error_message=desc.error_message
                or f"PubChem assay description failed for AID {aid}",
            )
        aid_key = str(aid)
        descriptions[aid_key] = desc.data
        summary = summarize_assay_description(desc.data)
        summary["aid"] = summary.get("aid") or aid
        description_summaries.append(summary)

    csv_payloads: dict[str, Any] = {}
    for aid in csv_aids or []:
        csv_res = assay_csv(aid, gene_symbol=gene_symbol, settings=cfg)
        last_url = csv_res.request_url
        last_params = csv_res.request_params
        last_status = csv_res.status_code
        if not csv_res.success:
            return _tool_result(
                endpoint_name="fetch_bioassays",
                gene_symbol=gene_symbol or gid,
                request_url=csv_res.request_url,
                request_params=csv_res.request_params,
                success=False,
                status_code=csv_res.status_code,
                data={
                    "gene_id": gid,
                    "aids": aids,
                    "aids_raw": aids_res.data,
                    "descriptions": descriptions,
                    "description_summaries": description_summaries,
                    "assay_csv": csv_payloads,
                    "failed_csv_aid": aid,
                },
                error_type=csv_res.error_type or "csv_failed",
                error_message=csv_res.error_message
                or f"PubChem assay CSV failed for AID {aid}",
            )
        raw_csv = ""
        if isinstance(csv_res.data, dict):
            raw_csv = str(csv_res.data.get("raw_csv") or "")
        rows = parse_assay_csv(raw_csv)
        csv_payloads[str(aid)] = {
            "aid": str(aid),
            "raw_csv": raw_csv,
            "content_type": (
                csv_res.data.get("content_type")
                if isinstance(csv_res.data, dict)
                else None
            ),
            "rows": rows,
            "activity_summaries": [summarize_activity_row(r) for r in rows],
            "row_count": len(rows),
        }

    return _tool_result(
        endpoint_name="fetch_bioassays",
        gene_symbol=gene_symbol or gid,
        request_url=last_url,
        request_params={
            "gene_id": gid,
            "max_descriptions": max_descriptions,
            "csv_aids": [str(a) for a in (csv_aids or [])],
            **last_params,
        },
        success=True,
        status_code=last_status,
        data={
            "gene_id": gid,
            "gene_symbol": gene_symbol or None,
            "aids": aids,
            "aid_count": len(aids),
            "aids_raw": aids_res.data,
            "descriptions": descriptions,
            "description_summaries": description_summaries,
            "assay_csv": csv_payloads,
        },
    )


__all__ = [
    "SOURCE_NAME",
    "PUG_BASE",
    "DEFAULT_MAX_DESCRIPTIONS",
    "DEFAULT_GENE_ID_SREBF2",
    "extract_aid_list",
    "summarize_assay_description",
    "parse_assay_csv",
    "summarize_activity_row",
    "aids_by_geneid",
    "assay_description",
    "assay_csv",
    "fetch_bioassays",
]
