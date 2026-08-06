"""Section 4a Harmonizome association parsing and selection.

Membership allowlists are sets. Ordered tuples alone drive selection,
presentation, workbook grouping, evidence order, and tests.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from typing import Any
from urllib.parse import quote

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

PARSER_VERSION = "section_4a_harmonizome_v1"
SUPPLEMENTARY_SCOPE = "exact Section 4a allowlisted transcription-factor datasets"
HARMONIZOME_SITE_BASE = "https://maayanlab.cloud/Harmonizome"
GENE_PAGE_TMPL = HARMONIZOME_SITE_BASE + "/gene/{gene}"

CURATED_TF_DATASET_ORDER = (
    "ENCODE Transcription Factor Binding Site Profiles",
    "ENCODE Transcription Factor Targets",
    "ChEA Transcription Factor Binding Site Profiles",
    "ChEA Transcription Factor Targets",
)
PREDICTED_TF_DATASET_ORDER = (
    "JASPAR Predicted Transcription Factor Targets",
    "MotifMap Predicted Transcription Factor Targets",
)

CURATED_TF_DATASETS = set(CURATED_TF_DATASET_ORDER)
PREDICTED_TF_DATASETS = set(PREDICTED_TF_DATASET_ORDER)
SECTION_4A_TF_DATASETS = CURATED_TF_DATASETS | PREDICTED_TF_DATASETS

PARSE_COMPLETE = "parsed_complete"
PARSE_PARTIAL = "parsed_partial"
PARSE_UNPARSED = "unparsed_attribute"
PARSE_UNSUPPORTED = "unsupported_format"

_ENCODE_BUILD_RE = re.compile(
    r"_(?P<build>(?:hg|mm)\d+(?:_\d+)?)$",
    re.IGNORECASE,
)
_CHEA_PMID_RE = re.compile(
    r"^(?P<tf>.+?)-(?P<pmid>\d+)-(?P<rest>.+)$",
)

_ORGANISM_MAP = {
    "human": "Human",
    "homo sapiens": "Human",
    "mouse": "Mouse",
    "mus musculus": "Mouse",
    "rat": "Rat",
    "rattus norvegicus": "Rat",
}


def absolute_harmonizome_href(href: str | None) -> str | None:
    if not href:
        return None
    text = str(href).strip()
    if not text:
        return None
    if text.startswith("http://") or text.startswith("https://"):
        return text
    if text.startswith("/"):
        return HARMONIZOME_SITE_BASE + text
    return HARMONIZOME_SITE_BASE + "/" + text.lstrip("/")


def gene_page_url(gene_symbol: str) -> str:
    return GENE_PAGE_TMPL.format(gene=quote((gene_symbol or "").strip(), safe=""))


def display_dataset_label(dataset_name: str) -> str:
    """Polished table label: ChEA -> CHEA; others unchanged."""
    name = dataset_name or ""
    if name.startswith("ChEA "):
        return "CHEA " + name[len("ChEA ") :]
    if name == "ChEA":
        return "CHEA"
    return name


def split_gene_set_name(gene_set_name: str) -> tuple[str | None, str | None]:
    """Split ``ATTRIBUTE/DATASET`` on the last slash."""
    text = (gene_set_name or "").strip()
    if not text or "/" not in text:
        return None, None
    attribute, dataset = text.rsplit("/", 1)
    attribute = attribute.strip()
    dataset = dataset.strip()
    if not attribute or not dataset:
        return None, None
    return attribute, dataset


def _organism_from_build(build: str | None) -> tuple[str | None, str | None]:
    if not build:
        return None, None
    lower = build.lower()
    if lower.startswith("hg"):
        return "Human", "genome_build_prefix"
    if lower.startswith("mm"):
        return "Mouse", "genome_build_prefix"
    return None, None


def _map_organism(token: str | None) -> str | None:
    """Return polished organism only for explicitly recognized taxa tokens."""
    if not token:
        return None
    return _ORGANISM_MAP.get(token.strip().lower())


def parse_encode_binding(attribute: str) -> dict[str, Any]:
    """ENCODE Binding Site Profiles: TF_tissue_build with right-anchored build."""
    raw = (attribute or "").strip()
    m = _ENCODE_BUILD_RE.search(raw)
    if not m:
        return {
            "association": raw or None,
            "tissue_cells": None,
            "organism": None,
            "organism_audit": None,
            "organism_derivation": None,
            "genome_build": None,
            "parse_status": PARSE_UNPARSED if raw else PARSE_UNSUPPORTED,
        }
    build = m.group("build")
    head = raw[: m.start()]
    if not head:
        return {
            "association": raw,
            "tissue_cells": None,
            "organism": None,
            "organism_audit": None,
            "organism_derivation": None,
            "genome_build": build,
            "parse_status": PARSE_PARTIAL,
        }
    parts = head.split("_")
    tf = parts[0]
    tissue = "_".join(parts[1:]) if len(parts) > 1 else None
    organism_audit, derivation = _organism_from_build(build)
    status = PARSE_COMPLETE if tf and tissue and build else PARSE_PARTIAL
    return {
        "association": tf or raw,
        "tissue_cells": tissue or None,
        "organism": None,  # polished ENCODE organism blank
        "organism_audit": organism_audit,
        "organism_derivation": derivation,
        "genome_build": build,
        "parse_status": status,
    }


def parse_chea_binding(attribute: str) -> dict[str, Any]:
    """ChEA Binding: TF-PMID-CONTEXT-ORGANISM with PMID anchor."""
    raw = (attribute or "").strip()
    m = _CHEA_PMID_RE.match(raw)
    if not m:
        return {
            "association": raw or None,
            "tissue_cells": None,
            "organism": None,
            "organism_audit": None,
            "organism_derivation": None,
            "genome_build": None,
            "parse_status": PARSE_UNPARSED if raw else PARSE_UNSUPPORTED,
            "pubmed_id": None,
        }
    tf = m.group("tf").strip()
    pmid = m.group("pmid")
    rest = m.group("rest").strip()
    # Last hyphen-separated token is organism; context may contain hyphens.
    if "-" not in rest:
        return {
            "association": tf or raw,
            "tissue_cells": rest or None,
            "organism": None,
            "organism_audit": None,
            "organism_derivation": None,
            "genome_build": None,
            "parse_status": PARSE_PARTIAL,
            "pubmed_id": pmid,
        }
    context, organism_token = rest.rsplit("-", 1)
    organism = _map_organism(organism_token)
    if organism is None:
        return {
            "association": tf or raw,
            "tissue_cells": context or None,
            "organism": None,
            "organism_audit": None,
            "organism_derivation": None,
            "organism_token_unparsed": organism_token,
            "genome_build": None,
            "parse_status": PARSE_PARTIAL,
            "pubmed_id": pmid,
        }
    status = PARSE_COMPLETE if tf and context and organism else PARSE_PARTIAL
    return {
        "association": tf or raw,
        "tissue_cells": context or None,
        "organism": organism,
        "organism_audit": organism,
        "organism_derivation": "attribute_token",
        "organism_token_unparsed": None,
        "genome_build": None,
        "parse_status": status,
        "pubmed_id": pmid,
    }


def parse_association_only(attribute: str) -> dict[str, Any]:
    raw = (attribute or "").strip()
    if not raw:
        return {
            "association": None,
            "tissue_cells": None,
            "organism": None,
            "organism_audit": None,
            "organism_derivation": None,
            "genome_build": None,
            "parse_status": PARSE_UNSUPPORTED,
        }
    return {
        "association": raw,
        "tissue_cells": None,
        "organism": None,
        "organism_audit": None,
        "organism_derivation": None,
        "genome_build": None,
        "parse_status": PARSE_COMPLETE,
    }


def parse_attribute_for_dataset(dataset_name: str, attribute: str) -> dict[str, Any]:
    if dataset_name == "ENCODE Transcription Factor Binding Site Profiles":
        return parse_encode_binding(attribute)
    if dataset_name == "ChEA Transcription Factor Binding Site Profiles":
        return parse_chea_binding(attribute)
    if dataset_name in SECTION_4A_TF_DATASETS:
        return parse_association_only(attribute)
    return {
        "association": (attribute or "").strip() or None,
        "tissue_cells": None,
        "organism": None,
        "organism_audit": None,
        "organism_derivation": None,
        "genome_build": None,
        "parse_status": PARSE_UNSUPPORTED,
    }


_PARSE_RANK = {
    PARSE_COMPLETE: 0,
    PARSE_PARTIAL: 1,
    PARSE_UNPARSED: 2,
    PARSE_UNSUPPORTED: 3,
}


def _blank(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"none", "null"}:
        return ""
    return text


def extract_association_fields(row: dict[str, Any]) -> dict[str, Any]:
    """Extract attribute/dataset/href from live geneSet.name or legacy shape."""
    gene_set = row.get("geneSet") if isinstance(row.get("geneSet"), dict) else {}
    gene_set_name = str(gene_set.get("name") or "").strip()
    href = gene_set.get("href")
    attribute, dataset = split_gene_set_name(gene_set_name)

    if not dataset:
        dataset_obj = row.get("dataset") if isinstance(row.get("dataset"), dict) else {}
        attribute_obj = row.get("attribute") if isinstance(row.get("attribute"), dict) else {}
        dataset = dataset_obj.get("name")
        attribute = attribute_obj.get("name")
        href = href or attribute_obj.get("href")
        if attribute and dataset and not gene_set_name:
            gene_set_name = f"{attribute}/{dataset}"

    return {
        "gene_set_name": gene_set_name or None,
        "attribute_name": attribute,
        "dataset_name": dataset,
        "href": absolute_harmonizome_href(href if isinstance(href, str) else None),
        "threshold_value": row.get("thresholdValue"),
        "standardized_value": row.get("standardizedValue"),
    }


def build_parsed_record(
    row: dict[str, Any],
    *,
    source_order: int,
    query_gene: str,
) -> dict[str, Any] | None:
    fields = extract_association_fields(row)
    dataset = fields.get("dataset_name")
    if not dataset or dataset not in SECTION_4A_TF_DATASETS:
        return None
    attribute = str(fields.get("attribute_name") or "")
    parsed = parse_attribute_for_dataset(str(dataset), attribute)
    identity = (
        query_gene,
        str(dataset),
        str(fields.get("gene_set_name") or f"{attribute}/{dataset}"),
        str(fields.get("href") or ""),
        int(source_order),
    )
    return {
        "query_gene": query_gene,
        "dataset_name": dataset,
        "dataset_display": display_dataset_label(str(dataset)),
        "gene_set_name": fields.get("gene_set_name"),
        "attribute_name": attribute,
        "href": fields.get("href"),
        "association": parsed.get("association"),
        "tissue_cells": parsed.get("tissue_cells"),
        "organism": parsed.get("organism"),
        "organism_audit": parsed.get("organism_audit"),
        "organism_derivation": parsed.get("organism_derivation"),
        "organism_token_unparsed": parsed.get("organism_token_unparsed"),
        "genome_build": parsed.get("genome_build"),
        "parse_status": parsed.get("parse_status"),
        "pubmed_id": parsed.get("pubmed_id"),
        "threshold_value": fields.get("threshold_value"),
        "standardized_value": fields.get("standardized_value"),
        "source_order": int(source_order),
        "identity_key": identity,
        "category": "curated" if dataset in CURATED_TF_DATASETS else "predicted",
        "displayed": False,
        "displayed_rank": None,
        "selection_reason": None,
    }


def presentation_row(record: dict[str, Any], *, curated: bool) -> dict[str, Any]:
    if curated:
        return {
            "association": _blank(record.get("association")),
            "dataset": _blank(record.get("dataset_display") or record.get("dataset_name")),
            "tissue_cells": _blank(record.get("tissue_cells")),
            "organism": _blank(record.get("organism")),
            "genome_build": _blank(record.get("genome_build")),
            "href": record.get("href"),
            "gene_set_name": record.get("gene_set_name"),
            "parse_status": record.get("parse_status"),
            "source_order": record.get("source_order"),
            "dataset_name": record.get("dataset_name"),
        }
    return {
        "predicted_association": _blank(record.get("association")),
        "dataset": _blank(record.get("dataset_display") or record.get("dataset_name")),
        "href": record.get("href"),
        "gene_set_name": record.get("gene_set_name"),
        "parse_status": record.get("parse_status"),
        "source_order": record.get("source_order"),
        "dataset_name": record.get("dataset_name"),
    }


def _display_dedupe_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        record.get("dataset_name"),
        record.get("association"),
        record.get("tissue_cells"),
        record.get("organism"),
        record.get("genome_build"),
        record.get("gene_set_name"),
    )


def round_robin_select(
    records_by_dataset: dict[str, list[dict[str, Any]]],
    dataset_order: tuple[str, ...],
    *,
    max_rows: int,
) -> list[dict[str, Any]]:
    """Deterministic round-robin across ordered datasets; prefer better parse."""
    queues: dict[str, list[dict[str, Any]]] = {}
    for name in dataset_order:
        rows = list(records_by_dataset.get(name) or [])
        rows.sort(
            key=lambda r: (
                _PARSE_RANK.get(str(r.get("parse_status")), 99),
                int(r.get("source_order") or 0),
            )
        )
        queues[name] = rows

    selected: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    indices = {name: 0 for name in dataset_order}
    while len(selected) < max_rows:
        progressed = False
        for name in dataset_order:
            if len(selected) >= max_rows:
                break
            idx = indices[name]
            queue = queues[name]
            while idx < len(queue):
                candidate = queue[idx]
                idx += 1
                key = _display_dedupe_key(candidate)
                if key in seen:
                    continue
                seen.add(key)
                selected.append(candidate)
                progressed = True
                break
            indices[name] = idx
        if not progressed:
            break
    return selected


def collect_section_4a_from_payload(
    payload: dict[str, Any],
    *,
    query_gene: str,
    max_displayed_curated: int = 14,
    max_displayed_predicted: int = 25,
) -> dict[str, Any]:
    associations = payload.get("associations") if isinstance(payload, dict) else None
    if not isinstance(associations, list):
        associations = []

    official_symbol = str(payload.get("symbol") or "").strip()
    curated_records: list[dict[str, Any]] = []
    predicted_records: list[dict[str, Any]] = []
    out_of_scope: Counter[str] = Counter()
    curated_by_ds: dict[str, list[dict[str, Any]]] = defaultdict(list)
    predicted_by_ds: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for i, row in enumerate(associations):
        if not isinstance(row, dict):
            continue
        fields = extract_association_fields(row)
        dataset = fields.get("dataset_name")
        if not dataset:
            out_of_scope["__missing_dataset__"] += 1
            continue
        if dataset not in SECTION_4A_TF_DATASETS:
            out_of_scope[str(dataset)] += 1
            continue
        record = build_parsed_record(row, source_order=i, query_gene=query_gene)
        if record is None:
            continue
        if record["category"] == "curated":
            curated_records.append(record)
            curated_by_ds[str(record["dataset_name"])].append(record)
        else:
            predicted_records.append(record)
            predicted_by_ds[str(record["dataset_name"])].append(record)

    curated_counts = {
        name: len(curated_by_ds.get(name) or []) for name in CURATED_TF_DATASET_ORDER
    }
    predicted_counts = {
        name: len(predicted_by_ds.get(name) or []) for name in PREDICTED_TF_DATASET_ORDER
    }

    curated_total = sum(curated_counts.values())
    predicted_total = sum(predicted_counts.values())

    curated_display_raw = round_robin_select(
        curated_by_ds,
        CURATED_TF_DATASET_ORDER,
        max_rows=max(1, int(max_displayed_curated)),
    )
    if predicted_total <= int(max_displayed_predicted):
        predicted_display_raw: list[dict[str, Any]] = []
        for name in PREDICTED_TF_DATASET_ORDER:
            predicted_display_raw.extend(
                sorted(
                    predicted_by_ds.get(name) or [],
                    key=lambda r: int(r.get("source_order") or 0),
                )
            )
    else:
        predicted_display_raw = round_robin_select(
            predicted_by_ds,
            PREDICTED_TF_DATASET_ORDER,
            max_rows=int(max_displayed_predicted),
        )

    curated_display = [presentation_row(r, curated=True) for r in curated_display_raw]
    predicted_display = [presentation_row(r, curated=False) for r in predicted_display_raw]

    curated_rank = {
        int(r.get("source_order") or -1): idx + 1
        for idx, r in enumerate(curated_display_raw)
    }
    predicted_rank = {
        int(r.get("source_order") or -1): idx + 1
        for idx, r in enumerate(predicted_display_raw)
    }
    predicted_reason = (
        "all_predicted_within_cap"
        if predicted_total <= int(max_displayed_predicted)
        else "round_robin_predicted_display"
    )
    for record in curated_records:
        order = int(record.get("source_order") or -1)
        rank = curated_rank.get(order)
        record["displayed"] = rank is not None
        record["displayed_rank"] = rank
        record["selection_reason"] = (
            "round_robin_curated_display" if rank is not None else "allowlisted_not_displayed"
        )
    for record in predicted_records:
        order = int(record.get("source_order") or -1)
        rank = predicted_rank.get(order)
        record["displayed"] = rank is not None
        record["displayed_rank"] = rank
        record["selection_reason"] = (
            predicted_reason if rank is not None else "allowlisted_not_displayed"
        )

    out_rows = [
        {"dataset_name": name, "association_count": count}
        for name, count in sorted(out_of_scope.items(), key=lambda kv: (-kv[1], kv[0]))
    ]

    return {
        "query_gene": query_gene,
        "official_symbol": official_symbol or None,
        "parser_version": PARSER_VERSION,
        "supplementary_scope": SUPPLEMENTARY_SCOPE,
        "total_association_count": int(len(associations)),
        "in_scope_association_count": curated_total + predicted_total,
        "curated_total": curated_total,
        "predicted_total": predicted_total,
        "displayed_curated_count": len(curated_display),
        "displayed_predicted_count": len(predicted_display),
        "curated_counts": curated_counts,
        "predicted_counts": predicted_counts,
        "curated_records": curated_records,
        "predicted_records": predicted_records,
        "curated_display": curated_display,
        "predicted_display": predicted_display,
        "out_of_scope_summary": out_rows,
        "out_of_scope_total": int(sum(out_of_scope.values())),
        "max_displayed_curated": int(max_displayed_curated),
        "max_displayed_predicted": int(max_displayed_predicted),
        "gene_page_url": gene_page_url(official_symbol or query_gene),
    }


def collect_section_4a_harmonizome(
    gene_symbol: str,
    *,
    max_displayed_curated: int = 14,
    max_displayed_predicted: int = 25,
    settings: Settings | None = None,
    transient: Any | None = None,
    tool_result: ToolResult | None = None,
    gene_associations_fn: Any | None = None,
) -> dict[str, Any]:
    """Production collector: one gene associations GET; never calls gene_set."""
    # Lazy import avoids circular import with tools.harmonizome re-exports.
    from gene_dossier.tools import harmonizome as hz

    cfg = settings or get_settings()
    symbol = (gene_symbol or "").strip()
    retrieved_at = datetime.now(timezone.utc).isoformat()
    fetch_fn = gene_associations_fn or hz.gene_associations
    result = tool_result or fetch_fn(
        symbol,
        show_associations=True,
        settings=cfg,
        transient=transient,
    )

    base = {
        "gene_symbol": symbol,
        "retrieved_at": retrieved_at,
        "tool_result": result,
        "parser_version": PARSER_VERSION,
        "supplementary_scope": SUPPLEMENTARY_SCOPE,
        "gene_page_url": gene_page_url(symbol),
    }
    if not result.success:
        return {
            **base,
            "scientific_status": "source_unavailable",
            "presentation_status": "failed",
            "payload": result.data if isinstance(result.data, dict) else {},
            "collection": None,
        }

    payload = result.data if isinstance(result.data, dict) else {}
    # Drop internal meta from scientific payload views while preserving on ToolResult.
    official = str(payload.get("symbol") or "").strip()
    if official and official.upper() != symbol.upper():
        return {
            **base,
            "scientific_status": "gene_mismatch",
            "presentation_status": "failed",
            "payload": payload,
            "official_symbol": official,
            "collection": None,
        }

    collection = collect_section_4a_from_payload(
        payload,
        query_gene=official or symbol,
        max_displayed_curated=max_displayed_curated,
        max_displayed_predicted=max_displayed_predicted,
    )
    curated_total = int(collection["curated_total"])
    predicted_total = int(collection["predicted_total"])
    if curated_total == 0 and predicted_total == 0:
        scientific = "no_associations"
        presentation = "failed"
    else:
        scientific = "success"
        presentation = "success"
        if not collection["curated_display"] and curated_total > 0:
            presentation = "partial"
        if not collection["predicted_display"] and predicted_total > 0:
            presentation = "partial"

    return {
        **base,
        "scientific_status": scientific,
        "presentation_status": presentation,
        "payload": payload,
        "official_symbol": official or symbol,
        "collection": collection,
    }


__all__ = [
    "CURATED_TF_DATASET_ORDER",
    "CURATED_TF_DATASETS",
    "GENE_PAGE_TMPL",
    "HARMONIZOME_SITE_BASE",
    "PARSER_VERSION",
    "PARSE_COMPLETE",
    "PARSE_PARTIAL",
    "PARSE_UNPARSED",
    "PARSE_UNSUPPORTED",
    "PREDICTED_TF_DATASET_ORDER",
    "PREDICTED_TF_DATASETS",
    "SECTION_4A_TF_DATASETS",
    "SUPPLEMENTARY_SCOPE",
    "absolute_harmonizome_href",
    "build_parsed_record",
    "collect_section_4a_from_payload",
    "collect_section_4a_harmonizome",
    "display_dataset_label",
    "extract_association_fields",
    "gene_page_url",
    "parse_attribute_for_dataset",
    "parse_chea_binding",
    "parse_encode_binding",
    "presentation_row",
    "round_robin_select",
    "split_gene_set_name",
]
