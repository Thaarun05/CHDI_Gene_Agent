"""Read-only scientific evidence overlap audits for agent preparation runs."""

from __future__ import annotations

import json
import re
from collections import Counter
from typing import Any

from gene_dossier.models import EvidenceRecord


_SPACE_RE = re.compile(r"\s+")


def _norm_text(value: Any) -> str:
    return _SPACE_RE.sub(" ", str(value or "").strip().casefold())


def _stable_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            str(key): _stable_value(val)
            for key, val in sorted(value.items())
            if key
            not in {
                "api_run_id",
                "raw_artifact_id",
                "dossier_run_id",
                "retrieved_at",
                "created_at",
                "timestamp",
            }
        }
    if isinstance(value, list):
        return [_stable_value(item) for item in value]
    if isinstance(value, str):
        return _norm_text(value)
    return value


def scientific_fingerprint(record: EvidenceRecord) -> str:
    """Return a deterministic content fingerprint excluding provenance IDs."""
    payload = {
        "gene": _norm_text(record.gene_symbol),
        "source": _norm_text(record.source_name),
        "source_id": _norm_text(record.source_id),
        "assertion_type": _norm_text(getattr(record.assertion_type, "value", record.assertion_type)),
        "fact_type": _norm_text(record.fact_type),
        "organism": _norm_text(record.organism),
        "species": _norm_text(record.species),
        "taxon_id": record.taxon_id,
        "value": _stable_value(record.value or {}),
        "display_text": _norm_text(record.display_text),
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def audit_evidence_overlap(
    previous: list[EvidenceRecord],
    new: list[EvidenceRecord],
) -> dict[str, Any]:
    previous_by_fp: dict[str, list[EvidenceRecord]] = {}
    new_by_fp: dict[str, list[EvidenceRecord]] = {}
    for record in previous:
        previous_by_fp.setdefault(scientific_fingerprint(record), []).append(record)
    for record in new:
        new_by_fp.setdefault(scientific_fingerprint(record), []).append(record)

    previous_fps = set(previous_by_fp)
    new_fps = set(new_by_fp)
    duplicate_fps = sorted(previous_fps & new_fps)

    previous_by_source_id: dict[str, set[str]] = {}
    new_by_source_id: dict[str, set[str]] = {}
    for fp, records in previous_by_fp.items():
        for record in records:
            previous_by_source_id.setdefault(record.source_id, set()).add(fp)
    for fp, records in new_by_fp.items():
        for record in records:
            new_by_source_id.setdefault(record.source_id, set()).add(fp)

    distinct_same_source_id = []
    for source_id in sorted(set(previous_by_source_id) & set(new_by_source_id)):
        if previous_by_source_id[source_id] != new_by_source_id[source_id]:
            distinct_same_source_id.append(source_id)

    def _count(records: list[EvidenceRecord], attr: str) -> dict[str, int]:
        return dict(Counter(str(getattr(record, attr)) for record in records))

    return {
        "exactDuplicateCount": sum(len(new_by_fp[fp]) for fp in duplicate_fps),
        "exactDuplicateSourceIds": sorted(
            {record.source_id for fp in duplicate_fps for record in new_by_fp[fp]}
        ),
        "semanticallyDistinctSharedSourceIds": distinct_same_source_id,
        "uniquePreviousCount": sum(len(previous_by_fp[fp]) for fp in sorted(previous_fps - new_fps)),
        "uniqueNewCount": sum(len(new_by_fp[fp]) for fp in sorted(new_fps - previous_fps)),
        "overlapBySource": {
            source: {
                "previous": _count(previous, "source_name").get(source, 0),
                "new": _count(new, "source_name").get(source, 0),
            }
            for source in sorted(set(_count(previous, "source_name")) | set(_count(new, "source_name")))
        },
        "overlapByAssertionType": {
            assertion: {
                "previous": _count(previous, "assertion_type").get(assertion, 0),
                "new": _count(new, "assertion_type").get(assertion, 0),
            }
            for assertion in sorted(set(_count(previous, "assertion_type")) | set(_count(new, "assertion_type")))
        },
    }


__all__ = ["audit_evidence_overlap", "scientific_fingerprint"]
