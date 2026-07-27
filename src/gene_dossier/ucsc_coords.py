"""UCSC coordinate conversion and assembly display labels.

Display (browser/search) positions are 1-based inclusive.
REST /getData/track intervals use 0-based half-open [start, end).

Conversion (mandatory)::

    api_start_0_based = display_start_1_based - 1
    api_end_exclusive = display_end_1_based
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# Version-controlled UCSC assembly registry (outside presentation builders).
_ASSEMBLY_REGISTRY: dict[str, dict[str, str]] = {
    "hg38": {
        "assembly_accession_name": "GRCh38",
        "ucsc_database": "hg38",
        "ucsc_display_name": "Human Dec. 2013 (GRCh38/hg38)",
    },
    "hg19": {
        "assembly_accession_name": "GRCh37",
        "ucsc_database": "hg19",
        "ucsc_display_name": "Human Feb. 2009 (GRCh37/hg19)",
    },
}


@dataclass(frozen=True)
class GenomicInterval:
    """Paired API and display coordinates for one locus."""

    chrom: str
    api_start_0_based: int
    api_end_exclusive: int
    display_start_1_based: int
    display_end_1_based: int
    genome: str = "hg38"
    coordinate_system: str = "ucsc_hg38_0based_half_open"

    @property
    def display_position(self) -> str:
        return f"{self.chrom}:{self.display_start_1_based}-{self.display_end_1_based}"

    def as_dict(self) -> dict[str, Any]:
        return {
            "chromosome": self.chrom,
            "chrom": self.chrom,
            "api_start_0_based": self.api_start_0_based,
            "api_end_exclusive": self.api_end_exclusive,
            "display_start_1_based": self.display_start_1_based,
            "display_end_1_based": self.display_end_1_based,
            "display_position": self.display_position,
            "coordinate_system": self.coordinate_system,
            "genome": self.genome,
        }


def display_to_api_start(display_start_1_based: int) -> int:
    """Convert a 1-based display start to a 0-based API start."""
    return int(display_start_1_based) - 1


def display_to_api_end(display_end_1_based: int) -> int:
    """Convert a 1-based inclusive display end to an exclusive API end."""
    return int(display_end_1_based)


def api_to_display_start(api_start_0_based: int) -> int:
    """Convert a 0-based API start to a 1-based display start."""
    return int(api_start_0_based) + 1


def api_to_display_end(api_end_exclusive: int) -> int:
    """Convert an exclusive API end to a 1-based inclusive display end."""
    return int(api_end_exclusive)


def interval_from_display(
    chrom: str,
    display_start_1_based: int,
    display_end_1_based: int,
    *,
    genome: str = "hg38",
) -> GenomicInterval:
    """Build an interval from a human-facing display range."""
    start = int(display_start_1_based)
    end = int(display_end_1_based)
    if end < start:
        raise ValueError(f"invalid display interval {chrom}:{start}-{end}")
    return GenomicInterval(
        chrom=chrom,
        api_start_0_based=display_to_api_start(start),
        api_end_exclusive=display_to_api_end(end),
        display_start_1_based=start,
        display_end_1_based=end,
        genome=genome,
    )


def interval_from_api(
    chrom: str,
    api_start_0_based: int,
    api_end_exclusive: int,
    *,
    genome: str = "hg38",
) -> GenomicInterval:
    """Build an interval from UCSC REST / bigGenePred coordinates."""
    start = int(api_start_0_based)
    end = int(api_end_exclusive)
    if end < start:
        raise ValueError(f"invalid API interval {chrom}:{start}-{end}")
    return GenomicInterval(
        chrom=chrom,
        api_start_0_based=start,
        api_end_exclusive=end,
        display_start_1_based=api_to_display_start(start),
        display_end_1_based=api_to_display_end(end),
        genome=genome,
    )


def assembly_label(genome: str) -> dict[str, str]:
    """Return assembly display metadata from the versioned registry.

    When unknown, falls back to a short ``GRCh38/hg38``-style label without inventing
    the long historical display string.
    """
    key = (genome or "").strip()
    entry = _ASSEMBLY_REGISTRY.get(key)
    if entry:
        return {
            "assembly_display_name": entry["ucsc_display_name"],
            "assembly_accession_name": entry["assembly_accession_name"],
            "ucsc_database": entry["ucsc_database"],
            "assembly_label_source": "versioned_registry",
        }
    short = key or "unknown"
    if short == "hg38":
        short_label = "GRCh38/hg38"
    elif short == "hg19":
        short_label = "GRCh37/hg19"
    else:
        short_label = short
    return {
        "assembly_display_name": short_label,
        "assembly_accession_name": short_label.split("/")[0] if "/" in short_label else short,
        "ucsc_database": short,
        "assembly_label_source": "versioned_registry",
    }


__all__ = [
    "GenomicInterval",
    "display_to_api_start",
    "display_to_api_end",
    "api_to_display_start",
    "api_to_display_end",
    "interval_from_display",
    "interval_from_api",
    "assembly_label",
]
