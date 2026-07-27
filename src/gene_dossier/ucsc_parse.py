"""UCSC search / knownGene parsing for Section 1b conservation evidence.

Parses nested ``positionMatches`` groups and bigGenePred track rows.
Does not emit EvidenceRecords (see ``normalize.ucsc_conservation``).
"""

from __future__ import annotations

import html
import re
from dataclasses import dataclass, field
from typing import Any, Sequence

from gene_dossier.ucsc_coords import (
    GenomicInterval,
    assembly_label,
    interval_from_api,
    interval_from_display,
)

_POSITION_RE = re.compile(
    r"^(?P<chrom>chr[\w.]+):(?P<start>\d+)-(?P<end>\d+)$",
    re.IGNORECASE,
)
_GENCODE_RELEASE_RE = re.compile(
    r"GENCODE\s+Version\s+(\d+)|gencodeV(\d+)|BasicV(\d+)|CompV(\d+)",
    re.IGNORECASE,
)
_ENST_RE = re.compile(r"\b(ENST\d+(?:\.\d+)?)\b", re.IGNORECASE)
_COMMA_ARRAY_RE = re.compile(r",\s*")


@dataclass
class ParseDiagnostic:
    code: str
    message: str
    severity: str = "warning"


@dataclass
class PositionalSearchMatch:
    genome: str
    group_name: str
    track_name: str
    track_description: str
    position: str
    chromosome: str
    display_start: int
    display_end: int
    matched_identifier: str | None
    transcript_id: str | None
    displayed_gene_symbol: str
    canonical_flag: bool
    release_version: int | None
    group_kind: str  # mane|knownGene|gencode_basic|gencode_comp|hgnc|refseq|other


@dataclass
class SearchInventory:
    genome: str
    matches: list[PositionalSearchMatch] = field(default_factory=list)
    basic_release: int | None = None
    basic_exact_gene_count: int = 0
    comprehensive_release: int | None = None
    comprehensive_exact_gene_count: int = 0
    mane_exact_gene_count: int = 0
    known_gene_exact_gene_count: int = 0
    refseq_exact_gene_count: int = 0
    diagnostics: list[ParseDiagnostic] = field(default_factory=list)

    @property
    def selected_display_interval(self) -> GenomicInterval | None:
        """Prefer MANE positional match, else knownGene canonical, else first match."""
        mane = [m for m in self.matches if m.group_kind == "mane"]
        if mane:
            m = mane[0]
            return interval_from_display(
                m.chromosome, m.display_start, m.display_end, genome=self.genome
            )
        known = [m for m in self.matches if m.group_kind == "knownGene" and m.canonical_flag]
        if known:
            m = known[0]
            return interval_from_display(
                m.chromosome, m.display_start, m.display_end, genome=self.genome
            )
        if self.matches:
            m = self.matches[0]
            return interval_from_display(
                m.chromosome, m.display_start, m.display_end, genome=self.genome
            )
        return None


@dataclass
class ExonInterval:
    start: int
    end: int


@dataclass
class TranscriptRow:
    raw: dict[str, Any]
    chrom: str
    chrom_start: int
    chrom_end: int
    name: str
    strand: str | None
    gene_name: str | None
    gene_name2: str | None
    gene_type: str | None
    transcript_class: str | None
    transcript_type: str | None
    source: str | None
    tag_raw: str | None
    tags: frozenset[str]
    tier_raw: str | None
    tiers: frozenset[str]
    rank: int | None
    block_count: int | None
    exons: list[ExonInterval]
    exon_ok: bool
    exon_diagnostics: list[ParseDiagnostic] = field(default_factory=list)

    @property
    def is_mane_select(self) -> bool:
        return "MANE_Select" in self.tags

    @property
    def is_ensembl_canonical(self) -> bool:
        return "Ensembl_canonical" in self.tags

    @property
    def is_gencode_primary(self) -> bool:
        return "GENCODE_Primary" in self.tags

    @property
    def is_canonical_tier(self) -> bool:
        return any(t.lower() == "canonical" for t in self.tiers)

    @property
    def is_protein_coding(self) -> bool:
        for val in (self.transcript_class, self.transcript_type, self.gene_type):
            if val and "protein" in str(val).lower() and "coding" in str(val).lower():
                return True
            if val and str(val).lower() in {"coding", "protein_coding"}:
                return True
        return False


@dataclass
class RegionalInventory:
    genome: str
    track_big_data_url: str | None
    track_data_time: str | None
    current_gencode_release: int | None
    release_source: str | None
    rows: list[TranscriptRow]
    exact_rows: list[TranscriptRow]
    excluded_neighbor_symbols: list[str]
    regional_total_count: int
    exact_gene_transcript_count: int
    excluded_neighbor_count: int
    malformed_row_count: int
    diagnostics: list[ParseDiagnostic] = field(default_factory=list)


def parse_position(position: str) -> dict[str, Any] | None:
    """Parse ``chrom:start-end`` display coordinates."""
    text = html.unescape((position or "").strip()).replace(",", "")
    match = _POSITION_RE.match(text)
    if not match:
        return None
    start = int(match.group("start"))
    end = int(match.group("end"))
    if end < start:
        return None
    return {
        "chrom": match.group("chrom"),
        "start": start,
        "end": end,
        "position": f"{match.group('chrom')}:{start}-{end}",
    }


def _decode(text: Any) -> str:
    return html.unescape(str(text or "")).strip()


def _leading_symbol(pos_name: str) -> str:
    text = _decode(pos_name)
    if " (" in text:
        return text.split(" (", 1)[0].strip()
    return text.strip()


def _parse_gencode_release(*texts: str) -> int | None:
    for text in texts:
        m = _GENCODE_RELEASE_RE.search(text or "")
        if not m:
            continue
        for g in m.groups():
            if g:
                return int(g)
    return None


def _classify_group(name: str, track_name: str, description: str) -> str:
    blob = f"{name}|{track_name}|{description}".lower()
    if name.lower() == "mane" or "mane" == track_name.lower():
        return "mane"
    if name.lower() == "knowngene" or track_name.lower() == "knowngene":
        return "knownGene"
    if "trackdb" in blob or "publichub" in blob or "public hub" in blob:
        return "reject"
    if "gencode" in blob and "basic" in blob:
        return "gencode_basic"
    if "gencode" in blob and ("comp" in blob or "comprehensive" in blob):
        return "gencode_comp"
    if name.lower() == "hgnc" or "hgnc" in blob:
        return "hgnc"
    if "refseq" in blob or name.lower() in {"refgene", "ncbirefseqcurated", "ncbirefseqpredicted"}:
        return "refseq"
    return "other"


def _exact_match_for_group(
    *,
    group_kind: str,
    gene_symbol: str,
    pos_name: str,
    hg_find: str,
) -> tuple[bool, str | None]:
    """Return (is_exact, displayed_symbol). Never uses description text."""
    target = gene_symbol.strip().upper()
    pos = _decode(pos_name)
    find = _decode(hg_find)
    leading = _leading_symbol(pos)

    if group_kind in {"gencode_basic", "gencode_comp", "hgnc"}:
        # posName must equal the requested symbol after HTML decode/trim.
        if pos.upper() == target:
            return True, leading or pos
        return False, None

    if group_kind == "knownGene":
        if leading.upper() == target:
            return True, leading
        return False, None

    if group_kind == "mane":
        if leading.upper() != target:
            return False, None
        if not (_ENST_RE.search(find) or _ENST_RE.search(pos)):
            return False, None
        return True, leading

    if group_kind == "refseq":
        # Do not infer from NM_/NR_ alone without a validated gene symbol.
        if leading.upper() == target:
            return True, leading
        return False, None

    return False, None


def parse_search_response(
    payload: dict[str, Any] | None,
    *,
    gene_symbol: str,
    genome: str = "hg38",
) -> SearchInventory:
    """Parse nested UCSC /search ``positionMatches`` into exact-gene positional hits."""
    inv = SearchInventory(genome=genome)
    if not isinstance(payload, dict):
        inv.diagnostics.append(
            ParseDiagnostic("missing_search_summary", "Search payload missing or invalid")
        )
        return inv

    groups = payload.get("positionMatches") or []
    if not isinstance(groups, list):
        inv.diagnostics.append(
            ParseDiagnostic("missing_search_summary", "positionMatches is not a list")
        )
        return inv

    basic_by_release: dict[int, int] = {}
    comp_by_release: dict[int, int] = {}

    for group in groups:
        if not isinstance(group, dict):
            continue
        name = _decode(group.get("name"))
        track_name = _decode(group.get("trackName") or name)
        description = _decode(group.get("description"))
        kind = _classify_group(name, track_name, description)
        if kind == "reject":
            continue
        release = _parse_gencode_release(name, track_name, description)
        matches = group.get("matches") or []
        if not isinstance(matches, list):
            continue

        exact_in_group = 0
        for row in matches:
            if not isinstance(row, dict):
                continue
            pos_raw = row.get("position")
            parsed = parse_position(str(pos_raw) if pos_raw is not None else "")
            if not parsed:
                if pos_raw:
                    inv.diagnostics.append(
                        ParseDiagnostic(
                            "malformed_position",
                            f"Rejected malformed position in {name}: {pos_raw!r}",
                            severity="info",
                        )
                    )
                continue
            pos_name = _decode(row.get("posName"))
            hg_find = _decode(row.get("hgFindMatches"))
            ok, displayed = _exact_match_for_group(
                group_kind=kind,
                gene_symbol=gene_symbol,
                pos_name=pos_name,
                hg_find=hg_find,
            )
            if not ok or not displayed:
                continue
            # Neighbor guard: symbol-AS1 etc.
            if displayed.upper() != gene_symbol.strip().upper():
                continue
            exact_in_group += 1
            enst = None
            m_enst = _ENST_RE.search(hg_find) or _ENST_RE.search(pos_name)
            if m_enst:
                enst = m_enst.group(1)
            inv.matches.append(
                PositionalSearchMatch(
                    genome=genome,
                    group_name=name,
                    track_name=track_name,
                    track_description=description,
                    position=parsed["position"],
                    chromosome=parsed["chrom"],
                    display_start=parsed["start"],
                    display_end=parsed["end"],
                    matched_identifier=hg_find or None,
                    transcript_id=enst,
                    displayed_gene_symbol=displayed,
                    canonical_flag=bool(row.get("canonical")),
                    release_version=release,
                    group_kind=kind,
                )
            )

        if kind == "gencode_basic" and release is not None:
            basic_by_release[release] = exact_in_group
        elif kind == "gencode_comp" and release is not None:
            comp_by_release[release] = exact_in_group
        elif kind == "mane":
            inv.mane_exact_gene_count += exact_in_group
        elif kind == "knownGene":
            inv.known_gene_exact_gene_count += exact_in_group
        elif kind == "refseq":
            inv.refseq_exact_gene_count += exact_in_group

    if basic_by_release:
        inv.basic_release = max(basic_by_release)
        inv.basic_exact_gene_count = basic_by_release[inv.basic_release]
    if comp_by_release:
        inv.comprehensive_release = max(comp_by_release)
        inv.comprehensive_exact_gene_count = comp_by_release[inv.comprehensive_release]
    if (
        inv.basic_release is not None
        and inv.comprehensive_release is not None
        and inv.basic_release != inv.comprehensive_release
    ):
        inv.diagnostics.append(
            ParseDiagnostic(
                "gencode_release_disagreement",
                f"Basic GENCODE V{inv.basic_release} vs Comprehensive V{inv.comprehensive_release}",
            )
        )
    return inv


def parse_ucsc_int_array(raw: Any) -> list[int]:
    """Parse comma-terminated UCSC integer arrays safely."""
    if raw is None:
        return []
    if isinstance(raw, list):
        out: list[int] = []
        for item in raw:
            try:
                out.append(int(item))
            except (TypeError, ValueError):
                continue
        return out
    text = str(raw).strip()
    if not text:
        return []
    parts = [p for p in _COMMA_ARRAY_RE.split(text) if p != ""]
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except ValueError:
            continue
    return out


def _parse_tag_set(raw: Any) -> tuple[str | None, frozenset[str]]:
    if raw is None:
        return None, frozenset()
    text = str(raw).strip()
    if not text:
        return None, frozenset()
    parts = {p.strip() for p in text.split(",") if p.strip()}
    return text, frozenset(parts)


def reconstruct_exons(
    *,
    chrom_start: int,
    chrom_end: int,
    block_count: Any,
    block_sizes: Any,
    chrom_starts: Any,
) -> tuple[list[ExonInterval], bool, list[ParseDiagnostic]]:
    """Reconstruct absolute exon intervals from bigGenePred blocks."""
    diags: list[ParseDiagnostic] = []
    sizes = parse_ucsc_int_array(block_sizes)
    starts = parse_ucsc_int_array(chrom_starts)
    try:
        count = int(block_count) if block_count is not None else None
    except (TypeError, ValueError):
        count = None
        diags.append(ParseDiagnostic("malformed_exon_arrays", "blockCount is not an integer"))

    if count is None:
        return [], False, diags
    if count != len(sizes) or count != len(starts):
        diags.append(
            ParseDiagnostic(
                "malformed_exon_arrays",
                f"blockCount={count} sizes={len(sizes)} starts={len(starts)}",
            )
        )
        return [], False, diags

    exons: list[ExonInterval] = []
    prev_end = None
    for rel_start, size in zip(starts, sizes, strict=True):
        if size <= 0:
            diags.append(ParseDiagnostic("malformed_exon_arrays", f"non-positive block size {size}"))
            return [], False, diags
        if rel_start < 0:
            diags.append(ParseDiagnostic("malformed_exon_arrays", f"negative relative start {rel_start}"))
            return [], False, diags
        abs_start = chrom_start + rel_start
        abs_end = abs_start + size
        if abs_start < chrom_start or abs_end > chrom_end:
            diags.append(
                ParseDiagnostic(
                    "exon_out_of_range",
                    f"exon {abs_start}-{abs_end} outside {chrom_start}-{chrom_end}",
                )
            )
            return [], False, diags
        if prev_end is not None and abs_start < prev_end:
            diags.append(ParseDiagnostic("malformed_exon_arrays", "exons not in genomic order"))
            return [], False, diags
        exons.append(ExonInterval(start=abs_start, end=abs_end))
        prev_end = abs_end
    return exons, True, diags


def parse_transcript_row(row: dict[str, Any]) -> TranscriptRow | None:
    """Parse one knownGene / bigGenePred or legacy genePred row."""
    if not isinstance(row, dict):
        return None
    chrom = row.get("chrom")
    if row.get("chromStart") is not None:
        try:
            chrom_start = int(row["chromStart"])
            chrom_end = int(row["chromEnd"])
        except (TypeError, ValueError, KeyError):
            return None
        block_count = row.get("blockCount")
        block_sizes = row.get("blockSizes")
        chrom_starts = row.get("chromStarts")
    elif row.get("txStart") is not None:
        try:
            chrom_start = int(row["txStart"])
            chrom_end = int(row["txEnd"])
        except (TypeError, ValueError, KeyError):
            return None
        # Legacy genePred: absolute exonStarts/exonEnds
        exon_starts = parse_ucsc_int_array(row.get("exonStarts"))
        exon_ends = parse_ucsc_int_array(row.get("exonEnds"))
        block_count = len(exon_starts)
        block_sizes = [e - s for s, e in zip(exon_starts, exon_ends, strict=False)]
        chrom_starts = [s - chrom_start for s in exon_starts]
    else:
        return None

    exons, ok, ediags = reconstruct_exons(
        chrom_start=chrom_start,
        chrom_end=chrom_end,
        block_count=block_count,
        block_sizes=block_sizes,
        chrom_starts=chrom_starts,
    )
    tag_raw, tags = _parse_tag_set(row.get("tag"))
    tier_raw, tiers = _parse_tag_set(row.get("tier"))
    rank_val = row.get("rank")
    try:
        rank = int(rank_val) if rank_val is not None and str(rank_val).strip() != "" else None
    except (TypeError, ValueError):
        rank = None

    return TranscriptRow(
        raw=dict(row),
        chrom=str(chrom),
        chrom_start=chrom_start,
        chrom_end=chrom_end,
        name=str(row.get("name") or ""),
        strand=str(row.get("strand")) if row.get("strand") is not None else None,
        gene_name=str(row["geneName"]) if row.get("geneName") is not None else None,
        gene_name2=str(row["geneName2"]) if row.get("geneName2") is not None else None,
        gene_type=str(row["geneType"]) if row.get("geneType") is not None else None,
        transcript_class=(
            str(row["transcriptClass"]) if row.get("transcriptClass") is not None else None
        ),
        transcript_type=(
            str(row["transcriptType"]) if row.get("transcriptType") is not None else None
        ),
        source=str(row["source"]) if row.get("source") is not None else None,
        tag_raw=tag_raw,
        tags=tags,
        tier_raw=tier_raw,
        tiers=tiers,
        rank=rank,
        block_count=int(block_count) if block_count is not None else None,
        exons=exons,
        exon_ok=ok,
        exon_diagnostics=ediags,
    )


def parse_gencode_release_from_big_data_url(url: str | None) -> int | None:
    if not url:
        return None
    m = re.search(r"gencodeV(\d+)", str(url), re.IGNORECASE)
    return int(m.group(1)) if m else None


def filter_exact_gene_rows(
    rows: Sequence[TranscriptRow],
    gene_symbol: str,
) -> tuple[list[TranscriptRow], list[TranscriptRow], list[str]]:
    """Split regional rows into exact-gene vs excluded neighbors."""
    target = gene_symbol.strip().upper()
    exact: list[TranscriptRow] = []
    excluded: list[TranscriptRow] = []
    symbols: list[str] = []
    for row in rows:
        gn = (row.gene_name or "").strip()
        if not gn:
            excluded.append(row)
            continue
        if gn.upper() == target:
            exact.append(row)
        else:
            excluded.append(row)
            if gn not in symbols:
                symbols.append(gn)
    return exact, excluded, symbols


def select_canonical_transcript(
    candidates: Sequence[TranscriptRow],
) -> tuple[TranscriptRow | None, list[ParseDiagnostic]]:
    """Deterministic MANE/canonical selection. Ties emit a warning and return None."""
    diags: list[ParseDiagnostic] = []
    if not candidates:
        diags.append(ParseDiagnostic("missing_exact_gene_transcripts", "No exact-gene rows"))
        return None, diags

    def score(row: TranscriptRow) -> tuple:
        rank = row.rank if row.rank is not None and row.rank > 0 else 10**9
        return (
            0 if row.is_mane_select else 1,
            0 if row.is_ensembl_canonical else 1,
            0 if row.is_gencode_primary else 1,
            0 if row.is_canonical_tier else 1,
            0 if row.is_protein_coding else 1,
            rank,
            row.name,
        )

    ranked = sorted(candidates, key=score)
    best = ranked[0]
    best_key = score(best)[:-1]  # exclude name
    tied = [r for r in ranked if score(r)[:-1] == best_key]
    if len(tied) > 1:
        diags.append(
            ParseDiagnostic(
                "ambiguous_canonical_transcript",
                "Multiple transcripts tied on all selection criteria: "
                + ", ".join(t.name for t in tied),
            )
        )
        return None, diags
    return best, diags


def parse_known_gene_region(
    payload: dict[str, Any] | None,
    *,
    gene_symbol: str,
    genome: str = "hg38",
    search_inventory: SearchInventory | None = None,
) -> RegionalInventory:
    """Parse a knownGene track response and build a regional inventory."""
    diags: list[ParseDiagnostic] = []
    if not isinstance(payload, dict):
        return RegionalInventory(
            genome=genome,
            track_big_data_url=None,
            track_data_time=None,
            current_gencode_release=None,
            release_source=None,
            rows=[],
            exact_rows=[],
            excluded_neighbor_symbols=[],
            regional_total_count=0,
            exact_gene_transcript_count=0,
            excluded_neighbor_count=0,
            malformed_row_count=0,
            diagnostics=[ParseDiagnostic("missing_track", "knownGene payload missing")],
        )

    big_url = payload.get("bigDataUrl")
    data_time = payload.get("dataTime")
    release = parse_gencode_release_from_big_data_url(
        str(big_url) if big_url is not None else None
    )
    release_source = "track_big_data_url" if release is not None else None

    search_cross: dict[str, Any] = {}
    if search_inventory is not None:
        search_cross = {
            "basic_release": search_inventory.basic_release,
            "comprehensive_release": search_inventory.comprehensive_release,
            "basic_exact_gene_count": search_inventory.basic_exact_gene_count,
            "comprehensive_exact_gene_count": search_inventory.comprehensive_exact_gene_count,
        }
        search_newest = None
        for val in (search_inventory.comprehensive_release, search_inventory.basic_release):
            if val is None:
                continue
            search_newest = val if search_newest is None else max(search_newest, val)
        if release is not None and search_newest is not None and release != search_newest:
            diags.append(
                ParseDiagnostic(
                    "gencode_release_crosscheck_disagreement",
                    f"Track bigDataUrl V{release} vs search newest V{search_newest}",
                )
            )
        if release is None and search_newest is not None:
            release = search_newest
            release_source = "search_crosscheck"
            diags.append(
                ParseDiagnostic(
                    "release_from_search_fallback",
                    "Used search GENCODE release because bigDataUrl release was absent",
                    severity="info",
                )
            )

    raw_rows = payload.get("knownGene")
    if raw_rows is None:
        # Sometimes keyed by track name only
        raw_rows = payload.get("track") if isinstance(payload.get("track"), list) else []
    if not isinstance(raw_rows, list):
        raw_rows = []

    parsed: list[TranscriptRow] = []
    malformed = 0
    for row in raw_rows:
        tr = parse_transcript_row(row) if isinstance(row, dict) else None
        if tr is None:
            malformed += 1
            continue
        if not tr.exon_ok:
            malformed += 1
            diags.extend(tr.exon_diagnostics)
        parsed.append(tr)

    exact, excluded, symbols = filter_exact_gene_rows(parsed, gene_symbol)
    return RegionalInventory(
        genome=genome,
        track_big_data_url=str(big_url) if big_url is not None else None,
        track_data_time=str(data_time) if data_time is not None else None,
        current_gencode_release=release,
        release_source=release_source,
        rows=parsed,
        exact_rows=exact,
        excluded_neighbor_symbols=symbols,
        regional_total_count=len(parsed),
        exact_gene_transcript_count=len(exact),
        excluded_neighbor_count=len(excluded),
        malformed_row_count=malformed,
        diagnostics=diags,
    )


def selection_reasons(row: TranscriptRow) -> list[str]:
    reasons: list[str] = []
    if row.is_mane_select:
        reasons.append("MANE_Select")
    if row.is_ensembl_canonical:
        reasons.append("Ensembl_canonical")
    if row.is_gencode_primary:
        reasons.append("GENCODE_Primary")
    if row.is_canonical_tier:
        reasons.append("canonical_tier")
    if row.is_protein_coding:
        reasons.append("protein_coding")
    if row.rank is not None:
        reasons.append(f"rank={row.rank}")
    return reasons


__all__ = [
    "ParseDiagnostic",
    "PositionalSearchMatch",
    "SearchInventory",
    "ExonInterval",
    "TranscriptRow",
    "RegionalInventory",
    "parse_position",
    "parse_search_response",
    "parse_ucsc_int_array",
    "reconstruct_exons",
    "parse_transcript_row",
    "parse_gencode_release_from_big_data_url",
    "filter_exact_gene_rows",
    "select_canonical_transcript",
    "parse_known_gene_region",
    "selection_reasons",
    "assembly_label",
    "interval_from_api",
    "interval_from_display",
]
