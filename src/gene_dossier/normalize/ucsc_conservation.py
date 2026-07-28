"""Normalize UCSC conservation evidence for Section 1b.

Primary UCSC conservation normalizer. Gene-identity keeps only a thin
compatibility wrapper that delegates here when needed.
"""

from __future__ import annotations

from typing import Any

from gene_dossier.models import (
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    SourceType,
    ToolResult,
)
from gene_dossier.source_ids import make_source_id
from gene_dossier.ucsc_coords import assembly_label, interval_from_api
from gene_dossier.ucsc_parse import (
    parse_known_gene_region,
    parse_search_response,
    select_canonical_transcript,
    selection_reasons,
)
from gene_dossier.ucsc_figure import build_safe_hgtracks_url

SOURCE_NAME = "UCSC"
SECTION_GENERAL = "General Gene Information"
SUBSECTION = "UCSC genome browser"


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    assertion_type: AssertionType = AssertionType.gene_identity,
    evidence_grade: EvidenceGrade = EvidenceGrade.C,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    source_id = make_source_id(SOURCE_NAME, gene_symbol, assertion_type, fact_type, key)
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=gene_symbol,
        section=SECTION_GENERAL,
        subsection=SUBSECTION,
        source_name=SOURCE_NAME,
        source_type=SourceType.curated_database,
        assertion_type=assertion_type,
        fact_type=fact_type,
        organism="Homo sapiens",
        species="human",
        taxon_id=9606,
        evidence_grade=evidence_grade,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def build_conservation_evidence(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    genome: str,
    search_payload: dict[str, Any] | None,
    track_payload: dict[str, Any] | None,
    figure_value: dict[str, Any] | None = None,
    search_api_run_id: str | None = None,
    search_raw_artifact_id: str | None = None,
    track_api_run_id: str | None = None,
    track_raw_artifact_id: str | None = None,
    figure_api_run_id: str | None = None,
    figure_raw_artifact_id: str | None = None,
) -> tuple[list[EvidenceRecord], list[dict[str, str]]]:
    """Build the four UCSC Section 1b EvidenceRecords from parsed payloads."""
    diagnostics: list[dict[str, str]] = []
    records: list[EvidenceRecord] = []

    search_inv = parse_search_response(search_payload, gene_symbol=gene_symbol, genome=genome)
    for d in search_inv.diagnostics:
        diagnostics.append({"code": d.code, "message": d.message, "severity": d.severity})

    region = parse_known_gene_region(
        track_payload,
        gene_symbol=gene_symbol,
        genome=genome,
        search_inventory=search_inv,
    )
    for d in region.diagnostics:
        diagnostics.append({"code": d.code, "message": d.message, "severity": d.severity})

    selected, sel_diags = select_canonical_transcript(region.exact_rows)
    for d in sel_diags:
        diagnostics.append({"code": d.code, "message": d.message, "severity": d.severity})

    asm = assembly_label(genome)
    display_interval = search_inv.selected_display_interval
    if selected is not None:
        api_interval = interval_from_api(
            selected.chrom,
            selected.chrom_start,
            selected.chrom_end,
            genome=genome,
        )
    elif display_interval is not None:
        api_interval = display_interval
    else:
        api_interval = None
        diagnostics.append(
            {
                "code": "missing_assembly_locus",
                "message": "Could not resolve a UCSC locus interval",
                "severity": "warning",
            }
        )

    if api_interval is not None:
        locus_value = {
            **api_interval.as_dict(),
            "assembly_display_name": asm["assembly_display_name"],
            "assembly_label_source": asm["assembly_label_source"],
            "requested_gene_symbol": gene_symbol,
        }
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                fact_type="ucsc_gene_locus",
                key=f"{genome}:{api_interval.display_position}",
                value=locus_value,
                display_text=(
                    f"{gene_symbol} UCSC/{genome} locus "
                    f"{api_interval.display_position}."
                ),
                api_run_id=track_api_run_id or search_api_run_id,
                raw_artifact_id=track_raw_artifact_id or search_raw_artifact_id,
            )
        )

    # Count disagreement warning
    if (
        region.exact_gene_transcript_count
        and search_inv.comprehensive_exact_gene_count
        and region.exact_gene_transcript_count != search_inv.comprehensive_exact_gene_count
    ):
        diagnostics.append(
            {
                "code": "transcript_count_disagreement",
                "message": (
                    f"Regional exact count {region.exact_gene_transcript_count} vs "
                    f"Comprehensive search count {search_inv.comprehensive_exact_gene_count}"
                ),
                "severity": "warning",
            }
        )

    release_label = (
        f"V{region.current_gencode_release}"
        if region.current_gencode_release is not None
        else None
    )
    inventory_value = {
        "genome": genome,
        "current_gencode_release": release_label,
        "current_gencode_release_number": region.current_gencode_release,
        "release_source": region.release_source,
        "track_big_data_url": region.track_big_data_url,
        "track_data_time": region.track_data_time,
        "search_release_crosscheck": {
            "basic_release": search_inv.basic_release,
            "basic_exact_gene_count": search_inv.basic_exact_gene_count,
            "comprehensive_release": search_inv.comprehensive_release,
            "comprehensive_exact_gene_count": search_inv.comprehensive_exact_gene_count,
        },
        "basic_release": (
            f"V{search_inv.basic_release}" if search_inv.basic_release is not None else None
        ),
        "basic_exact_gene_count": search_inv.basic_exact_gene_count,
        "comprehensive_release": (
            f"V{search_inv.comprehensive_release}"
            if search_inv.comprehensive_release is not None
            else None
        ),
        "comprehensive_exact_gene_count": search_inv.comprehensive_exact_gene_count,
        "regional_total_count": region.regional_total_count,
        "exact_gene_transcript_count": region.exact_gene_transcript_count,
        "excluded_neighbor_count": region.excluded_neighbor_count,
        "excluded_neighbor_symbols": list(region.excluded_neighbor_symbols),
        "malformed_row_count": region.malformed_row_count,
    }
    records.append(
        _record(
            dossier_run_id=dossier_run_id,
            gene_symbol=gene_symbol,
            fact_type="ucsc_transcript_inventory",
            key=f"{genome}:{release_label or 'unknown'}",
            value=inventory_value,
            display_text=(
                f"{gene_symbol} has {region.exact_gene_transcript_count} exact-gene "
                f"transcript models in the selected UCSC locus "
                f"({release_label or 'unknown GENCODE'})."
            ),
            api_run_id=track_api_run_id or search_api_run_id,
            raw_artifact_id=track_raw_artifact_id or search_raw_artifact_id,
        )
    )

    if selected is not None and api_interval is not None:
        browser_url = build_safe_hgtracks_url(
            genome=genome,
            display_position=api_interval.display_position,
            transcript_id=selected.name or None,
        )
        transcript_value = {
            "transcript_id": selected.name,
            "requested_gene_symbol": gene_symbol,
            "source_gene_symbol": selected.gene_name,
            "chromosome": selected.chrom,
            "strand": selected.strand,
            "genome": genome,
            "assembly_display_name": asm["assembly_display_name"],
            "assembly_label_source": asm["assembly_label_source"],
            "api_start_0_based": selected.chrom_start,
            "api_end_exclusive": selected.chrom_end,
            "display_start_1_based": api_interval.display_start_1_based,
            "display_end_1_based": api_interval.display_end_1_based,
            "display_position": api_interval.display_position,
            "coordinate_system": api_interval.coordinate_system,
            "exon_count": len(selected.exons) if selected.exon_ok else selected.block_count,
            "exons": [{"start": e.start, "end": e.end} for e in selected.exons],
            "exon_validation_ok": selected.exon_ok,
            "gene_type": selected.gene_type,
            "transcript_class": selected.transcript_class,
            "transcript_type": selected.transcript_type,
            "source": selected.source,
            "tag": selected.tag_raw,
            "tags": sorted(selected.tags),
            "tier": selected.tier_raw,
            "tiers": sorted(selected.tiers),
            "rank": selected.rank,
            "is_mane_select": selected.is_mane_select,
            "is_ensembl_canonical": selected.is_ensembl_canonical,
            "is_gencode_primary": selected.is_gencode_primary,
            "is_canonical_tier": selected.is_canonical_tier,
            "selection_reasons": selection_reasons(selected),
            "browser_url": browser_url,
            "raw_row": selected.raw,
        }
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                fact_type="ucsc_canonical_transcript",
                key=f"{genome}:{selected.name}",
                value=transcript_value,
                display_text=(
                    f"{gene_symbol} canonical transcript {selected.name} at "
                    f"{api_interval.display_position}."
                ),
                api_run_id=track_api_run_id,
                raw_artifact_id=track_raw_artifact_id,
            )
        )
    elif not region.exact_rows:
        diagnostics.append(
            {
                "code": "missing_exact_gene_transcript_rows",
                "message": "No exact-gene knownGene rows after filtering",
                "severity": "warning",
            }
        )

    if figure_value:
        fig_hash = str(figure_value.get("sha256") or figure_value.get("content_hash") or "figure")
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                fact_type="ucsc_conservation_figure",
                key=f"{genome}:{fig_hash[:16]}",
                value=dict(figure_value),
                display_text=(
                    f"UCSC conservation figure for {gene_symbol} "
                    f"({figure_value.get('display_position') or 'locus'})."
                ),
                assertion_type=AssertionType.visual_observation,
                evidence_grade=EvidenceGrade.C,
                api_run_id=figure_api_run_id,
                raw_artifact_id=figure_raw_artifact_id,
            )
        )
    else:
        diagnostics.append(
            {
                "code": "missing_figure",
                "message": "No validated UCSC conservation figure provided",
                "severity": "warning",
            }
        )

    return records, diagnostics


def normalize_ucsc_conservation(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize a combined UCSC conservation ToolResult payload."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    genome = str(data.get("genome") or "hg38")
    search_payload = data.get("search") if isinstance(data.get("search"), dict) else None
    track_payload = data.get("track_data") if isinstance(data.get("track_data"), dict) else None
    figure_value = data.get("figure") if isinstance(data.get("figure"), dict) else None
    figure_api_run_id = (
        str(figure_value.get("figure_api_run_id"))
        if figure_value and figure_value.get("figure_api_run_id")
        else api_run_id
    )
    figure_raw_artifact_id = (
        str(figure_value.get("figure_raw_artifact_id"))
        if figure_value and figure_value.get("figure_raw_artifact_id")
        else raw_artifact_id
    )
    records, _diags = build_conservation_evidence(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        genome=genome,
        search_payload=search_payload,
        track_payload=track_payload,
        figure_value=figure_value,
        search_api_run_id=api_run_id,
        search_raw_artifact_id=raw_artifact_id,
        track_api_run_id=api_run_id,
        track_raw_artifact_id=raw_artifact_id,
        figure_api_run_id=figure_api_run_id,
        figure_raw_artifact_id=figure_raw_artifact_id,
    )
    return records


__all__ = [
    "build_conservation_evidence",
    "normalize_ucsc_conservation",
]
