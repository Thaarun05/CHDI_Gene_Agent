"""Normalize protein / structure ToolResults into EvidenceRecords.

Consumes successful client payloads from UniProt, AlphaFold, PDBe, and CDD.
Does **not** call the network.

Rules:
- Prefer the UniProt client's ``selected_accession`` entry only
- Extract function / location / disease / domain facts from payload fields only
- Do not invent structures, domains, or disease relationships
- AlphaFold predictions are graded as computational (E)
- PDBe experimental mappings are curated structure evidence (C)
- CDD hits are computational domain matches (E)
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

SECTION_GENERAL = "General gene information"
SECTION_STRUCTURE = "Known structure / domains"
SECTION_STRUCTURE_BUNDLE = "AlphaFold / PDBe / CDD"

# UniProt commentType / feature type strings observed in REST JSON.
_FUNCTION_TYPES = {"FUNCTION"}
_LOCATION_TYPES = {"SUBCELLULAR LOCATION", "SUBCELLULAR_LOCATION"}
_DISEASE_TYPES = {"DISEASE"}
_DOMAIN_FEATURE_TYPES = {"Domain", "Repeat", "Region", "Motif", "Zinc finger"}


def _as_dict(data: Any) -> dict[str, Any]:
    return data if isinstance(data, dict) else {}


def _record(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    source_name: str,
    source_type: SourceType,
    assertion_type: AssertionType,
    fact_type: str,
    key: str,
    value: dict[str, Any],
    display_text: str,
    evidence_grade: EvidenceGrade,
    section: str,
    organism: str | None = None,
    taxon_id: int | None = None,
    official_symbol: str | None = None,
    subsection: str | None = None,
    confidence_notes: str | None = None,
    manual_review_required: bool = False,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> EvidenceRecord:
    """Build one EvidenceRecord with a deterministic source_id."""
    source_id = make_source_id(source_name, gene_symbol, assertion_type, key)
    return EvidenceRecord(
        source_id=source_id,
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=official_symbol or gene_symbol,
        section=section,
        subsection=subsection,
        source_name=source_name,
        source_type=source_type,
        assertion_type=assertion_type,
        fact_type=fact_type,
        organism=organism,
        taxon_id=taxon_id,
        evidence_grade=evidence_grade,
        manual_review_required=manual_review_required,
        confidence_notes=confidence_notes,
        value=value,
        display_text=display_text,
        api_run_id=api_run_id,
        raw_artifact_id=raw_artifact_id,
    )


def _comment_type(comment: dict[str, Any]) -> str:
    return str(comment.get("commentType") or comment.get("type") or "").strip().upper()


def _comment_texts(comment: dict[str, Any]) -> list[str]:
    """Collect free-text values from a UniProt comment block."""
    texts: list[str] = []
    for item in comment.get("texts") or []:
        if isinstance(item, dict) and item.get("value"):
            texts.append(str(item["value"]).strip())
        elif isinstance(item, str) and item.strip():
            texts.append(item.strip())
    # Some comment shapes use a single ``value`` / ``note``.
    if not texts and comment.get("value"):
        texts.append(str(comment["value"]).strip())
    note = comment.get("note")
    if isinstance(note, dict) and note.get("texts"):
        for item in note["texts"]:
            if isinstance(item, dict) and item.get("value"):
                texts.append(str(item["value"]).strip())
    return [t for t in texts if t]


def _subcellular_locations(comment: dict[str, Any]) -> list[str]:
    """Extract location labels from a SUBCELLULAR LOCATION comment."""
    locations: list[str] = []
    for row in comment.get("subcellularLocations") or []:
        if not isinstance(row, dict):
            continue
        loc = row.get("location") or {}
        if isinstance(loc, dict) and loc.get("value"):
            locations.append(str(loc["value"]).strip())
        elif isinstance(loc, str) and loc.strip():
            locations.append(loc.strip())
        topo = row.get("topology") or {}
        if isinstance(topo, dict) and topo.get("value"):
            locations.append(f"topology:{topo['value']}")
    if not locations:
        locations.extend(_comment_texts(comment))
    return [x for x in locations if x]


def _disease_texts(comment: dict[str, Any]) -> list[str]:
    """Extract disease annotation text without inventing phenotype links."""
    texts = _comment_texts(comment)
    disease = comment.get("disease") or {}
    if isinstance(disease, dict):
        if disease.get("diseaseId"):
            texts.insert(0, f"diseaseId:{disease['diseaseId']}")
        desc = disease.get("description")
        if desc:
            texts.append(str(desc).strip())
        name = disease.get("diseaseName") or disease.get("diseaseAcronym")
        if name:
            texts.insert(0, str(name).strip())
    return [t for t in texts if t]


def _feature_span(feature: dict[str, Any]) -> tuple[Any, Any]:
    location = feature.get("location") or {}
    if not isinstance(location, dict):
        return None, None
    start = location.get("start")
    end = location.get("end")
    start_v = start.get("value") if isinstance(start, dict) else start
    end_v = end.get("value") if isinstance(end, dict) else end
    return start_v, end_v


def _entry_accession(entry: dict[str, Any]) -> str | None:
    accession = entry.get("accession") or entry.get("primaryAccession")
    return str(accession) if accession else None


def _selected_uniprot_entry(data: dict[str, Any]) -> dict[str, Any] | None:
    """Return the selected summarized UniProt entry, if available."""
    selected = data.get("selected_accession")
    selected_str = str(selected) if selected else None

    selected_entry = data.get("selected_entry")
    if isinstance(selected_entry, dict):
        entry_acc = _entry_accession(selected_entry)
        if selected_str is None or entry_acc == selected_str:
            return selected_entry

    entries = data.get("entries") or []
    if not isinstance(entries, list):
        entries = []

    if selected_str:
        for entry in entries:
            if isinstance(entry, dict) and _entry_accession(entry) == selected_str:
                return entry
    elif len(entries) == 1 and isinstance(entries[0], dict):
        return entries[0]

    # Top-level single-entry fallback (accession / primaryAccession on data itself).
    top_acc = _entry_accession(data)
    if top_acc and (selected_str is None or top_acc == selected_str):
        return data

    return None


def normalize_uniprot_protein(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize UniProt search payloads into protein function/structure records.

    Gene-identity accession records belong in ``gene_identity.py``; this module
    emits function, localization, curated disease comments, and domain features.
    """
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    entry = _selected_uniprot_entry(data)
    if not entry:
        return []

    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    accession = _entry_accession(entry) or ""
    if not gene_symbol or not accession:
        return []

    protein_name = entry.get("protein_name")
    organism = entry.get("organism_name")
    organism_id = entry.get("organism_id")
    try:
        taxon_id = int(organism_id) if organism_id is not None else None
    except (TypeError, ValueError):
        taxon_id = None

    gene_names = entry.get("gene_names") or []
    official = gene_symbol
    if isinstance(gene_names, list) and gene_names and isinstance(gene_names[0], str):
        official = gene_names[0]

    common_kwargs = {
        "dossier_run_id": dossier_run_id,
        "gene_symbol": gene_symbol,
        "official_symbol": str(official),
        "source_name": "UniProt",
        "source_type": SourceType.curated_database,
        "organism": str(organism) if organism else None,
        "taxon_id": taxon_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }

    records: list[EvidenceRecord] = []

    if protein_name:
        records.append(
            _record(
                **common_kwargs,
                assertion_type=AssertionType.protein_function,
                fact_type="protein_name",
                key=f"{accession}-name",
                value={
                    "uniprot_accession": accession,
                    "protein_name": protein_name,
                },
                display_text=f"{official} encodes {protein_name} (UniProt {accession}).",
                evidence_grade=EvidenceGrade.C,
                section=SECTION_GENERAL,
                subsection="UniProt protein",
            )
        )

    comments = entry.get("comments") if isinstance(entry.get("comments"), list) else []
    function_idx = 0
    location_idx = 0
    disease_idx = 0
    for comment in comments:
        if not isinstance(comment, dict):
            continue
        ctype = _comment_type(comment)
        if ctype in _FUNCTION_TYPES:
            for text in _comment_texts(comment):
                function_idx += 1
                records.append(
                    _record(
                        **common_kwargs,
                        assertion_type=AssertionType.protein_function,
                        fact_type="protein_function",
                        key=f"{accession}-function-{function_idx}",
                        value={
                            "uniprot_accession": accession,
                            "comment_type": "FUNCTION",
                            "text": text,
                        },
                        display_text=f"{official} UniProt function: {text}",
                        evidence_grade=EvidenceGrade.C,
                        section=SECTION_GENERAL,
                        subsection="UniProt function",
                    )
                )
        elif ctype in _LOCATION_TYPES:
            for loc in _subcellular_locations(comment):
                location_idx += 1
                records.append(
                    _record(
                        **common_kwargs,
                        assertion_type=AssertionType.protein_function,
                        fact_type="subcellular_location",
                        key=f"{accession}-location-{location_idx}",
                        value={
                            "uniprot_accession": accession,
                            "comment_type": "SUBCELLULAR LOCATION",
                            "location": loc,
                        },
                        display_text=(
                            f"{official} UniProt subcellular location: {loc}."
                        ),
                        evidence_grade=EvidenceGrade.C,
                        section=SECTION_GENERAL,
                        subsection="UniProt localization",
                    )
                )
        elif ctype in _DISEASE_TYPES:
            for text in _disease_texts(comment):
                disease_idx += 1
                records.append(
                    _record(
                        **common_kwargs,
                        assertion_type=AssertionType.disease_association,
                        fact_type="disease_annotation",
                        key=f"{accession}-disease-{disease_idx}",
                        value={
                            "uniprot_accession": accession,
                            "comment_type": "DISEASE",
                            "text": text,
                            "caveat": (
                                "UniProt curated disease comment; not proof of "
                                "HD relevance."
                            ),
                        },
                        display_text=(
                            f"{official} UniProt disease annotation: {text}"
                        ),
                        evidence_grade=EvidenceGrade.C,
                        section=SECTION_GENERAL,
                        subsection="UniProt disease annotation",
                        confidence_notes=(
                            "Curated UniProt disease comment; not proof of HD relevance."
                        ),
                        manual_review_required=True,
                    )
                )

    features = entry.get("features") if isinstance(entry.get("features"), list) else []
    for feature in features:
        if not isinstance(feature, dict):
            continue
        ftype = str(feature.get("type") or "").strip()
        if ftype not in _DOMAIN_FEATURE_TYPES:
            continue
        start, end = _feature_span(feature)
        description = feature.get("description") or ftype
        span = ""
        if start is not None and end is not None:
            span = f"{start}-{end}"
        key_bits = [accession, ftype.lower().replace(" ", "-"), span or "na", str(description)]
        key = "-".join(str(b) for b in key_bits if b)
        records.append(
            _record(
                **common_kwargs,
                assertion_type=AssertionType.protein_structure,
                fact_type="sequence_feature",
                key=key,
                value={
                    "uniprot_accession": accession,
                    "feature_type": ftype,
                    "description": description,
                    "start": start,
                    "end": end,
                },
                display_text=(
                    f"{official} UniProt {ftype.lower()}"
                    + (f" {description}" if description else "")
                    + (f" at residues {span}." if span else ".")
                ),
                evidence_grade=EvidenceGrade.C,
                section=SECTION_STRUCTURE,
                subsection="UniProt domains / features",
            )
        )

    return records


def normalize_alphafold(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize AlphaFold ``fetch_prediction`` payloads (predicted structures)."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    accession = data.get("uniprot_accession")
    summaries = data.get("prediction_summaries") or []
    if not isinstance(summaries, list):
        summaries = []
    if not gene_symbol or not summaries:
        return []

    entry_url = data.get("entry_url")
    records: list[EvidenceRecord] = []
    for summary in summaries:
        if not isinstance(summary, dict):
            continue
        entry_id = summary.get("entry_id") or summary.get("model_entity_id")
        uniprot_acc = summary.get("uniprot_accession") or accession
        if not entry_id and not uniprot_acc:
            continue
        key = str(entry_id or uniprot_acc)
        organism = summary.get("organism_scientific_name")
        tax_id = summary.get("tax_id")
        try:
            taxon_id = int(tax_id) if tax_id is not None else None
        except (TypeError, ValueError):
            taxon_id = None
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="AlphaFold",
                source_type=SourceType.structure_database,
                assertion_type=AssertionType.protein_structure,
                fact_type="alphafold_prediction",
                key=key,
                value={
                    "uniprot_accession": uniprot_acc,
                    "entry_id": summary.get("entry_id"),
                    "model_entity_id": summary.get("model_entity_id"),
                    "global_metric_value": summary.get("global_metric_value"),
                    "fraction_plddt_very_high": summary.get("fraction_plddt_very_high"),
                    "fraction_plddt_confident": summary.get("fraction_plddt_confident"),
                    "latest_version": summary.get("latest_version"),
                    "pdb_url": summary.get("pdb_url"),
                    "cif_url": summary.get("cif_url"),
                    "entry_url": summary.get("entry_url") or entry_url,
                    "caveat": "AlphaFold structure is a prediction, not an experimental model.",
                },
                display_text=(
                    f"{gene_symbol} has an AlphaFold predicted structure"
                    + (f" for UniProt {uniprot_acc}" if uniprot_acc else "")
                    + (f" (entry {entry_id})." if entry_id else ".")
                ),
                evidence_grade=EvidenceGrade.E,
                section=SECTION_STRUCTURE_BUNDLE,
                subsection="AlphaFold prediction",
                organism=str(organism) if organism else None,
                taxon_id=taxon_id,
                confidence_notes=(
                    "Predicted structure from AlphaFold DB; not an experimental model."
                ),
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_pdbe(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize PDBe ``fetch_structures`` best-structure summaries."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    accession = data.get("uniprot_accession")
    summaries = data.get("structure_summaries") or []
    if not isinstance(summaries, list):
        summaries = []
    if not gene_symbol or not summaries:
        return []

    entry_by_pdb: dict[str, dict[str, Any]] = {}
    for entry in data.get("entry_summaries") or []:
        if isinstance(entry, dict) and entry.get("pdb_id"):
            entry_by_pdb[str(entry["pdb_id"]).lower()] = entry

    records: list[EvidenceRecord] = []
    for row in summaries:
        if not isinstance(row, dict):
            continue
        pdb_id = row.get("pdb_id")
        if not pdb_id:
            continue
        chain = row.get("chain_id")
        key = f"{pdb_id}-{chain}" if chain else str(pdb_id)
        entry_meta = entry_by_pdb.get(str(pdb_id).lower()) or {}
        title = entry_meta.get("title")
        method = row.get("experimental_method") or entry_meta.get("experimental_method")
        resolution = row.get("resolution")
        display = (
            f"{gene_symbol} maps to PDBe entry {pdb_id}"
            + (f" chain {chain}" if chain else "")
        )
        if method:
            display += f" ({method}"
            if resolution is not None:
                display += f", {resolution} Å"
            display += ")"
        display += "."

        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="PDBe",
                source_type=SourceType.structure_database,
                assertion_type=AssertionType.protein_structure,
                fact_type="pdb_structure",
                key=key,
                value={
                    "uniprot_accession": accession,
                    "pdb_id": pdb_id,
                    "chain_id": chain,
                    "unp_start": row.get("unp_start"),
                    "unp_end": row.get("unp_end"),
                    "coverage": row.get("coverage"),
                    "resolution": resolution,
                    "experimental_method": method,
                    "title": title,
                },
                display_text=display,
                evidence_grade=EvidenceGrade.C,
                section=SECTION_STRUCTURE_BUNDLE,
                subsection="PDBe experimental structure",
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_cdd(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Normalize CDD ``fetch_domains`` hit summaries (computational matches)."""
    if not tool_result.success:
        return []
    data = _as_dict(tool_result.data)
    gene_symbol = tool_result.gene_symbol or str(data.get("gene_symbol") or "")
    summaries = data.get("hit_summaries") or []
    if not isinstance(summaries, list):
        summaries = []
    if not gene_symbol or not summaries:
        return []

    queries = data.get("queries")
    cdsid = data.get("cdsid")
    records: list[EvidenceRecord] = []
    for idx, hit in enumerate(summaries, start=1):
        if not isinstance(hit, dict):
            continue
        domain_acc = hit.get("domain_accession")
        short_name = hit.get("domain_short_name")
        if not domain_acc and not short_name:
            continue
        start = hit.get("from_residue")
        end = hit.get("to_residue")
        key = "-".join(
            str(x)
            for x in (domain_acc or short_name, start, end, idx)
            if x is not None and str(x) != ""
        )
        label = short_name or domain_acc
        span = None
        if start is not None and end is not None:
            span = f"{start}-{end}"
        records.append(
            _record(
                dossier_run_id=dossier_run_id,
                gene_symbol=gene_symbol,
                source_name="CDD",
                source_type=SourceType.structure_database,
                assertion_type=AssertionType.protein_structure,
                fact_type="conserved_domain_hit",
                key=key,
                value={
                    "queries": queries,
                    "cdsid": cdsid,
                    "query_accession": hit.get("query_accession"),
                    "domain_accession": domain_acc,
                    "domain_short_name": short_name,
                    "domain_description": hit.get("domain_description"),
                    "from_residue": start,
                    "to_residue": end,
                    "evalue": hit.get("evalue"),
                    "bitscore": hit.get("bitscore"),
                    "superfamily": hit.get("superfamily"),
                    "caveat": "CDD hit from computational domain search; verify before citing as curated domain.",
                },
                display_text=(
                    f"{gene_symbol} CDD hit {label}"
                    + (f" ({domain_acc})" if domain_acc and short_name else "")
                    + (f" at residues {span}" if span else "")
                    + "."
                ),
                evidence_grade=EvidenceGrade.E,
                section=SECTION_STRUCTURE_BUNDLE,
                subsection="CDD conserved domains",
                confidence_notes=(
                    "Computational CDD domain hit; not a curated UniProt feature."
                ),
                manual_review_required=True,
                api_run_id=api_run_id,
                raw_artifact_id=raw_artifact_id,
            )
        )
    return records


def normalize_protein(
    tool_result: ToolResult,
    *,
    dossier_run_id: str,
    api_run_id: str | None = None,
    raw_artifact_id: str | None = None,
) -> list[EvidenceRecord]:
    """Dispatch protein/structure normalization by ``tool_result.source_name``."""
    source = (tool_result.source_name or "").strip()
    kwargs = {
        "dossier_run_id": dossier_run_id,
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    if source == "UniProt":
        return normalize_uniprot_protein(tool_result, **kwargs)
    if source == "AlphaFold":
        return normalize_alphafold(tool_result, **kwargs)
    if source == "PDBe":
        return normalize_pdbe(tool_result, **kwargs)
    if source == "CDD":
        return normalize_cdd(tool_result, **kwargs)
    return []


__all__ = [
    "normalize_uniprot_protein",
    "normalize_alphafold",
    "normalize_pdbe",
    "normalize_cdd",
    "normalize_protein",
]
