"""Tests for species-identity EvidenceRecord hygiene and Human Ensembl taxon gating."""

from __future__ import annotations

from datetime import datetime, timezone

from gene_dossier.db import (
    delete_evidence_record,
    init_db,
    list_evidence_for_run,
    save_api_run,
    save_evidence_record,
    save_raw_artifact,
    session_scope,
)
from gene_dossier.identity_hygiene import (
    assert_dossier_gene_matches,
    canonical_primary_identifier,
    dedupe_species_identity_records,
    identity_dedup_key,
    normalize_ensembl_gene_id,
    normalize_entrez_id,
    normalize_uniprot_accession,
    resolve_taxon_id,
    select_preferred_identity_record,
)
from gene_dossier.models import (
    ApiRun,
    AssertionType,
    EvidenceGrade,
    EvidenceRecord,
    RawArtifact,
    SourceType,
)
from gene_dossier.report_presentation import (
    NOT_AVAILABLE,
    build_gene_aliases_blocks,
)
from gene_dossier.source_ids import make_source_id


def _identity(
    *,
    evidence_id: str,
    gene_symbol: str,
    official_symbol: str,
    taxon_id: int,
    entrez: str,
    query_gene_symbol: str | None = None,
    species_gene_symbol: str | None = None,
    raw_artifact_id: str | None = None,
    api_run_id: str | None = None,
    created_at: datetime | None = None,
    dossier_run_id: str = "hygiene-run",
) -> EvidenceRecord:
    value = {
        "entrez_gene_id": entrez,
        "nomenclaturesymbol": official_symbol,
    }
    if query_gene_symbol is not None:
        value["query_gene_symbol"] = query_gene_symbol
    if species_gene_symbol is not None:
        value["species_gene_symbol"] = species_gene_symbol
    return EvidenceRecord(
        id=evidence_id,
        source_id=make_source_id(
            "NCBI Gene", gene_symbol, AssertionType.gene_identity, entrez
        ),
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        official_symbol=official_symbol,
        section="General gene information",
        source_name="NCBI Gene",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="entrez_gene_id",
        organism="Mus musculus" if taxon_id == 10090 else "Homo sapiens",
        taxon_id=taxon_id,
        evidence_grade=EvidenceGrade.C,
        value=value,
        display_text=f"{official_symbol} Entrez Gene ID is {entrez}.",
        raw_artifact_id=raw_artifact_id,
        api_run_id=api_run_id,
        created_at=created_at or datetime(2026, 1, 1, tzinfo=timezone.utc),
    )


def test_dedupe_removes_stale_species_symbol_identity_record():
    stale = _identity(
        evidence_id="stale-mouse",
        gene_symbol="Srebf2",
        official_symbol="Srebf2",
        taxon_id=10090,
        entrez="20788",
        created_at=datetime(2026, 1, 1, tzinfo=timezone.utc),
    )
    corrected = _identity(
        evidence_id="corrected-mouse",
        gene_symbol="SREBF2",
        official_symbol="Srebf2",
        taxon_id=10090,
        entrez="20788",
        query_gene_symbol="SREBF2",
        species_gene_symbol="Srebf2",
        raw_artifact_id="raw-1",
        api_run_id="api-1",
        created_at=datetime(2026, 2, 1, tzinfo=timezone.utc),
    )
    result = dedupe_species_identity_records(
        [stale, corrected], query_symbol="SREBF2"
    )
    assert len(result.removed) == 1
    assert result.removed[0].id == "stale-mouse"
    retained_ids = {r.id for r in result.retained}
    assert "corrected-mouse" in retained_ids
    assert "stale-mouse" not in retained_ids


def test_dedupe_retains_corrected_query_symbol_record():
    stale = _identity(
        evidence_id="old",
        gene_symbol="Srebf2",
        official_symbol="Srebf2",
        taxon_id=10090,
        entrez="20788",
    )
    corrected = _identity(
        evidence_id="new",
        gene_symbol="SREBF2",
        official_symbol="Srebf2",
        taxon_id=10090,
        entrez="20788",
        query_gene_symbol="SREBF2",
        species_gene_symbol="Srebf2",
        raw_artifact_id="raw-x",
        api_run_id="api-x",
    )
    preferred = select_preferred_identity_record(
        [stale, corrected], query_symbol="SREBF2"
    )
    assert preferred.id == "new"
    assert preferred.gene_symbol == "SREBF2"
    assert preferred.official_symbol == "Srebf2"
    assert preferred.value["query_gene_symbol"] == "SREBF2"
    assert preferred.value["species_gene_symbol"] == "Srebf2"


def test_dedupe_preserves_raw_artifacts_and_api_runs(tmp_path):
    from gene_dossier.db import ApiRunRow, RawArtifactRow, get_engine, save_dossier_run
    from gene_dossier.models import DossierRun
    from sqlmodel import select

    engine = get_engine("sqlite://")
    init_db(engine)

    run_id = "hygiene-persist-run"
    run = DossierRun(id=run_id, gene_symbol="SREBF2")
    api = ApiRun(
        id="api-keep",
        dossier_run_id=run_id,
        gene_symbol="SREBF2",
        source_name="NCBI Gene",
        endpoint_name="lookup_gene_mouse",
        request_url="https://example.test/ncbi",
        success=True,
    )
    artifact = RawArtifact(
        id="raw-keep",
        dossier_run_id=run_id,
        source_name="NCBI Gene",
        artifact_type="json",
        file_path=str(tmp_path / "raw.json"),
        content_hash="abc",
        captured_at=datetime.now(timezone.utc),
        api_run_id=api.id,
    )
    stale = _identity(
        evidence_id="stale-row",
        gene_symbol="Srebf2",
        official_symbol="Srebf2",
        taxon_id=10090,
        entrez="20788",
        dossier_run_id=run_id,
        raw_artifact_id=artifact.id,
        api_run_id=api.id,
    )
    corrected = _identity(
        evidence_id="corrected-row",
        gene_symbol="SREBF2",
        official_symbol="Srebf2",
        taxon_id=10090,
        entrez="20788",
        query_gene_symbol="SREBF2",
        species_gene_symbol="Srebf2",
        dossier_run_id=run_id,
        raw_artifact_id=artifact.id,
        api_run_id=api.id,
    )

    with session_scope(engine) as session:
        save_dossier_run(session, run)
        save_api_run(session, api)
        save_raw_artifact(session, artifact)
        save_evidence_record(session, stale)
        save_evidence_record(session, corrected)

    with session_scope(engine) as session:
        evidence = list_evidence_for_run(session, run_id)
        result = dedupe_species_identity_records(evidence, query_symbol="SREBF2")
        assert {r.id for r in result.removed} == {"stale-row"}
        for rec in result.removed:
            assert delete_evidence_record(session, rec.id)

    with session_scope(engine) as session:
        remaining = list_evidence_for_run(session, run_id)
        assert [r.id for r in remaining] == ["corrected-row"]
        assert session.get(ApiRunRow, "api-keep") is not None
        assert session.get(RawArtifactRow, "raw-keep") is not None
        assert session.exec(select(ApiRunRow)).first() is not None
        assert session.exec(select(RawArtifactRow)).first() is not None


def test_rejects_taxon_unspecified_ensembl_for_human_column():
    records = [
        EvidenceRecord(
            id="ens-unspec",
            source_id="ensembl:srebf2:gene_identity:ensg00000198911",
            dossier_run_id="r",
            gene_symbol="SREBF2",
            official_symbol="SREBF2",
            section="General gene information",
            source_name="Ensembl",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ensembl_gene_id",
            organism=None,
            taxon_id=None,
            evidence_grade=EvidenceGrade.C,
            value={"ensembl_gene_id": "ENSG00000198911", "display_name": "SREBF2"},
            display_text="SREBF2 Ensembl gene ID is ENSG00000198911.",
        )
    ]
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2", evidence_records=records
    ).blocks[0]
    assert block.table_rows[3][1] == NOT_AVAILABLE


def test_accepts_explicit_homo_sapiens_legacy_ensembl_evidence():
    records = [
        EvidenceRecord(
            id="ens-legacy",
            source_id="ensembl:srebf2:gene_identity:ensg00000198911",
            dossier_run_id="r",
            gene_symbol="SREBF2",
            official_symbol="SREBF2",
            section="General gene information",
            source_name="Ensembl",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ensembl_gene_id",
            organism="Homo sapiens",
            taxon_id=None,
            evidence_grade=EvidenceGrade.C,
            value={"ensembl_gene_id": "ENSG00000198911", "display_name": "SREBF2"},
            display_text="SREBF2 Ensembl gene ID is ENSG00000198911.",
        )
    ]
    block = build_gene_aliases_blocks(
        gene_symbol="SREBF2", evidence_records=records
    ).blocks[0]
    assert "ENSG00000198911" in block.table_rows[3][1]


def test_exact_final_diagnostic_behavior_for_rat_ensembl_xref_conflict():
    records = [
        EvidenceRecord(
            id="ens-rat",
            source_id="ensembl:srebf2:gene_identity:ensrnog00000007400",
            dossier_run_id="r",
            gene_symbol="SREBF2",
            official_symbol="Srebf2",
            section="General gene information",
            source_name="Ensembl",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ensembl_gene_id",
            organism="Rattus norvegicus",
            taxon_id=10116,
            evidence_grade=EvidenceGrade.C,
            value={
                "ensembl_gene_id": "ENSRNOG00000007400",
                "display_name": "Srebf2",
                "query_gene_symbol": "SREBF2",
                "species_gene_symbol": "Srebf2",
            },
            display_text="Srebf2 Ensembl gene ID is ENSRNOG00000007400.",
        ),
        EvidenceRecord(
            id="up-xref",
            source_id="uniprot:srebf2:gene_identity:q3t1i5-ensrnog00055029130",
            dossier_run_id="r",
            gene_symbol="SREBF2",
            official_symbol="Srebf2",
            section="General gene information",
            source_name="UniProt",
            source_type=SourceType.curated_database,
            assertion_type=AssertionType.gene_identity,
            fact_type="ensembl_xref",
            organism="Rattus norvegicus",
            taxon_id=10116,
            evidence_grade=EvidenceGrade.C,
            value={
                "ensembl_gene_id": "ENSRNOG00055029130",
                "uniprot_accession": "Q3T1I5",
                "primary": False,
                "query_gene_symbol": "SREBF2",
                "species_gene_symbol": "Srebf2",
            },
            display_text=(
                "Srebf2 UniProt Q3T1I5 cross-references Ensembl ENSRNOG00055029130."
            ),
        ),
    ]
    result = build_gene_aliases_blocks(
        gene_symbol="SREBF2", evidence_records=records
    )
    assert "ENSRNOG00000007400" in result.blocks[0].table_rows[3][3]
    assert "ENSRNOG00055029130" not in result.blocks[0].table_rows[3][3]
    warnings = [
        d
        for d in result.diagnostics
        if d.severity == "warning" and d.field == "rat.ensembl"
    ]
    assert len(warnings) == 1
    diag = warnings[0]
    assert "kept direct Ensembl='ENSRNOG00000007400'" in diag.reason
    assert "ENSRNOG00055029130" in diag.reason
    # Must not claim a clean diagnostic set when a source conflict exists.
    assert any(d.severity == "warning" for d in result.diagnostics)


def test_normalize_entrez_strips_leading_zeros():
    assert normalize_entrez_id("006721") == "6721"
    assert normalize_entrez_id(6721) == "6721"


def test_normalize_ensembl_drops_version_and_uppercases():
    assert normalize_ensembl_gene_id("ensg00000198911.14") == "ENSG00000198911"
    assert normalize_ensembl_gene_id("ENSG00000198911") == "ENSG00000198911"
    assert normalize_ensembl_gene_id("ENST00000389809.8") is None


def test_normalize_uniprot_uppercases():
    assert normalize_uniprot_accession("q12772") == "Q12772"


def test_versioned_and_unversioned_ensembl_ids_dedupe_together():
    older = EvidenceRecord(
        id="ens-unversioned",
        source_id="ensembl:srebf2:gene_identity:ensg00000198911",
        dossier_run_id="hygiene-run",
        gene_symbol="Srebf2",
        official_symbol="SREBF2",
        section="General gene information",
        source_name="Ensembl",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="ensembl_gene_id",
        organism="Homo sapiens",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.C,
        value={"ensembl_gene_id": "ENSG00000198911"},
        display_text="SREBF2 Ensembl gene ID is ENSG00000198911.",
        created_at=datetime(2026, 1, 1, tzinfo=timezone.utc),
    )
    newer = EvidenceRecord(
        id="ens-versioned",
        source_id="ensembl:srebf2:gene_identity:ensg00000198911-v",
        dossier_run_id="hygiene-run",
        gene_symbol="SREBF2",
        official_symbol="SREBF2",
        section="General gene information",
        source_name="Ensembl",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="ensembl_gene_id",
        organism="Homo sapiens",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.C,
        value={
            "ensembl_gene_id": "ENSG00000198911.14",
            "query_gene_symbol": "SREBF2",
            "species_gene_symbol": "SREBF2",
        },
        display_text="SREBF2 Ensembl gene ID is ENSG00000198911.14.",
        raw_artifact_id="raw-ens",
        api_run_id="api-ens",
        created_at=datetime(2026, 2, 1, tzinfo=timezone.utc),
    )
    assert canonical_primary_identifier(older) == canonical_primary_identifier(newer)
    assert identity_dedup_key(older) == identity_dedup_key(newer)
    result = dedupe_species_identity_records(
        [older, newer], query_symbol="SREBF2"
    )
    assert {r.id for r in result.removed} == {"ens-unversioned"}
    retained = [r for r in result.retained if r.id == "ens-versioned"]
    assert len(retained) == 1


def test_legacy_taxon_only_in_value_groups_with_explicit_taxon():
    legacy = EvidenceRecord(
        id="ens-legacy-tax",
        source_id="ensembl:srebf2:gene_identity:ensg-legacy",
        dossier_run_id="hygiene-run",
        gene_symbol="Srebf2",
        official_symbol="SREBF2",
        section="General gene information",
        source_name="Ensembl",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="ensembl_gene_id",
        organism=None,
        taxon_id=None,
        evidence_grade=EvidenceGrade.C,
        value={
            "ensembl_gene_id": "ENSG00000198911",
            "taxon_id": 9606,
        },
        display_text="SREBF2 Ensembl gene ID is ENSG00000198911.",
    )
    current = EvidenceRecord(
        id="ens-current-tax",
        source_id="ensembl:srebf2:gene_identity:ensg-current",
        dossier_run_id="hygiene-run",
        gene_symbol="SREBF2",
        official_symbol="SREBF2",
        section="General gene information",
        source_name="Ensembl",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="ensembl_gene_id",
        organism="Homo sapiens",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.C,
        value={
            "ensembl_gene_id": "ENSG00000198911",
            "query_gene_symbol": "SREBF2",
            "species_gene_symbol": "SREBF2",
        },
        display_text="SREBF2 Ensembl gene ID is ENSG00000198911.",
        raw_artifact_id="raw-tax",
        api_run_id="api-tax",
    )
    assert resolve_taxon_id(legacy) == 9606
    assert identity_dedup_key(legacy) == identity_dedup_key(current)
    result = dedupe_species_identity_records(
        [legacy, current], query_symbol="SREBF2"
    )
    assert {r.id for r in result.removed} == {"ens-legacy-tax"}


def test_genomic_location_key_includes_assembly():
    loc_a = EvidenceRecord(
        id="loc-a",
        source_id="ensembl:srebf2:gene_identity:loc-a",
        dossier_run_id="hygiene-run",
        gene_symbol="SREBF2",
        official_symbol="SREBF2",
        section="General gene information",
        source_name="Ensembl",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="genomic_location",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.C,
        value={
            "seq_region_name": "22",
            "start": 41833354,
            "end": 41903307,
            "assembly_name": "GRCh38",
        },
        display_text="loc a",
    )
    loc_b = EvidenceRecord(
        id="loc-b",
        source_id="ensembl:srebf2:gene_identity:loc-b",
        dossier_run_id="hygiene-run",
        gene_symbol="SREBF2",
        official_symbol="SREBF2",
        section="General gene information",
        source_name="Ensembl",
        source_type=SourceType.curated_database,
        assertion_type=AssertionType.gene_identity,
        fact_type="genomic_location",
        taxon_id=9606,
        evidence_grade=EvidenceGrade.C,
        value={
            "seq_region_name": "22",
            "start": 41833354,
            "end": 41903307,
            "assembly_name": "GRCh37",
        },
        display_text="loc b",
    )
    assert canonical_primary_identifier(loc_a) == "GRCh38:22:41833354-41903307"
    assert canonical_primary_identifier(loc_b) == "GRCh37:22:41833354-41903307"
    assert identity_dedup_key(loc_a) != identity_dedup_key(loc_b)
    result = dedupe_species_identity_records(
        [loc_a, loc_b], query_symbol="SREBF2"
    )
    assert result.removed == ()


def test_dedupe_script_aborts_on_gene_mismatch():
    from gene_dossier.identity_hygiene import assert_dossier_gene_matches

    assert_dossier_gene_matches("SREBF2", "SREBF2")
    assert_dossier_gene_matches("srebf2", "SREBF2")
    try:
        assert_dossier_gene_matches("SREBF2", "CDH10")
        raise AssertionError("expected ValueError on gene mismatch")
    except ValueError as exc:
        assert "Gene mismatch" in str(exc)
        assert "SREBF2" in str(exc)
        assert "CDH10" in str(exc)
