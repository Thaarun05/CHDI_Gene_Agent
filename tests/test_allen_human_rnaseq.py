"""Offline tests for Allen human RNA-seq donor ZIP client."""

from __future__ import annotations

import io
import zipfile

from gene_dossier.tools import allen_human_rnaseq as ahr


def _zip_bytes(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, content in members.items():
            zf.writestr(f"package/{name}", content)
    return buf.getvalue()


def _sample_annot(*, brain: str = "1", n_bio: int = 3, n_tech: int = 1) -> str:
    rows = [
        "RNAseq_sample_name,replicate_sample,sample_name,well_id,brain",
    ]
    for i in range(n_bio):
        rows.append(f"S_BIO_{i},No,S{i},{1000 + i},{brain}")
    for i in range(n_tech):
        rows.append(f"S_TECH_{i},Yes,T{i},{2000 + i},{brain}")
    return "\n".join(rows) + "\n"


def _genes_csv() -> str:
    return (
        "gene_symbol,gene_id,entrez_id\n"
        "SREBF2,6681,6721\n"
        "OTHER,1,2\n"
    )


def _tpm_csv(*, n_samples: int = 4, gene_values: list[float] | None = None) -> str:
    header = ["gene_symbol"] + [f"S_BIO_{i}" for i in range(3)] + ["S_TECH_0"]
    assert len(header) - 1 == n_samples
    # Keep values in (0, 1) so the filler gene can make each column sum to 1.0.
    vals = gene_values or [0.10, 0.20, 0.30, 0.05]
    other = [1.0 - v for v in vals]
    assert all(v >= 0 for v in other)
    lines = [
        ",".join(header),
        "SREBF2," + ",".join(str(v) for v in vals),
        "OTHER," + ",".join(str(v) for v in other),
    ]
    return "\n".join(lines) + "\n"


def _full_members(**overrides: bytes) -> dict[str, bytes]:
    members = {
        "Contents.txt": b"Allen Human Brain Atlas\n",
        "Genes.csv": _genes_csv().encode(),
        "SampleAnnot.csv": _sample_annot().encode(),
        "RNAseqTPM.csv": _tpm_csv().encode(),
        "RNAseqCounts.csv": b"gene_symbol,S_BIO_0\nSREBF2,1\n",
        "Ontology.csv": b"id,name\n1,brain\n",
    }
    members.update(overrides)
    return members


def _mean(vals: list[float]) -> float:
    return sum(vals) / len(vals)


def test_download_url_for_donor():
    assert ahr.download_url_for_donor(278447594).endswith("/278447594")


def test_unpack_rejects_missing_and_duplicate_members():
    incomplete = _zip_bytes({"Contents.txt": b"x", "Genes.csv": b"a"})
    bad = ahr.unpack_zip_members(incomplete)
    assert bad["success"] is False
    assert bad["error_type"] == "missing_zip_member"

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, content in _full_members().items():
            zf.writestr(f"a/{name}", content)
            if name == "Genes.csv":
                zf.writestr(f"b/{name}", content)
    dup = ahr.unpack_zip_members(buf.getvalue())
    assert dup["success"] is False
    assert dup["error_type"] == "duplicate_zip_member"


def test_unpack_rejects_html_and_bad_magic():
    assert ahr.unpack_zip_members(b"<html>not zip</html>")["error_type"] == (
        "invalid_zip_magic"
    )
    assert ahr.unpack_zip_members(b"notazip")["success"] is False


def test_detect_tpm_scale_fraction_and_conventional():
    frac = ahr.detect_tpm_scale([0.98, 1.01, 1.00])
    assert frac["success"] is True
    assert frac["scale_mode"] == "fraction_of_million"
    assert frac["scale_multiplier"] == 1_000_000.0
    assert "tolerances" in frac

    conv = ahr.detect_tpm_scale([1_000_000.0, 990_000.0, 1_010_000.0])
    assert conv["success"] is True
    assert conv["scale_mode"] == "conventional_tpm"

    bad = ahr.detect_tpm_scale([50.0, 60.0, 70.0])
    assert bad["success"] is False


def test_filter_biological_rejects_unknown_replicate_and_dup_well():
    rows = [
        {
            "RNAseq_sample_name": "A",
            "replicate_sample": "Maybe",
            "well_id": "1",
            "brain": "1",
        }
    ]
    out = ahr.filter_biological_samples(rows)
    assert out["success"] is False
    assert out["error_type"] == "unrecognized_replicate_sample"

    rows = [
        {
            "RNAseq_sample_name": "A",
            "replicate_sample": "No",
            "well_id": "1",
            "brain": "1",
        },
        {
            "RNAseq_sample_name": "B",
            "replicate_sample": "No",
            "well_id": "1",
            "brain": "1",
        },
    ]
    out = ahr.filter_biological_samples(rows)
    assert out["success"] is False
    assert out["error_type"] == "duplicate_well_id"


def test_parse_donor_package_and_pooled_mean():
    members = _full_members()
    zip_bytes = _zip_bytes(members)
    parsed = ahr.parse_donor_package(
        zip_bytes,
        gene_symbol="SREBF2",
        entrez_gene_id="6721",
        well_known_file_id=278447594,
    )
    assert parsed["success"] is True
    assert parsed["brain_identity"] == "1"
    assert parsed["retained_sample_count"] == 3
    assert parsed["technical_replicate_count"] == 1
    assert parsed["scale_detection"]["scale_mode"] == "fraction_of_million"
    assert abs(parsed["mean_tpm"] - _mean([0.10e6, 0.20e6, 0.30e6])) < 1e-6
    assert set(parsed["derived_members"]) == set(ahr.PERSISTED_DERIVED_MEMBERS)
    assert parsed["validated_only_members_present"]["RNAseqCounts.csv"] is True

    members2 = _full_members()
    members2["SampleAnnot.csv"] = _sample_annot(brain="2").encode()
    parsed2 = ahr.parse_donor_package(
        _zip_bytes(members2),
        gene_symbol="SREBF2",
        entrez_gene_id="6721",
        well_known_file_id=278448166,
    )
    assert parsed2["success"] is True
    pooled = ahr.pooled_mean_across_donors([parsed, parsed2])
    assert pooled["success"] is True
    assert pooled["retained_sample_count"] == 6


def test_pooled_requires_distinct_brains():
    members = _full_members()
    p1 = ahr.parse_donor_package(
        _zip_bytes(members),
        gene_symbol="SREBF2",
        entrez_gene_id="6721",
        well_known_file_id=1,
    )
    p2 = ahr.parse_donor_package(
        _zip_bytes(members),
        gene_symbol="SREBF2",
        entrez_gene_id="6721",
        well_known_file_id=2,
    )
    pooled = ahr.pooled_mean_across_donors([p1, p2])
    assert pooled["success"] is False
    assert pooled["error_type"] == "brain_identity_not_distinct"


def test_download_donor_zip_never_raises(monkeypatch):
    class _Resp:
        status_code = 200
        content = b"<html>nope</html>"
        url = "https://example.test/x"
        headers = {"content-type": "text/html"}

        @property
        def is_success(self):
            return True

    class _Client:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, *a, **k):
            return _Resp()

    monkeypatch.setattr(ahr.httpx, "Client", _Client)
    result = ahr.download_donor_zip(278447594, gene_symbol="SREBF2")
    assert result.success is False
    assert result.error_type == "html_masquerading_as_zip"


def test_gene_entrez_mismatch_fails():
    members = _full_members()
    parsed = ahr.parse_donor_package(
        _zip_bytes(members),
        gene_symbol="SREBF2",
        entrez_gene_id="999999",
        well_known_file_id=1,
    )
    assert parsed["success"] is False
    assert parsed["error_type"] == "gene_identity_mismatch"
