"""Protected dataset-level and gene-level layout for Section 6a (CTD).

Layout::

    data/outputs/section_6a/
      attempts/sources/{timestamp}_ctd_chem_gene_ixns/
      accepted/sources/ctd_chem_gene_ixns.json
      attempts/genes/{timestamp}_{GENE}/   # or {run_id}/{GENE}/
      accepted/genes/{GENE}.json

The shared bulk gzip is pinned once. Gene attempts reuse the accepted pointer
and must cite the original ApiRun / RawArtifact / source_attempt_id.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from gene_dossier.tools.ctd import CHEM_GENE_IXNS_BULK_URL, CHEM_GENE_IXNS_SOURCE_KEY

SECTION_6A_OUTPUT_DIRNAME = "section_6a"
ATTEMPTS_DIRNAME = "attempts"
ACCEPTED_DIRNAME = "accepted"
SOURCES_DIRNAME = "sources"
GENES_DIRNAME = "genes"
MANIFEST_FILENAME = "manifest.json"
SOURCE_KEY = CHEM_GENE_IXNS_SOURCE_KEY
OFFICIAL_URL = CHEM_GENE_IXNS_BULK_URL


def utc_stamp(moment: datetime | None = None) -> str:
    now = moment or datetime.now(timezone.utc)
    return now.strftime("%Y%m%dT%H%M%SZ")


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", (value or "").strip())
    return token.strip("._-") or "unknown"


@dataclass(frozen=True)
class Section6aPaths:
    root: Path

    @property
    def source_attempts(self) -> Path:
        return self.root / ATTEMPTS_DIRNAME / SOURCES_DIRNAME

    @property
    def gene_attempts(self) -> Path:
        return self.root / ATTEMPTS_DIRNAME / GENES_DIRNAME

    @property
    def accepted_sources(self) -> Path:
        return self.root / ACCEPTED_DIRNAME / SOURCES_DIRNAME

    @property
    def accepted_genes(self) -> Path:
        return self.root / ACCEPTED_DIRNAME / GENES_DIRNAME

    def ensure(self) -> None:
        for path in (
            self.source_attempts,
            self.gene_attempts,
            self.accepted_sources,
            self.accepted_genes,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def new_source_attempt(self, source_key: str = SOURCE_KEY, *, moment: datetime | None = None) -> Path:
        self.ensure()
        stamp = utc_stamp(moment)
        path = self.source_attempts / f"{stamp}_{_safe_token(source_key)}"
        if path.exists():
            suffix = 2
            while True:
                candidate = self.source_attempts / f"{stamp}_{_safe_token(source_key)}_{suffix}"
                if not candidate.exists():
                    path = candidate
                    break
                suffix += 1
        path.mkdir(parents=True, exist_ok=False)
        return path

    def new_gene_attempt(
        self,
        gene_symbol: str,
        *,
        run_id: str | None = None,
        moment: datetime | None = None,
    ) -> Path:
        self.ensure()
        run_token = _safe_token(run_id or utc_stamp(moment))
        gene_token = _safe_token(gene_symbol.upper())
        base = self.gene_attempts / run_token / gene_token
        if not base.exists():
            base.mkdir(parents=True, exist_ok=False)
            return base
        suffix = 2
        while True:
            candidate = self.gene_attempts / f"{run_token}_{suffix}" / gene_token
            if not candidate.exists():
                candidate.mkdir(parents=True, exist_ok=False)
                return candidate
            suffix += 1

    def accepted_source_pointer(self, source_key: str = SOURCE_KEY) -> Path:
        return self.accepted_sources / f"{_safe_token(source_key)}.json"

    def accepted_gene_pointer(self, gene_symbol: str) -> Path:
        return self.accepted_genes / f"{_safe_token(gene_symbol.upper())}.json"


def paths_for(output_root: str | Path) -> Section6aPaths:
    root = Path(output_root)
    if root.name != SECTION_6A_OUTPUT_DIRNAME:
        root = root / SECTION_6A_OUTPUT_DIRNAME
    paths = Section6aPaths(root=root)
    paths.ensure()
    return paths


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(
                payload,
                fh,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                default=str,
            )
            fh.write("\n")
        os.replace(tmp_name, path)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise


def accept_source(
    paths: Section6aPaths,
    *,
    source_key: str,
    attempt_dir: Path,
    artifact_path: Path,
    official_url: str,
    sha256: str,
    byte_size: int,
    validation: dict[str, Any],
    extra: dict[str, Any] | None = None,
) -> Path:
    """Atomically pin a validated CTD bulk attempt.

    ``extra`` must carry resolvable ORIGINAL download provenance:
    api_run_id, raw_artifact_id, source_attempt_id, ctd_report_created,
    retrieval_timestamp, parser_version, dossier_run_id.

    Refuses to pin when ``api_run_id`` or ``raw_artifact_id`` is missing.
    """
    extra = dict(extra or {})
    api_run_id = str(extra.get("api_run_id") or "").strip()
    raw_artifact_id = str(extra.get("raw_artifact_id") or "").strip()
    if not api_run_id or not raw_artifact_id:
        raise ValueError(
            "refusing to accept CTD bulk source without real api_run_id and "
            "raw_artifact_id"
        )
    pointer = paths.accepted_source_pointer(source_key)
    payload: dict[str, Any] = {
        "source_key": source_key,
        "official_url": official_url,
        "artifact_path": str(artifact_path),
        "attempt_dir": str(attempt_dir),
        "source_attempt_id": str(extra.get("source_attempt_id") or Path(attempt_dir).name),
        "sha256": sha256,
        "byte_size": byte_size,
        "validation_status": "accepted",
        "validation": validation,
        "accepted_at": datetime.now(timezone.utc).isoformat(),
        "api_run_id": api_run_id,
        "raw_artifact_id": raw_artifact_id,
    }
    payload.update(extra)
    payload["api_run_id"] = api_run_id
    payload["raw_artifact_id"] = raw_artifact_id
    write_json_atomic(pointer, payload)
    return pointer


def load_accepted_source(
    paths: Section6aPaths,
    *,
    source_key: str = SOURCE_KEY,
    official_url: str | None = OFFICIAL_URL,
) -> dict[str, Any] | None:
    """Return accepted CTD bulk record only on a genuine cache hit."""
    pointer = paths.accepted_source_pointer(source_key)
    if not pointer.is_file():
        return None
    try:
        record = json.loads(pointer.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(record, dict):
        return None
    if record.get("source_key") != source_key:
        return None
    if record.get("validation_status") != "accepted":
        return None
    if official_url and record.get("official_url") != official_url:
        return None
    artifact = Path(str(record.get("artifact_path") or ""))
    if not artifact.is_file():
        return None
    stored = str(record.get("sha256") or "")
    if not stored or sha256_file(artifact) != stored:
        return None
    # Genuine cache hit requires resolvable ORIGINAL download provenance.
    if not str(record.get("api_run_id") or "").strip():
        return None
    if not str(record.get("raw_artifact_id") or "").strip():
        return None
    return record


def accept_gene_report(
    paths: Section6aPaths,
    *,
    gene_symbol: str,
    attempt_dir: Path,
    acceptance: dict[str, Any],
    artifacts: dict[str, Any] | None = None,
) -> Path:
    pointer = paths.accepted_gene_pointer(gene_symbol)
    payload = {
        "gene_symbol": gene_symbol,
        "attempt_dir": str(attempt_dir),
        "acceptance": acceptance,
        "artifacts": artifacts or {},
        "accepted_at": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(pointer, payload)
    return pointer


__all__ = [
    "MANIFEST_FILENAME",
    "OFFICIAL_URL",
    "SECTION_6A_OUTPUT_DIRNAME",
    "SOURCE_KEY",
    "Section6aPaths",
    "accept_gene_report",
    "accept_source",
    "load_accepted_source",
    "paths_for",
    "sha256_bytes",
    "sha256_file",
    "utc_stamp",
    "write_json_atomic",
]
