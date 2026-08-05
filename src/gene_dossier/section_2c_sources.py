"""Protected acquisition of Section 2c dataset-level sources.

Human M1 and the DropViz GEO matrix are *dataset-level*: they are shared by every
gene, so they are acquired once, validated, and pinned by an accepted pointer.
No gene run may claim to own them, and a failed or partial retry must never
replace a previously accepted source.

Layout::

    data/outputs/section_2c/
      attempts/sources/{timestamp}_{source_key}/   # immutable, one per attempt
      accepted/sources/{source_key}.json           # atomic pointer to the accepted attempt

Acceptance requires a validated payload. A cache hit requires the accepted
pointer, a matching source key and official URL, bytes that still exist on disk,
and a matching SHA-256 - a pointer to a missing file is not a hit.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SECTION_2C_OUTPUT_DIRNAME = "section_2c"
ATTEMPTS_DIRNAME = "attempts"
ACCEPTED_DIRNAME = "accepted"
SOURCES_DIRNAME = "sources"
GENES_DIRNAME = "genes"

MANIFEST_FILENAME = "manifest.json"


def utc_stamp(moment: datetime | None = None) -> str:
    """Compact UTC timestamp used for attempt directory names."""
    now = moment or datetime.now(timezone.utc)
    return now.strftime("%Y%m%dT%H%M%SZ")


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


@dataclass(frozen=True)
class Section2cPaths:
    """Resolved directory layout for protected Section 2c outputs."""

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

    def new_source_attempt(self, source_key: str, *, moment: datetime | None = None) -> Path:
        self.ensure()
        stamp = utc_stamp(moment)
        path = self.source_attempts / f"{stamp}_{source_key}"
        if path.exists():
            suffix = 2
            while True:
                candidate = self.source_attempts / f"{stamp}_{source_key}_{suffix}"
                if not candidate.exists():
                    path = candidate
                    break
                suffix += 1
        path.mkdir(parents=True, exist_ok=False)
        return path

    def new_gene_attempt(self, gene_symbol: str, *, moment: datetime | None = None) -> Path:
        self.ensure()
        stamp = utc_stamp(moment)
        path = self.gene_attempts / f"{stamp}_{gene_symbol}"
        if path.exists():
            # Same-second retries (offline tests, rapid CLI re-runs) must not
            # collide with an earlier attempt directory.
            suffix = 2
            while True:
                candidate = self.gene_attempts / f"{stamp}_{gene_symbol}_{suffix}"
                if not candidate.exists():
                    path = candidate
                    break
                suffix += 1
        path.mkdir(parents=True, exist_ok=False)
        return path

    def accepted_source_pointer(self, source_key: str) -> Path:
        return self.accepted_sources / f"{source_key}.json"

    def accepted_gene_pointer(self, gene_symbol: str) -> Path:
        return self.accepted_genes / f"{gene_symbol}.json"


def paths_for(output_root: str | Path) -> Section2cPaths:
    """Return the protected layout rooted at ``output_root/section_2c``."""
    root = Path(output_root)
    if root.name != SECTION_2C_OUTPUT_DIRNAME:
        root = root / SECTION_2C_OUTPUT_DIRNAME
    return Section2cPaths(root=root)


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON via a temp file + atomic replace so readers never see a partial file."""
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
    paths: Section2cPaths,
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
    """Atomically pin ``attempt_dir`` as the accepted artifact for ``source_key``."""
    pointer = paths.accepted_source_pointer(source_key)
    payload: dict[str, Any] = {
        "source_key": source_key,
        "official_url": official_url,
        "artifact_path": str(artifact_path),
        "attempt_dir": str(attempt_dir),
        "sha256": sha256,
        "byte_size": byte_size,
        "validation_status": "accepted",
        "validation": validation,
        "accepted_at": datetime.now(timezone.utc).isoformat(),
    }
    if extra:
        payload.update(extra)
    write_json_atomic(pointer, payload)
    return pointer


def load_accepted_source(
    paths: Section2cPaths,
    *,
    source_key: str,
    official_url: str | None = None,
) -> dict[str, Any] | None:
    """Return the accepted source record only when it is a genuine cache hit.

    All of these must hold: the pointer exists and is accepted, the source key
    matches, the official URL matches when one is supplied, the artifact bytes
    still exist, and the on-disk SHA-256 still matches the recorded digest.
    """
    pointer = paths.accepted_source_pointer(source_key)
    if not pointer.is_file():
        return None
    try:
        record = json.loads(pointer.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(record, dict):
        return None
    if record.get("validation_status") != "accepted":
        return None
    if record.get("source_key") != source_key:
        return None
    if official_url is not None and record.get("official_url") != official_url:
        return None

    artifact = Path(str(record.get("artifact_path") or ""))
    if not artifact.is_file():
        return None
    expected = str(record.get("sha256") or "")
    if not expected:
        return None
    if sha256_bytes(artifact.read_bytes()) != expected:
        return None
    return record


def accept_gene_report(
    paths: Section2cPaths,
    *,
    gene_symbol: str,
    attempt_dir: Path,
    acceptance: dict[str, Any],
    artifacts: dict[str, Any] | None = None,
) -> Path:
    """Pin a gene attempt as the accepted report for ``gene_symbol``."""
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
    "Section2cPaths",
    "accept_gene_report",
    "accept_source",
    "load_accepted_source",
    "paths_for",
    "sha256_bytes",
    "utc_stamp",
    "write_json_atomic",
]
