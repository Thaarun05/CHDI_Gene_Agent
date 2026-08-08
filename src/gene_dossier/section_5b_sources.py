"""Protected attempt layout and accepted gene pointers for Section 5b.

Layout::

    data/outputs/section_5b/
      attempts/{RUN_ID}/{GENE}/
      accepted/genes/{GENE}.json
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

SECTION_5B_OUTPUT_DIRNAME = "section_5b"
ATTEMPTS_DIRNAME = "attempts"
ACCEPTED_DIRNAME = "accepted"
GENES_DIRNAME = "genes"
MANIFEST_FILENAME = "manifest.json"


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
class Section5bPaths:
    root: Path

    @property
    def attempts(self) -> Path:
        return self.root / ATTEMPTS_DIRNAME

    @property
    def accepted_genes(self) -> Path:
        return self.root / ACCEPTED_DIRNAME / GENES_DIRNAME

    def ensure(self) -> None:
        self.attempts.mkdir(parents=True, exist_ok=True)
        self.accepted_genes.mkdir(parents=True, exist_ok=True)

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
        base = self.attempts / run_token / gene_token
        if not base.exists():
            base.mkdir(parents=True, exist_ok=False)
            return base
        suffix = 2
        while True:
            candidate = self.attempts / f"{run_token}_{suffix}" / gene_token
            if not candidate.exists():
                candidate.mkdir(parents=True, exist_ok=False)
                return candidate
            suffix += 1

    def accepted_gene_pointer(self, gene_symbol: str) -> Path:
        return self.accepted_genes / f"{_safe_token(gene_symbol.upper())}.json"


def paths_for(output_root: str | Path) -> Section5bPaths:
    root = Path(output_root)
    if root.name != SECTION_5B_OUTPUT_DIRNAME:
        root = root / SECTION_5B_OUTPUT_DIRNAME
    paths = Section5bPaths(root=root)
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


def accept_gene_report(
    paths: Section5bPaths,
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
    "SECTION_5B_OUTPUT_DIRNAME",
    "Section5bPaths",
    "accept_gene_report",
    "paths_for",
    "sha256_bytes",
    "sha256_file",
    "utc_stamp",
    "write_json_atomic",
]
