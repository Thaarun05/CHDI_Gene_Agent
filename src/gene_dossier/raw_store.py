"""Raw artifact storage.

Persists the *original source material* (API responses, TSV dumps, screenshots, manual
notes) to disk and returns a :class:`~gene_dossier.models.RawArtifact` describing it.
This is the bottom of the provenance chain: every :class:`EvidenceRecord` should point
back to a stored raw artifact via ``raw_artifact_id``.

Layout on disk::

    {base_dir}/{dossier_run_id}/{source_slug}/{hint-}{hash12}.{ext}

Files are named by their content hash, so writing identical content is idempotent
(re-running a source overwrites the byte-identical file instead of duplicating it).
JSON is serialized with sorted keys so equal data always yields the same hash.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from .config import get_settings
from .models import RawArtifact
from .source_ids import slugify

# Number of hex chars from the content hash used in filenames.
_FILENAME_HASH_LEN = 12


def compute_hash(content: bytes) -> str:
    """Return the SHA-256 hex digest of ``content``."""
    return hashlib.sha256(content).hexdigest()


class RawStore:
    """Content-addressed store for raw source artifacts."""

    def __init__(self, base_dir: str | Path | None = None) -> None:
        """Create a store rooted at ``base_dir`` (defaults to ``settings.raw_data_path``)."""
        self.base_dir = Path(base_dir) if base_dir is not None else get_settings().raw_data_path

    # -- internal helpers ---------------------------------------------------------------
    def _dir_for(self, dossier_run_id: str, source_name: str) -> Path:
        """Return (creating if needed) the directory for a run + source."""
        run_part = slugify(dossier_run_id) or "misc"
        source_part = slugify(source_name, allow_underscore=True) or "unknown-source"
        target = self.base_dir / run_part / source_part
        target.mkdir(parents=True, exist_ok=True)
        return target

    def _write(
        self,
        *,
        dossier_run_id: str,
        source_name: str,
        content: bytes,
        artifact_type: str,
        extension: str,
        api_run_id: str | None,
        original_url: str | None,
        notes: str | None,
        filename_hint: str | None,
    ) -> RawArtifact:
        """Write ``content`` to disk and return the describing :class:`RawArtifact`."""
        content_hash = compute_hash(content)
        directory = self._dir_for(dossier_run_id, source_name)
        hint = slugify(filename_hint) if filename_hint else ""
        stem = f"{hint}-" if hint else ""
        path = directory / f"{stem}{content_hash[:_FILENAME_HASH_LEN]}.{extension.lstrip('.')}"
        path.write_bytes(content)
        return RawArtifact(
            dossier_run_id=dossier_run_id,
            api_run_id=api_run_id,
            source_name=source_name,
            artifact_type=artifact_type,
            file_path=str(path),
            original_url=original_url,
            content_hash=content_hash,
            notes=notes,
        )

    # -- public save API ----------------------------------------------------------------
    def save_json(
        self,
        dossier_run_id: str,
        source_name: str,
        data: Any,
        *,
        api_run_id: str | None = None,
        original_url: str | None = None,
        notes: str | None = None,
        filename_hint: str | None = None,
    ) -> RawArtifact:
        """Serialize ``data`` to canonical JSON and store it.

        Keys are sorted so equal data hashes identically (idempotent storage).
        """
        content = json.dumps(
            data, sort_keys=True, ensure_ascii=False, indent=2
        ).encode("utf-8")
        return self._write(
            dossier_run_id=dossier_run_id,
            source_name=source_name,
            content=content,
            artifact_type="json",
            extension="json",
            api_run_id=api_run_id,
            original_url=original_url,
            notes=notes,
            filename_hint=filename_hint,
        )

    def save_text(
        self,
        dossier_run_id: str,
        source_name: str,
        text: str,
        *,
        artifact_type: str = "text",
        extension: str = "txt",
        api_run_id: str | None = None,
        original_url: str | None = None,
        notes: str | None = None,
        filename_hint: str | None = None,
    ) -> RawArtifact:
        """Store a text artifact (e.g. TSV: pass ``extension="tsv"``)."""
        return self._write(
            dossier_run_id=dossier_run_id,
            source_name=source_name,
            content=text.encode("utf-8"),
            artifact_type=artifact_type,
            extension=extension,
            api_run_id=api_run_id,
            original_url=original_url,
            notes=notes,
            filename_hint=filename_hint,
        )

    def save_bytes(
        self,
        dossier_run_id: str,
        source_name: str,
        content: bytes,
        *,
        extension: str,
        artifact_type: str = "binary",
        api_run_id: str | None = None,
        original_url: str | None = None,
        notes: str | None = None,
        filename_hint: str | None = None,
    ) -> RawArtifact:
        """Store an arbitrary binary artifact (e.g. a screenshot or PDF)."""
        return self._write(
            dossier_run_id=dossier_run_id,
            source_name=source_name,
            content=content,
            artifact_type=artifact_type,
            extension=extension,
            api_run_id=api_run_id,
            original_url=original_url,
            notes=notes,
            filename_hint=filename_hint,
        )

    # -- public load API ----------------------------------------------------------------
    def load_bytes(self, artifact: RawArtifact) -> bytes:
        """Read the raw bytes for ``artifact`` back from disk."""
        return Path(artifact.file_path).read_bytes()

    def load_text(self, artifact: RawArtifact) -> str:
        """Read a text artifact back as a string."""
        return Path(artifact.file_path).read_text(encoding="utf-8")

    def load_json(self, artifact: RawArtifact) -> Any:
        """Read a JSON artifact back as a Python object."""
        return json.loads(self.load_text(artifact))

    def verify(self, artifact: RawArtifact) -> bool:
        """Return True if the on-disk content still matches ``artifact.content_hash``."""
        path = Path(artifact.file_path)
        if not path.exists():
            return False
        return compute_hash(path.read_bytes()) == artifact.content_hash


__all__ = ["RawStore", "compute_hash"]
