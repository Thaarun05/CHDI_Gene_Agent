"""Normalizers: turn successful :class:`~gene_dossier.models.ToolResult` payloads
into :class:`~gene_dossier.models.EvidenceRecord` lists.

Rules:
- No network I/O
- Do not invent facts beyond what the client payload contains
- Prefer curated, safely selected identifiers (skip ambiguous selections)
"""

from __future__ import annotations

from gene_dossier.normalize.gene_identity import normalize_gene_identity

__all__ = ["normalize_gene_identity"]
