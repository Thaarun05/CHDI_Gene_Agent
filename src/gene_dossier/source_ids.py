"""Deterministic ``source_id`` generation.

A ``source_id`` is the stable, human-readable key that ties a report claim back to a
specific normalized :class:`~gene_dossier.models.EvidenceRecord`. It MUST be
deterministic: the same source + gene + assertion + record identity always produces the
same id, so re-running a dossier is idempotent and claims can reliably cite evidence.

Format::

    {source}:{gene}:{assertion}:{key}

Each component is slugified (lowercased, non-alphanumerics collapsed to ``-``). When the
key material is long or contains characters that would be lost by slugification, a short
deterministic hash suffix is appended so ids stay unique and stable.

Examples::

    make_source_id("NCBI Gene", "SREBF2", "gene_identity", "6721")
        -> "ncbi-gene:srebf2:gene_identity:6721"
    make_source_id("PubMed", "SREBF2", "literature_summary", "PMID:12345678")
        -> "pubmed:srebf2:literature_summary:pmid-12345678"
"""

from __future__ import annotations

import hashlib
import re
from enum import Enum

# Component separator and the maximum length a key slug may reach before we switch to a
# hashed representation (keeps ids readable and filesystem/URL friendly).
SEP = ":"
_MAX_KEY_SLUG_LEN = 60
_HASH_LEN = 10

# Number of ``SEP``-delimited components in a well-formed source id.
_NUM_COMPONENTS = 4

_slug_re = re.compile(r"[^a-z0-9]+")


def _to_str(value: object) -> str:
    """Return the string form of ``value``, using ``.value`` for enums."""
    if isinstance(value, Enum):
        return str(value.value)
    return str(value)


def slugify(text: object, *, allow_underscore: bool = False) -> str:
    """Lowercase ``text`` and collapse runs of non-alphanumerics into single ``-``.

    Args:
        text: Value to slugify (enums use their ``.value``).
        allow_underscore: If True, underscores are preserved (used for source names and
            assertion types where ``_`` is meaningful, e.g. ``gene_identity``).

    Returns:
        A slug containing ``[a-z0-9]`` plus ``-`` (and ``_`` when allowed). Empty input
        yields ``""``.
    """
    raw = _to_str(text).strip().lower()
    if not raw:
        return ""
    if allow_underscore:
        # Treat underscore as a keepable character, everything else non-alnum -> "-".
        cleaned = re.sub(r"[^a-z0-9_]+", "-", raw)
    else:
        cleaned = _slug_re.sub("-", raw)
    return cleaned.strip("-_") if allow_underscore else cleaned.strip("-")


def _short_hash(text: str) -> str:
    """Return a short, stable hex digest of ``text``."""
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:_HASH_LEN]


def _key_slug(key_parts: tuple[object, ...]) -> str:
    """Build the key component from one or more parts.

    Parts are slugified and joined with ``-``. If the result is empty or exceeds
    ``_MAX_KEY_SLUG_LEN``, a deterministic hash of the raw joined parts is used instead
    (prefixed so it is recognizable).
    """
    raw_joined = "-".join(_to_str(p) for p in key_parts if _to_str(p).strip())
    slug = "-".join(s for s in (slugify(p) for p in key_parts) if s)
    if not slug:
        # Nothing survived slugification (e.g. all symbols) but there was raw content.
        return f"h-{_short_hash(raw_joined)}" if raw_joined else ""
    if len(slug) > _MAX_KEY_SLUG_LEN:
        return f"h-{_short_hash(raw_joined)}"
    return slug


def make_source_id(
    source_name: str,
    gene_symbol: str,
    assertion: object,
    *key_parts: object,
) -> str:
    """Build a deterministic ``source_id``.

    Args:
        source_name: Human name of the source (e.g. ``"NCBI Gene"``).
        gene_symbol: Gene symbol (e.g. ``"SREBF2"``).
        assertion: Assertion type (``AssertionType`` enum or string).
        *key_parts: One or more values uniquely identifying the record within the
            source (e.g. an accession, PMID, pathway id, interactor symbol).

    Returns:
        A stable id of the form ``{source}:{gene}:{assertion}:{key}``.

    Raises:
        ValueError: If ``source_name``, ``gene_symbol``, ``assertion``, or the resulting
            key slug is empty.
    """
    source = slugify(source_name, allow_underscore=True)
    gene = slugify(gene_symbol)
    assertion_slug = slugify(assertion, allow_underscore=True)
    key = _key_slug(key_parts)

    missing = [
        name
        for name, val in (
            ("source_name", source),
            ("gene_symbol", gene),
            ("assertion", assertion_slug),
            ("key", key),
        )
        if not val
    ]
    if missing:
        raise ValueError(f"cannot build source_id; empty component(s): {', '.join(missing)}")

    return SEP.join((source, gene, assertion_slug, key))


def is_valid_source_id(source_id: str) -> bool:
    """Return True if ``source_id`` has the expected 4-component shape with no blanks."""
    if not isinstance(source_id, str):
        return False
    parts = source_id.split(SEP)
    return len(parts) == _NUM_COMPONENTS and all(p.strip() for p in parts)


def parse_source_id(source_id: str) -> dict[str, str]:
    """Split a ``source_id`` into its named components.

    Args:
        source_id: An id produced by :func:`make_source_id`.

    Returns:
        Dict with keys ``source``, ``gene``, ``assertion``, ``key``.

    Raises:
        ValueError: If ``source_id`` is not a valid 4-component id.
    """
    if not is_valid_source_id(source_id):
        raise ValueError(f"invalid source_id: {source_id!r}")
    source, gene, assertion, key = source_id.split(SEP)
    return {"source": source, "gene": gene, "assertion": assertion, "key": key}


__all__ = [
    "SEP",
    "slugify",
    "make_source_id",
    "is_valid_source_id",
    "parse_source_id",
]
