"""Gene Dossier Platform.

A provenance-first platform for generating CHDI-style gene dossiers for Huntington's
disease research. Every reported fact traces back to a ``source_id``, and every
``source_id`` traces back to a raw API response or artifact. The LLM is never the
source of truth.

See ``IMPLEMENTATION_PLAN.md`` for the ordered build plan and architecture.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
