"""Biomedical API clients for the Gene Dossier Platform.

Each module in this package wraps one external source. Clients:

- build request URLs and parameters
- call APIs with timeouts
- never raise unhandled exceptions
- return :class:`~gene_dossier.models.ToolResult`
- do **not** normalize into evidence records (that happens in ``normalize/``)

Clients are built one file at a time (Priority A → B → C).
"""

from __future__ import annotations

__all__: list[str] = []
