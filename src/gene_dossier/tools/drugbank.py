"""DrugBank stub client (license / credentials constrained).

No scraping. Until DrugBank credentials exist on :class:`~gene_dossier.config.Settings`,
all status fetches return ``unavailable_not_configured``.

Never raises: failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "DrugBank"

# Candidate Settings field / env names. None exist on Settings today.
_DRUGBANK_KEY_CANDIDATES = (
    "DRUGBANK_API_KEY",
    "DRUGBANK_ACCESS_KEY",
    "DRUGBANK_USERNAME",
    "DRUGBANK_PASSWORD",
    "DRUGBANK",
)


def _tool_result(
    *,
    endpoint_name: str,
    gene_symbol: str,
    request_url: str,
    request_params: dict[str, Any],
    success: bool,
    status_code: int | None = None,
    data: Any | None = None,
    error_type: str | None = None,
    error_message: str | None = None,
) -> ToolResult:
    return ToolResult(
        source_name=SOURCE_NAME,
        endpoint_name=endpoint_name,
        success=success,
        gene_symbol=gene_symbol,
        request_url=request_url,
        request_params=request_params,
        status_code=status_code,
        data=data,
        error_type=error_type,
        error_message=error_message,
    )


def check_configured(settings: Settings | None = None) -> bool:
    """Return True only when a DrugBank credential field exists and is non-empty.

    Uses :meth:`Settings.has_key` for known ``DRUGBANK*`` names when the
    corresponding Settings attribute exists; otherwise False.
    """
    cfg = settings or get_settings()
    for name in _DRUGBANK_KEY_CANDIDATES:
        field = name.lower()
        if not hasattr(cfg, field):
            continue
        if cfg.has_key(name):
            return True
    return False


def fetch_status(
    gene_symbol: str = "",
    *,
    settings: Settings | None = None,
) -> ToolResult:
    """Report DrugBank availability. Never scrapes DrugBank HTML/APIs.

    When not configured, always returns ``error_type="unavailable_not_configured"``.
    """
    cfg = settings or get_settings()
    if not check_configured(cfg):
        return _tool_result(
            endpoint_name="fetch_status",
            gene_symbol=gene_symbol,
            request_url="",
            request_params={"gene_symbol": gene_symbol},
            success=False,
            data={
                "status": "unavailable_not_configured",
                "message": "DrugBank API access unavailable",
            },
            error_type="unavailable_not_configured",
            error_message="DrugBank API access unavailable (not configured)",
        )
    # Credentials exist but no licensed client is wired — still no scrape.
    return _tool_result(
        endpoint_name="fetch_status",
        gene_symbol=gene_symbol,
        request_url="",
        request_params={"gene_symbol": gene_symbol},
        success=False,
        data={
            "status": "unavailable_not_implemented",
            "message": "DrugBank credentials present but client not implemented",
        },
        error_type="unavailable_not_implemented",
        error_message="DrugBank client not implemented (no scrape)",
    )


__all__ = [
    "SOURCE_NAME",
    "check_configured",
    "fetch_status",
]
