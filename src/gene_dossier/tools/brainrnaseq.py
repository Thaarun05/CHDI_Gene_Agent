"""BrainRNASeq client (CSV bulk download + local gene filter).

Downloads published human/mouse cell-type expression CSVs and filters rows for a
gene. Priority C scaffold: CSV download, not a JSON API. Does **not** normalize
into evidence records — that belongs in ``normalize/expression.py``.

Key endpoints (validated)::

    GET https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-124.csv
        (human)
    GET https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-120.csv
        (mouse)

NOTE: Parse rows where ``gene_id`` / ``id`` matches the requested symbol
(case-insensitive). Prefer exact matches by default to avoid short-symbol
false positives. Downloads send polite ``User-Agent`` / ``Accept`` / ``Referer``
headers only; HTTP 403 is reported as ``access_forbidden`` (no access bypass).
When direct download is forbidden, :func:`download_csv_via_browser` may retrieve
the same published CSV bytes via a genuine Playwright browser context.

Never raises: all failures return :class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

import csv
import hashlib
import io
import logging
import tempfile
from pathlib import Path
from typing import Any, Literal
from urllib.parse import urlparse

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

logger = logging.getLogger(__name__)

SOURCE_NAME = "BrainRNASeq"
HUMAN_CSV_URL = (
    "https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-124.csv"
)
MOUSE_CSV_URL = (
    "https://brainrnaseq.org/wp-content/uploads/2022/09/fe-wp-dataset-120.csv"
)

# Polite identification only — not a browser spoof and not an access bypass.
REQUEST_HEADERS = {
    "User-Agent": "GeneDossier/0.1.0 (research; provenance-first gene dossier client)",
    "Accept": "text/csv,text/plain,*/*;q=0.8",
    "Referer": "https://brainrnaseq.org/",
}
SAFE_RESPONSE_HEADER_KEYS = ("content-type", "server", "cf-ray")
RAW_TEXT_PREVIEW_LIMIT = 500
MIN_CSV_BYTE_SIZE = 64
ID_COLUMN_CANDIDATES = ("gene_id", "id", "Gene", "gene", "symbol")

SpeciesChoice = Literal["human", "mouse"]

# Cell-type column prefixes noted in the validated API map.
# Section 2b also accepts a documented ``microglla_*`` typo alias when grouping.
HUMAN_CELLTYPE_PREFIXES = (
    "astrocytes_fetal_",
    "astrocytes_mature_",
    "endothelial_",
    "microglia_",
    "neurons_",
    "oligodendrocytes_",
)
MOUSE_CELLTYPE_PREFIXES = (
    "astrocytes_",
    "endothelial_",
    "microglia_macrophage_",
    "myelinating_oligodendrocyte_",
    "neurons_",
    "newly_formed_oligodendrocyte_",
    "opc_",
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
    """Build a uniform :class:`ToolResult` for this source."""
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


def csv_url_for_species(species: SpeciesChoice) -> str:
    """Return the validated BrainRNASeq CSV URL for ``species``."""
    if species == "mouse":
        return MOUSE_CSV_URL
    return HUMAN_CSV_URL


def published_csv_urls() -> frozenset[str]:
    """Exact published CSV URLs allowed for direct and browser retrieval."""
    return frozenset({HUMAN_CSV_URL, MOUSE_CSV_URL})


def normalize_published_csv_url(url: str) -> str | None:
    """Return the canonical published CSV URL, or None if not allowlisted.

    ``www.brainrnaseq.org`` is accepted only when the path matches a published
    CSV exactly (same path as the non-www canonical URL).
    """
    text = (url or "").strip()
    if not text:
        return None
    try:
        parsed = urlparse(text)
    except Exception:  # noqa: BLE001
        return None
    host = (parsed.hostname or "").lower()
    if host not in {"brainrnaseq.org", "www.brainrnaseq.org"}:
        return None
    path = parsed.path or ""
    for canonical in (HUMAN_CSV_URL, MOUSE_CSV_URL):
        canon = urlparse(canonical)
        if path == canon.path:
            return canonical
    return None


def is_html_payload(text: str) -> bool:
    """True when body looks like an HTML challenge / soft-200 page."""
    head = (text or "").lstrip("\ufeff").lstrip()[:200].lower()
    return head.startswith("<!doctype") or head.startswith("<html")


def celltype_prefixes_for_species(species: SpeciesChoice) -> tuple[str, ...]:
    """Return expected replicate column prefixes for ``species``."""
    if species == "mouse":
        return MOUSE_CELLTYPE_PREFIXES
    return HUMAN_CELLTYPE_PREFIXES


def validate_brainrnaseq_csv_bytes(
    content: bytes | str,
    *,
    species: SpeciesChoice,
) -> dict[str, Any]:
    """Validate BrainRNASeq CSV bytes used by direct and browser success paths.

    Requires nontrivial size, non-HTML body, parseable header with an id
    column, at least one expected species replicate prefix, and ≥1 data row.
    """
    if isinstance(content, bytes):
        raw_bytes = content
        try:
            text = content.decode("utf-8-sig")
        except UnicodeDecodeError as exc:
            return {
                "ok": False,
                "error_type": "invalid_or_html_response",
                "error_message": f"CSV bytes are not valid UTF-8: {exc}",
                "byte_size": len(content),
            }
    else:
        text = str(content or "")
        raw_bytes = text.encode("utf-8")

    byte_size = len(raw_bytes)
    if byte_size < MIN_CSV_BYTE_SIZE:
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": f"CSV payload too small ({byte_size} bytes)",
            "byte_size": byte_size,
        }
    if is_html_payload(text):
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": "Response body looks like HTML, not CSV",
            "byte_size": byte_size,
        }

    stripped = text.lstrip("\ufeff").strip()
    if not stripped:
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": "CSV body is empty",
            "byte_size": byte_size,
        }

    try:
        reader = csv.reader(io.StringIO(stripped))
        header = next(reader, None)
    except csv.Error as exc:
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": f"CSV header is not parseable: {exc}",
            "byte_size": byte_size,
        }
    if not header:
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": "CSV has no header row",
            "byte_size": byte_size,
        }

    header_norm = [str(h or "").strip() for h in header]
    header_lower = {h.lower() for h in header_norm if h}
    id_cols = [c for c in ID_COLUMN_CANDIDATES if c.lower() in header_lower]
    if not id_cols:
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": "CSV header lacks a gene identifier column",
            "byte_size": byte_size,
            "header": header_norm[:20],
        }

    prefixes = celltype_prefixes_for_species(species)
    prefix_hits = [
        col
        for col in header_norm
        if any(col.lower().startswith(prefix) for prefix in prefixes)
    ]
    if not prefix_hits:
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": f"CSV header lacks {species} cell-type replicate prefixes",
            "byte_size": byte_size,
            "header": header_norm[:20],
        }

    rows = parse_csv(stripped)
    if not rows:
        return {
            "ok": False,
            "error_type": "invalid_or_html_response",
            "error_message": "CSV has no data rows",
            "byte_size": byte_size,
        }

    return {
        "ok": True,
        "byte_size": byte_size,
        "sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "row_count": len(rows),
        "id_columns": id_cols,
        "replicate_column_count": len(prefix_hits),
        "content_type_hint": "text/csv",
    }


def parse_csv(raw_csv: str) -> list[dict[str, str]]:
    """Parse BrainRNASeq CSV text into row dicts."""
    text = (raw_csv or "").lstrip("\ufeff").strip()
    if not text:
        return []
    reader = csv.DictReader(io.StringIO(text))
    rows: list[dict[str, str]] = []
    for row in reader:
        if not isinstance(row, dict):
            continue
        cleaned = {
            str(k).strip(): ("" if v is None else str(v).strip())
            for k, v in row.items()
            if k is not None
        }
        if any(v for v in cleaned.values()):
            rows.append(cleaned)
    return rows


# Exact trailing species suffixes accepted in published BrainRNASeq gene_id values.
_SPECIES_ID_SUFFIXES: dict[SpeciesChoice, str] = {
    "human": " - Homo sapiens",
    "mouse": " - Mus musculus",
}
_ALL_SPECIES_ID_SUFFIXES: tuple[str, ...] = tuple(_SPECIES_ID_SUFFIXES.values())


def normalize_gene_identifier(
    value: str,
    *,
    species: SpeciesChoice | None = None,
) -> str | None:
    """Return the bare gene symbol for an accepted BrainRNASeq identifier.

    Accepts:
    - bare symbols (``GENEX``, ``Genex``)
    - exact trailing suffixes `` - Homo sapiens`` / `` - Mus musculus``

    Rejects arbitrary annotations such as ``GENEX - pseudogene`` or
    ``GENEX - antisense`` (any other `` - …`` suffix). When ``species`` is set,
    only that species' exact suffix is accepted in addition to bare symbols.
    """
    text = str(value or "").strip()
    if not text:
        return None
    lower = text.lower()
    allowed_suffixes = (
        (_SPECIES_ID_SUFFIXES[species],)
        if species is not None
        else _ALL_SPECIES_ID_SUFFIXES
    )
    for suffix in allowed_suffixes:
        if lower.endswith(suffix.lower()):
            return text[: -len(suffix)].strip() or None
    if " - " in text:
        return None
    return text


def row_matches_gene(
    row: dict[str, Any],
    gene_symbol: str,
    *,
    allow_substring_match: bool = False,
    species: SpeciesChoice | None = None,
) -> bool:
    """True if ``gene_id`` / ``id`` matches ``gene_symbol``.

    Exact case-insensitive match by default, including published
    ``SYMBOL - Homo sapiens`` / ``SYMBOL - Mus musculus`` forms only.
    Optional substring match is off by default to avoid short-symbol
    false positives.
    """
    target = gene_symbol.strip().lower()
    if not target:
        return False
    for key in ID_COLUMN_CANDIDATES:
        raw = str(row.get(key) or "").strip()
        if not raw:
            continue
        value = raw.lower()
        token = normalize_gene_identifier(raw, species=species)
        if token is None:
            continue
        token_l = token.lower()
        if value == target or token_l == target:
            return True
        if allow_substring_match and (target in value or target in token_l):
            return True
    return False


def filter_gene_rows(
    rows: list[dict[str, Any]],
    gene_symbol: str,
    *,
    allow_substring_match: bool = False,
    species: SpeciesChoice | None = None,
) -> list[dict[str, Any]]:
    """Filter CSV rows for ``gene_symbol``."""
    return [
        row
        for row in rows
        if row_matches_gene(
            row,
            gene_symbol,
            allow_substring_match=allow_substring_match,
            species=species,
        )
    ]


def summarize_expression_row(
    row: dict[str, Any],
    *,
    species: SpeciesChoice = "human",
) -> dict[str, Any]:
    """Extract gene identifiers plus cell-type expression columns (not evidence)."""
    prefixes = celltype_prefixes_for_species(species)
    cell_types: dict[str, Any] = {}
    for key, value in row.items():
        key_s = str(key)
        if any(key_s.startswith(prefix) for prefix in prefixes):
            cell_types[key_s] = value
    return {
        "gene_id": row.get("gene_id") or row.get("Gene") or row.get("gene"),
        "id": row.get("id"),
        "species": species,
        "cell_type_values": cell_types,
        "cell_type_count": len(cell_types),
    }


def _raw_text_preview(text: str, *, limit: int = RAW_TEXT_PREVIEW_LIMIT) -> str:
    """Short whitespace-collapsed preview for soft-fail diagnostics."""
    collapsed = " ".join((text or "").split())
    if len(collapsed) <= limit:
        return collapsed
    return collapsed[: limit - 3] + "..."


def _safe_response_headers(headers: Any) -> dict[str, str]:
    """Keep a small, non-sensitive response-header subset for debugging."""
    out: dict[str, str] = {}
    if headers is None:
        return out
    for key in SAFE_RESPONSE_HEADER_KEYS:
        value = headers.get(key)
        if value is not None and str(value).strip():
            out[key] = str(value)
    return out


def _should_fallback_to_browser(result: ToolResult) -> bool:
    """True when direct download failed in a way browser retrieval may fix."""
    if result.success:
        # Soft-200 HTML / invalid CSV still needs fallback.
        if isinstance(result.data, dict) and result.data.get("validation_failed"):
            return True
        return False
    if result.status_code == 403 or result.error_type == "access_forbidden":
        return True
    if result.error_type in {
        "invalid_or_html_response",
        "http_error",
        "timeout",
    }:
        return True
    return False


def download_csv(
    species: SpeciesChoice = "human",
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Download a BrainRNASeq CSV bulk file (raw text preserved on success)."""
    cfg = settings or get_settings()
    url = csv_url_for_species(species)
    request_params = {
        "species": species,
        "url": url,
        "request_headers": dict(REQUEST_HEADERS),
        "retrieval_method": "httpx_direct",
    }
    try:
        with httpx.Client(timeout=cfg.http_timeout_seconds) as client:
            response = client.get(url, headers=REQUEST_HEADERS)
        content_type = response.headers.get("content-type")
        response_headers = _safe_response_headers(response.headers)
        if response.is_success:
            validation = validate_brainrnaseq_csv_bytes(response.text, species=species)
            if not validation.get("ok"):
                preview = _raw_text_preview(response.text)
                return _tool_result(
                    endpoint_name="download_csv",
                    gene_symbol=gene_symbol or species,
                    request_url=url,
                    request_params=request_params,
                    success=False,
                    status_code=response.status_code,
                    data={
                        "species": species,
                        "content_type": content_type,
                        "response_headers": response_headers,
                        "raw_text_preview": preview,
                        "validation": validation,
                        "validation_failed": True,
                        "retrieval_method": "httpx_direct",
                    },
                    error_type=str(
                        validation.get("error_type") or "invalid_or_html_response"
                    ),
                    error_message=str(
                        validation.get("error_message")
                        or "BrainRNASeq response failed CSV validation"
                    ),
                )
            return _tool_result(
                endpoint_name="download_csv",
                gene_symbol=gene_symbol or species,
                request_url=url,
                request_params=request_params,
                success=True,
                status_code=response.status_code,
                data={
                    "raw_csv": response.text,
                    "content_type": content_type,
                    "response_headers": response_headers,
                    "species": species,
                    "retrieval_method": "httpx_direct",
                    "final_url": url,
                    "byte_size": validation.get("byte_size"),
                    "sha256": validation.get("sha256"),
                    "validation": validation,
                },
            )

        preview = _raw_text_preview(response.text)
        fail_data = {
            "species": species,
            "content_type": content_type,
            "response_headers": response_headers,
            "raw_text_preview": preview,
            "retrieval_method": "httpx_direct",
        }
        if response.status_code == 403:
            return _tool_result(
                endpoint_name="download_csv",
                gene_symbol=gene_symbol or species,
                request_url=url,
                request_params=request_params,
                success=False,
                status_code=403,
                data=fail_data,
                error_type="access_forbidden",
                error_message=(
                    "BrainRNASeq public CSV endpoint returned HTTP 403 Forbidden; "
                    "the published wp-content CSV appears inaccessible to automated "
                    "clients (access forbidden / unavailable)."
                ),
            )
        return _tool_result(
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params=request_params,
            success=False,
            status_code=response.status_code,
            data=fail_data,
            error_type="http_error",
            error_message=f"HTTP {response.status_code}",
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        return _tool_result(
            endpoint_name="download_csv",
            gene_symbol=gene_symbol or species,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


def _launch_playwright_chromium(pw: Any) -> tuple[Any, list[dict[str, Any]]]:
    """Launch Chrome channel first, then bundled Chromium; keep both attempts."""
    attempts: list[dict[str, Any]] = []
    try:
        browser = pw.chromium.launch(headless=True, channel="chrome")
        attempts.append({"channel": "chrome", "success": True})
        return browser, attempts
    except Exception as chrome_exc:  # noqa: BLE001
        attempts.append(
            {
                "channel": "chrome",
                "success": False,
                "error_type": type(chrome_exc).__name__,
                "error_message": str(chrome_exc)[:400],
            }
        )
    try:
        browser = pw.chromium.launch(headless=True)
        attempts.append({"channel": "chromium", "success": True})
        return browser, attempts
    except Exception as chromium_exc:  # noqa: BLE001
        attempts.append(
            {
                "channel": "chromium",
                "success": False,
                "error_type": type(chromium_exc).__name__,
                "error_message": str(chromium_exc)[:400],
            }
        )
        raise


def _url_matches_published_csv(url: str, expected: str) -> bool:
    """True when ``url`` resolves to the same published CSV as ``expected``."""
    canonical = normalize_published_csv_url(url)
    return canonical is not None and canonical == expected


def download_csv_via_browser(
    species: SpeciesChoice = "human",
    *,
    gene_symbol: str = "",
    settings: Settings | None = None,
) -> ToolResult:
    """Download published BrainRNASeq CSV via Playwright BrowserContext.

    Captures exact CSV bytes from a matching Response body or Download event.
    Never scrapes rendered page text. Never raises.
    """
    _ = settings  # reserved for future timeout/config wiring
    url = csv_url_for_species(species)
    request_params: dict[str, Any] = {
        "species": species,
        "url": url,
        "retrieval_method": "official_browser_download",
        "allowlisted_urls": sorted(published_csv_urls()),
    }
    symbol = gene_symbol or species

    try:
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name="download_csv_browser",
            gene_symbol=symbol,
            request_url=url,
            request_params=request_params,
            success=False,
            error_type="playwright_unavailable",
            error_message=(
                "Playwright is unavailable. Install with "
                ".venv/bin/python -m playwright install chromium. "
                f"{type(exc).__name__}: {exc}"
            ),
            data={"species": species, "browser_launch_attempts": []},
        )

    launch_attempts: list[dict[str, Any]] = []
    browser = None
    try:
        with sync_playwright() as pw:
            try:
                browser, launch_attempts = _launch_playwright_chromium(pw)
            except Exception as launch_exc:  # noqa: BLE001
                return _tool_result(
                    endpoint_name="download_csv_browser",
                    gene_symbol=symbol,
                    request_url=url,
                    request_params={**request_params, "browser_launch_attempts": launch_attempts},
                    success=False,
                    error_type="browser_launch_failed",
                    error_message=str(launch_exc)[:500],
                    data={
                        "species": species,
                        "browser_launch_attempts": launch_attempts,
                    },
                )

            browser_channel = "chromium"
            for attempt in reversed(launch_attempts):
                if attempt.get("success"):
                    browser_channel = str(attempt.get("channel") or "chromium")
                    break

            context = browser.new_context(accept_downloads=True)
            page = context.new_page()
            captured: dict[str, Any] = {
                "bytes": None,
                "final_url": None,
                "content_type": None,
                "capture_via": None,
                "status_code": None,
            }

            def _maybe_capture_response(response: Any) -> None:
                if captured["bytes"] is not None:
                    return
                resp_url = str(getattr(response, "url", "") or "")
                if not _url_matches_published_csv(resp_url, url):
                    return
                try:
                    body = response.body()
                except Exception:  # noqa: BLE001 — body often unavailable for downloads
                    return
                if not isinstance(body, (bytes, bytearray)) or not body:
                    return
                headers = getattr(response, "headers", {}) or {}
                content_type = None
                if hasattr(headers, "get"):
                    content_type = headers.get("content-type")
                observed_status = getattr(response, "status", None)
                try:
                    observed_status = (
                        int(observed_status) if observed_status is not None else None
                    )
                except (TypeError, ValueError):
                    observed_status = None
                captured["bytes"] = bytes(body)
                captured["final_url"] = resp_url
                captured["content_type"] = content_type
                captured["capture_via"] = "response_body"
                captured["status_code"] = observed_status
                captured["http_status_observed"] = observed_status is not None

            page.on("response", _maybe_capture_response)

            # Chrome serves these CSVs as file downloads. page.goto raises
            # "Download is starting"; capture via expect_download + save_as.
            try:
                from playwright.sync_api import Error as PlaywrightError
            except Exception:  # noqa: BLE001
                PlaywrightError = Exception  # type: ignore[misc, assignment]

            try:
                with page.expect_download(timeout=60_000) as download_info:
                    try:
                        page.goto(url, wait_until="commit", timeout=60_000)
                    except PlaywrightError as goto_exc:
                        # Expected when the browser starts a file download.
                        if "Download is starting" not in str(goto_exc):
                            raise
                download = download_info.value
                download_url = str(getattr(download, "url", None) or url)
                # Prefer save_as — download.path() can be canceled under races.
                with tempfile.TemporaryDirectory(prefix="brs_csv_") as tmp:
                    dest = Path(tmp) / "dataset.csv"
                    download.save_as(str(dest))
                    download_bytes = dest.read_bytes()
                if download_bytes and _url_matches_published_csv(download_url, url):
                    captured["bytes"] = download_bytes
                    captured["final_url"] = download_url
                    captured["content_type"] = "text/csv"
                    captured["capture_via"] = "download_event"
                    # Download events do not expose a reliable HTTP status.
                    captured["status_code"] = None
                    captured["http_status_observed"] = False
            except Exception as download_exc:  # noqa: BLE001
                logger.info(
                    "BrainRNASeq download-event capture failed (%s); "
                    "trying navigation response body",
                    download_exc,
                )
                if captured["bytes"] is None:
                    try:
                        response = page.goto(url, wait_until="commit", timeout=60_000)
                        if response is not None and _url_matches_published_csv(
                            str(response.url or ""), url
                        ):
                            body = response.body()
                            if body:
                                captured["bytes"] = bytes(body)
                                captured["final_url"] = str(response.url)
                                headers = getattr(response, "headers", {}) or {}
                                if hasattr(headers, "get"):
                                    captured["content_type"] = headers.get(
                                        "content-type"
                                    )
                                captured["capture_via"] = "navigation_response_body"
                                observed_status = getattr(response, "status", None)
                                try:
                                    captured["status_code"] = (
                                        int(observed_status)
                                        if observed_status is not None
                                        else None
                                    )
                                except (TypeError, ValueError):
                                    captured["status_code"] = None
                                captured["http_status_observed"] = (
                                    captured["status_code"] is not None
                                )
                    except Exception:  # noqa: BLE001
                        pass

            final_url = str(captured.get("final_url") or page.url or "")
            context.close()
            browser.close()
            browser = None

            raw = captured.get("bytes")
            if not isinstance(raw, (bytes, bytearray)) or not raw:
                return _tool_result(
                    endpoint_name="download_csv_browser",
                    gene_symbol=symbol,
                    request_url=url,
                    request_params={
                        **request_params,
                        "browser_launch_attempts": launch_attempts,
                        "browser_channel": browser_channel,
                    },
                    success=False,
                    error_type="browser_download_failed",
                    error_message=(
                        "Playwright did not capture CSV bytes from Response "
                        "body or Download event"
                    ),
                    data={
                        "species": species,
                        "browser_channel": browser_channel,
                        "browser_launch_attempts": launch_attempts,
                        "final_url": final_url,
                        "retrieval_method": "official_browser_download",
                    },
                )

            canonical_final = normalize_published_csv_url(final_url) or (
                normalize_published_csv_url(url) if _url_matches_published_csv(final_url or url, url) else None
            )
            # Prefer the response URL if it matched; otherwise require page URL.
            if canonical_final is None:
                # Still accept when capture matched the expected published URL.
                if _url_matches_published_csv(str(captured.get("final_url") or ""), url):
                    canonical_final = url
                else:
                    return _tool_result(
                        endpoint_name="download_csv_browser",
                        gene_symbol=symbol,
                        request_url=url,
                        request_params={
                            **request_params,
                            "browser_launch_attempts": launch_attempts,
                            "browser_channel": browser_channel,
                        },
                        success=False,
                        error_type="url_not_allowlisted",
                        error_message=(
                            f"Browser final URL is not an allowlisted CSV: {final_url!r}"
                        ),
                        data={
                            "species": species,
                            "browser_channel": browser_channel,
                            "browser_launch_attempts": launch_attempts,
                            "final_url": final_url,
                            "retrieval_method": "official_browser_download",
                        },
                    )

            validation = validate_brainrnaseq_csv_bytes(bytes(raw), species=species)
            if not validation.get("ok"):
                preview = _raw_text_preview(
                    bytes(raw).decode("utf-8", errors="replace")
                )
                return _tool_result(
                    endpoint_name="download_csv_browser",
                    gene_symbol=symbol,
                    request_url=url,
                    request_params={
                        **request_params,
                        "browser_launch_attempts": launch_attempts,
                        "browser_channel": browser_channel,
                    },
                    success=False,
                    error_type=str(
                        validation.get("error_type") or "invalid_or_html_response"
                    ),
                    error_message=str(
                        validation.get("error_message")
                        or "Browser CSV failed validation"
                    ),
                    data={
                        "species": species,
                        "browser_channel": browser_channel,
                        "browser_launch_attempts": launch_attempts,
                        "final_url": canonical_final,
                        "content_type": captured.get("content_type"),
                        "raw_text_preview": preview,
                        "validation": validation,
                        "retrieval_method": "official_browser_download",
                        "capture_via": captured.get("capture_via"),
                    },
                )

            raw_csv = bytes(raw).decode("utf-8-sig")
            capture_via = str(captured.get("capture_via") or "")
            if capture_via == "download_event":
                observed_status: int | None = None
                http_status_observed = False
            else:
                observed_status = captured.get("status_code")
                try:
                    observed_status = (
                        int(observed_status) if observed_status is not None else None
                    )
                except (TypeError, ValueError):
                    observed_status = None
                http_status_observed = observed_status is not None
            return _tool_result(
                endpoint_name="download_csv_browser",
                gene_symbol=symbol,
                request_url=url,
                request_params={
                    **request_params,
                    "browser_launch_attempts": launch_attempts,
                    "browser_channel": browser_channel,
                },
                success=True,
                status_code=observed_status,
                data={
                    "raw_csv": raw_csv,
                    "content_type": captured.get("content_type") or "text/csv",
                    "species": species,
                    "retrieval_method": "official_browser_download",
                    "source_url": url,
                    "final_url": canonical_final,
                    "browser_channel": browser_channel,
                    "browser_launch_attempts": launch_attempts,
                    "byte_size": validation.get("byte_size"),
                    "sha256": validation.get("sha256"),
                    "validation": validation,
                    "capture_via": capture_via,
                    "http_status_observed": http_status_observed,
                },
            )
    except Exception as exc:  # noqa: BLE001 — clients must never raise
        logger.warning("BrainRNASeq browser download failed: %s", exc)
        return _tool_result(
            endpoint_name="download_csv_browser",
            gene_symbol=symbol,
            request_url=url,
            request_params={**request_params, "browser_launch_attempts": launch_attempts},
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc)[:500],
            data={
                "species": species,
                "browser_launch_attempts": launch_attempts,
                "retrieval_method": "official_browser_download",
            },
        )
    finally:
        if browser is not None:
            try:
                browser.close()
            except Exception:  # noqa: BLE001
                pass


def expression_from_raw_csv(
    raw_csv: str,
    gene_symbol: str,
    *,
    species: SpeciesChoice = "human",
    allow_substring_match: bool = False,
    request_url: str | None = None,
    content_type: str | None = None,
    status_code: int | None = 200,
) -> ToolResult:
    """Filter an already-downloaded BrainRNASeq CSV (no HTTP).

    Success ``data`` matches :func:`fetch_gene_expression` shape so Section 2b
    can download once per species via :func:`download_csv`, persist the raw
    artifact, then filter locally without re-downloading.
    """
    symbol = gene_symbol.strip()
    url = request_url or csv_url_for_species(species)
    if not symbol:
        return _tool_result(
            endpoint_name="expression_from_raw_csv",
            gene_symbol=gene_symbol,
            request_url=url,
            request_params={"species": species},
            success=False,
            error_type="invalid_request",
            error_message="gene_symbol is required",
        )
    rows = parse_csv(raw_csv)
    matched = filter_gene_rows(
        rows,
        symbol,
        allow_substring_match=allow_substring_match,
        species=species,
    )
    summaries = [
        summarize_expression_row(row, species=species) for row in matched
    ]
    return _tool_result(
        endpoint_name="expression_from_raw_csv",
        gene_symbol=symbol,
        request_url=url,
        request_params={
            "species": species,
            "gene_symbol": symbol,
            "allow_substring_match": allow_substring_match,
            "from_raw_csv": True,
        },
        success=True,
        status_code=status_code,
        data={
            "gene_symbol": symbol,
            "species": species,
            "raw_csv": raw_csv,
            "content_type": content_type,
            "matched_rows": matched,
            "expression_summaries": summaries,
            "match_count": len(matched),
            "row_count_total": len(rows),
        },
    )


def fetch_gene_expression(
    gene_symbol: str,
    *,
    species: SpeciesChoice = "human",
    allow_substring_match: bool = False,
    settings: Settings | None = None,
) -> ToolResult:
    """Download CSV and return rows matching ``gene_symbol``.

    On success, ``data`` includes raw CSV, matched rows, and light summaries.
    Never raises.
    """
    cfg = settings or get_settings()
    symbol = gene_symbol.strip()
    if not symbol:
        return _tool_result(
            endpoint_name="fetch_gene_expression",
            gene_symbol=gene_symbol,
            request_url=csv_url_for_species(species),
            request_params={"species": species},
            success=False,
            error_type="invalid_request",
            error_message="gene_symbol is required",
        )

    downloaded = download_csv(species, gene_symbol=symbol, settings=cfg)
    if not downloaded.success:
        return _tool_result(
            endpoint_name="fetch_gene_expression",
            gene_symbol=symbol,
            request_url=downloaded.request_url,
            request_params=downloaded.request_params,
            success=False,
            status_code=downloaded.status_code,
            data=downloaded.data,
            error_type=downloaded.error_type or "download_failed",
            error_message=downloaded.error_message or "BrainRNASeq CSV download failed",
        )

    raw_csv = ""
    content_type = None
    if isinstance(downloaded.data, dict):
        raw_csv = str(downloaded.data.get("raw_csv") or "")
        content_type = downloaded.data.get("content_type")
    filtered = expression_from_raw_csv(
        raw_csv,
        symbol,
        species=species,
        allow_substring_match=allow_substring_match,
        request_url=downloaded.request_url,
        content_type=content_type if isinstance(content_type, str) else None,
        status_code=downloaded.status_code,
    )
    # Preserve fetch_gene_expression endpoint_name for existing callers/tests.
    return filtered.model_copy(
        update={
            "endpoint_name": "fetch_gene_expression",
            "request_params": {
                "species": species,
                "gene_symbol": symbol,
                "allow_substring_match": allow_substring_match,
            },
        }
    )


__all__ = [
    "SOURCE_NAME",
    "HUMAN_CSV_URL",
    "MOUSE_CSV_URL",
    "REQUEST_HEADERS",
    "HUMAN_CELLTYPE_PREFIXES",
    "MOUSE_CELLTYPE_PREFIXES",
    "MIN_CSV_BYTE_SIZE",
    "csv_url_for_species",
    "published_csv_urls",
    "normalize_published_csv_url",
    "is_html_payload",
    "validate_brainrnaseq_csv_bytes",
    "parse_csv",
    "normalize_gene_identifier",
    "row_matches_gene",
    "filter_gene_rows",
    "summarize_expression_row",
    "download_csv",
    "download_csv_via_browser",
    "_should_fallback_to_browser",
    "expression_from_raw_csv",
    "fetch_gene_expression",
]
