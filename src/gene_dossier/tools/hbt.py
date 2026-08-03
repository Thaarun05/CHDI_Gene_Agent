"""Human Brain Transcriptome (HBT) whole-brain PDF client.

Fetches the official gene-specific developmental-expression PDF from::

    https://hbatlas.org/hbtd/images/wholeBrain/{GENE}.pdf

Does **not** normalize into evidence records. Never raises: failures return
:class:`~gene_dossier.models.ToolResult`.
"""

from __future__ import annotations

from typing import Any
from urllib.parse import quote

import httpx

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ToolResult

SOURCE_NAME = "Human Brain Transcriptome"
HBT_BASE = "https://hbatlas.org/hbtd/images/wholeBrain"
HBT_HOME = "https://hbatlas.org/"

REQUEST_HEADERS = {
    "User-Agent": "GeneDossier/0.1.0 (research; provenance-first gene dossier client)",
    "Accept": "application/pdf,*/*",
}

_MIN_PDF_BYTES = 500


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


def whole_brain_pdf_url(gene_symbol: str) -> str:
    """Build the official HBT whole-brain PDF URL for ``gene_symbol``."""
    symbol = (gene_symbol or "").strip()
    return f"{HBT_BASE}/{quote(symbol, safe='')}.pdf"


def extract_pdf_text(pdf_bytes: bytes, *, max_pages: int | None = None) -> list[str]:
    """Extract plain text per page via PyMuPDF when available. Never uses OCR."""
    try:
        import fitz  # type: ignore
    except ImportError:
        return []
    texts: list[str] = []
    try:
        doc = fitz.open(stream=pdf_bytes, filetype="pdf")
    except Exception:  # noqa: BLE001
        return []
    try:
        limit = len(doc) if max_pages is None else min(len(doc), max_pages)
        for index in range(limit):
            try:
                texts.append(doc.load_page(index).get_text("text") or "")
            except Exception:  # noqa: BLE001
                texts.append("")
    finally:
        doc.close()
    return texts


def select_plot_page(
    page_texts: list[str],
    *,
    gene_symbol: str,
) -> int:
    """Pick the page most likely to contain the gene plot (0-based)."""
    if not page_texts:
        return 0
    target = (gene_symbol or "").strip().upper()
    scored: list[tuple[int, int]] = []
    for index, text in enumerate(page_texts):
        upper = (text or "").upper()
        score = 0
        if target and target in upper:
            score += 10
        for token in ("NCX", "HIP", "AMY", "STR", "CBC", "PERIOD", "AGE"):
            if token in upper:
                score += 1
        scored.append((score, index))
    scored.sort(key=lambda item: (-item[0], item[1]))
    best_score, best_index = scored[0]
    if best_score <= 0:
        return 0
    return best_index


def rasterize_pdf_page(
    pdf_bytes: bytes,
    *,
    page_index: int = 0,
    dpi: int = 180,
) -> tuple[bytes, dict[str, Any]] | None:
    """Rasterize one PDF page to PNG bytes. Returns ``None`` on failure."""
    try:
        import fitz  # type: ignore
    except ImportError:
        return None
    try:
        doc = fitz.open(stream=pdf_bytes, filetype="pdf")
    except Exception:  # noqa: BLE001
        return None
    try:
        if page_index < 0 or page_index >= len(doc):
            return None
        page = doc.load_page(page_index)
        zoom = max(dpi, 72) / 72.0
        matrix = fitz.Matrix(zoom, zoom)
        pix = page.get_pixmap(matrix=matrix, alpha=False)
        png = pix.tobytes("png")
        return png, {
            "width": int(pix.width),
            "height": int(pix.height),
            "dpi": int(dpi),
            "page_index": int(page_index),
            "page_count": int(len(doc)),
            "media_type": "image/png",
            "byte_size": len(png),
        }
    except Exception:  # noqa: BLE001
        return None
    finally:
        doc.close()


def fetch_whole_brain_pdf(
    gene_symbol: str,
    *,
    settings: Settings | None = None,
) -> ToolResult:
    """HTTP GET the official HBT whole-brain PDF for ``gene_symbol``."""
    cfg = settings or get_settings()
    symbol = (gene_symbol or "").strip()
    url = whole_brain_pdf_url(symbol)
    params = {"gene_symbol": symbol}
    try:
        with httpx.Client(
            timeout=cfg.http_timeout_seconds,
            follow_redirects=True,
            headers=REQUEST_HEADERS,
        ) as client:
            response = client.get(url)
        content = response.content or b""
        content_type = (response.headers.get("content-type") or "").lower()
        if not response.is_success:
            return _tool_result(
                endpoint_name="fetch_whole_brain_pdf",
                gene_symbol=symbol,
                request_url=str(response.url),
                request_params=params,
                success=False,
                status_code=response.status_code,
                data={
                    "content_type": content_type,
                    "byte_size": len(content),
                    "source_url": url,
                },
                error_type="http_error",
                error_message=f"HTTP {response.status_code}",
            )
        if not content.startswith(b"%PDF"):
            return _tool_result(
                endpoint_name="fetch_whole_brain_pdf",
                gene_symbol=symbol,
                request_url=str(response.url),
                request_params=params,
                success=False,
                status_code=response.status_code,
                data={
                    "content_type": content_type,
                    "byte_size": len(content),
                    "source_url": url,
                    "magic_prefix": content[:8].decode("latin-1", errors="replace"),
                },
                error_type="invalid_pdf",
                error_message="Response does not begin with %PDF",
            )
        if len(content) < _MIN_PDF_BYTES:
            return _tool_result(
                endpoint_name="fetch_whole_brain_pdf",
                gene_symbol=symbol,
                request_url=str(response.url),
                request_params=params,
                success=False,
                status_code=response.status_code,
                data={"byte_size": len(content), "source_url": url},
                error_type="invalid_pdf",
                error_message=f"PDF too small ({len(content)} bytes)",
            )

        page_texts = extract_pdf_text(content)
        gene_text_found: bool | None = None
        if page_texts:
            joined = "\n".join(page_texts).upper()
            gene_text_found = bool(symbol) and symbol.upper() in joined
            if gene_text_found is False:
                return _tool_result(
                    endpoint_name="fetch_whole_brain_pdf",
                    gene_symbol=symbol,
                    request_url=str(response.url),
                    request_params=params,
                    success=False,
                    status_code=response.status_code,
                    data={
                        "byte_size": len(content),
                        "source_url": url,
                        "page_count": len(page_texts),
                        "gene_text_found": False,
                    },
                    error_type="gene_text_mismatch",
                    error_message=f"Gene symbol {symbol!r} not found in PDF text",
                )

        plot_page = select_plot_page(page_texts, gene_symbol=symbol)
        return _tool_result(
            endpoint_name="fetch_whole_brain_pdf",
            gene_symbol=symbol,
            request_url=str(response.url),
            request_params=params,
            success=True,
            status_code=response.status_code,
            data={
                "gene_symbol": symbol,
                "source_url": url,
                "content_type": content_type or "application/pdf",
                "byte_size": len(content),
                "pdf_bytes": content,
                "page_count": len(page_texts) if page_texts else None,
                "page_texts_available": bool(page_texts),
                "gene_text_found": gene_text_found,
                "selected_page_index": plot_page,
            },
        )
    except httpx.TimeoutException as exc:
        return _tool_result(
            endpoint_name="fetch_whole_brain_pdf",
            gene_symbol=symbol,
            request_url=url,
            request_params=params,
            success=False,
            error_type="timeout",
            error_message=str(exc),
        )
    except httpx.HTTPError as exc:
        return _tool_result(
            endpoint_name="fetch_whole_brain_pdf",
            gene_symbol=symbol,
            request_url=url,
            request_params=params,
            success=False,
            error_type="http_error",
            error_message=str(exc),
        )
    except Exception as exc:  # noqa: BLE001
        return _tool_result(
            endpoint_name="fetch_whole_brain_pdf",
            gene_symbol=symbol,
            request_url=url,
            request_params=params,
            success=False,
            error_type=type(exc).__name__,
            error_message=str(exc),
        )


__all__ = [
    "SOURCE_NAME",
    "HBT_BASE",
    "HBT_HOME",
    "whole_brain_pdf_url",
    "extract_pdf_text",
    "select_plot_page",
    "rasterize_pdf_page",
    "fetch_whole_brain_pdf",
]
