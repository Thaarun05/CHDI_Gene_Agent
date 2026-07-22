"""Rancho BioSciences / CHDI polished dossier renderer.

Consumes a :class:`~gene_dossier.report_schema.ReportDocument` and writes a
visual report that follows ``SREBF2_report.pdf``:

- cover page (gene + CHR, Gene Report, prepared-for, curator/date)
- table of contents
- 15 major sections (green) with lettered subsections (orange)
- tables with light-green header rows
- embedded figures when ``figure_path`` is set on a block
- References + Compiled List of Relevant Databases
- optional provenance endnotes (``source_id`` citations)

This is the **final** report format. ``rendering.py`` remains a provenance/debug
markdown view only.

HTML is always produced. PDF export uses PyMuPDF ``Story`` when available and
soft-fails to HTML-only otherwise (pymupdf is optional at runtime).
"""

from __future__ import annotations

import base64
import html
import json
import logging
import mimetypes
from pathlib import Path
from typing import Any, Iterable

from gene_dossier.config import Settings, get_settings
from gene_dossier.report_schema import (
    REPORT_STYLE,
    ReportContentBlock,
    ReportDocument,
    ReportMajorSection,
    ReportSubsection,
    cover_lines,
    iter_toc_entries,
)

logger = logging.getLogger(__name__)

_ASSETS_DIR = Path(__file__).resolve().parent / "assets"


def _escape(text: str | None) -> str:
    return html.escape(text or "", quote=True)


def _asset_data_uri(name: str) -> str | None:
    """Load a packaged asset as a data URI for self-contained HTML."""
    path = _ASSETS_DIR / name
    if not path.is_file():
        return None
    mime, _ = mimetypes.guess_type(path.name)
    mime = mime or "application/octet-stream"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def _rancho_css() -> str:
    s = REPORT_STYLE
    # Prefer literal hex colors (not only CSS variables) so PyMuPDF Story PDF
    # export keeps Rancho branding when var() support is limited.
    return f"""
* {{ box-sizing: border-box; }}
html, body {{
  margin: 0;
  padding: 0;
  color: {s.brown_body};
  font-family: Arial, Helvetica, sans-serif;
  font-size: {s.body_pt}pt;
  line-height: 1.45;
  background: #ffffff;
}}
a {{ color: {s.orange_link}; text-decoration: none; }}
.page-header {{
  padding: 12pt 0 6pt 0;
  border-bottom: 1.5pt solid {s.rule_green};
  margin-bottom: 18pt;
}}
.page-header img.rancho {{ height: 42px; }}
.page-header img.chdi {{ height: 64px; float: right; }}
.page-footer {{
  margin-top: 28pt;
  padding-top: 8pt;
  border-top: 1pt solid {s.rule_green};
  font-size: 9pt;
  color: {s.green_major};
}}
.page-footer img {{ height: 28px; vertical-align: middle; }}
.cover {{
  min-height: 90vh;
  page-break-after: always;
}}
.cover-logos {{
  margin-bottom: 48pt;
}}
.cover-logos img.rancho {{ height: 72px; }}
.cover-logos img.chdi {{ height: 110px; float: right; }}
.cover-title-block {{
  text-align: center;
  margin: 72pt 0 72pt 0;
  clear: both;
}}
.cover-title {{
  color: {s.green_major};
  font-size: {s.cover_title_pt}pt;
  font-weight: 700;
  margin: 0 0 8pt 0;
}}
.cover-subtitle {{
  color: {s.green_major};
  font-size: {s.cover_subtitle_pt}pt;
  font-weight: 500;
  margin: 0;
}}
.cover-meta {{
  color: {s.brown_body};
  font-size: 16pt;
  line-height: 1.5;
  clear: both;
  margin-top: 96pt;
}}
.cover-meta .date {{
  color: {s.orange_sub};
  margin-top: 6pt;
}}
.toc {{
  page-break-after: always;
  clear: both;
}}
.toc h1 {{
  color: {s.brown_body};
  font-size: 22pt;
  margin: 0 0 16pt 0;
}}
table.toc-table {{
  width: 100%;
  border-collapse: collapse;
}}
table.toc-table td {{
  padding: 2pt 0;
  vertical-align: bottom;
  border: none;
}}
td.toc-leader {{
  width: 100%;
  border-bottom: 1px dotted #8a8a8a;
}}
td.toc-page {{
  text-align: right;
  padding-left: 8pt;
  white-space: nowrap;
  font-variant-numeric: tabular-nums;
}}
tr.toc-entry.major td {{
  color: {s.green_major};
  font-weight: 700;
  text-transform: uppercase;
  padding-top: 8pt;
}}
tr.toc-entry.major a {{ color: {s.green_major}; }}
tr.toc-entry.sub td.toc-title {{
  padding-left: 18pt;
}}
tr.toc-entry.sub td {{
  color: {s.orange_sub};
  font-weight: 600;
  text-transform: uppercase;
  font-size: 10.5pt;
}}
tr.toc-entry.sub a {{ color: {s.orange_sub}; }}
tr.toc-entry.back td {{
  color: {s.brown_body};
  font-weight: 700;
  text-transform: uppercase;
  padding-top: 10pt;
}}
tr.toc-entry.back a {{ color: {s.brown_body}; }}
.major-heading {{
  color: {s.green_major};
  font-size: {s.major_header_pt}pt;
  font-weight: 700;
  margin: 22pt 0 10pt 0;
}}
.sub-heading {{
  color: {s.orange_sub};
  font-size: {s.subsection_header_pt}pt;
  font-weight: 600;
  margin: 14pt 0 8pt 0;
}}
.block {{ margin: 0 0 10pt 0; }}
.empty-note {{
  font-style: italic;
  color: #7a5a45;
}}
table.rancho {{
  border-collapse: collapse;
  width: 100%;
  margin: 8pt 0 12pt 0;
  font-size: 10.5pt;
}}
table.rancho th {{
  background-color: {s.table_header_bg};
  color: {s.brown_body};
  font-weight: 700;
  text-align: left;
  padding: 6pt 8pt;
  border: 1px solid #b0b0b0;
}}
table.rancho td {{
  padding: 5pt 8pt;
  border: 1px solid #c8c8c8;
  vertical-align: top;
}}
figure.rancho-figure {{
  margin: 10pt 0 14pt 0;
  text-align: center;
}}
figure.rancho-figure img {{
  max-width: 100%;
  height: auto;
}}
.back-matter h2, .endnotes h2 {{
  color: {s.green_major};
  font-size: 18pt;
}}
.db-list, .ref-list, .endnote-list {{ padding-left: 18pt; }}
.provenance-muted {{ color: #8a7a6a; font-size: 8.5pt; }}
""".strip()


def _img_tag(data_uri: str | None, *, cls: str, alt: str) -> str:
    if not data_uri:
        return f'<span class="{_escape(cls)}-fallback">{_escape(alt)}</span>'
    return f'<img class="{_escape(cls)}" src="{data_uri}" alt="{_escape(alt)}" />'


def _render_block(block: ReportContentBlock) -> str:
    parts: list[str] = ['<div class="block">']
    if block.title:
        parts.append(
            f'<p><strong style="color:{REPORT_STYLE.orange_link};">'
            f"{_escape(block.title)}</strong></p>"
        )

    if block.kind == "table" and block.table_headers:
        parts.append('<table class="rancho"><thead><tr>')
        for h in block.table_headers:
            parts.append(
                f'<th bgcolor="{REPORT_STYLE.table_header_bg}" '
                f'style="background-color:{REPORT_STYLE.table_header_bg};">'
                f"{_escape(h)}</th>"
            )
        parts.append("</tr></thead><tbody>")
        for row in block.table_rows:
            parts.append("<tr>")
            for cell in row:
                parts.append(f"<td>{_escape(cell)}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table>")
    elif block.kind == "figure" and block.figure_path:
        fig_uri = None
        fig_path = Path(block.figure_path)
        if fig_path.is_file():
            mime, _ = mimetypes.guess_type(fig_path.name)
            mime = mime or "image/png"
            fig_uri = (
                f"data:{mime};base64,"
                + base64.b64encode(fig_path.read_bytes()).decode("ascii")
            )
        caption = _escape(block.figure_caption or "")
        if fig_uri:
            parts.append('<figure class="rancho-figure">')
            parts.append(f'<img src="{fig_uri}" alt="{caption or "figure"}" />')
            if caption:
                parts.append(f"<figcaption>{caption}</figcaption>")
            parts.append("</figure>")
        elif block.text:
            parts.append(f'<div class="narrative"><p>{_escape(block.text)}</p></div>')
    elif block.kind == "list" and block.table_rows:
        parts.append("<ul>")
        for row in block.table_rows:
            parts.append(f"<li>{_escape(' — '.join(row))}</li>")
        parts.append("</ul>")
    elif block.kind == "link" and block.links:
        parts.append("<ul>")
        for link in block.links:
            label = _escape(link.get("label") or link.get("url") or "link")
            url = _escape(link.get("url") or "#")
            parts.append(f'<li><a href="{url}">{label}</a></li>')
        parts.append("</ul>")
    else:
        text = (block.text or "").strip()
        if text:
            paras = text.split("\n\n")
            parts.append('<div class="narrative">')
            for para in paras:
                parts.append(f"<p>{_escape(para.strip())}</p>")
            parts.append("</div>")

    if block.links and block.kind not in {"link"}:
        for link in block.links:
            label = _escape(link.get("label") or "Link")
            url = _escape(link.get("url") or "#")
            parts.append(f'<p><a href="{url}">{label}</a></p>')

    parts.append("</div>")
    return "\n".join(parts)


def _render_subsection(sub: ReportSubsection) -> str:
    heading = f"{sub.key}. {sub.title}"
    parts = [
        (
            f'<h3 class="sub-heading" style="color:{REPORT_STYLE.orange_sub};">'
            f"{_escape(heading)}</h3>"
        ),
    ]
    if not sub.blocks:
        parts.append(
            '<p class="empty-note">No evidence available for this subsection in the '
            "current dossier run.</p>"
        )
    else:
        for block in sub.blocks:
            parts.append(_render_block(block))
    return "\n".join(parts)


def _render_major(section: ReportMajorSection) -> str:
    heading = f"{section.number}. {section.title}"
    parts = [
        f'<section id="section-{section.number}">',
        (
            f'<h2 class="major-heading" style="color:{REPORT_STYLE.green_major};">'
            f"{_escape(heading)}</h2>"
        ),
    ]
    if section.subsections:
        for sub in section.subsections:
            parts.append(_render_subsection(sub))
    elif section.blocks:
        for block in section.blocks:
            parts.append(_render_block(block))
    else:
        parts.append(
            '<p class="empty-note">No evidence available for this section in the '
            "current dossier run.</p>"
        )
    parts.append("</section>")
    return "\n".join(parts)


def _collect_endnotes(doc: ReportDocument) -> list[dict[str, str]]:
    """Build provenance endnotes from block source_ids (order preserved)."""
    notes: list[dict[str, str]] = []
    seen: set[str] = set()
    for section in doc.sections:
        blocks: list[ReportContentBlock] = list(section.blocks)
        for sub in section.subsections:
            blocks.extend(sub.blocks)
        for block in blocks:
            for sid in block.source_ids:
                if not sid or sid in seen:
                    continue
                seen.add(sid)
                text = (block.text or block.title or "").strip()
                if len(text) > 120:
                    text = text[:117] + "…"
                notes.append({"source_id": sid, "excerpt": text})
    return notes


def _render_cover(
    doc: ReportDocument,
    rancho_uri: str | None,
    chdi_uri: str | None,
    *,
    show_cover_logos: bool = False,
) -> str:
    cover = doc.cover
    lines = cover_lines(cover)
    meta_html: list[str] = []
    for line in lines[2:]:
        if not line:
            continue
        cls = "date" if line == (cover.report_date or "") else ""
        style = (
            f' style="color:{REPORT_STYLE.orange_sub};"'
            if cls == "date"
            else ""
        )
        meta_html.append(f'<div class="{cls}"{style}>{_escape(line)}</div>')
    logos = ""
    if show_cover_logos:
        logos = (
            '<div class="cover-logos">'
            f'{_img_tag(rancho_uri, cls="rancho", alt="Rancho BioSciences")}'
            f'{_img_tag(chdi_uri, cls="chdi", alt="CHDI Foundation")}'
            "</div>"
        )
    return f"""
<section class="cover">
  {logos}
  <div class="cover-title-block">
    <h1 class="cover-title" style="color:{REPORT_STYLE.green_major};text-align:center;">
      {_escape(cover.gene_line)}
    </h1>
    <p class="cover-subtitle" style="color:{REPORT_STYLE.green_major};text-align:center;">
      {_escape(cover.title_line)}
    </p>
  </div>
  <div class="cover-meta" style="color:{REPORT_STYLE.brown_body};">
    {"".join(meta_html)}
  </div>
</section>
""".strip()


def _toc_entry_key(entry: dict[str, Any]) -> str:
    """Stable lookup key for ``toc_page_numbers`` per TOC entry."""
    kind = entry["kind"]
    if kind == "major":
        return entry["key"]
    if kind == "subsection":
        return entry["slot"]
    return entry["display_title"].lower().replace(" ", "-")


def _render_toc(
    _doc: ReportDocument,
    toc_page_numbers: dict[str, int] | None = None,
) -> str:
    """TOC with dotted leaders and page-number slots.

    ``toc_page_numbers`` maps a TOC key (major key ``"1"``, subsection slot
    ``"1a"``, or back-matter slug) to a page number. Missing entries render a
    blank page slot while keeping the dotted-leader visual format.
    """
    pages = toc_page_numbers or {}
    entries = iter_toc_entries()
    rows: list[str] = []
    for entry in entries:
        kind = entry["kind"]
        if kind == "major":
            title = f'{entry["number"]}. {entry["title"]}'
            href = f'#section-{entry["number"]}'
            css = "major"
        elif kind == "subsection":
            title = f'{entry["key"].upper()}. {entry["title"]}'
            href = f'#section-{entry["number"]}'
            css = "sub"
        else:
            title = entry["title"]
            slug = entry["display_title"].lower().replace(" ", "-")
            href = f"#{slug}"
            css = "back"

        page_val = pages.get(_toc_entry_key(entry))
        page_str = _escape(str(page_val)) if page_val is not None else ""
        rows.append(
            f'<tr class="toc-entry {css}">'
            f'<td class="toc-title"><a href="{href}">{_escape(title)}</a></td>'
            f'<td class="toc-leader"></td>'
            f'<td class="toc-page">{page_str}</td>'
            f"</tr>"
        )
    return (
        '<section class="toc">'
        "<h1>Table of Contents</h1>"
        '<table class="toc-table">'
        + "\n".join(rows)
        + "</table>"
        + "</section>"
    )


def _render_back_matter(doc: ReportDocument) -> str:
    parts: list[str] = ['<section class="back-matter">']
    parts.append('<h2 id="references">References</h2>')
    if doc.references:
        parts.append('<ol class="ref-list">')
        for ref in doc.references:
            parts.append(f"<li>{_escape(ref)}</li>")
        parts.append("</ol>")
    else:
        parts.append(
            '<p class="empty-note">References will be listed here when literature '
            "citations are attached to the dossier run.</p>"
        )

    parts.append(
        '<h2 id="compiled-list-of-relevant-databases">'
        "Compiled List of Relevant Databases</h2>"
    )
    parts.append('<ol class="db-list">')
    for item in doc.compiled_databases:
        name = _escape(item.get("name") or "")
        url = _escape(item.get("url") or "")
        parts.append(f'<li>{name}: <a href="{url}">{url}</a></li>')
    parts.append("</ol>")
    parts.append("</section>")
    return "\n".join(parts)


def _render_endnotes(doc: ReportDocument) -> str:
    notes = _collect_endnotes(doc)
    parts = [
        '<section class="endnotes" id="provenance-endnotes">',
        "<h2>Provenance endnotes</h2>",
        '<p class="provenance-muted">Internal source_id citations retained for '
        "traceability. These are not shown as inline debug markers in the body.</p>",
    ]
    if not notes:
        parts.append('<p class="empty-note">No cited source_ids in this report.</p>')
    else:
        parts.append('<ol class="endnote-list">')
        for note in notes:
            parts.append(
                f"<li><code>{_escape(note['source_id'])}</code>"
                f" — {_escape(note['excerpt'])}</li>"
            )
        parts.append("</ol>")
    if doc.unmapped_source_ids:
        parts.append("<h3>Unmapped evidence source_ids</h3><ul>")
        for sid in doc.unmapped_source_ids:
            parts.append(f"<li><code>{_escape(sid)}</code></li>")
        parts.append("</ul>")
    parts.append("</section>")
    return "\n".join(parts)


def render_rancho_html(
    doc: ReportDocument,
    *,
    include_endnotes: bool = False,
    include_page_chrome: bool = True,
    toc_page_numbers: dict[str, int] | None = None,
    show_cover_logos: bool = False,
) -> str:
    """Render a full Rancho/CHDI-style HTML dossier from ``doc``.

    Provenance endnotes are off by default so the polished report matches the
    reference PDF (ends at References + Compiled Databases). ``source_id``s
    remain in the JSON sidecar regardless.
    """
    rancho = _asset_data_uri("rancho_wordmark.png")
    chdi = _asset_data_uri("chdi_wordmark.png")
    rancho_header = _asset_data_uri("rancho_header_bar.png") or rancho
    rancho_footer = _asset_data_uri("rancho_footer.png")

    header = ""
    footer = ""
    if include_page_chrome:
        header = (
            '<div class="page-header">'
            f'{_img_tag(rancho_header, cls="rancho", alt="Rancho BioSciences")}'
            f'{_img_tag(chdi, cls="chdi", alt="CHDI Foundation")}'
            "</div>"
        )
        footer = (
            '<div class="page-footer">'
            f'{_img_tag(rancho_footer, cls="rancho", alt="Rancho BioSciences")}'
            f"<span>{_escape(REPORT_STYLE.footer_url)}</span>"
            "</div>"
        )

    body_sections = "\n".join(_render_major(sec) for sec in doc.sections)
    endnotes = _render_endnotes(doc) if include_endnotes else ""
    title = _escape(f"{doc.cover.gene_line} Gene Report")
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title}</title>
  <style>
{_rancho_css()}
  </style>
</head>
<body>
{_render_cover(doc, rancho, chdi, show_cover_logos=show_cover_logos)}
{header}
{_render_toc(doc, toc_page_numbers)}
{body_sections}
{_render_back_matter(doc)}
{endnotes}
{footer}
</body>
</html>
"""


def _hex_to_rgb01(value: str) -> tuple[float, float, float]:
    """Convert ``#RRGGBB`` to a 0..1 RGB tuple for PyMuPDF drawing calls."""
    v = (value or "").lstrip("#")
    if len(v) != 6:
        return (0.0, 0.0, 0.0)
    return tuple(int(v[i : i + 2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore[return-value]


def _stamp_pdf_chrome(
    pdf_path: str | Path,
    *,
    stamp_cover: bool = False,
    start_page: int = 1,
) -> bool:
    """Reopen a rendered PDF and stamp repeated Rancho page chrome.

    For every page at/after ``start_page`` (and the cover only when
    ``stamp_cover`` is True) this stamps the Rancho logo/header, a green rule,
    and a footer with the site URL and page number, matching the repeated page
    chrome in ``SREBF2_report.pdf``. Soft-fails to the unstamped PDF.
    """
    try:
        import fitz
    except Exception as exc:  # noqa: BLE001
        logger.warning("Page stamping skipped (pymupdf not importable): %s", exc)
        return False

    dest = Path(pdf_path)
    if not dest.is_file():
        return False

    rancho_logo = _ASSETS_DIR / "rancho_wordmark.png"
    chdi_logo = _ASSETS_DIR / "chdi_wordmark.png"
    green = _hex_to_rgb01(REPORT_STYLE.rule_green)
    footer_color = _hex_to_rgb01(REPORT_STYLE.green_major)
    url = REPORT_STYLE.footer_url

    def _logo_rect(path: Path, x: float, y: float, height: float) -> "fitz.Rect | None":
        if not path.is_file():
            return None
        try:
            pix = fitz.Pixmap(str(path))
            ratio = (pix.width / pix.height) if pix.height else 3.0
        except Exception:  # noqa: BLE001
            ratio = 3.0
        return fitz.Rect(x, y, x + height * ratio, y + height)

    doc = None
    try:
        doc = fitz.open(str(dest))
        for index, page in enumerate(doc):
            if index < start_page and not (stamp_cover and index == 0):
                continue
            rect = page.rect
            left = rect.x0 + 36
            right = rect.x1 - 36

            # Header: Rancho logo (left) + CHDI logo (right) + green rule.
            r_rect = _logo_rect(rancho_logo, left, rect.y0 + 18, 20)
            if r_rect is not None:
                page.insert_image(r_rect, filename=str(rancho_logo))
            c_rect = _logo_rect(chdi_logo, 0, rect.y0 + 12, 28)
            if c_rect is not None:
                c_rect = fitz.Rect(right - c_rect.width, rect.y0 + 12, right,
                                   rect.y0 + 12 + c_rect.height)
                page.insert_image(c_rect, filename=str(chdi_logo))
            rule_y = rect.y0 + 46
            page.draw_line(
                fitz.Point(left, rule_y), fitz.Point(right, rule_y),
                color=green, width=1.5,
            )

            # Footer: green rule + site URL (left) + page number (right).
            foot_y = rect.y1 - 30
            page.draw_line(
                fitz.Point(left, foot_y), fitz.Point(right, foot_y),
                color=green, width=1.0,
            )
            page.insert_text(
                fitz.Point(left, foot_y + 12), url,
                fontsize=8, color=footer_color,
            )
            page_label = str(index + 1)
            page.insert_text(
                fitz.Point(right - 4 * len(page_label) - 6, foot_y + 12),
                page_label, fontsize=8, color=footer_color,
            )

        tmp = dest.with_suffix(".stamped.pdf")
        doc.save(str(tmp))
        doc.close()
        doc = None
        tmp.replace(dest)
        return True
    except Exception as exc:  # noqa: BLE001
        logger.warning("Page stamping failed; unstamped PDF retained: %s", exc)
        return False
    finally:
        if doc is not None:
            try:
                doc.close()
            except Exception:  # noqa: BLE001
                pass


def render_rancho_pdf(
    html_document: str,
    pdf_path: str | Path,
    *,
    page_size: str = "letter",
    stamp_page_chrome: bool = True,
    stamp_cover: bool = False,
) -> Path | None:
    """Write PDF via PyMuPDF Story. Returns path on success, else ``None``.

    When ``stamp_page_chrome`` is True the rendered PDF is reopened and repeated
    Rancho page chrome (header logos, green rule, footer URL, page number) is
    stamped on body pages. ``stamp_cover`` controls whether the cover page is
    also stamped, so output can match ``SREBF2_report.pdf``.
    """
    try:
        import fitz
    except Exception as exc:  # noqa: BLE001
        logger.warning("PDF export unavailable (pymupdf not importable): %s", exc)
        return None

    dest = Path(pdf_path)
    dest.parent.mkdir(parents=True, exist_ok=True)

    try:
        mediabox = fitz.paper_rect("a4" if page_size.lower() == "a4" else "letter")
        story = fitz.Story(html=html_document)
        writer = fitz.DocumentWriter(str(dest))
        more = True
        while more:
            device = writer.begin_page(mediabox)
            where = mediabox + (36, 56, -36, -48)
            more, _ = story.place(where)
            story.draw(device)
            writer.end_page()
        writer.close()
    except Exception as exc:  # noqa: BLE001
        logger.warning("PDF export failed; HTML report remains available: %s", exc)
        if dest.exists():
            try:
                dest.unlink()
            except OSError:
                pass
        return None

    if stamp_page_chrome:
        _stamp_pdf_chrome(dest, stamp_cover=stamp_cover, start_page=1)
    return dest


def write_rancho_report(
    doc: ReportDocument,
    *,
    output_dir: str | Path | None = None,
    settings: Settings | None = None,
    include_endnotes: bool = False,
    write_pdf: bool = True,
    toc_page_numbers: dict[str, int] | None = None,
    stamp_cover: bool = False,
    show_cover_logos: bool = False,
) -> dict[str, Path]:
    """Write ``{run_id}_rancho_report.html`` (+ json, optional pdf)."""
    cfg = settings or get_settings()
    out = Path(output_dir) if output_dir is not None else cfg.output_path
    out.mkdir(parents=True, exist_ok=True)

    html_doc = render_rancho_html(
        doc,
        include_endnotes=include_endnotes,
        toc_page_numbers=toc_page_numbers,
        show_cover_logos=show_cover_logos,
    )
    html_path = out / f"{doc.dossier_run_id}_rancho_report.html"
    json_path = out / f"{doc.dossier_run_id}_rancho_report.json"
    html_path.write_text(html_doc, encoding="utf-8")
    json_path.write_text(
        json.dumps(doc.model_dump(mode="json"), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    paths: dict[str, Path] = {"html": html_path, "json": json_path}
    if write_pdf:
        pdf_path = out / f"{doc.dossier_run_id}_rancho_report.pdf"
        # Page chrome is stamped onto every PDF page, so render the PDF source
        # HTML without the single in-flow header/footer to avoid duplication.
        pdf_html = render_rancho_html(
            doc,
            include_endnotes=include_endnotes,
            include_page_chrome=False,
            toc_page_numbers=toc_page_numbers,
            show_cover_logos=show_cover_logos,
        )
        written = render_rancho_pdf(pdf_html, pdf_path, stamp_cover=stamp_cover)
        if written is not None:
            paths["pdf"] = written
    return paths


def build_and_write_rancho_report(
    *,
    dossier_run_id: str,
    gene_symbol: str,
    evidence_records: Iterable[Any],
    curator: str | None = None,
    report_date: str | None = None,
    chromosome: str | None = None,
    references: Iterable[str] | None = None,
    output_dir: str | Path | None = None,
    settings: Settings | None = None,
    include_endnotes: bool = False,
    write_pdf: bool = True,
    toc_page_numbers: dict[str, int] | None = None,
    stamp_cover: bool = False,
    show_cover_logos: bool = False,
) -> tuple[ReportDocument, dict[str, Path]]:
    """Build a ReportDocument from evidence and write the Rancho report."""
    from gene_dossier.report_schema import build_report_document

    doc = build_report_document(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        evidence_records=evidence_records,
        curator=curator,
        report_date=report_date,
        chromosome=chromosome,
        references=references,
    )
    paths = write_rancho_report(
        doc,
        output_dir=output_dir,
        settings=settings,
        include_endnotes=include_endnotes,
        write_pdf=write_pdf,
        toc_page_numbers=toc_page_numbers,
        stamp_cover=stamp_cover,
        show_cover_logos=show_cover_logos,
    )
    return doc, paths


__all__ = [
    "render_rancho_html",
    "render_rancho_pdf",
    "write_rancho_report",
    "build_and_write_rancho_report",
]
