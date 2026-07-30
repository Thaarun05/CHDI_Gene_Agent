"""Rancho BioSciences / CHDI polished dossier renderer.

Consumes a :class:`~gene_dossier.report_schema.ReportDocument` and writes a
visual report that follows ``SREBF2_report.pdf``:

- cover page (gene + CHR, Gene Report, prepared-for, curator/date)
- table of contents
- 15 major sections (green) with lettered subsections (orange)
- optional synthesized ``narrative_markdown`` as preferred prose
- supporting EvidenceRecord blocks (tables / display_text)
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
import re
from pathlib import Path
from typing import Any, Iterable

from gene_dossier.config import Settings, get_settings
from gene_dossier.models import ReportSection
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

# Page-break sentinel for section-bundle output. It is a plain HTML comment
# (invisible to browsers) that :func:`render_rancho_pdf` uses to split a document
# into consecutive ``fitz.Story`` segments, each starting on a fresh page. Only
# ``render_section_bundle_html`` emits it; the full dossier renderer neither
# emits nor depends on it.
SECTION_1C_PDF_PAGE_BREAK = "<!--RANCHO_PDF_PAGE_BREAK-->"

# Known synthesis section titles (duplicate-heading suppression only — not a slot map).
_SYNTHESIS_TITLE_ALLOWLIST: frozenset[str] = frozenset(
    {
        "general gene information",
        "gene aliases and identifiers",
        "conservation / orthologs",
        "known structure / domains",
        "alphafold / pdbe / cdd",
        "homologues",
        "tissue and cell expression",
        "geo perturbations",
        "transcription factors",
        "protein-protein interactions",
        "ctd perturbations",
        "chemical tools",
        "eqtls",
        "clinvar / omim / open targets / snps",
        "pathways",
        "knockouts / model phenotypes",
        "major labs / literature",
        "antibodies",
        "patents",
        "nih/erc grants",
    }
)


def _escape(text: str | None) -> str:
    return html.escape(text or "", quote=True)


_SOURCE_ID_TOKEN_RE = re.compile(
    r"\[source_id\s*=\s*[^\]]+\]",
    flags=re.IGNORECASE,
)


def sanitize_polished_citation_tokens(text: str | None) -> str:
    """Remove internal ``[source_id=...]`` tokens from polished visible prose."""
    if not text:
        return ""
    cleaned = _SOURCE_ID_TOKEN_RE.sub("", text)
    return re.sub(r"[ \t]{2,}", " ", cleaned)


def _evidence_attr(block: ReportContentBlock) -> str:
    """Nonvisual opaque evidence attributes for polished HTML."""
    attrs = []
    if block.evidence_ref:
        attrs.append(f' data-evidence-ref="{_escape(block.evidence_ref)}"')
        attrs.append(' data-evidence-supported="true"')
    if block.presentation_item_key:
        attrs.append(f' data-item-key="{_escape(block.presentation_item_key)}"')
    return "".join(attrs)


def _norm_title(text: str) -> str:
    return " ".join((text or "").strip().lower().split())


def _truncate_excerpt(text: str, *, limit: int = 120) -> str:
    text = (text or "").strip()
    if len(text) > limit:
        return text[: limit - 3] + "…"
    return text


def _strip_leading_synthesis_title(markdown: str) -> str:
    """Drop a leading bold line only when it matches a known synthesis title."""
    lines = markdown.splitlines()
    idx = next((i for i, line in enumerate(lines) if line.strip()), None)
    if idx is None:
        return markdown
    stripped = lines[idx].strip()
    match = re.match(r"^\*\*(.+?)\*\*(?:\s*\([^)]*\))?\s*$", stripped)
    if not match:
        return markdown
    title_norm = _norm_title(match.group(1))
    if title_norm in {"key findings", "limitations"}:
        return markdown
    if title_norm not in _SYNTHESIS_TITLE_ALLOWLIST:
        return markdown
    rest = lines[:idx] + lines[idx + 1 :]
    while idx < len(rest) and not rest[idx].strip():
        rest = rest[:idx] + rest[idx + 1 :]
    return "\n".join(rest)


def _plain_excerpt_from_markdown(markdown: str | None) -> str:
    """Plain-text excerpt for endnotes; empty when no usable narrative."""
    if not (markdown or "").strip():
        return ""
    text = _strip_leading_synthesis_title(markdown or "")
    text = re.sub(r"\*\*", "", text)
    text = re.sub(r"^#+\s*", "", text, flags=re.MULTILINE)
    text = re.sub(r"^[-*]\s+", "", text, flags=re.MULTILINE)
    text = re.sub(r"^\d+\.\s+", "", text, flags=re.MULTILINE)
    return _truncate_excerpt(" ".join(text.split()))


def _inline_format(text: str) -> str:
    """Escape text and convert simple ``**bold**`` spans."""
    parts: list[str] = []
    remaining = sanitize_polished_citation_tokens(text)
    while True:
        start = remaining.find("**")
        if start < 0:
            parts.append(_escape(remaining))
            break
        parts.append(_escape(remaining[:start]))
        end = remaining.find("**", start + 2)
        if end < 0:
            parts.append(_escape(remaining[start:]))
            break
        parts.append(f"<strong>{_escape(remaining[start + 2 : end])}</strong>")
        remaining = remaining[end + 2 :]
    return "".join(parts)


def _render_narrative_markdown(
    markdown_text: str | None,
    *,
    synthesis_status: str | None = None,
) -> str:
    """Convert limited synthesis markdown into escaped HTML (no Markdown library)."""
    if not (markdown_text or "").strip():
        return ""
    text = sanitize_polished_citation_tokens(
        _strip_leading_synthesis_title(markdown_text or "")
    )
    if not text.strip():
        return ""

    body: list[str] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if not stripped:
            i += 1
            continue

        if stripped.startswith("### "):
            body.append(
                f'<h4 class="narrative-subheading">'
                f"{_inline_format(stripped[4:].strip())}</h4>"
            )
            i += 1
            continue
        if stripped.startswith("## "):
            body.append(
                f'<h4 class="narrative-subheading">'
                f"{_inline_format(stripped[3:].strip())}</h4>"
            )
            i += 1
            continue

        bold_only = re.match(r"^\*\*(.+?)\*\*\s*$", stripped)
        if bold_only:
            inner = bold_only.group(1).strip()
            inner_norm = _norm_title(inner)
            if inner_norm == "key findings":
                body.append(
                    f'<h4 class="narrative-subheading key-findings">'
                    f"{_escape(inner)}</h4>"
                )
                i += 1
                continue
            if inner_norm == "limitations":
                body.append(
                    f'<h4 class="narrative-subheading limitations">'
                    f"{_escape(inner)}</h4>"
                )
                i += 1
                continue
            # Other standalone bold lines stay as emphasis, not headings.
            body.append(f"<p><strong>{_escape(inner)}</strong></p>")
            i += 1
            continue

        if stripped.startswith("- "):
            items: list[str] = []
            while i < len(lines) and lines[i].strip().startswith("- "):
                items.append(f"<li>{_inline_format(lines[i].strip()[2:])}</li>")
                i += 1
            body.append("<ul>" + "".join(items) + "</ul>")
            continue

        if re.match(r"^\d+\.\s", stripped):
            items = []
            while i < len(lines) and re.match(r"^\d+\.\s", lines[i].strip()):
                item_text = re.sub(r"^\d+\.\s+", "", lines[i].strip())
                items.append(f"<li>{_inline_format(item_text)}</li>")
                i += 1
            body.append("<ol>" + "".join(items) + "</ol>")
            continue

        para_lines: list[str] = []
        while i < len(lines):
            candidate = lines[i].strip()
            if not candidate:
                break
            if candidate.startswith("##") or candidate.startswith("- "):
                break
            if re.match(r"^\d+\.\s", candidate):
                break
            if re.match(r"^\*\*(.+?)\*\*\s*$", candidate):
                break
            para_lines.append(candidate)
            i += 1
        if para_lines:
            body.append(f"<p>{_inline_format(' '.join(para_lines))}</p>")
        else:
            body.append(f"<p>{_escape(stripped)}</p>")
            i += 1

    if not body:
        return ""

    attrs = ' class="synthesized-narrative"'
    if synthesis_status:
        attrs += f' data-synthesis-status="{_escape(synthesis_status)}"'
    return f"<div{attrs}>\n" + "\n".join(body) + "\n</div>"


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
  break-after: avoid;
  page-break-after: avoid;
}}
.section-bundle-body .report-subsection + .subsection-c,
.section-bundle-body .report-subsection + .subsection-1c {{
  break-before: page;
  page-break-before: always;
}}
.block {{ margin: 0 0 10pt 0; }}
.empty-note {{
  font-style: italic;
  color: #7a5a45;
}}
.synthesized-narrative {{
  color: {s.brown_body};
  margin: 0 0 12pt 0;
}}
.synthesized-narrative p {{
  margin: 0 0 8pt 0;
}}
.synthesized-narrative ul,
.synthesized-narrative ol {{
  margin: 0 0 10pt 0;
  padding-left: 18pt;
}}
.narrative-subheading {{
  color: {s.orange_sub};
  font-size: {s.subsection_header_pt}pt;
  font-weight: 600;
  margin: 12pt 0 6pt 0;
}}
.key-findings, .limitations {{
  color: {s.orange_sub};
}}
.supporting-evidence {{
  margin: 12pt 0 8pt 0;
  padding-top: 8pt;
  border-top: 1pt solid {s.rule_green};
}}
.supporting-evidence-heading {{
  color: {s.orange_sub};
  font-size: 12pt;
  font-weight: 600;
  margin: 0 0 8pt 0;
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
  border: 1px solid #6e5a4a;
}}
table.rancho tr {{
  break-inside: avoid;
  page-break-inside: avoid;
}}
table.rancho tr.header-row th {{
  background-color: {s.table_header_bg};
}}
table.rancho td {{
  padding: 5pt 8pt;
  border: 1px solid #6e5a4a;
  vertical-align: top;
  text-align: left;
  word-wrap: break-word;
  overflow-wrap: break-word;
}}
table.gene-aliases-table {{
  width: 100%;
  table-layout: fixed;
}}
table.gene-aliases-table col.col-label {{ width: 12%; }}
table.gene-aliases-table th:first-child,
table.gene-aliases-table td:first-child {{
  width: 12%;
  font-weight: 600;
  white-space: normal;
  overflow-wrap: break-word;
  word-break: normal;
  hyphens: manual;
}}
table.gene-aliases-table col.col-human {{ width: 29%; }}
table.gene-aliases-table col.col-mouse {{ width: 32%; }}
table.gene-aliases-table col.col-rat {{ width: 27%; }}
table.section-1c-pdb-table {{
  table-layout: fixed;
  font-size: 10pt;
}}
table.section-1c-pdb-table col.col-pdb {{ width: 10%; }}
table.section-1c-pdb-table col.col-chains {{ width: 14%; }}
table.section-1c-pdb-table col.col-method {{ width: 17%; }}
table.section-1c-pdb-table col.col-resolution {{ width: 13%; }}
table.section-1c-pdb-table col.col-span {{ width: 14%; }}
table.section-1c-pdb-table col.col-coverage {{ width: 18%; }}
table.section-1c-pdb-table col.col-selection {{ width: 14%; }}
table.section-1c-pdb-table th:first-child,
table.section-1c-pdb-table td:first-child {{
  white-space: nowrap;
  overflow-wrap: normal;
  word-break: keep-all;
}}
a.id-link {{
  color: {s.orange_link};
  text-decoration: underline;
}}
a.id-link-uniprot {{
  color: {s.brown_body};
  text-decoration: none;
}}
a.id-link-uniprot:hover,
a.id-link-uniprot:focus {{
  text-decoration: underline;
}}
a.id-link:focus {{
  outline: 2px solid {s.orange_link};
  outline-offset: 1px;
}}
figure.rancho-figure {{
  margin: 10pt 0 14pt 0;
  text-align: center;
}}
figure.rancho-figure img {{
  max-width: 100%;
  height: auto;
}}
figure.rancho-figure.ucsc-conservation-figure {{
  break-inside: avoid;
  page-break-inside: avoid;
  margin: 6pt 0 12pt 0;
}}
figure.rancho-figure.ucsc-conservation-figure img {{
  width: 100%;
  max-width: 100%;
  height: auto;
  display: block;
}}
.section-bundle-body figure.rancho-figure.section-1c-domain-architecture-figure {{
  text-align: left;
  margin: 6pt 0 10pt 0;
  break-inside: avoid;
  page-break-inside: avoid;
}}
.section-bundle-body figure.rancho-figure.section-1c-domain-architecture-figure img {{
  width: 100%;
  max-width: 7.3in;
  max-height: 0.95in;
  object-fit: contain;
  display: block;
}}
.section-bundle-body figure.rancho-figure.section-1c-pdb-domain-focus-figure {{
  text-align: left;
  margin: 8pt 0 12pt 0;
  break-inside: avoid;
  page-break-inside: avoid;
}}
.section-bundle-body figure.rancho-figure.section-1c-pdb-domain-focus-figure img {{
  width: 2.1in;
  max-width: 50%;
  display: block;
}}
.section-bundle-body figure.rancho-figure.section-1c-pdb-assembly-figure {{
  text-align: left;
  margin: 10pt 0 12pt 0;
  break-inside: avoid;
  page-break-inside: avoid;
}}
.section-bundle-body figure.rancho-figure.section-1c-pdb-assembly-figure img {{
  width: 3.2in;
  max-width: 70%;
  display: block;
}}
a.ucsc-transcript-link {{
  color: {s.orange_link};
  text-decoration: underline;
}}
.ucsc-transcript-line {{
  margin: 8pt 0 6pt 0;
  break-after: avoid;
  page-break-after: avoid;
}}
.section-1c-link {{
  color: {s.orange_link};
  text-decoration: underline;
}}
.section-1c-link-line {{
  margin: 7pt 0 7pt 0;
  break-after: avoid;
  page-break-after: avoid;
}}
.back-matter h2, .endnotes h2 {{
  color: {s.green_major};
  font-size: 18pt;
}}
.db-list, .ref-list, .endnote-list {{ padding-left: 18pt; }}
.provenance-muted {{ color: #8a7a6a; font-size: 8.5pt; }}
.section-preview-body, .report-page {{
  max-width: 8.5in;
  margin: 24pt auto;
  padding: 0 24pt 24pt 24pt;
  background: #ffffff;
}}
.report-page.report-chrome {{
  padding-bottom: 0;
  padding-top: 0;
}}
.section-bundle-body .block {{
  margin: 6pt 0 8pt 0;
}}
.section-bundle-body .narrative p {{
  background-color: #ffffff;
}}
.section-bundle-body h3.sub-heading {{
  margin: 10pt 0 6pt 0;
}}
.section-bundle-body h2.major-heading {{
  margin: 0 0 8pt 0;
}}
.section-bundle-body .subsection-c,
.section-bundle-body .subsection-1c {{
  max-width: 6.9in;
  font-size: 9.5pt;
  line-height: 1.18;
  overflow-wrap: anywhere;
  word-break: normal;
}}
.section-bundle-body .subsection-c .block,
.section-bundle-body .subsection-1c .block {{
  margin: 4pt 0 6pt 0;
}}
.section-bundle-body .section-1c-domain-group,
.section-bundle-body .section-1c-feature-group,
.section-bundle-body .section-1c-pdb-group {{
  break-inside: avoid;
  page-break-inside: avoid;
  margin: 4pt 0 8pt 0;
}}
.section-bundle-body .section-1c-pdb-group {{
  display: block;
  break-before: page;
  page-break-before: always;
  margin-top: 0;
}}
@media screen {{
  .section-bundle-body .section-1c-pdb-group {{
    margin-top: 72pt;
    padding-top: 24pt;
    border-top: 1px solid transparent;
  }}
}}
/* Section-bundle continuation pages already start a new page/Story, so the
   CSS break would otherwise add an empty leading page. */
.section-bundle-body.section-1c-continuation .section-1c-pdb-group {{
  break-before: auto;
  page-break-before: auto;
  margin-top: 0;
  padding-top: 0;
}}
@media screen {{
  .section-bundle-body.section-1c-continuation .section-1c-pdb-group {{
    margin-top: 0;
    padding-top: 0;
  }}
}}
.section-bundle-body .section-1c-domain-group p,
.section-bundle-body .section-1c-feature-group p {{
  margin-bottom: 6pt;
}}
.section-bundle-body figure.rancho-figure.section-1c-domain-thumbnail,
.section-bundle-body figure.rancho-figure.section-1c-feature-thumbnail {{
  text-align: left;
  margin: 5pt 0 8pt 0;
}}
/* Absolute widths: the PyMuPDF Story engine honours ``width`` but ignores
   ``max-width``/``max-height``, so thumbnails need an explicit size to match the
   original report scale in PDF and PNG as well as in the browser. */
.section-bundle-body figure.rancho-figure.section-1c-domain-thumbnail img {{
  width: 1.45in;
  max-width: 1.55in;
  max-height: 1.75in;
  height: auto;
  object-fit: contain;
}}
.section-bundle-body figure.rancho-figure.section-1c-feature-thumbnail img {{
  width: 2.15in;
  max-width: 2.25in;
  max-height: 2.45in;
  height: auto;
  object-fit: contain;
}}
.section-bundle-body figure.rancho-figure.section-1c-pdb-official-image {{
  text-align: left;
  margin: 8pt 0 4pt 0;
}}
.section-bundle-body figure.rancho-figure.section-1c-pdb-official-image img {{
  width: 3.4in;
  max-width: 55%;
  max-height: 4.8in;
  height: auto;
  object-fit: contain;
}}
.section-bundle-body figure.rancho-figure.section-1c-pdb-official-image figcaption,
.section-bundle-body .section-1c-image-attribution {{
  color: #8a7a6a;
  font-size: 8.5pt;
  text-align: left;
}}
.section-bundle-body .subsection-d {{
  margin-top: 10pt;
}}
.section-bundle-body .section-1d-link-line {{
  margin: 4pt 0 6pt 0;
  font-size: 13pt;
}}
.section-bundle-body .section-1d-link {{
  text-decoration: underline;
}}
.section-bundle-body .section-1d-status-line {{
  margin: 4pt 0 6pt 0;
  font-size: 12.5pt;
  color: #8a7a6a;
}}
.section-bundle-body .section-1d-visual-table {{
  width: 100%;
  border-collapse: collapse;
  border: none;
  background: transparent;
  margin: 2pt 0 14pt 0;
  break-inside: avoid;
  page-break-inside: avoid;
}}
.section-bundle-body .section-1d-visual-table td {{
  border: none;
  background: transparent;
  vertical-align: top;
  padding: 0 6pt 0 0;
}}
.section-bundle-body .section-1d-visual-table col.section-1d-model-col {{
  width: 47%;
}}
.section-bundle-body .section-1d-visual-table col.section-1d-confidence-col {{
  width: 24%;
}}
.section-bundle-body .section-1d-visual-table col.section-1d-blurb-col {{
  width: 29%;
}}
.section-bundle-body figure.rancho-figure.section-1d-human-structure-capture {{
  margin: 0;
  max-width: 100%;
}}
.section-bundle-body figure.rancho-figure.section-1d-human-structure-capture img {{
  width: 100%;
  max-width: 4.0in;
  max-height: 3.25in;
  height: auto;
  object-fit: contain;
  display: block;
}}
.section-bundle-body .section-1d-confidence-legend {{
  font-size: 10pt;
  line-height: 1.3;
}}
.section-bundle-body .section-1d-confidence-legend .legend-title {{
  font-weight: 700;
  font-size: 10.5pt;
  margin-bottom: 4pt;
}}
.section-bundle-body .section-1d-confidence-legend .legend-row {{
  margin: 1.5pt 0;
  font-size: 10pt;
}}
.section-bundle-body .section-1d-confidence-legend .swatch {{
  display: inline-block;
  width: 0.55em;
  height: 0.55em;
  margin-right: 0.35em;
  vertical-align: middle;
}}
.section-bundle-body .section-1d-confidence-legend .swatch-very-high {{ background: #0053D6; }}
.section-bundle-body .section-1d-confidence-legend .swatch-high {{ background: #65CBF3; }}
.section-bundle-body .section-1d-confidence-legend .swatch-low {{ background: #FFDB13; }}
.section-bundle-body .section-1d-confidence-legend .swatch-very-low {{ background: #FF7D45; }}
.section-bundle-body .section-1d-blurb {{
  font-size: 9.5pt;
  line-height: 1.35;
  color: {REPORT_STYLE.brown_body};
}}
@media (max-width: 900px) {{
  .section-preview-body, .report-page {{
    max-width: 100%;
    margin: 12pt auto;
    padding: 0 12pt 16pt 12pt;
  }}
}}
@page {{
  size: Letter;
  margin: 0.5in;
}}
@media print {{
  .report-page {{
    max-width: none;
    margin: 0;
    padding-left: 0;
    padding-right: 0;
  }}
  .page-header {{
    margin-bottom: 12pt;
  }}
  .page-footer {{
    margin-top: 12pt;
  }}
}}
""".strip()


def _img_tag(data_uri: str | None, *, cls: str, alt: str) -> str:
    if not data_uri:
        return f'<span class="{_escape(cls)}-fallback">{_escape(alt)}</span>'
    return f'<img class="{_escape(cls)}" src="{data_uri}" alt="{_escape(alt)}" />'


def _section_1c_emphasized_paragraph(text: str) -> str:
    """Escape a Section 1c paragraph, bolding a known specific-domain lead phrase."""
    from gene_dossier.report_presentation import SECTION_1C_BOLD_LEAD_PHRASES

    for phrase in SECTION_1C_BOLD_LEAD_PHRASES:
        if text.startswith(phrase):
            return (
                f"<strong>{_escape(phrase)}</strong>"
                f"{_escape(text[len(phrase):])}"
            )
    return _escape(text)


def _render_block(block: ReportContentBlock) -> str:
    from gene_dossier.report_presentation import format_safe_table_cell_html

    parts: list[str] = [f'<div class="block"{_evidence_attr(block)}>']
    inline_section_1c_link_rendered = False
    if block.title:
        parts.append(
            f'<p><strong style="color:{REPORT_STYLE.orange_link};">'
            f"{_escape(block.title)}</strong></p>"
        )

    if block.kind == "table" and block.table_headers:
        classes = ["rancho"]
        if block.presentation_role == "gene_aliases_table":
            classes.append("gene-aliases-table")
        elif block.presentation_role == "section_1c_pdb_table":
            classes.append("section-1c-pdb-table")
        class_attr = " ".join(classes)
        parts.append(f'<table class="{class_attr}">')
        if block.presentation_role == "gene_aliases_table":
            parts.append(
                "<colgroup>"
                '<col class="col-label" />'
                '<col class="col-human" />'
                '<col class="col-mouse" />'
                '<col class="col-rat" />'
                "</colgroup>"
            )
        elif block.presentation_role == "section_1c_pdb_table":
            parts.append(
                "<colgroup>"
                '<col class="col-pdb" />'
                '<col class="col-chains" />'
                '<col class="col-method" />'
                '<col class="col-resolution" />'
                '<col class="col-span" />'
                '<col class="col-coverage" />'
                '<col class="col-selection" />'
                "</colgroup>"
            )
        parts.append('<tbody><tr class="header-row">')
        for h in block.table_headers:
            parts.append(
                f'<th style="background-color:{REPORT_STYLE.table_header_bg};">'
                f"{_escape(h)}</th>"
            )
        parts.append("</tr>")
        for row in block.table_rows:
            parts.append("<tr>")
            for cell in row:
                parts.append(f"<td>{format_safe_table_cell_html(cell)}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table>")
    elif block.kind == "figure" and block.figure_path:
        fig_uri = None
        fig_path = Path(block.figure_path)
        if not fig_path.is_file():
            try:
                from gene_dossier.ucsc_figure import resolve_artifact_path

                fig_path = resolve_artifact_path(block.figure_path)
            except Exception:  # noqa: BLE001
                fig_path = Path(block.figure_path)
        if fig_path.is_file():
            mime, _ = mimetypes.guess_type(fig_path.name)
            mime = mime or "image/png"
            fig_uri = (
                f"data:{mime};base64,"
                + base64.b64encode(fig_path.read_bytes()).decode("ascii")
            )
        caption = _escape(block.figure_caption or "")
        alt = caption or "figure"
        if fig_uri:
            fig_classes = ["rancho-figure"]
            if block.presentation_role:
                fig_classes.append(
                    re.sub(
                        r"[^a-z0-9]+",
                        "-",
                        str(block.presentation_role).lower(),
                    ).strip("-")
                )
            if block.presentation_role == "ucsc_conservation_figure":
                fig_classes.append("ucsc-conservation-figure")
            parts.append(f'<figure class="{" ".join(fig_classes)}">')
            parts.append(f'<img src="{fig_uri}" alt="{alt}" />')
            if (
                caption
                and block.presentation_role != "ucsc_conservation_figure"
                and not str(block.presentation_role or "").startswith("section_1d_")
                and (
                    not str(block.presentation_role or "").startswith("section_1c_")
                    or block.presentation_role == "section_1c_pdb_official_image"
                )
            ):
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
        # Single UCSC transcript line: render as one underlined orange link.
        if len(block.links) == 1:
            link = block.links[0]
            label = _escape(link.get("label") or block.text or link.get("url") or "link")
            url = _escape(link.get("url") or "#")
            role = str(block.presentation_role or "")
            if role.startswith("section_1d_"):
                prefix = _escape(block.text or "")
                parts.append(
                    f'<p class="section-1d-link-line">'
                    f"{prefix}"
                    f'<a class="section-1d-link" style="color:{REPORT_STYLE.orange_link};" '
                    f'href="{url}">{label}</a></p>'
                )
            else:
                line_class = "ucsc-transcript-line"
                link_class = "ucsc-transcript-link"
                if role.startswith("section_1c_"):
                    line_class = "section-1c-link-line"
                    link_class = "section-1c-link"
                parts.append(
                    f'<p class="{line_class}">'
                    f'<a class="{link_class}" style="color:{REPORT_STYLE.orange_link};" '
                    f'href="{url}">{label}</a></p>'
                )
        else:
            parts.append("<ul>")
            for link in block.links:
                label = _escape(link.get("label") or link.get("url") or "link")
                url = _escape(link.get("url") or "#")
                parts.append(f'<li><a href="{url}">{label}</a></li>')
            parts.append("</ul>")
    else:
        text = sanitize_polished_citation_tokens((block.text or "").strip())
        if text:
            parts.append('<div class="narrative">')
            if (
                block.presentation_role == "section_1c_domain_summary"
                and len(block.links) == 1
                and block.links[0].get("label")
                and text.startswith(str(block.links[0].get("label")))
            ):
                link = block.links[0]
                label_raw = str(link.get("label") or "")
                label = _escape(label_raw)
                url = _escape(link.get("url") or "#")
                parts.append(
                    f'<p class="section-1c-link-line">'
                    f'<a class="section-1c-link" '
                    f'style="color:{REPORT_STYLE.orange_link};" '
                    f'href="{url}">{label}</a></p>'
                )
                rest = text[len(label_raw) :].lstrip()
                for para in rest.split("\n\n"):
                    para = para.strip()
                    if para:
                        parts.append(f"<p>{_escape(para)}</p>")
                inline_section_1c_link_rendered = True
            else:
                paras = [para.strip() for para in text.split("\n\n") if para.strip()]
                # Rancho body style closes a domain paragraph with its CDD link.
                trailing_link = (
                    block.links[0]
                    if (
                        block.presentation_role == "section_1c_domain_summary"
                        and len(block.links) == 1
                        and paras
                    )
                    else None
                )
                for index, para in enumerate(paras):
                    body = (
                        _section_1c_emphasized_paragraph(para)
                        if block.presentation_role == "section_1c_domain_summary"
                        else _escape(para)
                    )
                    if trailing_link is not None and index == len(paras) - 1:
                        label = _escape(trailing_link.get("label") or "Link")
                        url = _escape(trailing_link.get("url") or "#")
                        body += (
                            f' <a class="section-1c-link" '
                            f'style="color:{REPORT_STYLE.orange_link};" '
                            f'href="{url}">{label}</a>'
                        )
                    parts.append(f"<p>{body}</p>")
                if trailing_link is not None:
                    inline_section_1c_link_rendered = True
            parts.append("</div>")

    if block.links and block.kind not in {"link"} and not inline_section_1c_link_rendered:
        for link in block.links:
            label = _escape(link.get("label") or "Link")
            url = _escape(link.get("url") or "#")
            if str(block.presentation_role or "").startswith("section_1c_"):
                parts.append(
                    f'<p class="section-1c-link-line">'
                    f'<a class="section-1c-link" '
                    f'style="color:{REPORT_STYLE.orange_link};" '
                    f'href="{url}">{label}</a></p>'
                )
            else:
                parts.append(f'<p><a href="{url}">{label}</a></p>')

    parts.append("</div>")
    return "\n".join(parts)


def _render_supporting_evidence(blocks: list[ReportContentBlock]) -> str:
    parts = [
        '<div class="supporting-evidence">',
        (
            f'<h4 class="supporting-evidence-heading" '
            f'style="color:{REPORT_STYLE.orange_sub};">Supporting evidence</h4>'
        ),
    ]
    for block in blocks:
        parts.append(_render_block(block))
    parts.append("</div>")
    return "\n".join(parts)


def _section_1c_group_class(block: ReportContentBlock) -> str | None:
    role = str(block.presentation_role or "")
    if role in {
        "section_1c_domain_summary",
        "section_1c_domain_thumbnail",
    }:
        return "section-1c-domain-group"
    if role in {
        "section_1c_feature_summary",
        "section_1c_feature_thumbnail",
    }:
        return "section-1c-feature-group"
    if role in {
        "section_1c_pdb_link",
        "section_1c_pdb_official_image",
        "section_1c_image_attribution",
    }:
        return "section-1c-pdb-group"
    return None


def _render_section_1c_grouped_blocks(blocks: list[ReportContentBlock]) -> str:
    parts: list[str] = []
    current_class: str | None = None
    current_key: str | None = None
    current_blocks: list[ReportContentBlock] = []

    def flush() -> None:
        nonlocal current_class, current_key, current_blocks
        if not current_blocks:
            return
        if current_class and current_key:
            parts.append(
                f'<div class="{_escape(current_class)}" data-item-key="{_escape(current_key)}">'
            )
            parts.extend(_render_block(block) for block in current_blocks)
            parts.append("</div>")
        else:
            parts.extend(_render_block(block) for block in current_blocks)
        current_class = None
        current_key = None
        current_blocks = []

    for block in blocks:
        group_class = _section_1c_group_class(block)
        group_key = block.presentation_item_key if group_class else None
        if group_class and group_key:
            if current_blocks and (group_class != current_class or group_key != current_key):
                flush()
            current_class = group_class
            current_key = group_key
            current_blocks.append(block)
        else:
            flush()
            parts.append(_render_block(block))
    flush()
    return "\n".join(parts)


def split_section_1c_page_segments(
    blocks: list[ReportContentBlock],
) -> list[list[ReportContentBlock]]:
    """Group Section 1c blocks into one list per rendered page.

    A new segment starts at every block flagged ``presentation_page_break_before``
    by the presentation builder, so page breaks are decided by content roles
    rather than by CSS pagination.
    """
    segments: list[list[ReportContentBlock]] = []
    current: list[ReportContentBlock] = []
    for block in blocks:
        if block.presentation_page_break_before and current:
            segments.append(current)
            current = []
        current.append(block)
    if current:
        segments.append(current)
    return segments


def render_section_1c_subsection_segments(sub: ReportSubsection) -> list[str]:
    """Render Section 1c as one subsection HTML string per page segment.

    Only the first segment carries the subsection heading; the section-bundle
    renderer wraps later segments in their own ``report-page`` so the C-terminal
    family block and the PDB group each start on a real new page.
    """
    segments = split_section_1c_page_segments(list(sub.presentation_blocks or []))
    heading = f"{sub.key}. {sub.title}"
    subsection_class = re.sub(r"[^a-z0-9]+", "-", sub.key.lower()).strip("-")
    rendered: list[str] = []
    for index, segment in enumerate(segments):
        parts = [
            f'<section class="report-subsection subsection-{_escape(subsection_class)}">'
        ]
        if index == 0:
            parts.append(
                f'<h3 class="sub-heading" style="color:{REPORT_STYLE.orange_sub};">'
                f"{_escape(heading)}</h3>"
            )
        parts.append(_render_section_1c_grouped_blocks(segment))
        parts.append("</section>")
        rendered.append("\n".join(parts))
    return rendered


def _render_section_1d_confidence_legend(block: ReportContentBlock) -> str:
    """Deterministic report-side pLDDT color key (title + four rows)."""
    _ = block
    return (
        '<div class="section-1d-confidence-legend">'
        '<div class="legend-title">Model Confidence</div>'
        '<div class="legend-row">'
        '<span class="swatch swatch-very-high"></span>'
        "Very high (pLDDT &gt; 90)</div>"
        '<div class="legend-row">'
        '<span class="swatch swatch-high"></span>'
        "High (90 &gt; pLDDT &gt; 70)</div>"
        '<div class="legend-row">'
        '<span class="swatch swatch-low"></span>'
        "Low (70 &gt; pLDDT &gt; 50)</div>"
        '<div class="legend-row">'
        '<span class="swatch swatch-very-low"></span>'
        "Very low (pLDDT &lt; 50)</div>"
        "</div>"
    )


def _render_section_1d_confidence_blurb() -> str:
    return (
        '<div class="section-1d-blurb">'
        "AlphaFold produces a per-residue model confidence score (pLDDT) "
        "between 0 and 100. Some regions below 50 pLDDT may be unstructured "
        "in isolation."
        "</div>"
    )


def _render_section_1d_blocks(blocks: list[ReportContentBlock]) -> str:
    """Render 1d blocks; pair capture + legend into one PDF-stable visual table."""
    parts: list[str] = []
    index = 0
    while index < len(blocks):
        block = blocks[index]
        role = str(block.presentation_role or "")
        if role == "section_1d_human_structure_capture":
            legend = None
            if index + 1 < len(blocks) and (
                blocks[index + 1].presentation_role == "section_1d_confidence_legend"
            ):
                legend = blocks[index + 1]
                index += 1
            if legend is None:
                # Capture without legend is incomplete; still emit the image alone.
                parts.append(_render_block(block))
            else:
                parts.append(
                    '<table class="section-1d-visual-table">'
                    "<colgroup>"
                    '<col class="section-1d-model-col" />'
                    '<col class="section-1d-confidence-col" />'
                    '<col class="section-1d-blurb-col" />'
                    "</colgroup>"
                    "<tr>"
                    f'<td class="section-1d-model-cell">{_render_block(block)}</td>'
                    f'<td class="section-1d-confidence-cell">'
                    f"{_render_section_1d_confidence_legend(legend)}</td>"
                    f'<td class="section-1d-blurb-cell">'
                    f"{_render_section_1d_confidence_blurb()}</td>"
                    "</tr></table>"
                )
        elif role == "section_1d_confidence_legend":
            # Legend-only visuals are not accepted; skip standalone emission.
            pass
        elif role == "section_1d_species_status":
            text = _escape((block.text or "").strip())
            parts.append(
                f'<p class="section-1d-status-line" {_evidence_attr(block)}>'
                f"{text}</p>"
            )
        else:
            parts.append(_render_block(block))
        index += 1
    return "\n".join(parts)


def _render_subsection(sub: ReportSubsection) -> str:
    heading = f"{sub.key}. {sub.title}"
    subsection_class = re.sub(r"[^a-z0-9]+", "-", sub.key.lower()).strip("-")
    parts = [
        f'<section class="report-subsection subsection-{_escape(subsection_class)}">',
        (
            f'<h3 class="sub-heading" style="color:{REPORT_STYLE.orange_sub};">'
            f"{_escape(heading)}</h3>"
        ),
    ]
    # Polished presentation_blocks are the complete human-facing subsection body.
    if sub.presentation_blocks:
        if sub.key == "c":
            parts.append(_render_section_1c_grouped_blocks(list(sub.presentation_blocks)))
        elif sub.key == "d":
            parts.append(_render_section_1d_blocks(list(sub.presentation_blocks)))
        else:
            for block in sub.presentation_blocks:
                parts.append(_render_block(block))
        parts.append("</section>")
        return "\n".join(parts)

    narrative_html = _render_narrative_markdown(
        sub.narrative_markdown,
        synthesis_status=sub.synthesis_status,
    )
    if narrative_html:
        parts.append(narrative_html)
    if sub.blocks:
        if narrative_html:
            parts.append(_render_supporting_evidence(list(sub.blocks)))
        else:
            for block in sub.blocks:
                parts.append(_render_block(block))
    elif not narrative_html:
        parts.append(
            '<p class="empty-note">No evidence available for this subsection in the '
            "current dossier run.</p>"
        )
    parts.append("</section>")
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
    narrative_html = _render_narrative_markdown(
        section.narrative_markdown,
        synthesis_status=section.synthesis_status,
    )
    if narrative_html:
        parts.append(narrative_html)
    if section.blocks:
        if narrative_html:
            parts.append(_render_supporting_evidence(list(section.blocks)))
        else:
            for block in section.blocks:
                parts.append(_render_block(block))
    if section.subsections:
        for sub in section.subsections:
            parts.append(_render_subsection(sub))
    elif not narrative_html and not section.blocks:
        parts.append(
            '<p class="empty-note">No evidence available for this section in the '
            "current dossier run.</p>"
        )
    parts.append("</section>")
    return "\n".join(parts)


def _block_excerpt(block: ReportContentBlock) -> str:
    return _truncate_excerpt((block.text or block.title or "").strip())


def _neutral_audit_excerpt(*, kind: str, title: str) -> str:
    if kind == "major":
        label = f"Cited by synthesized section: {title}"
    else:
        label = f"Cited by synthesized subsection: {title}"
    return _truncate_excerpt(label)


def _append_slot_endnotes(
    notes: list[dict[str, str]],
    seen: set[str],
    *,
    source_ids: list[str],
    blocks: list[ReportContentBlock],
    narrative_markdown: str | None,
    audit_kind: str,
    audit_title: str,
) -> None:
    evidence_excerpts: dict[str, str] = {}
    for block in blocks:
        excerpt = _block_excerpt(block)
        for sid in block.source_ids:
            if sid and sid not in evidence_excerpts and excerpt:
                evidence_excerpts[sid] = excerpt

    narrative_excerpt = _plain_excerpt_from_markdown(narrative_markdown)
    neutral = _neutral_audit_excerpt(kind=audit_kind, title=audit_title)

    for sid in source_ids:
        if not sid or sid in seen:
            continue
        seen.add(sid)
        if sid in evidence_excerpts:
            excerpt = evidence_excerpts[sid]
        elif narrative_excerpt:
            excerpt = narrative_excerpt
        else:
            excerpt = neutral
        notes.append({"source_id": sid, "excerpt": excerpt})

    for block in blocks:
        excerpt = _block_excerpt(block) or neutral
        for sid in block.source_ids:
            if not sid or sid in seen:
                continue
            seen.add(sid)
            notes.append({"source_id": sid, "excerpt": excerpt})


def _major_owned_source_ids(section: ReportMajorSection) -> list[str]:
    """Return major-level source_ids not exclusively inherited from subsections.

    ``ReportMajorSection.source_ids`` aggregates subsection IDs. Endnotes must
    process subsection-owned IDs under their subsection so precise evidence
    excerpts are preferred over major narrative / neutral audit text.
    """
    subsection_source_ids = {
        sid
        for sub in section.subsections
        for sid in (sub.source_ids or [])
        if sid
    }
    major_block_source_ids = {
        sid
        for block in section.blocks
        for sid in (block.source_ids or [])
        if sid
    }
    return [
        sid
        for sid in (section.source_ids or [])
        if sid
        and (sid not in subsection_source_ids or sid in major_block_source_ids)
    ]


def _collect_endnotes(doc: ReportDocument) -> list[dict[str, str]]:
    """Build provenance endnotes from evidence blocks and synthesis source_ids."""
    notes: list[dict[str, str]] = []
    seen: set[str] = set()
    for section in doc.sections:
        _append_slot_endnotes(
            notes,
            seen,
            source_ids=_major_owned_source_ids(section),
            blocks=list(section.blocks or []),
            narrative_markdown=section.narrative_markdown,
            audit_kind="major",
            audit_title=section.title,
        )
        for sub in section.subsections:
            _append_slot_endnotes(
                notes,
                seen,
                source_ids=list(sub.source_ids or []),
                blocks=list(sub.blocks or []),
                narrative_markdown=sub.narrative_markdown,
                audit_kind="subsection",
                audit_title=sub.title,
            )
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


def render_rancho_section_fragment(
    *,
    document: ReportDocument,
    section_number: int,
    subsection_key: str,
    show_cover_logos: bool = False,
) -> str:
    """Render one major section heading + one subsection with production helpers.

    Uses the same CSS/page chrome and block/subsection renderers as the full
    report. Does not render unrelated sections or hide them with CSS.
    """
    major = next(
        (sec for sec in document.sections if sec.number == section_number),
        None,
    )
    if major is None:
        raise ValueError(f"Unknown section number: {section_number}")
    sub = next((s for s in major.subsections if s.key == subsection_key), None)
    if sub is None:
        raise ValueError(
            f"Unknown subsection {subsection_key!r} under section {section_number}"
        )

    heading = f"{major.number}. {major.title}"
    body_parts = [
        f'<section id="section-{major.number}" class="section-preview-body">',
        (
            f'<h2 class="major-heading" style="color:{REPORT_STYLE.green_major};">'
            f"{_escape(heading)}</h2>"
        ),
        _render_subsection(sub),
        "</section>",
    ]
    body = "\n".join(body_parts)

    rancho = _asset_data_uri("rancho_wordmark.png")
    chdi = _asset_data_uri("chdi_wordmark.png")
    rancho_header = _asset_data_uri("rancho_header_bar.png") or rancho
    rancho_footer = _asset_data_uri("rancho_footer.png")
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
    # Keep logos optional for light previews.
    if not show_cover_logos:
        header = '<div class="page-header"></div>'
        footer = (
            f'<div class="page-footer"><span>{_escape(REPORT_STYLE.footer_url)}</span></div>'
        )

    title = _escape(
        f"{document.cover.gene_line} — "
        f"Section {section_number}{subsection_key} preview"
    )
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
{header}
{body}
{footer}
</body>
</html>
"""


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


def _split_pdf_page_segments(html_document: str) -> list[str]:
    """Split a document at page-break sentinels into standalone documents.

    Each segment keeps the original ``<head>``, and therefore the report
    stylesheet, so every ``fitz.Story`` renders with identical styling.
    Documents without a sentinel are returned unchanged as a single segment.
    """
    if SECTION_1C_PDF_PAGE_BREAK not in html_document:
        return [html_document]
    body_open = html_document.find("<body>")
    body_close = html_document.rfind("</body>")
    if body_open == -1 or body_close == -1 or body_close < body_open:
        return [html_document]
    body_start = body_open + len("<body>")
    prefix = html_document[:body_start]
    suffix = html_document[body_close:]
    body = html_document[body_start:body_close]
    parts = [
        part for part in body.split(SECTION_1C_PDF_PAGE_BREAK) if part.strip()
    ]
    return [f"{prefix}{part}{suffix}" for part in parts] or [html_document]


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
        writer = fitz.DocumentWriter(str(dest))
        # One Story per segment: a new Story always begins on a new page, which
        # makes the section-bundle page splits deterministic.
        for segment in _split_pdf_page_segments(html_document):
            story = fitz.Story(html=segment)
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


def rasterize_pdf_page_to_png(
    pdf_path: str | Path,
    png_path: str | Path,
    *,
    page_index: int = 0,
    dpi: int = 150,
) -> Path | None:
    """Rasterize one PDF page to PNG at a fixed DPI. Returns path or ``None``."""
    try:
        import fitz
    except Exception as exc:  # noqa: BLE001
        logger.warning("PNG rasterization unavailable (pymupdf not importable): %s", exc)
        return None

    src = Path(pdf_path)
    dest = Path(png_path)
    if not src.is_file():
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with fitz.open(str(src)) as doc:
            if page_index < 0 or page_index >= doc.page_count:
                return None
            page = doc[page_index]
            pix = page.get_pixmap(dpi=dpi, alpha=False)
            pix.save(str(dest))
        return dest
    except Exception as exc:  # noqa: BLE001
        logger.warning("PNG rasterization failed: %s", exc)
        if dest.exists():
            try:
                dest.unlink()
            except OSError:
                pass
        return None


def clear_stale_bundle_pngs(output_dir: str | Path) -> None:
    """Remove prior section-bundle PNG stems from ``output_dir`` only."""
    out = Path(output_dir)
    if not out.is_dir():
        return
    for path in out.glob("section_1.png"):
        path.unlink(missing_ok=True)
    for path in out.glob("section_1_page_*.png"):
        path.unlink(missing_ok=True)
    for path in out.glob("section_1_contact_sheet.png"):
        path.unlink(missing_ok=True)


def rasterize_pdf_pages_to_pngs(
    pdf_path: str | Path,
    output_dir: str | Path,
    *,
    stem: str = "section_1",
    dpi: int = 150,
) -> list[Path]:
    """Rasterize every PDF page.

    One page → ``{stem}.png``. Multiple pages → ``{stem}_page_1.png`` …
    ``{stem}_page_<n>.png``. Clears prior matching PNGs in ``output_dir`` first.
    """
    clear_stale_bundle_pngs(output_dir)
    try:
        import fitz
    except Exception as exc:  # noqa: BLE001
        logger.warning("PNG rasterization unavailable (pymupdf not importable): %s", exc)
        return []

    src = Path(pdf_path)
    out = Path(output_dir)
    if not src.is_file():
        return []
    out.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    try:
        with fitz.open(str(src)) as doc:
            count = doc.page_count
            for index in range(count):
                if count == 1:
                    dest = out / f"{stem}.png"
                else:
                    dest = out / f"{stem}_page_{index + 1}.png"
                page = doc[index]
                pix = page.get_pixmap(dpi=dpi, alpha=False)
                pix.save(str(dest))
                written.append(dest)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Multi-page PNG rasterization failed: %s", exc)
        return written
    return written


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
    report_sections: Iterable[ReportSection] | None = None,
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
    """Build a ReportDocument from evidence (and optional synthesis) and write it.

    When ``report_sections`` is omitted, behavior matches the evidence-only path.
    Synthesized prose is preferred narrative; evidence blocks remain supporting.
    """
    from gene_dossier.report_schema import build_report_document

    doc = build_report_document(
        dossier_run_id=dossier_run_id,
        gene_symbol=gene_symbol,
        evidence_records=evidence_records,
        report_sections=report_sections,
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
    "rasterize_pdf_page_to_png",
    "rasterize_pdf_pages_to_pngs",
    "clear_stale_bundle_pngs",
    "sanitize_polished_citation_tokens",
    "render_rancho_section_fragment",
    "write_rancho_report",
    "build_and_write_rancho_report",
]
