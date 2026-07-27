#!/usr/bin/env python3
"""Offline Rancho section preview (production builders only).

Produces HTML plus US Letter PDF and a fixed-DPI PNG raster when PyMuPDF is
available.

Example::

    .venv/bin/python scripts/render_rancho_section_preview.py \\
      --gene SREBF2 \\
      --section 1a \\
      --dossier-run-id 9923bf6d326246a4a9abb6d56053c898
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from gene_dossier.db import list_evidence_for_run, session_scope  # noqa: E402
from gene_dossier.rancho_report import (  # noqa: E402
    rasterize_pdf_page_to_png,
    render_rancho_pdf,
    render_rancho_section_fragment,
)
from gene_dossier.report_presentation import (  # noqa: E402
    build_section_presentation,
)
from gene_dossier.report_schema import build_report_document  # noqa: E402


def _parse_section_key(section: str) -> tuple[int, str]:
    text = (section or "").strip().lower()
    if text in {"1a", "1.a"}:
        return 1, "a"
    if len(text) >= 2 and text[0].isdigit() and text[-1].isalpha():
        return int(text[:-1]), text[-1]
    raise SystemExit(
        f"Unsupported section key {section!r}. Use e.g. '1a'."
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Render a single Rancho report section preview from stored evidence."
    )
    parser.add_argument("--gene", required=True, help="Gene symbol (e.g. SREBF2)")
    parser.add_argument("--section", required=True, help="Section key (e.g. 1a)")
    parser.add_argument(
        "--dossier-run-id",
        required=True,
        help="Existing dossier run id whose evidence should be loaded",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_REPO_ROOT / "data" / "outputs" / "section_previews",
        help="Preview output directory",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="PNG rasterization DPI (default: 150)",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Skip PDF/PNG generation",
    )
    args = parser.parse_args(argv)

    section_number, subsection_key = _parse_section_key(args.section)
    gene = args.gene.strip()
    run_id = args.dossier_run_id.strip()

    with session_scope() as session:
        evidence = list_evidence_for_run(session, run_id)

    if not evidence:
        raise SystemExit(
            f"No evidence records found for dossier_run_id={run_id!r}. "
            "Confirm the run exists in the configured database."
        )

    # Diagnostics from the presentation builder (not rendered into HTML).
    presentation = build_section_presentation(
        section_key=args.section,
        gene_symbol=gene,
        evidence_records=evidence,
    )

    doc = build_report_document(
        dossier_run_id=run_id,
        gene_symbol=gene,
        evidence_records=evidence,
        report_sections=None,
        curator="Gene Dossier Platform",
    )
    html = render_rancho_section_fragment(
        document=doc,
        section_number=section_number,
        subsection_key=subsection_key,
        show_cover_logos=False,
    )

    out_dir: Path = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{gene}_{args.section}_gene_aliases" if args.section.lower() in {
        "1a",
        "1.a",
    } else (
        f"{gene}_{args.section}_conservation"
        if args.section.lower() in {"1b", "1.b"}
        else f"{gene}_{args.section}"
    )
    html_path = out_dir / f"{stem}.html"
    html_path.write_text(html, encoding="utf-8")
    print(f"Wrote {html_path}")

    pdf_path = out_dir / f"{stem}.pdf"
    png_path = out_dir / f"{stem}.png"
    if not args.no_pdf:
        written_pdf = render_rancho_pdf(
            html,
            pdf_path,
            page_size="letter",
            stamp_page_chrome=True,
            stamp_cover=False,
        )
        if written_pdf is not None:
            print(f"Wrote {written_pdf}")
            written_png = rasterize_pdf_page_to_png(
                written_pdf, png_path, page_index=0, dpi=args.dpi
            )
            if written_png is not None:
                print(f"Wrote {written_png} ({args.dpi} DPI)")
            else:
                print("PNG rasterization skipped (PyMuPDF unavailable or failed).")
        else:
            print("PDF export skipped (PyMuPDF unavailable or failed).")

    print(f"Evidence records loaded: {len(evidence)}")
    print(f"Presentation blocks: {len(presentation.blocks)}")
    print("Diagnostics:")
    if not presentation.diagnostics:
        print("  (none)")
    for diag in presentation.diagnostics:
        print(f"  [{diag.severity}] {diag.field}: {diag.reason}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
