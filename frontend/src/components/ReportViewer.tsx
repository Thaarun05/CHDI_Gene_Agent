import { Download, FileArchive, Search } from 'lucide-react'
import type { ReportArtifact } from '@/api/types'
import { cn } from '@/lib/utils'

const SECTION_NAV = [
  'General Gene Information',
  'Structure',
  'Expression',
  'GEO Perturbations',
  'Transcription Factors',
  'Protein Interactions',
  'Chemical Perturbations',
  'Chemical Tools',
]

export function ReportViewer({
  report,
  className,
}: {
  report: ReportArtifact
  className?: string
}) {
  return (
    <div className={cn('flex min-h-[70vh] flex-col gap-4 lg:flex-row', className)}>
      <nav className="surface-card w-full shrink-0 p-3 lg:w-56">
        <p className="mb-2 px-2 text-xs font-medium tracking-wide text-text-muted uppercase">
          Sections
        </p>
        <ul className="space-y-0.5">
          {SECTION_NAV.map((label) => (
            <li key={label}>
              <button
                type="button"
                className="w-full rounded-lg px-2.5 py-2 text-left text-sm text-text-secondary transition hover:bg-white/5 hover:text-text"
              >
                {label}
              </button>
            </li>
          ))}
        </ul>
      </nav>

      <div className="flex min-w-0 flex-1 flex-col gap-3">
        <div className="flex flex-wrap items-center gap-2">
          <a
            href={report.pdfUrl || '#'}
            className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3 py-2 text-sm text-text transition hover:border-accent/40"
          >
            <Download className="size-3.5" />
            Download PDF
          </a>
          <button
            type="button"
            className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3 py-2 text-sm text-text-secondary"
          >
            <FileArchive className="size-3.5" />
            Supplementary Files
          </button>
          <button
            type="button"
            className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3 py-2 text-sm text-text-secondary"
          >
            <Search className="size-3.5" />
            Inspect Evidence
          </button>
        </div>

        <div className="surface-card flex-1 overflow-hidden p-1">
          {report.htmlUrl ? (
            <iframe
              title={`${report.geneSymbol} report`}
              src={report.htmlUrl}
              className="h-[min(80vh,900px)] w-full rounded-[14px] bg-white"
            />
          ) : (
            <div className="flex h-80 items-center justify-center text-sm text-text-muted">
              Report HTML URL not available in mock data.
            </div>
          )}
        </div>
        <p className="text-xs text-text-muted">
          Embedded Rancho HTML artifact — not re-rendered in React. Demo iframe may require the
          local artifact server.
        </p>
      </div>
    </div>
  )
}
