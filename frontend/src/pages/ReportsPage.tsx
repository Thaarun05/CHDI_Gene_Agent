import { useEffect, useState } from 'react'
import { Link } from 'react-router-dom'
import { listReports } from '@/api/client'
import type { ReportArtifact } from '@/api/types'
import { StatusBadge } from '@/components/StatusBadge'
import { ErrorState, LoadingSkeleton } from '@/components/EmptyState'

export function ReportsPage() {
  const [reports, setReports] = useState<ReportArtifact[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const acceptedReports = reports.filter((report) => report.reportOrigin === 'accepted')
  const generatedReports = reports.filter((report) => report.reportOrigin === 'generated')

  useEffect(() => {
    listReports()
      .then(setReports)
      .catch((e: Error) => setError(e.message))
      .finally(() => setLoading(false))
  }, [])

  return (
    <div className="fade-in space-y-6">
      <div>
        <h1 className="text-2xl font-medium tracking-tight text-text">Reports</h1>
        <p className="mt-2 text-sm text-text-secondary">
          HD-focused gene dossiers and related artifacts.
        </p>
      </div>

      {loading && <LoadingSkeleton rows={4} />}
      {error && <ErrorState description={error} />}

      <ReportGroup title="Accepted" reports={acceptedReports} />
      <ReportGroup title="Generated" reports={generatedReports} />
    </div>
  )
}

function ReportGroup({ title, reports }: { title: string; reports: ReportArtifact[] }) {
  if (reports.length === 0 && title === 'Generated') {
    return (
      <section>
        <h2 className="mb-3 text-sm font-medium text-text-secondary">{title}</h2>
        <p className="text-sm text-text-muted">No generated dossiers yet.</p>
      </section>
    )
  }

  return (
    <section>
      <h2 className="mb-3 text-sm font-medium text-text-secondary">{title}</h2>
      <div className="grid gap-4 sm:grid-cols-2">
        {reports.map((report) => (
          <article key={report.id} className="surface-card p-5">
            <div className="flex items-start justify-between gap-3">
              <div>
                <p className="text-xs text-text-muted">{report.geneSymbol}</p>
                <h3 className="mt-1 text-base font-medium text-text">{report.title}</h3>
              </div>
              <StatusBadge status={report.status} />
            </div>
            <p className="mt-3 text-xs text-text-muted">
              Created {new Date(report.createdAt).toLocaleString()} · {report.sections.length}{' '}
              sections
            </p>
            <div className="mt-4 flex flex-wrap gap-2">
              {report.htmlUrl ? (
                <Link
                  to={`/reports/${report.id}`}
                  className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text transition hover:border-accent/40"
                >
                  Open
                </Link>
              ) : (
                <span className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text-muted">
                  Report unavailable
                </span>
              )}
              {report.pdfUrl ? (
                <a
                  href={report.pdfUrl}
                  className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text-secondary"
                >
                  Download PDF
                </a>
              ) : (
                <span className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text-muted">
                  {report.reportOrigin === 'generated' ? 'PDF not generated' : 'PDF unavailable'}
                </span>
              )}
              <Link
                to={`/evidence?gene=${encodeURIComponent(report.geneSymbol)}`}
                className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text-secondary"
              >
                Evidence
              </Link>
            </div>
          </article>
        ))}
      </div>
    </section>
  )
}
