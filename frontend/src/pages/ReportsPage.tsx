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

      <div className="grid gap-4 sm:grid-cols-2">
        {reports.map((r) => (
          <article key={r.id} className="surface-card p-5">
            <div className="flex items-start justify-between gap-3">
              <div>
                <p className="text-xs text-text-muted">{r.geneSymbol}</p>
                <h2 className="mt-1 text-base font-medium text-text">{r.title}</h2>
              </div>
              <StatusBadge status={r.status} />
            </div>
            <p className="mt-3 text-xs text-text-muted">
              Created {new Date(r.createdAt).toLocaleString()} · {r.sections.length} sections
            </p>
            <div className="mt-4 flex flex-wrap gap-2">
              <Link
                to={`/reports/${r.id}`}
                className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text transition hover:border-accent/40"
              >
                Open
              </Link>
              <a
                href={r.pdfUrl || '#'}
                className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text-secondary"
              >
                Download PDF
              </a>
              <Link
                to={`/evidence`}
                className="rounded-xl border border-border bg-bg-secondary px-3 py-2 text-sm text-text-secondary"
              >
                Evidence
              </Link>
            </div>
          </article>
        ))}
      </div>
    </div>
  )
}
