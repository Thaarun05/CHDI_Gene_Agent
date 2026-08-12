import { useEffect, useState } from 'react'
import { useParams } from 'react-router-dom'
import { getReport } from '@/api/client'
import type { ReportArtifact } from '@/api/types'
import { ReportViewer } from '@/components/ReportViewer'
import { ErrorState, LoadingSkeleton } from '@/components/EmptyState'
import { StatusBadge } from '@/components/StatusBadge'

export function ReportViewerPage() {
  const { id = '' } = useParams()
  const [report, setReport] = useState<ReportArtifact | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    setLoading(true)
    getReport(id)
      .then(setReport)
      .catch((e: Error) => setError(e.message))
      .finally(() => setLoading(false))
  }, [id])

  if (loading) return <LoadingSkeleton rows={6} />
  if (error) return <ErrorState title="Report unavailable" description={error} />
  if (!report) return <ErrorState title="Report not found" />

  return (
    <div className="fade-in space-y-6">
      <div className="flex flex-wrap items-end justify-between gap-3">
        <div>
          <p className="text-xs text-text-muted">
            {report.geneSymbol} · {report.reportOrigin === 'accepted' ? 'Accepted' : 'Generated'}
          </p>
          <h1 className="mt-1 text-2xl font-medium text-text">{report.title}</h1>
        </div>
        <StatusBadge status={report.status} />
      </div>
      <ReportViewer report={report} />
    </div>
  )
}
