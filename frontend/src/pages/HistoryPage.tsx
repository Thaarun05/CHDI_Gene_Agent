import { useEffect, useState } from 'react'
import { listHistory } from '@/api/client'
import type { HistoryItem } from '@/api/types'
import { StatusBadge } from '@/components/StatusBadge'
import { ErrorState, LoadingSkeleton } from '@/components/EmptyState'

export function HistoryPage() {
  const [items, setItems] = useState<HistoryItem[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    listHistory()
      .then(setItems)
      .catch((e: Error) => setError(e.message))
      .finally(() => setLoading(false))
  }, [])

  return (
    <div className="fade-in space-y-6">
      <div>
        <h1 className="text-2xl font-medium tracking-tight text-text">History</h1>
        <p className="mt-2 text-sm text-text-secondary">Workflow activity across the research workspace.</p>
      </div>

      {loading && <LoadingSkeleton rows={5} />}
      {error && <ErrorState description={error} />}

      <ul className="surface-card divide-y divide-border overflow-hidden">
        {items.map((item) => (
          <li
            key={item.id}
            className="flex flex-wrap items-center justify-between gap-3 px-5 py-4"
          >
            <div>
              <p className="text-sm font-medium text-text">{item.geneLabel}</p>
              <p className="mt-0.5 text-sm text-text-secondary">{item.workflow}</p>
              <p className="mt-1 text-xs text-text-muted">
                {new Date(item.createdAt).toLocaleString()}
              </p>
            </div>
            <StatusBadge status={item.status} />
          </li>
        ))}
      </ul>
    </div>
  )
}
