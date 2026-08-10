import type { EvidenceCoverageRow } from '@/api/types'
import { StatusBadge } from '@/components/StatusBadge'

export function EvidenceCoverage({ rows }: { rows: EvidenceCoverageRow[] }) {
  return (
    <div className="surface-card overflow-hidden">
      <div className="border-b border-border px-4 py-3">
        <h3 className="text-sm font-medium text-text">Evidence Coverage</h3>
        <p className="mt-0.5 text-xs text-text-muted">
          Status vocabulary distinguishes No Results, Unavailable, Missing, and Source Error.
        </p>
      </div>
      <ul className="divide-y divide-border">
        {rows.map((row) => (
          <li key={row.category} className="flex items-center justify-between gap-3 px-4 py-3">
            <div className="min-w-0">
              <p className="text-sm text-text">{row.category}</p>
              {row.detail && <p className="truncate text-xs text-text-muted">{row.detail}</p>}
            </div>
            <StatusBadge status={row.status} />
          </li>
        ))}
      </ul>
    </div>
  )
}
