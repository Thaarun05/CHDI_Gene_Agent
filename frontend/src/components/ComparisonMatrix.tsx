import type { ComparisonResponse } from '@/api/types'
import { StatusBadge } from '@/components/StatusBadge'
import { useEvidenceDrawer } from '@/context/EvidenceDrawerContext'
import { cn } from '@/lib/utils'

export function ComparisonMatrix({
  data,
  className,
}: {
  data: ComparisonResponse
  className?: string
}) {
  const { openEvidence } = useEvidenceDrawer()

  return (
    <div className={cn('surface-card overflow-hidden', className)}>
      <div className="overflow-x-auto">
        <table className="w-full min-w-[640px] border-collapse text-left text-sm">
          <thead>
            <tr className="border-b border-border bg-bg-secondary/60">
              <th className="px-4 py-3 font-medium text-text-secondary">Evidence Dimension</th>
              {data.genes.map((g) => (
                <th key={g} className="px-4 py-3 font-medium text-text">
                  {g}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {data.matrix.map((row) => (
              <tr key={row.dimension} className="border-b border-border/80 last:border-0">
                <td className="px-4 py-3.5 align-top text-text-secondary">{row.dimension}</td>
                {data.genes.map((g) => {
                  const cell = row.cells[g]
                  if (!cell) {
                    return (
                      <td key={g} className="px-4 py-3.5 text-text-muted">
                        —
                      </td>
                    )
                  }
                  return (
                    <td key={g} className="px-4 py-3.5 align-top">
                      <button
                        type="button"
                        className="w-full rounded-xl border border-transparent p-2 text-left transition hover:border-border hover:bg-white/[0.03]"
                        onClick={() => {
                          const id = cell.evidenceRecordIds[0]
                          if (id) void openEvidence(id)
                        }}
                      >
                        <StatusBadge status={cell.status} />
                        <p className="mt-2 text-sm text-text">{cell.summary}</p>
                        <p className="mt-1 text-xs text-text-muted">
                          {cell.evidenceCount} evidence
                        </p>
                      </button>
                    </td>
                  )
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  )
}
