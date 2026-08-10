import type { EvidenceRecord } from '@/api/types'
import { StatusBadge } from '@/components/StatusBadge'
import { useEvidenceDrawer } from '@/context/EvidenceDrawerContext'
import { cn } from '@/lib/utils'

export function EvidenceCard({
  record,
  className,
}: {
  record: EvidenceRecord
  className?: string
}) {
  const { openEvidence } = useEvidenceDrawer()
  return (
    <button
      type="button"
      onClick={() => void openEvidence(record)}
      className={cn(
        'surface-card w-full p-4 text-left transition duration-150 hover:border-accent/30',
        className,
      )}
    >
      <div className="mb-2 flex items-center justify-between gap-2">
        <span className="text-xs font-medium text-accent-secondary">{record.sourceName}</span>
        {record.status && <StatusBadge status={record.status} />}
      </div>
      <p className="text-sm leading-relaxed text-text">{record.displayText}</p>
      <p className="mt-2 text-xs text-text-muted">
        {record.section} · {record.id}
      </p>
    </button>
  )
}
