import type { EvidenceStatus, JobStatus } from '@/api/types'
import { cn } from '@/lib/utils'

const statusStyles: Record<EvidenceStatus | JobStatus, string> = {
  Available: 'bg-accent/15 text-accent border-accent/25',
  Limited: 'bg-accent-secondary/15 text-accent-secondary border-accent-secondary/25',
  Missing: 'bg-white/5 text-text-muted border-white/10',
  'No Results': 'bg-white/[0.04] text-text-secondary border-white/10',
  Unavailable: 'bg-warning/10 text-warning border-warning/20',
  'Source Error': 'bg-danger/10 text-danger border-danger/25',
  Stale: 'bg-warning/10 text-warning border-warning/20',
  Running: 'bg-accent-secondary/15 text-accent-secondary border-accent-secondary/25',
  Queued: 'bg-white/5 text-text-muted border-white/10',
  Completed: 'bg-accent/15 text-accent border-accent/25',
  Partial: 'bg-warning/10 text-warning border-warning/20',
  Failed: 'bg-danger/10 text-danger border-danger/25',
}

const dotStyles: Record<EvidenceStatus | JobStatus, string> = {
  Available: 'bg-accent',
  Limited: 'bg-accent-secondary',
  Missing: 'bg-text-muted',
  'No Results': 'bg-text-secondary',
  Unavailable: 'bg-warning',
  'Source Error': 'bg-danger',
  Stale: 'bg-warning',
  Running: 'bg-accent-secondary animate-pulse',
  Queued: 'bg-text-muted',
  Completed: 'bg-accent',
  Partial: 'bg-warning',
  Failed: 'bg-danger',
}

export function StatusBadge({
  status,
  className,
}: {
  status: EvidenceStatus | JobStatus
  className?: string
}) {
  return (
    <span
      className={cn(
        'inline-flex items-center gap-1.5 rounded-full border px-2.5 py-0.5 text-xs font-medium',
        statusStyles[status],
        className,
      )}
    >
      <span className={cn('size-1.5 rounded-full', dotStyles[status])} />
      {status}
    </span>
  )
}
