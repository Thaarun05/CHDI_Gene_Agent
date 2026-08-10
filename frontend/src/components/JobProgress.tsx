import type { WorkflowJob } from '@/api/types'
import { StatusBadge } from '@/components/StatusBadge'
import { cn } from '@/lib/utils'

const stageDot: Record<string, string> = {
  Complete: 'bg-accent',
  Running: 'bg-accent-secondary animate-pulse',
  Queued: 'bg-text-muted',
  Waiting: 'bg-text-muted/60',
  Failed: 'bg-danger',
}

export function JobProgress({ job, className }: { job: WorkflowJob; className?: string }) {
  return (
    <div className={cn('surface-card p-5', className)}>
      <div className="mb-4 flex items-center justify-between gap-3">
        <div>
          <h3 className="text-sm font-medium text-text">Dossier job progress</h3>
          <p className="mt-0.5 text-xs text-text-muted">
            {job.geneSymbol} · {job.id}
          </p>
        </div>
        <StatusBadge status={job.status} />
      </div>
      <ol className="space-y-2.5">
        {job.stages.map((stage) => (
          <li
            key={stage.id}
            className="flex items-center justify-between gap-3 rounded-xl border border-border/70 bg-bg-secondary/50 px-3.5 py-2.5"
          >
            <div className="flex items-center gap-3">
              <span className={cn('size-2 rounded-full', stageDot[stage.status])} />
              <span className="text-sm text-text">{stage.label}</span>
            </div>
            <span className="text-xs text-text-muted">{stage.status}</span>
          </li>
        ))}
      </ol>
    </div>
  )
}
