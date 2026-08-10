import { CheckCircle2, Circle, Loader2 } from 'lucide-react'
import { cn } from '@/lib/utils'

export function AgentActivity({
  steps,
  className,
}: {
  steps: string[]
  className?: string
}) {
  return (
    <div className={cn('surface-card p-4', className)}>
      <p className="mb-3 text-xs font-medium tracking-wide text-text-muted uppercase">
        Agent activity
      </p>
      <ul className="space-y-2">
        {steps.map((step, i) => {
          const isLast = i === steps.length - 1
          return (
            <li key={step} className="flex items-center gap-2.5 text-sm">
              {isLast ? (
                <Loader2 className="size-3.5 shrink-0 animate-spin text-accent-secondary" />
              ) : (
                <CheckCircle2 className="size-3.5 shrink-0 text-accent" />
              )}
              <span className={cn(isLast ? 'text-text' : 'text-text-secondary')}>{step}</span>
            </li>
          )
        })}
      </ul>
    </div>
  )
}

export function AgentActivityIdle({ className }: { className?: string }) {
  return (
    <div className={cn('surface-card p-4', className)}>
      <p className="mb-2 text-xs font-medium tracking-wide text-text-muted uppercase">
        Agent activity
      </p>
      <div className="flex items-center gap-2 text-sm text-text-muted">
        <Circle className="size-3.5" />
        Waiting for a question
      </div>
    </div>
  )
}
