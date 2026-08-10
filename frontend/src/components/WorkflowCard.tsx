import type { LucideIcon } from 'lucide-react'
import { Link } from 'react-router-dom'
import { cn } from '@/lib/utils'

export function WorkflowCard({
  title,
  description,
  icon: Icon,
  to,
  comingNext,
  className,
}: {
  title: string
  description: string
  icon: LucideIcon
  to?: string
  comingNext?: boolean
  className?: string
}) {
  const inner = (
    <>
      <div className="mb-4 flex items-center justify-between">
        <div className="flex size-9 items-center justify-center rounded-xl border border-border bg-white/[0.03]">
          <Icon className="size-4 text-accent" strokeWidth={1.75} />
        </div>
        {comingNext && (
          <span className="rounded-full border border-border px-2 py-0.5 text-[10px] tracking-wide text-text-muted uppercase">
            Coming next
          </span>
        )}
      </div>
      <h3 className="text-[15px] font-medium text-text">{title}</h3>
      <p className="mt-2 text-sm leading-relaxed text-text-secondary">{description}</p>
    </>
  )

  const classes = cn(
    'surface-card block p-5 transition duration-200 hover:-translate-y-0.5 hover:border-accent/30 hover:shadow-[0_12px_40px_-24px_rgba(140,203,94,0.45)]',
    comingNext && 'opacity-80',
    className,
  )

  if (to && !comingNext) {
    return (
      <Link to={to} className={classes}>
        {inner}
      </Link>
    )
  }
  return <div className={classes}>{inner}</div>
}
