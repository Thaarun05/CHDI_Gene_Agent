import { cn } from '@/lib/utils'

export function EmptyState({
  title,
  description,
  className,
}: {
  title: string
  description?: string
  className?: string
}) {
  return (
    <div
      className={cn(
        'flex flex-col items-center justify-center rounded-2xl border border-dashed border-border bg-card/40 px-6 py-12 text-center',
        className,
      )}
    >
      <p className="text-sm font-medium text-text">{title}</p>
      {description && <p className="mt-2 max-w-sm text-sm text-text-muted">{description}</p>}
    </div>
  )
}

export function ErrorState({
  title = 'Something went wrong',
  description,
  className,
}: {
  title?: string
  description?: string
  className?: string
}) {
  return (
    <div
      className={cn(
        'rounded-2xl border border-danger/30 bg-danger/5 px-5 py-4 text-sm',
        className,
      )}
    >
      <p className="font-medium text-danger">{title}</p>
      {description && <p className="mt-1 text-text-secondary">{description}</p>}
    </div>
  )
}

export function LoadingSkeleton({
  rows = 4,
  className,
}: {
  rows?: number
  className?: string
}) {
  return (
    <div className={cn('space-y-3', className)}>
      {Array.from({ length: rows }).map((_, i) => (
        <div
          key={i}
          className="h-10 animate-pulse rounded-xl bg-white/[0.04]"
          style={{ width: `${88 - i * 8}%` }}
        />
      ))}
    </div>
  )
}
