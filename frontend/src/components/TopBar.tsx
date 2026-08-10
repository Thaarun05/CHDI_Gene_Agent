import { MoreHorizontal, ChevronDown } from 'lucide-react'
import { useState } from 'react'
import { cn } from '@/lib/utils'

const MODES = ['Gene Intelligence', 'Evidence Assistant', 'Target Discovery'] as const

export function TopBar({ className }: { className?: string }) {
  const [mode, setMode] = useState<(typeof MODES)[number]>('Gene Intelligence')
  const [open, setOpen] = useState(false)

  return (
    <div className={cn('flex items-center justify-between', className)}>
      <div className="relative">
        <button
          type="button"
          onClick={() => setOpen((o) => !o)}
          className="inline-flex items-center gap-1.5 rounded-xl border border-border bg-card/60 px-3 py-2 text-sm text-text transition hover:bg-card"
        >
          {mode}
          <ChevronDown className="size-3.5 text-text-muted" />
        </button>
        {open && (
          <div className="absolute top-full left-0 z-20 mt-1 min-w-[200px] overflow-hidden rounded-xl border border-border bg-card-elevated shadow-xl">
            {MODES.map((m) => (
              <button
                key={m}
                type="button"
                className={cn(
                  'block w-full px-3.5 py-2.5 text-left text-sm transition hover:bg-white/5',
                  m === mode ? 'text-accent' : 'text-text-secondary',
                )}
                onClick={() => {
                  setMode(m)
                  setOpen(false)
                }}
              >
                {m}
              </button>
            ))}
          </div>
        )}
      </div>
      <button
        type="button"
        className="rounded-lg p-2 text-text-muted transition hover:bg-white/5 hover:text-text"
        aria-label="Menu"
      >
        <MoreHorizontal className="size-4" />
      </button>
    </div>
  )
}
