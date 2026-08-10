import { useState, type FormEvent, type KeyboardEvent } from 'react'
import { Paperclip, Send, ChevronDown } from 'lucide-react'
import { cn } from '@/lib/utils'

const SUGGESTIONS = ['SREBF2', 'CDH10', 'Chemical tools', 'Protein interactions']

export function SearchComposer({
  value,
  onChange,
  onSubmit,
  placeholder = 'Ask about a gene, target, pathway, compound, or evidence question...',
  className,
}: {
  value: string
  onChange: (v: string) => void
  onSubmit: () => void
  placeholder?: string
  className?: string
}) {
  const [geneOpen, setGeneOpen] = useState(false)

  function handleSubmit(e?: FormEvent) {
    e?.preventDefault()
    if (!value.trim()) return
    onSubmit()
  }

  function onKeyDown(e: KeyboardEvent<HTMLTextAreaElement>) {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault()
      handleSubmit()
    }
  }

  return (
    <div className={cn('w-full', className)}>
      <div className="mb-3 flex justify-center">
        <span className="inline-flex items-center gap-1.5 rounded-full border border-accent/25 bg-accent/10 px-3 py-1 text-xs text-accent">
          <span className="size-1.5 rounded-full bg-accent" />
          Evidence-grounded
        </span>
      </div>

      <form
        onSubmit={handleSubmit}
        className="relative rounded-[22px] border border-border bg-card shadow-[0_24px_80px_-40px_rgba(0,0,0,0.8)]"
      >
        <textarea
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={onKeyDown}
          rows={3}
          placeholder={placeholder}
          className="w-full resize-none bg-transparent px-5 pt-5 pb-16 text-[15px] leading-relaxed text-text outline-none placeholder:text-text-muted"
        />

        <div className="absolute right-3 bottom-3 left-3 flex items-center justify-between gap-2">
          <div className="flex flex-wrap items-center gap-1.5">
            <button
              type="button"
              className="rounded-lg p-2 text-text-muted transition hover:bg-white/5 hover:text-text"
              aria-label="Attach"
            >
              <Paperclip className="size-4" />
            </button>
            <div className="relative">
              <button
                type="button"
                onClick={() => setGeneOpen((o) => !o)}
                className="inline-flex items-center gap-1 rounded-lg border border-border px-2.5 py-1.5 text-xs text-text-secondary transition hover:bg-white/5"
              >
                Gene / Target
                <ChevronDown className="size-3" />
              </button>
              {geneOpen && (
                <div className="absolute bottom-full left-0 mb-1 w-36 overflow-hidden rounded-xl border border-border bg-card-elevated shadow-xl">
                  {['SREBF2', 'CDH10'].map((g) => (
                    <button
                      key={g}
                      type="button"
                      className="block w-full px-3 py-2 text-left text-xs text-text-secondary hover:bg-white/5 hover:text-text"
                      onClick={() => {
                        onChange(value ? `${value} ${g}` : g)
                        setGeneOpen(false)
                      }}
                    >
                      {g}
                    </button>
                  ))}
                </div>
              )}
            </div>
            <button
              type="button"
              className="rounded-lg border border-border px-2.5 py-1.5 text-xs text-text-secondary transition hover:bg-white/5"
            >
              Evidence Sources
            </button>
            <button
              type="button"
              className="rounded-lg border border-border px-2.5 py-1.5 text-xs text-text-secondary transition hover:bg-white/5"
            >
              Research Mode
            </button>
          </div>
          <button
            type="submit"
            disabled={!value.trim()}
            className="flex size-10 items-center justify-center rounded-full bg-accent text-bg transition hover:brightness-110 disabled:cursor-not-allowed disabled:opacity-40"
            aria-label="Send"
          >
            <Send className="size-4" />
          </button>
        </div>
      </form>

      <div className="mt-4 flex flex-wrap justify-center gap-2">
        {SUGGESTIONS.map((s) => (
          <button
            key={s}
            type="button"
            onClick={() => onChange(s)}
            className="rounded-full border border-border bg-bg-secondary/80 px-3.5 py-1.5 text-xs text-text-secondary transition hover:border-accent/35 hover:text-text"
          >
            {s}
          </button>
        ))}
      </div>
    </div>
  )
}
