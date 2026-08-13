import { useEffect, useRef, useState, type FormEvent, type KeyboardEvent } from 'react'
import { ChevronDown, Paperclip, Send, Square } from 'lucide-react'
import type { ResearchMode } from '@/api/types'
import { cn } from '@/lib/utils'

const SUGGESTIONS = ['SREBF2', 'CDH10', 'Chemical tools', 'Protein interactions']
const CONTEXT_OPTIONS = ['SREBF2', 'CDH10'] as const
const RESEARCH_MODES: Array<{ value: ResearchMode; label: string }> = [
  { value: 'auto', label: 'Auto' },
  { value: 'deep_research', label: 'Deep Research' },
  { value: 'stored_only', label: 'Stored Evidence Only' },
]

export function SearchComposer({
  value,
  onChange,
  onSubmit,
  selectedGene = null,
  onSelectGene,
  researchMode = 'auto',
  onSelectResearchMode,
  submitting = false,
  onCancel,
  placeholder = 'Ask about a gene, target, pathway, compound, or evidence question...',
  className,
}: {
  value: string
  onChange: (v: string) => void
  onSubmit: () => void
  selectedGene?: string | null
  onSelectGene?: (gene: string | null) => void
  researchMode?: ResearchMode
  onSelectResearchMode?: (mode: ResearchMode) => void
  submitting?: boolean
  onCancel?: () => void
  placeholder?: string
  className?: string
}) {
  const [geneOpen, setGeneOpen] = useState(false)
  const [modeOpen, setModeOpen] = useState(false)
  const textareaRef = useRef<HTMLTextAreaElement>(null)

  useEffect(() => {
    const textarea = textareaRef.current
    if (!textarea) return
    textarea.style.height = '0px'
    textarea.style.height = `${Math.min(Math.max(textarea.scrollHeight, 92), 240)}px`
    textarea.style.overflowY = textarea.scrollHeight > 240 ? 'auto' : 'hidden'
  }, [value])

  function handleSubmit(e?: FormEvent) {
    e?.preventDefault()
    if (!value.trim() || submitting) return
    onSubmit()
  }

  function onKeyDown(e: KeyboardEvent<HTMLTextAreaElement>) {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault()
      handleSubmit()
    }
  }

  const modeLabel = RESEARCH_MODES.find((item) => item.value === researchMode)?.label ?? 'Auto'

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
        className="rounded-[22px] border border-border bg-card shadow-[0_24px_80px_-40px_rgba(0,0,0,0.8)]"
      >
        <textarea
          ref={textareaRef}
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={onKeyDown}
          rows={3}
          placeholder={placeholder}
          className="block min-h-[92px] max-h-60 w-full resize-none bg-transparent px-5 pt-5 pb-3 text-[15px] leading-relaxed text-text outline-none placeholder:text-text-muted"
        />

        <div className="flex flex-wrap items-end justify-between gap-3 border-t border-border/60 px-3 py-3">
          <div className="flex min-w-0 flex-wrap items-center gap-1.5">
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
                onClick={() => setGeneOpen((open) => !open)}
                className="inline-flex max-w-full items-center gap-1 rounded-lg border border-border px-2.5 py-1.5 text-xs text-text-secondary transition hover:bg-white/5"
              >
                {selectedGene ? `Context Gene: ${selectedGene}` : 'No context gene'}
                <ChevronDown className="size-3 shrink-0" />
              </button>
              {geneOpen && (
                <div className="absolute bottom-full left-0 z-20 mb-1 w-44 overflow-hidden rounded-lg border border-border bg-card-elevated shadow-xl">
                  <button
                    type="button"
                    className="block w-full px-3 py-2 text-left text-xs text-text-secondary hover:bg-white/5 hover:text-text"
                    onClick={() => {
                      onSelectGene?.(null)
                      setGeneOpen(false)
                    }}
                  >
                    No context gene
                  </button>
                  {CONTEXT_OPTIONS.map((gene) => (
                    <button
                      key={gene}
                      type="button"
                      className="block w-full px-3 py-2 text-left text-xs text-text-secondary hover:bg-white/5 hover:text-text"
                      onClick={() => {
                        onSelectGene?.(gene)
                        setGeneOpen(false)
                      }}
                    >
                      {gene}
                    </button>
                  ))}
                </div>
              )}
            </div>
            <div className="relative">
              <button
                type="button"
                onClick={() => setModeOpen((open) => !open)}
                className="inline-flex items-center gap-1 rounded-lg border border-border px-2.5 py-1.5 text-xs text-text-secondary transition hover:bg-white/5"
              >
                {modeLabel}
                <ChevronDown className="size-3" />
              </button>
              {modeOpen && (
                <div className="absolute bottom-full left-0 z-20 mb-1 w-48 overflow-hidden rounded-lg border border-border bg-card-elevated shadow-xl">
                  {RESEARCH_MODES.map((mode) => (
                    <button
                      key={mode.value}
                      type="button"
                      className="block w-full px-3 py-2 text-left text-xs text-text-secondary hover:bg-white/5 hover:text-text"
                      onClick={() => {
                        onSelectResearchMode?.(mode.value)
                        setModeOpen(false)
                      }}
                    >
                      {mode.label}
                    </button>
                  ))}
                </div>
              )}
            </div>
          </div>
          {submitting ? (
            <button
              type="button"
              onClick={onCancel}
              className="flex size-10 shrink-0 items-center justify-center rounded-full border border-accent/40 text-accent transition hover:bg-accent/10"
              aria-label="Cancel request"
            >
              <Square className="size-3.5 fill-current" />
            </button>
          ) : (
            <button
              type="submit"
              disabled={!value.trim()}
              className="flex size-10 shrink-0 items-center justify-center rounded-full bg-accent text-bg transition hover:brightness-110 disabled:cursor-not-allowed disabled:opacity-40"
              aria-label="Send"
            >
              <Send className="size-4" />
            </button>
          )}
        </div>
      </form>

      <div className="mt-4 flex flex-wrap justify-center gap-2">
        {SUGGESTIONS.map((suggestion) => (
          <button
            key={suggestion}
            type="button"
            onClick={() => onChange(suggestion)}
            className="rounded-full border border-border bg-bg-secondary/80 px-3.5 py-1.5 text-xs text-text-secondary transition hover:border-accent/35 hover:text-text"
          >
            {suggestion}
          </button>
        ))}
      </div>
    </div>
  )
}
