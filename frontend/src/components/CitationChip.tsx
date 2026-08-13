import { useEvidenceDrawer } from '@/context/EvidenceDrawerContext'
import { cn } from '@/lib/utils'

export function CitationChip({
  label,
  evidenceReference,
  className,
}: {
  label: string
  evidenceReference: string
  className?: string
}) {
  const { openEvidence } = useEvidenceDrawer()
  return (
    <button
      type="button"
      onClick={() => void openEvidence(evidenceReference)}
      className={cn(
        'inline-flex items-center rounded-md border border-accent/30 bg-accent/10 px-2 py-0.5 text-xs font-medium text-accent transition hover:bg-accent/20',
        className,
      )}
    >
      [{label}]
    </button>
  )
}
