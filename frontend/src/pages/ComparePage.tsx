import { useEffect, useState } from 'react'
import { useSearchParams } from 'react-router-dom'
import { compareGenes } from '@/api/client'
import type { ComparisonResponse } from '@/api/types'
import { ComparisonMatrix } from '@/components/ComparisonMatrix'
import { ErrorState, LoadingSkeleton } from '@/components/EmptyState'
import { Plus } from 'lucide-react'

export function ComparePage() {
  const [params] = useSearchParams()
  const initial = (params.get('genes') || 'SREBF2,CDH10')
    .split(',')
    .map((g) => g.trim().toUpperCase())
    .filter(Boolean)

  const [genes, setGenes] = useState<string[]>(initial.length ? initial : ['SREBF2', 'CDH10'])
  const [data, setData] = useState<ComparisonResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    setLoading(true)
    compareGenes(genes)
      .then(setData)
      .catch((e: Error) => setError(e.message))
      .finally(() => setLoading(false))
  }, [genes])

  function addGene() {
    const next = genes.includes('CDH10') && genes.includes('SREBF2') ? null : 'CDH10'
    if (!next || genes.includes(next)) return
    setGenes([...genes, next])
  }

  return (
    <div className="fade-in space-y-8">
      <div className="flex flex-wrap items-end justify-between gap-4">
        <div>
          <h1 className="text-2xl font-medium tracking-tight text-text">Compare Genes</h1>
          <p className="mt-2 text-sm text-text-secondary">
            Evidence matrix across biological dimensions. No AI scores or unexplained rankings.
          </p>
        </div>
        <button
          type="button"
          onClick={addGene}
          className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3.5 py-2 text-sm text-text-secondary transition hover:text-text"
        >
          <Plus className="size-3.5" />
          Add Gene
        </button>
      </div>

      <div className="flex flex-wrap gap-2">
        {genes.map((g) => (
          <span
            key={g}
            className="rounded-full border border-accent/30 bg-accent/10 px-3 py-1 text-sm text-accent"
          >
            {g}
          </span>
        ))}
      </div>

      {loading && <LoadingSkeleton rows={6} />}
      {error && <ErrorState description={error} />}
      {data && !loading && (
        <>
          <ComparisonMatrix data={data} />
          <section className="surface-card p-6">
            <h2 className="text-sm font-medium text-text">Evidence-Grounded Comparison</h2>
            <p className="mt-3 text-[15px] leading-relaxed text-text-secondary">{data.narrative}</p>
            <p className="mt-3 text-xs text-text-muted">
              Comparison based on selected provenance-backed evidence runs. No AI score.
            </p>
            <div className="mt-3 flex flex-wrap gap-2 text-xs text-text-muted">
              {Object.entries(data.evidenceUniverses).map(([geneSymbol, universe]) => (
                <span key={geneSymbol} className="rounded-full border border-border px-2.5 py-1">
                  {geneSymbol}: {universe.dossierRunIds.join(', ')}
                </span>
              ))}
            </div>
          </section>
        </>
      )}
    </div>
  )
}
