import { useEffect, useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { getEvidenceCoverage, getGene, listGeneEvidence } from '@/api/client'
import type { EvidenceCoverageRow, EvidenceRecord, Gene } from '@/api/types'
import { GeneHeader } from '@/components/GeneHeader'
import { EvidenceCoverage } from '@/components/EvidenceCoverage'
import { EvidenceCard } from '@/components/EvidenceCard'
import { ErrorState, LoadingSkeleton } from '@/components/EmptyState'
import { cn } from '@/lib/utils'

const TABS = [
  'Overview',
  'Evidence',
  'Expression',
  'Perturbations',
  'Interactions',
  'Chemical Tools',
  'Publications',
  'Reports',
  'Provenance',
] as const

type Tab = (typeof TABS)[number]

export function GeneWorkspacePage() {
  const { symbol = 'SREBF2' } = useParams()
  const [gene, setGene] = useState<Gene | null>(null)
  const [coverage, setCoverage] = useState<EvidenceCoverageRow[]>([])
  const [evidence, setEvidence] = useState<EvidenceRecord[]>([])
  const [tab, setTab] = useState<Tab>('Overview')
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    let cancelled = false
    setLoading(true)
    setError(null)
    Promise.all([
      getGene(symbol),
      getEvidenceCoverage(symbol),
      listGeneEvidence(symbol),
    ])
      .then(([g, c, e]) => {
        if (cancelled) return
        setGene(g)
        setCoverage(c.rows)
        setEvidence(e.records)
      })
      .catch((err: Error) => {
        if (!cancelled) setError(err.message)
      })
      .finally(() => {
        if (!cancelled) setLoading(false)
      })
    return () => {
      cancelled = true
    }
  }, [symbol])

  if (loading) return <LoadingSkeleton rows={8} className="mt-8" />
  if (error) return <ErrorState title="Unable to load gene" description={error} />
  if (!gene) return <ErrorState title="Gene not found" />

  return (
    <div className="fade-in space-y-8">
      <div className="flex flex-wrap gap-2 text-sm items-center">
        <Link to="/genes" className="text-text-secondary hover:text-text">
          Genes
        </Link>
        <span className="text-text-muted">/</span>
        <span className="text-text font-medium">{gene.symbol}</span>
      </div>

      <GeneHeader gene={gene} />

      <div className="flex gap-1 overflow-x-auto border-b border-border pb-px">
        {TABS.map((t) => (
          <button
            key={t}
            type="button"
            onClick={() => setTab(t)}
            className={cn(
              'shrink-0 border-b-2 px-3 py-2.5 text-sm transition',
              tab === t
                ? 'border-accent text-text'
                : 'border-transparent text-text-muted hover:text-text-secondary',
            )}
          >
            {t}
          </button>
        ))}
      </div>

      {tab === 'Overview' && (
        <div className="grid gap-4 lg:grid-cols-2">
          <div className="surface-card p-5">
            <h3 className="text-sm font-medium text-text">Target Summary</h3>
            <p className="mt-3 text-sm leading-relaxed text-text-secondary">{gene.summary}</p>
            <p className="mt-3 text-xs text-text-muted">Live FastAPI provenance database data</p>
          </div>
          <EvidenceCoverage rows={coverage} />
          <OverviewTile
            title="Expression"
            body="Tissue and cell expression sections available from GTEx / Human Brain Transcriptome artifacts in the stored dossier."
          />
          <OverviewTile
            title="Protein Interactions"
            body="STRING and BioGRID partner tables present for this gene in the demo coverage set."
          />
          <OverviewTile
            title="Chemical Perturbations"
            body="Baseline CTD chemical-gene interaction evidence is not present in the accepted demo evidence run."
          />
          <OverviewTile
            title="Chemical Tools"
            body={
              gene.symbol === 'SREBF2'
                ? 'ChEMBL workbook, PubMed tools, PubChem assays, and NCATS rows available.'
                : 'Limited chemical tools: no authoritative ChEMBL target; indirect literature only.'
            }
          />
          <div className="surface-card p-5 lg:col-span-2">
            <h3 className="mb-3 text-sm font-medium text-text">Recent Evidence</h3>
            <div className="grid gap-3 sm:grid-cols-2">
              {evidence.slice(0, 4).map((rec) => (
                <EvidenceCard key={rec.id} record={rec} />
              ))}
            </div>
          </div>
        </div>
      )}

      {tab === 'Evidence' && (
        <div className="grid gap-3 sm:grid-cols-2">
          {evidence.map((rec) => (
            <EvidenceCard key={rec.id} record={rec} />
          ))}
        </div>
      )}

      {tab !== 'Overview' && tab !== 'Evidence' && (
        <div className="surface-card p-8 text-center">
          <p className="text-sm text-text-secondary">
            {tab} workspace panel — layout shell ready for FastAPI evidence payloads.
          </p>
          <p className="mt-2 text-xs text-text-muted">Demo placeholder for {gene.symbol}</p>
        </div>
      )}
    </div>
  )
}

function OverviewTile({ title, body }: { title: string; body: string }) {
  return (
    <div className="surface-card p-5">
      <h3 className="text-sm font-medium text-text">{title}</h3>
      <p className="mt-3 text-sm leading-relaxed text-text-secondary">{body}</p>
    </div>
  )
}
