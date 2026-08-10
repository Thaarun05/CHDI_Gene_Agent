import { Link } from 'react-router-dom'
import { FileText, GitCompare, MessageSquareText, RefreshCw } from 'lucide-react'
import type { Gene } from '@/api/types'
import { cn } from '@/lib/utils'

export function GeneHeader({ gene, className }: { gene: Gene; className?: string }) {
  return (
    <div className={cn('space-y-4', className)}>
      <div>
        <h1 className="text-3xl font-medium tracking-tight text-text">{gene.symbol}</h1>
        <p className="mt-1 text-text-secondary">{gene.name}</p>
      </div>
      <div className="flex flex-wrap gap-2">
        {[gene.organism, `NCBI Gene ${gene.entrezGeneId}`, `UniProt ${gene.uniprotAccession}`].map(
          (chip) => (
            <span
              key={chip}
              className="rounded-full border border-border bg-card px-3 py-1 text-xs text-text-secondary"
            >
              {chip}
            </span>
          ),
        )}
      </div>
      <div className="flex flex-wrap gap-2">
        <Link
          to={`/ask?gene=${gene.symbol}`}
          className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3.5 py-2 text-sm text-text transition hover:border-accent/40"
        >
          <MessageSquareText className="size-3.5 text-accent" />
          Ask about {gene.symbol}
        </Link>
        <Link
          to={`/generate?gene=${gene.symbol}`}
          className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3.5 py-2 text-sm text-text transition hover:border-accent/40"
        >
          <FileText className="size-3.5 text-accent" />
          Generate Dossier
        </Link>
        <Link
          to={`/compare?genes=${gene.symbol},CDH10`}
          className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3.5 py-2 text-sm text-text transition hover:border-accent/40"
        >
          <GitCompare className="size-3.5 text-accent" />
          Compare
        </Link>
        <button
          type="button"
          className="inline-flex items-center gap-2 rounded-xl border border-border bg-card px-3.5 py-2 text-sm text-text-secondary"
        >
          <RefreshCw className="size-3.5" />
          Refresh Evidence
        </button>
      </div>
    </div>
  )
}
