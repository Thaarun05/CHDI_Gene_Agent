import { useEffect, useMemo, useState } from 'react'
import { listAllEvidence } from '@/api/client'
import type { EvidenceRecord } from '@/api/types'
import { StatusBadge } from '@/components/StatusBadge'
import { useEvidenceDrawer } from '@/context/EvidenceDrawerContext'
import { ErrorState, LoadingSkeleton } from '@/components/EmptyState'

export function EvidencePage() {
  const { openEvidence } = useEvidenceDrawer()
  const [rows, setRows] = useState<EvidenceRecord[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [gene, setGene] = useState('')
  const [source, setSource] = useState('')
  const [type, setType] = useState('')
  const [section, setSection] = useState('')

  useEffect(() => {
    setLoading(true)
    listAllEvidence({
      gene: gene || undefined,
      source: source || undefined,
      type: type || undefined,
      section: section || undefined,
    })
      .then(setRows)
      .catch((e: Error) => setError(e.message))
      .finally(() => setLoading(false))
  }, [gene, source, type, section])

  const sources = useMemo(
    () => Array.from(new Set(rows.map((r) => r.sourceName))).sort(),
    [rows],
  )
  const types = useMemo(
    () => Array.from(new Set(rows.map((r) => r.evidenceType))).sort(),
    [rows],
  )

  return (
    <div className="fade-in space-y-6">
      <div>
        <h1 className="text-2xl font-medium tracking-tight text-text">Evidence</h1>
        <p className="mt-2 text-sm text-text-secondary">
          Searchable EvidenceRecords (demo/mock). Click a row to open the provenance drawer.
        </p>
      </div>

      <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
        <FilterSelect
          label="Gene"
          value={gene}
          onChange={setGene}
          options={['', 'SREBF2', 'CDH10']}
        />
        <FilterSelect
          label="Source"
          value={source}
          onChange={setSource}
          options={['', ...sources]}
        />
        <FilterSelect label="Evidence type" value={type} onChange={setType} options={['', ...types]} />
        <label className="block">
          <span className="mb-1.5 block text-xs text-text-muted">Section</span>
          <input
            value={section}
            onChange={(e) => setSection(e.target.value)}
            placeholder="Filter section…"
            className="w-full rounded-xl border border-border bg-card px-3 py-2.5 text-sm text-text outline-none"
          />
        </label>
      </div>

      {loading && <LoadingSkeleton rows={6} />}
      {error && <ErrorState description={error} />}

      {!loading && !error && (
        <div className="surface-card overflow-hidden">
          <div className="overflow-x-auto">
            <table className="w-full min-w-[800px] border-collapse text-left text-sm">
              <thead>
                <tr className="border-b border-border bg-bg-secondary/50">
                  {['Gene', 'Evidence', 'Source', 'Type', 'Section', 'Retrieved', 'Status'].map(
                    (h) => (
                      <th key={h} className="px-4 py-3 font-medium text-text-secondary">
                        {h}
                      </th>
                    ),
                  )}
                </tr>
              </thead>
              <tbody>
                {rows.map((r) => (
                  <tr
                    key={r.id}
                    className="cursor-pointer border-b border-border/70 last:border-0 hover:bg-white/[0.02]"
                    onClick={() => void openEvidence(r)}
                  >
                    <td className="px-4 py-3 text-text">{r.geneSymbol}</td>
                    <td className="max-w-xs truncate px-4 py-3 text-text-secondary">
                      {r.displayText}
                    </td>
                    <td className="px-4 py-3 text-text-secondary">{r.sourceName}</td>
                    <td className="px-4 py-3 text-text-muted">{r.evidenceType}</td>
                    <td className="px-4 py-3 text-text-muted">{r.section}</td>
                    <td className="px-4 py-3 text-text-muted">
                      {new Date(r.retrievedAt).toLocaleDateString()}
                    </td>
                    <td className="px-4 py-3">
                      {r.status ? <StatusBadge status={r.status} /> : '—'}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  )
}

function FilterSelect({
  label,
  value,
  onChange,
  options,
}: {
  label: string
  value: string
  onChange: (v: string) => void
  options: string[]
}) {
  return (
    <label className="block">
      <span className="mb-1.5 block text-xs text-text-muted">{label}</span>
      <select
        value={value}
        onChange={(e) => onChange(e.target.value)}
        className="w-full rounded-xl border border-border bg-card px-3 py-2.5 text-sm text-text outline-none"
      >
        {options.map((o) => (
          <option key={o || 'all'} value={o}>
            {o || 'All'}
          </option>
        ))}
      </select>
    </label>
  )
}
