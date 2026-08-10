import { X, ExternalLink, Copy, FileJson } from 'lucide-react'
import { useEvidenceDrawer } from '@/context/EvidenceDrawerContext'
import { cn } from '@/lib/utils'

function MetaRow({ label, value }: { label: string; value?: string }) {
  return (
    <div className="flex items-start justify-between gap-4 border-b border-border/60 py-2.5 last:border-0">
      <dt className="text-xs text-text-muted">{label}</dt>
      <dd className="max-w-[60%] text-right text-sm text-text break-all">{value || '—'}</dd>
    </div>
  )
}

export function EvidenceDrawer() {
  const { open, record, loading, closeEvidence } = useEvidenceDrawer()

  return (
    <>
      <div
        className={cn(
          'fixed inset-0 z-40 bg-black/50 transition-opacity duration-200',
          open ? 'opacity-100' : 'pointer-events-none opacity-0',
        )}
        onClick={closeEvidence}
      />
      <aside
        className={cn(
          'fixed top-0 right-0 z-50 flex h-full w-full max-w-md flex-col border-l border-border bg-sidebar shadow-2xl transition-transform duration-200',
          open ? 'translate-x-0' : 'translate-x-full',
        )}
      >
        <div className="flex items-center justify-between border-b border-border px-5 py-4">
          <div>
            <p className="text-xs tracking-wide text-text-muted uppercase">Evidence Record</p>
            <h2 className="mt-1 text-lg font-medium text-text">Provenance drawer</h2>
          </div>
          <button
            type="button"
            onClick={closeEvidence}
            className="rounded-lg p-2 text-text-secondary transition hover:bg-white/5 hover:text-text"
            aria-label="Close"
          >
            <X className="size-4" />
          </button>
        </div>

        <div className="flex-1 overflow-y-auto px-5 py-4">
          {loading && (
            <div className="space-y-3">
              {Array.from({ length: 6 }).map((_, i) => (
                <div key={i} className="h-8 animate-pulse rounded bg-white/5" />
              ))}
            </div>
          )}

          {!loading && record && (
            <div className="fade-in space-y-6">
              <p className="rounded-xl border border-border bg-card p-4 text-sm leading-relaxed text-text-secondary">
                {record.displayText}
              </p>

              <section>
                <h3 className="mb-2 text-xs font-medium tracking-wide text-accent uppercase">
                  Evidence
                </h3>
                <dl>
                  <MetaRow label="Gene" value={record.geneSymbol} />
                  <MetaRow label="Source" value={record.sourceName} />
                  <MetaRow label="Evidence type" value={record.evidenceType} />
                  <MetaRow label="Fact type" value={record.factType} />
                  <MetaRow label="Evidence class" value={record.evidenceClass} />
                  <MetaRow label="Source identifier" value={record.sourceIdentifier} />
                  <MetaRow
                    label="Retrieved date"
                    value={new Date(record.retrievedAt).toLocaleString()}
                  />
                </dl>
              </section>

              <section>
                <h3 className="mb-3 text-xs font-medium tracking-wide text-accent-secondary uppercase">
                  Provenance
                </h3>
                <ol className="space-y-2">
                  {['API Request', 'Raw Artifact', 'Evidence Record', 'Current Answer / Report'].map(
                    (step, i) => (
                      <li
                        key={step}
                        className="flex items-center gap-3 rounded-lg border border-border bg-card px-3 py-2 text-sm"
                      >
                        <span className="flex size-6 items-center justify-center rounded-full bg-accent/15 text-xs text-accent">
                          {i + 1}
                        </span>
                        <span className="text-text-secondary">{step}</span>
                      </li>
                    ),
                  )}
                </ol>
                <dl className="mt-4">
                  <MetaRow label="EvidenceRecord ID" value={record.id} />
                  <MetaRow label="ApiRun ID" value={record.apiRunId} />
                  <MetaRow label="RawArtifact ID" value={record.rawArtifactId} />
                </dl>
              </section>
            </div>
          )}

          {!loading && !record && (
            <p className="text-sm text-text-muted">No evidence record selected.</p>
          )}
        </div>

        <div className="flex gap-2 border-t border-border p-4">
          {record?.sourceUrl && (
            <a
              href={record.sourceUrl}
              target="_blank"
              rel="noreferrer"
              className="inline-flex flex-1 items-center justify-center gap-2 rounded-xl border border-border bg-card px-3 py-2.5 text-sm text-text transition hover:border-accent/40"
            >
              <ExternalLink className="size-3.5" />
              View Source
            </a>
          )}
          <button
            type="button"
            className="inline-flex flex-1 items-center justify-center gap-2 rounded-xl border border-border bg-card px-3 py-2.5 text-sm text-text transition hover:border-accent/40"
            onClick={() => {
              if (record) void navigator.clipboard.writeText(record.id)
            }}
          >
            <Copy className="size-3.5" />
            Copy Citation
          </button>
          <button
            type="button"
            className="inline-flex items-center justify-center gap-2 rounded-xl border border-border bg-card px-3 py-2.5 text-sm text-text-secondary"
            title="Raw evidence viewer connects when FastAPI is ready"
          >
            <FileJson className="size-3.5" />
          </button>
        </div>
      </aside>
    </>
  )
}
