import { useEffect, useMemo, useRef, useState } from 'react'
import { useSearchParams } from 'react-router-dom'
import { askEvidenceQuestion } from '@/api/client'
import type { AskResponse } from '@/api/types'
import { SearchComposer } from '@/components/SearchComposer'
import { AgentActivity, AgentActivityIdle } from '@/components/AgentActivity'
import { CitationChip } from '@/components/CitationChip'
import { LoadingSkeleton } from '@/components/EmptyState'
import { TopBar } from '@/components/TopBar'

type AskGene = 'SREBF2' | 'CDH10'

function parseAskGene(value: string | null): AskGene {
  return value?.trim().toUpperCase() === 'CDH10' ? 'CDH10' : 'SREBF2'
}

export function AskPage() {
  const [params, setParams] = useSearchParams()
  const initialQ = params.get('q') ?? ''
  const initialGene = parseAskGene(params.get('gene'))
  const refreshIfAvailable =
    params.get('refresh_if_available') === 'true' || params.get('refresh') === 'true'

  const [query, setQuery] = useState(initialQ)
  const [selectedGene, setSelectedGene] = useState<AskGene>(initialGene)
  const [loading, setLoading] = useState(false)
  const [response, setResponse] = useState<AskResponse | null>(null)
  const [detail, setDetail] = useState<'evidence' | 'sources' | 'tools' | null>(null)
  const requestGeneration = useRef(0)

  const citationMap = useMemo(() => {
    const m = new Map<string, AskResponse['citations'][number]>()
    response?.citations.forEach((c) => m.set(c.id, c))
    return m
  }, [response])

  async function submit(q = query) {
    const generation = ++requestGeneration.current
    const requestGene = selectedGene
    setLoading(true)
    setDetail(null)
    try {
      const res = await askEvidenceQuestion(q, requestGene, {
        refreshIfAvailable,
      })
      if (generation === requestGeneration.current) {
        setResponse(res)
      }
    } finally {
      if (generation === requestGeneration.current) {
        setLoading(false)
      }
    }
  }

  function selectGene(value: string) {
    const gene = parseAskGene(value)
    requestGeneration.current += 1
    setSelectedGene(gene)
    setResponse(null)
    setDetail(null)
    setLoading(false)
    setParams((current) => {
      const next = new URLSearchParams(current)
      next.set('gene', gene)
      return next
    }, { replace: true })
  }

  useEffect(() => {
    if (initialQ.trim()) {
      void submit(initialQ)
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  return (
    <div className="fade-in mx-auto flex max-w-3xl flex-col gap-8">
      <TopBar />

      <div className="text-center">
        <p className="text-sm text-accent-secondary">Ask Evidence Question</p>
        <h1 className="mt-2 text-2xl font-medium tracking-tight text-text">
          Evidence-grounded scientific Q&A
        </h1>
        <p className="mt-2 text-sm text-text-secondary">
          Answers cite stored EvidenceRecords. The LLM is not the scientific source of truth.
        </p>
      </div>

      <SearchComposer
        value={query}
        onChange={setQuery}
        onSubmit={() => void submit()}
        selectedGene={selectedGene}
        onSelectGene={selectGene}
      />

      {loading && (
        <div className="space-y-4">
          <AgentActivity
            loading={loading}
            steps={[
              `Resolved ${selectedGene}`,
              'Searching stored evidence',
              'Checking chemical-tool evidence',
              'Building grounded answer',
            ]}
          />
          <LoadingSkeleton rows={5} />
        </div>
      )}

      {!loading && !response && <AgentActivityIdle />}

      {!loading && response && (
        <div className="space-y-5">
          <AgentActivity steps={response.agentActivity} />

          <article className="surface-card space-y-6 p-6">
            <section className="flex flex-wrap gap-2 text-xs text-text-secondary">
              <span className="rounded-full border border-border px-2.5 py-1">
                {response.retrievalMethod === 'semantic' ? 'Semantic Retrieval' : 'Keyword Retrieval'}
              </span>
              <span className="rounded-full border border-border px-2.5 py-1">
                {response.generationMethod === 'grounded_llm'
                  ? 'Grounded LLM'
                  : 'Deterministic Grounded Summary'}
              </span>
              <span className="rounded-full border border-border px-2.5 py-1">
                Embedding: {response.embeddingBackend}
              </span>
            </section>

            <section>
              <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
                Summary
              </h2>
              <p className="mt-2 text-[15px] leading-relaxed text-text">
                <InlineCitationSummary response={response} />
              </p>
            </section>

            <section>
              <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
                Evidence
              </h2>
              <div className="mt-3 space-y-4">
                {response.evidenceBlocks.map((block) => (
                  <div key={block.sourceGroup}>
                    <h3 className="text-sm font-medium text-accent-secondary">
                      {block.sourceGroup}
                    </h3>
                    <ul className="mt-2 space-y-2">
                      {block.items.map((item, i) => (
                        <li key={i} className="text-sm leading-relaxed text-text-secondary">
                          {item.text}{' '}
                          {item.citationIds.map((cid) => {
                            const c = citationMap.get(cid)
                            if (!c) return null
                            return (
                              <CitationChip
                                key={cid}
                                label={c.label}
                                evidenceRecordId={c.evidenceRecordId}
                                className="ml-1"
                              />
                            )
                          })}
                        </li>
                      ))}
                    </ul>
                  </div>
                ))}
              </div>
            </section>

            <section>
              <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
                Limitations
              </h2>
              <ul className="mt-2 list-disc space-y-1 pl-5 text-sm text-text-secondary">
                {response.limitations.map((l) => (
                  <li key={l}>{l}</li>
                ))}
              </ul>
            </section>
          </article>

          <div className="flex flex-wrap gap-2">
            <MetaButton
              active={detail === 'evidence'}
              onClick={() => setDetail(detail === 'evidence' ? null : 'evidence')}
              label={`Evidence used · ${response.evidenceUsedCount}`}
            />
            <MetaButton
              active={detail === 'sources'}
              onClick={() => setDetail(detail === 'sources' ? null : 'sources')}
              label={`Sources · ${response.sourcesCount}`}
            />
            <MetaButton
              active={detail === 'tools'}
              onClick={() => setDetail(detail === 'tools' ? null : 'tools')}
              label={`Tools invoked · ${response.toolsInvokedCount}`}
            />
          </div>

          {detail && (
            <div className="surface-elevated p-4 text-sm text-text-secondary">
              {detail === 'evidence' && (
                <ul className="space-y-1">
                  {response.citations.map((c) => (
                    <li key={c.id}>
                      <CitationChip label={c.label} evidenceRecordId={c.evidenceRecordId} /> —{' '}
                      {c.sourceName}
                    </li>
                  ))}
                </ul>
              )}
              {detail === 'sources' && (
                <p>{response.sourcesUsed.length ? response.sourcesUsed.join(', ') : 'No sources used.'}</p>
              )}
              {detail === 'tools' && (
                <div className="space-y-2">
                  {response.toolActivity.length === 0 && <p>No tools invoked.</p>}
                  {response.toolActivity.map((tool, i) => (
                    <div key={i}>
                      <p>Selected tool: {String(tool.toolName ?? 'unknown')}</p>
                      <p>Sections: {Array.isArray(tool.sectionKeys) ? tool.sectionKeys.join(', ') : 'none'}</p>
                      <p>New evidence run: {String(tool.dossierRunId ?? 'none')}</p>
                      <p>Evidence re-indexed: {String(tool.indexedRecords ?? 0)}</p>
                    </div>
                  ))}
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  )
}

function InlineCitationSummary({ response }: { response: AskResponse }) {
  return response.summary.split(/(\[\[\d+\]\])/g).map((part, index) => {
    const marker = /^\[\[(\d+)\]\]$/.exec(part)
    if (!marker) return part

    const ordinal = Number(marker[1])
    const citation = response.citations[ordinal - 1]
    if (!citation) return part

    return (
      <CitationChip
        key={`${ordinal}-${index}`}
        label={`[${ordinal}]`}
        evidenceRecordId={citation.evidenceRecordId}
        className="mx-0.5"
      />
    )
  })
}

function MetaButton({
  label,
  onClick,
  active,
}: {
  label: string
  onClick: () => void
  active: boolean
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      className={
        active
          ? 'rounded-xl border border-accent/40 bg-accent/10 px-3.5 py-2 text-sm text-accent'
          : 'rounded-xl border border-border bg-card px-3.5 py-2 text-sm text-text-secondary transition hover:text-text'
      }
    >
      {label}
    </button>
  )
}
