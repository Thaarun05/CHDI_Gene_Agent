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
        allowToolAcquisition: true,
        evidenceSelection: 'accepted_or_latest_generated',
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
              `Using ${selectedGene} as context`,
              'Planning scientific evidence needs',
              'Retrieving and assessing evidence',
              'Building a provenance-grounded answer',
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
              {response.status && (
                <span className="rounded-full border border-border px-2.5 py-1">
                  {response.status.replaceAll('_', ' ')}
                </span>
              )}
              {response.plannerMethod && (
                <span className="rounded-full border border-border px-2.5 py-1">
                  Plan: {response.plannerMethod.replaceAll('_', ' ')}
                </span>
              )}
              <span className="rounded-full border border-border px-2.5 py-1">
                {retrievalLabel(response.retrievalMethod)}
              </span>
              <span className="rounded-full border border-border px-2.5 py-1">
                {generationLabel(response.generationMethod)}
              </span>
              <span className="rounded-full border border-border px-2.5 py-1">
                Embedding: {response.embeddingBackend}
              </span>
            </section>

            <ResolvedEntities response={response} />

            <section>
              <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
                Summary
              </h2>
              <p className="mt-2 text-[15px] leading-relaxed text-text">
                <InlineCitationSummary response={response} />
              </p>
            </section>

            <RequirementSummary response={response} />

            <HdModifierMatrix response={response} />

            <FailureSummary response={response} />

            <EvidenceCategorySummary response={response} />

            {!!response.evidenceGaps?.length && (
              <section>
                <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
                  Evidence gaps
                </h2>
                <ul className="mt-2 list-disc space-y-1 pl-5 text-sm text-text-secondary">
                  {response.evidenceGaps.map((gap) => (
                    <li key={gap}>{gap}</li>
                  ))}
                </ul>
              </section>
            )}

            <RecommendationSummary response={response} />

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
                      <p>Gene: {String(tool.geneSymbol ?? 'unknown')}</p>
                      <p>Sections: {Array.isArray(tool.sectionKeys) ? tool.sectionKeys.join(', ') : 'none'}</p>
                      <p>{tool.reused ? 'Reused evidence run' : 'New evidence run'}: {String(tool.dossierRunId ?? 'none')}</p>
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

function retrievalLabel(method: string) {
  if (method === 'semantic') return 'Semantic Retrieval'
  if (method === 'keyword') return 'Keyword Retrieval'
  if (method === 'metadata') return 'Metadata Retrieval'
  return 'No Retrieval'
}

function generationLabel(method: string) {
  if (method === 'grounded_llm') return 'Grounded LLM'
  if (method === 'hybrid') return 'Hybrid Grounded Synthesis'
  if (method === 'deterministic') return 'Deterministic Grounded Summary'
  return 'No Synthesis'
}

function ResolvedEntities({ response }: { response: AskResponse }) {
  const groups = response.resolvedEntities
    ? Object.entries(response.resolvedEntities).filter(([, values]) => values?.length)
    : []
  if (!groups.length) return null

  return (
    <section>
      <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
        Resolved entities
      </h2>
      <div className="mt-2 flex flex-wrap gap-2">
        {groups.flatMap(([kind, values]) =>
          (values ?? []).map((value) => (
            <span key={`${kind}-${value}`} className="rounded-full border border-border px-2.5 py-1 text-xs text-text-secondary">
              {value}
            </span>
          )),
        )}
      </div>
    </section>
  )
}

function RequirementSummary({ response }: { response: AskResponse }) {
  const assessments = response.requirementAssessments ?? []
  if (!assessments.length) return null

  return (
    <section>
      <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
        Evidence requirements
      </h2>
      <div className="mt-3 overflow-x-auto">
        <table className="w-full min-w-[540px] text-left text-xs">
          <thead className="text-text-muted">
            <tr className="border-b border-border">
              <th className="px-2 py-2 font-medium">Gene</th>
              <th className="px-2 py-2 font-medium">Requirement</th>
              <th className="px-2 py-2 font-medium">Role</th>
              <th className="px-2 py-2 font-medium">Status</th>
              <th className="px-2 py-2 font-medium">Support</th>
            </tr>
          </thead>
          <tbody>
            {assessments.map((item) => (
              <tr key={`${item.requirement_id}-${item.gene_symbol}`} className="border-b border-border/60 text-text-secondary">
                <td className="px-2 py-2 text-text">{item.gene_symbol}</td>
                <td className="px-2 py-2">{item.evidence_need.replaceAll('_', ' ')}</td>
                <td className="px-2 py-2">{item.required ? 'Required' : 'Supporting'}</td>
                <td className="px-2 py-2">{item.status.replaceAll('_', ' ')}</td>
                <td className="px-2 py-2">{item.qualifying_count}/{item.minimum_support}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  )
}

function HdModifierMatrix({ response }: { response: AskResponse }) {
  const rows = response.comparisonMatrix ?? []
  const genes = response.resolvedEntities?.genes ?? []
  if (!rows.length || !genes.length) return null

  return (
    <section>
      <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
        HD modifier evidence matrix
      </h2>
      <div className="mt-3 overflow-x-auto">
        <table className="w-full min-w-[680px] text-left text-xs">
          <thead className="text-text-muted">
            <tr className="border-b border-border">
              <th className="px-2 py-2 font-medium">Dimension</th>
              {genes.map((gene) => <th key={gene} className="px-2 py-2 font-medium">{gene}</th>)}
            </tr>
          </thead>
          <tbody>
            {rows.map((row) => (
              <tr key={row.dimension} className="border-b border-border/60 align-top">
                <th className="px-2 py-3 font-medium text-text">{row.dimension}</th>
                {genes.map((gene) => {
                  const cell = row.cells[gene]
                  return (
                    <td key={gene} className="px-2 py-3 text-text-secondary">
                      <p className="font-medium text-text">{cell?.status ?? 'Missing'}</p>
                      <p className="mt-1">{cell?.summary ?? 'No qualifying evidence.'}</p>
                      {!!cell?.evidenceRecordIds.length && (
                        <div className="mt-2 flex flex-wrap gap-1">
                          {cell.evidenceRecordIds.map((recordId) => (
                            <CitationChip key={recordId} label="Evidence" evidenceRecordId={recordId} />
                          ))}
                        </div>
                      )}
                    </td>
                  )
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  )
}

function FailureSummary({ response }: { response: AskResponse }) {
  const failures = response.failures ?? []
  if (!failures.length) return null

  return (
    <section>
      <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
        Generation status
      </h2>
      <ul className="mt-2 list-disc space-y-1 pl-5 text-sm text-text-secondary">
        {failures.map((failure, index) => (
          <li key={index}>
            {String(failure.failure_type ?? 'failure').replaceAll('_', ' ')}: {String(failure.message ?? 'No detail available')}
          </li>
        ))}
      </ul>
    </section>
  )
}

function EvidenceCategorySummary({ response }: { response: AskResponse }) {
  const categories = response.evidenceCategories ?? []
  if (!categories.length) return null

  return (
    <section>
      <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
        Evidence categories
      </h2>
      <div className="mt-3 overflow-x-auto">
        <table className="w-full min-w-[640px] text-left text-xs">
          <thead className="text-text-muted">
            <tr className="border-b border-border">
              <th className="px-2 py-2 font-medium">Gene</th>
              <th className="px-2 py-2 font-medium">Category</th>
              <th className="px-2 py-2 font-medium">Evidence system</th>
              <th className="px-2 py-2 font-medium">Claim type</th>
              <th className="px-2 py-2 font-medium">Status</th>
            </tr>
          </thead>
          <tbody>
            {categories.map((item, index) => (
              <tr key={index} className="border-b border-border/60 text-text-secondary">
                <td className="px-2 py-2 text-text">{String(item.gene_symbol ?? '')}</td>
                <td className="px-2 py-2">{String(item.category ?? '').replaceAll('_', ' ')}</td>
                <td className="px-2 py-2">{String(item.evidence_system ?? 'not specified')}</td>
                <td className="px-2 py-2">{String(item.claim_type ?? 'not specified')}</td>
                <td className="px-2 py-2">{String(item.status ?? '').replaceAll('_', ' ')}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  )
}

function RecommendationSummary({ response }: { response: AskResponse }) {
  const recommendations = response.recommendations ?? []
  if (!recommendations.length) return null

  return (
    <section>
      <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">
        Recommendations
      </h2>
      <ul className="mt-2 list-disc space-y-2 pl-5 text-sm text-text-secondary">
        {recommendations.map((item, index) => (
          <li key={index}>
            <span className="text-text">{String(item.label ?? 'Recommendation')}:</span>{' '}
            {String(item.description ?? '')}
            {Array.isArray(item.gap_ids) && item.gap_ids.length ? (
              <span className="block text-xs text-text-muted">
                Gaps: {item.gap_ids.map(String).join(', ')}
              </span>
            ) : null}
          </li>
        ))}
      </ul>
    </section>
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
