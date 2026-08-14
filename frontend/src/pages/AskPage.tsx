import { useEffect, useRef, useState } from 'react'
import { AlertCircle, RefreshCw } from 'lucide-react'
import { useSearchParams } from 'react-router-dom'
import { askEvidenceQuestion } from '@/api/client'
import type { AskResponse, PublicEvidenceItem, ResearchMode } from '@/api/types'
import { AgentActivity, AgentActivityIdle } from '@/components/AgentActivity'
import { CitationChip } from '@/components/CitationChip'
import { LoadingSkeleton } from '@/components/EmptyState'
import { SearchComposer } from '@/components/SearchComposer'
import { TopBar } from '@/components/TopBar'

const CONTEXT_GENE_RE = /^[A-Z0-9][A-Z0-9-]{1,14}$/

function parseContextGene(value: string | null): string | null {
  const normalized = value?.trim().toUpperCase() ?? ''
  return CONTEXT_GENE_RE.test(normalized) ? normalized : null
}

function parseResearchMode(value: string | null): ResearchMode {
  if (value === 'deep_research' || value === 'stored_only') return value
  return 'auto'
}

function humanize(value: string) {
  return value.replaceAll('_', ' ').replace(/\b\w/g, (character) => character.toUpperCase())
}

type RequestError = {
  kind: 'cancelled' | 'backend_unavailable' | 'provider_failure' | 'request_failed'
  message: string
}

export function AskPage() {
  const [params, setParams] = useSearchParams()
  const initialQuestion = params.get('q') ?? ''
  const initialContext = parseContextGene(params.get('gene'))
  const initialMode = parseResearchMode(params.get('research_mode'))
  const refreshIfAvailable =
    params.get('refresh_if_available') === 'true' || params.get('refresh') === 'true'

  const [query, setQuery] = useState(initialQuestion)
  const [contextGene, setContextGene] = useState<string | null>(initialContext)
  const [researchMode, setResearchMode] = useState<ResearchMode>(initialMode)
  const [loading, setLoading] = useState(false)
  const [response, setResponse] = useState<AskResponse | null>(null)
  const [error, setError] = useState<RequestError | null>(null)
  const [detail, setDetail] = useState<'evidence' | 'contextual' | 'sources' | 'tools' | null>(null)
  const requestGeneration = useRef(0)
  const activeController = useRef<AbortController | null>(null)
  const cancelledByUser = useRef(false)
  const autoSubmitted = useRef(false)

  async function submit(
    question = query,
    selectedContext = contextGene,
    selectedMode = researchMode,
  ) {
    if (!question.trim() || loading) return
    const generation = ++requestGeneration.current
    activeController.current?.abort()
    const controller = new AbortController()
    activeController.current = controller
    cancelledByUser.current = false

    setLoading(true)
    setError(null)
    setResponse(null)
    setDetail(null)
    try {
      const result = await askEvidenceQuestion(question.trim(), selectedContext, {
        refreshIfAvailable,
        allowToolAcquisition: selectedMode !== 'stored_only',
        evidenceSelection: 'accepted_or_latest_generated',
        researchMode: selectedMode,
        signal: controller.signal,
      })
      if (generation !== requestGeneration.current) return
      setResponse(result)
      if (result.failures?.some((failure) => failure.failure_type === 'provider_failure')) {
        setError({
          kind: 'provider_failure',
          message: 'The answer provider failed. A deterministic grounded response is shown where available.',
        })
      }
    } catch (caught) {
      if (generation !== requestGeneration.current) return
      if (controller.signal.aborted) {
        // React StrictMode remounts abort in-flight requests; only surface explicit user cancels.
        if (!cancelledByUser.current) return
        setError({
          kind: 'cancelled',
          message: 'The request was cancelled.',
        })
      } else {
        const message = caught instanceof Error ? caught.message : 'Unknown request failure'
        setError({
          kind: caught instanceof TypeError ? 'backend_unavailable' : 'request_failed',
          message:
            caught instanceof TypeError
              ? 'The scientific backend is unavailable. Check the connection and retry.'
              : `The request failed: ${message}`,
        })
      }
    } finally {
      if (generation === requestGeneration.current) {
        activeController.current = null
        setLoading(false)
      }
    }
  }

  function cancelRequest() {
    cancelledByUser.current = true
    activeController.current?.abort()
  }

  function selectContextGene(gene: string | null) {
    requestGeneration.current += 1
    activeController.current?.abort()
    activeController.current = null
    setContextGene(gene)
    setResponse(null)
    setError(null)
    setDetail(null)
    setLoading(false)
    setParams((current) => {
      const next = new URLSearchParams(current)
      if (gene) next.set('gene', gene)
      else next.delete('gene')
      return next
    }, { replace: true })
  }

  function selectResearchMode(mode: ResearchMode) {
    setResearchMode(mode)
    setParams((current) => {
      const next = new URLSearchParams(current)
      next.set('research_mode', mode)
      return next
    }, { replace: true })
  }

  useEffect(() => {
    if (autoSubmitted.current || !initialQuestion.trim()) return
    autoSubmitted.current = true
    void submit(initialQuestion, initialContext, initialMode)
    // Route state is intentionally captured once for initial auto-submission.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  useEffect(() => () => activeController.current?.abort(), [])

  return (
    <div className="fade-in mx-auto flex max-w-4xl flex-col gap-8">
      <TopBar />

      <div className="text-center">
        <p className="text-sm text-accent-secondary">Ask Evidence Question</p>
        <h1 className="mt-2 text-2xl font-medium tracking-tight text-text">
          Evidence-grounded scientific Q&amp;A
        </h1>
        <p className="mt-2 text-sm text-text-secondary">
          Answers cite persisted scientific evidence. The backend validates every supporting record.
        </p>
      </div>

      <SearchComposer
        value={query}
        onChange={setQuery}
        onSubmit={() => void submit()}
        selectedGene={contextGene}
        onSelectGene={selectContextGene}
        researchMode={researchMode}
        onSelectResearchMode={selectResearchMode}
        submitting={loading}
        onCancel={cancelRequest}
      />

      {loading && (
        <div className="space-y-4">
          <AgentActivity
            loading
            steps={[
              'Planning requirements',
              'Searching stored evidence',
              'Deciding whether acquisition is needed',
              'Executing approved tools',
              'Validating retrieved evidence',
              'Generating the grounded answer',
            ]}
          />
          <p className="text-center text-xs text-text-muted">
            This can take several minutes. You can cancel anytime; the request stops automatically after 6 minutes.
          </p>
          <LoadingSkeleton rows={4} />
        </div>
      )}

      {!loading && error && (
        <div className="flex items-start justify-between gap-4 rounded-lg border border-danger/35 bg-danger/5 p-4">
          <div className="flex gap-3">
            <AlertCircle className="mt-0.5 size-4 shrink-0 text-danger" />
            <div>
              <p className="text-sm font-medium text-text">{humanize(error.kind)}</p>
              <p className="mt-1 text-sm text-text-secondary">{error.message}</p>
            </div>
          </div>
          {error.kind !== 'cancelled' && (
            <button
              type="button"
              onClick={() => void submit()}
              className="inline-flex shrink-0 items-center gap-1.5 rounded-lg border border-border px-3 py-2 text-xs text-text-secondary hover:text-text"
            >
              <RefreshCw className="size-3.5" />
              Retry
            </button>
          )}
        </div>
      )}

      {!loading && !response && !error && <AgentActivityIdle />}

      {!loading && response && (
        <div className="space-y-5">
          <AgentActivity steps={response.agentActivity} />

          <article className="surface-card space-y-7 p-6">
            <ResponseBadges response={response} researchMode={researchMode} />
            <ResolvedEntities response={response} />
            <StructuredAnswer response={response} />
            <RequirementSummary response={response} />
            <ComparisonDecisionSummary response={response} />
            <HdModifierMatrix response={response} onViewEvidence={() => setDetail('evidence')} />
            <EvidenceGapSummary response={response} />
            <RecommendationSummary response={response} />
            <Limitations response={response} />
            <TechnicalDiagnostics response={response} />
          </article>

          <div className="flex flex-wrap gap-2">
            <MetaButton
              active={detail === 'evidence'}
              onClick={() => setDetail(detail === 'evidence' ? null : 'evidence')}
              label={`Qualifying evidence · ${response.activitySummary?.accepted_evidence ?? response.evidenceUsedCount}`}
            />
            <MetaButton
              active={detail === 'contextual'}
              onClick={() => setDetail(detail === 'contextual' ? null : 'contextual')}
              label={`Contextual / excluded · ${response.activitySummary?.rejected_evidence ?? 0}`}
            />
            <MetaButton
              active={detail === 'sources'}
              onClick={() => setDetail(detail === 'sources' ? null : 'sources')}
              label={`Sources · ${response.sourcesCount}`}
            />
            <MetaButton
              active={detail === 'tools'}
              onClick={() => setDetail(detail === 'tools' ? null : 'tools')}
              label={`Tools executed · ${response.activitySummary?.tools_executed ?? response.toolsInvokedCount}`}
            />
          </div>

          {detail && (
            <div className="surface-elevated p-4 text-sm text-text-secondary">
              {detail === 'evidence' && <EvidenceList items={response.evidenceItems ?? []} />}
              {detail === 'contextual' && (
                <EvidenceList items={response.contextualEvidence ?? []} contextual />
              )}
              {detail === 'sources' && (
                <p>{response.sourcesUsed.length ? response.sourcesUsed.join(', ') : 'No sources used.'}</p>
              )}
              {detail === 'tools' && <ToolSummary response={response} />}
            </div>
          )}
        </div>
      )}
    </div>
  )
}

function ResponseBadges({ response, researchMode }: { response: AskResponse; researchMode: ResearchMode }) {
  return (
    <section className="flex flex-wrap gap-2 text-xs text-text-secondary">
      <span className="rounded-full border border-border px-2.5 py-1">{humanize(response.status)}</span>
      {response.plannerMethod && (
        <span className="rounded-full border border-border px-2.5 py-1">
          Plan: {humanize(response.plannerMethod)}
        </span>
      )}
      <span className="rounded-full border border-border px-2.5 py-1">
        {retrievalLabel(response.retrievalMethod)}
      </span>
      <span className="rounded-full border border-border px-2.5 py-1">
        {generationLabel(response.generationMethod)}
      </span>
      <span className="rounded-full border border-border px-2.5 py-1">
        {researchMode === 'stored_only' ? 'Stored Evidence Only' : researchMode === 'deep_research' ? 'Deep Research' : 'Auto'}
      </span>
    </section>
  )
}

function StructuredAnswer({ response }: { response: AskResponse }) {
  const sections = response.answerSections ?? []
  if (!sections.length) {
    return (
      <section>
        <h2 className="text-sm font-medium text-text">Direct answer</h2>
        <p className="mt-2 text-[15px] leading-relaxed text-text-secondary">
          <InlineCitations text={response.summary} response={response} />
        </p>
      </section>
    )
  }
  return (
    <div className="space-y-6">
      {sections.map((section) => (
        <section key={section.key}>
          <h2 className="text-sm font-medium text-text">{section.title}</h2>
          <div className="mt-2 space-y-2 text-[15px] leading-relaxed text-text-secondary">
            {section.paragraphs.map((paragraph, index) => (
              <p key={`${section.key}-${index}`}>
                <InlineCitations text={paragraph} response={response} />
              </p>
            ))}
          </div>
        </section>
      ))}
    </div>
  )
}

function ResolvedEntities({ response }: { response: AskResponse }) {
  const groups = response.resolvedEntities
    ? Object.entries(response.resolvedEntities).filter(([, values]) => values?.length)
    : []
  if (!groups.length) return null
  return (
    <section>
      <h2 className="text-xs font-medium tracking-wide text-text-muted uppercase">Resolved entities</h2>
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
      <h2 className="text-sm font-medium text-text">Evidence requirements</h2>
      <div className="mt-3 overflow-x-auto">
        <table className="w-full min-w-[620px] text-left text-xs">
          <thead className="text-text-muted">
            <tr className="border-b border-border">
              <th className="px-2 py-2 font-medium">Gene</th>
              <th className="px-2 py-2 font-medium">Requirement</th>
              <th className="px-2 py-2 font-medium">Role</th>
              <th className="px-2 py-2 font-medium">Status</th>
              <th className="px-2 py-2 font-medium">Unique support</th>
              <th className="px-2 py-2 font-medium">Sources</th>
            </tr>
          </thead>
          <tbody>
            {assessments.map((item) => (
              <tr key={`${item.requirement_id}-${item.gene_symbol}`} className="border-b border-border/60 text-text-secondary">
                <td className="px-2 py-2 text-text">{item.gene_symbol}</td>
                <td className="px-2 py-2">{humanize(item.evidence_need)}</td>
                <td className="px-2 py-2">{item.required ? 'Required' : 'Supporting'}</td>
                <td className="px-2 py-2">{humanize(item.status)}</td>
                <td className="px-2 py-2">{item.qualifying_count}/{item.minimum_support}</td>
                <td className="px-2 py-2">{item.distinct_source_count ?? 0}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  )
}

function ComparisonDecisionSummary({ response }: { response: AskResponse }) {
  const decision = response.comparisonDecision
  if (!decision) return null
  return (
    <section>
      <h2 className="text-sm font-medium text-text">Comparison conclusion</h2>
      <p className="mt-2 text-sm leading-relaxed text-text-secondary">{decision.summary}</p>
      {!!decision.limitations.length && (
        <ul className="mt-2 list-disc space-y-1 pl-5 text-xs text-text-muted">
          {decision.limitations.map((limitation) => <li key={limitation}>{limitation}</li>)}
        </ul>
      )}
    </section>
  )
}

function HdModifierMatrix({
  response,
  onViewEvidence,
}: {
  response: AskResponse
  onViewEvidence: () => void
}) {
  const rows = response.comparisonMatrix ?? []
  const genes = response.resolvedEntities?.genes ?? []
  if (!rows.length || !genes.length) return null
  return (
    <section>
      <h2 className="text-sm font-medium text-text">HD modifier evidence matrix</h2>
      <div className="mt-3 overflow-x-auto">
        <table className="w-full min-w-[720px] text-left text-xs">
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
                      <p className="mt-1">{cell?.evidenceCount ?? 0} unique item(s) · {cell?.distinctSourceCount ?? 0} source(s)</p>
                      <p className="mt-1">{cell?.summary ?? 'No qualifying evidence.'}</p>
                      {!!cell?.citationOrdinals?.length && (
                        <div className="mt-2 flex flex-wrap gap-1">
                          {cell.citationOrdinals.slice(0, 3).map((ordinal) => {
                            const citation = response.citations[ordinal - 1]
                            return citation ? (
                              <CitationChip
                                key={`${gene}-${row.dimension}-${ordinal}`}
                                label={String(ordinal)}
                                evidenceReference={citation.publicEvidenceRef}
                              />
                            ) : null
                          })}
                        </div>
                      )}
                      {(cell?.evidenceCount ?? 0) > 3 && (
                        <button type="button" onClick={onViewEvidence} className="mt-2 text-xs text-accent hover:underline">
                          View all {cell.evidenceCount}
                        </button>
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

function EvidenceGapSummary({ response }: { response: AskResponse }) {
  if (!response.evidenceGaps?.length) return null
  return (
    <section>
      <h2 className="text-sm font-medium text-text">Evidence gaps</h2>
      <ul className="mt-2 list-disc space-y-1 pl-5 text-sm text-text-secondary">
        {response.evidenceGaps.map((gap) => <li key={gap}>{humanize(gap)}</li>)}
      </ul>
    </section>
  )
}

function RecommendationSummary({ response }: { response: AskResponse }) {
  const recommendations = response.recommendations ?? []
  if (!recommendations.length) return null
  return (
    <section>
      <h2 className="text-sm font-medium text-text">Recommended next studies</h2>
      <ul className="mt-2 space-y-3 text-sm text-text-secondary">
        {recommendations.map((item, index) => (
          <li key={index} className="border-l-2 border-accent/30 pl-3">
            <span className="text-text">{String(item.label ?? 'Recommendation')}:</span>{' '}
            {String(item.description ?? '')}
            {Array.isArray(item.gap_labels) && item.gap_labels.length ? (
              <span className="mt-1 block text-xs text-text-muted">
                Addresses: {item.gap_labels.map(String).join(', ')}
              </span>
            ) : null}
          </li>
        ))}
      </ul>
    </section>
  )
}

function Limitations({ response }: { response: AskResponse }) {
  if (!response.limitations.length) return null
  return (
    <section>
      <h2 className="text-sm font-medium text-text">Limitations</h2>
      <ul className="mt-2 list-disc space-y-1 pl-5 text-sm text-text-secondary">
        {response.limitations.map((limitation) => <li key={limitation}>{limitation}</li>)}
      </ul>
    </section>
  )
}

function TechnicalDiagnostics({ response }: { response: AskResponse }) {
  return (
    <details className="border-t border-border pt-4">
      <summary className="cursor-pointer text-xs font-medium tracking-wide text-text-muted uppercase">
        Technical diagnostics
      </summary>
      <div className="mt-3 grid gap-2 text-xs text-text-secondary sm:grid-cols-2">
        <p>Embedding: {response.embeddingBackend}</p>
        <p>Evidence universe: {humanize(response.evidenceUniverse)}</p>
        <p>Accepted evidence: {response.activitySummary?.accepted_evidence ?? response.evidenceUsedCount}</p>
        <p>Rejected evidence: {response.activitySummary?.rejected_evidence ?? 0}</p>
        <p>Failed tool executions: {response.activitySummary?.tools_failed ?? 0}</p>
        <p>Reused runs: {response.activitySummary?.runs_reused ?? 0}</p>
        <p>Tools skipped: {response.activitySummary?.tools_skipped ?? 0}</p>
      </div>
    </details>
  )
}

function EvidenceList({ items, contextual = false }: { items: PublicEvidenceItem[]; contextual?: boolean }) {
  if (!items.length) return <p>No {contextual ? 'contextual or excluded' : 'qualifying'} evidence items.</p>
  return (
    <ul className="space-y-3">
      {items.map((item) => (
        <li key={`${item.public_evidence_ref}-${item.exclusion_reason ?? 'accepted'}`} className="border-b border-border/60 pb-3 last:border-0">
          <div className="flex flex-wrap items-center gap-2">
            <CitationChip label="Evidence" evidenceReference={item.public_evidence_ref} />
            <span className="text-text">{item.title || item.public_identifier || item.source_name}</span>
            <span className="text-xs text-text-muted">{item.source_name}</span>
          </div>
          <p className="mt-1 text-sm leading-relaxed">{item.display_text}</p>
          <p className="mt-1 text-xs text-text-muted">
            {humanize(item.evidence_need)} · {humanize(item.designation)}
            {item.public_identifier ? ` · ${item.public_identifier}` : ''}
            {item.backing_record_count > 1 ? ` · ${item.backing_record_count} deduplicated rows` : ''}
          </p>
          {item.exclusion_reason && (
            <p className="mt-1 text-xs text-warning">Excluded from support: {humanize(item.exclusion_reason)}</p>
          )}
        </li>
      ))}
    </ul>
  )
}

function ToolSummary({ response }: { response: AskResponse }) {
  if (!response.toolActivity.length) return <p>No tools executed or reused.</p>
  return (
    <div className="space-y-3">
      {response.toolActivity.map((tool, index) => (
        <div key={`${String(tool.publicRunRef ?? index)}-${index}`} className="border-b border-border/60 pb-3 last:border-0">
          <p className="text-text">{tool.reused ? 'Reused run' : 'Executed tool'}: {String(tool.toolName ?? 'Approved capability')}</p>
          <p>Gene: {String(tool.geneSymbol ?? 'Not available')}</p>
          <p>Execution: {tool.executionSucceeded ? 'Succeeded' : humanize(String(tool.status ?? 'unknown'))}</p>
          <p>Accepted evidence: {String(tool.qualifyingEvidenceCount ?? 0)}</p>
          <p>Rejected evidence: {String(tool.rejectedEvidenceCount ?? 0)}</p>
          <p>Scientific retrieval: {tool.scientificRetrievalSucceeded ? 'Qualifying evidence found' : 'No qualifying evidence found'}</p>
        </div>
      ))}
    </div>
  )
}

function InlineCitations({ text, response }: { text: string; response: AskResponse }) {
  return text.split(/(\[\[\d+\]\])/g).map((part, index) => {
    const marker = /^\[\[(\d+)\]\]$/.exec(part)
    if (!marker) return part
    const ordinal = Number(marker[1])
    const citation = response.citations[ordinal - 1]
    return citation ? (
      <CitationChip
        key={`${ordinal}-${index}`}
        label={String(ordinal)}
        evidenceReference={citation.publicEvidenceRef}
        className="mx-0.5"
      />
    ) : part
  })
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

function MetaButton({ label, onClick, active }: { label: string; onClick: () => void; active: boolean }) {
  return (
    <button
      type="button"
      onClick={onClick}
      className={
        active
          ? 'rounded-lg border border-accent/40 bg-accent/10 px-3.5 py-2 text-sm text-accent'
          : 'rounded-lg border border-border bg-card px-3.5 py-2 text-sm text-text-secondary transition hover:text-text'
      }
    >
      {label}
    </button>
  )
}
