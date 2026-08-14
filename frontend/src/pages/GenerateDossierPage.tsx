import { useEffect, useState } from 'react'
import { Link, useSearchParams } from 'react-router-dom'
import { getJob, getJobArtifacts, startDossierJob } from '@/api/client'
import type { ReportArtifact, WorkflowJob } from '@/api/types'
import { JobProgress } from '@/components/JobProgress'
import { cn } from '@/lib/utils'

const SECTIONS = [
  'General Gene Information',
  'Structure',
  'Expression',
  'GEO Perturbations',
  'Transcription Factors',
  'Protein Interactions',
  'Chemical Perturbations',
  'Chemical Tools',
]

const ACCEPTED_DEMO_GENES = ['SREBF2', 'CDH10'] as const

const SECTION_KEY_MAP: Record<string, string[]> = {
  'General Gene Information': ['1a', '1b', '1c', '1d', '1e'],
  Structure: ['1b', '1c', '1d'],
  Expression: ['2a', '2b', '2c'],
  'GEO Perturbations': ['3a'],
  'Transcription Factors': ['4a'],
  'Protein Interactions': ['5a', '5b'],
  'Chemical Perturbations': ['6a'],
  'Chemical Tools': ['7a'],
}

function sectionKeysFor(selectedSections: string[]) {
  return Array.from(
    new Set(selectedSections.flatMap((section) => SECTION_KEY_MAP[section] ?? [])),
  )
}

function isAcceptedDemoGene(symbol: string) {
  return (ACCEPTED_DEMO_GENES as readonly string[]).includes(symbol)
}

export function GenerateDossierPage() {
  const [params, setParams] = useSearchParams()
  const initialJobId = params.get('job')
  const initialGene = (params.get('gene') || 'SREBF2').trim().toUpperCase() || 'SREBF2'
  const [gene, setGene] = useState(initialGene)
  const [selected, setSelected] = useState<string[]>([...SECTIONS])
  const [useAccepted, setUseAccepted] = useState(isAcceptedDemoGene(initialGene))
  const acceptedAvailable = isAcceptedDemoGene(gene)
  const [job, setJob] = useState<WorkflowJob | null>(null)
  const [artifact, setArtifact] = useState<ReportArtifact | null>(null)
  const [artifactResolved, setArtifactResolved] = useState(false)
  const [starting, setStarting] = useState(false)
  const [sessionMessage, setSessionMessage] = useState<string | null>(null)

  async function loadArtifacts(jobId: string) {
    setArtifactResolved(false)
    try {
      const artifacts = await getJobArtifacts(jobId)
      setArtifact(artifacts.report)
      if (artifacts.report?.reportOrigin === 'generated') {
        window.dispatchEvent(new Event('generated-report-updated'))
      }
    } catch {
      setArtifact(null)
    } finally {
      setArtifactResolved(true)
    }
  }

  useEffect(() => {
    if (!initialJobId) return
    let active = true
    void getJob(initialJobId)
      .then(async (restoredJob) => {
        if (!active) return
        setJob(restoredJob)
        setSessionMessage(null)
        if (restoredJob.status === 'Completed' || restoredJob.status === 'Partial') {
          await loadArtifacts(restoredJob.id)
        }
      })
      .catch((error: Error) => {
        if (!active) return
        setJob(null)
        setArtifact(null)
        setArtifactResolved(false)
        if (error.message.startsWith('API 404:')) {
          setSessionMessage(
            'This job session has expired. Completed generated dossiers remain available in Reports.',
          )
          const next = new URLSearchParams(params)
          next.delete('job')
          setParams(next, { replace: true })
        } else {
          setSessionMessage(error.message)
        }
      })
    return () => {
      active = false
    }
  }, [initialJobId, params, setParams])

  useEffect(() => {
    if (
      !job ||
      job.status === 'Completed' ||
      job.status === 'Partial' ||
      job.status === 'Failed'
    )
      return
    const t = setInterval(() => {
      void getJob(job.id).then((j) => {
        setJob(j)
        if (j.status === 'Completed' || j.status === 'Partial') {
          void loadArtifacts(j.id)
        }
      })
    }, 700)
    return () => clearInterval(t)
  }, [job])

  async function onGenerate() {
    const normalized = gene.trim().toUpperCase()
    if (!normalized) {
      setSessionMessage('Enter a gene symbol to generate a dossier.')
      return
    }
    if (normalized !== gene) {
      setGene(normalized)
    }
    const useExistingAccepted = isAcceptedDemoGene(normalized) && useAccepted
    setStarting(true)
    setArtifact(null)
    setArtifactResolved(false)
    setSessionMessage(null)
    try {
      const j = await startDossierJob(normalized, {
        sectionKeys: sectionKeysFor(selected),
        useExistingAccepted,
      })
      setJob(j)
      const next = new URLSearchParams(params)
      next.set('gene', normalized)
      next.set('job', j.id)
      setParams(next, { replace: true })
      if (j.status === 'Completed' || j.status === 'Partial') {
        await loadArtifacts(j.id)
      }
    } finally {
      setStarting(false)
    }
  }

  function changeGene(raw: string) {
    const nextGene = raw.trim().toUpperCase()
    setGene(nextGene)
    if (!isAcceptedDemoGene(nextGene)) {
      setUseAccepted(false)
    }
    setJob(null)
    setArtifact(null)
    setArtifactResolved(false)
    setSessionMessage(null)
    const next = new URLSearchParams(params)
    if (nextGene) {
      next.set('gene', nextGene)
    } else {
      next.delete('gene')
    }
    next.delete('job')
    setParams(next, { replace: true })
  }

  function toggle(section: string) {
    setSelected((prev) =>
      prev.includes(section) ? prev.filter((s) => s !== section) : [...prev, section],
    )
  }

  return (
    <div className="fade-in mx-auto max-w-3xl space-y-8">
      <div>
        <h1 className="text-2xl font-medium tracking-tight text-text">
          Generate HD-Focused Gene Dossier
        </h1>
        <p className="mt-2 text-sm text-text-secondary">
          Generate a dossier from validated deterministic workflows, or open the accepted demo
          dossier instantly.
        </p>
      </div>

      <div className="surface-card space-y-5 p-5">
        <label className="block">
          <span className="mb-1.5 block text-xs text-text-muted">Gene</span>
          <input
            list="accepted-demo-genes"
            value={gene}
            onChange={(e) => changeGene(e.target.value)}
            placeholder="e.g. LRPAP1"
            spellCheck={false}
            autoCapitalize="characters"
            className="w-full rounded-xl border border-border bg-bg-secondary px-3 py-2.5 text-sm uppercase text-text outline-none"
          />
          <datalist id="accepted-demo-genes">
            {ACCEPTED_DEMO_GENES.map((symbol) => (
              <option key={symbol} value={symbol} />
            ))}
          </datalist>
        </label>

        <div>
          <p className="mb-2 text-xs text-text-muted">Sections</p>
          <div className="grid gap-2 sm:grid-cols-2">
            {SECTIONS.map((section) => (
              <label
                key={section}
                className={cn(
                  'flex cursor-pointer items-center gap-2.5 rounded-xl border px-3 py-2.5 text-sm transition',
                  selected.includes(section)
                    ? 'border-accent/35 bg-accent/5 text-text'
                    : 'border-border text-text-secondary',
                )}
              >
                <input
                  type="checkbox"
                  checked={selected.includes(section)}
                  onChange={() => toggle(section)}
                  className="accent-[var(--color-accent)]"
                />
                {section}
              </label>
            ))}
          </div>
        </div>

        <div className="space-y-1.5">
          <label
            className={cn(
              'flex items-center justify-between gap-4 rounded-xl border border-border bg-bg-secondary px-3 py-2.5 text-sm',
              acceptedAvailable ? 'text-text-secondary' : 'text-text-muted',
            )}
          >
            <span>Use accepted dossier (instant)</span>
            <input
              type="checkbox"
              checked={acceptedAvailable && useAccepted}
              disabled={!acceptedAvailable}
              onChange={(e) => setUseAccepted(e.target.checked)}
              className="accent-[var(--color-accent)]"
            />
          </label>
          {!acceptedAvailable && (
            <p className="px-1 text-xs text-text-muted">
              Accepted dossiers are available only for SREBF2 and CDH10. Other genes run a fresh
              job.
            </p>
          )}
        </div>

        <button
          type="button"
          onClick={() => void onGenerate()}
          disabled={starting || selected.length === 0 || !gene.trim()}
          className="rounded-xl bg-accent px-5 py-2.5 text-sm font-medium text-bg transition hover:brightness-110 disabled:opacity-40"
        >
          {starting ? 'Starting…' : 'Generate Dossier'}
        </button>
      </div>

      {sessionMessage && (
        <div className="surface-card flex flex-wrap items-center justify-between gap-3 p-4 text-sm text-text-secondary">
          <span>{sessionMessage}</span>
          <Link to="/reports" className="text-accent hover:underline">
            Open Reports
          </Link>
        </div>
      )}

      {job && <JobProgress job={job} />}

      {(job?.status === 'Completed' || job?.status === 'Partial') && (
        <div className="surface-card flex flex-wrap gap-2 p-5">
          {!artifactResolved ? (
            <span className="px-1 py-2.5 text-sm text-text-muted">
              Loading report artifact…
            </span>
          ) : artifact?.id ? (
            <Link
              to={`/reports/${artifact.id}`}
              className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text transition hover:border-accent/40"
            >
              View Interactive Report
            </Link>
          ) : (
            <span className="px-1 py-2.5 text-sm text-text-muted">
              Report artifact is not available.
            </span>
          )}
          {artifact?.pdfUrl && (
            <a
              href={artifact.pdfUrl}
              className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text-secondary"
            >
              Download PDF
            </a>
          )}
          {artifactResolved && artifact && !artifact.pdfUrl && (
            <span className="px-1 py-2.5 text-sm text-text-muted">
              {artifact.reportOrigin === 'generated' ? 'PDF not generated' : 'PDF unavailable'}
            </span>
          )}
          <button
            type="button"
            className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text-secondary"
          >
            Supplementary Files
          </button>
          <Link
            to={`/evidence?gene=${encodeURIComponent(job.geneSymbol)}`}
            className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text-secondary"
          >
            Inspect Provenance
          </Link>
        </div>
      )}
    </div>
  )
}
