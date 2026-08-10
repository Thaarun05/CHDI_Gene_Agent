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

export function GenerateDossierPage() {
  const [params] = useSearchParams()
  const [gene, setGene] = useState(params.get('gene')?.toUpperCase() || 'SREBF2')
  const [selected, setSelected] = useState<string[]>([...SECTIONS])
  const [useAccepted, setUseAccepted] = useState(true)
  const [job, setJob] = useState<WorkflowJob | null>(null)
  const [artifact, setArtifact] = useState<ReportArtifact | null>(null)
  const [starting, setStarting] = useState(false)

  useEffect(() => {
    if (!job || job.status === 'Completed' || job.status === 'Failed') return
    const t = setInterval(() => {
      void getJob(job.id).then((j) => {
        setJob(j)
        if (j.status === 'Completed') {
          void getJobArtifacts(j.id).then((artifacts) => setArtifact(artifacts.report))
        }
      })
    }, 700)
    return () => clearInterval(t)
  }, [job])

  async function onGenerate() {
    setStarting(true)
    setArtifact(null)
    try {
      const j = await startDossierJob(gene, {
        sectionKeys: sectionKeysFor(selected),
        useExistingAccepted: useAccepted,
      })
      setJob(j)
      if (j.status === 'Completed') {
        const artifacts = await getJobArtifacts(j.id)
        setArtifact(artifacts.report)
      }
    } finally {
      setStarting(false)
    }
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
          <select
            value={gene}
            onChange={(e) => setGene(e.target.value)}
            className="w-full rounded-xl border border-border bg-bg-secondary px-3 py-2.5 text-sm text-text outline-none"
          >
            <option value="SREBF2">SREBF2</option>
            <option value="CDH10">CDH10</option>
          </select>
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

        <label className="flex items-center justify-between gap-4 rounded-xl border border-border bg-bg-secondary px-3 py-2.5 text-sm text-text-secondary">
          <span>Use accepted dossier (instant)</span>
          <input
            type="checkbox"
            checked={useAccepted}
            onChange={(e) => setUseAccepted(e.target.checked)}
            className="accent-[var(--color-accent)]"
          />
        </label>

        <button
          type="button"
          onClick={() => void onGenerate()}
          disabled={starting || selected.length === 0}
          className="rounded-xl bg-accent px-5 py-2.5 text-sm font-medium text-bg transition hover:brightness-110 disabled:opacity-40"
        >
          {starting ? 'Starting…' : 'Generate Dossier'}
        </button>
      </div>

      {job && <JobProgress job={job} />}

      {job?.status === 'Completed' && (
        <div className="surface-card flex flex-wrap gap-2 p-5">
          <Link
            to={`/reports/${artifact?.id ?? 'rep-srebf2'}`}
            className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text transition hover:border-accent/40"
          >
            View Interactive Report
          </Link>
          <a
            href={artifact?.pdfUrl || '#'}
            className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text-secondary"
          >
            Download PDF
          </a>
          <button
            type="button"
            className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text-secondary"
          >
            Supplementary Files
          </button>
          <Link
            to="/evidence"
            className="rounded-xl border border-border bg-bg-secondary px-4 py-2.5 text-sm text-text-secondary"
          >
            Inspect Provenance
          </Link>
        </div>
      )}
    </div>
  )
}
