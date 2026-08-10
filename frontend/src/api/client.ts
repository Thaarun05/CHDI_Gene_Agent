/**
 * API service layer.
 *
 * Pages import these functions only — never call fetch() directly.
 * Flip USE_MOCKS to false when FastAPI endpoints are ready; keep signatures stable.
 */

import type {
  AskResponse,
  ComparisonResponse,
  EvidenceCoverageRow,
  EvidenceRecord,
  Gene,
  HistoryItem,
  RecentWorkItem,
  ReportArtifact,
  WorkflowJob,
} from '@/api/types'
import {
  askResponseSrebf2,
  compareResponse,
  completedJob,
  coverageByGene,
  createMockJob,
  evidenceRecords,
  genes,
  history,
  recentWork,
  reports,
} from '@/mocks/data'

const USE_MOCKS = true
const API_BASE = import.meta.env.VITE_API_BASE ?? '/api'

async function http<T>(path: string, init?: RequestInit): Promise<T> {
  const res = await fetch(`${API_BASE}${path}`, {
    headers: { 'Content-Type': 'application/json', ...(init?.headers ?? {}) },
    ...init,
  })
  if (!res.ok) {
    throw new Error(`API ${res.status}: ${path}`)
  }
  return res.json() as Promise<T>
}

function delay(ms = 350) {
  return new Promise((r) => setTimeout(r, ms))
}

const jobStore = new Map<string, WorkflowJob>()

export async function getGene(symbol: string): Promise<Gene> {
  if (USE_MOCKS) {
    await delay()
    const gene = genes[symbol.toUpperCase()]
    if (!gene) throw new Error(`Gene not found: ${symbol}`)
    return gene
  }
  return http<Gene>(`/genes/${encodeURIComponent(symbol)}`)
}

export async function getEvidenceCoverage(symbol: string): Promise<EvidenceCoverageRow[]> {
  if (USE_MOCKS) {
    await delay(200)
    return coverageByGene[symbol.toUpperCase()] ?? []
  }
  return http(`/genes/${encodeURIComponent(symbol)}/coverage`)
}

export async function listGeneEvidence(symbol: string): Promise<EvidenceRecord[]> {
  if (USE_MOCKS) {
    await delay()
    return evidenceRecords.filter(
      (e) => e.geneSymbol.toUpperCase() === symbol.toUpperCase(),
    )
  }
  return http(`/genes/${encodeURIComponent(symbol)}/evidence`)
}

export async function getEvidenceRecord(id: string): Promise<EvidenceRecord> {
  if (USE_MOCKS) {
    await delay(150)
    const rec = evidenceRecords.find((e) => e.id === id)
    if (!rec) throw new Error(`Evidence not found: ${id}`)
    return rec
  }
  return http(`/evidence/${encodeURIComponent(id)}`)
}

export async function listAllEvidence(filters?: {
  gene?: string
  source?: string
  type?: string
  section?: string
}): Promise<EvidenceRecord[]> {
  if (USE_MOCKS) {
    await delay()
    return evidenceRecords.filter((e) => {
      if (filters?.gene && e.geneSymbol !== filters.gene) return false
      if (filters?.source && e.sourceName !== filters.source) return false
      if (filters?.type && e.evidenceType !== filters.type) return false
      if (filters?.section && !e.section.toLowerCase().includes(filters.section.toLowerCase()))
        return false
      return true
    })
  }
  const qs = new URLSearchParams()
  if (filters?.gene) qs.set('gene', filters.gene)
  if (filters?.source) qs.set('source', filters.source)
  if (filters?.type) qs.set('type', filters.type)
  if (filters?.section) qs.set('section', filters.section)
  const q = qs.toString()
  return http(`/evidence${q ? `?${q}` : ''}`)
}

export async function startDossierJob(geneSymbol: string): Promise<WorkflowJob> {
  if (USE_MOCKS) {
    await delay(400)
    const job = createMockJob(geneSymbol)
    jobStore.set(job.id, job)
    // Simulate progress in background for demo UX.
    void simulateJob(job.id)
    return job
  }
  return http('/jobs', {
    method: 'POST',
    body: JSON.stringify({ gene_symbol: geneSymbol, job_type: 'hd_dossier' }),
  })
}

async function simulateJob(jobId: string) {
  const order = [2, 3, 4, 5] // stage indices to advance
  for (const idx of order) {
    await delay(900)
    const job = jobStore.get(jobId)
    if (!job) return
    job.stages = job.stages.map((s, i) => {
      if (i < idx) return { ...s, status: 'Complete' }
      if (i === idx) return { ...s, status: 'Running' }
      if (i === idx + 1) return { ...s, status: 'Queued' }
      return s
    })
    job.status = 'Running'
    jobStore.set(jobId, { ...job })
  }
  await delay(800)
  const done = jobStore.get(jobId)
  if (!done) return
  done.stages = done.stages.map((s) => ({ ...s, status: 'Complete' }))
  done.status = 'Completed'
  done.completedAt = new Date().toISOString()
  done.artifactIds = [
    reports.find((r) => r.geneSymbol === done.geneSymbol)?.id ?? 'rep-srebf2',
  ]
  jobStore.set(jobId, { ...done })
}

export async function getJob(jobId: string): Promise<WorkflowJob> {
  if (USE_MOCKS) {
    await delay(120)
    const job = jobStore.get(jobId) ?? completedJob
    return { ...job, stages: job.stages.map((s) => ({ ...s })) }
  }
  return http(`/jobs/${encodeURIComponent(jobId)}`)
}

export async function getJobArtifacts(jobId: string): Promise<ReportArtifact | null> {
  if (USE_MOCKS) {
    await delay(150)
    const job = jobStore.get(jobId) ?? completedJob
    const id = job.artifactIds?.[0]
    return reports.find((r) => r.id === id) ?? reports[0] ?? null
  }
  return http(`/jobs/${encodeURIComponent(jobId)}/artifacts`)
}

export async function listReports(): Promise<ReportArtifact[]> {
  if (USE_MOCKS) {
    await delay()
    return reports
  }
  return http('/reports')
}

export async function getReport(id: string): Promise<ReportArtifact> {
  if (USE_MOCKS) {
    await delay()
    const report = reports.find((r) => r.id === id)
    if (!report) throw new Error(`Report not found: ${id}`)
    return report
  }
  return http(`/reports/${encodeURIComponent(id)}`)
}

export async function askEvidenceQuestion(
  question: string,
  geneSymbol = 'SREBF2',
): Promise<AskResponse> {
  if (USE_MOCKS) {
    await delay(700)
    return {
      ...askResponseSrebf2,
      question,
      geneSymbol: geneSymbol.toUpperCase(),
    }
  }
  return http('/ask', {
    method: 'POST',
    body: JSON.stringify({ question, gene_symbol: geneSymbol }),
  })
}

export async function compareGenes(geneSymbols: string[]): Promise<ComparisonResponse> {
  if (USE_MOCKS) {
    await delay(500)
    return {
      ...compareResponse,
      genes: geneSymbols.length ? geneSymbols : compareResponse.genes,
    }
  }
  return http('/compare', {
    method: 'POST',
    body: JSON.stringify({ genes: geneSymbols }),
  })
}

export async function listHistory(): Promise<HistoryItem[]> {
  if (USE_MOCKS) {
    await delay()
    return history
  }
  return http('/history')
}

export async function listRecentWork(): Promise<RecentWorkItem[]> {
  if (USE_MOCKS) {
    await delay(100)
    return recentWork
  }
  return http('/recent')
}

export { USE_MOCKS }
