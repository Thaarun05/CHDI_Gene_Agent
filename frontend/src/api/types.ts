/** Shared TypeScript contracts aligned with planned FastAPI routes. */

export type EvidenceStatus =
  | 'Available'
  | 'Limited'
  | 'Missing'
  | 'No Results'
  | 'Unavailable'
  | 'Source Error'
  | 'Stale'
  | 'Running'

export type JobStatus = 'Queued' | 'Running' | 'Completed' | 'Partial' | 'Failed'

export type JobStageStatus = 'Complete' | 'Running' | 'Queued' | 'Waiting' | 'Failed'

export interface Gene {
  symbol: string
  name: string
  organism: string
  entrezGeneId: string
  uniprotAccession: string
  summary: string
}

export interface EvidenceCoverageRow {
  category: string
  status: EvidenceStatus
  detail?: string
}

export interface EvidenceRecord {
  id: string
  geneSymbol: string
  sourceName: string
  evidenceType: string
  factType: string
  evidenceClass?: string
  section: string
  subsection?: string
  sourceIdentifier?: string
  retrievedAt: string
  displayText: string
  status?: EvidenceStatus
  apiRunId?: string
  rawArtifactId?: string
  sourceUrl?: string
}

export interface ApiRun {
  id: string
  sourceName: string
  endpointName: string
  requestUrl: string
  success: boolean
  retrievedAt: string
}

export interface RawArtifact {
  id: string
  filePath: string
  contentHash: string
  artifactType: string
  apiRunId: string
}

export interface ReportArtifact {
  id: string
  geneSymbol: string
  title: string
  status: JobStatus
  createdAt: string
  sections: string[]
  htmlUrl?: string
  pdfUrl?: string
}

export interface WorkflowJobStage {
  id: string
  label: string
  status: JobStageStatus
}

export interface WorkflowJob {
  id: string
  geneSymbol: string
  jobType: string
  status: JobStatus
  stages: WorkflowJobStage[]
  createdAt: string
  completedAt?: string
  artifactIds?: string[]
}

export interface Citation {
  id: string
  label: string
  evidenceRecordId: string
  sourceName: string
}

export interface AskResponse {
  question: string
  geneSymbol: string
  summary: string
  evidenceBlocks: Array<{
    sourceGroup: string
    items: Array<{ text: string; citationIds: string[] }>
  }>
  limitations: string[]
  citations: Citation[]
  evidenceUsedCount: number
  sourcesCount: number
  toolsInvokedCount: number
  agentActivity: string[]
}

export interface ComparisonCell {
  status: EvidenceStatus
  summary: string
  evidenceCount: number
  evidenceRecordIds: string[]
}

export interface ComparisonResponse {
  genes: string[]
  dimensions: string[]
  matrix: Array<{
    dimension: string
    cells: Record<string, ComparisonCell>
  }>
  narrative: string
}

export interface HistoryItem {
  id: string
  geneLabel: string
  workflow: string
  status: JobStatus
  createdAt: string
}

export interface RecentWorkItem {
  id: string
  label: string
  href: string
}
