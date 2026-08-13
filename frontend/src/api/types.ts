/** Shared TypeScript contracts aligned with planned FastAPI routes. */

export type EvidenceStatus =
  | 'Available'
  | 'Limited'
  | 'Missing'
  | 'Not available'
  | 'No Results'
  | 'Unavailable'
  | 'Source Error'
  | 'Stale'
  | 'Running'

export type ComparisonStrength =
  | 'Strong'
  | 'Moderate'
  | 'Limited'
  | 'Weak'
  | 'Missing'

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

export type EvidenceUniverseName =
  | 'accepted_demo'
  | 'explicit_run'
  | 'accepted_demo_with_tool_overlay'
  | 'explicit_run_with_tool_overlay'
  | 'latest_generated'
  | 'latest_generated_with_tool_overlay'
  | 'no_base_evidence'
  | 'tool_overlay_only'
  | 'multi_gene'

export interface EvidenceUniverseMeta {
  baseEvidenceRunId?: string | null
  reusedToolRunIds?: string[]
  createdToolRunIds?: string[]
  toolRunIds: string[]
  dossierRunIds: string[]
  evidenceUniverse: EvidenceUniverseName
}

export interface EvidenceCoverageResponse extends EvidenceUniverseMeta {
  geneSymbol: string
  rows: EvidenceCoverageRow[]
}

export interface EvidenceRecord {
  id: string
  geneSymbol: string
  sourceName: string
  evidenceType: string
  factType: string
  evidenceClass?: string
  evidenceGrade?: string
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

export interface EvidenceListResponse extends EvidenceUniverseMeta {
  geneSymbol: string
  records: EvidenceRecord[]
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
  reportOrigin: 'accepted' | 'generated'
  dossierRunId: string
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
  dossierRunId?: string
  sectionKeys?: string[]
  errors?: string[]
}

export interface Citation {
  id: string
  label: string
  evidenceRecordId: string
  sourceName: string
}

export interface AskResponse {
  status: string
  question: string
  geneSymbol: string
  contextGene?: string | null
  summary: string
  retrievalMethod: string
  generationMethod: string
  embeddingBackend: string
  baseEvidenceRunId?: string | null
  reusedToolRunIds?: string[]
  createdToolRunIds?: string[]
  toolRunIds: string[]
  dossierRunIds: string[]
  evidenceUniverse: EvidenceUniverseName
  evidenceBlocks: Array<{
    sourceGroup: string
    items: Array<{ text: string; citationIds: string[] }>
  }>
  limitations: string[]
  citations: Citation[]
  evidenceUsedCount: number
  sourcesCount: number
  sourcesUsed: string[]
  toolsInvokedCount: number
  toolActivity: Array<Record<string, unknown>>
  agentActivity: string[]
  plannerMethod?: string | null
  intent?: string | null
  resolvedEntities?: {
    genes?: string[]
    diseases?: string[]
    biological_processes?: string[]
    pathways?: string[]
    chemicals?: string[]
  }
  evidenceRequirements?: Array<{
    id: string
    label: string
    genes: string[]
    evidence_need: string
    required: boolean
    minimum_support: number
  }>
  requirementAssessments?: Array<{
    requirement_id: string
    gene_symbol: string
    evidence_need: string
    required: boolean
    minimum_support: number
    status: string
    qualifying_count: number
    detail: string
  }>
  evidenceUniverses?: Record<string, EvidenceUniverseMeta>
  evidenceGaps?: string[]
  comparisonDimensions?: string[]
  comparisonMatrix?: Array<{
    dimension: string
    cells: Record<string, AgentComparisonCell>
  }>
  evidenceCategories?: Array<Record<string, unknown>>
  structuredGaps?: Array<Record<string, unknown>>
  recommendations?: Array<Record<string, unknown>>
  citationRegistry?: Array<Record<string, unknown>>
  sourceAttempts?: Array<Record<string, unknown>>
  retrievalTimestamps?: string[]
  failures?: Array<Record<string, unknown>>
  metadata?: Record<string, unknown>
}

export interface AgentComparisonCell {
  status: ComparisonStrength
  summary: string
  evidenceCount: number
  evidenceRecordIds: string[]
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
  evidenceUniverses: Record<string, EvidenceUniverseMeta>
}

export interface Artifact {
  id: string
  path: string
  artifactType: string
  exists: boolean
}

export interface JobArtifactsResponse {
  jobId: string
  dossierRunId?: string | null
  report: ReportArtifact | null
  supplementaryArtifacts: Artifact[]
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
