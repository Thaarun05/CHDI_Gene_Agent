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

export type ResearchMode = 'auto' | 'deep_research' | 'stored_only'

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
  publicEvidenceRef: string
  sourceName: string
}

export interface AskResponse {
  status: string
  question: string
  geneSymbol: string | null
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
    distinct_source_count?: number
    direct_count?: number
    supporting_count?: number
    contextual_count?: number
    excluded_count?: number
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
  answerSections?: AnswerSection[]
  comparisonDecision?: ComparisonDecision | null
  evidenceItems?: PublicEvidenceItem[]
  contextualEvidence?: PublicEvidenceItem[]
  activitySummary?: ActivitySummary
  costSummary?: CostSummary
  baseEvidenceRef?: string | null
  reusedToolRunRefs?: string[]
  createdToolRunRefs?: string[]
  toolRunRefs?: string[]
  dossierRunRefs?: string[]
  evidenceUniverseRefs?: Record<string, PublicEvidenceUniverse>
}

export interface AnswerSection {
  key: 'status' | 'direct_answer' | 'conditional_conclusion' | 'key_findings' | 'evidence_by_dimension'
  title: string
  paragraphs: string[]
}

export interface ComparisonDecision {
  outcome: 'not_rankable' | 'dimension_specific_difference' | 'conditional_preference' | 'supported_preference'
  summary: string
  preferred_gene?: string | null
  criterion?: string | null
  limitations: string[]
}

export interface PublicEvidenceItem {
  public_evidence_ref: string
  gene_symbol: string
  source_name: string
  public_identifier?: string | null
  title?: string | null
  source_url?: string | null
  evidence_need: string
  designation: 'direct' | 'supporting' | 'contextual' | 'excluded'
  display_text: string
  retrieved_at?: string | null
  backing_record_count: number
  exclusion_reason?: string | null
}

export interface ActivitySummary {
  requirements_planned: number
  persisted_retrieval_completed: boolean
  tools_executed: number
  tools_failed: number
  runs_reused: number
  tools_skipped: number
  accepted_evidence: number
  rejected_evidence: number
  skip_reasons: Record<string, number>
  rejection_reasons: Record<string, number>
}

export interface CostSummary {
  estimated_model_cost_usd?: number | null
  external_tool_cost_usd?: number | null
  actual_billed_cost_usd?: number | null
  cost_basis: string[]
  provider_reported_usage: Record<string, unknown>
}

export interface PublicEvidenceUniverse {
  geneSymbol: string
  baseEvidenceRef?: string | null
  explicitRunRefs: string[]
  reusedToolRunRefs: string[]
  createdToolRunRefs: string[]
  toolRunRefs: string[]
  dossierRunRefs: string[]
  evidenceUniverse: EvidenceUniverseName
}

export interface AgentComparisonCell {
  status: ComparisonStrength
  summary: string
  evidenceCount: number
  publicEvidenceRefs: string[]
  citationOrdinals: number[]
  distinctSourceCount: number
  directCount: number
  supportingCount: number
  excludedCount: number
  directionalityKnown: boolean
  hasConflict: boolean
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
