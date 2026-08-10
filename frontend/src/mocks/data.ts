/**
 * MOCK / DEMO DATA ONLY — not live backend results.
 * Replace via api/client.ts when FastAPI endpoints are ready.
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

export const MOCK_NOTE = 'demo/mock data — not live FastAPI results'

export const genes: Record<string, Gene> = {
  SREBF2: {
    symbol: 'SREBF2',
    name: 'Sterol regulatory element-binding transcription factor 2',
    organism: 'Human',
    entrezGeneId: '6721',
    uniprotAccession: 'Q12772',
    summary:
      'Transcription factor that regulates cholesterol biosynthesis genes. Frequently assessed as an HD-relevant metabolic target in dossier workflows.',
  },
  CDH10: {
    symbol: 'CDH10',
    name: 'Cadherin-10',
    organism: 'Human',
    entrezGeneId: '1008',
    uniprotAccession: 'Q9Y6N8',
    summary:
      'Type II classical cadherin involved in cell–cell adhesion. Used as a contrast gene with limited chemical-tool tractability signals.',
  },
}

export const coverageByGene: Record<string, EvidenceCoverageRow[]> = {
  SREBF2: [
    { category: 'Gene Identity', status: 'Available', detail: 'NCBI + UniProt resolved' },
    { category: 'Structure', status: 'Available' },
    { category: 'Expression', status: 'Available', detail: 'GTEx / HBT figures present' },
    { category: 'GEO Perturbations', status: 'Available' },
    { category: 'Transcription Factors', status: 'Available' },
    { category: 'Protein Interactions', status: 'Available', detail: 'STRING + BioGRID' },
    { category: 'Chemical Perturbations', status: 'Not available', detail: 'No CTD perturbation records' },
    { category: 'Chemical Tools', status: 'Available', detail: 'ChEMBL workbook + literature tools' },
  ],
  CDH10: [
    { category: 'Gene Identity', status: 'Available' },
    { category: 'Structure', status: 'Available' },
    { category: 'Expression', status: 'Available' },
    { category: 'GEO Perturbations', status: 'Limited' },
    { category: 'Transcription Factors', status: 'Available' },
    { category: 'Protein Interactions', status: 'Available' },
    { category: 'Chemical Perturbations', status: 'Not available', detail: 'No CTD perturbation records' },
    {
      category: 'Chemical Tools',
      status: 'Limited',
      detail: 'ChEMBL no authoritative target; literature UNC0642 indirect',
    },
  ],
}

export const evidenceRecords: EvidenceRecord[] = [
  {
    id: 'EV-101',
    geneSymbol: 'SREBF2',
    sourceName: 'ChEMBL',
    evidenceType: 'chemical_tool',
    factType: 'section_7a_chembl_workbook',
    evidenceClass: 'target_linked_activity',
    section: 'Chemical tools',
    sourceIdentifier: 'CHEMBL1795166',
    retrievedAt: '2026-08-10T14:20:00Z',
    displayText:
      'Exact-target ChEMBL activities for SREBF2 (CHEMBL1795166) with assay relationship metadata in supplementary workbook.',
    status: 'Available',
    apiRunId: 'API-701',
    rawArtifactId: 'RAW-701',
    sourceUrl: 'https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL1795166/',
  },
  {
    id: 'EV-102',
    geneSymbol: 'SREBF2',
    sourceName: 'PubMed',
    evidenceType: 'chemical_tool',
    factType: 'section_7a_pubmed_tool',
    evidenceClass: 'literature_negative_effect',
    section: 'Chemical tools',
    sourceIdentifier: 'PMID 27756839',
    retrievedAt: '2026-08-10T14:22:00Z',
    displayText:
      'Betulin, a small-molecule inhibitor, abrogates SREBP-2 activation in the cited study.',
    status: 'Available',
    apiRunId: 'API-702',
    rawArtifactId: 'RAW-702',
    sourceUrl: 'https://pubmed.ncbi.nlm.nih.gov/27756839/',
  },
  {
    id: 'EV-103',
    geneSymbol: 'SREBF2',
    sourceName: 'PubChem',
    evidenceType: 'chemical_tool',
    factType: 'section_7a_pubchem_assay',
    section: 'Chemical tools',
    sourceIdentifier: 'AID focused set',
    retrievedAt: '2026-08-10T14:23:00Z',
    displayText: 'Focused PubChem bioassays linked to SREBF2 gene target (polished representatives).',
    status: 'Available',
    apiRunId: 'API-703',
    rawArtifactId: 'RAW-703',
  },
  {
    id: 'EV-201',
    geneSymbol: 'CDH10',
    sourceName: 'ChEMBL',
    evidenceType: 'chemical_tool',
    factType: 'section_7a_summary',
    evidenceClass: 'no_authoritative_target',
    section: 'Chemical tools',
    retrievedAt: '2026-08-10T14:25:00Z',
    displayText: 'No authoritative ChEMBL single-protein target resolved for CDH10 UniProt accession.',
    status: 'No Results',
    apiRunId: 'API-711',
    rawArtifactId: 'RAW-711',
  },
  {
    id: 'EV-202',
    geneSymbol: 'CDH10',
    sourceName: 'PubMed',
    evidenceType: 'chemical_tool',
    factType: 'section_7a_pubmed_tool',
    evidenceClass: 'indirect_pathway_effect',
    section: 'Chemical tools',
    sourceIdentifier: 'PMID 32292512',
    retrievedAt: '2026-08-10T14:26:00Z',
    displayText:
      'UNC0642 inhibits G9a/EHMT2 with downstream effects on CDH10 expression (indirect pathway effect).',
    status: 'Limited',
    apiRunId: 'API-712',
    rawArtifactId: 'RAW-712',
    sourceUrl: 'https://pubmed.ncbi.nlm.nih.gov/32292512/',
  },
  {
    id: 'EV-105',
    geneSymbol: 'SREBF2',
    sourceName: 'STRING',
    evidenceType: 'ppi',
    factType: 'section_5a_direct_partner',
    section: 'Protein-protein interaction (PPI) partners',
    retrievedAt: '2026-08-10T14:15:00Z',
    displayText: 'STRING direct functional partners for SREBF2 (demo coverage row).',
    status: 'Available',
    apiRunId: 'API-501',
    rawArtifactId: 'RAW-501',
  },
  {
    id: 'EV-106',
    geneSymbol: 'SREBF2',
    sourceName: 'GTEx',
    evidenceType: 'expression',
    factType: 'gtex_tissue_expression_summary',
    section: 'Tissue and cell expression',
    retrievedAt: '2026-08-10T14:10:00Z',
    displayText: 'GTEx tissue expression summary for SREBF2 across adult tissues.',
    status: 'Available',
    apiRunId: 'API-201',
    rawArtifactId: 'RAW-201',
  },
]

export const recentWork: RecentWorkItem[] = [
  { id: 'rw1', label: 'SREBF2 Gene Dossier', href: '/reports/rep-srebf2' },
  { id: 'rw2', label: 'CDH10 Chemical Tools', href: '/genes/CDH10' },
  { id: 'rw3', label: 'SREBF2 vs CDH10', href: '/compare' },
  { id: 'rw4', label: 'SREBF2 Pharmacology Question', href: '/ask' },
]

export const reports: ReportArtifact[] = [
  {
    id: 'rep-srebf2',
    geneSymbol: 'SREBF2',
    title: 'HD-Focused Gene Dossier',
    status: 'Completed',
    createdAt: '2026-08-10T14:36:00Z',
    sections: ['1a', '1b', '1c', '1d', '1e', '2a', '2b', '2c', '3a', '4a', '5a', '5b', '6a', '7a'],
    htmlUrl:
      'http://127.0.0.1:8901/data/outputs/section_validation/SREBF2_full_1a7a/407e1a4293c6424e8b6b830a1f0a7c60/section_1.html',
    pdfUrl:
      '/Users/thaarun/Desktop/Gene_Dossier/CHDI_Gene_Agent/data/outputs/section_validation/full_1a7a_delivery/SREBF2_sections_1a-7a.pdf',
  },
  {
    id: 'rep-cdh10',
    geneSymbol: 'CDH10',
    title: 'HD-Focused Gene Dossier',
    status: 'Completed',
    createdAt: '2026-08-10T14:33:00Z',
    sections: ['1a', '1b', '1c', '1d', '1e', '2a', '2b', '2c', '3a', '4a', '5a', '5b', '6a', '7a'],
    htmlUrl:
      'http://127.0.0.1:8901/data/outputs/section_validation/CDH10_full_1a7a/d94f392f4a3941d5a59f697f58d18234/section_1.html',
    pdfUrl:
      '/Users/thaarun/Desktop/Gene_Dossier/CHDI_Gene_Agent/data/outputs/section_validation/full_1a7a_delivery/CDH10_sections_1a-7a.pdf',
  },
]

export const history: HistoryItem[] = [
  {
    id: 'h1',
    geneLabel: 'SREBF2',
    workflow: 'HD Dossier',
    status: 'Completed',
    createdAt: '2026-08-10T14:36:00Z',
  },
  {
    id: 'h2',
    geneLabel: 'SREBF2',
    workflow: 'Evidence Question',
    status: 'Completed',
    createdAt: '2026-08-10T15:02:00Z',
  },
  {
    id: 'h3',
    geneLabel: 'SREBF2 vs CDH10',
    workflow: 'Gene Comparison',
    status: 'Completed',
    createdAt: '2026-08-10T15:10:00Z',
  },
  {
    id: 'h4',
    geneLabel: 'CDH10',
    workflow: 'HD Dossier',
    status: 'Completed',
    createdAt: '2026-08-10T14:33:00Z',
  },
]

const jobStagesTemplate = [
  { id: 's1', label: 'Resolving gene identity', status: 'Complete' as const },
  { id: 's2', label: 'Retrieving expression', status: 'Complete' as const },
  { id: 's3', label: 'Retrieving perturbations', status: 'Running' as const },
  { id: 's4', label: 'Retrieving protein interactions', status: 'Queued' as const },
  { id: 's5', label: 'Retrieving chemical evidence', status: 'Queued' as const },
  { id: 's6', label: 'Rendering report', status: 'Waiting' as const },
]

export function createMockJob(geneSymbol: string): WorkflowJob {
  return {
    id: `job-${geneSymbol.toLowerCase()}-${Date.now()}`,
    geneSymbol,
    jobType: 'hd_dossier',
    status: 'Running',
    stages: jobStagesTemplate.map((s) => ({ ...s })),
    createdAt: new Date().toISOString(),
  }
}

export const completedJob: WorkflowJob = {
  id: 'job-srebf2-demo',
  geneSymbol: 'SREBF2',
  jobType: 'hd_dossier',
  status: 'Completed',
  stages: jobStagesTemplate.map((s) => ({ ...s, status: 'Complete' })),
  createdAt: '2026-08-10T14:14:00Z',
  completedAt: '2026-08-10T14:36:00Z',
  artifactIds: ['rep-srebf2'],
  dossierRunId: 'cb9030ab81dc42db80b81dd15d48e653',
}

export const askResponseSrebf2: AskResponse = {
  status: 'answered',
  question: 'What evidence suggests SREBF2 can be pharmacologically manipulated?',
  geneSymbol: 'SREBF2',
  retrievalMethod: 'semantic',
  generationMethod: 'deterministic',
  embeddingBackend: 'local_minilm',
  baseEvidenceRunId: '407e1a4293c6424e8b6b830a1f0a7c60',
  toolRunIds: [],
  dossierRunIds: ['407e1a4293c6424e8b6b830a1f0a7c60'],
  evidenceUniverse: 'accepted_demo',
  summary:
    'Stored dossier evidence indicates SREBF2 has small-molecule and chemical-tool signals across ChEMBL, PubMed tool literature, and PubChem focused assays. These records support pharmacological manipulation as a research hypothesis, with source-specific limitations.',
  evidenceBlocks: [
    {
      sourceGroup: 'ChEMBL',
      items: [
        {
          text: 'Exact-target activity inventory exists for ChEMBL target CHEMBL1795166 with assay relationship metadata.',
          citationIds: ['c1'],
        },
      ],
    },
    {
      sourceGroup: 'PubMed',
      items: [
        {
          text: 'Literature tools such as betulin are reported to interfere with SREBP-2 activation.',
          citationIds: ['c2'],
        },
      ],
    },
    {
      sourceGroup: 'PubChem',
      items: [
        {
          text: 'Focused gene-linked bioassays are present in the polished PubChem table.',
          citationIds: ['c3'],
        },
      ],
    },
  ],
  limitations: [
    'DrugBank API access was unavailable in the stored run.',
    'Literature tool eligibility depends on local evidence spans; endogenous lipids are not treated as tools by default.',
    'Open Targets small-molecule tractability (section 7b) is not implemented.',
  ],
  citations: [
    { id: 'c1', label: 'ChEMBL', evidenceRecordId: 'EV-101', sourceName: 'ChEMBL' },
    {
      id: 'c2',
      label: 'PubMed PMID 27756839',
      evidenceRecordId: 'EV-102',
      sourceName: 'PubMed',
    },
    { id: 'c3', label: 'PubChem', evidenceRecordId: 'EV-103', sourceName: 'PubChem' },
  ],
  evidenceUsedCount: 3,
  sourcesCount: 3,
  sourcesUsed: ['ChEMBL', 'PubChem', 'PubMed'],
  toolsInvokedCount: 0,
  toolActivity: [],
  agentActivity: [
    'Resolved SREBF2',
    'Semantic retrieval attempted first',
    'Checking chemical-tool evidence',
    'Validating sources',
    'Building grounded answer',
  ],
}

export const compareResponse: ComparisonResponse = {
  genes: ['SREBF2', 'CDH10'],
  dimensions: [
    'Gene Identity',
    'Expression',
    'GEO Perturbations',
    'Protein Interactions',
    'Chemical Perturbations',
    'Chemical Tools',
  ],
  matrix: [
    {
      dimension: 'Gene Identity',
      cells: {
        SREBF2: {
          status: 'Available',
          summary: 'Entrez 6721 · UniProt Q12772',
          evidenceCount: 4,
          evidenceRecordIds: [],
        },
        CDH10: {
          status: 'Available',
          summary: 'Entrez 1008 · UniProt Q9Y6N8',
          evidenceCount: 4,
          evidenceRecordIds: [],
        },
      },
    },
    {
      dimension: 'Expression',
      cells: {
        SREBF2: {
          status: 'Available',
          summary: 'GTEx/HBT tissue profiles present',
          evidenceCount: 8,
          evidenceRecordIds: ['EV-106'],
        },
        CDH10: {
          status: 'Available',
          summary: 'Tissue expression dossier section complete',
          evidenceCount: 6,
          evidenceRecordIds: [],
        },
      },
    },
    {
      dimension: 'GEO Perturbations',
      cells: {
        SREBF2: {
          status: 'Available',
          summary: 'GEO profile charts captured',
          evidenceCount: 10,
          evidenceRecordIds: [],
        },
        CDH10: {
          status: 'Limited',
          summary: 'Fewer high-confidence GEO charts',
          evidenceCount: 3,
          evidenceRecordIds: [],
        },
      },
    },
    {
      dimension: 'Protein Interactions',
      cells: {
        SREBF2: {
          status: 'Available',
          summary: 'STRING + BioGRID partners',
          evidenceCount: 24,
          evidenceRecordIds: ['EV-105'],
        },
        CDH10: {
          status: 'Available',
          summary: 'STRING network available',
          evidenceCount: 18,
          evidenceRecordIds: [],
        },
      },
    },
    {
      dimension: 'Chemical Perturbations',
      cells: {
        SREBF2: {
          status: 'Not available',
          summary: 'No CTD perturbation records',
          evidenceCount: 0,
          evidenceRecordIds: [],
        },
        CDH10: {
          status: 'Not available',
          summary: 'No CTD perturbation records',
          evidenceCount: 0,
          evidenceRecordIds: [],
        },
      },
    },
    {
      dimension: 'Chemical Tools',
      cells: {
        SREBF2: {
          status: 'Available',
          summary: 'ChEMBL workbook + literature tools',
          evidenceCount: 20,
          evidenceRecordIds: ['EV-101', 'EV-102', 'EV-103'],
        },
        CDH10: {
          status: 'Limited',
          summary: 'No ChEMBL target; UNC0642 indirect only',
          evidenceCount: 2,
          evidenceRecordIds: ['EV-201', 'EV-202'],
        },
      },
    },
  ],
  narrative:
    'SREBF2 shows denser chemical-tool evidence in the accepted baseline, while CDH10 remains stronger on identity/expression/PPI but limited for direct chemical tools. Neither accepted baseline contains persisted CTD chemical-perturbation evidence. This comparison is evidence-matrix based, not a scored ranking.',
  evidenceUniverses: {
    SREBF2: {
      baseEvidenceRunId: '407e1a4293c6424e8b6b830a1f0a7c60',
      toolRunIds: [],
      dossierRunIds: ['407e1a4293c6424e8b6b830a1f0a7c60'],
      evidenceUniverse: 'accepted_demo',
    },
    CDH10: {
      baseEvidenceRunId: 'd94f392f4a3941d5a59f697f58d18234',
      toolRunIds: [],
      dossierRunIds: ['d94f392f4a3941d5a59f697f58d18234'],
      evidenceUniverse: 'accepted_demo',
    },
  },
}
