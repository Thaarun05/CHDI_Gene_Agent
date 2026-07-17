# Gene Dossier Platform - Implementation Plan

Provenance-first platform that generates CHDI-style gene dossiers for Huntington's disease
research genes (SREBF2, HTT, MSH3, FAN1, PMS2, MLH1, and others).

---

## 1. Project Goal

Build a **provenance-first** system where every reported fact traces back to a `source_id`,
and every `source_id` traces back to a raw API response, artifact, screenshot, or manual note.

The LLM is **never** the source of truth. The source of truth is retrieved evidence and raw
API responses. The system must run retrieval, normalization, and reporting **even with no LLM
API key configured**.

Primary deliverable at this stage:

```
python scripts/run_srebf2_full_api_pass.py
```

...which creates a dossier run for SREBF2, attempts every configured source, stores raw
artifacts + evidence records in SQLite, generates a source coverage report, and produces a
partial CHDI-style markdown report with verified, source-cited claims.

---

## 2. Architecture Summary

```
Validated biomedical APIs
  -> raw source responses            (raw_store: files on disk + content hash)
  -> normalized evidence records     (normalize/*: source-level factual units)
  -> provenance database             (SQLite via SQLModel: runs, api_runs, artifacts, evidence)
  -> report tables / figures         (rendering: presentation layer)
  -> LLM-written report sections     (synthesis: optional, evidence-constrained)
  -> claim verification              (verification: rule-based first)
  -> final gene dossier              (data/outputs/{run_id}_report.md)
```

Layer separation is mandatory:

- **API clients** (`tools/`): build request, call with timeout, never raise, return `ToolResult`.
  They do NOT normalize.
- **Normalizers** (`normalize/`): turn raw responses into `EvidenceRecord`s. No network calls.
- **Raw store** (`raw_store.py`): persists raw artifacts with content hashing.
- **Database** (`db.py`): provenance store for runs, api_runs, artifacts, evidence.
- **Coverage** (`coverage.py`): never silently omit a source; report status for every source.
- **Verification** (`verification.py`): rule-based checks that every claim cites existing sources.
- **Synthesis** (`synthesis.py`): optional LLM section writing, constrained to evidence records.

### Reference assets (already available, external to this repo)

- `../CHDI_Data_APIs_Gene_Report_SREBF2.md`: validated endpoint formats, chained identifiers,
  gotchas, and the canonical CHDI section list. Use as the endpoint spec of record.
- `../Gene-Agent/`: prior simpler implementation for pattern reference only; superseded here.

---

## 3. Core Development Rule: One File At A Time

Do NOT generate the whole project at once. For each step:

1. State the single file being created/modified.
2. Explain why it is needed.
3. Create/modify only that one file.
4. Run the smallest relevant check or test.
5. Show the result.
6. State what the next file should be.
7. STOP and wait for explicit approval ("continue" / "go to the next file").

Only advance when the current file's check passes and approval is given.

---

## 4. Full Implementation Order

Progress legend: [ ] pending, [~] in progress, [x] done, [-] deferred.

### Phase 0 - Planning and setup
- [~] 1. `IMPLEMENTATION_PLAN.md` (this file)
- [ ] 2. `pyproject.toml` - deps + package config (target Python >=3.11)
- [ ] 3. `.env.example` - env var template (no real keys)
- [ ] 4. `README.md` - overview + quickstart
- [ ] 5. `src/gene_dossier/__init__.py`
- [ ] 6. `src/gene_dossier/config.py` - settings, paths, key loading

### Phase 1 - Core models and tests
- [ ] 7. `src/gene_dossier/models.py` - enums + Pydantic v2 models
- [ ] 8. `tests/test_models.py`

### Phase 2 - Source IDs
- [ ] 9. `src/gene_dossier/source_ids.py` - deterministic source ID generation
- [ ] 10. `tests/test_source_ids.py`

### Phase 3 - Raw artifact storage
- [ ] 11. `src/gene_dossier/raw_store.py` - save/load raw artifacts + content hash
- [ ] 12. `tests/test_raw_store.py`

### Phase 4 - Source registry
- [ ] 13. `src/gene_dossier/source_registry.py` - full source map (A/B/C)
- [ ] 14. `tests/test_source_registry.py`

### Phase 5 - Database
- [ ] 15. `src/gene_dossier/db.py` - SQLModel engine + tables + helpers
- [ ] 16. `tests/test_db.py` (basic smoke test)

### Phase 6 - Coverage reporting
- [ ] 17. `src/gene_dossier/coverage.py` - build + write coverage report (md + json)
- [ ] 18. `tests/test_coverage.py`

### Phase 7 - Priority A API clients (one file at a time)
- [ ] 19. `tools/ncbi_gene.py`
- [ ] 20. `tools/pubmed.py`
- [ ] 21. `tools/ensembl.py`
- [ ] 22. `tools/uniprot.py`
- [ ] 23. `tools/gtex.py`
- [ ] 24. `tools/string_db.py`
- [ ] 25. `tools/reactome.py`
- [ ] 26. `tools/clinvar.py`
- [ ] 27. `tools/opentargets.py`
- [ ] 28. `tools/mousemine.py`
- [ ] 29. `tools/ctd.py`
- [ ] 30. `tools/chembl.py`
- [ ] 31. `tools/pubchem.py`
- [ ] 32. `tools/nih_reporter.py`

### Phase 8 - Priority B API clients (one file at a time)
- [ ] 33. `tools/geo.py`
- [ ] 34. `tools/harmonizome.py`
- [ ] 35. `tools/biogrid.py`
- [ ] 36. `tools/wikipathways.py`
- [ ] 37. `tools/alphafold.py`
- [ ] 38. `tools/pdbe.py`
- [ ] 39. `tools/cdd.py`
- [ ] 40. `tools/ncbi_datasets.py`
- [ ] 41. `tools/ucsc.py`

### Phase 9 - Priority C scaffolds
- [ ] 42. `tools/allen_brain.py`
- [ ] 43. `tools/brainrnaseq.py`
- [ ] 44. `tools/patents.py`
- [ ] 45. `tools/antibodies.py`
- [ ] 46. `tools/omim.py`

### Phase 10 - Normalizers (one file at a time)
- [ ] 47. `normalize/gene_identity.py`
- [ ] 48. `normalize/literature.py`
- [ ] 49. `normalize/protein.py`
- [ ] 50. `normalize/expression.py`
- [ ] 51. `normalize/ppi.py`
- [ ] 52. `normalize/pathways.py`
- [ ] 53. `normalize/chemicals.py`
- [ ] 54. `normalize/variants.py`
- [ ] 55. `normalize/model_organisms.py`
- [ ] 56. `normalize/grants.py`
- (also: `normalize/perturbation.py`, `normalize/transcription_factors.py` per final structure)

### Phase 11 - Verification and report generation
- [ ] 57. `src/gene_dossier/verification.py`
- [ ] 58. `tests/test_verification.py`
- [ ] 59. `src/gene_dossier/synthesis.py`
- [ ] 60. `src/gene_dossier/rendering.py`

### Phase 12 - Workflow
- [ ] 61. `src/gene_dossier/workflow.py` - `run_gene_dossier_full_api_pass(...)`
- [ ] 62. `scripts/run_srebf2_full_api_pass.py`

### Phase 13 - Retrieval skeleton
- [ ] 63. `src/gene_dossier/retrieval.py` - keyword search + metadata filters

### Phase 14 - FastAPI
- [ ] 64. `src/gene_dossier/api/main.py`

### Phase 15 - Utility scripts
- [ ] 65. `scripts/run_source_smoke_tests.py`
- [ ] 66. `scripts/print_source_coverage_report.py`

---

## 5. Source Priority Levels

**Priority A** - full client + deep normalizer target:
NCBI Gene, PubMed, UniProt, Ensembl, GTEx, STRING, Reactome, ClinVar, Open Targets,
MouseMine/MGI, CTD, ChEMBL, PubChem, NIH RePORTER.

**Priority B** - client + raw artifact storage first, basic normalizer OK:
GEO, Harmonizome, BioGRID, WikiPathways, AlphaFold, PDBe, CDD, NCBI Datasets, UCSC.

**Priority C** - scaffold; mark manual / semi-structured / requires_key / deferred:
Allen Brain Atlas, BrainRNASeq, patents, antibodies, OMIM, DrugBank, NCATS, ERC grants.

### Key / access notes
- `NCBI_API_KEY` optional (higher rate limits) for NCBI Gene, PubMed, ClinVar.
- `BIOGRID_ACCESSKEY` required -> else `requires_key`.
- `OMIM_API_KEY` required -> else `requires_key`.
- `SERPAPI_API_KEY` for patents -> else `manual` / `requires_key`.
- Antibodies, Allen Brain, BrainRNASeq: treat as semi-structured / manual if no clean API.

### Validated SREBF2 anchors (from reference report)
Entrez `6721` (mouse `20788`), Ensembl `ENSG00000198911`, UniProt `Q12772`,
GTEx GENCODE `ENSG00000198911.11`, MGI `MGI:107585`. Prefer exact official-symbol match;
never blindly trust the first search hit.

---

## 6. Data Model Summary (built in `models.py`)

Enums: `EvidenceGrade` (A-F), `SourceType`, `AssertionType`, `SourceStatus`.

Models: `DossierRun`, `ApiRun`, `RawArtifact`, `EvidenceRecord`, `ReportSection`, `Claim`,
`VerificationResult`, `ToolResult`, `SourceCoverageResult` (fields per spec).

`EvidenceGrade`: A=direct human genetic/curated causal; B=human expression/eQTL/disease;
C=curated protein/pathway/PPI; D=mouse/cell model; E=predicted/computational/text-mining;
F=weak mention / needs review.

`SourceStatus`: success, failed, deferred, manual, requires_key, partial, skipped,
not_implemented.

---

## 7. Report Structure (CHDI-style)

Sections rendered when evidence exists (else explicitly marked not available / failed /
requires key / deferred / manual review):

1. General gene information
2. Gene aliases and identifiers
3. Conservation / orthologs
4. Known structure / domains
5. AlphaFold / PDBe / CDD
6. Homologues
7. Tissue and cell expression
8. GEO perturbations
9. Transcription factors
10. Protein-protein interactions
11. CTD perturbations
12. Chemical tools
13. eQTLs
14. ClinVar / OMIM / Open Targets / SNPs
15. Pathways
16. Knockouts / model phenotypes
17. Major labs / literature
18. Antibodies
19. Patents
20. NIH/ERC grants
21. Missing / deferred / manual sources
22. Verification warnings

---

## 8. Verification Rules (MVP, rule-based)

- Every claim must cite >= 1 `source_id`.
- Every cited `source_id` must exist in `evidence_records`.
- Flag causal language: causes, proves, drives, therapeutic target, clinically validated,
  pathogenic, disease-modifying.
- Causal language + evidence grade below A (or no explicit disease evidence)
  -> `warning` / `human_review`.
- Flag unsupported claims. Return `VerificationResult` objects.

---

## 9. Acceptance Criteria

1. Codebase builds one file at a time; each file checked before advancing.
2. `python scripts/run_srebf2_full_api_pass.py` runs end to end.
3. Creates a dossier run for SREBF2 and attempts all configured sources.
4. Saves raw artifacts for every successful source.
5. Logs failed / deferred / manual / requires_key sources (never silent).
6. Creates evidence records for sources with implemented normalizers; stores in SQLite.
7. Generates a source coverage report (`.md` + `.json`).
8. Generates a partial CHDI-style markdown report that separates completed / partial /
   missing / manual-review / failed / requires-key.
9. Every generated claim cites `source_id`s; verification runs and flags weak/unsupported.
10. Code stays modular: one client per source, separate normalizers, raw store, evidence DB,
    and report generation.

---

## 10. Deferred (TODOs - not implemented yet)

- HDinHD MCP integration (leave architecture notes only).
- Chroma vector search, hybrid RAG reranking, evidence-sufficiency checks.
- Streamlit UI, React UI.
- Figure artifacts, richer report tables.
- Gene comparison workflow.
- LLM-based (beyond rule-based) verification.
- Human review UI.
- Report export to DOCX / PDF.

---

## 11. Testing Approach

- Each file gets a minimal check (import, instantiation, or unit test).
- Normalizer tests use mocked sample responses (no live network).
- Client tests confirm import + graceful failure returns `ToolResult`.
- `pytest` for the `tests/` suite; external-API-free where possible.
