# Gene Dossier Platform

A **provenance-first** platform that generates CHDI-style gene dossiers for Huntington's
disease research genes (SREBF2, HTT, MSH3, FAN1, PMS2, MLH1, and others).

This is **not** a chatbot-first project. Facts come from validated biomedical APIs and their
raw responses. The LLM is used only for section synthesis, Q&A, and claim verification, and
is **never** treated as the source of truth.

## Core principle: provenance first

- Every fact traces back to a `source_id`.
- Every `source_id` traces back to a raw API response, artifact, or manual note.
- No report claim exists without cited `source_id`s.
- The system runs retrieval, normalization, and reporting **even with no LLM API key**.



## Pipeline

```
Validated biomedical APIs
  -> raw source responses         (data/raw + content hash)
  -> normalized evidence records  (source-level factual units)
  -> provenance database          (SQLite)
  -> report tables / figures      (presentation layer)
  -> LLM-written report sections  (optional, evidence-constrained)
  -> claim verification           (rule-based first)
  -> final gene dossier           (data/outputs/{run_id}_report.md)
```



## Tech stack

Python 3.11+, FastAPI, Pydantic v2, SQLModel + SQLite, httpx/requests, pytest.
Optional (deferred): LangGraph/LangChain (LLM), Chroma (RAG), Streamlit (UI).

## Repository layout

```
gene_dossier/
  IMPLEMENTATION_PLAN.md      # build order + acceptance criteria (read this first)
  pyproject.toml
  .env.example
  data/{raw,outputs,indexes}/
  src/gene_dossier/           # config, models, db, raw_store, source registry, coverage,
                              # workflow, synthesis, verification, retrieval, rendering
    tools/                    # one API client per source (returns ToolResult, never raises)
    normalize/                # raw responses -> EvidenceRecords (no network calls)
    api/                      # FastAPI app
  scripts/                    # CLI entry points
  tests/                      # pytest suite (mocked responses, no live APIs)
```



## Quickstart

```bash
# 1. Create a virtual environment (Python 3.11+)
python3 -m venv .venv && source .venv/bin/activate

# 2. Install (editable) with dev extras
pip install -e ".[dev]"

# 3. Configure environment (all keys optional; missing keys degrade gracefully)
cp .env.example .env    # then fill in any keys you have

# 4. Run the SREBF2 full API pass (available once the workflow is built)
python scripts/run_srebf2_full_api_pass.py

# 5. Run tests
pytest
```

Outputs land in `data/outputs/`:

- `{dossier_run_id}_report.md` - partial CHDI-style dossier
- `{dossier_run_id}_source_coverage.md` / `.json` - per-source status



## Source coverage

The platform never silently omits a source. Every configured source reports one of:
`success`, `failed`, `deferred`, `manual`, `requires_key`, `partial`, `skipped`,
`not_implemented` - with raw artifact path, evidence count, and any error.

Priority levels:

- **A** (full client + deep normalizer): NCBI Gene, PubMed, UniProt, Ensembl, GTEx, STRING,
Reactome, ClinVar, Open Targets, MouseMine/MGI, CTD, ChEMBL, PubChem, NIH RePORTER.
- **B** (client + raw storage, basic normalizer): GEO, Harmonizome, BioGRID, WikiPathways,
AlphaFold, PDBe, CDD, NCBI Datasets, UCSC.
- **C** (scaffold; manual/semi-structured/requires_key): Allen Brain, BrainRNASeq, patents,
antibodies, OMIM, DrugBank, NCATS, ERC grants.



## Development rule

Built one file at a time. See `IMPLEMENTATION_PLAN.md` for the full ordered checklist,
data models, verification rules, acceptance criteria, and deferred work (HDinHD MCP, Chroma
vector search, hybrid RAG, Streamlit/React UI, DOCX/PDF export).



