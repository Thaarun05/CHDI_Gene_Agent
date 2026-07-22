# Gene Dossier Platform

A **provenance-first** platform that generates CHDI-style gene dossiers for Huntington's
disease research genes (SREBF2, HTT, MSH3, FAN1, PMS2, MLH1, and others).

This is **not** a chatbot-first project. Facts come from validated biomedical APIs and their
raw responses. The LLM is used only for optional section synthesis; it is **never** the
source of truth. Chroma is an optional index only — structured truth stays in the
provenance DB and on-disk raw artifacts.

## Core principle: provenance first

- Every fact traces back to a `source_id`.
- Every `source_id` traces back to a raw API response, artifact, or manual note.
- No report claim exists without cited `source_id`s.
- Retrieval, normalization, and reporting run **even with no LLM API key**.

## Pipeline

```
Validated biomedical APIs
  -> raw source responses         (data/raw + content hash)
  -> normalized evidence records  (source-level factual units)
  -> provenance database          (SQLite or Postgres / Supabase)
  -> optional Chroma index        (display_text only; never source of truth)
  -> report tables / sections     (deterministic; optional LLM synthesis)
  -> claim verification           (rule-based)
  -> gene dossier                 (markdown + Rancho HTML/PDF + coverage)
```

## Tech stack

Python 3.11+, FastAPI, Pydantic v2, SQLModel (SQLite / Postgres), httpx/requests,
LangGraph (workflow), LangChain (optional LLM), Chroma (optional index), pytest.
Deferred: Streamlit/React UI, hybrid RAG reranking, HDinHD MCP.

## Repository layout

```
CHDI_Gene_Agent/
  IMPLEMENTATION_PLAN.md
  pyproject.toml
  .env.example
  data/{raw,outputs,indexes}/
  src/gene_dossier/
    config, models, db, raw_store, source_registry, coverage
    workflow, synthesis, verification, retrieval, rendering
    report_schema, rancho_report
    tools/       # one API client per source (ToolResult; never raises)
    normalize/   # raw responses -> EvidenceRecords (no network)
    api/         # FastAPI app
    assets/      # Rancho report branding
  scripts/
    run_srebf2_full_api_pass.py
    run_source_smoke_tests.py
    print_source_coverage_report.py
  tests/
```

## Quickstart

```bash
# 1. Create a virtual environment (Python 3.11+)
python3 -m venv .venv && source .venv/bin/activate

# 2. Install (editable) with dev extras
pip install -e ".[dev]"

# 3. Configure environment (all keys optional; missing keys degrade gracefully)
cp .env.example .env

# 4. SREBF2 full API pass (live network; soft-fails per source)
python scripts/run_srebf2_full_api_pass.py
# Useful flags: --no-rancho --no-pdf --no-db --sources GTEx,STRING --allow-llm

# 5. Source smoke tests / coverage printer
python scripts/run_source_smoke_tests.py
python scripts/print_source_coverage_report.py --gene SREBF2

# 6. API
uvicorn gene_dossier.api.main:app --reload

# 7. Tests (mocked / offline; no live APIs required)
pytest
```

Outputs land in `data/outputs/` (or `--output-dir`):

- `{dossier_run_id}_report.md` — debug CHDI-style markdown
- `{dossier_run_id}_rancho.html` / `.pdf` — polished Rancho visual dossier
- `{dossier_run_id}_source_coverage.md` / `.json` — per-source status

## HTTP API (MVP)

| Method | Path | Role |
|--------|------|------|
| GET | `/health` | Liveness; `database=sqlite\|postgres\|other` (no raw URL) |
| GET | `/version` | Package version |
| GET | `/sources` | Full source registry |
| GET | `/sources/summary` | Counts by priority / status |
| GET | `/sources/{name}` | One source |
| POST | `/dossier/runs` | Start LangGraph pass (`wait=true` sync; default background) |
| GET | `/dossier/runs/{id}` | Persisted run status |
| GET | `/dossier/runs/{id}/evidence` | Evidence from DB |
| GET | `/dossier/runs/{id}/coverage` | Coverage rows |
| POST | `/dossier/runs/{id}/search` | Keyword/metadata evidence search |

`wait=false` + `persist_db=false` returns **422** (background runs must be pollable).

## Source coverage

The platform never silently omits a source. Every configured source reports one of:
`success`, `failed`, `deferred`, `manual`, `requires_key`, `partial`, `skipped`,
`not_implemented` — with raw artifact path, evidence count, and any error.

Priority levels:

- **A** (full client + deep normalizer): NCBI Gene, PubMed, UniProt, Ensembl, GTEx, STRING,
  Reactome, ClinVar, Open Targets, MouseMine/MGI, CTD, ChEMBL, PubChem, NIH RePORTER.
- **B** (client + raw storage, basic normalizer): GEO, Harmonizome, BioGRID, WikiPathways,
  AlphaFold, PDBe, CDD, NCBI Datasets, UCSC.
- **C** (scaffold; manual/semi-structured/requires_key): Allen Brain, BrainRNASeq, patents,
  antibodies, OMIM, DrugBank, NCATS, ERC grants.

## Status

Phases 0–15 of `IMPLEMENTATION_PLAN.md` (numbered steps 1–66) are implemented.
Still deferred: HDinHD MCP, hybrid RAG, Streamlit/React UI, richer figure artifacts,
gene comparison, LLM-based verification beyond rules, human review UI.
