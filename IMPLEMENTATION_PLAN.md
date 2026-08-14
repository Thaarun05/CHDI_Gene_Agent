# Gene Dossier Platform — Current Status and Remaining Work

Provenance-first platform that generates CHDI-style (Rancho-layout) gene dossiers for
Huntington's disease research. Demo genes with accepted artifacts: **SREBF2** and **CDH10**.
The architecture can attempt other symbols through the backend; arbitrary-gene UI is incomplete.

This document is the **current engineering status**, not a from-scratch build checklist.
Original phased checkboxes from the first implementation are in the appendix (all core
phases are done). Do not treat unchecked historical items as unfinished work.

See [README.md](README.md) for how to run the stack.

---

## 1. Project Goal

Build a **provenance-first** system where every reported fact traces back to a `source_id`,
and every `source_id` traces back to a raw API response, artifact, screenshot, or manual note.

The LLM is **never** the source of truth. The source of truth is retrieved evidence and raw
API responses. Retrieval, normalization, and reporting must run **even with no LLM API key**.

### Current primary paths

| Path | Entry | What it produces |
|---|---|---|
| Generate UI / `POST /api/jobs` | `run_section_bundle` | Deterministic Rancho HTML/PDF for selected `1a`–`7a` |
| CLI | `scripts/run_section_bundle.py` | Same bundle; default keys **1a–4a** unless `--sections` is passed |
| Full source pass | `scripts/run_srebf2_full_api_pass.py` / `POST /dossier/runs` | LangGraph: `CLIENT_DISPATCH` sources, coverage, optional Rancho |
| Ask | `POST /api/ask` | `ScientificAgentService` over stored evidence (optional grounded LLM) |
| Compare | `POST /api/compare` | EvidenceRecord **coverage** matrix (not a scientific ranking) |

Default `DATABASE_URL` in `.env.example` is **SQLite**
(`sqlite:///data/gene_dossier.db`). Postgres (including Supabase) is supported via the
same setting; it is not the implied local default.

Polished interactive output is Rancho **HTML/PDF** (`rancho_report.py`), not
`data/outputs/{run_id}_report.md`. Markdown from `rendering.py` is a provenance/debug view.

---

## 2. Architecture Summary

```
Validated biomedical APIs / browser captures
  -> raw source responses            (raw_store: local files + content hash)
  -> normalized evidence records     (normalize/* and section-owned nodes)
  -> provenance database             (SQLModel → SQLite or Postgres via DATABASE_URL)
  -> optional Chroma index           (semantic retrieval; NOT source of truth)
  -> report_presentation             (deterministic Rancho blocks; no LLM in the bundle)
  -> rancho_report HTML/PDF          (final polished dossier)
```

Generate / section bundle (no LLM):

```
identity → section-owned acquisition → EvidenceRecords
  → report_presentation → rancho HTML/PDF
```

Full LangGraph pass (`workflow.py`):

```
create_dossier_run
  -> resolve_gene_identity
  -> call_source_clients
  -> save_raw_artifacts
  -> normalize_evidence
  -> (optional) index_evidence_in_chroma
  -> build_report_sections          # synthesis; force_deterministic=True by default
  -> verify_claims
  -> render_outputs
```

Ask (`src/gene_dossier/agent/`):

```
question → plan_scientific_question → evidence universe
  → retrieve (Chroma then keyword) → sufficiency
  → optional allowlisted acquisition (section_bundle or source_workflow)
  → grounded synthesis (LLM fills slots only; citations are server-rendered)
```

### Provenance truth hierarchy

1. **Raw artifacts** (files on disk; object storage not implemented) = original source material.
   Postgres/SQLite stores only **metadata** (path, hash, timestamps, source, urls) — never
   the large raw response body.
2. **SQLModel database** = structured source of truth for runs, evidence, claims, coverage,
   verification, and generated-report pointers. SQLite locally; Postgres when `DATABASE_URL`
   is set that way.
3. **Chroma** = semantic search index over evidence records only. Not the source of truth.

### Database: dual-backend via `DATABASE_URL`

| Mode | Example URL | When to use |
|------|-------------|-------------|
| Local SQLite | `sqlite:///data/gene_dossier.db` | Default local / offline |
| Postgres (e.g. Supabase) | `postgresql+psycopg://...@host:5432/postgres` | Shared/dev, multi-machine |
| In-memory SQLite | `sqlite://` or `sqlite:///:memory:` | Unit tests |

Driver: **`psycopg[binary]`**. `db.py` uses SQLite pragmas only for SQLite URLs.
Never hardcode credentials; use `.env` only.

Tables: `dossier_runs`, `api_runs`, `raw_artifacts` (pointers only), `evidence_records`,
`report_sections`, `claims`, `verification_results`, `source_coverage_results`,
`generated_reports`. Schema is created with `SQLModel.metadata.create_all` — **no Alembic
migrations**.

### Layer separation (mandatory)

- **API clients** (`tools/`): timeout, never raise, return `ToolResult`. They do NOT normalize.
- **Normalizers** (`normalize/`): `ToolResult` → `EvidenceRecord`. No network calls.
- **Section nodes** (`section_*.py`): polished 1c–7a acquisition/presentation; many own HTTP
  instead of the generic workflow client for that source.
- **Raw store** (`raw_store.py`): artifact **bytes** on disk.
- **Database** (`db.py`): SQLModel provenance store.
- **Coverage** (`coverage.py`): never silently omit a registered source.
- **Retrieval** (`retrieval.py`): keyword/metadata first; Chroma second.
- **Verification** (`verification.py`): rule-based citation and causal-language checks.
- **Bundle presentation** (`report_presentation.py` + `rancho_report.py`): deterministic Rancho UI.
- **Workflow synthesis** (`synthesis.py`): optional LLM section markdown for the LangGraph path;
  default `force_deterministic=True`.
- **Ask synthesis** (`agent/synthesis.py`): grounded prose slots; not used to build section HTML.

### Runtime processes (verified)

- **Backend framework:** FastAPI (`src/gene_dossier/api/main.py`).
- **App server:** Uvicorn — `uvicorn gene_dossier.api.main:app --host 0.0.0.0 --port 8001`.
- **Frontend dev:** Vite on port 5173.
- **Frontend demo host:** Vercel static SPA (`https://chdi-gene-agent.vercel.app`). No `vercel.json` in repo.
- **Demo API exposure:** local Uvicorn + Cloudflare Quick Tunnel. **Not production.**
- **Jobs:** in-process `threading.Thread` + `_JOB_STORE` (lost on restart). No Celery/Redis.

---

## 3. What Is Implemented

### Polished Rancho sections (bundle-supported)

`SUPPORTED_SECTION_BUNDLE_KEYS`: `1a` `1b` `1c` `1d` `1e` `2a` `2b` `2c` `3a` `4a` `5a` `5b` `6a` `7a`.

| Default | Keys |
|---|---|
| CLI `scripts/run_section_bundle.py` (omit `--sections`) | **1a–4a** (`DEFAULT_SECTION_BUNDLE_KEYS`) |
| API Generate / `_HD_DOSSIER_DEFAULT_SECTIONS` | **1a–7a** |

1a and 1b have no `section_1a.py` / `section_1b.py`; they are assembled from gene identity
and UCSC normalization plus `report_presentation.py`.

### Frontend

React 19 + TypeScript + Vite 8 + Tailwind 4. Routes: `/`, `/generate`, `/ask`, `/compare`,
`/genes/:symbol`, `/reports`, `/reports/:id`, `/evidence`, `/history`.

### Tests

Offline/mocked `pytest` under `tests/` (section modules, clients, workflow, API frontend
routes, scientific agent). No frontend unit/e2e suite and no GitHub Actions CI.

### Source clients wired in `CLIENT_DISPATCH`

NCBI Gene, Ensembl, UniProt, PubMed, GTEx, STRING, Reactome, ClinVar, Open Targets,
MouseMine, CTD, ChEMBL, PubChem, NIH RePORTER, GEO, Harmonizome, BioGRID, WikiPathways,
AlphaFold, PDBe, CDD, NCBI Datasets, UCSC, Allen Brain Atlas, BrainRNASeq, Patents,
Antibodies, OMIM.

Section-owned modules also call OrthoDB, PubTator3, HBT, DropViz/GEO Profiles, NCATS, etc.
**DrugBank** is a stub; **ERC Grants** and **HDinHD** have no client.

`source_registry.py` still leaves `client_implemented` / `normalizer_implemented` at the
default **False**. Treat `CLIENT_DISPATCH` and section nodes as the real inventory.

---

## 4. Rancho TOC vs code

Canonical layout: `report_schema.py` `REPORT_SECTIONS` (15 majors + lettered subsections,
copied from `SREBF2_report.pdf` wording). Older `CHDI_REPORT_SECTIONS` in `synthesis.py`
is a parallel 22-name list used by LangGraph markdown synthesis — do not confuse the two.

| Rancho slot | Status |
|---|---|
| 1a–1e, 2a–2c, 3a, 4a, 5a, 5b, 6a, 7a | Polished generators implemented |
| 7b Tractability | Schema only |
| 8a GTEx eQTL | Schema; GTEx client exists; no section module |
| 9a–c ClinVar / OMIM / Open Targets | Clients exist; no polished section |
| 9d SNPs3D | Schema only; no client |
| 10a Pathways | Reactome/WikiPathways clients; Ask can fetch Reactome; no Rancho section |
| 11a–b Knockouts | MouseMine client; no polished section |
| 12 Labs | Schema only |
| 13 Antibodies / 14 Patents | Scaffold clients; SREBF2-oriented query maps |
| 15a NIH RePORTER | Client + grants normalizer; no polished section |
| 15b ERC | Registry deferred; no client |
| HDinHD | Deferred |

---

## 5. Source access notes

- `NCBI_API_KEY` optional (rate limits) for NCBI Gene, PubMed, ClinVar, and related.
- `BIOGRID_ACCESSKEY` required or BioGRID is `requires_key`.
- `OMIM_API_KEY` required or OMIM is `requires_key`.
- `SERPAPI_API_KEY` for patents (and antibodies search); else `requires_key` / manual.
- `UCSC_BROWSER_API_KEY` optional for hgRenderTracks.
- LLM keys optional (OpenAI, Anthropic, NVIDIA NIM, Google Gemini). Without them, Ask and
  reports still run deterministically.
- Playwright Chromium required for several figure captures (CDD, AlphaFold viewer, BioGRID,
  GEO Profiles, Allen, DropViz).
- PyMuPDF (`fitz`) required for PDF/PNG; **not declared in `pyproject.toml`**. HTML still writes.

Validated SREBF2 anchors (from the reference report): Entrez `6721`, Ensembl `ENSG00000198911`,
UniProt `Q12772`, GTEx GENCODE `ENSG00000198911.11`. Prefer exact official-symbol match.

Licensing and terms of use for GTEx, BioGRID, OMIM, DrugBank, SerpAPI, and browser captures
**must be confirmed with CHDI**. Not specified in this repository.

---

## 6. Data model (`models.py`)

Enums: `EvidenceGrade` (A–F), `SourceType`, `AssertionType`, `SourceStatus`.

Models: `DossierRun`, `ApiRun`, `RawArtifact`, `EvidenceRecord`, `ReportSection`, `Claim`,
`VerificationResult`, `ToolResult`, `SourceCoverageResult`.

Provenance chain: `DossierRun → ApiRun → RawArtifact → EvidenceRecord → report / retrieval / answer`.

---

## 7. Verification rules (rule-based; no LLM)

Implemented in `verification.py`:

- Every claim must cite >= 1 `source_id`.
- Every cited `source_id` must exist in evidence.
- Flag causal language (causes, proves, drives, therapeutic target, pathogenic, …).
- Causal language + weak grade / no disease evidence → `warning` / `human_review`.

This does **not** certify biological correctness. CHDI scientific review remains mandatory.

---

## 8. Remaining work

### P0 — before real CHDI users

- Authentication / authorization (none today; CORS is `allow_origins=["*"]`).
- Replace Cloudflare-tunnel-to-laptop as the API.
- CI: pytest, frontend build, secret scan.
- Confirm secrets were never committed; rotate if needed.
- Production `VITE_USE_MOCKS=false` and a stable `VITE_API_BASE`.
- Declare PyMuPDF in packaging **or** document HTML-only as supported.
- Legal/ToS review of scraped and browser-captured sources.

### P1 — production launch / full Rancho

- Durable jobs (replace `_JOB_STORE`); object storage for artifacts; Postgres + migrations.
- Align ports: README/Uvicorn **8001** vs Vite proxy and `frontend/.env.example` **8000**.
- Fix stale `client_implemented` flags; Pydantic body for `POST /jobs` (currently a raw `dict`).
- Polished sections **7b** and **8–15** following the 4a/5a pattern, with SME review.
- Observability, backups, rate limits, request timeouts, worker isolation for Playwright.
- Frontend tests and accessibility.

### Still deferred (and still true)

- HDinHD MCP integration.
- Object store (S3 / Supabase Storage) for raw artifact bytes.
- LLM-based verification (beyond rules).
- DOCX export.
- Human review workflow productization.

### No longer deferred (do not plan these as greenfield)

- React UI (exists).
- PDF export (optional runtime via PyMuPDF).
- Gene comparison (coverage matrix, not ranking).
- Chroma Ask path and grounded LLM slots.
- Figure-rich Rancho sections 1a–7a.

### CHDI decisions required (not engineering alone)

Hosting, SSO, data classification, API licenses, LLM vendors, required report scope
(HD 1a–7a vs full 15-section Rancho TOC), who signs scientific acceptance.

---

## 9. Current acceptance criteria (what “done” means now)

1. `pip install -e ".[dev]"` and `uvicorn gene_dossier.api.main:app --port 8001` serve `/health`.
2. Accepted SREBF2/CDH10 reports open without live APIs (`use_existing_accepted`).
3. Fresh `POST /api/jobs` runs `run_section_bundle`, persists evidence, writes HTML
   (and PDF if PyMuPDF is present).
4. Sources soft-fail with coverage status; they are not silently omitted.
5. Ask returns cited evidence or abstains; it does not invent `source_id`s.
6. Compare reports EvidenceRecord counts, not druggability.
7. Without LLM keys, acquisition + deterministic reports + deterministic Ask still run.
8. `python -m pytest` passes offline (mocked HTTP).

Historical criterion “codebase builds one file at a time” is **retired**.

---

## 10. Testing approach

- `pytest` with `pythonpath = ["src"]`; fixtures and `respx` mocks; no live network unless a
  test says otherwise.
- DB tests use in-memory SQLite.
- Frontend: `npm run build` and `npm run lint` only (no component/e2e tests).
- No CI in `.github/`.

---

## 11. Dependency and database setup

```bash
pip install -e ".[dev]"
python -m playwright install chromium   # figure captures
pip install pymupdf                     # PDF until added to pyproject.toml
```

```bash
# Local / offline (default in .env.example)
DATABASE_URL=sqlite:///data/gene_dossier.db

# Optional Postgres
# DATABASE_URL=postgresql+psycopg://USER:PASSWORD@HOST:5432/postgres
```

LLM keys remain optional. Credentials remain in `.env` only.

---

## Appendix A — Original build process (historical)

The first implementation used a one-file-at-a-time approval loop. That process is **complete**
for the core platform. Do not wait on “continue / go to the next file” to treat the repo as
unbuilt.

## Appendix B — Original phase checklist (all core items done)

Legend: [x] done in tree, [~] partial, [-] deferred.

- [x] Package, config, `.env.example` (SQLite + Postgres examples), README
- [x] `models.py`, `source_ids.py`, `raw_store.py`, `source_registry.py` (+ tests)
- [x] `db.py` (SQLite + Postgres, coverage + generated_reports) and `tests/test_db.py`
- [x] `coverage.py`, `verification.py`, `synthesis.py`, `rendering.py`, `rancho_report.py`
- [x] Priority A/B clients in `CLIENT_DISPATCH`; Priority C scaffolds (Allen, BrainRNASeq,
      patents, antibodies, OMIM); DrugBank/ERC/HDinHD still deferred or stubbed
- [x] Normalizers under `normalize/`
- [x] `workflow.py`, `scripts/run_srebf2_full_api_pass.py`, `retrieval.py`
- [x] FastAPI `api/main.py`, React frontend, section bundle 1a–7a, Ask agent, Compare
- [x] `scripts/run_source_smoke_tests.py`, `scripts/print_source_coverage_report.py`
- [~] Registry `client_implemented` flags (always False; stale)
- [-] Streamlit extra in `pyproject.toml` (unused; React is the UI)
- [-] Alembic migrations, CI, Docker, auth, object storage
