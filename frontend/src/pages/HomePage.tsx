import { useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { FileText, MessageSquareText, GitCompare, FlaskConical, ScanSearch, ListOrdered } from 'lucide-react'
import { TopBar } from '@/components/TopBar'
import { MolecularHero } from '@/components/MolecularHero'
import { SearchComposer } from '@/components/SearchComposer'
import { WorkflowCard } from '@/components/WorkflowCard'

export function HomePage() {
  const [query, setQuery] = useState('')
  const navigate = useNavigate()

  function handleSubmit() {
    const q = query.trim()
    if (!q) return
    const upper = q.toUpperCase()
    if (upper === 'SREBF2' || upper === 'CDH10') {
      navigate(`/genes/${upper}`)
      return
    }
    navigate(`/ask?q=${encodeURIComponent(q)}`)
  }

  return (
    <div className="fade-in flex min-h-[calc(100vh-3rem)] flex-col">
      <TopBar className="mb-8" />

      <div className="flex flex-1 flex-col items-center justify-center pb-10">
        <MolecularHero className="mb-8" />

        <p className="text-sm tracking-wide text-accent-secondary">CHDI Gene Intelligence</p>
        <h1 className="mt-3 max-w-2xl text-center text-3xl font-medium tracking-tight text-balance text-text sm:text-4xl">
          What would you like to investigate?
        </h1>
        <p className="mt-3 max-w-xl text-center text-sm leading-relaxed text-text-secondary sm:text-[15px]">
          Explore therapeutic targets using deterministic biological evidence, provenance-aware
          retrieval, and agentic scientific workflows.
        </p>

        <div className="mt-10 w-full max-w-2xl">
          <SearchComposer value={query} onChange={setQuery} onSubmit={handleSubmit} />
        </div>

        <div className="mt-10 grid w-full max-w-4xl gap-4 sm:grid-cols-3">
          <WorkflowCard
            title="Generate HD Gene Dossier"
            description="Build a provenance-backed Huntington's disease target dossier from validated biological sources."
            icon={FileText}
            to="/generate"
          />
          <WorkflowCard
            title="Ask Evidence Question"
            description="Ask a scientific question and retrieve only the evidence required to answer it."
            icon={MessageSquareText}
            to="/ask"
          />
          <WorkflowCard
            title="Compare Genes"
            description="Compare therapeutic targets across biological evidence dimensions."
            icon={GitCompare}
            to="/compare"
          />
        </div>

        <div className="mt-4 grid w-full max-w-4xl gap-3 sm:grid-cols-3">
          <WorkflowCard
            title="Therapeutic Target Assessment"
            description="Structured assessment of target suitability for HD programs."
            icon={FlaskConical}
            comingNext
          />
          <WorkflowCard
            title="Evidence Gap Analysis"
            description="Identify missing or weak evidence categories for a target."
            icon={ScanSearch}
            comingNext
          />
          <WorkflowCard
            title="Target Prioritization"
            description="Rank candidates using evidence coverage — not AI confidence scores."
            icon={ListOrdered}
            comingNext
          />
        </div>
      </div>
    </div>
  )
}
