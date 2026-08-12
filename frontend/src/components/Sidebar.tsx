import { useEffect, useState } from 'react'
import { NavLink, Link } from 'react-router-dom'
import {
  Home,
  Dna,
  FileText,
  GitCompare,
  Database,
  History,
  Search,
  X,
  type LucideIcon,
} from 'lucide-react'
import { listRecentWork } from '@/api/client'
import type { RecentWorkItem } from '@/api/types'
import { cn } from '@/lib/utils'

const NAV: Array<{ to: string; label: string; icon: LucideIcon }> = [
  { to: '/', label: 'Home', icon: Home },
  { to: '/genes/SREBF2', label: 'Genes', icon: Dna },
  { to: '/reports', label: 'Reports', icon: FileText },
  { to: '/compare', label: 'Compare', icon: GitCompare },
  { to: '/evidence', label: 'Evidence', icon: Database },
  { to: '/history', label: 'History', icon: History },
]

export function Sidebar({
  onClose,
  className,
}: {
  onClose?: () => void
  className?: string
}) {
  const [recent, setRecent] = useState<RecentWorkItem[]>([])
  const [query, setQuery] = useState('')

  useEffect(() => {
    const refresh = () => void listRecentWork().then(setRecent)
    refresh()
    window.addEventListener('generated-report-updated', refresh)
    return () => window.removeEventListener('generated-report-updated', refresh)
  }, [])

  return (
    <aside
      className={cn(
        'flex h-full w-[250px] flex-col border-r border-border bg-sidebar',
        className,
      )}
    >
      <div className="flex items-start justify-between gap-2 px-4 pt-5 pb-4">
        <Link to="/" className="block" onClick={onClose}>
          <p className="text-[13px] leading-snug font-medium tracking-tight text-text">
            CHDI Gene Intelligence
          </p>
          <p className="mt-1 text-[11px] text-text-muted">Target intelligence</p>
        </Link>
        {onClose && (
          <button
            type="button"
            className="rounded-lg p-1.5 text-text-muted lg:hidden"
            onClick={onClose}
            aria-label="Close sidebar"
          >
            <X className="size-4" />
          </button>
        )}
      </div>

      <div className="px-3 pb-4">
        <div className="flex items-center gap-2 rounded-xl border border-border bg-bg-secondary px-3 py-2">
          <Search className="size-3.5 shrink-0 text-text-muted" />
          <input
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Search genes, reports, or evidence"
            className="w-full bg-transparent text-xs text-text outline-none placeholder:text-text-muted"
          />
        </div>
      </div>

      <nav className="flex-1 space-y-0.5 overflow-y-auto px-2">
        {NAV.map(({ to, label, icon: Icon }) => (
          <NavLink
            key={to}
            to={to}
            end={to === '/'}
            onClick={onClose}
            className={({ isActive }) =>
              cn(
                'flex items-center gap-2.5 rounded-xl px-3 py-2 text-sm transition',
                isActive
                  ? 'bg-white/[0.06] text-text'
                  : 'text-text-secondary hover:bg-white/[0.03] hover:text-text',
              )
            }
          >
            <Icon className="size-4 opacity-80" strokeWidth={1.75} />
            {label}
          </NavLink>
        ))}

        <div className="mt-6 px-3">
          <p className="mb-2 text-[10px] font-medium tracking-wider text-text-muted uppercase">
            Recent Work
          </p>
          <ul className="space-y-1">
            {recent
              .filter((r) => !query || r.label.toLowerCase().includes(query.toLowerCase()))
              .map((item) => (
                <li key={item.id}>
                  <Link
                    to={item.href}
                    onClick={onClose}
                    className="block truncate rounded-lg px-2 py-1.5 text-xs text-text-secondary transition hover:bg-white/[0.03] hover:text-text"
                  >
                    {item.label}
                  </Link>
                </li>
              ))}
          </ul>
        </div>
      </nav>

      <div className="border-t border-border px-4 py-4">
        <p className="text-xs text-text-muted">Research Workspace</p>
      </div>
    </aside>
  )
}
