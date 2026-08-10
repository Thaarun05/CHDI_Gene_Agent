import { useState, type ReactNode } from 'react'
import { Menu } from 'lucide-react'
import { Sidebar } from '@/components/Sidebar'
import { EvidenceDrawer } from '@/components/EvidenceDrawer'
import { EvidenceDrawerProvider } from '@/context/EvidenceDrawerContext'
import { cn } from '@/lib/utils'

export function AppShell({
  children,
  className,
}: {
  children: ReactNode
  className?: string
}) {
  const [mobileOpen, setMobileOpen] = useState(false)
  const [collapsed, setCollapsed] = useState(false)

  return (
    <EvidenceDrawerProvider>
      <div className="flex min-h-screen bg-bg">
        {/* Desktop / tablet sidebar */}
        <div
          className={cn(
            'sticky top-0 hidden h-screen shrink-0 transition-[width] duration-200 md:block',
            collapsed ? 'w-0 overflow-hidden' : 'w-[250px]',
          )}
        >
          <Sidebar />
        </div>

        {/* Mobile drawer */}
        <div
          className={cn(
            'fixed inset-0 z-40 bg-black/50 transition-opacity md:hidden',
            mobileOpen ? 'opacity-100' : 'pointer-events-none opacity-0',
          )}
          onClick={() => setMobileOpen(false)}
        />
        <div
          className={cn(
            'fixed top-0 left-0 z-50 h-full transition-transform duration-200 md:hidden',
            mobileOpen ? 'translate-x-0' : '-translate-x-full',
          )}
        >
          <Sidebar onClose={() => setMobileOpen(false)} />
        </div>

        <div className="flex min-w-0 flex-1 flex-col">
          <header className="sticky top-0 z-30 flex items-center gap-3 border-b border-border bg-bg/90 px-4 py-3 backdrop-blur md:hidden">
            <button
              type="button"
              onClick={() => setMobileOpen(true)}
              className="rounded-lg p-2 text-text-secondary"
              aria-label="Open menu"
            >
              <Menu className="size-5" />
            </button>
            <span className="text-sm font-medium text-text">CHDI Gene Intelligence</span>
          </header>

          {/* Tablet collapse control */}
          <div className="hidden border-b border-border px-4 py-2 md:block lg:hidden">
            <button
              type="button"
              onClick={() => setCollapsed((c) => !c)}
              className="text-xs text-text-muted transition hover:text-text"
            >
              {collapsed ? 'Show sidebar' : 'Collapse sidebar'}
            </button>
          </div>

          <main className={cn('mx-auto w-full max-w-[1280px] flex-1 px-4 py-6 sm:px-6 lg:px-8', className)}>
            {children}
          </main>
        </div>
      </div>
      <EvidenceDrawer />
    </EvidenceDrawerProvider>
  )
}
