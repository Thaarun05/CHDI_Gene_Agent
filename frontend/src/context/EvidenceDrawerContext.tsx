import {
  createContext,
  useCallback,
  useContext,
  useMemo,
  useState,
  type ReactNode,
} from 'react'
import type { EvidenceRecord } from '@/api/types'
import { getEvidenceRecord } from '@/api/client'

interface EvidenceDrawerContextValue {
  open: boolean
  record: EvidenceRecord | null
  loading: boolean
  openEvidence: (idOrRecord: string | EvidenceRecord) => Promise<void>
  closeEvidence: () => void
}

const EvidenceDrawerContext = createContext<EvidenceDrawerContextValue | null>(null)

export function EvidenceDrawerProvider({ children }: { children: ReactNode }) {
  const [open, setOpen] = useState(false)
  const [record, setRecord] = useState<EvidenceRecord | null>(null)
  const [loading, setLoading] = useState(false)

  const openEvidence = useCallback(async (idOrRecord: string | EvidenceRecord) => {
    setOpen(true)
    if (typeof idOrRecord !== 'string') {
      setRecord(idOrRecord)
      setLoading(false)
      return
    }
    setLoading(true)
    try {
      const rec = await getEvidenceRecord(idOrRecord)
      setRecord(rec)
    } catch {
      setRecord(null)
    } finally {
      setLoading(false)
    }
  }, [])

  const closeEvidence = useCallback(() => {
    setOpen(false)
  }, [])

  const value = useMemo(
    () => ({ open, record, loading, openEvidence, closeEvidence }),
    [open, record, loading, openEvidence, closeEvidence],
  )

  return (
    <EvidenceDrawerContext.Provider value={value}>{children}</EvidenceDrawerContext.Provider>
  )
}

export function useEvidenceDrawer() {
  const ctx = useContext(EvidenceDrawerContext)
  if (!ctx) throw new Error('useEvidenceDrawer must be used within EvidenceDrawerProvider')
  return ctx
}
