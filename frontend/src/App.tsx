import { BrowserRouter, Navigate, Route, Routes } from 'react-router-dom'
import { AppShell } from '@/components/AppShell'
import { HomePage } from '@/pages/HomePage'
import { GeneWorkspacePage } from '@/pages/GeneWorkspacePage'
import { GenerateDossierPage } from '@/pages/GenerateDossierPage'
import { ReportViewerPage } from '@/pages/ReportViewerPage'
import { AskPage } from '@/pages/AskPage'
import { ComparePage } from '@/pages/ComparePage'
import { EvidencePage } from '@/pages/EvidencePage'
import { ReportsPage } from '@/pages/ReportsPage'
import { HistoryPage } from '@/pages/HistoryPage'

export default function App() {
  return (
    <BrowserRouter>
      <AppShell>
        <Routes>
          <Route path="/" element={<HomePage />} />
          <Route path="/genes/:symbol" element={<GeneWorkspacePage />} />
          <Route path="/generate" element={<GenerateDossierPage />} />
          <Route path="/reports" element={<ReportsPage />} />
          <Route path="/reports/:id" element={<ReportViewerPage />} />
          <Route path="/ask" element={<AskPage />} />
          <Route path="/compare" element={<ComparePage />} />
          <Route path="/evidence" element={<EvidencePage />} />
          <Route path="/history" element={<HistoryPage />} />
          <Route path="*" element={<Navigate to="/" replace />} />
        </Routes>
      </AppShell>
    </BrowserRouter>
  )
}
