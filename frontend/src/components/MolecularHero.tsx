import { cn } from '@/lib/utils'

/** CSS-only scientific network sphere for the home hero. */
export function MolecularHero({ className }: { className?: string }) {
  return (
    <div className={cn('relative mx-auto h-36 w-36', className)} aria-hidden>
      <div className="absolute inset-0 rounded-full bg-[radial-gradient(circle_at_30%_30%,rgba(140,203,94,0.22),transparent_55%),radial-gradient(circle_at_70%_65%,rgba(89,195,182,0.18),transparent_50%)] blur-sm" />
      <div className="absolute inset-3 rounded-full border border-white/10 bg-[radial-gradient(circle_at_40%_35%,rgba(255,255,255,0.08),transparent_60%)] shadow-[inset_0_0_40px_rgba(140,203,94,0.12)]" />
      <div className="absolute inset-0">
        {[
          [18, 42],
          [72, 28],
          [84, 68],
          [36, 78],
          [52, 48],
          [22, 64],
          [68, 52],
        ].map(([x, y], i) => (
          <span
            key={i}
            className="absolute size-1.5 rounded-full bg-accent/70"
            style={{ left: `${x}%`, top: `${y}%` }}
          />
        ))}
        <svg className="absolute inset-0 h-full w-full opacity-40" viewBox="0 0 100 100">
          <line x1="22" y1="46" x2="52" y2="50" stroke="#8CCB5E" strokeWidth="0.4" />
          <line x1="52" y1="50" x2="74" y2="32" stroke="#59C3B6" strokeWidth="0.4" />
          <line x1="52" y1="50" x2="84" y2="70" stroke="#8CCB5E" strokeWidth="0.35" />
          <line x1="52" y1="50" x2="40" y2="78" stroke="#59C3B6" strokeWidth="0.35" />
          <line x1="26" y1="66" x2="52" y2="50" stroke="#8CCB5E" strokeWidth="0.3" />
        </svg>
      </div>
      <div className="absolute inset-[-12%] rounded-full bg-accent/5 blur-2xl" />
    </div>
  )
}
