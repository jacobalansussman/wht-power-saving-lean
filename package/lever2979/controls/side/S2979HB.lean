import Work.XHostT.Inst.S2979Facts
import Work.XHostT.HF.GeomInst

/-!
# (key: heads-first) #2979's head bins: `headBins sched2979` has 16,086 elements

p-hfirst's `headBins` (A4 bins with remaining climb 48 in which no slot is hosted) on `sched2979`: with positive
residuals (`resK2979_pos`) p-hfirst's `headBins_eq` reduces it to the A4 bins with remainder 48, counted by
`sched2979_heads`.  This is the `hb.card` of the bridge's `hLlen`.
x-tables (wht45/lever2979, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.XHostT.Inst
open OAI.PowerSaving.XHostT OAI.PowerSaving.XHostT.K.DH2979

/-- **the head bins of #2979**: 16,086. -/
theorem headBins2979_card : (headBins sched2979).card = 16086 := by
  rw [headBins_eq sched2979 (resK ct resT 74 408) sched2979_fits resK2979_pos]
  exact sched2979_heads

end OAI.PowerSaving.XHostT.Inst

#print axioms OAI.PowerSaving.XHostT.Inst.headBins2979_card
