import Work.XHostT.Inst.S5Facts
import Work.XHostT.HF.GeomInst

/-!
# (key: heads-first) #3945's head bins: `headBins sched3945` has 57,312 elements

p-hfirst's `headBins` (A4 bins with remaining climb 48 in which no slot is hosted) on `sched3945`: with positive
residuals (`resK3945_pos`) p-hfirst's `headBins_eq` reduces it to the A4 bins with remainder 48, counted by
`sched3945_heads` (x-tables' S2979HB.lean with the instance changed).
i-final (wht45/lever, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.XHostT.Inst
open OAI.PowerSaving.XHostT OAI.PowerSaving.XHostT.K.DS5

/-- **the head bins of #3945**: 57,312. -/
theorem headBins3945_card : (headBins sched3945).card = 57312 := by
  rw [headBins_eq sched3945 (resK ct resT 360 408) sched3945_fits resK3945_pos]
  exact sched3945_heads

end OAI.PowerSaving.XHostT.Inst

#print axioms OAI.PowerSaving.XHostT.Inst.headBins3945_card
