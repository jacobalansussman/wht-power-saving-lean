import Work.XHostT.Inst.S5HFacts
import Work.XHostT.HF.GeomInst

/-!
# (key: heads-first) design B's head bins: `headBins sched3945H` has 80,640 elements

p-hfirst's `headBins` (A4 bins with remaining climb 48 in which no slot is hosted) on `sched3945H`: with positive
residuals (`resK3945H_pos`) p-hfirst's `headBins_eq` reduces it to the A4 bins with remainder 48, counted by
`sched3945H_heads` (i-final's S5HB.lean with the instance and the count changed).
bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.XHostT.Inst
open OAI.PowerSaving.XHostT OAI.PowerSaving.XHostT.K.DS5H

/-- **the head bins of design B**: 80,640. -/
theorem headBins3945H_card : (headBins sched3945H).card = 80640 := by
  rw [headBins_eq sched3945H (resK ct resT 360 408) sched3945H_fits resK3945H_pos]
  exact sched3945H_heads

end OAI.PowerSaving.XHostT.Inst

#print axioms OAI.PowerSaving.XHostT.Inst.headBins3945H_card
