import Work.XHostT.KHeadC
import Work.XHostT.Data.S5H.Calls
/-! bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head): per-kind data facts of the host table S5H (p-hostck's generator on the all-heads fill), decided by the kernel with
x-tables' counters (Work/XHostT/KHeadC.lean): the number of A4 bins (kind 4) with remainder 48 (the 80,640 literal
heads: 40,320 of child 8 + 40,320 of child 10), the F1 bins with remainder 48 (57,312, not heads), and the residual
bounds of the two providers (every residual is at least 5; the two vectors are equal).  i-final's Data/S5/Heads.lean
with the instance and the counts changed.  No `sorry`. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
/-- the A4 bins with remainder 48 (all of them unhosted: the 80,640 literal heads). -/
theorem heads_cnt : cntK rt 48 4 360 = 80640 := by decide +kernel
/-- the F1 bins with remainder 48 (unhosted, not heads). -/
theorem f1_cnt : cntK rt 48 6 360 = 57312 := by decide +kernel
theorem res_ge0 : resGE resT 0 5 408 = true := by decide +kernel
theorem res_ge1 : resGE resT 1 5 408 = true := by decide +kernel
theorem res_eq01 : resEqP resT 408 = true := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
#print axioms OAI.PowerSaving.XHostT.K.DS5H.heads_cnt
#print axioms OAI.PowerSaving.XHostT.K.DS5H.f1_cnt
#print axioms OAI.PowerSaving.XHostT.K.DS5H.res_eq01
