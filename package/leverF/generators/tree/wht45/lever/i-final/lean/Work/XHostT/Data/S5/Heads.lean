import Work.XHostT.KHeadC
import Work.XHostT.Data.S5.Calls
/-! i-final (wht45/lever): per-kind data facts of #3945's host table S5 (p-hostck's data), decided by the kernel with
x-tables' counters (Work/XHostT/KHeadC.lean): the number of A4 bins (kind 4) with remainder 48 (the 57,312 literal
heads: 28,672 of child 8 + 28,640 of child 10), the F1 bins with remainder 48 (80,640, not heads), and the residual
bounds of the two providers (every residual is at least 5; the two vectors are equal).  Same text as x-tables'
Data/H2979/Heads.lean with the instance changed.  No `sorry`. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
/-- the A4 bins with remainder 48 (all of them unhosted: the 57,312 literal heads). -/
theorem heads_cnt : cntK rt 48 4 360 = 57312 := by decide +kernel
/-- the F1 bins with remainder 48 (unhosted, not heads). -/
theorem f1_cnt : cntK rt 48 6 360 = 80640 := by decide +kernel
theorem res_ge0 : resGE resT 0 5 408 = true := by decide +kernel
theorem res_ge1 : resGE resT 1 5 408 = true := by decide +kernel
theorem res_eq01 : resEqP resT 408 = true := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
#print axioms OAI.PowerSaving.XHostT.K.DS5.heads_cnt
#print axioms OAI.PowerSaving.XHostT.K.DS5.f1_cnt
#print axioms OAI.PowerSaving.XHostT.K.DS5.res_eq01
