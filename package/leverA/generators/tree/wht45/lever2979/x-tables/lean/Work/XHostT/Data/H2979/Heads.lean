import Work.XHostT.KHeadC
import Work.XHostT.Data.H2979.Calls
/-! x-tables (wht45/lever2979): per-kind data facts of #2979's host table H2979 (decided by the kernel):
the number of A4 bins (kind 4) with remainder 48 (the literal heads: 8,064 of child 8 + 8,022 of child 10), and the
residual bounds of the two providers (every residual is at least 4; the two vectors are equal). -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
/-- the A4 bins with remainder 48 (all of them unhosted: the 16,086 literal heads). -/
theorem heads_cnt : cntK rt 48 4 74 = 16086 := by decide +kernel
/-- the F1 bins with remainder 48 (unhosted, not heads: 6,512 + 5,592). -/
theorem f1_cnt : cntK rt 48 6 74 = 12104 := by decide +kernel
theorem res_ge0 : resGE resT 0 4 408 = true := by decide +kernel
theorem res_ge1 : resGE resT 1 4 408 = true := by decide +kernel
theorem res_eq01 : resEqP resT 408 = true := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979
#print axioms OAI.PowerSaving.XHostT.K.DH2979.heads_cnt
