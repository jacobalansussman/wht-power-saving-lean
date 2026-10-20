import Work.XHostT.KHeadC
import Work.XHostT.Data.H2979.Calls
/-! x-final (wht45/lever2979) Lean control, compile form: on #2979's host table (x-tables' data) the kernel decides that
the number of A4 bins with remainder 48 is NOT 16,087 (one head more), and that NOT every residual is at least 5 (the
smallest width is 4).  These are the negations of the claims one step too far; `Work/XHostT/Data/H2979/Heads.lean`
proves 16,086 and 4.  i-final's Heads3945Neg.lean with the #2979 data.  No `sorry`. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
theorem heads_cnt_ne_succ : cntK rt 48 4 74 ≠ 16087 := by decide +kernel
theorem res_ge0_5_false : resGE resT 0 5 408 = false := by decide +kernel
theorem res_ge1_5_false : resGE resT 1 5 408 = false := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979
#print axioms OAI.PowerSaving.XHostT.K.DH2979.heads_cnt_ne_succ
#print axioms OAI.PowerSaving.XHostT.K.DH2979.res_ge0_5_false
#print axioms OAI.PowerSaving.XHostT.K.DH2979.res_ge1_5_false
