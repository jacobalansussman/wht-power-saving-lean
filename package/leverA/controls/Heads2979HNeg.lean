import Work.XHostT.KHeadC
import Work.XHostT.Data.H2979H.Calls
/-! aA-final (wht45/allheads) Lean control, compile form: on design A's host table (#2979 with every head kept,
aA-data's data) the kernel decides that
the number of A4 bins with remainder 48 is NOT 16,577 (one head more), and that NOT every residual is at least 5 (the
smallest width is 4).  These are the negations of the claims one step too far; `Work/XHostT/Data/H2979H/Heads.lean`
proves 16,576 and 4.  x-final's Heads2979Neg.lean with the design-A data.  No `sorry`. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979H
open OAI.PowerSaving.XHostT.K
theorem heads_cnt_ne_succ : cntK rt 48 4 74 ≠ 16577 := by decide +kernel
theorem res_ge0_5_false : resGE resT 0 5 408 = false := by decide +kernel
theorem res_ge1_5_false : resGE resT 1 5 408 = false := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979H
#print axioms OAI.PowerSaving.XHostT.K.DH2979H.heads_cnt_ne_succ
#print axioms OAI.PowerSaving.XHostT.K.DH2979H.res_ge0_5_false
#print axioms OAI.PowerSaving.XHostT.K.DH2979H.res_ge1_5_false
