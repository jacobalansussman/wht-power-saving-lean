import Work.XHostT.KHeadC
import Work.XHostT.Data.S5.Calls
/-! i-final (wht45/lever) Lean control, compile form: on #3945's host table (p-hostck's data) the kernel decides that the
number of A4 bins with remainder 48 is NOT 57,313 (one head more), and that NOT every residual is at least 6 (the
smallest width is 5).  These are the negations of the claims one step too far; `Work/XHostT/Data/S5/Heads.lean` proves
57,312 and 5.  No `sorry`. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem heads_cnt_ne_succ : cntK rt 48 4 360 ≠ 57313 := by decide +kernel
theorem res_ge0_6_false : resGE resT 0 6 408 = false := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
#print axioms OAI.PowerSaving.XHostT.K.DS5.heads_cnt_ne_succ
#print axioms OAI.PowerSaving.XHostT.K.DS5.res_ge0_6_false
