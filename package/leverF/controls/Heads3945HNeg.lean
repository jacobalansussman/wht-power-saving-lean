import Work.XHostT.KHeadC
import Work.XHostT.Data.S5H.Calls
/-! bB-final (wht45/allheads, design B) Lean control, compile form (i-final's Heads3945Neg.lean on the design-B table):
on the host table of #3945 with every head kept (bB-data's S5H data) the kernel decides that the number of A4 bins with
remainder 48 is NOT 80,641 (one head more), and that NOT every residual is at least 6 (the smallest width is 5).  These
are the negations of the claims one step too far; `Work/XHostT/Data/S5H/Heads.lean` proves 80,640 and 5.  No `sorry`. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem heads_cnt_ne_succ : cntK rt 48 4 360 ≠ 80641 := by decide +kernel
theorem res_ge0_6_false : resGE resT 0 6 408 = false := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
#print axioms OAI.PowerSaving.XHostT.K.DS5H.heads_cnt_ne_succ
#print axioms OAI.PowerSaving.XHostT.K.DS5H.res_ge0_6_false
