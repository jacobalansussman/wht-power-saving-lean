import Work.XHostT.KHeadC
import Work.XHostT.Data.H2979Fc.Calls
/-! afc-final (wht45/afc) Lean control, compile form: on A-fc's host table (#2979's words with the re-solved fill
u1-2979sep, afc-data's data) the kernel decides that
the number of A4 bins with remainder 48 is NOT 16,577 (one head more), and that NOT every residual is at least 5 (the
smallest width is 4).  These are the negations of the claims one step too far; `Work/XHostT/Data/H2979Fc/Heads.lean`
proves 16,576 and 4.  aA-final's Heads2979HNeg.lean with the A-fc data.  No `sorry`. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979Fc
open OAI.PowerSaving.XHostT.K
theorem heads_cnt_ne_succ : cntK rt 48 4 74 ≠ 16577 := by decide +kernel
theorem res_ge0_5_false : resGE resT 0 5 408 = false := by decide +kernel
theorem res_ge1_5_false : resGE resT 1 5 408 = false := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979Fc
#print axioms OAI.PowerSaving.XHostT.K.DH2979Fc.heads_cnt_ne_succ
#print axioms OAI.PowerSaving.XHostT.K.DH2979Fc.res_ge0_5_false
#print axioms OAI.PowerSaving.XHostT.K.DH2979Fc.res_ge1_5_false
