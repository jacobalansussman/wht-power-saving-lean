import Work.XHostT.KHeadC
import Work.XHostT.Data.S5H.Calls
/-! bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head): NEGATIVE control (must FAIL to compile): the head count of #3945's OWN fill (57,312) on the S5H table. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem heads_cnt_old : cntK rt 48 4 360 = 57312 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
