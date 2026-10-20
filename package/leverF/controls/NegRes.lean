import Work.XHostT.KHeadC
import Work.XHostT.Data.S5H.Calls
/-! bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head): NEGATIVE control (must FAIL to compile): a residual bound one too high (the smallest width is 5). -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem res_ge0_6 : resGE resT 0 6 408 = true := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
