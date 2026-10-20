import Work.XHostT.KHeadC
import Work.XHostT.Data.S5.Calls
/-! i-final NEGATIVE control (must FAIL to compile): a residual bound one too high (#3945's smallest width is 5). -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem res_ge0_6 : resGE resT 0 6 408 = true := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
