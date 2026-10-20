import Work.XHostT.KHeadC
import Work.XHostT.Data.S5.Calls
/-! i-final NEGATIVE control (must FAIL to compile): one head more than #3945's table has. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem heads_cnt_plus1 : cntK rt 48 4 360 = 57313 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
