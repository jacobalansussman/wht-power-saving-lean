import Work.XHostT.KHeadC
import Work.XHostT.Data.H2979.Calls
/-! x-tables NEGATIVE control: the head count one too high; this module must FAIL to compile. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
theorem heads_bad : cntK rt 48 4 74 = 16087 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979
