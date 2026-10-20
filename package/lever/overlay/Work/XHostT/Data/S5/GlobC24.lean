import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 24 of S5: host kinds 2880..2999. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin24 : inR rt [0, 3, 5, 48] 2880 120 = true := by decide +kernel
theorem gc0_24 : cntR rt 0 2880 120 = 21056 := by decide +kernel
theorem gc3_24 : cntR rt 3 2880 120 = 0 := by decide +kernel
theorem gc5_24 : cntR rt 5 2880 120 = 0 := by decide +kernel
theorem gc48_24 : cntR rt 48 2880 120 = 5824 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
