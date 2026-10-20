import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 21 of S5H: host kinds 2520..2639. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin21 : inR rt [0, 3, 5, 48] 2520 120 = true := by decide +kernel
theorem gc0_21 : cntR rt 0 2520 120 = 21056 := by decide +kernel
theorem gc3_21 : cntR rt 3 2520 120 = 2240 := by decide +kernel
theorem gc5_21 : cntR rt 5 2520 120 = 0 := by decide +kernel
theorem gc48_21 : cntR rt 48 2520 120 = 3584 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
