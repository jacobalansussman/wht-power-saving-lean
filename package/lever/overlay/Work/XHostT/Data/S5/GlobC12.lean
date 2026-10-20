import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 12 of S5: host kinds 1440..1559. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin12 : inR rt [0, 3, 5, 48] 1440 120 = true := by decide +kernel
theorem gc0_12 : cntR rt 0 1440 120 = 21056 := by decide +kernel
theorem gc3_12 : cntR rt 3 1440 120 = 1344 := by decide +kernel
theorem gc5_12 : cntR rt 5 1440 120 = 0 := by decide +kernel
theorem gc48_12 : cntR rt 48 1440 120 = 4480 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
