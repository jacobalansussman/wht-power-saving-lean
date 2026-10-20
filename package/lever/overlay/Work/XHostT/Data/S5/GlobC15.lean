import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 15 of S5: host kinds 1800..1919. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin15 : inR rt [0, 3, 5, 48] 1800 120 = true := by decide +kernel
theorem gc0_15 : cntR rt 0 1800 120 = 21056 := by decide +kernel
theorem gc3_15 : cntR rt 3 1800 120 = 224 := by decide +kernel
theorem gc5_15 : cntR rt 5 1800 120 = 0 := by decide +kernel
theorem gc48_15 : cntR rt 48 1800 120 = 5600 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
