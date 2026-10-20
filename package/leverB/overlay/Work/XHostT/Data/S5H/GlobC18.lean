import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 18 of S5H: host kinds 2160..2279. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin18 : inR rt [0, 3, 5, 48] 2160 120 = true := by decide +kernel
theorem gc0_18 : cntR rt 0 2160 120 = 21056 := by decide +kernel
theorem gc3_18 : cntR rt 3 2160 120 = 1120 := by decide +kernel
theorem gc5_18 : cntR rt 5 2160 120 = 0 := by decide +kernel
theorem gc48_18 : cntR rt 48 2160 120 = 4704 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
