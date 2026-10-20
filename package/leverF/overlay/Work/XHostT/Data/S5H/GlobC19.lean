import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 19 of S5H: host kinds 2280..2399. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin19 : inR rt [0, 3, 5, 48] 2280 120 = true := by decide +kernel
theorem gc0_19 : cntR rt 0 2280 120 = 20832 := by decide +kernel
theorem gc3_19 : cntR rt 3 2280 120 = 448 := by decide +kernel
theorem gc5_19 : cntR rt 5 2280 120 = 0 := by decide +kernel
theorem gc48_19 : cntR rt 48 2280 120 = 5600 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
