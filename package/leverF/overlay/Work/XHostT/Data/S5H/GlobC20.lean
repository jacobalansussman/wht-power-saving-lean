import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 20 of S5H: host kinds 2400..2519. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin20 : inR rt [0, 3, 5, 48] 2400 120 = true := by decide +kernel
theorem gc0_20 : cntR rt 0 2400 120 = 20384 := by decide +kernel
theorem gc3_20 : cntR rt 3 2400 120 = 0 := by decide +kernel
theorem gc5_20 : cntR rt 5 2400 120 = 448 := by decide +kernel
theorem gc48_20 : cntR rt 48 2400 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
