import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 1 of S5: host kinds 120..239. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin1 : inR rt [0, 3, 5, 48] 120 120 = true := by decide +kernel
theorem gc0_1 : cntR rt 0 120 120 = 20832 := by decide +kernel
theorem gc3_1 : cntR rt 3 120 120 = 0 := by decide +kernel
theorem gc5_1 : cntR rt 5 120 120 = 0 := by decide +kernel
theorem gc48_1 : cntR rt 48 120 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
