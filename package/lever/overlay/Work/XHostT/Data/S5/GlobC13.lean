import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 13 of S5: host kinds 1560..1679. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin13 : inR rt [0, 3, 5, 48] 1560 120 = true := by decide +kernel
theorem gc0_13 : cntR rt 0 1560 120 = 20384 := by decide +kernel
theorem gc3_13 : cntR rt 3 1560 120 = 1568 := by decide +kernel
theorem gc5_13 : cntR rt 5 1560 120 = 448 := by decide +kernel
theorem gc48_13 : cntR rt 48 1560 120 = 4480 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
