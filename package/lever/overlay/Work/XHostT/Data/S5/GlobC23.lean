import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 23 of S5: host kinds 2760..2879. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin23 : inR rt [0, 3, 5, 48] 2760 120 = true := by decide +kernel
theorem gc0_23 : cntR rt 0 2760 120 = 20384 := by decide +kernel
theorem gc3_23 : cntR rt 3 2760 120 = 1792 := by decide +kernel
theorem gc5_23 : cntR rt 5 2760 120 = 448 := by decide +kernel
theorem gc48_23 : cntR rt 48 2760 120 = 4256 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
