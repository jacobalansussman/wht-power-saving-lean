import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 3 of S5: host kinds 360..479. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin3 : inR rt [0, 3, 5, 48] 360 120 = true := by decide +kernel
theorem gc0_3 : cntR rt 0 360 120 = 20384 := by decide +kernel
theorem gc3_3 : cntR rt 3 360 120 = 1120 := by decide +kernel
theorem gc5_3 : cntR rt 5 360 120 = 672 := by decide +kernel
theorem gc48_3 : cntR rt 48 360 120 = 4704 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
