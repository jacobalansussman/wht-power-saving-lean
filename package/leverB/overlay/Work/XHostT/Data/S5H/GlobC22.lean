import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 22 of S5H: host kinds 2640..2759. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin22 : inR rt [0, 3, 5, 48] 2640 120 = true := by decide +kernel
theorem gc0_22 : cntR rt 0 2640 120 = 20832 := by decide +kernel
theorem gc3_22 : cntR rt 3 2640 120 = 0 := by decide +kernel
theorem gc5_22 : cntR rt 5 2640 120 = 0 := by decide +kernel
theorem gc48_22 : cntR rt 48 2640 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
