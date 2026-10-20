import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 2 of S5H: host kinds 240..359. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin2 : inR rt [0, 3, 5, 48] 240 120 = true := by decide +kernel
theorem gc0_2 : cntR rt 0 240 120 = 20788 := by decide +kernel
theorem gc3_2 : cntR rt 3 240 120 = 224 := by decide +kernel
theorem gc5_2 : cntR rt 5 240 120 = 44 := by decide +kernel
theorem gc48_2 : cntR rt 48 240 120 = 5824 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
