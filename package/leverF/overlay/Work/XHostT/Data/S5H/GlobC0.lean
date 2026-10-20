import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 0 of S5H: host kinds 0..119. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin0 : inR rt [0, 3, 5, 48] 0 120 = true := by decide +kernel
theorem gc0_0 : cntR rt 0 0 120 = 20384 := by decide +kernel
theorem gc3_0 : cntR rt 3 0 120 = 672 := by decide +kernel
theorem gc5_0 : cntR rt 5 0 120 = 672 := by decide +kernel
theorem gc48_0 : cntR rt 48 0 120 = 5152 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
