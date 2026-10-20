import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 6 of S5H: host kinds 720..839. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin6 : inR rt [0, 3, 5, 48] 720 120 = true := by decide +kernel
theorem gc0_6 : cntR rt 0 720 120 = 20384 := by decide +kernel
theorem gc3_6 : cntR rt 3 720 120 = 672 := by decide +kernel
theorem gc5_6 : cntR rt 5 720 120 = 672 := by decide +kernel
theorem gc48_6 : cntR rt 48 720 120 = 5152 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
