import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 4 of S5: host kinds 480..599. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin4 : inR rt [0, 3, 5, 48] 480 120 = true := by decide +kernel
theorem gc0_4 : cntR rt 0 480 120 = 20832 := by decide +kernel
theorem gc3_4 : cntR rt 3 480 120 = 0 := by decide +kernel
theorem gc5_4 : cntR rt 5 480 120 = 0 := by decide +kernel
theorem gc48_4 : cntR rt 48 480 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
