import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 25 of S5: host kinds 3000..3119. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin25 : inR rt [0, 3, 5, 48] 3000 120 = true := by decide +kernel
theorem gc0_25 : cntR rt 0 3000 120 = 20832 := by decide +kernel
theorem gc3_25 : cntR rt 3 3000 120 = 0 := by decide +kernel
theorem gc5_25 : cntR rt 5 3000 120 = 0 := by decide +kernel
theorem gc48_25 : cntR rt 48 3000 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
