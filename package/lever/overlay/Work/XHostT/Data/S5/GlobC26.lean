import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 26 of S5: host kinds 3120..3239. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin26 : inR rt [0, 3, 5, 48] 3120 120 = true := by decide +kernel
theorem gc0_26 : cntR rt 0 3120 120 = 20832 := by decide +kernel
theorem gc3_26 : cntR rt 3 3120 120 = 2912 := by decide +kernel
theorem gc5_26 : cntR rt 5 3120 120 = 0 := by decide +kernel
theorem gc48_26 : cntR rt 48 3120 120 = 3136 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
