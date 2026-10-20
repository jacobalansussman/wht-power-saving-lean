import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 5 of S5: host kinds 600..719. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin5 : inR rt [0, 3, 5, 48] 600 120 = true := by decide +kernel
theorem gc0_5 : cntR rt 0 600 120 = 20832 := by decide +kernel
theorem gc3_5 : cntR rt 3 600 120 = 480 := by decide +kernel
theorem gc5_5 : cntR rt 5 600 120 = 0 := by decide +kernel
theorem gc48_5 : cntR rt 48 600 120 = 5568 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
