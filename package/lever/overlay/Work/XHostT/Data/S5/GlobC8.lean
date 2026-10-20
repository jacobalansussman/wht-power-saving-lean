import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 8 of S5: host kinds 960..1079. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin8 : inR rt [0, 3, 5, 48] 960 120 = true := by decide +kernel
theorem gc0_8 : cntR rt 0 960 120 = 20832 := by decide +kernel
theorem gc3_8 : cntR rt 3 960 120 = 0 := by decide +kernel
theorem gc5_8 : cntR rt 5 960 120 = 0 := by decide +kernel
theorem gc48_8 : cntR rt 48 960 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
