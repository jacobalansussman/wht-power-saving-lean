import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 17 of S5: host kinds 2040..2159. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin17 : inR rt [0, 3, 5, 48] 2040 120 = true := by decide +kernel
theorem gc0_17 : cntR rt 0 2040 120 = 20832 := by decide +kernel
theorem gc3_17 : cntR rt 3 2040 120 = 0 := by decide +kernel
theorem gc5_17 : cntR rt 5 2040 120 = 0 := by decide +kernel
theorem gc48_17 : cntR rt 48 2040 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
