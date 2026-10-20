import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 11 of S5: host kinds 1320..1439. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin11 : inR rt [0, 3, 5, 48] 1320 120 = true := by decide +kernel
theorem gc0_11 : cntR rt 0 1320 120 = 20832 := by decide +kernel
theorem gc3_11 : cntR rt 3 1320 120 = 0 := by decide +kernel
theorem gc5_11 : cntR rt 5 1320 120 = 0 := by decide +kernel
theorem gc48_11 : cntR rt 48 1320 120 = 6048 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
