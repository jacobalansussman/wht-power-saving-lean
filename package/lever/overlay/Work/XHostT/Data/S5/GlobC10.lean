import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 10 of S5: host kinds 1200..1319. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin10 : inR rt [0, 3, 5, 48] 1200 120 = true := by decide +kernel
theorem gc0_10 : cntR rt 0 1200 120 = 20320 := by decide +kernel
theorem gc3_10 : cntR rt 3 1200 120 = 448 := by decide +kernel
theorem gc5_10 : cntR rt 5 1200 120 = 512 := by decide +kernel
theorem gc48_10 : cntR rt 48 1200 120 = 5600 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
