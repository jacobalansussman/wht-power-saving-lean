import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 9 of S5: host kinds 1080..1199. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin9 : inR rt [0, 3, 5, 48] 1080 120 = true := by decide +kernel
theorem gc0_9 : cntR rt 0 1080 120 = 21056 := by decide +kernel
theorem gc3_9 : cntR rt 3 1080 120 = 2464 := by decide +kernel
theorem gc5_9 : cntR rt 5 1080 120 = 0 := by decide +kernel
theorem gc48_9 : cntR rt 48 1080 120 = 3360 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
