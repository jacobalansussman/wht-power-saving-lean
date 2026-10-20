import Work.XHostT.KHistC
import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: global chunk 16 of S5: host kinds 1920..2039. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem gin16 : inR rt [0, 3, 5, 48] 1920 120 = true := by decide +kernel
theorem gc0_16 : cntR rt 0 1920 120 = 20384 := by decide +kernel
theorem gc3_16 : cntR rt 3 1920 120 = 2688 := by decide +kernel
theorem gc5_16 : cntR rt 5 1920 120 = 448 := by decide +kernel
theorem gc48_16 : cntR rt 48 1920 120 = 3360 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
