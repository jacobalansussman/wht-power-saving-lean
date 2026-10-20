import Work.XHostT.KHistC
import Work.XHostT.Data.S5H.Calls
/-! p-hostck mk_parts.py: global chunk 7 of S5H: host kinds 840..959. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
theorem gin7 : inR rt [0, 3, 5, 48] 840 120 = true := by decide +kernel
theorem gc0_7 : cntR rt 0 840 120 = 20832 := by decide +kernel
theorem gc3_7 : cntR rt 3 840 120 = 1120 := by decide +kernel
theorem gc5_7 : cntR rt 5 840 120 = 0 := by decide +kernel
theorem gc48_7 : cntR rt 48 840 120 = 4928 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
