import Work.XHostT.KHistC
import Work.XHostT.Data.H2979H.Calls
/-! p-hostck mk_parts.py: global chunk 1 of H2979H: host kinds 111..221. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979H
open OAI.PowerSaving.XHostT.K
theorem gin1 : inR rt [0, 1, 3, 5, 48] 111 111 = true := by decide +kernel
theorem gc0_1 : cntR rt 0 111 111 = 16382 := by decide +kernel
theorem gc1_1 : cntR rt 1 111 111 = 2883 := by decide +kernel
theorem gc3_1 : cntR rt 3 111 111 = 499 := by decide +kernel
theorem gc5_1 : cntR rt 5 111 111 = 0 := by decide +kernel
theorem gc48_1 : cntR rt 48 111 111 = 5100 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979H
