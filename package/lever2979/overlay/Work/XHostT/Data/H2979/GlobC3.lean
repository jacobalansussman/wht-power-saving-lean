import Work.XHostT.KHistC
import Work.XHostT.Data.H2979.Calls
/-! p-hostck mk_parts.py: global chunk 3 of H2979: host kinds 333..443. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
theorem gin3 : inR rt [0, 1, 3, 5, 48] 333 111 = true := by decide +kernel
theorem gc0_3 : cntR rt 0 333 111 = 15158 := by decide +kernel
theorem gc1_3 : cntR rt 1 333 111 = 3504 := by decide +kernel
theorem gc3_3 : cntR rt 3 333 111 = 384 := by decide +kernel
theorem gc5_3 : cntR rt 5 333 111 = 826 := by decide +kernel
theorem gc48_3 : cntR rt 48 333 111 = 4992 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979
