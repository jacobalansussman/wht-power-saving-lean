import Work.XHostT.KHistC
import Work.XHostT.Data.H2979.Calls
/-! p-hostck mk_parts.py: global chunk 5 of H2979: host kinds 555..665. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
theorem gin5 : inR rt [0, 1, 3, 5, 48] 555 111 = true := by decide +kernel
theorem gc0_5 : cntR rt 0 555 111 = 16695 := by decide +kernel
theorem gc1_5 : cntR rt 1 555 111 = 2570 := by decide +kernel
theorem gc3_5 : cntR rt 3 555 111 = 1561 := by decide +kernel
theorem gc5_5 : cntR rt 5 555 111 = 1 := by decide +kernel
theorem gc48_5 : cntR rt 48 555 111 = 4037 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979
