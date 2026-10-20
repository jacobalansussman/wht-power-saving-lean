import Work.XHostT.KHistC
import Work.XHostT.Data.H2979.Calls
/-! p-hostck mk_parts.py NEGATIVE control: one chunk's remainder count one too high; this module must FAIL to compile. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
theorem gbad : cntR rt 1 0 111 = 3505 := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979
