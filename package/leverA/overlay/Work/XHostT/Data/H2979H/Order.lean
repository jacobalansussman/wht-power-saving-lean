import Work.XHostT.Data.H2979H.Calls
/-! p-hostck mk_parts.py: the global order of H2979H is valid (times distinct, below NC, increasing along each replica). -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979H
open OAI.PowerSaving.XHostT.K
theorem order_ok : chkOrder ct NC NC 0 0 = true := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979H
