import Work.XHostT.Data.S5.Calls
/-! p-hostck mk_parts.py: the global order of S5 is valid (times distinct, below NC, increasing along each replica). -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT.K
theorem order_ok : chkOrder ct NC NC 0 0 = true := by decide +kernel
end OAI.PowerSaving.XHostT.K.DS5
