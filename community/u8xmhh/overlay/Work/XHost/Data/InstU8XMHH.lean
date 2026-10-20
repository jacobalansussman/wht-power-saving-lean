import Work.GCert.Data.InstU8XM
import Work.XHost.Data.RateU8XMHHx
namespace OAI.PowerSaving.BridgeRate.U8XMHH
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
/-- the remaining climbs of ONE run (2,016 hosts, ranks `rem`). -/
def remPrice (φ : ℕ → ℝ) : ℝ := 290 * bcost φ 1 + 110 * bcost φ 2 + 32 * bcost φ 3 + 112 * bcost φ 4 + 8 * bcost φ 12 + 326 * bcost φ 48
def padPrice (φ : ℕ → ℝ) : ℝ := 1 * bcost φ 1 + 1 * bcost φ 61
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * BridgeRate.U8XM.invPrice φ + 8 * remPrice φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocksU8XMHHx := by
  unfold unitPrice BridgeRate.U8XM.invPrice remPrice padPrice
  simp only [listCost, BlockAccounting.blocksU8XMHHx, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
end
end OAI.PowerSaving.BridgeRate.U8XMHH
#print axioms OAI.PowerSaving.BridgeRate.U8XMHH.unit_cost_x
