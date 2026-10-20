import Work.GCert.Data.InstU8X0
import Work.XHost.Data.RateU8X0HQx
namespace OAI.PowerSaving.BridgeRate.U8X0HQ
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
/-- the remaining climbs of ONE run (2,016 hosts, ranks `rem`). -/
def remPrice (φ : ℕ → ℝ) : ℝ := 290 * bcost φ 1 + 108 * bcost φ 2 + 32 * bcost φ 3 + 112 * bcost φ 4 + 12 * bcost φ 12 + 336 * bcost φ 48
def padPrice (φ : ℕ → ℝ) : ℝ := 1 * bcost φ 59
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * BridgeRate.U8X0.invPrice φ + 8 * remPrice φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocksU8X0HQx := by
  unfold unitPrice BridgeRate.U8X0.invPrice remPrice padPrice
  simp only [listCost, BlockAccounting.blocksU8X0HQx, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
end
end OAI.PowerSaving.BridgeRate.U8X0HQ
#print axioms OAI.PowerSaving.BridgeRate.U8X0HQ.unit_cost_x
