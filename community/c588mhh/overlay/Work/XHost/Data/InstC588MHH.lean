import Work.GCert.Data.InstC588M
import Work.XHost.Data.RateC588MHHx
namespace OAI.PowerSaving.BridgeRate.C588MHH
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
/-- the remaining climbs of ONE run (2,016 hosts, ranks `rem`). -/
def remPrice (φ : ℕ → ℝ) : ℝ := 224 * bcost φ 1 + 8 * bcost φ 3 + 44 * bcost φ 4 + 212 * bcost φ 5 + 40 * bcost φ 7 + 126 * bcost φ 8 + 210 * bcost φ 48
def padPrice (φ : ℕ → ℝ) : ℝ := 1 * bcost φ 39
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * BridgeRate.C588M.invPrice φ + 8 * remPrice φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocksC588MHHx := by
  unfold unitPrice BridgeRate.C588M.invPrice remPrice padPrice
  simp only [listCost, BlockAccounting.blocksC588MHHx, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
end
end OAI.PowerSaving.BridgeRate.C588MHH
#print axioms OAI.PowerSaving.BridgeRate.C588MHH.unit_cost_x
