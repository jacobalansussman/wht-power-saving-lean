import Work.Bridge377.Bank
import Work.FoldRate.Group
import Work.GCert.Data.RateC588Mx
namespace OAI.PowerSaving.BridgeRate.C588M
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
def invPrice (φ : ℕ → ℝ) : ℝ := 2484 * φ 1 + 492 * φ 2 + 280 * φ 3 + 232 * φ 4 + 56 * φ 5 + 56 * φ 6 + 8 * φ 8
def padPrice (φ : ℕ → ℝ) : ℝ := 72 * bcost φ 28
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * invPrice φ + ((112:ℕ):ℝ) * GX.bankPrice377 φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocksC588Mx := by
  unfold unitPrice invPrice padPrice
  simp only [GX.bankPrice377, listCost, BlockAccounting.blocksC588Mx, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
theorem unit_tally : unitPrice (fun r => (r:ℝ)) = ((11200:ℕ):ℝ) * 73 - 2560 := by
  rw [unit_cost_x]
  exact BlockAccounting.blocksC588Mx_tally
end
end OAI.PowerSaving.BridgeRate.C588M
#print axioms OAI.PowerSaving.BridgeRate.C588M.unit_cost_x
