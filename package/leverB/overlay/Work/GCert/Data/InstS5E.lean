import Work.Bridge377.Bank
import Work.FoldRate.Group
import Work.GCert.Data.RateS5Ex
namespace OAI.PowerSaving.BridgeRate.S5E
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
def invPrice (φ : ℕ → ℝ) : ℝ := 2482 * φ 1 + 712 * φ 2 + 218 * φ 3 + 109 * φ 4 + 20 * φ 5 + 8 * φ 8
def padPrice (φ : ℕ → ℝ) : ℝ := 1 * bcost φ 19
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * invPrice φ + ((112:ℕ):ℝ) * GX.bankPrice377 φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocksS5Ex := by
  unfold unitPrice invPrice padPrice
  simp only [GX.bankPrice377, listCost, BlockAccounting.blocksS5Ex, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
theorem unit_tally : unitPrice (fun r => (r:ℝ)) = ((10427:ℕ):ℝ) * 73 - 2560 := by
  rw [unit_cost_x]
  exact BlockAccounting.blocksS5Ex_tally
end
end OAI.PowerSaving.BridgeRate.S5E
#print axioms OAI.PowerSaving.BridgeRate.S5E.unit_cost_x
