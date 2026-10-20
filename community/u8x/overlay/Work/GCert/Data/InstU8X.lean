import Work.Bridge377.Bank
import Work.FoldRate.Group
import Work.GCert.Data.RateU8Xx
namespace OAI.PowerSaving.BridgeRate.U8X
open Finset Binary Matrix RAM CB FoldRate
noncomputable section
def invPrice (φ : ℕ → ℝ) : ℝ := 2660 * φ 1 + 338 * φ 2 + 446 * φ 3 + 162 * φ 4 + 72 * φ 5 + 120 * φ 8
def padPrice (φ : ℕ → ℝ) : ℝ := 1 * bcost φ 41
def unitPrice (φ : ℕ → ℝ) : ℝ := 72 * invPrice φ + ((112:ℕ):ℝ) * GX.bankPrice377 φ + padPrice φ
theorem unit_cost_x (φ : ℕ → ℝ) : unitPrice φ = listCost φ BlockAccounting.blocksU8Xx := by
  unfold unitPrice invPrice padPrice
  simp only [GX.bankPrice377, listCost, BlockAccounting.blocksU8Xx, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  norm_num [bcost]
  try ring
theorem unit_tally : unitPrice (fun r => (r:ℝ)) = ((11889:ℕ):ℝ) * 73 - 2560 := by
  rw [unit_cost_x]
  exact BlockAccounting.blocksU8Xx_tally
end
end OAI.PowerSaving.BridgeRate.U8X
#print axioms OAI.PowerSaving.BridgeRate.U8X.unit_cost_x
