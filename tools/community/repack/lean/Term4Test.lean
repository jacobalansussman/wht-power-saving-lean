import Work.BlockAccounting.Moment4
namespace OAI.PowerSaving.BlockAccounting
/-! x-pack: drop-in test of term_le4 on rows b0 (rank 1, the largest), b9 (pad 37) and b11/b12 (idle terms) of the generated quartic certificate. Not part of any certificate. -/
example : T 73 (1 - 12383292/(10:ℝ)^10) 1 187636 ≤ ((2584048812127:ℝ)/1000000000) :=
  term_le4 73 1 187636 (12383292/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
      (by norm_num [Finset.sum_range_succ, Nat.factorial]))
    (by norm_num) (by norm_num)
example : T 73 (1 - 12383292/(10:ℝ)^10) 37 1 ≤ ((507276007:ℝ)/1000000000) :=
  term_le4 73 37 1 (12383292/(10:ℝ)^10) ((339770764253:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (log_le_of_le_expsum _ _ 13 (by norm_num) (by norm_num)
      (by norm_num [Finset.sum_range_succ, Nat.factorial]))
    (by norm_num) (by norm_num)
example : T 73 (1 - 12383292/(10:ℝ)^10) 72 1 ≤ ((986318217:ℝ)/1000000000) :=
  term_le4 73 72 1 (12383292/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
      (by norm_num [Finset.sum_range_succ, Nat.factorial]))
    (by norm_num) (by norm_num)
example : T 73 (1 - 12383292/(10:ℝ)^10) 1 1 ≤ ((2754321:ℝ)/200000000) :=
  term_le4 73 1 1 (12383292/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
      (by norm_num [Finset.sum_range_succ, Nat.factorial]))
    (by norm_num) (by norm_num)
end OAI.PowerSaving.BlockAccounting
