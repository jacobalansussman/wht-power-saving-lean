import Work.BlockAccounting.Moment
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of PR588 atm3 WITHOUT source copies (mutable sources), #401 banks for the helpers that are not hosted, same-run idle-axis hosting (DaysSky #442) of the others; h-L3-build. -/
def blocksC588MHHx : List (ℕ × ℕ) := [(1, 180640), (2, 35424), (3, 20224), (4, 17056), (5, 5728), (6, 4032), (7, 320), (8, 1584), (39, 1), (48, 1680)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksC588MHHx_rank : (blocksC588MHHx.map fun b => b.1 * b.2).sum = 7279 * 73 - 2560 := by
  norm_num [blocksC588MHHx]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksC588MHHx_tally : (blocksC588MHHx.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((7279:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksC588MHHx]

theorem blocksC588MHHx_count : (blocksC588MHHx.map fun b => b.2).sum = 266689 := by
  norm_num [blocksC588MHHx]

theorem blocksC588MHHx_proper : ∀ b ∈ blocksC588MHHx, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksC588MHHx` at exponent `1 - 15699306/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 7279).  Relative slack of the rational certificate: 9.86e-11. -/
theorem fold_C588MHHx :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 15699306/(10:ℝ)^10) blocksC588MHHx
      + ((7279:ℕ):ℝ) * (T 73 (1 - 15699306/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 15699306/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((7279:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 15699306/(10:ℝ)^10) blocksC588MHHx
      = T 73 (1 - 15699306/(10:ℝ)^10) 1 180640 + (T 73 (1 - 15699306/(10:ℝ)^10) 2 35424 + (T 73 (1 - 15699306/(10:ℝ)^10) 3 20224 + (T 73 (1 - 15699306/(10:ℝ)^10) 4 17056 + (T 73 (1 - 15699306/(10:ℝ)^10) 5 5728 + (T 73 (1 - 15699306/(10:ℝ)^10) 6 4032 + (T 73 (1 - 15699306/(10:ℝ)^10) 7 320 + (T 73 (1 - 15699306/(10:ℝ)^10) 8 1584 + (T 73 (1 - 15699306/(10:ℝ)^10) 39 1 + (T 73 (1 - 15699306/(10:ℝ)^10) 48 1680 + 0))))))))) := rfl
  have b0 : T 73 (1 - 15699306/(10:ℝ)^10) 1 180640 ≤ ((2491244536823:ℝ)/1000000000) :=
    term_le 73 1 180640 (15699306/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 15699306/(10:ℝ)^10) 2 35424 ≤ ((195203421693:ℝ)/200000000) :=
    term_le 73 2 35424 (15699306/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 15699306/(10:ℝ)^10) 3 20224 ≤ ((835298486563:ℝ)/1000000000) :=
    term_le 73 3 20224 (15699306/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 15699306/(10:ℝ)^10) 4 17056 ≤ ((234711530111:ℝ)/250000000) :=
    term_le 73 4 17056 (15699306/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 15699306/(10:ℝ)^10) 5 5728 ≤ ((393983567579:ℝ)/1000000000) :=
    term_le 73 5 5728 (15699306/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 15699306/(10:ℝ)^10) 6 4032 ≤ ((166349907457:ℝ)/500000000) :=
    term_le 73 6 4032 (15699306/(10:ℝ)^10) ((1249349985961:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 15699306/(10:ℝ)^10) 7 320 ≤ ((15399042091:ℝ)/500000000) :=
    term_le 73 7 320 (15699306/(10:ℝ)^10) ((468909858419:ℝ)/200000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 15699306/(10:ℝ)^10) 8 1584 ≤ ((43548160289:ℝ)/250000000) :=
    term_le 73 8 1584 (15699306/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 15699306/(10:ℝ)^10) 39 1 ≤ ((534772633:ℝ)/1000000000) :=
    term_le 73 39 1 (15699306/(10:ℝ)^10) ((31344889751:ℝ)/50000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 13 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 15699306/(10:ℝ)^10) 48 1680 ≤ ((1105384866529:ℝ)/1000000000) :=
    term_le 73 48 1680 (15699306/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 15699306/(10:ℝ)^10) 72 1 ≤ ((986322729:ℝ)/1000000000) :=
    term_le 73 72 1 (15699306/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b11 : T 73 (1 - 15699306/(10:ℝ)^10) 1 1 ≤ ((13791213:ℝ)/1000000000) :=
    term_le 73 1 1 (15699306/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((7279:ℕ):ℝ) = 7279 := by norm_num
  have e4 : T 73 (1 - 15699306/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 15699306/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_C588MHHx
