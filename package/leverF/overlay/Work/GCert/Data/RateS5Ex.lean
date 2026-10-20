import Work.BlockAccounting.Moment
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of the warped #3945 word S5E (gcert/1w, p-lower's lowering; 408 helpers, 91 zero buffers), #401 dedicated banks, and an explicit x-exit-words pad-minimal (remapped) bank schedule. -/
def blocksS5Ex : List (ℕ × ℕ) := [(1, 178704), (2, 51264), (3, 15696), (4, 7848), (5, 1440), (8, 5952), (16, 3584), (19, 1), (32, 3584), (48, 3584)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksS5Ex_rank : (blocksS5Ex.map fun b => b.1 * b.2).sum = 10427 * 73 - 2560 := by
  norm_num [blocksS5Ex]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksS5Ex_tally : (blocksS5Ex.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((10427:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksS5Ex]

theorem blocksS5Ex_count : (blocksS5Ex.map fun b => b.2).sum = 271657 := by
  norm_num [blocksS5Ex]

theorem blocksS5Ex_proper : ∀ b ∈ blocksS5Ex, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksS5Ex` at exponent `1 - 14548216/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 10427).  Relative slack of the rational certificate: 1.87e-10. -/
theorem fold_S5Ex :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 14548216/(10:ℝ)^10) blocksS5Ex
      + ((10427:ℕ):ℝ) * (T 73 (1 - 14548216/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 14548216/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((10427:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 14548216/(10:ℝ)^10) blocksS5Ex
      = T 73 (1 - 14548216/(10:ℝ)^10) 1 178704 + (T 73 (1 - 14548216/(10:ℝ)^10) 2 51264 + (T 73 (1 - 14548216/(10:ℝ)^10) 3 15696 + (T 73 (1 - 14548216/(10:ℝ)^10) 4 7848 + (T 73 (1 - 14548216/(10:ℝ)^10) 5 1440 + (T 73 (1 - 14548216/(10:ℝ)^10) 8 5952 + (T 73 (1 - 14548216/(10:ℝ)^10) 16 3584 + (T 73 (1 - 14548216/(10:ℝ)^10) 19 1 + (T 73 (1 - 14548216/(10:ℝ)^10) 32 3584 + (T 73 (1 - 14548216/(10:ℝ)^10) 48 3584 + 0))))))))) := rfl
  have b0 : T 73 (1 - 14548216/(10:ℝ)^10) 1 178704 ≤ ((123166393827:ℝ)/50000000) :=
    term_le 73 1 178704 (14548216/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 14548216/(10:ℝ)^10) 2 51264 ≤ ((705931385257:ℝ)/500000000) :=
    term_le 73 2 51264 (14548216/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 14548216/(10:ℝ)^10) 3 15696 ≤ ((324021678487:ℝ)/500000000) :=
    term_le 73 3 15696 (14548216/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 14548216/(10:ℝ)^10) 4 7848 ≤ ((215924063263:ℝ)/500000000) :=
    term_le 73 4 7848 (14548216/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 14548216/(10:ℝ)^10) 5 1440 ≤ ((19803117261:ℝ)/200000000) :=
    term_le 73 5 1440 (14548216/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 14548216/(10:ℝ)^10) 8 5952 ≤ ((327187740109:ℝ)/500000000) :=
    term_le 73 8 5952 (14548216/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 14548216/(10:ℝ)^10) 16 3584 ≤ ((787270804859:ℝ)/1000000000) :=
    term_le 73 16 3584 (14548216/(10:ℝ)^10) ((151787071891:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 17 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 14548216/(10:ℝ)^10) 19 1 ≤ ((130392073:ℝ)/500000000) :=
    term_le 73 19 1 (14548216/(10:ℝ)^10) ((1346020461983:ℝ)/1000000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 17 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 14548216/(10:ℝ)^10) 32 3584 ≤ ((78647731641:ℝ)/50000000) :=
    term_le 73 32 3584 (14548216/(10:ℝ)^10) ((16494470767:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 14 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 14548216/(10:ℝ)^10) 48 3584 ≤ ((589510144787:ℝ)/250000000) :=
    term_le 73 48 3584 (14548216/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 14548216/(10:ℝ)^10) 72 1 ≤ ((493160581:ℝ)/500000000) :=
    term_le 73 72 1 (14548216/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b11 : T 73 (1 - 14548216/(10:ℝ)^10) 1 1 ≤ ((13784403:ℝ)/1000000000) :=
    term_le 73 1 1 (14548216/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((10427:ℕ):ℝ) = 10427 := by norm_num
  have e4 : T 73 (1 - 14548216/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 14548216/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_S5Ex
