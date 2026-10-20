import Work.BlockAccounting.Moment
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of PR588 atm3 WITHOUT source copies (mutable sources), #401 dedicated banks, and an explicit packed bank schedule. -/
def blocksC588Mx : List (ℕ × ℕ) := [(1, 178848), (2, 35424), (3, 20160), (4, 16704), (5, 4032), (6, 4032), (8, 5952), (16, 3584), (32, 3584), (45, 1), (48, 3584)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksC588Mx_rank : (blocksC588Mx.map fun b => b.1 * b.2).sum = 11173 * 73 - 2560 := by
  norm_num [blocksC588Mx]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksC588Mx_tally : (blocksC588Mx.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((11173:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksC588Mx]

theorem blocksC588Mx_count : (blocksC588Mx.map fun b => b.2).sum = 275905 := by
  norm_num [blocksC588Mx]

theorem blocksC588Mx_proper : ∀ b ∈ blocksC588Mx, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksC588Mx` at exponent `1 - 13566303/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 11173).  Relative slack of the rational certificate: 1.71e-10. -/
theorem fold_C588Mx :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 13566303/(10:ℝ)^10) blocksC588Mx
      + ((11173:ℕ):ℝ) * (T 73 (1 - 13566303/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 13566303/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((11173:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 13566303/(10:ℝ)^10) blocksC588Mx
      = T 73 (1 - 13566303/(10:ℝ)^10) 1 178848 + (T 73 (1 - 13566303/(10:ℝ)^10) 2 35424 + (T 73 (1 - 13566303/(10:ℝ)^10) 3 20160 + (T 73 (1 - 13566303/(10:ℝ)^10) 4 16704 + (T 73 (1 - 13566303/(10:ℝ)^10) 5 4032 + (T 73 (1 - 13566303/(10:ℝ)^10) 6 4032 + (T 73 (1 - 13566303/(10:ℝ)^10) 8 5952 + (T 73 (1 - 13566303/(10:ℝ)^10) 16 3584 + (T 73 (1 - 13566303/(10:ℝ)^10) 32 3584 + (T 73 (1 - 13566303/(10:ℝ)^10) 45 1 + (T 73 (1 - 13566303/(10:ℝ)^10) 48 3584 + 0)))))))))) := rfl
  have b0 : T 73 (1 - 13566303/(10:ℝ)^10) 1 178848 ≤ ((2464274441781:ℝ)/1000000000) :=
    term_le 73 1 178848 (13566303/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 13566303/(10:ℝ)^10) 2 35424 ≤ ((60954280423:ℝ)/62500000) :=
    term_le 73 2 35424 (13566303/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 13566303/(10:ℝ)^10) 3 20160 ≤ ((416044218957:ℝ)/500000000) :=
    term_le 73 3 20160 (13566303/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 13566303/(10:ℝ)^10) 4 16704 ≤ ((229725227397:ℝ)/250000000) :=
    term_le 73 4 16704 (13566303/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 13566303/(10:ℝ)^10) 5 4032 ≤ ((27717066587:ℝ)/100000000) :=
    term_le 73 5 4032 (13566303/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 13566303/(10:ℝ)^10) 6 4032 ≤ ((332522541597:ℝ)/1000000000) :=
    term_le 73 6 4032 (13566303/(10:ℝ)^10) ((1249349985961:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 13566303/(10:ℝ)^10) 8 5952 ≤ ((654233428713:ℝ)/1000000000) :=
    term_le 73 8 5952 (13566303/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 13566303/(10:ℝ)^10) 16 3584 ≤ ((98394184667:ℝ)/125000000) :=
    term_le 73 16 3584 (13566303/(10:ℝ)^10) ((151787071891:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 17 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 13566303/(10:ℝ)^10) 32 3584 ≤ ((1572827259019:ℝ)/1000000000) :=
    term_le 73 32 3584 (13566303/(10:ℝ)^10) ((16494470767:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 14 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 13566303/(10:ℝ)^10) 45 1 ≤ ((616843079:ℝ)/1000000000) :=
    term_le 73 45 1 (13566303/(10:ℝ)^10) ((24189847569:ℝ)/50000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 12 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 13566303/(10:ℝ)^10) 48 3584 ≤ ((2357943506433:ℝ)/1000000000) :=
    term_le 73 48 3584 (13566303/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b11 : T 73 (1 - 13566303/(10:ℝ)^10) 72 1 ≤ ((986319827:ℝ)/1000000000) :=
    term_le 73 72 1 (13566303/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b12 : T 73 (1 - 13566303/(10:ℝ)^10) 1 1 ≤ ((13778597:ℝ)/1000000000) :=
    term_le 73 1 1 (13566303/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((11173:ℕ):ℝ) = 11173 := by norm_num
  have e4 : T 73 (1 - 13566303/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 13566303/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_C588Mx
