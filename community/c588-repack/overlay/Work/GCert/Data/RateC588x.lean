import Work.BlockAccounting.Moment4
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of PR588 atm3 with source copies, #401 dedicated banks, and an explicit 5,009-bank schedule (pad-optimal; x-pack). -/
def blocksC588x : List (ℕ × ℕ) := [(1, 187636), (2, 35424), (3, 20160), (4, 16704), (5, 4032), (6, 4032), (8, 14016), (16, 3584), (32, 3584), (37, 1), (48, 3584)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksC588x_rank : (blocksC588x.map fun b => b.1 * b.2).sum = 12177 * 73 - 2560 := by
  norm_num [blocksC588x]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksC588x_tally : (blocksC588x.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((12177:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksC588x]

theorem blocksC588x_count : (blocksC588x.map fun b => b.2).sum = 292757 := by
  norm_num [blocksC588x]

theorem blocksC588x_proper : ∀ b ∈ blocksC588x, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksC588x` at exponent `1 - 12383292/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 12177).  Relative slack of the rational certificate: 1.03e-10. -/
theorem fold_C588x :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 12383292/(10:ℝ)^10) blocksC588x
      + ((12177:ℕ):ℝ) * (T 73 (1 - 12383292/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 12383292/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((12177:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 12383292/(10:ℝ)^10) blocksC588x
      = T 73 (1 - 12383292/(10:ℝ)^10) 1 187636 + (T 73 (1 - 12383292/(10:ℝ)^10) 2 35424 + (T 73 (1 - 12383292/(10:ℝ)^10) 3 20160 + (T 73 (1 - 12383292/(10:ℝ)^10) 4 16704 + (T 73 (1 - 12383292/(10:ℝ)^10) 5 4032 + (T 73 (1 - 12383292/(10:ℝ)^10) 6 4032 + (T 73 (1 - 12383292/(10:ℝ)^10) 8 14016 + (T 73 (1 - 12383292/(10:ℝ)^10) 16 3584 + (T 73 (1 - 12383292/(10:ℝ)^10) 32 3584 + (T 73 (1 - 12383292/(10:ℝ)^10) 37 1 + (T 73 (1 - 12383292/(10:ℝ)^10) 48 3584 + 0)))))))))) := rfl
  have b0 : T 73 (1 - 12383292/(10:ℝ)^10) 1 187636 ≤ ((2584048812127:ℝ)/1000000000) :=
    term_le4 73 1 187636 (12383292/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 12383292/(10:ℝ)^10) 2 35424 ≤ ((38994141109:ℝ)/40000000) :=
    term_le4 73 2 35424 (12383292/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 12383292/(10:ℝ)^10) 3 20160 ≤ ((831774297719:ℝ)/1000000000) :=
    term_le4 73 3 20160 (12383292/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 12383292/(10:ℝ)^10) 4 16704 ≤ ((918585257677:ℝ)/1000000000) :=
    term_le4 73 4 16704 (12383292/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 12383292/(10:ℝ)^10) 5 4032 ≤ ((55416553893:ℝ)/200000000) :=
    term_le4 73 5 4032 (12383292/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 12383292/(10:ℝ)^10) 6 4032 ≤ ((332424262093:ℝ)/1000000000) :=
    term_le4 73 6 4032 (12383292/(10:ℝ)^10) ((1249349985961:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 12383292/(10:ℝ)^10) 8 14016 ≤ ((385052820359:ℝ)/250000000) :=
    term_le4 73 8 14016 (12383292/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 12383292/(10:ℝ)^10) 16 3584 ≤ ((787012143821:ℝ)/1000000000) :=
    term_le4 73 16 3584 (12383292/(10:ℝ)^10) ((151787071891:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 17 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 12383292/(10:ℝ)^10) 32 3584 ≤ ((1572673812387:ℝ)/1000000000) :=
    term_le4 73 32 3584 (12383292/(10:ℝ)^10) ((16494470767:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 14 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 12383292/(10:ℝ)^10) 37 1 ≤ ((507276007:ℝ)/1000000000) :=
    term_le4 73 37 1 (12383292/(10:ℝ)^10) ((339770764253:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 13 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 12383292/(10:ℝ)^10) 48 3584 ≤ ((1178913279149:ℝ)/500000000) :=
    term_le4 73 48 3584 (12383292/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b11 : T 73 (1 - 12383292/(10:ℝ)^10) 72 1 ≤ ((986318217:ℝ)/1000000000) :=
    term_le4 73 72 1 (12383292/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b12 : T 73 (1 - 12383292/(10:ℝ)^10) 1 1 ≤ ((2754321:ℝ)/200000000) :=
    term_le4 73 1 1 (12383292/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((12177:ℕ):ℝ) = 12177 := by norm_num
  have e4 : T 73 (1 - 12383292/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 12383292/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_C588x
