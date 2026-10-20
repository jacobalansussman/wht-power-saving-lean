import Work.BlockAccounting.Moment
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of the exit-frame certificate U8X (hub tail in phase B), #401 dedicated banks, and an explicit x-exit-words pad-minimal (remapped) bank schedule. -/
def blocksU8Xx : List (ℕ × ℕ) := [(1, 191520), (2, 24336), (3, 32112), (4, 11664), (5, 5184), (8, 14016), (16, 3584), (32, 3584), (41, 1), (48, 3584)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksU8Xx_rank : (blocksU8Xx.map fun b => b.1 * b.2).sum = 11889 * 73 - 2560 := by
  norm_num [blocksU8Xx]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksU8Xx_tally : (blocksU8Xx.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((11889:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksU8Xx]

theorem blocksU8Xx_count : (blocksU8Xx.map fun b => b.2).sum = 289585 := by
  norm_num [blocksU8Xx]

theorem blocksU8Xx_proper : ∀ b ∈ blocksU8Xx, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksU8Xx` at exponent `1 - 12704024/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 11889).  Relative slack of the rational certificate: 6.49e-12. -/
theorem fold_U8Xx :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 12704024/(10:ℝ)^10) blocksU8Xx
      + ((11889:ℕ):ℝ) * (T 73 (1 - 12704024/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 12704024/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((11889:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 12704024/(10:ℝ)^10) blocksU8Xx
      = T 73 (1 - 12704024/(10:ℝ)^10) 1 191520 + (T 73 (1 - 12704024/(10:ℝ)^10) 2 24336 + (T 73 (1 - 12704024/(10:ℝ)^10) 3 32112 + (T 73 (1 - 12704024/(10:ℝ)^10) 4 11664 + (T 73 (1 - 12704024/(10:ℝ)^10) 5 5184 + (T 73 (1 - 12704024/(10:ℝ)^10) 8 14016 + (T 73 (1 - 12704024/(10:ℝ)^10) 16 3584 + (T 73 (1 - 12704024/(10:ℝ)^10) 32 3584 + (T 73 (1 - 12704024/(10:ℝ)^10) 41 1 + (T 73 (1 - 12704024/(10:ℝ)^10) 48 3584 + 0))))))))) := rfl
  have b0 : T 73 (1 - 12704024/(10:ℝ)^10) 1 191520 ≤ ((13189503606963:ℝ)/5000000000) :=
    term_le 73 1 191520 (12704024/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 12704024/(10:ℝ)^10) 2 24336 ≤ ((418621078733:ℝ)/625000000) :=
    term_le 73 2 24336 (12704024/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 12704024/(10:ℝ)^10) 3 32112 ≤ ((2650066553191:ℝ)/2000000000) :=
    term_le 73 3 32112 (12704024/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 12704024/(10:ℝ)^10) 4 11664 ≤ ((1282971326841:ℝ)/2000000000) :=
    term_le 73 4 11664 (12704024/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 12704024/(10:ℝ)^10) 5 5184 ≤ ((445349888271:ℝ)/1250000000) :=
    term_le 73 5 5184 (12704024/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 12704024/(10:ℝ)^10) 8 14016 ≤ ((15403205103837:ℝ)/10000000000) :=
    term_le 73 8 14016 (12704024/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 12704024/(10:ℝ)^10) 16 3584 ≤ ((3935252295773:ℝ)/5000000000) :=
    term_le 73 16 3584 (12704024/(10:ℝ)^10) ((151787071891:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 17 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 12704024/(10:ℝ)^10) 32 3584 ≤ ((7863577063273:ℝ)/5000000000) :=
    term_le 73 32 3584 (12704024/(10:ℝ)^10) ((16494470767:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 14 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 12704024/(10:ℝ)^10) 41 1 ≤ ((1405139009:ℝ)/2500000000) :=
    term_le 73 41 1 (12704024/(10:ℝ)^10) ((288443687223:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 12 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 12704024/(10:ℝ)^10) 48 3584 ≤ ((23578582641287:ℝ)/10000000000) :=
    term_le 73 48 3584 (12704024/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 12704024/(10:ℝ)^10) 72 1 ≤ ((9863186531:ℝ)/10000000000) :=
    term_le 73 72 1 (12704024/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b11 : T 73 (1 - 12704024/(10:ℝ)^10) 1 1 ≤ ((137735001:ℝ)/10000000000) :=
    term_le 73 1 1 (12704024/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((11889:ℕ):ℝ) = 11889 := by norm_num
  have e4 : T 73 (1 - 12704024/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 12704024/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_U8Xx
