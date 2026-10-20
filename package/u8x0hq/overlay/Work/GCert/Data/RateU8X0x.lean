import Work.BlockAccounting.Moment
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of the exit-frame certificate U8X0 with mutable sources (hub tail in phase B), #401 dedicated banks, and an explicit x-exit-words pad-minimal (remapped) bank schedule. -/
def blocksU8X0x : List (ℕ × ℕ) := [(1, 183456), (2, 24336), (3, 32112), (4, 11664), (5, 5184), (8, 5952), (16, 3584), (32, 3584), (48, 3584), (55, 1)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksU8X0x_rank : (blocksU8X0x.map fun b => b.1 * b.2).sum = 10895 * 73 - 2560 := by
  norm_num [blocksU8X0x]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksU8X0x_tally : (blocksU8X0x.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((10895:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksU8X0x]

theorem blocksU8X0x_count : (blocksU8X0x.map fun b => b.2).sum = 273457 := by
  norm_num [blocksU8X0x]

theorem blocksU8X0x_proper : ∀ b ∈ blocksU8X0x, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksU8X0x` at exponent `1 - 13928619/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 10895).  Relative slack of the rational certificate: 1.34e-11. -/
theorem fold_U8X0x :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 13928619/(10:ℝ)^10) blocksU8X0x
      + ((10895:ℕ):ℝ) * (T 73 (1 - 13928619/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 13928619/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((10895:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 13928619/(10:ℝ)^10) blocksU8X0x
      = T 73 (1 - 13928619/(10:ℝ)^10) 1 183456 + (T 73 (1 - 13928619/(10:ℝ)^10) 2 24336 + (T 73 (1 - 13928619/(10:ℝ)^10) 3 32112 + (T 73 (1 - 13928619/(10:ℝ)^10) 4 11664 + (T 73 (1 - 13928619/(10:ℝ)^10) 5 5184 + (T 73 (1 - 13928619/(10:ℝ)^10) 8 5952 + (T 73 (1 - 13928619/(10:ℝ)^10) 16 3584 + (T 73 (1 - 13928619/(10:ℝ)^10) 32 3584 + (T 73 (1 - 13928619/(10:ℝ)^10) 48 3584 + (T 73 (1 - 13928619/(10:ℝ)^10) 55 1 + 0))))))))) := rfl
  have b0 : T 73 (1 - 13928619/(10:ℝ)^10) 1 183456 ≤ ((3160198986773:ℝ)/1250000000) :=
    term_le 73 1 183456 (13928619/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 13928619/(10:ℝ)^10) 2 24336 ≤ ((3350444265061:ℝ)/5000000000) :=
    term_le 73 2 24336 (13928619/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 13928619/(10:ℝ)^10) 3 32112 ≤ ((1656939122237:ℝ)/1250000000) :=
    term_le 73 3 32112 (13928619/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 13928619/(10:ℝ)^10) 4 11664 ≤ ((641713844189:ℝ)/1000000000) :=
    term_le 73 4 11664 (13928619/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 13928619/(10:ℝ)^10) 5 5184 ≤ ((222748064163:ℝ)/625000000) :=
    term_le 73 5 5184 (13928619/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 13928619/(10:ℝ)^10) 8 5952 ≤ ((6542858406923:ℝ)/10000000000) :=
    term_le 73 8 5952 (13928619/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 13928619/(10:ℝ)^10) 16 3584 ≤ ((7871967679719:ℝ)/10000000000) :=
    term_le 73 16 3584 (13928619/(10:ℝ)^10) ((151787071891:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 17 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 13928619/(10:ℝ)^10) 32 3584 ≤ ((15728742574659:ℝ)/10000000000) :=
    term_le 73 32 3584 (13928619/(10:ℝ)^10) ((16494470767:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 14 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 13928619/(10:ℝ)^10) 48 3584 ≤ ((1473737078013:ℝ)/625000000) :=
    term_le 73 48 3584 (13928619/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 13928619/(10:ℝ)^10) 55 1 ≤ ((1507443667:ℝ)/2000000000) :=
    term_le 73 55 1 (13928619/(10:ℝ)^10) ((283126255917:ℝ)/1000000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 10 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 13928619/(10:ℝ)^10) 72 1 ≤ ((9863203191:ℝ)/10000000000) :=
    term_le 73 72 1 (13928619/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b11 : T 73 (1 - 13928619/(10:ℝ)^10) 1 1 ≤ ((137807387:ℝ)/10000000000) :=
    term_le 73 1 1 (13928619/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((10895:ℕ):ℝ) = 10895 := by norm_num
  have e4 : T 73 (1 - 13928619/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 13928619/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_U8X0x
