import Work.BlockAccounting.Moment
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of the exit-frame certificate U8XM with mutable sources (hub tail in phase B), #401 dedicated banks, and an explicit x-exit-words pad-minimal (remapped) bank schedule. -/
def blocksU8XMx : List (ℕ × ℕ) := [(1, 183816), (2, 24336), (3, 32112), (4, 11664), (5, 5184), (8, 6312), (16, 3584), (27, 1), (32, 3584), (48, 3584)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksU8XMx_rank : (blocksU8XMx.map fun b => b.1 * b.2).sum = 10939 * 73 - 2560 := by
  norm_num [blocksU8XMx]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksU8XMx_tally : (blocksU8XMx.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((10939:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksU8XMx]

theorem blocksU8XMx_count : (blocksU8XMx.map fun b => b.2).sum = 274177 := by
  norm_num [blocksU8XMx]

theorem blocksU8XMx_proper : ∀ b ∈ blocksU8XMx, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksU8XMx` at exponent `1 - 13868857/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 10939).  Relative slack of the rational certificate: 3.58e-12. -/
theorem fold_U8XMx :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 13868857/(10:ℝ)^10) blocksU8XMx
      + ((10939:ℕ):ℝ) * (T 73 (1 - 13868857/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 13868857/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((10939:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 13868857/(10:ℝ)^10) blocksU8XMx
      = T 73 (1 - 13868857/(10:ℝ)^10) 1 183816 + (T 73 (1 - 13868857/(10:ℝ)^10) 2 24336 + (T 73 (1 - 13868857/(10:ℝ)^10) 3 32112 + (T 73 (1 - 13868857/(10:ℝ)^10) 4 11664 + (T 73 (1 - 13868857/(10:ℝ)^10) 5 5184 + (T 73 (1 - 13868857/(10:ℝ)^10) 8 6312 + (T 73 (1 - 13868857/(10:ℝ)^10) 16 3584 + (T 73 (1 - 13868857/(10:ℝ)^10) 27 1 + (T 73 (1 - 13868857/(10:ℝ)^10) 32 3584 + (T 73 (1 - 13868857/(10:ℝ)^10) 48 3584 + 0))))))))) := rfl
  have b0 : T 73 (1 - 13868857/(10:ℝ)^10) 1 183816 ≤ ((10132221219803:ℝ)/4000000000) :=
    term_le 73 1 183816 (13868857/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 13868857/(10:ℝ)^10) 2 24336 ≤ ((33503722368247:ℝ)/50000000000) :=
    term_le 73 2 24336 (13868857/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 13868857/(10:ℝ)^10) 3 32112 ≤ ((6627630064513:ℝ)/5000000000) :=
    term_le 73 3 32112 (13868857/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 13868857/(10:ℝ)^10) 4 11664 ≤ ((32085135337637:ℝ)/50000000000) :=
    term_le 73 4 11664 (13868857/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 13868857/(10:ℝ)^10) 5 5184 ≤ ((7127823847757:ℝ)/20000000000) :=
    term_le 73 5 5184 (13868857/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 13868857/(10:ℝ)^10) 8 6312 ≤ ((69385041279969:ℝ)/100000000000) :=
    term_le 73 8 6312 (13868857/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 13868857/(10:ℝ)^10) 16 3584 ≤ ((39359481362957:ℝ)/50000000000) :=
    term_le 73 16 3584 (13868857/(10:ℝ)^10) ((151787071891:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 17 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 13868857/(10:ℝ)^10) 27 1 ≤ ((18518678257:ℝ)/50000000000) :=
    term_le 73 27 1 (13868857/(10:ℝ)^10) ((497311287573:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 15 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 13868857/(10:ℝ)^10) 32 3584 ≤ ((157286650523773:ℝ)/100000000000) :=
    term_le 73 32 3584 (13868857/(10:ℝ)^10) ((16494470767:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 14 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 13868857/(10:ℝ)^10) 48 3584 ≤ ((235797341674031:ℝ)/100000000000) :=
    term_le 73 48 3584 (13868857/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 13868857/(10:ℝ)^10) 72 1 ≤ ((3082250743:ℝ)/3125000000) :=
    term_le 73 72 1 (13868857/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b11 : T 73 (1 - 13868857/(10:ℝ)^10) 1 1 ≤ ((1378038531:ℝ)/100000000000) :=
    term_le 73 1 1 (13868857/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((10939:ℕ):ℝ) = 10939 := by norm_num
  have e4 : T 73 (1 - 13868857/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 13868857/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_U8XMx
