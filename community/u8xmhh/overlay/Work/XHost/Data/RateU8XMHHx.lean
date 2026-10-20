import Work.BlockAccounting.Moment
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of the U8 word (#869/#870/#890/#942) with hub exit frames, mutable sources and 5 centre-source copies, #401 banks for the helpers that are not hosted, same-run idle-axis hosting (DaysSky #442) of the others; fs-price. -/
def blocksU8XMHHx : List (ℕ × ℕ) := [(1, 186137), (2, 25216), (3, 32368), (4, 12560), (5, 5184), (8, 936), (12, 64), (48, 2608), (61, 1)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksU8XMHHx_rank : (blocksU8XMHHx.map fun b => b.1 * b.2).sum = 7478 * 73 - 2560 := by
  norm_num [blocksU8XMHHx]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksU8XMHHx_tally : (blocksU8XMHHx.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((7478:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksU8XMHHx]

theorem blocksU8XMHHx_count : (blocksU8XMHHx.map fun b => b.2).sum = 265074 := by
  norm_num [blocksU8XMHHx]

theorem blocksU8XMHHx_proper : ∀ b ∈ blocksU8XMHHx, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksU8XMHHx` at exponent `1 - 16197706/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 7478).  Relative slack of the rational certificate: 5.96e-11. -/
theorem fold_U8XMHHx :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 16197706/(10:ℝ)^10) blocksU8XMHHx
      + ((7478:ℕ):ℝ) * (T 73 (1 - 16197706/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 16197706/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((7478:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 16197706/(10:ℝ)^10) blocksU8XMHHx
      = T 73 (1 - 16197706/(10:ℝ)^10) 1 186137 + (T 73 (1 - 16197706/(10:ℝ)^10) 2 25216 + (T 73 (1 - 16197706/(10:ℝ)^10) 3 32368 + (T 73 (1 - 16197706/(10:ℝ)^10) 4 12560 + (T 73 (1 - 16197706/(10:ℝ)^10) 5 5184 + (T 73 (1 - 16197706/(10:ℝ)^10) 8 936 + (T 73 (1 - 16197706/(10:ℝ)^10) 12 64 + (T 73 (1 - 16197706/(10:ℝ)^10) 48 2608 + (T 73 (1 - 16197706/(10:ℝ)^10) 61 1 + 0)))))))) := rfl
  have b0 : T 73 (1 - 16197706/(10:ℝ)^10) 1 186137 ≤ ((25676038221493:ℝ)/10000000000) :=
    term_le 73 1 186137 (16197706/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 16197706/(10:ℝ)^10) 2 25216 ≤ ((434304080181:ℝ)/625000000) :=
    term_le 73 2 25216 (16197706/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 16197706/(10:ℝ)^10) 3 32368 ≤ ((3342716927303:ℝ)/2500000000) :=
    term_le 73 3 32368 (16197706/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 16197706/(10:ℝ)^10) 4 12560 ≤ ((6914642474967:ℝ)/10000000000) :=
    term_le 73 4 12560 (16197706/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 16197706/(10:ℝ)^10) 5 5184 ≤ ((1783068911169:ℝ)/5000000000) :=
    term_le 73 5 5184 (16197706/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 16197706/(10:ℝ)^10) 8 936 ≤ ((205886717379:ℝ)/2000000000) :=
    term_le 73 8 936 (16197706/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 16197706/(10:ℝ)^10) 12 64 ≤ ((105513611939:ℝ)/10000000000) :=
    term_le 73 12 64 (16197706/(10:ℝ)^10) ((902776395681:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 18 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 16197706/(10:ℝ)^10) 48 2608 ≤ ((17160142690363:ℝ)/10000000000) :=
    term_le 73 48 2608 (16197706/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 16197706/(10:ℝ)^10) 61 1 ≤ ((8358595441:ℝ)/10000000000) :=
    term_le 73 61 1 (16197706/(10:ℝ)^10) ((179585576977:ℝ)/1000000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 9 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 16197706/(10:ℝ)^10) 72 1 ≤ ((9863234061:ℝ)/10000000000) :=
    term_le 73 72 1 (16197706/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 16197706/(10:ℝ)^10) 1 1 ≤ ((27588323:ℝ)/2000000000) :=
    term_le 73 1 1 (16197706/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((7478:ℕ):ℝ) = 7478 := by norm_num
  have e4 : T 73 (1 - 16197706/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 16197706/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_U8XMHHx
