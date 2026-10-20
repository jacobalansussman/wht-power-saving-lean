import Work.BlockAccounting.Moment4
set_option maxRecDepth 4000
namespace OAI.PowerSaving.BlockAccounting
open Finset
/-- 72 copies of the U8 word (#869/#870/#890/#942) with hub exit frames, mutable sources (no copies: split centres on sources, fsp-cc), #401 banks for the helpers that are not hosted, same-run idle-axis hosting (DaysSky #442) of the others; fs-price. -/
def blocksU8X0HQx : List (ℕ × ℕ) := [(1, 185776), (2, 25200), (3, 32368), (4, 12560), (5, 5184), (8, 576), (12, 96), (48, 2688), (59, 1)]

/-- total rank of one unit (= number of unit moves it stands for): `W * m - D` -/
theorem blocksU8X0HQx_rank : (blocksU8X0HQx.map fun b => b.1 * b.2).sum = 7491 * 73 - 2560 := by
  norm_num [blocksU8X0HQx]

/-- the same as a price: the unit list at the price list `r ↦ r` (shape of `hPid` in `RAM.engine_program_group_perRank`) -/
theorem blocksU8X0HQx_tally : (blocksU8X0HQx.map fun b => (b.2:ℝ) * ((b.1:ℕ):ℝ)).sum = ((7491:ℕ):ℝ) * ((73:ℕ):ℝ) - ((2560:ℕ):ℝ) := by
  norm_num [blocksU8X0HQx]

theorem blocksU8X0HQx_count : (blocksU8X0HQx.map fun b => b.2).sum = 264449 := by
  norm_num [blocksU8X0HQx]

theorem blocksU8X0HQx_proper : ∀ b ∈ blocksU8X0HQx, 1 ≤ b.1 ∧ b.1 < 73 := by
  decide

/-- **whole-block fact** of `blocksU8X0HQx` at exponent `1 - 16256630/10^10`, fill parameter `s = 40`: the hypothesis `hnum` of `RAM.engine_program_group_list` (u = 73, W = 7491).  Relative slack of the rational certificate: 2.58e-10. -/
theorem fold_U8X0HQx :
    ((2:ℝ)^40 - 1) * moment 73 (1 - 16256630/(10:ℝ)^10) blocksU8X0HQx
      + ((7491:ℕ):ℝ) * (T 73 (1 - 16256630/(10:ℝ)^10) (73 - 1) 1 + T 73 (1 - 16256630/(10:ℝ)^10) 1 1)
      < (2:ℝ)^40 * ((7491:ℕ):ℝ) := by
  have hsum : moment 73 (1 - 16256630/(10:ℝ)^10) blocksU8X0HQx
      = T 73 (1 - 16256630/(10:ℝ)^10) 1 185776 + (T 73 (1 - 16256630/(10:ℝ)^10) 2 25200 + (T 73 (1 - 16256630/(10:ℝ)^10) 3 32368 + (T 73 (1 - 16256630/(10:ℝ)^10) 4 12560 + (T 73 (1 - 16256630/(10:ℝ)^10) 5 5184 + (T 73 (1 - 16256630/(10:ℝ)^10) 8 576 + (T 73 (1 - 16256630/(10:ℝ)^10) 12 96 + (T 73 (1 - 16256630/(10:ℝ)^10) 48 2688 + (T 73 (1 - 16256630/(10:ℝ)^10) 59 1 + 0)))))))) := rfl
  have b0 : T 73 (1 - 16256630/(10:ℝ)^10) 1 185776 ≤ ((2562688869521:ℝ)/1000000000) :=
    term_le4 73 1 185776 (16256630/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b1 : T 73 (1 - 16256630/(10:ℝ)^10) 2 25200 ≤ ((347230161539:ℝ)/500000000) :=
    term_le4 73 2 25200 (16256630/(10:ℝ)^10) ((359731226059:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 24 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b2 : T 73 (1 - 16256630/(10:ℝ)^10) 3 32368 ≤ ((668555954221:ℝ)/500000000) :=
    term_le4 73 3 32368 (16256630/(10:ℝ)^10) ((1595923576241:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 23 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b3 : T 73 (1 - 16256630/(10:ℝ)^10) 4 12560 ≤ ((4321725477:ℝ)/6250000) :=
    term_le4 73 4 12560 (16256630/(10:ℝ)^10) ((290416508003:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 22 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b4 : T 73 (1 - 16256630/(10:ℝ)^10) 5 5184 ≤ ((35661941433:ℝ)/100000000) :=
    term_le4 73 5 5184 (16256630/(10:ℝ)^10) ((670255382179:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 21 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b5 : T 73 (1 - 16256630/(10:ℝ)^10) 8 576 ≤ ((15837646093:ℝ)/250000000) :=
    term_le4 73 8 576 (16256630/(10:ℝ)^10) ((221101789947:ℝ)/100000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 20 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b6 : T 73 (1 - 16256630/(10:ℝ)^10) 12 96 ≤ ((3165442031:ℝ)/200000000) :=
    term_le4 73 12 96 (16256630/(10:ℝ)^10) ((902776395681:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 18 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b7 : T 73 (1 - 16256630/(10:ℝ)^10) 48 2688 ≤ ((1768657112869:ℝ)/1000000000) :=
    term_le4 73 48 2688 (16256630/(10:ℝ)^10) ((209629215121:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 11 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b8 : T 73 (1 - 16256630/(10:ℝ)^10) 59 1 ≤ ((101062373:ℝ)/125000000) :=
    term_le4 73 59 1 (16256630/(10:ℝ)^10) ((53230499311:ℝ)/250000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 10 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b9 : T 73 (1 - 16256630/(10:ℝ)^10) 72 1 ≤ ((986323487:ℝ)/1000000000) :=
    term_le4 73 72 1 (16256630/(10:ℝ)^10) ((6896661067:ℝ)/500000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 6 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have b10 : T 73 (1 - 16256630/(10:ℝ)^10) 1 1 ≤ ((1379451:ℝ)/100000000) :=
    term_le4 73 1 1 (16256630/(10:ℝ)^10) ((85809188823:ℝ)/20000000000) _ (by norm_num) (by norm_num) (by norm_num) (by norm_num)
      (log_le_of_le_expsum _ _ 26 (by norm_num) (by norm_num)
        (by norm_num [Finset.sum_range_succ, Nat.factorial]))
      (by norm_num) (by norm_num)
  have e1 : ((2:ℝ)^40 - 1) = 1099511627775 := by norm_num
  have e2 : (2:ℝ)^40 = 1099511627776 := by norm_num
  have e3 : ((7491:ℕ):ℝ) = 7491 := by norm_num
  have e4 : T 73 (1 - 16256630/(10:ℝ)^10) (73 - 1) 1 = T 73 (1 - 16256630/(10:ℝ)^10) 72 1 := rfl
  rw [e1, e2, e3, e4, hsum]
  linarith

end OAI.PowerSaving.BlockAccounting
#print axioms OAI.PowerSaving.BlockAccounting.fold_U8X0HQx
