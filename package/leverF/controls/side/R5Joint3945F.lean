import Work.Mixed.TwoTarget
import Work.Mixed.Rate3945F

/-!
# R5 x R1/R2 for the FULL cut of design B (suffix 3945F): the mixed rule, its rate fact and the R2 bound

fc-rate (wht45/fullcut, 19 Oct 2026): e-r45m1's `R5Joint.lean` / bB-data's `R5Joint3945H.lean` with the full-cut
numbers.  Two targets: `W - H = 241920` arrays keep the full target `m = 73`; the `H = 80640` head arrays get the cut
target `kc = 25 = 73 - 48` (each literal rank-48 head cut by all 48 axes).  Every paid rank lies in `[1, 24]`
(the 57,312 other idle rank-48 blocks split 24 + 24), so the rule is separated: `rmax = 24 < 25 = kc`.

For every enumeration `rk` of the paid rounds whose histogram is `cut3945F`:
* `rk_bounds3945F`: the ranks lie in `[1, 24]` (from the histogram alone);
* `rateQ3945F`: the rule `twoTarget 73 25 322560 80640 P 24 rk` satisfies `RateQ zRate3945F 1` (Theorem 1G's rate
  fact) at `zRate3945F = 1 - 1780693502/10^12`;
* `allotM3945F`: e-r12's R2 bound at the seam exponent `alpha = 1 - (100 N + 1)/10^12`, N = 17806935:
  `allotM n0 n ≤ C 2^n (n+1)^alpha` for every threshold `n0 ≥ 73`;
* `rateQ3945F'` / `allotM3945F'`: the same with the rank bounds discharged by `rk_bounds3945F`.
Prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.MX
open BlockAccounting

/-- Separation of the full cut: every paid rank is at most 24, below the cut target 25. -/
theorem sep3945F : (24 : ℕ) < 25 := by norm_num

/-- Every rank of an enumeration with the histogram of `cut3945F` lies in `[1, 24]`. -/
theorem rk_bounds3945F {P : ℕ} (rk : Fin P → ℕ)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut3945F r) (j : Fin P) :
    1 ≤ rk j ∧ rk j ≤ 24 := by
  have hpos : 0 < (Finset.univ.filter fun j' => rk j' = rk j).card :=
    Finset.card_pos.mpr ⟨j, by simp⟩
  rw [hcnt] at hpos
  have hk := mem_keys_of_histFn_ne_zero cut3945F (rk j) (Nat.pos_iff_ne_zero.mp hpos)
  unfold keys at hk
  rw [List.mem_toFinset, List.mem_map] at hk
  obtain ⟨b, hb, hbr⟩ := hk
  have := cut3945F_proper b hb
  rw [← hbr]
  exact this

theorem rateQ3945F {P : ℕ} (rk : Fin P → ℕ) (hrk1 : ∀ j, 1 ≤ rk j) (hrk : ∀ j, rk j ≤ 24)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut3945F r) :
    (Rule.twoTarget 73 25 322560 80640 P 24 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num) (by norm_num)).RateQ
      zRate3945F 1 :=
  Rule.twoTarget_rateQ_hist hrk1 hrk _ _ _ _ (rate3945F rk hcnt)

theorem allotM3945F {P : ℕ} (rk : Fin P → ℕ) (hrk1 : ∀ j, 1 ≤ rk j) (hrk : ∀ j, rk j ≤ 24)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut3945F r) {n0 : ℕ} (hn0 : 73 ≤ n0) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ n,
      ((Rule.twoTarget 73 25 322560 80640 P 24 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num)
        (by norm_num)).allotM n0 n : ℝ) ≤ C * 2 ^ n * ((n : ℝ) + 1) ^ (1 - (1780693501:ℝ)/(10:ℝ)^12) :=
  (Rule.twoTarget 73 25 322560 80640 P 24 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num)
    (by norm_num)).allotM_bound (rateQ3945F rk hrk1 hrk hcnt)
    (by norm_num [zRate3945F]) (by norm_num [zRate3945F]) (by norm_num) hn0

/-- `rateQ3945F` with the rank bounds taken from the histogram. -/
theorem rateQ3945F' {P : ℕ} (rk : Fin P → ℕ)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut3945F r) :
    (Rule.twoTarget 73 25 322560 80640 P 24 rk (fun j => (rk_bounds3945F rk hcnt j).1)
      (fun j => (rk_bounds3945F rk hcnt j).2) (by norm_num) (by norm_num) (by norm_num) (by norm_num)).RateQ
      zRate3945F 1 :=
  rateQ3945F rk _ _ hcnt

/-- `allotM3945F` with the rank bounds taken from the histogram. -/
theorem allotM3945F' {P : ℕ} (rk : Fin P → ℕ)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut3945F r) {n0 : ℕ} (hn0 : 73 ≤ n0) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ n,
      ((Rule.twoTarget 73 25 322560 80640 P 24 rk (fun j => (rk_bounds3945F rk hcnt j).1)
        (fun j => (rk_bounds3945F rk hcnt j).2) (by norm_num) (by norm_num) (by norm_num)
        (by norm_num)).allotM n0 n : ℝ) ≤ C * 2 ^ n * ((n : ℝ) + 1) ^ (1 - (1780693501:ℝ)/(10:ℝ)^12) :=
  allotM3945F rk _ _ hcnt hn0

end OAI.PowerSaving.MX

#print axioms OAI.PowerSaving.MX.sep3945F
#print axioms OAI.PowerSaving.MX.rk_bounds3945F
#print axioms OAI.PowerSaving.MX.rateQ3945F
#print axioms OAI.PowerSaving.MX.allotM3945F
#print axioms OAI.PowerSaving.MX.rateQ3945F'
#print axioms OAI.PowerSaving.MX.allotM3945F'
