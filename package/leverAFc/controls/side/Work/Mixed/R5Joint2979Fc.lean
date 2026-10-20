import Work.Mixed.TwoTarget
import Work.Mixed.Rate2979Fc

/-!
# R5 x R1/R2 for the FULL cut of A-fc (suffix 2979Fc): the mixed rule, its rate fact and the R2 bound

afc-cut (wht45/afc, 20 Oct 2026): fc-rate's `R5Joint3945F.lean` with A-fc's numbers.  A-fc is the hub's word of
PR #2979 (T74 clock) with the re-solved fill u1-2979sep (round wht45/fullcut).  Two targets: `W - H = 49728` arrays
keep the full target `m = 73`; the `H = 16576` head arrays get the cut target `kc = 25 = 73 - 48` (each literal
rank-48 head cut by all 48 axes).  There is no other rank-48 block, so nothing is split; every paid rank lies in
`[1, 23]`, so the rule is separated: `rmax = 23 < 25 = kc`.

For every enumeration `rk` of the paid rounds whose histogram is `cut2979Fc`:
* `rk_bounds2979Fc`: the ranks lie in `[1, 23]` (from the histogram alone);
* `rateQ2979Fc`: the rule `twoTarget 73 25 66304 16576 P 23 rk` satisfies `RateQ zRate2979Fc 1` (Theorem 1G's rate
  fact) at `zRate2979Fc = 1 - 1788704502/10^12`;
* `allotM2979Fc`: e-r12's R2 bound at the seam exponent `alpha = 1 - (100 N + 1)/10^12`, N = 17887045:
  `allotM n0 n ≤ C 2^n (n+1)^alpha` for every threshold `n0 ≥ 73`;
* `rateQ2979Fc'` / `allotM2979Fc'`: the same with the rank bounds discharged by `rk_bounds2979Fc`.
Prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.MX
open BlockAccounting

/-- Separation of the full cut: every paid rank is at most 23, below the cut target 25. -/
theorem sep2979Fc : (23 : ℕ) < 25 := by norm_num

/-- Every rank of an enumeration with the histogram of `cut2979Fc` lies in `[1, 23]`. -/
theorem rk_bounds2979Fc {P : ℕ} (rk : Fin P → ℕ)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut2979Fc r) (j : Fin P) :
    1 ≤ rk j ∧ rk j ≤ 23 := by
  have hpos : 0 < (Finset.univ.filter fun j' => rk j' = rk j).card :=
    Finset.card_pos.mpr ⟨j, by simp⟩
  rw [hcnt] at hpos
  have hk := mem_keys_of_histFn_ne_zero cut2979Fc (rk j) (Nat.pos_iff_ne_zero.mp hpos)
  unfold keys at hk
  rw [List.mem_toFinset, List.mem_map] at hk
  obtain ⟨b, hb, hbr⟩ := hk
  have := cut2979Fc_proper b hb
  rw [← hbr]
  exact this

theorem rateQ2979Fc {P : ℕ} (rk : Fin P → ℕ) (hrk1 : ∀ j, 1 ≤ rk j) (hrk : ∀ j, rk j ≤ 23)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut2979Fc r) :
    (Rule.twoTarget 73 25 66304 16576 P 23 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num) (by norm_num)).RateQ
      zRate2979Fc 1 :=
  Rule.twoTarget_rateQ_hist hrk1 hrk _ _ _ _ (rate2979Fc rk hcnt)

theorem allotM2979Fc {P : ℕ} (rk : Fin P → ℕ) (hrk1 : ∀ j, 1 ≤ rk j) (hrk : ∀ j, rk j ≤ 23)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut2979Fc r) {n0 : ℕ} (hn0 : 73 ≤ n0) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ n,
      ((Rule.twoTarget 73 25 66304 16576 P 23 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num)
        (by norm_num)).allotM n0 n : ℝ) ≤ C * 2 ^ n * ((n : ℝ) + 1) ^ (1 - (1788704501:ℝ)/(10:ℝ)^12) :=
  (Rule.twoTarget 73 25 66304 16576 P 23 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num)
    (by norm_num)).allotM_bound (rateQ2979Fc rk hrk1 hrk hcnt)
    (by norm_num [zRate2979Fc]) (by norm_num [zRate2979Fc]) (by norm_num) hn0

/-- `rateQ2979Fc` with the rank bounds taken from the histogram. -/
theorem rateQ2979Fc' {P : ℕ} (rk : Fin P → ℕ)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut2979Fc r) :
    (Rule.twoTarget 73 25 66304 16576 P 23 rk (fun j => (rk_bounds2979Fc rk hcnt j).1)
      (fun j => (rk_bounds2979Fc rk hcnt j).2) (by norm_num) (by norm_num) (by norm_num) (by norm_num)).RateQ
      zRate2979Fc 1 :=
  rateQ2979Fc rk _ _ hcnt

/-- `allotM2979Fc` with the rank bounds taken from the histogram. -/
theorem allotM2979Fc' {P : ℕ} (rk : Fin P → ℕ)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut2979Fc r) {n0 : ℕ} (hn0 : 73 ≤ n0) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ n,
      ((Rule.twoTarget 73 25 66304 16576 P 23 rk (fun j => (rk_bounds2979Fc rk hcnt j).1)
        (fun j => (rk_bounds2979Fc rk hcnt j).2) (by norm_num) (by norm_num) (by norm_num)
        (by norm_num)).allotM n0 n : ℝ) ≤ C * 2 ^ n * ((n : ℝ) + 1) ^ (1 - (1788704501:ℝ)/(10:ℝ)^12) :=
  allotM2979Fc rk _ _ hcnt hn0

end OAI.PowerSaving.MX

#print axioms OAI.PowerSaving.MX.sep2979Fc
#print axioms OAI.PowerSaving.MX.rk_bounds2979Fc
#print axioms OAI.PowerSaving.MX.rateQ2979Fc
#print axioms OAI.PowerSaving.MX.allotM2979Fc
#print axioms OAI.PowerSaving.MX.rateQ2979Fc'
#print axioms OAI.PowerSaving.MX.allotM2979Fc'
