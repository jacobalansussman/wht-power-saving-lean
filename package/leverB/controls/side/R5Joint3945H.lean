import Work.Mixed.TwoTarget
import Work.Mixed.Rate3945H

/-!
# R5 x R1/R2 joint check for design B (#3945 with every child-8/10 bin a literal head, suffix 3945H)

bB-data (wht45/allheads, 19 Oct 2026): e-r45m1's `R5Joint.lean` with the #3945H numbers.  For every enumeration `rk`
of the paid rounds of the half-cut word of #3945H (histogram `cut3945H`, ranks in `[1, 48]`):
* `rateQ3945H`: the rule `twoTarget 73 49 322560 80640 P 48 rk` satisfies `RateQ zRate3945H 1` (Theorem 1G's rate fact);
* `allotM3945H`: e-r12's R2 bound at the seam exponent `alpha = 1 - (100 N + 1)/10^12`, N = 17581520:
  `allotM n0 n ≤ C 2^n (n+1)^alpha` for every threshold `n0 ≥ 73`.
Prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.MX
open BlockAccounting

theorem rateQ3945H {P : ℕ} (rk : Fin P → ℕ) (hrk1 : ∀ j, 1 ≤ rk j) (hrk : ∀ j, rk j ≤ 48)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut3945H r) :
    (Rule.twoTarget 73 49 322560 80640 P 48 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num) (by norm_num)).RateQ
      zRate3945H 1 :=
  Rule.twoTarget_rateQ_hist hrk1 hrk _ _ _ _ (rate3945H rk hcnt)

theorem allotM3945H {P : ℕ} (rk : Fin P → ℕ) (hrk1 : ∀ j, 1 ≤ rk j) (hrk : ∀ j, rk j ≤ 48)
    (hcnt : ∀ r, (Finset.univ.filter fun j => rk j = r).card = histFn cut3945H r) {n0 : ℕ} (hn0 : 73 ≤ n0) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ n,
      ((Rule.twoTarget 73 49 322560 80640 P 48 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num)
        (by norm_num)).allotM n0 n : ℝ) ≤ C * 2 ^ n * ((n : ℝ) + 1) ^ (1 - (1758152001:ℝ)/(10:ℝ)^12) :=
  (Rule.twoTarget 73 49 322560 80640 P 48 rk hrk1 hrk (by norm_num) (by norm_num) (by norm_num)
    (by norm_num)).allotM_bound (rateQ3945H rk hrk1 hrk hcnt)
    (by norm_num [zRate3945H]) (by norm_num [zRate3945H]) (by norm_num) hn0

end OAI.PowerSaving.MX

#print axioms OAI.PowerSaving.MX.rateQ3945H
#print axioms OAI.PowerSaving.MX.allotM3945H
