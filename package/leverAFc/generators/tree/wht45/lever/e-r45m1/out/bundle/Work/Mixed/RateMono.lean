import Work.BlockAccounting.Moment

/-!
# The two-target rate fact fails below a failing exponent (e-r45m1, round wht45/lever, 19 Oct 2026)

`rate*_sharp` (generated, Work/Mixed/Rate*.lean) shows that the mixed rate fact `moment ≤ (W - H) + H (kc/m)^z` is
FALSE one grid step above the figure.  This file shows that it is then false at every smaller exponent too (every
larger saving), so no exponent of the chain can reach that step:

* `moment_eq_mul`          `moment m z L = (kc/m)^z * ∑ n (r/kc)^z`;
* `twoTarget_fail_mono`    `Wf + T m z' kc H < moment m z' L` and `z ≤ z'` give the same at `z`
                           (all ranks of `L` in `(0, kc]`, `0 < kc ≤ m`).
No `sorry`, no axioms.
-/

namespace OAI.PowerSaving.MX
open BlockAccounting

/-- the normalised children sum `∑ n (r/kc)^z` -/
noncomputable def kcSum (kc : ℕ) (z : ℝ) (L : List (ℕ × ℕ)) : ℝ := (L.map fun b => (b.2:ℝ) * ((b.1:ℝ)/kc)^z).sum

lemma moment_eq_mul (m kc : ℕ) (hm : 0 < m) (hkc : 0 < kc) (z : ℝ) (L : List (ℕ × ℕ)) :
    moment m z L = ((kc:ℝ)/m)^z * kcSum kc z L := by
  have hm' : (0:ℝ) < m := by exact_mod_cast hm
  have hk' : (0:ℝ) < kc := by exact_mod_cast hkc
  induction L with
  | nil => simp [moment, kcSum]
  | cons b L ih =>
    have hb : ((b.1:ℝ)/m)^z = ((kc:ℝ)/m)^z * ((b.1:ℝ)/kc)^z := by
      rw [← Real.mul_rpow (by positivity) (by positivity)]
      congr 1
      field_simp
    have e1 : moment m z (b :: L) = (b.2:ℝ) * ((b.1:ℝ)/m)^z + moment m z L := by
      simp [moment]
    have e2 : kcSum kc z (b :: L) = (b.2:ℝ) * ((b.1:ℝ)/kc)^z + kcSum kc z L := by
      simp [kcSum]
    rw [e1, e2, ih, hb]
    ring

lemma kcSum_anti (kc : ℕ) (hkc : 0 < kc) (L : List (ℕ × ℕ)) (hL : ∀ b ∈ L, 0 < b.1 ∧ b.1 ≤ kc)
    {z z' : ℝ} (hzz : z ≤ z') : kcSum kc z' L ≤ kcSum kc z L := by
  have hk' : (0:ℝ) < kc := by exact_mod_cast hkc
  induction L with
  | nil => simp [kcSum]
  | cons b L ih =>
    have hb := hL b (by simp)
    have hL' : ∀ c ∈ L, 0 < c.1 ∧ c.1 ≤ kc := fun c hc => hL c (by simp [hc])
    have e : ∀ w, kcSum kc w (b :: L) = (b.2:ℝ) * ((b.1:ℝ)/kc)^w + kcSum kc w L := by
      intro w; simp [kcSum]
    rw [e, e]
    have h0 : (0:ℝ) < (b.1:ℝ)/kc := div_pos (by exact_mod_cast hb.1) hk'
    have h1 : (b.1:ℝ)/kc ≤ 1 := (div_le_one hk').mpr (by exact_mod_cast hb.2)
    have ht : ((b.1:ℝ)/kc)^z' ≤ ((b.1:ℝ)/kc)^z := Real.rpow_le_rpow_of_exponent_ge h0 h1 hzz
    have := mul_le_mul_of_nonneg_left ht (Nat.cast_nonneg b.2 : (0:ℝ) ≤ b.2)
    linarith [ih hL']

/-- the comparison, normalised by `(kc/m)^z` -/
lemma fail_iff (m kc Wf H : ℕ) (hm : 0 < m) (hkc : 0 < kc) (L : List (ℕ × ℕ)) (z : ℝ) :
    (Wf:ℝ) + T m z kc H < moment m z L ↔ (Wf:ℝ) * ((m:ℝ)/kc)^z + H < kcSum kc z L := by
  have hm' : (0:ℝ) < m := by exact_mod_cast hm
  have hk' : (0:ℝ) < kc := by exact_mod_cast hkc
  have hc : 0 < ((kc:ℝ)/m)^z := Real.rpow_pos_of_pos (div_pos hk' hm') z
  have hinv : ((kc:ℝ)/m)^z * ((m:ℝ)/kc)^z = 1 := by
    rw [← Real.mul_rpow (by positivity) (by positivity)]
    rw [show (kc:ℝ)/m * ((m:ℝ)/kc) = 1 by field_simp]
    exact Real.one_rpow z
  have eW : (Wf:ℝ) + T m z kc H = ((kc:ℝ)/m)^z * ((Wf:ℝ) * ((m:ℝ)/kc)^z + H) := by
    unfold T
    have : ((kc:ℝ)/m)^z * ((Wf:ℝ) * ((m:ℝ)/kc)^z) = Wf := by
      rw [mul_comm (Wf:ℝ), ← mul_assoc, hinv, one_mul]
    rw [mul_add, this]
    ring
  rw [eW, moment_eq_mul m kc hm hkc z L]
  exact ⟨fun h => lt_of_mul_lt_mul_left h hc.le, fun h => mul_lt_mul_of_pos_left h hc⟩

/-- **The two-target rate fact fails at every exponent below a failing one.** -/
theorem twoTarget_fail_mono (m kc Wf H : ℕ) (L : List (ℕ × ℕ)) (hm : 0 < m) (hkc : 0 < kc) (hkm : kc ≤ m)
    (hL : ∀ b ∈ L, 0 < b.1 ∧ b.1 ≤ kc) {z z' : ℝ} (hzz : z ≤ z')
    (h : (Wf:ℝ) + T m z' kc H < moment m z' L) : (Wf:ℝ) + T m z kc H < moment m z L := by
  rw [fail_iff m kc Wf H hm hkc] at h ⊢
  have hk' : (0:ℝ) < kc := by exact_mod_cast hkc
  have h1 : (1:ℝ) ≤ (m:ℝ)/kc := (one_le_div hk').mpr (by exact_mod_cast hkm)
  have hp : ((m:ℝ)/kc)^z ≤ ((m:ℝ)/kc)^z' := Real.rpow_le_rpow_of_exponent_le h1 hzz
  have hW := mul_le_mul_of_nonneg_left hp (Nat.cast_nonneg Wf : (0:ℝ) ≤ Wf)
  have hs := kcSum_anti kc hkc L hL hzz
  linarith

end OAI.PowerSaving.MX

#print axioms OAI.PowerSaving.MX.twoTarget_fail_mono
