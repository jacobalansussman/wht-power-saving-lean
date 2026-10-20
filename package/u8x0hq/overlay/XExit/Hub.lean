import XExit.Halves

/-!
# (key: x-exit) The hub's tail-gate shape gives the end-gate factorisations of `GLabX`

The hub's hyperplane-ended helpers (#591/#593/#599): each exit helper `h` has ONE cleanup gate `h -= c * p`, moved
right after part B; the rest `R` of the cleanup never reads or writes `h`.  In matrix form: a tail matrix `K` whose
nonzero entries sit in exit-slot rows and in columns that are x roles or ordinary slots (`KOk`), with
`R := Mti * (1 + K)` acting as the identity on every exit slot (`hRrow`, `hRcol`).  Then, with `Rinv := (1 - K) * Mt`:
* forward  `(PX + PY + PS * R) * (1 - K) = g8 Mti`                 (`hub_gsplit`, tail `Gt = 1 - K`),
* backward `(PX + PY + Rinvᵀ * PS) * (1 + Kᵀ) = (h8 Mt)ᵀ`          (`hub_hsplit`, tail `Ht = 1 + Kᵀ`: "p += c h"),
* legality: the tails join a helper and its partner (`hub_Gt_legal`, `hub_Ht_legal`, given equal tail labels); the
  rests touch neither y roles nor exit slots (`hub_Gr_legal`, `hub_Hr_legal`).
So for a hub word the instance obligations are SPARSE: the support of `K`, the exit rows/columns of `Mti * (1 + K)`,
and equal tail labels of each helper and its partner.  Works for partners that are x roles or ordinary helpers.
Credits: the tail-gate shape and its commuting move are the hub's (#591 Kothe, #593/#599 Kothe / DreamingOfClouds);
the dirt gates `g8`, `h8` are Jacob Sussman's.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace XEX
open Binary Matrix Finset RAM SS CB BR CR GF RF GX B401
noncomputable section

section Ring
variable {A : Type*} [Ring A]

theorem ring_g (PX PY M K P : A) (hX : PX * K = 0) (hY : PY * K = 0) (hKK : K * K = 0) :
    (PX + PY + P * (M * (1 + K))) * (1 - K) = PX + PY + P * M := by
  have key : (PX + PY + P * (M * (1 + K))) * (1 - K)
      = PX + PY + P * M - PX * K - PY * K - P * M * (K * K) := by noncomm_ring
  rw [key, hX, hY, hKK]
  simp

theorem ring_h (PX PY Mt B P : A) (hA : Mt * PX = PX) (hPYB : PY * B = 0) (hBP : B * P = B)
    (hPB : P * B + PX * B = B) (hBB : B * B = 0) :
    (PX + PY + Mt * (1 - B) * P) * (1 + B) = PX + PY + Mt * P := by
  have key : (PX + PY + Mt * (1 - B) * P) * (1 + B)
      = PX + PY + Mt * P + (PX * B - Mt * PX * B) + PY * B + Mt * (P * B + PX * B - B) - Mt * (B * P - B)
        - Mt * (B * P * B) := by noncomm_ring
  rw [key, hA, hPYB, hPB, hBP, hBB]
  simp

end Ring

section Hub
variable {T Sl : Type} [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl]
variable (E : Finset Sl) (K : Matrix (L3 T Sl) (L3 T Sl) ℚ)

/-- the support of the tail matrix: exit-slot rows, columns that are x roles or ordinary slots. -/
def KOk : Prop := ∀ i j, K i j ≠ 0 → iE E i = 1 ∧ restOk E j

lemma iE_one {i : L3 T Sl} (h : iE E i = 1) : iX i = 0 ∧ iY i = 0 := by
  rcases i with ((t|t)|q)
  · simp [iE] at h
  · simp [iE] at h
  · exact ⟨by simp [iX], by simp [iY]⟩

variable {E K}

lemma PX_K (hK : KOk E K) : PX * K = 0 := by
  ext i j
  simp only [PX, Matrix.diagonal_mul, Matrix.zero_apply]
  by_cases h : K i j = 0
  · simp [h]
  · rw [(iE_one E (hK i j h).1).1, zero_mul]

lemma PY_K (hK : KOk E K) : PY * K = 0 := by
  ext i j
  simp only [PY, Matrix.diagonal_mul, Matrix.zero_apply]
  by_cases h : K i j = 0
  · simp [h]
  · rw [(iE_one E (hK i j h).1).2, zero_mul]

lemma K_PY (hK : KOk E K) : K * PY = 0 := by
  ext i j
  simp only [PY, Matrix.mul_diagonal, Matrix.zero_apply]
  by_cases h : K i j = 0
  · simp [h]
  · rw [(hK i j h).2.1, mul_zero]

lemma K_K (hK : KOk E K) : K * K = 0 := by
  ext i j
  rw [Matrix.mul_apply, Matrix.zero_apply]
  apply Finset.sum_eq_zero
  intro k _
  by_cases h1 : K i k = 0
  · simp [h1]
  by_cases h2 : K k j = 0
  · simp [h2]
  exfalso
  have a := (hK i k h1).2.2
  have b := (hK k j h2).1
  rw [a] at b
  norm_num at b

lemma PS_eq (E : Finset Sl) : (PE E + PR E : Matrix (L3 T Sl) (L3 T Sl) ℚ) = 1 - PX - PY := by
  rw [← sum_P E]; abel

lemma PS_K (hK : KOk E K) : (PE E + PR E) * K = K := by
  rw [PS_eq E, sub_mul, sub_mul, Matrix.one_mul, PX_K hK, PY_K hK, sub_zero, sub_zero]

lemma K_PSX (hK : KOk E K) : K * (PE E + PR E) + K * PX = K := by
  rw [PS_eq E, Matrix.mul_sub, Matrix.mul_sub, Matrix.mul_one, K_PY hK]
  abel

/-- **Forward split from the tail shape.** -/
theorem hub_gsplit (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) (hK : KOk E K) :
    (PX + PY + (PE E + PR E) * (M * (1 + K))) * (1 - K) = g8 M := by
  unfold g8
  rw [← PS_eq E]
  exact ring_g _ _ _ _ _ (PX_K hK) (PY_K hK) (K_K hK)

/-- **Backward (inverse-transpose) split from the tail shape.** -/
theorem hub_hsplit (Mt : Matrix (L3 T Sl) (L3 T Sl) ℚ) (hK : KOk E K) (hxm : PX * Mt = PX) :
    (PX + PY + ((1 - K) * Mt)ᵀ * (PE E + PR E)) * (1 + Kᵀ) = (h8 Mt)ᵀ := by
  have hA : Mtᵀ * PX = PX := by
    have h := congrArg Matrix.transpose hxm
    rwa [Matrix.transpose_mul, PX_T] at h
  have hPYB : PY * Kᵀ = 0 := by
    have h := congrArg Matrix.transpose (K_PY hK)
    rwa [Matrix.transpose_mul, PY_T, Matrix.transpose_zero] at h
  have hBP : Kᵀ * (PE E + PR E) = Kᵀ := by
    have h := congrArg Matrix.transpose (PS_K hK)
    rwa [Matrix.transpose_mul, Matrix.transpose_add, PE_T, PR_T] at h
  have hPB : (PE E + PR E) * Kᵀ + PX * Kᵀ = Kᵀ := by
    have h := congrArg Matrix.transpose (K_PSX hK)
    rwa [Matrix.transpose_add, Matrix.transpose_mul, Matrix.transpose_mul, Matrix.transpose_add, PE_T, PR_T,
      PX_T] at h
  have hBB : Kᵀ * Kᵀ = 0 := by
    rw [← Matrix.transpose_mul, K_K hK, Matrix.transpose_zero]
  have hT : (h8 Mt)ᵀ = PX + PY + Mtᵀ * (PE E + PR E) := by
    unfold h8
    rw [← PS_eq E, Matrix.transpose_add, Matrix.transpose_add, Matrix.transpose_mul, Matrix.transpose_add,
      PX_T, PY_T, PE_T, PR_T]
  rw [hT, Matrix.transpose_mul, Matrix.transpose_sub, Matrix.transpose_one]
  exact ring_h _ _ _ _ _ hA hPYB hBP hPB hBB

lemma hub_Gt_legal {β : Type} (lab : L3 T Sl → β) (hKlab : ∀ i j, K i j ≠ 0 → lab i = lab j) :
    ∀ i j, i ≠ j → (1 - K) i j ≠ 0 → lab i = lab j := by
  intro i j e h
  apply hKlab
  rw [Matrix.sub_apply, Matrix.one_apply_ne e, zero_sub, neg_ne_zero] at h
  exact h

lemma hub_Ht_legal {β : Type} (lab : L3 T Sl → β) (hKlab : ∀ i j, K i j ≠ 0 → lab i = lab j) :
    ∀ i j, i ≠ j → (1 + Kᵀ) i j ≠ 0 → lab i = lab j := by
  intro i j e h
  rw [Matrix.add_apply, Matrix.one_apply_ne e, zero_add, Matrix.transpose_apply] at h
  exact (hKlab j i h).symm

lemma gr_off (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) {i j : L3 T Sl} (e : i ≠ j) :
    (PX + PY + (PE E + PR E) * M : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j = (iE E i + iR E i) * M i j := by
  simp only [PX, PY, PE, PR, Matrix.add_mul, Matrix.add_apply, Matrix.diagonal_mul, Matrix.diagonal_apply_ne _ e]
  ring

lemma hr_off (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) {i j : L3 T Sl} (e : i ≠ j) :
    (PX + PY + Mᵀ * (PE E + PR E) : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j = M j i * (iE E j + iR E j) := by
  simp only [PX, PY, PE, PR, Matrix.mul_add, Matrix.add_apply, Matrix.mul_diagonal, Matrix.diagonal_apply_ne _ e,
    Matrix.transpose_apply]
  ring

lemma slot_of_iER {i : L3 T Sl} (h : iE E i + iR E i ≠ 0) : ∃ q, i = lS q := by
  rcases i with ((t|t)|q)
  · simp [iE, iR] at h
  · simp [iE, iR] at h
  · exact ⟨q, rfl⟩

lemma restOk_of {q : Sl} (h : iE E (lS q : L3 T Sl) ≠ 1) : restOk E (lS q : L3 T Sl) := by
  refine ⟨by simp [iY], ?_⟩
  by_cases hq : q ∈ E
  · exact absurd (by simp [iE, hq]) h
  · simp [iE, hq]

/-- the forward rest gate of the tail shape touches neither y roles nor exit slots. -/
theorem hub_Gr_legal (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) (hK : KOk E K) (hiy : M * PY = PY)
    (hRrow : ∀ i j, iE E i = 1 → (M * (1 + K)) i j = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j)
    (hRcol : ∀ i j, iE E j = 1 → (M * (1 + K)) i j = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j) :
    ∀ i j, i ≠ j → (PX + PY + (PE E + PR E) * (M * (1 + K)) : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j ≠ 0 → restOk E i ∧ restOk E j := by
  intro i j e h
  have hRy : (M * (1 + K)) * PY = PY := by
    rw [Matrix.mul_assoc, Matrix.add_mul, Matrix.one_mul, K_PY hK, add_zero, hiy]
  rw [gr_off (E := E) _ e] at h
  obtain ⟨h1, h2⟩ := mul_ne_zero_iff.mp h
  have hi_ne : iE E i ≠ 1 := fun hi => h2 (by rw [hRrow i j hi, Matrix.one_apply_ne e])
  have hj_ne : iE E j ≠ 1 := fun hj => h2 (by rw [hRcol i j hj, Matrix.one_apply_ne e])
  obtain ⟨q, rfl⟩ := slot_of_iER h1
  refine ⟨restOk_of hi_ne, ?_⟩
  rcases j with ((t|t)|q')
  · exact ⟨by simp [iY], by simp [iE]⟩
  · exact absurd (col_y_zero hRy q t) h2
  · exact restOk_of hj_ne

/-- the backward rest gate of the tail shape touches neither y roles nor exit slots. -/
theorem hub_Hr_legal (M Mt : Matrix (L3 T Sl) (L3 T Sl) ℚ) (hK : KOk E K) (hiy : M * PY = PY)
    (him : M * Mt = 1)
    (hRrow : ∀ i j, iE E i = 1 → (M * (1 + K)) i j = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j)
    (hRcol : ∀ i j, iE E j = 1 → (M * (1 + K)) i j = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j) :
    ∀ i j, i ≠ j → (PX + PY + ((1 - K) * Mt)ᵀ * (PE E + PR E) : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j ≠ 0 → restOk E i ∧ restOk E j := by
  intro i j e h
  set R := M * (1 + K) with hR
  set Ri := (1 - K) * Mt with hRi
  have hRRi : R * Ri = 1 := by
    have e1 : R * Ri = M * ((1 + K) * (1 - K)) * Mt := by
      rw [hR, hRi]; simp only [Matrix.mul_assoc]
    have e2 : (1 + K) * (1 - K) = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) := by
      rw [Matrix.add_mul, Matrix.one_mul, Matrix.mul_sub, Matrix.mul_one, K_K hK, sub_zero]; abel
    rw [e1, e2, Matrix.mul_one, him]
  have hRiR : Ri * R = 1 := mul_eq_one_comm.mp hRRi
  have hRy : R * PY = PY := by
    rw [hR, Matrix.mul_assoc, Matrix.add_mul, Matrix.one_mul, K_PY hK, add_zero, hiy]
  have hRiy : Ri * PY = PY := by
    calc Ri * PY = Ri * (R * PY) := by rw [hRy]
      _ = PY := by rw [← Matrix.mul_assoc, hRiR, Matrix.one_mul]
  -- rows and columns of `Ri` at exit slots
  have hrow : ∀ a b, iE E a = 1 → Ri a b = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) a b := by
    intro a b ha
    have h1 := congrFun (congrFun hRRi a) b
    rw [Matrix.mul_apply, Finset.sum_eq_single a] at h1
    · rw [hRrow a a ha, Matrix.one_apply_eq, one_mul] at h1; exact h1
    · intro k _ hk; rw [hRrow a k ha, Matrix.one_apply_ne (Ne.symm hk), zero_mul]
    · intro hn; exact absurd (Finset.mem_univ a) hn
  have hcol : ∀ a b, iE E b = 1 → Ri a b = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) a b := by
    intro a b hb
    have h1 := congrFun (congrFun hRiR a) b
    rw [Matrix.mul_apply, Finset.sum_eq_single b] at h1
    · rw [hRcol b b hb, Matrix.one_apply_eq, mul_one] at h1; exact h1
    · intro k _ hk; rw [hRcol k b hb, Matrix.one_apply_ne hk, mul_zero]
    · intro hn; exact absurd (Finset.mem_univ b) hn
  rw [hr_off (E := E) _ e] at h
  obtain ⟨h2, h1⟩ := mul_ne_zero_iff.mp h
  have hj_ne : iE E j ≠ 1 := fun hj => h2 (by rw [hrow j i hj, Matrix.one_apply_ne (Ne.symm e)])
  have hi_ne : iE E i ≠ 1 := fun hi => h2 (by rw [hcol j i hi, Matrix.one_apply_ne (Ne.symm e)])
  obtain ⟨q, rfl⟩ := slot_of_iER h1
  refine ⟨?_, restOk_of hj_ne⟩
  rcases i with ((t|t)|q')
  · exact ⟨by simp [iY], by simp [iE]⟩
  · exact absurd (col_y_zero hRiy q t) h2
  · exact restOk_of hi_ne

end Hub

section Tail
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

/-- **The tail data of a hub word**: the tail matrix `K` (one cleanup gate per exit helper), its support, equal tail
labels of each helper and its partner, and the rest `Mti * (1 + K)` acting as the identity on exit helpers. -/
structure HubTail (S : GScal T Sl C) (E : Finset Sl) (labT : L3 T Sl → SLbl H) where
  K : Matrix (L3 T Sl) (L3 T Sl) ℚ
  hK : KOk E K
  hKlab : ∀ i j, K i j ≠ 0 → labT i = labT j
  hRrow : ∀ i j, iE E i = 1 → (S.toScal.Mti * (1 + K)) i j = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j
  hRcol : ∀ i j, iE E j = 1 → (S.toScal.Mti * (1 + K)) i j = (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) i j

namespace HubTail
variable {S : GScal T Sl C} {E : Finset Sl} {labT : L3 T Sl → SLbl H} (τ : HubTail S E labT)

def Gt : Matrix (L3 T Sl) (L3 T Sl) ℚ := 1 - τ.K
def Gr : Matrix (L3 T Sl) (L3 T Sl) ℚ := PX + PY + (PE E + PR E) * (S.toScal.Mti * (1 + τ.K))
def Ht : Matrix (L3 T Sl) (L3 T Sl) ℚ := 1 + τ.Kᵀ
def Hr : Matrix (L3 T Sl) (L3 T Sl) ℚ := PX + PY + ((1 - τ.K) * S.toScal.Mt)ᵀ * (PE E + PR E)

/-- the six end-gate fields of `GLabX`, from the tail data. -/
theorem hGsplit : τ.Gr * τ.Gt = g8 S.toScal.Mti := hub_gsplit _ τ.hK
theorem hHsplit : τ.Hr * τ.Ht = (h8 S.toScal.Mt)ᵀ := hub_hsplit _ τ.hK S.toScal.hxm
theorem hGt : ∀ i j, i ≠ j → τ.Gt i j ≠ 0 → labT i = labT j := hub_Gt_legal labT τ.hKlab
theorem hHt : ∀ i j, i ≠ j → τ.Ht i j ≠ 0 → labT i = labT j := hub_Ht_legal labT τ.hKlab
theorem hGr : ∀ i j, i ≠ j → τ.Gr i j ≠ 0 → restOk E i ∧ restOk E j :=
  hub_Gr_legal _ τ.hK S.toScal.hiy τ.hRrow τ.hRcol
theorem hHr : ∀ i j, i ≠ j → τ.Hr i j ≠ 0 → restOk E i ∧ restOk E j :=
  hub_Hr_legal _ _ τ.hK S.toScal.hiy S.toScal.im τ.hRrow τ.hRcol

end HubTail
end Tail
end
end XEX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XEX.hub_gsplit
#print axioms OAI.PowerSaving.XEX.hub_hsplit
#print axioms OAI.PowerSaving.XEX.hub_Gr_legal
#print axioms OAI.PowerSaving.XEX.hub_Hr_legal
#print axioms OAI.PowerSaving.XEX.HubTail.hGsplit
#print axioms OAI.PowerSaving.XEX.HubTail.hHsplit
#print axioms OAI.PowerSaving.XEX.HubTail.hGr
#print axioms OAI.PowerSaving.XEX.HubTail.hHr
