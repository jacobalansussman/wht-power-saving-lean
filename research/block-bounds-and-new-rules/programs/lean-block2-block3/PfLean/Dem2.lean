import Work.PfLean.Hook

/-!
# The demand of the block-2 hook (pf-lean, 13 October 2026)

Block size 2 (`α = Fin 2`, `V = ℂ^4`), `f = 3` columns, certificate `hookN 1` (`V ⊗ Λ²V`,
repeated address).  We prove `raise N kernel ≥ 8`: eight vectors `Alt (δ (y, y, t))` of `N` whose
images under the kernel are linearly independent modulo `N`, witnessed by the values at eight
tuples `(x, a, b)` with `x ∉ {a, b}` (where every vector of `N` vanishes).  The `8 × 8` table of
values is computed entry by entry from the definitions; its inverse is given explicitly
(entries `±1, ±i`), so no rank computation is trusted.
-/

set_option linter.unusedSectionVars false
set_option maxHeartbeats 8000000

namespace OAI
namespace PowerSaving.PfLean
open Binary GF Matrix Module Submodule Complex

/-- The four addresses (bit `a` of the integer). -/
def ad : Fin 4 → Space (Fin 2) := ![![0, 0], ![1, 0], ![0, 1], ![1, 1]]

/-- Eight tuples `(x, a, b)`, `x ∉ {a, b}`. -/
def rowT : Fin 8 → Sky (Fin 2) 3 := ![![ad 0, ad 1, ad 2], ![ad 0, ad 1, ad 3], ![ad 0, ad 2, ad 3],
  ![ad 1, ad 0, ad 2], ![ad 1, ad 0, ad 3], ![ad 1, ad 2, ad 3], ![ad 2, ad 0, ad 3],
  ![ad 2, ad 1, ad 3]]

/-- Eight tuples `(y, y, t)`. -/
def colS : Fin 8 → Sky (Fin 2) 3 := ![![ad 0, ad 0, ad 1], ![ad 0, ad 0, ad 2], ![ad 0, ad 0, ad 3],
  ![ad 1, ad 1, ad 0], ![ad 1, ad 1, ad 2], ![ad 1, ad 1, ad 3], ![ad 2, ad 2, ad 0],
  ![ad 2, ad 2, ad 1]]

/-- The table of values `(kernel (Alt δ(colS c))) (rowT r)`. -/
noncomputable def Etab : Fin 8 → Fin 8 → ℂ := ![
  ![1/4, -1/4, 0, 1/4*I, 0, 1/4*I, -1/4*I, 0],
  ![0, 1/4*I, -1/4, 0, 1/4*I, -1/4, -1/4, -1/4*I],
  ![1/4*I, 0, -1/4, -1/4, -1/4*I, 0, 0, 1/4*I],
  ![0, -1/4, 1/4*I, 0, -1/4, 1/4*I, -1/4*I, -1/4],
  ![1/4*I, 1/4*I, 0, 1/4, 0, -1/4, -1/4, 0],
  ![1/4, 0, 1/4*I, -1/4*I, 1/4, 0, 0, 1/4],
  ![1/4*I, 1/4*I, 0, -1/4, 0, 1/4, 1/4, 0],
  ![0, 1/4, 1/4*I, 0, 1/4, 1/4*I, -1/4*I, 1/4]]

lemma I3 : I ^ 3 = -I := by rw [pow_succ, I_sq]; ring
lemma I4 : I ^ 4 = 1 := by rw [show (4:ℕ) = 2 * 2 from rfl, pow_mul, I_sq]; norm_num
lemma I5 : I ^ 5 = I := by rw [pow_succ, I4, one_mul]
lemma I6 : I ^ 6 = -1 := by rw [pow_succ, I5, I_mul_I]
lemma I7 : I ^ 7 = -I := by rw [pow_succ, I6]; ring

lemma alt1_apply (v : Xf (Fin 2) 3) (τ : Sky (Fin 2) 3) :
    Alt 1 v τ = v τ - v ![τ 0, τ 2, τ 1] := by
  have huniv : (Finset.univ : Finset (Equiv.Perm (Fin 2))) = {1, Equiv.swap 0 1} := by decide
  rw [Alt_apply, Finset.sum_apply, huniv, Finset.sum_pair (by decide)]
  have h1 : τ ∘ ⇑(lift1 (1 : Equiv.Perm (Fin 2))) = τ := by rw [lift1_one]; rfl
  have h2 : τ ∘ ⇑(lift1 (Equiv.swap (0 : Fin 2) 1)) = ![τ 0, τ 2, τ 1] := by
    funext l
    fin_cases l
    · simp
    · show τ (lift1 (Equiv.swap (0 : Fin 2) 1) (Fin.succ 0)) = τ 2
      rw [lift1_succ, Equiv.swap_apply_left]; rfl
    · show τ (lift1 (Equiv.swap (0 : Fin 2) 1) (Fin.succ 1)) = τ 1
      rw [lift1_succ, Equiv.swap_apply_right]; rfl
  simp only [Pi.smul_apply, permV_apply, smul_eq_mul, h1, h2, Equiv.Perm.sign_one,
    Equiv.Perm.sign_swap (show (0 : Fin 2) ≠ 1 by decide)]
  simp; ring

lemma matAct_del (M : CMat (Fin 2)) (σ τ : Sky (Fin 2) 3) :
    matAct 3 M (del σ) τ = ∏ l, M (τ l) (σ l) := by
  simp [matAct, del, Matrix.mulVec_single_one, Matrix.col_apply, digitProd_apply]

lemma kernel2_apply (x y : Space (Fin 2)) :
    kernel (Fin 2) x y = littleC (x 0) (y 0) * littleC (x 1) (y 1) := by
  simp [kernel, digitProd_apply, Fin.prod_univ_two]

lemma lc00 : littleC 0 0 = ca := by simp [littleC]
lemma lc01 : littleC 0 1 = cb := by simp [littleC]
lemma lc10 : littleC 1 0 = cb := by simp [littleC]
lemma lc11 : littleC 1 1 = ca := by simp [littleC]

/-- The `64` values, computed from the definitions. -/
lemma entry (r c : Fin 8) :
    matAct 3 (kernel (Fin 2)) (Alt 1 (del (colS c))) (rowT r) = Etab r c := by
  rw [matAct_Alt, alt1_apply, matAct_del, matAct_del]
  fin_cases r <;> fin_cases c <;>
    simp [rowT, colS, ad, Etab, Fin.prod_univ_three, kernel2_apply, lc00, lc01, lc10, lc11,
      ca, cb] <;>
    (try ring_nf) <;> (try simp only [I_sq, I3, I4, I5, I6, I7]) <;> (try ring_nf)

/-- **The demand of the block-2 hook is at least 8.** -/
theorem dem2 : 8 ≤ raise (hookN (α := Fin 2) 1) (matL 3 (kernel (Fin 2))) := by
  have hn : ∀ c, Alt 1 (del (colS c)) ∈ hookN (α := Fin 2) 1 := fun c =>
    ⟨_, del_mem_Dsub 1 _ (by revert c; decide), rfl⟩
  let y : Fin 8 → (Xf (Fin 2) 3 ⧸ hookN (α := Fin 2) 1) := fun c =>
    (hookN (α := Fin 2) 1).mkQ (matL 3 (kernel (Fin 2)) (Alt 1 (del (colS c))))
  have hli : LinearIndependent ℂ y := by
    rw [Fintype.linearIndependent_iff]
    intro a ha
    have hmem : matL 3 (kernel (Fin 2)) (∑ c, a c • Alt 1 (del (colS c))) ∈
        hookN (α := Fin 2) 1 := by
      rw [← Submodule.Quotient.mk_eq_zero, ← Submodule.mkQ_apply, map_sum, map_sum]
      simpa [y, map_smul] using ha
    have hr : ∀ r, ∑ c, a c * Etab r c = 0 := by
      intro r
      have h0 := hookN_vanish 1 hmem (rowT r) (by revert r; decide)
      rw [map_sum] at h0
      simpa only [map_smul, Finset.sum_apply, Pi.smul_apply, smul_eq_mul, matL_apply,
        entry] using h0
    have h0 := hr 0
    have h1 := hr 1
    have h2 := hr 2
    have h3 := hr 3
    have h4 := hr 4
    have h5 := hr 5
    have h6 := hr 6
    have h7 := hr 7
    simp only [Fin.sum_univ_eight, Etab, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.head_cons, Matrix.cons_val] at h0 h1 h2 h3 h4 h5 h6 h7
    have e0 : a 0 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        h0 - I * h4 + h5 - I * h6 - h7
    have e1 : a 1 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        -h0 - I * h4 - h5 - I * h6 + h7
    have e2 : a 2 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        I * h0 - h1 - h2 - I * h3 + h4 - I * h5
    have e3 : a 3 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        -I * h0 - h1 - h2 + I * h3 + h4 + I * h5
    have e4 : a 4 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        h0 - I * h1 + I * h2 - h3 + I * h4 + h5
    have e5 : a 5 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        -I * h0 - h4 + I * h5 + h6 - I * h7
    have e6 : a 6 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        -h1 - h2 + I * h3 + h6 + I * h7
    have e7 : a 7 = 0 := by
      linear_combination (norm := (ring_nf; try simp only [I_sq]; try ring_nf))
        I * h1 - I * h2 - h3 + I * h6 + h7
    intro c
    fin_cases c
    exacts [e0, e1, e2, e3, e4, e5, e6, e7]
  have hspan : span ℂ (Set.range y) ≤
      (hookN (α := Fin 2) 1).map ((hookN (α := Fin 2) 1).mkQ ∘ₗ matL 3 (kernel (Fin 2))) := by
    rw [Submodule.span_le]; rintro _ ⟨c, rfl⟩; exact ⟨_, hn c, rfl⟩
  have h := Submodule.finrank_mono hspan
  rw [finrank_span_eq_card hli, Fintype.card_fin] at h
  exact h

end PowerSaving.PfLean
end OAI

#print axioms OAI.PowerSaving.PfLean.dem2
