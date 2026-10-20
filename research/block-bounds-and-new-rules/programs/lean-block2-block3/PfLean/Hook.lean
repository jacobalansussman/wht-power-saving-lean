import Work.PfLean.TheoremR

/-!
# The hook certificate (pf-lean, 13 October 2026)

`X = Sky α (j+2) → ℂ` (functions of `j+2` digit columns).  `Alt` antisymmetrises the columns
`1 .. j+1`; `Dsub` is the coordinate subspace "column 0 = column 1"; the certificate is

  `hookN j = Alt (Dsub)`  (= the repeated-address part of `V ⊗ Λ^(j+1) V`, `V = ℂ^(Space α)`).

Proved here, for every `α` and `j`:
* `hook_kept`: every free letter of the engine (shift, address permutation, phase) keeps `hookN`;
* `hook_cap`: every direction `dir z`, `z ≠ 0`, raises by at most `(card/2) * C(card - 2, j)`,
  `card = 2 ^ |α|` (`80` at `|α| = 3, j = 3`; `4` at `|α| = 2, j = 1`).
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.PfLean
open Binary GF Matrix Module Submodule

section Hook
variable {α : Type*} [Fintype α] [DecidableEq α]

/-! ### Column permutations -/

/-- Extend a permutation of `Fin n` to `Fin (n+1)`, fixing `0`. -/
def lift1 {n : ℕ} (π : Equiv.Perm (Fin n)) : Equiv.Perm (Fin (n+1)) :=
  Equiv.Perm.decomposeFin.symm (0, π)

@[simp] lemma lift1_zero {n : ℕ} (π : Equiv.Perm (Fin n)) : lift1 π 0 = 0 :=
  Equiv.Perm.decomposeFin_symm_apply_zero _ _

@[simp] lemma lift1_succ {n : ℕ} (π : Equiv.Perm (Fin n)) (i : Fin n) :
    lift1 π i.succ = (π i).succ := by
  rw [lift1, Equiv.Perm.decomposeFin_symm_apply_succ, Equiv.swap_self, Equiv.refl_apply]

lemma sign_lift1 {n : ℕ} (π : Equiv.Perm (Fin n)) :
    Equiv.Perm.sign (lift1 π) = Equiv.Perm.sign π := by
  rw [lift1, Equiv.Perm.decomposeFin.symm_sign]; simp

lemma lift1_mul {n : ℕ} (π π' : Equiv.Perm (Fin n)) : lift1 (π * π') = lift1 π * lift1 π' := by
  ext i; cases i using Fin.cases <;> simp [Equiv.Perm.mul_apply]

lemma lift1_one {n : ℕ} : lift1 (1 : Equiv.Perm (Fin n)) = 1 := by
  ext i; cases i using Fin.cases <;> simp

lemma lift1_inv {n : ℕ} (π : Equiv.Perm (Fin n)) : lift1 π⁻¹ = (lift1 π)⁻¹ := by
  rw [eq_inv_iff_mul_eq_one, ← lift1_mul, inv_mul_cancel, lift1_one]

abbrev Xf (α : Type*) [Fintype α] [DecidableEq α] (f : ℕ) := Sky α f → ℂ

/-- Permute the columns: `(permV e v) σ = v (σ ∘ e)`. -/
noncomputable def permV {f : ℕ} (e : Equiv.Perm (Fin f)) : Xf α f →ₗ[ℂ] Xf α f :=
  LinearMap.funLeft ℂ ℂ (fun σ : Sky α f => σ ∘ e)

@[simp] lemma permV_apply {f : ℕ} (e : Equiv.Perm (Fin f)) (v : Xf α f) (σ : Sky α f) :
    permV e v σ = v (σ ∘ e) := rfl

lemma permV_permV {f : ℕ} (e e' : Equiv.Perm (Fin f)) (v : Xf α f) :
    permV e (permV e' v) = permV (e * e') v := by
  funext σ; simp only [permV_apply, Equiv.Perm.coe_mul]; rfl

/-- The engine's column-wise action commutes with permuting the columns. -/
lemma matAct_permV {f : ℕ} (M : CMat α) (e : Equiv.Perm (Fin f)) (v : Xf α f) :
    matAct f M (permV e v) = permV e (matAct f M v) := by
  funext τ
  simp only [permV_apply, matAct, Matrix.mulVec, dotProduct, digitProd_apply]
  refine Fintype.sum_equiv (Equiv.arrowCongr e.symm (Equiv.refl (Space α))) _ _ (fun σ => ?_)
  simp only [Equiv.arrowCongr_apply, Equiv.symm_symm, Equiv.coe_refl, Function.id_comp,
    Function.comp_apply]
  congr 1
  exact (Equiv.prod_comp e (fun l => M (τ l) (σ l))).symm

/-! ### The antisymmetriser -/

/-- Antisymmetrise the columns `1 .. j+1`. -/
noncomputable def Alt (j : ℕ) : Xf α (j+2) →ₗ[ℂ] Xf α (j+2) where
  toFun v := ∑ π : Equiv.Perm (Fin (j+1)), ((Equiv.Perm.sign π : ℤ) : ℂ) • permV (lift1 π) v
  map_add' v w := by simp [smul_add, Finset.sum_add_distrib]
  map_smul' c v := by simp [Finset.smul_sum, smul_smul, mul_comm]

lemma Alt_apply (j : ℕ) (v : Xf α (j+2)) :
    Alt j v = ∑ π : Equiv.Perm (Fin (j+1)), ((Equiv.Perm.sign π : ℤ) : ℂ) • permV (lift1 π) v :=
  rfl

lemma matAct_Alt (j : ℕ) (M : CMat α) (v : Xf α (j+2)) :
    matAct (j+2) M (Alt j v) = Alt j (matAct (j+2) M v) := by
  rw [← matL_apply, Alt_apply, map_sum, Alt_apply]
  simp only [map_smul, matL_apply, matAct_permV]

lemma sign_sq_C {n : ℕ} (t : Equiv.Perm (Fin n)) :
    ((Equiv.Perm.sign t : ℤ) : ℂ) * ((Equiv.Perm.sign t : ℤ) : ℂ) = 1 := by
  rw [← Int.cast_mul, ← Units.val_mul, Int.units_mul_self]; simp

lemma Alt_permV (j : ℕ) (t : Equiv.Perm (Fin (j+1))) (v : Xf α (j+2)) :
    Alt j (permV (lift1 t) v) = ((Equiv.Perm.sign t : ℤ) : ℂ) • Alt j v := by
  simp only [Alt_apply, permV_permV, ← lift1_mul, Finset.smul_sum, smul_smul]
  refine Fintype.sum_equiv (Equiv.mulRight t) _ _ (fun π => ?_)
  simp only [Equiv.coe_mulRight, Equiv.Perm.sign_mul, Units.val_mul, Int.cast_mul]
  congr 1
  calc ((Equiv.Perm.sign π : ℤ) : ℂ)
      = ((Equiv.Perm.sign π : ℤ) : ℂ) * (((Equiv.Perm.sign t : ℤ) : ℂ) *
          ((Equiv.Perm.sign t : ℤ) : ℂ)) := by rw [sign_sq_C, mul_one]
    _ = _ := by ring

lemma Alt_eq_zero_of_fixed (j : ℕ) (t : Equiv.Perm (Fin (j+1))) (ht : Equiv.Perm.sign t = -1)
    (v : Xf α (j+2)) (hv : permV (lift1 t) v = v) : Alt j v = 0 := by
  have h := Alt_permV j t v
  rw [hv, ht] at h
  funext σ
  have h2 := congrFun h σ
  simp only [Units.val_neg, Units.val_one, Int.cast_neg, Int.cast_one, neg_smul, one_smul,
    Pi.neg_apply] at h2
  simp only [Pi.zero_apply]
  linear_combination h2 / 2

/-- The unit vector at a tuple. -/
noncomputable def del {f : ℕ} (τ : Sky α f) : Xf α f := Pi.single τ 1

lemma permV_del {f : ℕ} (e : Equiv.Perm (Fin f)) (τ : Sky α f) :
    permV e (del τ) = del (τ ∘ e.symm) := by
  funext σ
  have h : (σ ∘ ⇑e = τ) ↔ (σ = τ ∘ ⇑e.symm) := by
    constructor
    · rintro rfl; funext l; simp
    · rintro rfl; funext l; simp
  simp only [permV_apply, del, Pi.single_apply, h]

lemma Alt_del_eq_zero (j : ℕ) (τ : Sky α (j+2)) (a b : Fin (j+1)) (hab : a ≠ b)
    (h : τ a.succ = τ b.succ) : Alt j (del τ) = 0 := by
  apply Alt_eq_zero_of_fixed j (Equiv.swap a b) (Equiv.Perm.sign_swap hab)
  have hs : (lift1 (Equiv.swap a b)).symm = lift1 (Equiv.swap a b) := by
    rw [← Equiv.Perm.inv_def, ← lift1_inv, Equiv.swap_inv]
  rw [permV_del, hs]
  congr 1
  funext i
  cases i using Fin.cases with
  | zero => simp
  | succ i =>
    simp only [Function.comp_apply, lift1_succ]
    rw [Equiv.swap_apply_def]
    split_ifs with h1 h2
    · subst h1; exact h.symm
    · subst h2; exact h
    · rfl

/-! ### The certificate -/

/-- Column 0 equals column 1. -/
def Dsub (j : ℕ) : Submodule ℂ (Xf α (j+2)) where
  carrier := {w | ∀ τ : Sky α (j+2), τ 0 ≠ τ 1 → w τ = 0}
  add_mem' ha hb τ h := by simp [ha τ h, hb τ h]
  zero_mem' τ h := rfl
  smul_mem' c w hw τ h := by simp [hw τ h]

/-- **The hook certificate.** -/
noncomputable def hookN (j : ℕ) : Submodule ℂ (Xf α (j+2)) := (Dsub j).map (Alt j)

lemma hookN_vanish (j : ℕ) {v : Xf α (j+2)} (hv : v ∈ hookN j) (τ : Sky α (j+2))
    (hτ : ∀ i : Fin (j+1), τ i.succ ≠ τ 0) : v τ = 0 := by
  obtain ⟨w, hw, rfl⟩ := hv
  simp only [Alt_apply, Finset.sum_apply, Pi.smul_apply, permV_apply, smul_eq_mul]
  apply Finset.sum_eq_zero; intro π _
  rw [hw _ ?_, mul_zero]
  simp only [Function.comp_apply, lift1_zero]
  rw [← Fin.succ_zero_eq_one, lift1_succ]; exact (hτ _).symm

lemma matAct_mono {f : ℕ} (M : CMat α) (g : Space α → Space α) (d : Space α → ℂ)
    (hM : ∀ x y, M x y = if g x = y then d x else 0) (v : Xf α f) (τ : Sky α f) :
    matAct f M v τ = (∏ l, d (τ l)) * v (g ∘ τ) := by
  simp only [matAct, Matrix.mulVec, dotProduct, digitProd_apply, hM]
  rw [Finset.sum_eq_single (g ∘ τ)]
  · congr 1; apply Finset.prod_congr rfl; intro l _; simp
  · intro σ _ hσ
    obtain ⟨l, hl⟩ : ∃ l, g (τ l) ≠ σ l := by
      by_contra hcon; push Not at hcon; exact hσ (funext fun l => (hcon l).symm)
    rw [Finset.prod_eq_zero (Finset.mem_univ l) (by simp [hl]), zero_mul]
  · intro h; exact absurd (Finset.mem_univ _) h

lemma Dsub_mono (j : ℕ) (M : CMat α) (g : Space α → Space α) (d : Space α → ℂ)
    (hM : ∀ x y, M x y = if g x = y then d x else 0) (hg : Function.Injective g)
    {w : Xf α (j+2)} (hw : w ∈ Dsub j) : matAct (j+2) M w ∈ Dsub (α := α) j := by
  intro τ hτ
  rw [matAct_mono M g d hM, hw (g ∘ τ) (fun h => hτ (hg h)), mul_zero]

lemma hookN_mono (j : ℕ) (M : CMat α) (g : Space α → Space α) (d : Space α → ℂ)
    (hM : ∀ x y, M x y = if g x = y then d x else 0) (hg : Function.Injective g)
    {v : Xf α (j+2)} (hv : v ∈ hookN j) : matAct (j+2) M v ∈ hookN (α := α) j := by
  obtain ⟨w, hw, rfl⟩ := hv
  exact ⟨_, Dsub_mono j M g d hM hg hw, (matAct_Alt j M w).symm⟩

/-- **Every free letter keeps the hook certificate.** -/
theorem hook_kept (j : ℕ) : Kept (hookN (α := α) j) where
  shift z v hv := hookN_mono j _ (· + z) (fun _ => 1) (fun x y => by simp [Binary.shift])
    (add_left_injective z) hv
  perm G v hv := hookN_mono j _ G.π (fun _ => 1) (fun x y => by simp [permMat])
    G.π.injective hv
  phase z c v hv := hookN_mono j _ id (fun x => tint (dot z x + c))
    (fun x y => by simp [phaseMat, Matrix.diagonal_apply]) Function.injective_id hv

/-! ### The raise of a direction -/

lemma dir_apply (z a y : Space α) :
    dir z a y = (if a = y then ca else 0) + (if a + z = y then cb else 0) := by
  simp [dir, Binary.shift, Matrix.one_apply]

lemma space_add_self (z : Space α) : z + z = 0 := by
  funext i
  show z i + z i = 0
  generalize z i = c
  revert c; decide

lemma dir_bracket (z : Space α) (hz : z ≠ 0) (x a b : Space α) (hab : a ≠ b) :
    dir z a x * dir z b x = dir z a (x + z) * dir z b (x + z) := by
  have hzz := space_add_self z
  have key : ∀ c y : Space α, (c + z = y ↔ c = y + z) := fun c y =>
    ⟨fun h => by rw [← h, add_assoc, hzz, add_zero], fun h => by rw [h, add_assoc, hzz, add_zero]⟩
  have hxz : x ≠ x + z := fun h => hz (add_left_cancel (h.symm.trans (add_zero x).symm))
  simp only [dir_apply, key, add_assoc, hzz, add_zero]
  by_cases h1 : a = x <;> by_cases h2 : a = x + z <;> by_cases h3 : b = x <;>
    by_cases h4 : b = x + z <;> simp_all [mul_comm]

/-- `(x, x, ρ)`. -/
def tup {j : ℕ} (x : Space α) (ρ : Fin j → Space α) : Sky α (j+2) :=
  Fin.cons x (Fin.cons x ρ : Fin (j+1) → Space α)

@[simp] lemma tup_zero {j : ℕ} (x : Space α) (ρ : Fin j → Space α) : tup x ρ 0 = x := rfl
@[simp] lemma tup_one {j : ℕ} (x : Space α) (ρ : Fin j → Space α) :
    tup x ρ (Fin.succ 0) = x := rfl
@[simp] lemma tup_one' {j : ℕ} (x : Space α) (ρ : Fin j → Space α) :
    tup x ρ 1 = x := by rw [← Fin.succ_zero_eq_one]; rfl
@[simp] lemma tup_ss {j : ℕ} (x : Space α) (ρ : Fin j → Space α) (i : Fin j) :
    tup x ρ i.succ.succ = ρ i := rfl

lemma del_mem_Dsub (j : ℕ) (τ : Sky α (j+2)) (h : τ 0 = τ 1) : del τ ∈ Dsub (α := α) j := by
  intro σ hσ
  simp only [del, Pi.single_apply]
  rw [ite_eq_right_iff.mpr]; rintro rfl; exact absurd h hσ

lemma tup_01 {j : ℕ} (x : Space α) (ρ : Fin j → Space α) : tup x ρ 0 = tup x ρ 1 := by
  rw [← Fin.succ_zero_eq_one]; rfl

/-- The pair vectors lie in `N ∩ (dir z)⁻¹ N`. -/
lemma pair_mem (j : ℕ) (z : Space α) (hz : z ≠ 0) (x : Space α) (ρ : Fin j → Space α) :
    Alt j (del (tup x ρ)) - Alt j (del (tup (x + z) ρ)) ∈
      hookN j ⊓ (hookN j).comap (matL (j+2) (dir z)) := by
  have hD : del (tup x ρ) - del (tup (x + z) ρ) ∈ Dsub (α := α) j :=
    Submodule.sub_mem _ (del_mem_Dsub j _ (tup_01 _ _)) (del_mem_Dsub j _ (tup_01 _ _))
  have hD' : matAct (j+2) (dir z) (del (tup x ρ) - del (tup (x + z) ρ)) ∈ Dsub (α := α) j := by
    intro σ hσ
    simp only [matAct, Matrix.mulVec_sub, Pi.sub_apply, del, Matrix.mulVec_single_one,
      Matrix.col_apply, digitProd_apply]
    simp only [Fin.prod_univ_succ, tup_zero, tup_one, tup_ss]
    have hb := dir_bracket z hz x (σ 0) (σ (Fin.succ 0)) (by rwa [Fin.succ_zero_eq_one])
    linear_combination (∏ i : Fin j, dir z (σ i.succ.succ) (ρ i)) * hb
  refine ⟨⟨_, hD, by rw [map_sub]⟩, ?_⟩
  show matL (j+2) (dir z) _ ∈ hookN j
  rw [matL_apply, ← map_sub, matAct_Alt]
  exact ⟨_, hD', rfl⟩

lemma tup_comp {j : ℕ} (a : Space α) (r : Fin j → Space α) (π : Equiv.Perm (Fin j)) :
    tup a (r ∘ π) = tup a r ∘ lift1 (lift1 π) := by
  funext l
  cases l using Fin.cases with
  | zero => simp
  | succ l =>
    cases l using Fin.cases with
    | zero =>
      simp only [Function.comp_apply, lift1_succ, lift1_zero]; rfl
    | succ l => simp

/-- A tuple listing a `j`-set (any fixed choice). -/
noncomputable def rep (j : ℕ) (T : Finset (Space α)) : Fin j → Space α :=
  if h : ∃ ρ : Fin j → Space α, Function.Injective ρ ∧ Finset.univ.image ρ = T then h.choose
  else fun _ => 0

lemma rep_perm (j : ℕ) (ρ : Fin j → Space α) (hρ : Function.Injective ρ) :
    ∃ π : Equiv.Perm (Fin j), ρ = rep j (Finset.univ.image ρ) ∘ π := by
  have hex : ∃ ρ' : Fin j → Space α, Function.Injective ρ' ∧
      Finset.univ.image ρ' = Finset.univ.image ρ := ⟨ρ, hρ, rfl⟩
  have hr : rep j (Finset.univ.image ρ) = hex.choose := dif_pos hex
  obtain ⟨hinj, himg⟩ := hex.choose_spec
  have hmem : ∀ i, ∃ i', hex.choose i' = ρ i := fun i => by
    have : ρ i ∈ Finset.univ.image hex.choose := by
      rw [himg]; exact Finset.mem_image_of_mem _ (Finset.mem_univ i)
    simpa using this
  choose π hπ using hmem
  have πinj : Function.Injective π := fun a b h => hρ (by rw [← hπ a, ← hπ b, h])
  refine ⟨Equiv.ofBijective π (Finite.injective_iff_bijective.mp πinj), ?_⟩
  rw [hr]; funext i; simp [hπ]

/-- The transversal of the pairs `{a, a + z}`. -/
def Tr (i0 : α) : Finset (Space α) := Finset.univ.filter (fun a => a i0 = 0)

/-- Index set of the `80` (resp. `4`) complement vectors. -/
def Idx (j : ℕ) (z : Space α) (i0 : α) : Finset (Σ _ : Space α, Finset (Space α)) :=
  (Tr i0).sigma (fun a => (Finset.univ \ {a, a + z}).powersetCard j)

noncomputable def Fv (j : ℕ) (p : Σ _ : Space α, Finset (Space α)) : Xf α (j+2) :=
  Alt j (del (tup p.1 (rep j p.2)))

noncomputable def Qsp (j : ℕ) (z : Space α) (i0 : α) : Submodule ℂ (Xf α (j+2)) :=
  span ℂ (Set.range (fun p : (Idx j z i0) => Fv j p.1))

lemma finrank_Qsp (j : ℕ) (z : Space α) (i0 : α) : finrank ℂ (Qsp j z i0) ≤ (Idx j z i0).card := by
  have h1 := finrank_range_le_card (R := ℂ) (fun p : (Idx j z i0) => Fv j p.1)
  exact h1.trans (Fintype.card_coe _).le

lemma zmod2_ne (c : ZMod 2) (h : c ≠ 0) : c = 1 := by revert c; decide
lemma zmod2_add_one (c : ZMod 2) (h : c ≠ 0) : c + 1 = 0 := by revert c; decide
lemma zmod2_add_one' (c : ZMod 2) (h : c = 0) : c + 1 ≠ 0 := by revert c; decide

/-- Every generator lies in `(N ∩ (dir z)⁻¹ N) + Q`. -/
lemma key_mem (j : ℕ) (z : Space α) (hz : z ≠ 0) (i0 : α) (hi0 : z i0 = 1)
    (x : Space α) (ρ : Fin j → Space α) :
    Alt j (del (tup x ρ)) ∈ (hookN j ⊓ (hookN j).comap (matL (j+2) (dir z))) ⊔ Qsp j z i0 := by
  have hzz := space_add_self z
  -- case A: ρ not injective
  by_cases hA : ¬ Function.Injective ρ
  · obtain ⟨a, b, hab, hne⟩ : ∃ a b, ρ a = ρ b ∧ a ≠ b := by
      unfold Function.Injective at hA; push Not at hA; exact hA
    rw [Alt_del_eq_zero j _ a.succ b.succ (fun h => hne (Fin.succ_injective _ h)) (by simpa using hab)]
    exact Submodule.zero_mem _
  push Not at hA
  -- case B: x occurs in ρ
  by_cases hB : ∃ i, ρ i = x
  · obtain ⟨i, hi⟩ := hB
    rw [Alt_del_eq_zero j _ 0 i.succ (Fin.succ_ne_zero i).symm (by simp [hi])]
    exact Submodule.zero_mem _
  push Not at hB
  -- case C: x + z occurs in ρ
  by_cases hC : ∃ i, ρ i = x + z
  · obtain ⟨i, hi⟩ := hC
    have h0 : Alt j (del (tup (x + z) ρ)) = 0 :=
      Alt_del_eq_zero j _ 0 i.succ (Fin.succ_ne_zero i).symm (by simp [hi])
    have := pair_mem j z hz x ρ
    rw [h0, sub_zero] at this
    exact Submodule.mem_sup_left this
  push Not at hC
  -- case D
  set a : Space α := if x i0 = 0 then x else x + z with ha_def
  have ha : a i0 = 0 := by
    rw [ha_def]; split_ifs with h
    · exact h
    · show x i0 + z i0 = 0; rw [hi0]; exact zmod2_add_one _ h
  have hpair : (a = x ∧ a + z = x + z) ∨ (a = x + z ∧ a + z = x) := by
    rw [ha_def]; split_ifs with h
    · exact Or.inl ⟨rfl, rfl⟩
    · exact Or.inr ⟨rfl, by rw [add_assoc, hzz, add_zero]⟩
  set T := Finset.univ.image ρ with hT
  have hidx : (⟨a, T⟩ : Σ _ : Space α, Finset (Space α)) ∈ Idx j z i0 := by
    simp only [Idx, Finset.mem_sigma, Tr, Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.mem_powersetCard]
    refine ⟨ha, ?_, ?_⟩
    · intro t ht
      rw [hT, Finset.mem_image] at ht
      obtain ⟨i, _, rfl⟩ := ht
      simp only [Finset.mem_sdiff, Finset.mem_univ, true_and, Finset.mem_insert,
        Finset.mem_singleton, not_or]
      rcases hpair with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · rw [h1]; exact ⟨hB i, hC i⟩
      · rw [h1]; exact ⟨hC i, by rw [add_assoc, hzz, add_zero]; exact hB i⟩
    · rw [hT, Finset.card_image_of_injective _ hA, Finset.card_univ, Fintype.card_fin]
  obtain ⟨π, hπ⟩ := rep_perm j ρ hA
  rw [← hT] at hπ
  have hQ : Alt j (del (tup a ρ)) ∈ Qsp j z i0 := by
    have e1 : del (tup a ρ) = permV (lift1 (lift1 π)⁻¹) (del (tup a (rep j T))) := by
      rw [permV_del, lift1_inv, ← Equiv.Perm.inv_def, inv_inv, ← tup_comp, ← hπ]
    rw [e1, Alt_permV]
    exact Submodule.smul_mem _ _ (Submodule.subset_span ⟨⟨⟨a, T⟩, hidx⟩, rfl⟩)
  have hU : Alt j (del (tup x ρ)) - Alt j (del (tup a ρ)) ∈
      hookN j ⊓ (hookN j).comap (matL (j+2) (dir z)) := by
    rcases hpair with ⟨h1, _⟩ | ⟨h1, _⟩
    · rw [h1, sub_self]; exact Submodule.zero_mem _
    · rw [h1]; exact pair_mem j z hz x ρ
  have := Submodule.add_mem_sup hU hQ
  rwa [sub_add_cancel] at this

lemma hookN_le (j : ℕ) (z : Space α) (hz : z ≠ 0) (i0 : α) (hi0 : z i0 = 1) :
    hookN j ≤ (hookN j ⊓ (hookN j).comap (matL (j+2) (dir z))) ⊔ Qsp j z i0 := by
  rintro _ ⟨w, hw, rfl⟩
  have hw' : w = ∑ τ, w τ • del τ := by
    funext σ; simp [del, Finset.sum_apply, Pi.single_apply]
  rw [hw', map_sum]
  apply Submodule.sum_mem; intro τ _
  rw [map_smul]
  by_cases h : τ 0 = τ 1
  · apply Submodule.smul_mem
    have hτ : τ = tup (τ 0) (fun i => τ i.succ.succ) := by
      funext l
      cases l using Fin.cases with
      | zero => rfl
      | succ l =>
        cases l using Fin.cases with
        | zero => show τ (Fin.succ 0) = τ 0; rw [Fin.succ_zero_eq_one, h]
        | succ l => rfl
    rw [hτ]; exact key_mem j z hz i0 hi0 _ _
  · rw [hw τ h, zero_smul]; exact Submodule.zero_mem _

lemma card_Tr (z : Space α) (i0 : α) (hi0 : z i0 = 1) :
    (Tr i0).card * 2 = Fintype.card (Space α) := by
  have hzz := space_add_self z
  have h1 := Finset.card_filter_add_card_filter_not (s := Finset.univ)
    (fun a : Space α => a i0 = 0)
  have h2 : (Finset.univ.filter (fun a : Space α => a i0 = 0)).card =
      (Finset.univ.filter (fun a : Space α => ¬ a i0 = 0)).card := by
    apply Finset.card_nbij' (· + z) (· + z)
    · intro a ha
      simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_setOf_eq] at ha ⊢
      show a i0 + z i0 ≠ 0; rw [hi0]; exact zmod2_add_one' _ ha
    · intro a ha
      simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_setOf_eq] at ha ⊢
      show a i0 + z i0 = 0; rw [hi0]; exact zmod2_add_one _ ha
    · intro a _; simp [add_assoc, hzz]
    · intro a _; simp [add_assoc, hzz]
  rw [Finset.card_univ] at h1
  unfold Tr; omega

lemma card_Idx (j : ℕ) (z : Space α) (hz : z ≠ 0) (i0 : α) :
    (Idx j z i0).card = (Tr i0).card * (Fintype.card (Space α) - 2).choose j := by
  rw [Idx, Finset.card_sigma, Finset.sum_const_nat]
  intro a _
  rw [Finset.card_powersetCard, ← Finset.compl_eq_univ_sdiff, Finset.card_compl,
    Finset.card_pair]
  intro h
  exact hz (add_left_cancel (h.symm.trans (add_zero a).symm))

/-- **The raise of a direction.** -/
theorem hook_cap (j : ℕ) (z : Space α) (hz : z ≠ 0) :
    raise (hookN (α := α) j) (matL (j+2) (dir z)) ≤
      (Fintype.card (Space α) / 2) * (Fintype.card (Space α) - 2).choose j := by
  obtain ⟨i0, hi0⟩ : ∃ i0, z i0 ≠ 0 := Function.ne_iff.mp hz
  have hi0' : z i0 = 1 := zmod2_ne _ hi0
  have hc := card_Tr z i0 hi0'
  have h1 := raise_le_of_le_sup (le_refl _) (hookN_le j z hz i0 hi0')
  have h2 := finrank_Qsp j z i0
  rw [card_Idx j z hz i0] at h2
  have h3 : Fintype.card (Space α) / 2 = (Tr i0).card := by omega
  rw [h3]; omega

end Hook
end PowerSaving.PfLean
end OAI

#print axioms OAI.PowerSaving.PfLean.hook_kept
#print axioms OAI.PowerSaving.PfLean.hook_cap
