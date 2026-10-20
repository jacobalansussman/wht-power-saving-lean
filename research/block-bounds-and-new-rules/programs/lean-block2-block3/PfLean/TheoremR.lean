import Work.GFrame.Engine.Letters

/-!
# Theorem R for the GFrame engine (pf-lean, 13 October 2026)

A subspace `N` of the one-role data space `X = Sky α f → ℂ` that every FREE letter keeps gives a
lower bound on the PAID letters of every correct word:

  `card σ * raise N (kernel) ≤ cap * (total block rank of w)`

whenever every direction `dir z` (`z ≠ 0`) has `raise N (dir z) ≤ cap`.  Here
`raise N M = finrank (N.map (N.mkQ ∘ M))` = `dim (N + M N) - dim N` for invertible `M`.

Proof: the potential `pot A = dim (A + Nᵖ) - dim Nᵖ` (`Nᵖ` = `N` in every role, `A` = image of the
inputs) starts at `0`, is not raised by free letters, is raised by at most `raise` by a paid
letter, and ends at least `card σ * raise N kernel` (the live roles hold `kernel` of their input).
No frame, no Clifford structure, no invertibility of the mixing matrices is used.
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.PfLean
open Module Submodule

/-! ## Linear algebra -/
section Abstract
variable {K V W : Type*} [Field K] [AddCommGroup V] [Module K V] [FiniteDimensional K V]
  [AddCommGroup W] [Module K W]

/-- Rank-nullity on a subspace. -/
lemma finrank_map_add_inf_ker (g : V →ₗ[K] W) (B : Submodule K V) :
    finrank K (B.map g) + finrank K ↥(B ⊓ LinearMap.ker g) = finrank K B := by
  have h := LinearMap.finrank_range_add_finrank_ker (g.domRestrict B)
  rw [LinearMap.range_domRestrict, LinearMap.ker_domRestrict] at h
  have h2 : (LinearMap.ker g).comap B.subtype = (B ⊓ LinearMap.ker g).comap B.subtype := by
    rw [Submodule.comap_inf, Submodule.comap_subtype_self, top_inf_eq]
  rw [h2, (Submodule.comapSubtypeEquivOfLe inf_le_left).finrank_eq] at h
  exact h

/-- The raise of `L` against `P`: the dimension of the image of `P` in `V ⧸ P` under `L`. -/
noncomputable def raise (P : Submodule K V) (L : V →ₗ[K] V) : ℕ :=
  finrank K (P.map (P.mkQ ∘ₗ L))

/-- The potential of a subspace `A`: the dimension of its image in `V ⧸ P`. -/
noncomputable def pot (P A : Submodule K V) : ℕ := finrank K (A.map P.mkQ)

lemma ker_mkQ_comp (P : Submodule K V) (L : V →ₗ[K] V) :
    LinearMap.ker (P.mkQ ∘ₗ L) = P.comap L := by
  rw [LinearMap.ker_comp, Submodule.ker_mkQ]

lemma raise_add (P : Submodule K V) (L : V →ₗ[K] V) :
    raise P L + finrank K ↥(P ⊓ P.comap L) = finrank K P := by
  rw [raise, ← ker_mkQ_comp]; exact finrank_map_add_inf_ker _ _

lemma pot_eq_zero_of_le {P A : Submodule K V} (h : A ≤ P) : pot P A = 0 := by
  have : A.map P.mkQ = ⊥ := by
    rw [eq_bot_iff]
    rintro _ ⟨a, ha, rfl⟩
    simp [(Submodule.Quotient.mk_eq_zero P).mpr (h ha)]
  rw [pot, this, finrank_bot]

/-- **One step of Theorem R.** -/
lemma pot_map_le (P A : Submodule K V) (L : V →ₗ[K] V) :
    pot P (A.map L) ≤ pot P A + raise P L := by
  have h1 : (A.map L).map P.mkQ ≤ (A ⊔ (P ⊓ P.comap L)).map (P.mkQ ∘ₗ L) := by
    rw [← Submodule.map_comp]; exact Submodule.map_mono le_sup_left
  have h2 := finrank_map_add_inf_ker (P.mkQ ∘ₗ L) (A ⊔ (P ⊓ P.comap L))
  have h3 : finrank K ↥(P ⊓ P.comap L) ≤
      finrank K ↥((A ⊔ (P ⊓ P.comap L)) ⊓ LinearMap.ker (P.mkQ ∘ₗ L)) := by
    apply Submodule.finrank_mono
    refine le_inf le_sup_right ?_
    rw [ker_mkQ_comp]; exact inf_le_right
  have h4 : finrank K ↥(A ⊔ (P ⊓ P.comap L)) ≤ finrank K ↥(A ⊔ P) :=
    Submodule.finrank_mono (sup_le_sup_left inf_le_left _)
  have h5 := finrank_map_add_inf_ker P.mkQ (A ⊔ P)
  have h6 : (A ⊔ P).map P.mkQ = A.map P.mkQ := by
    rw [Submodule.map_sup, Submodule.mkQ_map_self, sup_bot_eq]
  have h7 : (A ⊔ P) ⊓ LinearMap.ker P.mkQ = P := by
    rw [Submodule.ker_mkQ]; exact inf_eq_right.mpr le_sup_right
  have h8 := raise_add P L
  have h9 : finrank K ((A.map L).map P.mkQ) ≤
      finrank K ((A ⊔ (P ⊓ P.comap L)).map (P.mkQ ∘ₗ L)) :=
    Submodule.finrank_mono h1
  rw [h6, h7] at h5
  unfold pot
  omega

lemma raise_eq_pot (P : Submodule K V) (L : V →ₗ[K] V) : raise P L = pot P (P.map L) := by
  rw [raise, pot, Submodule.map_comp]

lemma raise_comp_le (P : Submodule K V) (L₁ L₂ : V →ₗ[K] V) :
    raise P (L₂ ∘ₗ L₁) ≤ raise P L₁ + raise P L₂ := by
  rw [raise_eq_pot, Submodule.map_comp, raise_eq_pot P L₁]
  exact pot_map_le _ _ _

lemma raise_eq_zero_of_mapsTo {P : Submodule K V} {L : V →ₗ[K] V} (h : ∀ x ∈ P, L x ∈ P) :
    raise P L = 0 := by
  rw [raise_eq_pot]; apply pot_eq_zero_of_le
  rintro _ ⟨x, hx, rfl⟩; exact h x hx

lemma raise_le_of_le_sup {P U Q : Submodule K V} {L : V →ₗ[K] V} (hU : U ≤ P ⊓ P.comap L)
    (hP : P ≤ U ⊔ Q) : raise P L ≤ finrank K Q := by
  have h1 := raise_add P L
  have h2 : finrank K P ≤ finrank K ↥(U ⊔ Q) := Submodule.finrank_mono hP
  have h3 := Submodule.finrank_add_le_finrank_add_finrank U Q
  have h4 : finrank K U ≤ finrank K ↥(P ⊓ P.comap L) := Submodule.finrank_mono hU
  omega

/-- The submodule `Π i, M i` has the sum of the dimensions. -/
def piEquiv {ι : Type*} (M : ι → Submodule K V) :
    (Submodule.pi Set.univ M) ≃ₗ[K] (∀ i, M i) where
  toFun x i := ⟨x.1 i, x.2 i (Set.mem_univ i)⟩
  invFun y := ⟨fun i => (y i).1, fun i _ => (y i).2⟩
  map_add' x y := rfl
  map_smul' c x := rfl
  left_inv x := rfl
  right_inv y := rfl

lemma finrank_pi_const {ι : Type*} [Fintype ι] (M : Submodule K V) :
    finrank K (Submodule.pi Set.univ (fun _ : ι => M)) = Fintype.card ι * finrank K M := by
  rw [(piEquiv (fun _ : ι => M)).finrank_eq, Module.finrank_pi_fintype]
  simp [Finset.sum_const, Finset.card_univ]

end Abstract

/-! ## The engine as linear maps -/
section Engine
open Binary GF Matrix
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- `matAct f M` as a linear map. -/
noncomputable def matL (f : ℕ) (M : CMat α) : (Sky α f → ℂ) →ₗ[ℂ] (Sky α f → ℂ) :=
  Matrix.mulVecLin (digitProd fun _ : Fin f => M)

@[simp] lemma matL_apply (f : ℕ) (M : CMat α) (x : Sky α f → ℂ) : matL f M x = matAct f M x := rfl

lemma matL_mul (f : ℕ) (M M' : CMat α) : matL f (M * M') = matL f M ∘ₗ matL f M' := by
  apply LinearMap.ext; intro x; simp [matAct_mul]

lemma matL_one (f : ℕ) : matL f (1 : CMat α) = LinearMap.id := by
  apply LinearMap.ext; intro x; simp

/-- `roleAct f l M` as a linear map. -/
noncomputable def roleL (f : ℕ) (l : ρ) (M : CMat α) : Data α ρ f →ₗ[ℂ] Data α ρ f where
  toFun := roleAct f l M
  map_add' x y := by
    funext r; by_cases h : r = l <;> simp [roleAct, h, matAct, Matrix.mulVec_add]
  map_smul' c x := by
    funext r; by_cases h : r = l <;> simp [roleAct, h, matAct, Matrix.mulVec_smul]

/-- The free mixing letter as a linear map. -/
noncomputable def pointL (f : ℕ) (g : Matrix ρ ρ ℚ) : Data α ρ f →ₗ[ℂ] Data α ρ f where
  toFun := point f (actPoint g)
  map_add' x y := by
    funext r i; simp [point, actPoint, mul_add, Finset.sum_add_distrib]
  map_smul' c x := by
    funext r i; simp [point, actPoint, Finset.mul_sum, mul_left_comm]

/-- Each letter of the engine as a linear map. -/
noncomputable def letterL (f : ℕ) : GMove α ρ → Data α ρ f →ₗ[ℂ] Data α ρ f
  | .old (.localM g) => pointL f g
  | .old (.shift l z) => roleL f l (Binary.shift z)
  | .old (.block l _ z _) => roleL f l (blockMat z)
  | .perm l G => roleL f l (permMat G)
  | .phase l z c => roleL f l (phaseMat z c)

lemma letterL_apply (f : ℕ) (mv : GMove α ρ) (x : Data α ρ f) : letterL f mv x = mv.act f x := by
  rcases mv with (g | ⟨l, z⟩ | ⟨l, rk, z, ind⟩) | ⟨l, G⟩ | ⟨l, z, c⟩ <;> rfl

/-- A word as a linear map. -/
noncomputable def wordL (f : ℕ) : GWord α ρ → Data α ρ f →ₗ[ℂ] Data α ρ f
  | [] => LinearMap.id
  | mv :: w => wordL f w ∘ₗ letterL f mv

lemma wordL_apply (f : ℕ) (w : GWord α ρ) (x : Data α ρ f) : wordL f w x = gwalk f w x := by
  induction w generalizing x with
  | nil => rfl
  | cons mv w ih =>
    simp only [wordL, LinearMap.comp_apply, ih, letterL_apply, gwalk_cons]

/-- **Theorem R along a word.** -/
lemma pot_wordL (f : ℕ) (P A : Submodule ℂ (Data α ρ f)) (w : GWord α ρ) :
    pot P (A.map (wordL f w)) ≤ pot P A + (w.map fun mv => raise P (letterL f mv)).sum := by
  induction w generalizing A with
  | nil => simp [wordL, Submodule.map_id]
  | cons mv w ih =>
    simp only [wordL, List.map_cons, List.sum_cons]
    rw [Submodule.map_comp]
    have h1 := ih (A.map (letterL f mv))
    have h2 := pot_map_le P A (letterL f mv)
    omega

/-- `N` in every role. -/
def Npi (f : ℕ) (N : Submodule ℂ (Sky α f → ℂ)) : Submodule ℂ (Data α ρ f) :=
  Submodule.pi Set.univ (fun _ => N)

lemma mem_Npi {f : ℕ} {N : Submodule ℂ (Sky α f → ℂ)} {x : Data α ρ f} :
    x ∈ (Npi f N : Submodule ℂ (Data α ρ f)) ↔ ∀ r, x r ∈ N := by
  simp [Npi, Submodule.mem_pi]

variable {f : ℕ} (N : Submodule ℂ (Sky α f → ℂ))

lemma raise_pointL (g : Matrix ρ ρ ℚ) : raise (Npi f N) (pointL (α := α) f g) = 0 := by
  apply raise_eq_zero_of_mapsTo
  intro x hx
  rw [mem_Npi] at hx ⊢
  intro r
  have : (pointL (α := α) f g x) r = ∑ j, ((g r j : ℚ) : ℂ) • x j := by
    funext i; simp [pointL, point, actPoint, Finset.sum_apply]
  rw [this]
  exact Submodule.sum_mem _ fun j _ => Submodule.smul_mem _ _ (hx j)

lemma raise_roleL_keep (l : ρ) (M : CMat α) (hM : ∀ v ∈ N, matAct f M v ∈ N) :
    raise (Npi f N) (roleL f l M) = 0 := by
  apply raise_eq_zero_of_mapsTo
  intro x hx
  rw [mem_Npi] at hx ⊢
  intro r
  by_cases h : r = l
  · subst h; simp [roleL, roleAct, hM _ (hx r)]
  · simp [roleL, roleAct, h, hx r]

lemma raise_roleL_le (l : ρ) (M : CMat α) :
    raise (Npi f N) (roleL f l M) ≤ raise N (matL f M) := by
  -- the inclusion of the role `l`, read modulo `N`
  have hle : N ≤ LinearMap.ker ((Npi f N : Submodule ℂ (Data α ρ f)).mkQ ∘ₗ
      LinearMap.single ℂ (fun _ : ρ => Sky α f → ℂ) l) := by
    intro v hv
    simp only [LinearMap.mem_ker, LinearMap.comp_apply, Submodule.mkQ_apply,
      Submodule.Quotient.mk_eq_zero, mem_Npi]
    intro r
    by_cases h : r = l
    · subst h; simpa using hv
    · simp [Pi.single_apply, h, Submodule.zero_mem]
  set ι := N.liftQ _ hle
  have hsub : (Npi f N).map ((Npi f N : Submodule ℂ (Data α ρ f)).mkQ ∘ₗ roleL f l M) ≤
      (N.map (N.mkQ ∘ₗ matL f M)).map ι := by
    rintro _ ⟨x, hx, rfl⟩
    replace hx := mem_Npi.mp hx
    refine ⟨N.mkQ (matL f M (x l)), ⟨x l, hx l, rfl⟩, ?_⟩
    simp only [ι, Submodule.liftQ_apply, LinearMap.comp_apply, Submodule.mkQ_apply,
      LinearMap.coe_single]
    rw [Submodule.Quotient.eq, mem_Npi]
    intro r
    by_cases h : r = l
    · subst h; simp [roleL, roleAct, Submodule.zero_mem]
    · simp [roleL, roleAct, h, Pi.single_apply, hx r]
  calc raise (Npi f N) (roleL f l M) ≤ finrank ℂ ((N.map (N.mkQ ∘ₗ matL f M)).map ι) :=
        Submodule.finrank_mono hsub
    _ ≤ raise N (matL f M) := Submodule.finrank_map_le _ _

lemma roleL_mul (l : ρ) (M M' : CMat α) :
    roleL (ρ := ρ) f l (M * M') = roleL f l M ∘ₗ roleL f l M' := by
  apply LinearMap.ext; intro x
  simp only [roleL, LinearMap.coe_mk, AddHom.coe_mk, LinearMap.comp_apply, roleAct_mul]

lemma raise_roleL_prod (l : ρ) (cap : ℕ) :
    ∀ Ms : List (CMat α), (∀ M ∈ Ms, raise N (matL f M) ≤ cap) →
      raise (Npi f N) (roleL (ρ := ρ) f l Ms.prod) ≤ Ms.length * cap
  | [], _ => by
    simp only [List.prod_nil, List.length_nil, zero_mul, Nat.le_zero]
    apply raise_roleL_keep; intro v hv; simpa using hv
  | M :: Ms, h => by
    rw [List.prod_cons, roleL_mul]
    have h1 := raise_comp_le (Npi f N) (roleL (ρ := ρ) f l Ms.prod) (roleL f l M)
    have h2 := raise_roleL_prod l cap Ms (fun M' hM' => h M' (List.mem_cons_of_mem _ hM'))
    have h3 := raise_roleL_le N l M
    have h4 := h M List.mem_cons_self
    simp only [List.length_cons]
    nlinarith

/-- The hypotheses on `N`: kept by the free letters. -/
structure Kept (N : Submodule ℂ (Sky α f → ℂ)) : Prop where
  shift : ∀ z : Space α, ∀ v ∈ N, matAct f (Binary.shift z) v ∈ N
  perm : ∀ G : APerm α, ∀ v ∈ N, matAct f (permMat G) v ∈ N
  phase : ∀ (z : Space α) (c : F), ∀ v ∈ N, matAct f (phaseMat z c) v ∈ N

lemma raise_letterL_le {N : Submodule ℂ (Sky α f → ℂ)} (hN : Kept N) (cap : ℕ)
    (hcap : ∀ z : Space α, z ≠ 0 → raise N (matL f (dir z)) ≤ cap) (mv : GMove α ρ) :
    raise (Npi f N) (letterL f mv) ≤ cap * mv.shadow.cost id := by
  rcases mv with (g | ⟨l, z⟩ | ⟨l, rk, z, ind⟩) | ⟨l, G⟩ | ⟨l, z, c⟩
  · simp [letterL, raise_pointL]
  · simp [letterL, raise_roleL_keep N l _ (hN.shift z)]
  · simp only [letterL, GMove.shadow, BMove.cost, BMove.rk?, id]
    have := raise_roleL_prod (ρ := ρ) N l cap (List.ofFn fun i => dir (z i)) (by
      intro M hM
      obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hM
      exact hcap _ (ind.ne_zero i))
    simp only [List.length_ofFn] at this
    unfold blockMat; linarith
  · simp [letterL, raise_roleL_keep N l _ (hN.perm G)]
  · simp [letterL, raise_roleL_keep N l _ (hN.phase z c)]

lemma sum_raise_le {N : Submodule ℂ (Sky α f → ℂ)} (hN : Kept N) (cap : ℕ)
    (hcap : ∀ z : Space α, z ≠ 0 → raise N (matL f (dir z)) ≤ cap) (w : GWord α ρ) :
    (w.map fun mv => raise (Npi f N) (letterL f mv)).sum ≤ cap * w.cost id := by
  induction w with
  | nil => simp
  | cons mv w ih =>
    simp only [List.map_cons, List.sum_cons, GWord.cost_cons, mul_add]
    have := raise_letterL_le hN cap hcap mv
    omega

/-- The inputs: `N` in the live roles, `0` in the scratch roles. -/
def inputs {σ : Type*} (e : σ → ρ) (f : ℕ) (N : Submodule ℂ (Sky α f → ℂ)) :
    Submodule ℂ (Data α ρ f) where
  carrier := {v | (∀ s, (∀ l, e l ≠ s) → v s = 0) ∧ ∀ r, v r ∈ N}
  add_mem' := by
    rintro a b ⟨ha1, ha2⟩ ⟨hb1, hb2⟩
    exact ⟨fun s hs => by simp [ha1 s hs, hb1 s hs], fun r => N.add_mem (ha2 r) (hb2 r)⟩
  zero_mem' := ⟨fun _ _ => rfl, fun _ => N.zero_mem⟩
  smul_mem' := by
    rintro c a ⟨ha1, ha2⟩
    exact ⟨fun s hs => by simp [ha1 s hs], fun r => N.smul_mem c (ha2 r)⟩

/-- Extension by zero from the live roles. -/
noncomputable def extL {σ : Type*} [Fintype σ] (e : σ → ρ) (f : ℕ) :
    (σ → Sky α f → ℂ) →ₗ[ℂ] Data α ρ f where
  toFun x r := ∑ l, if e l = r then x l else 0
  map_add' x y := by
    funext r; simp only [Pi.add_apply]; rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl; intro l _; split_ifs <;> simp
  map_smul' c x := by
    funext r; simp only [Pi.smul_apply, RingHom.id_apply]; rw [Finset.smul_sum]
    apply Finset.sum_congr rfl; intro l _; split_ifs <;> simp

lemma extL_apply_e {σ : Type*} [Fintype σ] {e : σ → ρ} (he : Function.Injective e)
    (x : σ → Sky α f → ℂ) (l : σ) : extL e f x (e l) = x l := by
  simp only [extL, LinearMap.coe_mk, AddHom.coe_mk]
  rw [Finset.sum_eq_single l]
  · simp
  · intro b _ hb; simp [he.ne hb]
  · simp

lemma extL_apply_off {σ : Type*} [Fintype σ] {e : σ → ρ}
    (x : σ → Sky α f → ℂ) (s : ρ) (hs : ∀ l, e l ≠ s) : extL e f x s = 0 := by
  simp [extL, hs]

/-- **The end count.** -/
lemma card_mul_raise_le_pot {σ : Type*} [Fintype σ] {e : σ → ρ} (he : Function.Injective e)
    (w : GWord α ρ) (hw : GLiveKernel e w) :
    Fintype.card σ * raise N (matL f (kernel α)) ≤
      pot (Npi f N) ((inputs e f N).map (wordL f w)) := by
  have hle : (Npi f N : Submodule ℂ (Data α ρ f)) ≤ LinearMap.ker
      (LinearMap.pi fun l : σ => N.mkQ ∘ₗ LinearMap.proj (e l)) := by
    intro x hx
    rw [mem_Npi] at hx
    simp only [LinearMap.mem_ker]
    funext l
    simp [(Submodule.Quotient.mk_eq_zero N).mpr (hx (e l))]
  have hsub : Submodule.pi Set.univ (fun _ : σ => N.map (N.mkQ ∘ₗ matL f (kernel α))) ≤
      (((inputs e f N).map (wordL f w)).map (Npi f N).mkQ).map
        ((Npi f N : Submodule ℂ (Data α ρ f)).liftQ _ hle) := by
    intro y hy
    rw [Submodule.mem_pi] at hy
    choose x hxN hxy using fun l => (hy l (Set.mem_univ l))
    refine ⟨(Npi f N).mkQ (wordL f w (extL e f x)), ⟨_, ⟨_, ⟨?_, ?_⟩, rfl⟩, rfl⟩, ?_⟩
    · intro s hs; exact extL_apply_off x s hs
    · intro r
      by_cases h : ∃ l, e l = r
      · obtain ⟨l, rfl⟩ := h; rw [extL_apply_e he]; exact hxN l
      · push Not at h; rw [extL_apply_off x r h]; exact N.zero_mem
    · funext l
      simp only [Submodule.mkQ_apply, Submodule.liftQ_apply, LinearMap.pi_apply,
        LinearMap.comp_apply, LinearMap.proj_apply, wordL_apply]
      rw [hw f (extL e f x) (fun s hs => extL_apply_off x s hs) l, extL_apply_e he, ← hxy l]
      rfl
  have h1 := Submodule.finrank_mono hsub
  have h2 := Submodule.finrank_map_le ((Npi f N : Submodule ℂ (Data α ρ f)).liftQ _ hle)
    (((inputs e f N).map (wordL f w)).map (Npi f N).mkQ)
  rw [finrank_pi_const] at h1
  unfold pot; unfold raise; omega

lemma pot_inputs (σ : Type*) (e : σ → ρ) : pot (Npi f N) (inputs e f N) = 0 :=
  pot_eq_zero_of_le fun _ hv => mem_Npi.mpr hv.2

/-- **Theorem R for the engine.**  If `N` is kept by every free letter and every direction raises
by at most `cap`, then every word that delivers the kernel to `card σ` distinct live roles (any
number of scratch roles, any free letters, blocks of any rank) has total block rank at least
`card σ * raise N kernel / cap`. -/
theorem engine_lower_bound {σ : Type*} [Fintype σ] {e : σ → ρ} (he : Function.Injective e)
    (hN : Kept N) (cap : ℕ) (hcap : ∀ z : Space α, z ≠ 0 → raise N (matL f (dir z)) ≤ cap)
    (w : GWord α ρ) (hw : GLiveKernel e w) :
    Fintype.card σ * raise N (matL f (kernel α)) ≤ cap * w.cost id := by
  have h1 := card_mul_raise_le_pot N he w hw
  have h2 := pot_wordL f (Npi f N) (inputs e f N) w
  have h3 := pot_inputs N σ e
  have h4 := sum_raise_le hN cap hcap w
  omega

end Engine
end PowerSaving.PfLean
end OAI

#print axioms OAI.PowerSaving.PfLean.engine_lower_bound
