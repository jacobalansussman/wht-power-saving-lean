import Work.Mixed.Prog.ClassBody

/-!
# Mixed engine, R3 program part 18: valid descriptors (agent e-r3prog)

The handler answers VALID descriptors: every top job, and every class of the lockstep tree
(`1 ≤ μ`, `mx * μ * kmin^h ≤ n`: rl-prove's region invariant; its subtree's consumption fits in the supply).

* `Eng.Valid`        the predicate;
* `Eng.kidsC_valid`  the calls of a valid class are valid;
* `Eng.class_body_valid`  the class body for any family of valid class descriptors.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

universe U

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)

/-- Valid descriptors. -/
def Valid (kmin : ℕ) : MDesc → Prop
  | .top _ _ _ => True
  | .cls n h μ o BS _ _ => 1 ≤ μ ∧ E.mx * μ * kmin ^ h ≤ n ∧ o + E.cons n h μ ≤ BS

/-- The smallest target rank (used by `Valid`). -/
structure KMin (kmin : ℕ) : Prop where
  two : 2 ≤ kmin
  att : ∃ l, E.tr l = kmin

lemma rk_lt_kmin {kmin : ℕ} (hk : E.KMin kmin) {mv : GMove α ρ} (hmv : mv ∈ E.wd) {r : ℕ}
    (hr : mv.rk? = some r) : r < kmin := by
  obtain ⟨l, hl⟩ := hk.att
  rw [← hl]; exact E.sep mv hmv r hr l

lemma accW_cons_len (n h μ : ℕ) (mv : GMove α ρ) (ws : GWord α ρ) :
    (E.accW n h μ ws).length ≤ (E.accW n h μ (mv :: ws)).length := by
  rw [accW_cons]; simp

/-- The children of a class along a word are valid. -/
lemma kidsUp_valid {kmin : ℕ} (hk : E.KMin kmin) (n h μ BS : ℕ) (sup : ℕ → Bits n → ℂ) (hμ : 1 ≤ μ)
    (hreg : E.mx * μ * kmin ^ (h+1) ≤ n) :
    ∀ (ws : GWord α ρ), (∀ mv ∈ ws, mv ∈ E.wd) → ∀ (V : E.CState n (h+1) μ) (o' : ℕ),
      o' + (E.accW n h μ ws).length ≤ BS → ∀ c ∈ E.kidsUp n h μ BS sup V o' ws, E.Valid kmin c := by
  intro ws
  induction ws with
  | nil => intro _ V o' _ c hc; simp [kidsUp] at hc
  | cons mv ws ih =>
    intro hsub V o' hfit c hc
    have hsub' : ∀ mv' ∈ ws, mv' ∈ E.wd := fun mv' hm => hsub mv' (List.mem_cons_of_mem _ hm)
    cases mv with
    | old mv =>
      cases mv with
      | localM g => exact ih hsub' _ o' (le_trans (by have := E.accW_cons_len n h μ (.old (.localM g)) ws; omega) hfit) c hc
      | shift l z => exact ih hsub' _ o' (le_trans (by have := E.accW_cons_len n h μ (.old (.shift l z)) ws; omega) hfit) c hc
      | block l r z ind =>
        have hrk : r < kmin := E.rk_lt_kmin hk (hsub _ List.mem_cons_self) rfl
        have hacc : E.accW n h μ (.old (.block l r z ind) :: ws) = E.accR n h μ r ++ E.accW n h μ ws := by
          rw [accW_cons]; rfl
        rw [hacc, List.length_append] at hfit
        simp only [kidsUp] at hc
        rcases List.mem_append.1 hc with hc | hc
        · split_ifs at hc with hr
          · rw [List.mem_singleton] at hc
            subst hc
            refine ⟨?_, ?_, ?_⟩
            · have h1 : (GMove.old (.block l r z ind) : GMove α ρ).Proper E.X.n := E.hP _ (hsub _ List.mem_cons_self)
              obtain ⟨h1, _⟩ := h1
              exact Nat.one_le_iff_ne_zero.2 (Nat.mul_ne_zero (by omega) (by omega))
            · calc E.mx * (μ * r) * kmin ^ h ≤ E.mx * μ * kmin * kmin ^ h := by
                    have := Nat.mul_le_mul_left (E.mx * μ) hrk.le
                    calc E.mx * (μ * r) * kmin ^ h = E.mx * μ * r * kmin ^ h := by ring
                      _ ≤ E.mx * μ * kmin * kmin ^ h := Nat.mul_le_mul_right _ this
                _ = E.mx * μ * kmin ^ (h + 1) := by ring
                _ ≤ n := hreg
            · have : (E.accR n h μ r).length = (E.ownSrc n h μ r).length + E.cons n h (μ * r) := by
                simp [accR, cons]
              omega
          · simp at hc
        · exact ih hsub' _ _ (by simp only [accR, cons, List.length_append] at hfit ⊢; omega) c hc
    | perm l G => exact ih hsub' _ o' (le_trans (by have := E.accW_cons_len n h μ (.perm l G) ws; omega) hfit) c hc
    | phase l z c' => exact ih hsub' _ o' (le_trans (by have := E.accW_cons_len n h μ (.phase l z c') ws; omega) hfit) c hc

/-- The exits of a leaf class are valid (top jobs). -/
lemma kidsLeaf_valid {kmin : ℕ} (n μ : ℕ) : ∀ (ws : GWord α ρ) (V : E.CState n 0 μ),
    ∀ c ∈ E.kidsLeaf n μ V ws, E.Valid kmin c := by
  intro ws
  induction ws with
  | nil => intro V c hc; simp [kidsLeaf] at hc
  | cons mv ws ih =>
    intro V c hc
    cases mv with
    | old mv =>
      cases mv with
      | localM g => exact ih _ c hc
      | shift l z => exact ih _ c hc
      | block l r z ind =>
        simp only [kidsLeaf] at hc
        rcases List.mem_append.1 hc with hc | hc
        · split_ifs at hc with hr
          · rw [List.mem_singleton] at hc; subst hc
            unfold exitC; split_ifs <;> trivial
          · simp at hc
        · exact ih _ c hc
    | perm l G => exact ih _ c hc
    | phase l z c' => exact ih _ c hc

/-- **The calls of a valid class are valid.** -/
theorem kidsC_valid {kmin : ℕ} (hk : E.KMin kmin) (d : MDesc) (hd : E.Valid kmin d) :
    ∀ c ∈ E.kidsC d, E.Valid kmin c := by
  cases d with
  | top n B v => intro c hc; simp [kidsC] at hc
  | cls n h μ o BS J sup =>
    obtain ⟨hμ, hreg, hfit⟩ := hd
    cases h with
    | zero => exact E.kidsLeaf_valid n μ _ _
    | succ h =>
      exact E.kidsUp_valid hk n h μ BS sup hμ hreg E.wd (fun _ hm => hm) _ o (by
        have e : (E.accW n h μ E.wd).length = E.cons n (h+1) μ := rfl
        omega)

lemma valid_facts {kmin : ℕ} (hk : E.KMin kmin) {n h μ o BS : ℕ} {J sup : ℕ → Bits n → ℂ}
    (hd : E.Valid kmin (MDesc.cls n h μ o BS J sup)) :
    (h = 0 → E.X.n * μ ≤ n) ∧ h ≤ n ∧ μ ≤ n ∧ o + E.cons n h μ ≤ BS := by
  obtain ⟨hμ, hreg, hfit⟩ := hd
  have hk1 : 1 ≤ kmin ^ h := Nat.one_le_pow _ _ (by have := hk.two; omega)
  have hm : 1 ≤ E.mx := le_trans (by norm_num) E.three
  refine ⟨fun h0 => by subst h0; simpa using hreg, ?_, ?_, hfit⟩
  · have h1 : h < 2 ^ h := Nat.lt_two_pow_self
    have h2 : 2 ^ h ≤ kmin ^ h := Nat.pow_le_pow_left hk.two _
    have h3 : kmin ^ h ≤ E.mx * μ * kmin ^ h := Nat.le_mul_of_pos_left _ (Nat.mul_pos hm hμ)
    omega
  · have : μ ≤ E.mx * μ * kmin ^ h := by
      calc μ ≤ E.mx * μ := Nat.le_mul_of_pos_left _ hm
        _ ≤ E.mx * μ * kmin ^ h := Nat.le_mul_of_pos_right _ hk1
    omega

/-- **The class body for valid descriptors.** -/
theorem class_body_valid {kmin : ℕ} (hk : E.KMin kmin) {ι : Type U} {m : Bool} {cl : Paint}
    {b : Bank (portM cl) ι} {n h μ o BS : ι → ℕ} {J sup : ∀ i, ℕ → Bits (n i) → ℂ} {τ : ι → MDesc → ℕ}
    (hB2 : Small b.B (fun i => 2^(n i))) (hBS : Small b.B BS)
    (hd : ∀ i, E.Valid kmin (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i)))
    (hkc : ∀ i c, c ∈ E.kidsC (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i)) →
      (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) :
    Able m b (portIn cl) (portOut cl)
      (fun i => E.enc cl (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.out cl (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.clsLoc (n i) (h i) (μ i))
      (fun i => ((E.kidsC (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i))).map (τ i)).sum) :=
  E.class_body hB2 hBS (fun i => (E.valid_facts hk (hd i)).1) (fun i => (E.valid_facts hk (hd i)).2.1)
    (fun i => (E.valid_facts hk (hd i)).2.2.1) (fun i => (E.valid_facts hk (hd i)).2.2.2) hkc

end Eng

end
end PowerSaving.RAM
end OAI
