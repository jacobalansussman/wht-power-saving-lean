import Work.Mixed.Prog.SweepUp

/-!
# Mixed engine, R3 program part 16: the sweep of a leaf class (agent e-r3prog)

`Eng.sweepLeaf`: every letter of the word on a leaf class; a block letter is ONE fresh top call
(`exit_letter`).  The charge is exactly the list of exits `kidsLeaf`.
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
variable {ι : Type U} {m : Bool} {cl : Paint} {b : Bank (portM cl) ι}
variable {n μ : ι → ℕ} {τ : ι → MDesc → ℕ}

lemma cU_leaf (n μ p : ℕ) : E.cUc n 0 μ p = μ := by
  simp [cUc, MXR.cU]

/-- Sweep cost of a word on a leaf class. -/
def lfCost (n μ : ℕ) : GWord α ρ → ℕ
  | [] => 0
  | mv :: ws => (match mv.rk? with
      | some _ => E.NUc n 0 μ * 2^n + 2 * n + 2
      | none => E.NUc n 0 μ * 2^n) + lfCost n μ ws

/-- **The sweep of a leaf class.** -/
theorem sweepLeaf (ws : GWord α ρ) :
    ∀ {s : Ty} {x : ι → T s} {V : ∀ i, E.CState (n i) 0 (μ i)},
    CIn m b.B s x n (fun _ => 0) μ →
    CCtx m b.B s x E.X.n n (fun i => E.NUc (n i) 0 (μ i)) (fun i p => E.cUc (n i) 0 (μ i) p)
      (fun i p => E.iUc (n i) 0 (μ i) p) →
    Has m b.B s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) 0 (μ i) (V i)) →
    ws.Proper E.X.n → (∀ i, E.X.n * μ i ≤ n i) →
    (∀ i c, c ∈ E.kidsLeaf (n i) (μ i) (V i) ws → (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) →
    Able m b s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) 0 (μ i) (E.runV (V i) ws))
      (fun i => E.lfCost (n i) (μ i) ws) (fun i => ((E.kidsLeaf (n i) (μ i) (V i) ws).map (τ i)).sum) := by
  induction ws with
  | nil =>
    intro s x V cin cc ct hP hfit hk
    refine Able.mono ((Able.of_can (b:=b) (ct.mono (W':=fun _ => 0) (fun _ => le_rfl))).cong (fun _ => rfl)
      (fun i => by rw [runV_nil])) (fun i => Nat.zero_le _) (fun i => Nat.zero_le _)
  | cons mv ws ih =>
    intro s x V cin cc ct hP hfit hk
    have cont : ∀ (V' : ∀ i, E.CState (n i) 0 (μ i)) (W1 Z1 : ι → ℕ),
        Able m b s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) 0 (μ i) (V' i)) W1 Z1 →
        (∀ i c, c ∈ E.kidsLeaf (n i) (μ i) (V' i) ws → (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) →
        Able m b s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) 0 (μ i) (E.runV (V' i) ws))
          (fun i => W1 i + E.lfCost (n i) (μ i) ws)
          (fun i => Z1 i + ((E.kidsLeaf (n i) (μ i) (V' i) ws).map (τ i)).sum) := by
      intro V' W1 Z1 H1 hk'
      let st (i : ι) : T (a (a (a (c cl)))) := E.ctabE (n i) 0 (μ i) (V' i)
      have H2 := ih (s:=p s (a (a (a (c cl))))) (x:=fun i => (x i, st i)) (V:=V') (cin.ext (t:=a (a (a (c cl)))) st)
        (cc.ext (t:=a (a (a (c cl)))) st)
        Can.second hP.tail hfit hk'
      exact H1.save H2
    have free : mv.rk? = none →
        (∀ i, E.kidsLeaf (n i) (μ i) (V i) (mv :: ws) = E.kidsLeaf (n i) (μ i) (E.stepV (V i) mv) ws) →
        Able m b s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) 0 (μ i) (E.runV (V i) (mv :: ws)))
          (fun i => E.lfCost (n i) (μ i) (mv :: ws)) (fun i => ((E.kidsLeaf (n i) (μ i) (V i) (mv :: ws)).map (τ i)).sum) := by
      intro hfree hkids
      have H1 : Able m b s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) 0 (μ i) (E.stepV (V i) mv))
          (fun i => E.NUc (n i) 0 (μ i) * 2^(n i)) (fun _ => 0) :=
        Able.of_can (b:=b) (cc.free_letter ct rfl mv hfree)
      have H := cont (fun i => E.stepV (V i) mv) _ _ H1 (fun i c hc => hk i c (by rw [hkids i]; exact hc))
      refine (H.cong (fun _ => rfl) (fun i => ?_)).mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
      · rw [runV_cons]
      · simp only [lfCost, hfree]
      · rw [hkids i, Nat.zero_add]
    cases mv with
    | old mv =>
      cases mv with
      | localM g => exact free rfl (fun i => rfl)
      | shift l z => exact free rfl (fun i => rfl)
      | block l r z ind =>
        have hpr : 1 ≤ r ∧ r < E.X.n := hP.head
        have hr : r ≤ E.X.n := hpr.2.le
        have hkn : ∀ i, r * μ i ≤ n i := fun i => le_trans (Nat.mul_le_mul_right _ hr) (hfit i)
        have hkids : ∀ i, E.kidsLeaf (n i) (μ i) (V i) (.old (.block l r z ind) :: ws) =
            E.exitC (n i) (μ i) (V i) l r z ind hr :: E.kidsLeaf (n i) (μ i) (E.stepV (V i) (.old (.block l r z ind))) ws :=
          fun i => by simp only [kidsLeaf, dif_pos hr, List.singleton_append]
        -- the powers for the exit
        have hke : Knows m b.B s x (fun i => r * μ i) := (Can.wc r).times cin.μK (Small.const _ _) cin.μS
        have hkes : Small b.B (fun i => r * μ i) := (Small.const _ r).mul cin.μS
        have hp1 : Can m b.B s w x (fun i => 2^(r * μ i)) (fun i => n i + 1) :=
          (hke.mono (fun _ => Nat.zero_le _)).wTwo hkes (cin.big.mono (fun i => Nat.pow_le_pow_right (by omega) (hkn i)))
            (Dom.of_le (fun i => by have := hkn i; omega))
        have hnk : Knows m b.B s x (fun i => n i - r * μ i) := cin.nK.minus hke cin.nS
        have hp2 : Can m b.B s w x (fun i => 2^(n i - r * μ i)) (fun i => n i + 1) :=
          (hnk.mono (fun _ => Nat.zero_le _)).wTwo (cin.nS.sub) (cin.big.mono (fun i => Nat.pow_le_pow_right (by omega)
            (Nat.sub_le _ _))) (Dom.of_le (fun i => by omega))
        let PW (i : ι) : T (p w (p w w)) := (2^(r * μ i), 2^(n i - r * μ i), r * μ i)
        have HPW : Can m b.B s (p w (p w w)) x PW (fun i => n i + 1 + (n i + 1 + 0)) := hp1.fork (hp2.fork hke)
        let x1 (i : ι) : T (p s (p w (p w w))) := (x i, PW i)
        have cc1 := cc.ext PW
        have ct1 : Has m b.B (p s (p w (p w w))) (a (a (a (c cl)))) x1 (fun i => E.ctabE (n i) 0 (μ i) (V i)) :=
          Can.first.then_do ct
        let S := bS z ind
        let X' := bX (α:=α) r
        have hX' : X'.n + r = E.X.n := E.bX_n hr
        have hk1 : ∀ i, (b.use i (mkIn cl (r * μ i) 0 0 0
            (jobsTape (r * μ i) (E.NUc (n i) 0 (μ i) * 2^(n i - r * μ i))
              (toVec (n i) (r * μ i) (hkn i) (exJobs E.X S X' hX' l (cc1.fit i) (V i))))
            (jobsTape (r * μ i) 0 (fun _ _ => 0)))).OK
          (jobsTape (r * μ i) (E.NUc (n i) 0 (μ i) * 2^(n i - r * μ i))
            (fun t => ripple (r * μ i) (toVec (n i) (r * μ i) (hkn i) (exJobs E.X S X' hX' l (cc1.fit i) (V i)) t)),
            jobsTape (r * μ i) 0 (fun _ _ => 0), (Layout.std 0).tape (fun _ => (0 : ℕ)))
          (τ i (E.exitC (n i) (μ i) (V i) l r z ind hr)) (b.ext i) := by
          intro i
          have := hk i (E.exitC (n i) (μ i) (V i) l r z ind hr) (by rw [hkids i]; exact List.mem_cons_self)
          have e : E.exitC (n i) (μ i) (V i) l r z ind hr = MDesc.top (r * μ i) (E.NUc (n i) 0 (μ i) * 2^(n i - r * μ i))
              (toVec (n i) (r * μ i) (hkn i) (exJobs E.X S X' hX' l (cc1.fit i) (V i))) := by
            unfold exitC; rw [dif_pos (hkn i)]
          rw [e] at this ⊢
          exact this
        have HE := cc1.exit_letter ct1 S X' hX' l μ (fun i p => E.cU_leaf _ _ _) hkn
          (Can.fst Can.second) (Can.fst (Can.snd Can.second)) (Can.snd (Can.snd Can.second)) hk1
        have H1 : Able m b s (a (a (a (c cl)))) x
            (fun i => E.ctabE (n i) 0 (μ i) (E.stepV (V i) (.old (.block l r z ind))))
            (fun i => n i + 1 + (n i + 1 + 0) + E.NUc (n i) 0 (μ i) * 2^(n i))
            (fun i => 0 + τ i (E.exitC (n i) (μ i) (V i) l r z ind hr)) := by
          have := (Able.of_can (b:=b) HPW).save (Able.fst HE)
          refine this.cong (fun _ => rfl) (fun i => ?_)
          rw [E.stepV_block]
        have H := cont (fun i => E.stepV (V i) (.old (.block l r z ind))) _ _ H1
          (fun i c hc => hk i c (by rw [hkids i]; exact List.mem_cons_of_mem _ hc))
        refine (H.cong (fun _ => rfl) (fun i => ?_)).mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
        · rw [runV_cons]
        · show _ = (E.NUc (n i) 0 (μ i) * 2^(n i) + 2 * n i + 2) + E.lfCost (n i) (μ i) ws
          ring
        · rw [hkids i, List.map_cons, List.sum_cons, Nat.zero_add]
    | perm l G => exact free rfl (fun i => rfl)
    | phase l z c' => exact free rfl (fun i => rfl)

end Eng

end
end PowerSaving.RAM
end OAI
