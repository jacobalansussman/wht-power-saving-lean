import Work.Mixed.Prog.Sweep

/-!
# Mixed engine, R3 program part 15: the sweep of a class of height `h+1` over a word (agent e-r3prog)

`Eng.sweepUp`: by induction on the word, every letter in order: a free letter is one pass, a block
letter is ONE child call (`block_step`).  The charge is exactly the list of child calls `kidsUp`.
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
variable {n h μ o BS : ι → ℕ} {sup : ∀ i, ℕ → Bits (n i) → ℂ} {τ : ι → MDesc → ℕ}

/-- **The sweep of a class of height `h+1`.** -/
theorem sweepUp (ws : GWord α ρ) :
    ∀ {s : Ty} {x : ι → T s} {V : ∀ i, E.CState (n i) (h i + 1) (μ i)} {A : ι → List ℕ},
    CIn m b.B s x n (fun i => h i + 1) μ →
    CCtx m b.B s x E.X.n n (fun i => E.NUc (n i) (h i + 1) (μ i)) (fun i p => E.cUc (n i) (h i + 1) (μ i) p)
      (fun i p => E.iUc (n i) (h i + 1) (μ i) p) →
    Knows m b.B s x (fun i => E.W^(h i + 1)) →
    Has m b.B s (a w) x (fun i => (Layout.std (E.W^(h i + 1))).tape (fun t => E.PCc (n i) (h i + 1) (μ i) t)) →
    Has m b.B s (a w) x (fun i => (Layout.std (E.NUc (n i) (h i + 1) (μ i))).tape
      (fun p => E.Uc (n i) (h i + 1) (μ i) p)) →
    Has m b.B s (jobsT cl) x (fun i => jobsTape (n i) (BS i) (sup i)) → Small b.B BS →
    Has m b.B s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) (h i + 1) (μ i) (V i)) →
    Has m b.B s (jobsT cl) x (fun i => accT (n i) (o i) (sup i) (A i)) →
    Has m b.B s (a w) x (fun i => szT (A i)) →
    Knows m b.B s x (fun i => o i + (A i).length) →
    ws.Proper E.X.n → (∀ mv ∈ ws, ∀ r, mv.rk? = some r → ∀ l, r < E.tr l) →
    (∀ i, o i + (A i).length + (E.accW (n i) (h i) (μ i) ws).length ≤ BS i) →
    (∀ i c, c ∈ E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) ws →
      (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) →
    Able m b s (StT cl) x
      (fun i => (E.ctabE (n i) (h i + 1) (μ i) (E.runV (V i) ws),
        accT (n i) (o i) (sup i) (A i ++ E.accW (n i) (h i) (μ i) ws),
        szT (A i ++ E.accW (n i) (h i) (μ i) ws), o i + (A i ++ E.accW (n i) (h i) (μ i) ws).length))
      (fun i => E.swCost (n i) (h i) (μ i) (A i).length ws)
      (fun i => ((E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) ws).map (τ i)).sum) := by
  induction ws with
  | nil =>
    intro s x V A cin cc hWp hPCp hUp tsup hBSs ct tacc tsz horun hP hsep hfit hk
    have H := ct.fork (tacc.fork (tsz.fork horun))
    have H2 := Able.of_can (b:=b) (H.mono (W':=fun _ => 0) (fun _ => by simp))
    refine Able.mono (H2.cong (fun _ => rfl) (fun i => ?_)) (fun i => Nat.zero_le _) (fun i => Nat.zero_le _)
    simp [runV_nil, accW, blockRanks]
  | cons mv ws ih =>
    intro s x V A cin cc hWp hPCp hUp tsup hBSs ct tacc tsz horun hP hsep hfit hk
    -- the continuation after one letter: the facts in the extended environment
    have cont : ∀ (V' : ∀ i, E.CState (n i) (h i + 1) (μ i)) (A' : ι → List ℕ) (W1 Z1 : ι → ℕ),
        Able m b s (StT cl) x (fun i => (E.ctabE (n i) (h i + 1) (μ i) (V' i), accT (n i) (o i) (sup i) (A' i),
          szT (A' i), o i + (A' i).length)) W1 Z1 →
        (∀ i, o i + (A' i).length + (E.accW (n i) (h i) (μ i) ws).length ≤ BS i) →
        (∀ i c, c ∈ E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (V' i) (o i + (A' i).length) ws →
          (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) →
        Able m b s (StT cl) x
          (fun i => (E.ctabE (n i) (h i + 1) (μ i) (E.runV (V' i) ws),
            accT (n i) (o i) (sup i) (A' i ++ E.accW (n i) (h i) (μ i) ws),
            szT (A' i ++ E.accW (n i) (h i) (μ i) ws), o i + (A' i ++ E.accW (n i) (h i) (μ i) ws).length))
          (fun i => W1 i + E.swCost (n i) (h i) (μ i) (A' i).length ws)
          (fun i => Z1 i + ((E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (V' i) (o i + (A' i).length) ws).map (τ i)).sum) := by
      intro V' A' W1 Z1 H1 hfit' hk'
      let st (i : ι) : T (StT cl) := (E.ctabE (n i) (h i + 1) (μ i) (V' i), accT (n i) (o i) (sup i) (A' i),
          szT (A' i), o i + (A' i).length)
      have F : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z → Has m b.B (p s (StT cl)) u (fun i => (x i, st i)) z :=
        fun g => Can.first.then_do g
      have hst : Has m b.B (p s (StT cl)) (StT cl) (fun i => (x i, st i)) st := Can.second
      have H2 := ih (s:=p s (StT cl)) (x:=fun i => (x i, st i)) (V:=V') (A:=A') (cin.ext st) (cc.ext st) (F hWp) (F hPCp) (F hUp)
        (F tsup) hBSs (Can.fst hst) (Can.fst (Can.snd hst)) (Can.fst (Can.snd (Can.snd hst)))
        (Can.snd (Can.snd (Can.snd hst))) hP.tail (fun mv' hm => hsep mv' (List.mem_cons_of_mem _ hm)) hfit' hk'
      exact H1.save H2
    -- a free letter
    have free : mv.rk? = none →
        (∀ i, E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) (mv :: ws) =
          E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (E.stepV (V i) mv) (o i + (A i).length) ws) →
        Able m b s (StT cl) x
          (fun i => (E.ctabE (n i) (h i + 1) (μ i) (E.runV (V i) (mv :: ws)),
            accT (n i) (o i) (sup i) (A i ++ E.accW (n i) (h i) (μ i) (mv :: ws)),
            szT (A i ++ E.accW (n i) (h i) (μ i) (mv :: ws)), o i + (A i ++ E.accW (n i) (h i) (μ i) (mv :: ws)).length))
          (fun i => E.swCost (n i) (h i) (μ i) (A i).length (mv :: ws))
          (fun i => ((E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) (mv :: ws)).map (τ i)).sum) := by
      intro hfree hkids
      have hacc : ∀ i, E.accW (n i) (h i) (μ i) (mv :: ws) = E.accW (n i) (h i) (μ i) ws := fun i => by
        rw [accW_cons, hfree]; rfl
      have H1 : Able m b s (StT cl) x (fun i => (E.ctabE (n i) (h i + 1) (μ i) (E.stepV (V i) mv),
          accT (n i) (o i) (sup i) (A i), szT (A i), o i + (A i).length))
          (fun i => E.NUc (n i) (h i + 1) (μ i) * 2^(n i)) (fun _ => 0) :=
        Able.of_can (b:=b) (((cc.free_letter ct rfl mv hfree).fork (tacc.fork (tsz.fork horun))).mono
          (fun i => by simp))
      have H := cont (fun i => E.stepV (V i) mv) A _ _ H1 (fun i => by rw [← hacc i]; exact hfit i)
        (fun i c hc => hk i c (by rw [hkids i]; exact hc))
      refine (H.cong (fun _ => rfl) (fun i => ?_)).mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
      · rw [runV_cons, hacc i]
      · simp only [swCost, hfree, Nat.add_zero]
      · rw [hkids i, Nat.zero_add]
    cases mv with
    | old mv =>
      cases mv with
      | localM g => exact free rfl (fun i => rfl)
      | shift l z => exact free rfl (fun i => rfl)
      | block l r z ind =>
        have hpr : 1 ≤ r ∧ r < E.X.n := hP.head
        have hr : r ≤ E.X.n := hpr.2.le
        have hsep' : ∀ l', r < E.tr l' := hsep _ List.mem_cons_self r rfl
        have hacc : ∀ i, E.accW (n i) (h i) (μ i) (.old (.block l r z ind) :: ws) =
            E.accR (n i) (h i) (μ i) r ++ E.accW (n i) (h i) (μ i) ws := fun i => by
          rw [accW_cons]; rfl
        have hlen : ∀ i, o i + (A i).length + (E.ownSrc (n i) (h i) (μ i) r).length + E.cons (n i) (h i) (μ i * r) =
            o i + (A i ++ E.accR (n i) (h i) (μ i) r).length := fun i => by
          simp [accR, cons]; omega
        have hkids : ∀ i, E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length)
            (.old (.block l r z ind) :: ws) =
            E.childC (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) l r z ind hr ::
              E.kidsUp (n i) (h i) (μ i) (BS i) (sup i) (E.stepV (V i) (.old (.block l r z ind)))
                (o i + (A i ++ E.accR (n i) (h i) (μ i) r).length) ws := fun i => by
          simp only [kidsUp, dif_pos hr, List.singleton_append]
          rw [hlen i]
        have hfit1 : ∀ i, o i + (A i).length + (E.accR (n i) (h i) (μ i) r).length ≤ BS i := fun i => by
          have := hfit i; rw [hacc i, List.length_append] at this; omega
        have H1 := E.block_step cin cc hWp hPCp hUp tsup hBSs ct tacc tsz horun l r z ind hr hsep' hfit1
          (τ:=fun i => τ i (E.childC (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) l r z ind hr))
          (fun i => hk i _ (by rw [hkids i]; exact List.mem_cons_self))
        have H := cont (fun i => E.stepV (V i) (.old (.block l r z ind))) (fun i => A i ++ E.accR (n i) (h i) (μ i) r)
          _ _ H1 (fun i => by have := hfit i; rw [hacc i] at this; simp only [List.length_append] at *; omega)
          (fun i c hc => hk i c (by rw [hkids i]; exact List.mem_cons_of_mem _ hc))
        refine (H.cong (fun _ => rfl) (fun i => ?_)).mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
        · rw [runV_cons, hacc i, List.append_assoc]
        · simp only [swCost, GMove.rk?, GMove.shadow, BMove.rk?, List.length_append]
        · rw [hkids i, List.map_cons, List.sum_cons]
    | perm l G => exact free rfl (fun i => rfl)
    | phase l z c' => exact free rfl (fun i => rfl)

end Eng

end
end PowerSaving.RAM
end OAI
