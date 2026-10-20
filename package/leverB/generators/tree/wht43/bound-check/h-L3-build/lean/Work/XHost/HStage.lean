import Work.XHost.HDef

/-!
# (key: bridge-host) The steps of one run of the HOSTED nine-stage word

Copy of Work/Bridge401/GNetStage.lean (Shea's complex588 tree; written by Nil Mamano's gen401.py) with same-run hosting:
the adapters also move every host of a window (`hHP`), a stage puts its slots on banks or hosts (`ghliftQ<l>`), the
idle climbs are the remaining climbs of the host chains (`hHE`), and the finals of pairs 1-3 are those of F1-F3.
The gate layers are Shea's, word for word.  Written by h-L3-build/gen/gen_hnet.py, prepared with assistance from
Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option linter.unusedSimpArgs false
set_option maxHeartbeats 4000000
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX
noncomputable section

section Net
variable {H T Sl C α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {nb : ℕ} (N : GHost401 H T Sl C α Γ nb)

/-- the free adapters before stage 0 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt0 (ρ : Fin 8) :
    GRoute (N.S0 ρ) (N.S1 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S0 ρ) (N.S1 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach ((N.G0 g).Xm.Φ (tt (N.I.tv t))) (N.Hp 0 0 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_A2 g t).1]
        exact N.hHP_A2_0 g (false, t)
      · show GReach ((N.G0 g).Xm.Φ 0) (N.Hp 0 0 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_A2 g t).2]
        exact N.hHP_A2_0 g (true, t)
      · show GReach ((N.G0 g).Xm.Φ (tt (N.I.tv t))) (N.Hp 2 0 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_A3 g t).1]
        exact N.hHP_A3_0 g (false, t)
      · show GReach ((N.G0 g).Xm.Φ 0) (N.Hp 2 0 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_A3 g t).2]
        exact N.hHP_A3_0 g (true, t)
      · show GReach ((N.G0 g).Xm.Φ (tt (N.I.tv t))) (N.Hp 4 0 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_A4 g t).1]
        exact N.hHP_A4_0 g (false, t)
      · show GReach ((N.G0 g).Xm.Φ 0) (N.Hp 4 0 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_A4 g t).2]
        exact N.hHP_A4_0 g (true, t)
      · exact N.hP0 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 0 of run `ρ` (pair 1): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage0 (ρ : Fin 8) :
    GRoute (N.S1 ρ) (N.S2 ρ) (liftQ0 (Γ := Γ) (T := T) (C := C) (N.bk ρ 0) GF)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn0 g q * N.S1 ρ (slotR (N.hostAt 0) (N.bk ρ 0) g q) = N.S2 ρ (slotR (N.hostAt 0) (N.bk ρ 0) g q) := by
    intro g q
    rcases hq : N.hostAt 0 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ0 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 0 q c (s, t) hq
      fin_cases c
      · cases s <;> exact N.hHZ_A2_0 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A3_0 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A4_0 g q _ hq
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
  have h := gistageE_fwd N.I (bid Γ T) N.G0 N.N0 N.Wn0 N.hN0 N.hW0
    (fun g q => N.S1 ρ (slotR (N.hostAt 0) (N.bk ρ 0) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv (bid Γ T) (fun d => (N.G0 d).Xm)
        (fun g q => N.S1 ρ (slotR (N.hostAt 0) (N.bk ρ 0) g q)) (cOne Γ C α))
      (mkS N.xi1 N.yi1 (fun g q => N.S2 ρ (slotR (N.hostAt 0) (N.bk ρ 0) g q)) (cOne Γ C α)) GF
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv (bid Γ T) (fun d => (N.G0 d).Xm)
        (fun g q => N.Wn0 g q * N.S1 ρ (slotR (N.hostAt 0) (N.bk ρ 0) g q)) (cOne Γ C α)
        = mkS N.xi1 N.yi1 (fun g q => N.S2 ρ (slotR (N.hostAt 0) (N.bk ρ 0) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · exact N.x01 g t
      · exact N.y01 g t
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ0 (N.hostAt 0) (N.bk ρ 0) (N.hok ρ 0) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GF_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hHZn_A2_0 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A2_0 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A3_0 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A3_0 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_0 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_0 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hZn0 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the free adapters before stage 1 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt1 (ρ : Fin 8) :
    GRoute (N.S2 ρ) (N.S3 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S2 ρ) (N.S3 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact N.hHP_A2_1 g (false, t)
      · exact N.hHP_A2_1 g (true, t)
      · exact N.hHP_A3_1 g (false, t)
      · exact N.hHP_A3_1 g (true, t)
      · exact N.hHP_A4_1 g (false, t)
      · exact N.hHP_A4_1 g (true, t)
      · exact N.hP1 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 1 of run `ρ` (pair 1): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage1 (ρ : Fin 8) :
    GRoute (N.S3 ρ) (N.S4 ρ) (liftQ1 (Γ := Γ) (T := T) (C := C) (N.bk ρ 1) GB)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn1 g q * N.S3 ρ (slotR (N.hostAt 1) (N.bk ρ 1) g q) = N.S4 ρ (slotR (N.hostAt 1) (N.bk ρ 1) g q) := by
    intro g q
    rcases hq : N.hostAt 1 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ1 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 1 q c (s, t) hq
      fin_cases c
      · cases s <;> exact N.hHZ_A2_1 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A3_1 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A4_1 g q _ hq
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
  have h := gistageE_bwd N.I N.r1 N.G1 N.N1 N.Wn1 N.hN1 N.hW1
    (fun g q => N.S3 ρ (slotR (N.hostAt 1) (N.bk ρ 1) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r1 (fun d => (N.G1 d).Xm)
        (fun g q => N.S3 ρ (slotR (N.hostAt 1) (N.bk ρ 1) g q)) (cOne Γ C α))
      (mkS N.xi2 N.yi2 (fun g q => N.S4 ρ (slotR (N.hostAt 1) (N.bk ρ 1) g q)) (cOne Γ C α)) GB
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r1 (fun d => (N.G1 d).Xm)
        (fun g q => N.Wn1 g q * N.S3 ρ (slotR (N.hostAt 1) (N.bk ρ 1) g q)) (cOne Γ C α)
        = mkS N.xi2 N.yi2 (fun g q => N.S4 ρ (slotR (N.hostAt 1) (N.bk ρ 1) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · exact N.x12 g t
      · exact N.y12 g t
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ1 (N.hostAt 1) (N.bk ρ 1) (N.hok ρ 1) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GB_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hHZn_A2_1 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A2_1 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A3_1 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A3_1 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_1 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_1 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hZn1 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the remaining climb of the host chains A2 (pair 2): ONE block of rank `hrem 0 j` per role. -/
lemma GHost401.climbA2 (ρ : Fin 8) :
    GRoute (N.S4 ρ) (N.S5 ρ) id (fun φ => (Fintype.card Γ : ℝ) * N.hcost 0 φ) := by
  have H1 := GRoute.reach_all (N.S4 ρ) (N.S5 ρ)
    (mkR4 (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ t φ => bcost φ (N.hrem 0 (false, t))) (fun _ t φ => bcost φ (N.hrem 0 (true, t))) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0)) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach (N.Hb 0 2 g (false, t)) ((N.G2 (N.r2 t g)).Xm.Φ (tt (N.I.tv t))) (fun φ => bcost φ (N.hrem 0 (false, t)))
        rw [← N.x12 g t]
        exact (N.hHE_A2 g t).1
      · show GReach (N.Hb 0 2 g (true, t)) ((N.G2 (N.r2 t g)).Xm.Φ 0) (fun φ => bcost φ (N.hrem 0 (true, t)))
        rw [← N.y12 g t]
        exact (N.hHE_A2 g t).2
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mkR4]
  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,
    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]
  ring

/-- gate layer A2 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateA2 (ρ : Fin 8) : XRoute (N.S5 ρ) (N.S5 ρ) (gscQ pABQ_2 kAQ_2) (fun _ => 0) := by
  refine gateP pABQ_2 kAQ_2 (N.S5 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- the free adapters before stage 2 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt2 (ρ : Fin 8) :
    GRoute (N.S5 ρ) (N.S6 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S5 ρ) (N.S6 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach ((N.G2 (N.r2 t g)).Xm.Φ (tt (N.I.tv t))) (N.Hp 1 2 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_B2 g t).1]
        exact N.hHP_B2_2 g (false, t)
      · show GReach ((N.G2 (N.r2 t g)).Xm.Φ 0) (N.Hp 1 2 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_B2 g t).2]
        exact N.hHP_B2_2 g (true, t)
      · exact N.hHP_A3_2 g (false, t)
      · exact N.hHP_A3_2 g (true, t)
      · exact N.hHP_A4_2 g (false, t)
      · exact N.hHP_A4_2 g (true, t)
      · exact N.hP2 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 2 of run `ρ` (pair 1): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage2 (ρ : Fin 8) :
    GRoute (N.S6 ρ) (N.S7 ρ) (liftQ2 (Γ := Γ) (T := T) (C := C) (N.bk ρ 2) GF)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn2 g q * N.S6 ρ (slotR (N.hostAt 2) (N.bk ρ 2) g q) = N.S7 ρ (slotR (N.hostAt 2) (N.bk ρ 2) g q) := by
    intro g q
    rcases hq : N.hostAt 2 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ2 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 2 q c (s, t) hq
      fin_cases c
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_B2_2 g q _ hq
      · cases s <;> exact N.hHZ_A3_2 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A4_2 g q _ hq
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
  have h := gistageE_fwd N.I N.r2 N.G2 N.N2 N.Wn2 N.hN2 N.hW2
    (fun g q => N.S6 ρ (slotR (N.hostAt 2) (N.bk ρ 2) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r2 (fun d => (N.G2 d).Xm)
        (fun g q => N.S6 ρ (slotR (N.hostAt 2) (N.bk ρ 2) g q)) (cOne Γ C α))
      (mkS N.xo2 N.yo2 (fun g q => N.S7 ρ (slotR (N.hostAt 2) (N.bk ρ 2) g q)) (cOne Γ C α)) GF
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r2 (fun d => (N.G2 d).Xm)
        (fun g q => N.Wn2 g q * N.S6 ρ (slotR (N.hostAt 2) (N.bk ρ 2) g q)) (cOne Γ C α)
        = mkS N.xo2 N.yo2 (fun g q => N.S7 ρ (slotR (N.hostAt 2) (N.bk ρ 2) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · rfl
      · rfl
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ2 (N.hostAt 2) (N.bk ρ 2) (N.hok ρ 2) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GF_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hHZn_B2_2 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_B2_2 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A3_2 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A3_2 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_2 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_2 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hZn2 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the remaining climb of the host chains B2 (pair 2): ONE block of rank `hrem 1 j` per role. -/
lemma GHost401.climbB2 (ρ : Fin 8) :
    GRoute (N.S7 ρ) (N.S8 ρ) id (fun φ => (Fintype.card Γ : ℝ) * N.hcost 1 φ) := by
  have H1 := GRoute.reach_all (N.S7 ρ) (N.S8 ρ)
    (mkR4 (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ t φ => bcost φ (N.hrem 1 (false, t))) (fun _ t φ => bcost φ (N.hrem 1 (true, t))) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0)) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach (N.Hb 1 3 g (false, t)) ((N.G3 (N.r3 t g)).Xm.Φ (tt (N.I.tv t))) (fun φ => bcost φ (N.hrem 1 (false, t)))
        rw [← N.x23 g t]
        exact (N.hHE_B2 g t).1
      · show GReach (N.Hb 1 3 g (true, t)) ((N.G3 (N.r3 t g)).Xm.Φ 0) (fun φ => bcost φ (N.hrem 1 (true, t)))
        rw [← N.y23 g t]
        exact (N.hHE_B2 g t).2
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mkR4]
  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,
    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]
  ring

/-- gate layer B2 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateB2 (ρ : Fin 8) : XRoute (N.S8 ρ) (N.S8 ρ) (gscQ pABQ_2 kBQ_2) (fun _ => 0) := by
  refine gateP pABQ_2 kBQ_2 (N.S8 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact N.x23 g t
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (N.y23 g t).symm
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- the free adapters before stage 3 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt3 (ρ : Fin 8) :
    GRoute (N.S8 ρ) (N.S9 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S8 ρ) (N.S9 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · show GReach ((N.G2 (N.r2 t g)).Xm.Φ 1) (N.Hp 6 3 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_F1 g t).1]
        exact N.hHP_F1_3 g (false, t)
      · show GReach ((N.G2 (N.r2 t g)).Xm.Φ (1 + tt (N.I.tv t))) (N.Hp 6 3 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_F1 g t).2]
        exact N.hHP_F1_3 g (true, t)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact N.hHP_A3_3 g (false, t)
      · exact N.hHP_A3_3 g (true, t)
      · exact N.hHP_A4_3 g (false, t)
      · exact N.hHP_A4_3 g (true, t)
      · exact N.hP3 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 3 of run `ρ` (pair 2): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage3 (ρ : Fin 8) :
    GRoute (N.S9 ρ) (N.S10 ρ) (liftQ3 (Γ := Γ) (T := T) (C := C) (N.bk ρ 3) GB)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn3 g q * N.S9 ρ (slotR (N.hostAt 3) (N.bk ρ 3) g q) = N.S10 ρ (slotR (N.hostAt 3) (N.bk ρ 3) g q) := by
    intro g q
    rcases hq : N.hostAt 3 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ3 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 3 q c (s, t) hq
      fin_cases c
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A3_3 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A4_3 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_F1_3 g q _ hq
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
  have h := gistageE_bwd N.I N.r3 N.G3 N.N3 N.Wn3 N.hN3 N.hW3
    (fun g q => N.S9 ρ (slotR (N.hostAt 3) (N.bk ρ 3) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r3 (fun d => (N.G3 d).Xm)
        (fun g q => N.S9 ρ (slotR (N.hostAt 3) (N.bk ρ 3) g q)) (cOne Γ C α))
      (mkS N.xi4 N.yi4 (fun g q => N.S10 ρ (slotR (N.hostAt 3) (N.bk ρ 3) g q)) (cOne Γ C α)) GB
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r3 (fun d => (N.G3 d).Xm)
        (fun g q => N.Wn3 g q * N.S9 ρ (slotR (N.hostAt 3) (N.bk ρ 3) g q)) (cOne Γ C α)
        = mkS N.xi4 N.yi4 (fun g q => N.S10 ρ (slotR (N.hostAt 3) (N.bk ρ 3) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · exact N.x34 g t
      · exact N.y34 g t
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ3 (N.hostAt 3) (N.bk ρ 3) (N.hok ρ 3) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GB_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact (N.hHZn_F1_3 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F1_3 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hHZn_A3_3 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A3_3 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_3 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_3 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hZn3 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the remaining climb of the host chains A3 (pair 3): ONE block of rank `hrem 2 j` per role. -/
lemma GHost401.climbA3 (ρ : Fin 8) :
    GRoute (N.S10 ρ) (N.S11 ρ) id (fun φ => (Fintype.card Γ : ℝ) * N.hcost 2 φ) := by
  have H1 := GRoute.reach_all (N.S10 ρ) (N.S11 ρ)
    (mkR4 (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ t φ => bcost φ (N.hrem 2 (false, t))) (fun _ t φ => bcost φ (N.hrem 2 (true, t))) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0)) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach (N.Hb 2 4 g (false, t)) ((N.G4 (N.r4 t g)).Xm.Φ (tt (N.I.tv t))) (fun φ => bcost φ (N.hrem 2 (false, t)))
        rw [← N.x34 g t]
        exact (N.hHE_A3 g t).1
      · show GReach (N.Hb 2 4 g (true, t)) ((N.G4 (N.r4 t g)).Xm.Φ 0) (fun φ => bcost φ (N.hrem 2 (true, t)))
        rw [← N.y34 g t]
        exact (N.hHE_A3 g t).2
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mkR4]
  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,
    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]
  ring

/-- gate layer A3 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateA3 (ρ : Fin 8) : XRoute (N.S11 ρ) (N.S11 ρ) (gscQ pABQ_3 kAQ_3) (fun _ => 0) := by
  refine gateP pABQ_3 kAQ_3 (N.S11 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- the free adapters before stage 4 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt4 (ρ : Fin 8) :
    GRoute (N.S11 ρ) (N.S12 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S11 ρ) (N.S12 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact N.hHP_F1_4 g (false, t)
      · exact N.hHP_F1_4 g (true, t)
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach ((N.G4 (N.r4 t g)).Xm.Φ (tt (N.I.tv t))) (N.Hp 3 4 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_B3 g t).1]
        exact N.hHP_B3_4 g (false, t)
      · show GReach ((N.G4 (N.r4 t g)).Xm.Φ 0) (N.Hp 3 4 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_B3 g t).2]
        exact N.hHP_B3_4 g (true, t)
      · exact N.hHP_A4_4 g (false, t)
      · exact N.hHP_A4_4 g (true, t)
      · exact N.hP4 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 4 of run `ρ` (pair 2): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage4 (ρ : Fin 8) :
    GRoute (N.S12 ρ) (N.S13 ρ) (liftQ4 (Γ := Γ) (T := T) (C := C) (N.bk ρ 4) GF)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn4 g q * N.S12 ρ (slotR (N.hostAt 4) (N.bk ρ 4) g q) = N.S13 ρ (slotR (N.hostAt 4) (N.bk ρ 4) g q) := by
    intro g q
    rcases hq : N.hostAt 4 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ4 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 4 q c (s, t) hq
      fin_cases c
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_B3_4 g q _ hq
      · cases s <;> exact N.hHZ_A4_4 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_F1_4 g q _ hq
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
  have h := gistageE_fwd N.I N.r4 N.G4 N.N4 N.Wn4 N.hN4 N.hW4
    (fun g q => N.S12 ρ (slotR (N.hostAt 4) (N.bk ρ 4) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r4 (fun d => (N.G4 d).Xm)
        (fun g q => N.S12 ρ (slotR (N.hostAt 4) (N.bk ρ 4) g q)) (cOne Γ C α))
      (mkS N.xo4 N.yo4 (fun g q => N.S13 ρ (slotR (N.hostAt 4) (N.bk ρ 4) g q)) (cOne Γ C α)) GF
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r4 (fun d => (N.G4 d).Xm)
        (fun g q => N.Wn4 g q * N.S12 ρ (slotR (N.hostAt 4) (N.bk ρ 4) g q)) (cOne Γ C α)
        = mkS N.xo4 N.yo4 (fun g q => N.S13 ρ (slotR (N.hostAt 4) (N.bk ρ 4) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · rfl
      · rfl
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ4 (N.hostAt 4) (N.bk ρ 4) (N.hok ρ 4) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GF_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact (N.hHZn_F1_4 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F1_4 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hHZn_B3_4 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_B3_4 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_4 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_4 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hZn4 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the remaining climb of the host chains B3 (pair 3): ONE block of rank `hrem 3 j` per role. -/
lemma GHost401.climbB3 (ρ : Fin 8) :
    GRoute (N.S13 ρ) (N.S14 ρ) id (fun φ => (Fintype.card Γ : ℝ) * N.hcost 3 φ) := by
  have H1 := GRoute.reach_all (N.S13 ρ) (N.S14 ρ)
    (mkR4 (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ t φ => bcost φ (N.hrem 3 (false, t))) (fun _ t φ => bcost φ (N.hrem 3 (true, t))) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0)) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach (N.Hb 3 5 g (false, t)) ((N.G5 (N.r5 t g)).Xm.Φ (tt (N.I.tv t))) (fun φ => bcost φ (N.hrem 3 (false, t)))
        rw [← N.x45 g t]
        exact (N.hHE_B3 g t).1
      · show GReach (N.Hb 3 5 g (true, t)) ((N.G5 (N.r5 t g)).Xm.Φ 0) (fun φ => bcost φ (N.hrem 3 (true, t)))
        rw [← N.y45 g t]
        exact (N.hHE_B3 g t).2
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mkR4]
  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,
    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]
  ring

/-- gate layer B3 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateB3 (ρ : Fin 8) : XRoute (N.S14 ρ) (N.S14 ρ) (gscQ pABQ_3 kBQ_3) (fun _ => 0) := by
  refine gateP pABQ_3 kBQ_3 (N.S14 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact N.x45 g t
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (N.y45 g t).symm
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- the free adapters before stage 5 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt5 (ρ : Fin 8) :
    GRoute (N.S14 ρ) (N.S15 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S14 ρ) (N.S15 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact N.hHP_F1_5 g (false, t)
      · exact N.hHP_F1_5 g (true, t)
      · show GReach ((N.G4 (N.r4 t g)).Xm.Φ 1) (N.Hp 7 5 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_F2 g t).1]
        exact N.hHP_F2_5 g (false, t)
      · show GReach ((N.G4 (N.r4 t g)).Xm.Φ (1 + tt (N.I.tv t))) (N.Hp 7 5 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_F2 g t).2]
        exact N.hHP_F2_5 g (true, t)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact N.hHP_A4_5 g (false, t)
      · exact N.hHP_A4_5 g (true, t)
      · exact N.hP5 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 5 of run `ρ` (pair 3): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage5 (ρ : Fin 8) :
    GRoute (N.S15 ρ) (N.S16 ρ) (liftQ5 (Γ := Γ) (T := T) (C := C) (N.bk ρ 5) GB)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn5 g q * N.S15 ρ (slotR (N.hostAt 5) (N.bk ρ 5) g q) = N.S16 ρ (slotR (N.hostAt 5) (N.bk ρ 5) g q) := by
    intro g q
    rcases hq : N.hostAt 5 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ5 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 5 q c (s, t) hq
      fin_cases c
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_A4_5 g q _ hq
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_F1_5 g q _ hq
      · cases s <;> exact N.hHZ_F2_5 g q _ hq
      · exact absurd hw (by decide)
  have h := gistageE_bwd N.I N.r5 N.G5 N.N5 N.Wn5 N.hN5 N.hW5
    (fun g q => N.S15 ρ (slotR (N.hostAt 5) (N.bk ρ 5) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r5 (fun d => (N.G5 d).Xm)
        (fun g q => N.S15 ρ (slotR (N.hostAt 5) (N.bk ρ 5) g q)) (cOne Γ C α))
      (mkS N.xi6 N.yi6 (fun g q => N.S16 ρ (slotR (N.hostAt 5) (N.bk ρ 5) g q)) (cOne Γ C α)) GB
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r5 (fun d => (N.G5 d).Xm)
        (fun g q => N.Wn5 g q * N.S15 ρ (slotR (N.hostAt 5) (N.bk ρ 5) g q)) (cOne Γ C α)
        = mkS N.xi6 N.yi6 (fun g q => N.S16 ρ (slotR (N.hostAt 5) (N.bk ρ 5) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · exact N.x56 g t
      · exact N.y56 g t
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ5 (N.hostAt 5) (N.bk ρ 5) (N.hok ρ 5) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GB_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact (N.hHZn_F1_5 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F1_5 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_5 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_5 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hHZn_A4_5 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_A4_5 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hZn5 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the remaining climb of the host chains A4 (pair 4): ONE block of rank `hrem 4 j` per role. -/
lemma GHost401.climbA4 (ρ : Fin 8) :
    GRoute (N.S16 ρ) (N.S17 ρ) id (fun φ => (Fintype.card Γ : ℝ) * N.hcost 4 φ) := by
  have H1 := GRoute.reach_all (N.S16 ρ) (N.S17 ρ)
    (mkR4 (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ t φ => bcost φ (N.hrem 4 (false, t))) (fun _ t φ => bcost φ (N.hrem 4 (true, t))) (fun _ _ _ => 0) (fun _ _ _ => 0)) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach (N.Hb 4 6 g (false, t)) ((N.G6 (N.r6 t g)).Xm.Φ (tt (N.I.tv t))) (fun φ => bcost φ (N.hrem 4 (false, t)))
        rw [← N.x56 g t]
        exact (N.hHE_A4 g t).1
      · show GReach (N.Hb 4 6 g (true, t)) ((N.G6 (N.r6 t g)).Xm.Φ 0) (fun φ => bcost φ (N.hrem 4 (true, t)))
        rw [← N.y56 g t]
        exact (N.hHE_A4 g t).2
      · exact GReach.refl _
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mkR4]
  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,
    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]
  ring

/-- gate layer A4 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateA4 (ρ : Fin 8) : XRoute (N.S17 ρ) (N.S17 ρ) (gscQ pABQ_4 kAQ_4) (fun _ => 0) := by
  refine gateP pABQ_4 kAQ_4 (N.S17 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- the free adapters before stage 6 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt6 (ρ : Fin 8) :
    GRoute (N.S17 ρ) (N.S18 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S17 ρ) (N.S18 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact N.hHP_F1_6 g (false, t)
      · exact N.hHP_F1_6 g (true, t)
      · exact N.hHP_F2_6 g (false, t)
      · exact N.hHP_F2_6 g (true, t)
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach ((N.G6 (N.r6 t g)).Xm.Φ (tt (N.I.tv t))) (N.Hp 5 6 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_B4 g t).1]
        exact N.hHP_B4_6 g (false, t)
      · show GReach ((N.G6 (N.r6 t g)).Xm.Φ 0) (N.Hp 5 6 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_B4 g t).2]
        exact N.hHP_B4_6 g (true, t)
      · exact N.hP6 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 6 of run `ρ` (pair 3): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage6 (ρ : Fin 8) :
    GRoute (N.S18 ρ) (N.S19 ρ) (liftQ6 (Γ := Γ) (T := T) (C := C) (N.bk ρ 6) GF)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn6 g q * N.S18 ρ (slotR (N.hostAt 6) (N.bk ρ 6) g q) = N.S19 ρ (slotR (N.hostAt 6) (N.bk ρ 6) g q) := by
    intro g q
    rcases hq : N.hostAt 6 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ6 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 6 q c (s, t) hq
      fin_cases c
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_B4_6 g q _ hq
      · cases s <;> exact N.hHZ_F1_6 g q _ hq
      · cases s <;> exact N.hHZ_F2_6 g q _ hq
      · exact absurd hw (by decide)
  have h := gistageE_fwd N.I N.r6 N.G6 N.N6 N.Wn6 N.hN6 N.hW6
    (fun g q => N.S18 ρ (slotR (N.hostAt 6) (N.bk ρ 6) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r6 (fun d => (N.G6 d).Xm)
        (fun g q => N.S18 ρ (slotR (N.hostAt 6) (N.bk ρ 6) g q)) (cOne Γ C α))
      (mkS N.xo6 N.yo6 (fun g q => N.S19 ρ (slotR (N.hostAt 6) (N.bk ρ 6) g q)) (cOne Γ C α)) GF
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r6 (fun d => (N.G6 d).Xm)
        (fun g q => N.Wn6 g q * N.S18 ρ (slotR (N.hostAt 6) (N.bk ρ 6) g q)) (cOne Γ C α)
        = mkS N.xo6 N.yo6 (fun g q => N.S19 ρ (slotR (N.hostAt 6) (N.bk ρ 6) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · rfl
      · rfl
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ6 (N.hostAt 6) (N.bk ρ 6) (N.hok ρ 6) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GF_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact (N.hHZn_F1_6 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F1_6 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_6 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_6 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hHZn_B4_6 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_B4_6 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hZn6 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the remaining climb of the host chains B4 (pair 4): ONE block of rank `hrem 5 j` per role. -/
lemma GHost401.climbB4 (ρ : Fin 8) :
    GRoute (N.S19 ρ) (N.S20 ρ) id (fun φ => (Fintype.card Γ : ℝ) * N.hcost 5 φ) := by
  have H1 := GRoute.reach_all (N.S19 ρ) (N.S20 ρ)
    (mkR4 (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ _ _ => 0) (fun _ t φ => bcost φ (N.hrem 5 (false, t))) (fun _ t φ => bcost φ (N.hrem 5 (true, t))) (fun _ _ _ => 0) (fun _ _ _ => 0)) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · exact GReach.refl _
      · show GReach (N.Hb 5 7 g (false, t)) ((N.G7 (N.r7 t g)).Xm.Φ (tt (N.I.tv t))) (fun φ => bcost φ (N.hrem 5 (false, t)))
        rw [← N.x67 g t]
        exact (N.hHE_B4 g t).1
      · show GReach (N.Hb 5 7 g (true, t)) ((N.G7 (N.r7 t g)).Xm.Φ 0) (fun φ => bcost φ (N.hrem 5 (true, t)))
        rw [← N.y67 g t]
        exact (N.hHE_B4 g t).2
      · exact GReach.refl _
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mkR4]
  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,
    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]
  ring

/-- gate layer B4 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateB4 (ρ : Fin 8) : XRoute (N.S20 ρ) (N.S20 ρ) (gscQ pABQ_4 kBQ_4) (fun _ => 0) := by
  refine gateP pABQ_4 kBQ_4 (N.S20 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact N.x67 g t
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (N.y67 g t).symm
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- the free adapters before stage 7 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt7 (ρ : Fin 8) :
    GRoute (N.S20 ρ) (N.S21 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S20 ρ) (N.S21 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact N.hHP_F1_7 g (false, t)
      · exact N.hHP_F1_7 g (true, t)
      · exact N.hHP_F2_7 g (false, t)
      · exact N.hHP_F2_7 g (true, t)
      · show GReach ((N.G6 (N.r6 t g)).Xm.Φ 1) (N.Hp 8 7 g (false, t)) (fun _ => 0)
        rw [← (N.hH0_F3 g t).1]
        exact N.hHP_F3_7 g (false, t)
      · show GReach ((N.G6 (N.r6 t g)).Xm.Φ (1 + tt (N.I.tv t))) (N.Hp 8 7 g (true, t)) (fun _ => 0)
        rw [← (N.hH0_F3 g t).2]
        exact N.hHP_F3_7 g (true, t)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact N.hP7 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 7 of run `ρ` (pair 4): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage7 (ρ : Fin 8) :
    GRoute (N.S21 ρ) (N.S22 ρ) (liftQ7 (Γ := Γ) (T := T) (C := C) (N.bk ρ 7) GB)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn7 g q * N.S21 ρ (slotR (N.hostAt 7) (N.bk ρ 7) g q) = N.S22 ρ (slotR (N.hostAt 7) (N.bk ρ 7) g q) := by
    intro g q
    rcases hq : N.hostAt 7 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ7 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 7 q c (s, t) hq
      fin_cases c
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_F1_7 g q _ hq
      · cases s <;> exact N.hHZ_F2_7 g q _ hq
      · cases s <;> exact N.hHZ_F3_7 g q _ hq
  have h := gistageE_bwd N.I N.r7 N.G7 N.N7 N.Wn7 N.hN7 N.hW7
    (fun g q => N.S21 ρ (slotR (N.hostAt 7) (N.bk ρ 7) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r7 (fun d => (N.G7 d).Xm)
        (fun g q => N.S21 ρ (slotR (N.hostAt 7) (N.bk ρ 7) g q)) (cOne Γ C α))
      (mkS N.xi8 N.yi8 (fun g q => N.S22 ρ (slotR (N.hostAt 7) (N.bk ρ 7) g q)) (cOne Γ C α)) GB
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r7 (fun d => (N.G7 d).Xm)
        (fun g q => N.Wn7 g q * N.S21 ρ (slotR (N.hostAt 7) (N.bk ρ 7) g q)) (cOne Γ C α)
        = mkS N.xi8 N.yi8 (fun g q => N.S22 ρ (slotR (N.hostAt 7) (N.bk ρ 7) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · exact N.x78 g t
      · exact N.y78 g t
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ7 (N.hostAt 7) (N.bk ρ 7) (N.hok ρ 7) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GB_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact (N.hHZn_F1_7 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F1_7 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_7 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_7 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F3_7 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F3_7 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hZn7 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the free adapters before stage 8 of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/
lemma GHost401.adapt8 (ρ : Fin 8) :
    GRoute (N.S22 ρ) (N.S23 ρ) id (fun _ => 0) := by
  have H1 := GRoute.reach_all (N.S22 ρ) (N.S23 ρ) (fun _ _ => 0) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact N.hHP_F1_8 g (false, t)
      · exact N.hHP_F1_8 g (true, t)
      · exact N.hHP_F2_8 g (false, t)
      · exact N.hHP_F2_8 g (true, t)
      · exact N.hHP_F3_8 g (false, t)
      · exact N.hHP_F3_8 g (true, t)
      · exact GReach.refl _
      · exact GReach.refl _
      · exact N.hP8 ρ g b
      · exact GReach.refl _)
  exact H1.cast (fun φ => by simp)

/-- stage 8 of run `ρ` (pair 4): the slot `q` of `g` on its bank or on its host. -/
lemma GHost401.stage8 (ρ : Fin 8) :
    GRoute (N.S23 ρ) (N.S24 ρ) (liftQ8 (Γ := Γ) (T := T) (C := C) (N.bk ρ 8) GF)
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
  have hUS : ∀ g q, N.Wn8 g q * N.S23 ρ (slotR (N.hostAt 8) (N.bk ρ 8) g q) = N.S24 ρ (slotR (N.hostAt 8) (N.bk ρ 8) g q) := by
    intro g q
    rcases hq : N.hostAt 8 q with _ | ⟨c, s, t⟩
    · rw [slotR_bank hq]
      exact N.hZ8 ρ g q hq
    · rw [slotR_host hq]
      have hw := N.hostWin 8 q c (s, t) hq
      fin_cases c
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · exact absurd hw (by decide)
      · cases s <;> exact N.hHZ_F1_8 g q _ hq
      · cases s <;> exact N.hHZ_F2_8 g q _ hq
      · cases s <;> exact N.hHZ_F3_8 g q _ hq
  have h := gistageE_fwd N.I N.r8 N.G8 N.N8 N.Wn8 N.hN8 N.hW8
    (fun g q => N.S23 ρ (slotR (N.hostAt 8) (N.bk ρ 8) g q)) (cOne Γ C α) (cOne Γ C α)
  have h' : GRoute (bIn N.I.tv N.r8 (fun d => (N.G8 d).Xm)
        (fun g q => N.S23 ρ (slotR (N.hostAt 8) (N.bk ρ 8) g q)) (cOne Γ C α))
      (mkS N.xo8 N.yo8 (fun g q => N.S24 ρ (slotR (N.hostAt 8) (N.bk ρ 8) g q)) (cOne Γ C α)) GF
      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by
    have e : bOut N.I.tv N.r8 (fun d => (N.G8 d).Xm)
        (fun g q => N.Wn8 g q * N.S23 ρ (slotR (N.hostAt 8) (N.bk ρ 8) g q)) (cOne Γ C α)
        = mkS N.xo8 N.yo8 (fun g q => N.S24 ρ (slotR (N.hostAt 8) (N.bk ρ 8) g q)) (cOne Γ C α) := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · rfl
      · rfl
      · exact hUS g q
      · rfl
    rw [← e]
    exact h
  refine ghliftQ8 (N.hostAt 8) (N.bk ρ 8) (N.hok ρ 8) h' (fun g t => rfl) (fun g t => rfl)
    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_
    (fun _ _ _ => rfl) GF_loc
  intro r hr
  rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
  · exact (N.hHZn_F1_8 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F1_8 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_8 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F2_8 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F3_8 g (false, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact (N.hHZn_F3_8 g (true, t) (fun q hq => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jX (g,t), rfl⟩) hr
  · exact absurd ((cover_iff _ _).mpr ⟨jY (g,t), rfl⟩) hr
  · exact (N.hZn8 ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr
      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm
  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr

/-- the finals of the replica of run `ρ`: the remaining climbs of F1-F3 and pair 4's final block. -/
lemma GHost401.fin (ρ : Fin 8) :
    GRoute (N.S24 ρ) (N.S25 ρ) id (fun φ => (Fintype.card Γ : ℝ) *
      (N.hcost 6 φ + N.hcost 7 φ + N.hcost 8 φ + (Fintype.card T : ℝ) * (N.cX4 φ + N.cY4 φ))) := by
  have H1 := GRoute.reach_all (N.S24 ρ) (N.S25 ρ)
    (mkR4 (fun _ t φ => bcost φ (N.hrem 6 (false, t))) (fun _ t φ => bcost φ (N.hrem 6 (true, t))) (fun _ t φ => bcost φ (N.hrem 7 (false, t))) (fun _ t φ => bcost φ (N.hrem 7 (true, t))) (fun _ t φ => bcost φ (N.hrem 8 (false, t))) (fun _ t φ => bcost φ (N.hrem 8 (true, t))) (fun _ _ => N.cX4) (fun _ _ => N.cY4) (fun _ _ _ => 0) (fun _ _ _ => 0)) (by
      rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
      · exact (N.hHE_F1 g t).1
      · exact (N.hHE_F1 g t).2
      · exact (N.hHE_F2 g t).1
      · exact (N.hHE_F2 g t).2
      · exact (N.hHE_F3 g t).1
      · exact (N.hHE_F3 g t).2
      · exact GReach.ofX (N.fX4 g t)
      · exact GReach.ofX (N.fY4 g t)
      · exact GReach.refl _
      · exact GReach.refl _)
  refine H1.cast (fun φ => ?_)
  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]
  rw [sum_mkR4]
  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,
    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]
  ring

/-- gate layer C4 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateC4 (ρ : Fin 8) : XRoute (N.S25 ρ) (N.S25 ρ) (gscQ pCQ_4 kCQ_4) (fun _ => 0) := by
  refine gateP pCQ_4 kCQ_4 (N.S25 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · rfl
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- gate layer C3 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateC3 (ρ : Fin 8) : XRoute (N.S25 ρ) (N.S25 ρ) (gscQ pCQ_3 kCQ_3) (fun _ => 0) := by
  refine gateP pCQ_3 kCQ_3 (N.S25 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · rfl
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim

/-- gate layer C2 of run `ρ` (Shea's, word for word): twins have equal frames. -/
lemma GHost401.gateC2 (ρ : Fin 8) : XRoute (N.S25 ρ) (N.S25 ρ) (gscQ pCQ_2 kCQ_2) (fun _ => 0) := by
  refine gateP pCQ_2 kCQ_2 (N.S25 ρ) ?_
  rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩) hk
  · exact (hk rfl).elim
  · rfl
  · rfl
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim
  · exact (hk rfl).elim

end Net
end
end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.GHost401.stage0
#print axioms OAI.PowerSaving.XHost.GHost401.stage1
#print axioms OAI.PowerSaving.XHost.GHost401.stage2
#print axioms OAI.PowerSaving.XHost.GHost401.stage3
#print axioms OAI.PowerSaving.XHost.GHost401.stage4
#print axioms OAI.PowerSaving.XHost.GHost401.stage5
#print axioms OAI.PowerSaving.XHost.GHost401.stage6
#print axioms OAI.PowerSaving.XHost.GHost401.stage7
#print axioms OAI.PowerSaving.XHost.GHost401.stage8
#print axioms OAI.PowerSaving.XHost.GHost401.fin
