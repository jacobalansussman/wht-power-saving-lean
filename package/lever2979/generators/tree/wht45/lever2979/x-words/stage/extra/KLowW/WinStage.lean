import Work.KLowW.Win
import Work.Bridge401.StageE

/-!
# (key: p-warpchk) A stage of windowed gauged invocations (`GInvEW`): drop-in for `gistageE_fwd/bwd`

Copy of `GX.gistageE_fwd` / `gistageE_bwd` / `GInvE.climbX` / `climbY` (Work/Bridge401/StageE.lean, DFT Bridge3 on
Jacob Sussman's `gistage_fwd`) for the windowed package `KLW.GInvEW`: the only change is the hypothesis
`hok : ∀ d, Ok (M d)` on the stage's windows, passed to the package's routes.  A network built on these asks `Ok` of
its windows and nothing else.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace KLW
open OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.RAM OAI.PowerSaving.SS OAI.PowerSaving.CB
  OAI.PowerSaving.RF OAI.PowerSaving.BR OAI.PowerSaving.GF OAI.PowerSaving.GX Matrix Finset
noncomputable section

section Pkg
variable {H T Sl C α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {Ok : ∀ {α : Type} [Fintype α] [DecidableEq α], GWin H α → Prop} (I : GInvEW H T Sl C Ok)

/-- the second idle climb of an x role (copy of `GInvE.climbX`). -/
lemma GInvEW.climbX (Xm : XMap H α) (t : T) :
    XReach (Xm.Φ (tt (I.tv t))) (Xm.Φ 1) (fun φ => bcost φ (Fintype.card H - 1)) := by
  obtain ⟨p, hp⟩ := I.hgap t
  obtain ⟨A, k, hA⟩ := exists_base_one (I.tv t) (I.hunit t) p hp
  have h := Xm.up A {k} univ (Finset.subset_univ _)
  rw [mat_single, OBase.mat_univ, hA] at h
  refine h.cast (fun φ => ?_)
  rw [Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_univ,
    Finset.card_singleton]

/-- the second idle climb of a y role (copy of `GInvE.climbY`). -/
lemma GInvEW.climbY (Xm : XMap H α) (t : T) :
    XReach (Xm.Φ 0) (Xm.Φ (1 + tt (I.tv t))) (fun φ => bcost φ (Fintype.card H - 1)) := by
  obtain ⟨p, hp⟩ := I.hgap t
  obtain ⟨A, k, hA⟩ := exists_base_one (I.tv t) (I.hunit t) p hp
  have h := Xm.up A ∅ (univ \ {k}) (Finset.empty_subset _)
  rw [mat_empty, mat_compl_single, hA] at h
  refine h.cast (fun φ => ?_)
  rw [Finset.sdiff_empty, Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_univ,
    Finset.card_singleton]

/-- **A forward stage of windowed gauged invocations.** -/
theorem gistageEW_fwd (εX : T → Γ ≃ Γ) (M : Γ → GWin H α) (hok : ∀ d, Ok (M d)) (N Wn : Γ → Sl → CMat α)
    (hN : ∀ d q, (M d).X.Φ (I.ent q) * N d q = 1)
    (hW : ∀ d q, (M d).X.Φ (I.ent q).top = Wn d q * (M d).X.Φ (I.ent q))
    (Z : Γ → Sl → CMat α) (c0 c1 : Γ → C → CMat α) :
    GRoute (bIn I.tv εX (fun d => (M d).Xm) Z c0)
      (bOut I.tv εX (fun d => (M d).Xm) (fun g q => Wn g q * Z g q) c1)
      GF (fun φ => (Fintype.card Γ : ℝ) * I.cost φ) := by
  refine gsstage εX (Equiv.refl Γ) (h := shearF) (cost := I.cost) (fun d => ?_) (fun d x => ?_)
  · rw [bIn_sem, bOut_sem]
    exact I.fwd (M d) (hok d) (N d) (Wn d) (hN d) (hW d) _ _ _
  · funext b
    rcases b with ((t|t)|(q|k)) <;> rfl

/-- **A backward stage of windowed gauged invocations.** -/
theorem gistageEW_bwd (εX : T → Γ ≃ Γ) (M : Γ → GWin H α) (hok : ∀ d, Ok (M d)) (N Wn : Γ → Sl → CMat α)
    (hN : ∀ d q, (M d).X.Φ (I.ent q) * N d q = 1)
    (hW : ∀ d q, (M d).X.Φ (I.ent q).top = Wn d q * (M d).X.Φ (I.ent q))
    (Z : Γ → Sl → CMat α) (c0 c1 : Γ → C → CMat α) :
    GRoute (bIn I.tv εX (fun d => (M d).Xm) Z c0)
      (bOut I.tv εX (fun d => (M d).Xm) (fun g q => Wn g q * Z g q) c1)
      GB (fun φ => (Fintype.card Γ : ℝ) * I.cost φ) := by
  refine gsstage εX (Equiv.refl Γ) (h := shearB) (cost := I.cost) (fun d => ?_) (fun d x => ?_)
  · rw [bIn_sem, bOut_sem]
    exact I.bwd (M d) (hok d) (N d) (Wn d) (hN d) (hW d) _ _ _
  · funext b
    rcases b with ((t|t)|(q|k)) <;> rfl

end Pkg
end
end KLW

#print axioms KLW.gistageEW_fwd
#print axioms KLW.gistageEW_bwd
