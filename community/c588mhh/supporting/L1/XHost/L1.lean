import Work.Gauge.Defs
import Work.GFrame.Labels.ConjB

/-!
# x-hosting, L1: the residual of a gauged slot is ONE block on the missing coordinates between free adapters

Hosting (DaysSky, CrocSwap/integer-mult-bounds PR #442, `research/idle-axis-hosting/PROOF.md` §1(a)) needs: a dirty
invocation acts on the stream owning a helper slot as an exact coordinate factor `P_g C_block P_g⁻¹` on `w`
coordinates, whatever the stream held.  In our engine the "whatever the stream held" half is already the `∀ Z` of
`GInvE.fwd`/`GInvE.bwd` (Nil Mamano / DFT Bridge3, `Work/Gauge/Defs.lean`, on Dugongue's #401 entrance gauges).  This
file supplies the other half for every slot `q` and every window `W`:

    W.X.Φ (ent q).top = C * blockMat z * C' * W.X.Φ (ent q),
    C = kg G', C' = kg G'⁻¹ (G' = W.X.liftG (ent q).G),  Free C, Free C',  C' * C = 1 = C * C',
    z = the unit vectors of the address coordinates `W.X.liftS univ \ W.X.liftS (ent q).s`   (rank |H| - dim (ent q)),

and, for the windows of the network (`BXMap.ofConj d q`: `U ↦ d (U ⊕ span q)`), those coordinates are exactly the
images `Sum.inl k` of the label coordinates `k ∉ (ent q).s` missing from the entrance label.

Main results (namespace `OAI.PowerSaving.XHost`):
* `slot_residual`            the identity above for ANY stage-B frame map and ANY label;
* `Wn_eq_resW`, `resW_pk`    ANY `Wn` meeting `GInvE.fwd`'s hypotheses equals it; its coordinate-kernel form `C * pk R * C'`;
* `ginvE_slot_residual`      L1 as stated in the plan, for every `GInvE` slot and window, with `z` indexed by
                             `Fin (card H - dim (ent q))`;
* `ofConj_missing`, `ofConj_resDir`   the missing coordinates of a network window are `inl` of the missing label coords;
* `ginvE_fwd_residual`, `ginvE_bwd_residual`   §1(a) in our engine: the invocation leaves ANY slot content `Z q` as
                             `C * blockMat z * C' * Z q` (the residual acts on whatever the slot held);
* `residual_step`            on a bank role the residual is ONE engine step of rank `card H - dim (ent q)`;
* `hosted_remaining_climb`   on a host whose idle climb (in the slot's adapters) is ONE block along `Fin.append z y`,
                             after the residual the rest is ONE block of rank `card y` (x-hosting's `Core.lean`
                             `host_remaining_climb_twisted`, restated here on L1's own adapters).

The mathematical content is NOT new: the pk form of L1 is already proved in Ryan Shea's complex588 tree (at 0e61af6)
as `B401.gauged_step` / `B401.gauged_Wn` (`Work/Bank401/Adapt.lean`), `GX.entW_spec` / `GX.entW_form`
(`Work/Bridge401/Sched.lean`) and `GX.gauge_step` (`Work/Gauge/Step.lean`; DFT Bridge3 / Nil Mamano), all on Dugongue's
#401 entrance gauges and on Jacob Sussman's `core_nested` (`Work/GFrame/Clifford/Core.lean`, Apache-2.0).  This file
restates it in the form the hosting lemmas consume: `blockMat` along explicit unit vectors, `C' * C = 1 = C * C'`,
the identification of the directions for the network windows, and the `GInvE`-level forward/backward routes.
No `sorry`.

Written for the wht43 hosting study (13 Oct 2026, agent h-L1), prepared with assistance from Claude.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB BR CR GF GX RF
noncomputable section

section Residual
variable {H α : Type} [Fintype H] [DecidableEq H] [Fintype α] [DecidableEq α]

/-- the address coordinates a slot with entrance label `U` still has to climb in the window `X`. -/
def missing (X : BXMap H α) (U : SLbl H) : Finset α := X.liftS univ \ X.liftS U.s

/-- the residual directions: the unit vectors of the missing coordinates. -/
def resDir (X : BXMap H α) (U : SLbl H) : Fin (missing X U).card → Space α :=
  (OBase.canonical α).fam (missing X U)

/-- the left adapter of the residual. -/
def resC (X : BXMap H α) (U : SLbl H) : CMat α := kg (X.liftG U.G)

/-- the right adapter of the residual. -/
def resC' (X : BXMap H α) (U : SLbl H) : CMat α := kg (apInv (X.liftG U.G))

/-- the residual of the slot (`Wn` of `GInvE.fwd`). -/
def resW (X : BXMap H α) (U : SLbl H) : CMat α := resC X U * blockMat (resDir X U) * resC' X U

lemma missing_sub (X : BXMap H α) (U : SLbl H) : X.liftS U.s ⊆ X.liftS univ :=
  X.mono (Finset.subset_univ _)

/-- its rank: the number of label coordinates missing from `U`. -/
theorem card_missing (X : BXMap H α) (U : SLbl H) : (missing X U).card = Fintype.card H - U.dim := by
  unfold missing
  rw [X.card (Finset.subset_univ U.s), Finset.card_sdiff_of_subset (Finset.subset_univ _), Finset.card_univ]
  rfl

theorem free_resC (X : BXMap H α) (U : SLbl H) : Free (resC X U) := free_kg _

theorem free_resC' (X : BXMap H α) (U : SLbl H) : Free (resC' X U) := free_kg _

theorem resC'_mul_resC (X : BXMap H α) (U : SLbl H) : resC' X U * resC X U = 1 := kg_inv_mul _

theorem resC_mul_resC' (X : BXMap H α) (U : SLbl H) : resC X U * resC' X U = 1 := kg_mul_inv _

theorem resDir_indep (X : BXMap H α) (U : SLbl H) : LinearIndependent F (resDir X U) :=
  (OBase.canonical α).fam_indep _

/-- every residual direction is the unit vector of a missing coordinate. -/
theorem resDir_unit (X : BXMap H α) (U : SLbl H) (i : Fin (missing X U).card) :
    ∃ k ∈ missing X U, resDir X U i = eu k :=
  ⟨((missing X U).equivFin.symm i).1, ((missing X U).equivFin.symm i).2, rfl⟩

/-- the residual block is the coordinate kernel on the missing coordinates (the hub's `C_block` on `S'`). -/
theorem blockMat_resDir (X : BXMap H α) (U : SLbl H) : blockMat (resDir X U) = pk (missing X U) :=
  (pk_eq_blockMat _).symm

/-- **L1, any frame map, any label.** -/
theorem slot_residual (X : BXMap H α) (U : SLbl H) :
    X.Φ U.top = resC X U * blockMat (resDir X U) * resC' X U * X.Φ U := by
  show core (X.liftG U.G) (X.liftS univ) = _ * core (X.liftG U.G) (X.liftS U.s)
  exact core_nested _ (missing_sub X U)

theorem slot_residual' (X : BXMap H α) (U : SLbl H) : X.Φ U.top = resW X U * X.Φ U :=
  slot_residual X U

/-- a right inverse of the entrance frame (the `N q` of `GInvE.fwd`; as Shea's `entN_spec`). -/
def resN (X : BXMap H α) (U : SLbl H) : CMat α :=
  kg (X.liftG U.G) * (pk (X.liftS U.s) * shift (indv (X.liftS U.s))) * kg (apInv (X.liftG U.G))

theorem resN_spec (X : BXMap H α) (U : SLbl H) : X.Φ U * resN X U = 1 := by
  show core (X.liftG U.G) (X.liftS U.s) * _ = 1
  unfold core resN
  have e : kg (X.liftG U.G) * pk (X.liftS U.s) * kg (apInv (X.liftG U.G))
        * (kg (X.liftG U.G) * (pk (X.liftS U.s) * shift (indv (X.liftS U.s))) * kg (apInv (X.liftG U.G)))
      = kg (X.liftG U.G) * ((pk (X.liftS U.s) * ((kg (apInv (X.liftG U.G)) * kg (X.liftG U.G))
          * pk (X.liftS U.s))) * shift (indv (X.liftS U.s))) * kg (apInv (X.liftG U.G)) := by
    simp only [Matrix.mul_assoc]
  rw [e, kg_inv_mul, Matrix.one_mul, pk_sq, shift_sq, Matrix.mul_one, kg_mul_inv]

/-- the residual in coordinate-kernel form (the form of Shea's `entW`, `B401.gauged_Wn` and h-L2's `HostOK`). -/
theorem resW_pk (X : BXMap H α) (U : SLbl H) : resW X U = resC X U * pk (missing X U) * resC' X U := by
  rw [resW, blockMat_resDir]

/-- **uniqueness**: ANY `Wn` meeting `GInvE.fwd`'s hypotheses for the slot (some right inverse `N` of the entrance
frame, `X.Φ U.top = Wn * X.Φ U`) IS the residual `C * blockMat z * C'` (so L1 applies to whatever `Wn` a generated
network uses, e.g. Shea's `entW`). -/
theorem Wn_eq_resW (X : BXMap H α) (U : SLbl H) {N Wn : CMat α} (hN : X.Φ U * N = 1)
    (hW : X.Φ U.top = Wn * X.Φ U) : Wn = resW X U := by
  calc Wn = Wn * (X.Φ U * N) := by rw [hN, Matrix.mul_one]
    _ = X.Φ U.top * N := by rw [hW, Matrix.mul_assoc]
    _ = resW X U * (X.Φ U * N) := by rw [slot_residual', Matrix.mul_assoc]
    _ = resW X U := by rw [hN, Matrix.mul_one]

/-- re-indexing a block along `Fin.cast` changes nothing. -/
lemma blockMat_cast {m n : ℕ} (h : m = n) (z : Fin m → Space α) :
    blockMat (fun i : Fin n => z (Fin.cast h.symm i)) = blockMat z := by
  subst h; rfl

end Residual

section Window
variable {H β : Type} [Fintype H] [DecidableEq H] [Fintype β] [DecidableEq β]

/-- **The missing coordinates of a network window** (`BXMap.ofConj d q`, the window of `gwinB`/`gwinC`): the images
`Sum.inl k` of the label coordinates `k` missing from `U`, whatever `d` and `q`. -/
theorem ofConj_missing (d : Matrix (H ⊕ β) (H ⊕ β) F) (hd : dᵀ * d = 1) (hd' : d * dᵀ = 1) (q : Finset β)
    (U : SLbl H) :
    missing (BXMap.ofConj d hd hd' q) U = (univ \ U.s).map Function.Embedding.inl := by
  ext k
  cases k with
  | inl i => simp [missing, BXMap.ofConj, Finset.inl_mem_disjSum]
  | inr b => simp [missing, BXMap.ofConj, Finset.inr_mem_disjSum]

/-- the residual directions of a network window are the unit vectors `eu (inl k)`, `k ∉ U.s`. -/
theorem ofConj_resDir (d : Matrix (H ⊕ β) (H ⊕ β) F) (hd : dᵀ * d = 1) (hd' : d * dᵀ = 1) (q : Finset β)
    (U : SLbl H) (i : Fin (missing (BXMap.ofConj d hd hd' q) U).card) :
    ∃ k : H, k ∉ U.s ∧ resDir (BXMap.ofConj d hd hd' q) U i = eu (Sum.inl k) := by
  obtain ⟨k, hk, e⟩ := resDir_unit (BXMap.ofConj d hd hd' q) U i
  rw [ofConj_missing] at hk
  obtain ⟨j, hj, rfl⟩ := Finset.mem_map.mp hk
  exact ⟨j, (Finset.mem_sdiff.mp hj).2, e⟩

end Window

section Slots
variable {H T Sl C α : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]

/-- **L1** (the plan's statement): for every slot `q` of a gauged invocation and every window `W`, the step of the slot
from its entrance frame to the top frame is `C * blockMat z * C' `, `C`, `C'` free and mutually inverse, `z` the
`card H - dim (ent q)` unit vectors of the missing address coordinates. -/
theorem ginvE_slot_residual (I : GInvE H T Sl C) (W : GWin H α) (q : Sl) :
    ∃ (Cm Cm' : CMat α) (z : Fin (Fintype.card H - (I.ent q).dim) → Space α),
      Free Cm ∧ Free Cm' ∧ Cm' * Cm = 1 ∧ Cm * Cm' = 1 ∧ LinearIndependent F z ∧
      (∀ i, ∃ k ∈ W.X.liftS univ \ W.X.liftS (I.ent q).s, z i = eu k) ∧
      W.X.Φ (I.ent q).top = Cm * blockMat z * Cm' * W.X.Φ (I.ent q) := by
  have hc := card_missing W.X (I.ent q)
  refine ⟨resC W.X (I.ent q), resC' W.X (I.ent q),
    fun i => resDir W.X (I.ent q) (Fin.cast hc.symm i),
    free_resC _ _, free_resC' _ _, resC'_mul_resC _ _, resC_mul_resC' _ _, ?_, ?_, ?_⟩
  · exact (resDir_indep W.X (I.ent q)).comp _ (Fin.cast_injective _)
  · intro i
    exact resDir_unit W.X (I.ent q) _
  · rw [blockMat_cast hc (resDir W.X (I.ent q))]
    exact slot_residual W.X (I.ent q)

/-- **§1(a) in our engine, forward**: the invocation leaves ANY slot content `Z q` as `resW * Z q`. -/
theorem ginvE_fwd_residual (I : GInvE H T Sl C) (W : GWin H α) (Z : Sl → CMat α) (c0 c1 : C → CMat α) :
    GRoute (iin W.Xm I.tv Z c0) (iout W.Xm I.tv (fun q => resW W.X (I.ent q) * Z q) c1) shearF I.cost :=
  I.fwd W (fun q => resN W.X (I.ent q)) (fun q => resW W.X (I.ent q))
    (fun q => resN_spec _ _) (fun q => slot_residual' _ _) Z c0 c1

/-- **§1(a) in our engine, backward** (the inverse-transpose program). -/
theorem ginvE_bwd_residual (I : GInvE H T Sl C) (W : GWin H α) (Z : Sl → CMat α) (c0 c1 : C → CMat α) :
    GRoute (iin W.Xm I.tv Z c0) (iout W.Xm I.tv (fun q => resW W.X (I.ent q) * Z q) c1) shearB I.cost :=
  I.bwd W (fun q => resN W.X (I.ent q)) (fun q => resW W.X (I.ent q))
    (fun q => resN_spec _ _) (fun q => slot_residual' _ _) Z c0 c1

end Slots

section Steps
variable {H α ρ : Type} [Fintype H] [DecidableEq H] [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- on a bank role the residual is ONE engine step of rank `card H - dim U` (what a bank pays for it). -/
theorem residual_step (l : ρ) (X : BXMap H α) (U : SLbl H) (Y : CMat α)
    (h1 : 1 ≤ (missing X U).card) (h2 : (missing X U).card < Fintype.card α) :
    GStep l Y (resW X U * Y) (fun φ => φ (Fintype.card H - U.dim)) := by
  have h := gstep_of_factor l Y (resDir X U) (resDir_indep X U) (free_resC X U) (free_resC' X U) h1 h2
  exact h.cast (fun φ => by rw [card_missing])

lemma blockMat_comm' {r s : ℕ} (z : Fin r → Space α) (y : Fin s → Space α) :
    blockMat z * blockMat y = blockMat y * blockMat z := by
  rw [blockMat_wrap, blockMat_wrap]
  apply wrap_comm

lemma blockMat_append' {r s : ℕ} (z : Fin r → Space α) (y : Fin s → Space α) :
    blockMat (Fin.append z y) = blockMat z * blockMat y := by
  unfold blockMat
  rw [List.ofFn_add, List.prod_append]
  simp [Fin.append_right]

/-- **the remaining climb of a host** that received the residual of a slot: if the host's idle climb, in the slot's
adapters, is ONE block along `Fin.append (resDir X U) y`, then after the residual the rest is ONE block of rank `s`. -/
theorem hosted_remaining_climb (l : ρ) (X : BXMap H α) (U : SLbl H) (Y : CMat α) {s : ℕ} (y : Fin s → Space α)
    (ind : LinearIndependent F y) (h1 : 1 ≤ s) (h2 : s < Fintype.card α) :
    GStep l (resW X U * Y)
      (resC X U * blockMat (Fin.append (resDir X U) y) * resC' X U * Y) (fun φ => φ s) := by
  have h := gstep_of_factor l (resW X U * Y) y ind (free_resC X U) (free_resC' X U) h1 h2
  have e : resC X U * blockMat y * resC' X U * (resW X U * Y)
      = resC X U * blockMat (Fin.append (resDir X U) y) * resC' X U * Y := by
    unfold resW
    calc resC X U * blockMat y * resC' X U * (resC X U * blockMat (resDir X U) * resC' X U * Y)
        = resC X U * blockMat y * (resC' X U * resC X U) * blockMat (resDir X U) * resC' X U * Y := by
          simp only [Matrix.mul_assoc]
      _ = resC X U * (blockMat y * blockMat (resDir X U)) * resC' X U * Y := by
          rw [resC'_mul_resC]; simp only [Matrix.mul_one, Matrix.mul_assoc]
      _ = resC X U * blockMat (Fin.append (resDir X U) y) * resC' X U * Y := by
          rw [blockMat_comm' y (resDir X U), ← blockMat_append']
  rw [e] at h
  exact h

end Steps

end
end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.slot_residual
#print axioms OAI.PowerSaving.XHost.Wn_eq_resW
#print axioms OAI.PowerSaving.XHost.resW_pk
#print axioms OAI.PowerSaving.XHost.ginvE_slot_residual
#print axioms OAI.PowerSaving.XHost.ofConj_missing
#print axioms OAI.PowerSaving.XHost.ofConj_resDir
#print axioms OAI.PowerSaving.XHost.ginvE_fwd_residual
#print axioms OAI.PowerSaving.XHost.ginvE_bwd_residual
#print axioms OAI.PowerSaving.XHost.residual_step
#print axioms OAI.PowerSaving.XHost.hosted_remaining_climb
