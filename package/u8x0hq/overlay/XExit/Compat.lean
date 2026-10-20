import XExit.Halves

/-!
# (key: x-exit) Conservativity: Shea's `GLab` is the exit-labelled `GLabX` with NO exit slots

`GLabX.ofLab L` (`E = ∅`, tail labels = `labB`, `rT = []`): same price (`ofLab_cost`), same exit labels
(`(ent q).top`, `ofLab_exq`), same bank residuals (`fitsX_ofLab ↔ S.fits L.ent`).  So everything that passes
through `GLab.inv` / `GInvE` today passes through `GLabX.inv` / `GInvX` with the same figures.
Credits as XExit/Halves.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace XEX
open Binary Matrix Finset RAM SS CB BR CR GF RF GX B401
noncomputable section

section Compat
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

lemma iE_empty (i : L3 T Sl) : iE (∅ : Finset Sl) i = 0 := by
  rcases i with (a|q) <;> simp [iE]

/-- a gauged label half is an exit-labelled one with no exit slots (tail gates `1`, rest gates Shea's end gates). -/
def GLabX.ofLab {S : GScal T Sl C} (L : GLab H S) : GLabX H S where
  tv := L.tv
  hunit := L.hunit
  hgap := L.hgap
  ent := L.ent
  E := ∅
  exl := fun q => (L.ent q).top
  hexG := fun q hq => by simp at hq
  hexs := fun q hq => by simp at hq
  lab0 := L.lab0
  labA := L.labA
  labB := L.labB
  labT := L.labB
  labF := L.labF
  cen := L.cen
  z0 := L.z0
  hz0 := L.hz0
  h0x := L.h0x
  h0y := L.h0y
  h0s' := L.h0s'
  hent := L.hent
  hAy := L.hAy
  hAc := L.hAc
  hcen := L.hcen
  hFx := L.hFx
  hFs := fun q _ => L.hFs q
  hFe := fun q hq => by simp at hq
  hFy := L.hFy
  Gt := 1
  Gr := g8 S.toScal.Mti
  Ht := 1
  Hr := (h8 S.toScal.Mt)ᵀ
  hGsplit := Matrix.mul_one _
  hHsplit := Matrix.mul_one _
  hGt := fun i j e h => absurd (Matrix.one_apply_ne e) h
  hHt := fun i j e h => absurd (Matrix.one_apply_ne e) h
  hGr := fun i j e h => by
    have f := g8_facts S.toScal.hiy
    obtain ⟨-, h2, h3⟩ := offdiag_s f.1 f.2.1 f.2.2 h e
    exact ⟨⟨h2, iE_empty i⟩, ⟨h3, iE_empty j⟩⟩
  hHr := fun i j e h => by
    have f := h8_facts S.toScal.hmy
    have h' : h8 S.toScal.Mt j i ≠ 0 := h
    obtain ⟨-, h2, h3⟩ := offdiag_s f.1 f.2.1 f.2.2 h' (Ne.symm e)
    exact ⟨⟨h3, iE_empty i⟩, ⟨h2, iE_empty j⟩⟩
  rA := L.rA
  rB := L.rB
  rT := []
  rF := L.rF
  pA := L.pA
  pB := L.pB
  pT := fun X _ => GRoute.refl.cast (fun φ => (rcost_nil φ).symm)
  pF := L.pF

lemma GLabX.ofLab_cost {S : GScal T Sl C} (L : GLab H S) (φ : ℕ → ℝ) :
    (GLabX.ofLab L).cost φ = L.cost φ := by
  simp only [GLabX.cost, GLab.cost, GLabX.ofLab, rcost_nil, add_zero]

lemma GLabX.ofLab_exq {S : GScal T Sl C} (L : GLab H S) (q : Sl) :
    (GLabX.ofLab L).exq q = (L.ent q).top := by
  simp [GLabX.exq, GLabX.ofLab]

/-- the exit residuals of `ofLab L` are Shea's residuals `card H - dim ent`. -/
theorem fitsX_ofLab {S : GScal T Sl C} (L : GLab H S) (Sc : BankSched Sl) :
    fitsX Sc (GLabX.ofLab L).inv ↔ Sc.fits L.ent := by
  have e : (fun q => ((GLabX.ofLab L).inv.ex q).dim - ((GLabX.ofLab L).inv.ent q).dim)
      = fun q => Fintype.card H - (L.ent q).dim := by
    funext q
    show ((GLabX.ofLab L).exq q).dim - (L.ent q).dim = _
    rw [GLabX.ofLab_exq, SLbl.top_dim]
  unfold fitsX BankSched.fits
  rw [e]

end Compat
end
end XEX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XEX.GLabX.ofLab
#print axioms OAI.PowerSaving.XEX.GLabX.ofLab_cost
#print axioms OAI.PowerSaving.XEX.fitsX_ofLab
