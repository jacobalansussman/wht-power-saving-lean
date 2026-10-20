import Work.GCert.Data.Glue
import XMut.ScalarZ
import XMut.ECheckZ
import XMut.Labels

/-!
# x-mutable (13 Oct 2026): Glue: a checked MUTABLE-SOURCE gauged certificate is a gauged invocation package

As Nil Mamano's `GXD.gscalOf` / `GXD.ginvOf` (`Work/Gauge/Glue.lean`, on Jacob Sussman's `Work/GCert/Data/Glue.lean`):
* `gscalOfZ`  the gauged mutable-source scalar half, from `GS.ScalOKZ c` (kinds `okAZK`/`okBZK`, no `hx`) and the
              E7 kernel checks `GS.eCheckZ`;
* `ginvOfZ`   the gauged invocation package `GX.GInvE` (the SAME package type `GX.engine_of_ginv_401` consumes), from
              the label check `GLab.ValidE` (unchanged), via `GLab.ValidE.glabZ` and `GX.GLabZ.inv`;
* `ginvOfZ_cost`  its price: `rcost rA + Σ_k bcost (dim of copy k) + rcost rB + rcost rF`, as `ginvOf_cost`.
No `sorry`.
-/

set_option linter.unusedVariables false
set_option linter.unusedSectionVars false

namespace GXD
open SSC OAI.PowerSaving OAI.PowerSaving.CR OAI.PowerSaving.CB Matrix Finset

section
variable (c : Raw) (K : GS.ScalOKZ c)
    (es : List (GS.SPar × List Nat)) (he : ∀ e ∈ es, GS.eCheckZ c e.1 e.2 = true)

/-- **the gauged mutable-source scalar half** -/
noncomputable def gscalOfZ : GX.GScalZ (Fin c.v) (Fin c.R) (Fin c.ret.length) where
  toScalZ := K.toScalZ
  gauged := univ.filter (fun q : Fin c.R => 2 * c.v + q.val ∈ es.flatMap Prod.snd)
  hE := fun q hq t => by
    obtain ⟨e, he1, he2⟩ := List.mem_flatMap.mp (Finset.mem_filter.mp hq).2
    exact GS.eCheckZ_sound K.hv (he e he1) q he2 t

theorem gscalOfZ_gauged (q : Fin c.R) :
    q ∈ (gscalOfZ c K es he).gauged ↔ 2 * c.v + q.val ∈ es.flatMap Prod.snd := by
  show q ∈ univ.filter (fun q : Fin c.R => 2 * c.v + q.val ∈ es.flatMap Prod.snd) ↔ _
  simp

theorem reg_slotZ (k : Fin c.ret.length) (q : Fin c.R) (h : (gscalOfZ c K es he).Cc k q ≠ 0) :
    reg c k = 2 * c.v + q.val := by
  have h' : (if 2 * c.v + q.val = (c.ret.getD k.val (0, 0)).1 then (1:ℚ) else 0) ≠ 0 := h
  by_cases e : 2 * c.v + q.val = (c.ret.getD k.val (0, 0)).1
  · exact e.symm
  · rw [if_neg e] at h'; exact absurd rfl h'

variable (a : GLab.Aux) (V : GLab.ValidE c a (es.flatMap Prod.snd))

/-- the gauged label half over the mutable-source scalar half -/
noncomputable def glabOfZ : GX.GLabZ (Fin a.p.h) (gscalOfZ c K es he) :=
  V.glabZ (gscalOfZ c K es he) (roles c K.hv) GS.toCoef rfl rfl rfl rfl
    (reg c) (reg_slotZ c K es he) (reg_mem c)
    (fun q hq => (gscalOfZ_gauged c K es he q).mpr hq)

/-- **a checked mutable-source gauged certificate is a gauged invocation package** -/
noncomputable def ginvOfZ : GX.GInvE (Fin a.p.h) (Fin c.v) (Fin c.R) (Fin c.ret.length) :=
  (glabOfZ c K es he a V).inv

theorem ginvOfZ_cost (φ : ℕ → ℝ) :
    (ginvOfZ c K es he a V).cost φ
      = rcost φ a.rA + (∑ k : Fin c.ret.length, bcost φ (GLab.cdim (a.p.tab.get (a.SA.get (reg c k)))))
        + rcost φ a.rB + rcost φ a.rF := by
  obtain ⟨h1, h2, h3, h4⟩ := V.glab_ranksZ (gscalOfZ c K es he) (roles c K.hv) GS.toCoef rfl rfl rfl rfl
    (reg c) (reg_slotZ c K es he) (reg_mem c)
    (fun q hq => (gscalOfZ_gauged c K es he q).mpr hq)
  show GX.GLabZ.cost (glabOfZ c K es he a V) φ = _
  unfold GX.GLabZ.cost
  unfold glabOfZ
  rw [h1, h2, h3]
  simp only [h4]

theorem ginvOfZ_ent_dim (q : Fin c.R) :
    ((ginvOfZ c K es he a V).ent q).dim = GLab.cdim (a.p.tab.get (a.S0.get (2 * c.v + q.val))) :=
  GLab.lab_dim (GLab.tabK_good _ _ V.tab _)

end
end GXD

#print axioms GXD.gscalOfZ
#print axioms GXD.ginvOfZ
#print axioms GXD.ginvOfZ_cost
#print axioms GXD.ginvOfZ_ent_dim
