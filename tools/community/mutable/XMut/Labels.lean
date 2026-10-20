import Work.Gauge.Labels
import XMut.GInv

/-!
# x-mutable (13 Oct 2026): the label half of a gauged certificate over a MUTABLE-SOURCE scalar half

`ValidE.glabZ` / `ValidE.glab_ranksZ` = Nil Mamano's `ValidE.glab` / `ValidE.glab_ranks` (`Work/Gauge/Labels.lean`,
on Jacob Sussman's label check), word for word, with the scalar half `GScalZ` (no `hx`) in place of `GScal`.
The label check never used `hx`: it reads the scalar half only through `MA`, `MAi`, `MB`, `MBi`, `Cc` and
`gauged`.  Generated from the vendored source by a text substitution (x-mutable/tools/make_labelsZ.sh).  No `sorry`.
-/

set_option linter.unusedVariables false
set_option linter.unusedSectionVars false
set_option linter.deprecated false

namespace GLab
open SSC GXD OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.GF OAI.PowerSaving.RAM
  OAI.PowerSaving.CB Finset Matrix

/-- **THE LABEL HALF of a gauged certificate, sources mutable**: `GX.GLabZ` for its gauged scalar half; the entrance label of slot
`q` is its label in the start state.  `hg`: every slot in the gauge list is gauged in `S`. -/
noncomputable def ValidE.glabZ {c : Raw} {a : Aux} {gl : List Nat} (V : ValidE c a gl) {C : Type} [Fintype C]
    [DecidableEq C] (S : OAI.PowerSaving.GX.GScalZ (Fin c.v) (Fin c.R) C)
    {emb : CR.L3 (Fin c.v) (Fin c.R) → Nat} {un : Nat → CR.L3 (Fin c.v) (Fin c.R)}
    (Ro : Roles c.v c.R emb un) (cf : Co → Coef)
    (hMA : S.MA = matP un (mic cf (c.gatesA.flatMap Gate.adds)))
    (hMAi : S.MAi = matPinv un (mic cf (c.gatesA.flatMap Gate.adds)))
    (hMB : S.MB = matP un (mic cf (c.gatesB.flatMap Gate.adds)))
    (hMBi : S.MBi = matPinv un (mic cf (c.gatesB.flatMap Gate.adds)))
    (reg : C → Nat) (hCc : ∀ k (q : Fin c.R), S.Cc k q ≠ 0 → reg k = 2 * c.v + q.val)
    (hreg : ∀ k, reg k ∈ c.ret.map Prod.fst)
    (hg : ∀ q : Fin c.R, 2 * c.v + q.val ∈ gl → q ∈ S.gauged) : OAI.PowerSaving.GX.GLabZ (Fin a.p.h) S :=
  have E := endKE_sound V.tab V.ends
  have C0 : Ctx a.p emb un := ⟨V.tab, E.hn, Ro.hun,
    fun r hr => Ro.hin r (by have := V.hn; omega),
    fun q => by have := V.hn; have := Ro.hlt q; omega⟩
  have wA := V.runA.wf E.w0
  have wB := V.runB.wf wA
  have G : ∀ k, (GLab.lab a.p.h (a.p.tab.get k)).dim = cdim (a.p.tab.get k) :=
    fun k => lab_dim (tabK_good _ _ V.tab k)
  { tv := fun t => tvec a.p.h (c.ports.getD t.val 0)
    hunit := fun t => (E.port t.val t.isLt).unit
    hgap := fun t => (E.port t.val t.isLt).gap
    ent := fun q => a.p.labS a.S0 (emb (CR.lS q))
    lab0 := fun q => a.p.labS a.S0 (emb q)
    labA := fun q => a.p.labS a.SA (emb q)
    labB := fun q => a.p.labS a.SB (emb q)
    labF := fun q => a.p.labS a.SF (emb q)
    cen := fun k => a.p.labS a.SA (reg k)
    z0 := GLab.lab a.p.h (a.p.tab.get 0)
    hz0 := by rw [G]; exact E.z0
    h0x := fun t => by
      show (GLab.lab a.p.h (a.p.tab.get (a.S0.get (emb (CR.lX t))))).dim = 1 ∧
        (GLab.lab a.p.h (a.p.tab.get (a.S0.get (emb (CR.lX t))))).Mem _
      rw [Ro.eX t, G]
      exact ⟨(E.port t.val t.isLt).x0d, (E.port t.val t.isLt).x0m⟩
    h0y := fun t => by
      show (GLab.lab a.p.h (a.p.tab.get (a.S0.get (emb (CR.lY t))))).dim = 0
      rw [Ro.eY t, (E.port t.val t.isLt).y0, G]; exact E.z0
    h0s' := fun q => rfl
    hent := fun q hq => by
      show (GLab.lab a.p.h (a.p.tab.get (a.S0.get (emb (CR.lS q))))).dim = 0
      rw [Ro.eS q]
      rcases (E.slot q.val (by have := V.hn; have := q.isLt; omega)).1 with h0 | hm
      · rw [h0, G]; exact E.z0
      · exact absurd (hg q hm) hq
    hAy := fun t => by
      show GLab.lab a.p.h (a.p.tab.get (a.SA.get (emb (CR.lY t)))) = GLab.lab a.p.h (a.p.tab.get 0)
      rw [Ro.eY t, (E.port t.val t.isLt).yA]
    hAc := fun k q h => by
      show GLab.lab a.p.h (a.p.tab.get (a.SA.get (reg k)))
        = GLab.lab a.p.h (a.p.tab.get (a.SA.get (emb (CR.lS q))))
      rw [Ro.eS q, hCc k q h]
    hcen := fun k => by
      show 0 < (GLab.lab a.p.h (a.p.tab.get (a.SA.get (reg k)))).dim
      rw [G]; exact E.ret _ (hreg k)
    hFx := fun t => by
      show (GLab.lab a.p.h (a.p.tab.get (a.SF.get (emb (CR.lX t))))).dim = Fintype.card (Fin a.p.h)
      rw [Ro.eX t, G, Fintype.card_fin]; exact (E.port t.val t.isLt).xF
    hFs := fun q => by
      show (GLab.lab a.p.h (a.p.tab.get (a.SF.get (emb (CR.lS q))))).dim = Fintype.card (Fin a.p.h)
      rw [Ro.eS q, G, Fintype.card_fin]
      exact (E.slot q.val (by have := V.hn; have := q.isLt; omega)).2
    hFy := fun t => by
      show (GLab.lab a.p.h (a.p.tab.get (a.SF.get (emb (CR.lY t))))).dim + 1
          = Fintype.card (Fin a.p.h) ∧
        ∀ x, (GLab.lab a.p.h (a.p.tab.get (a.SF.get (emb (CR.lY t))))).Mem x → dot _ x = 0
      rw [Ro.eY t, G, Fintype.card_fin]
      exact ⟨(E.port t.val t.isLt).yFd, (E.port t.val t.isLt).yFm⟩
    rA := a.rA
    rB := a.rB
    rF := a.rF
    pA := fun {α} _ _ X hα => by
      have hα' : a.p.h < Fintype.card α := by rwa [Fintype.card_fin] at hα
      rw [hMA, hMAi]
      exact (V.runA.route C0 X hα' cf E.w0).1
    pB := fun {α} _ _ X hα => by
      have hα' : a.p.h < Fintype.card α := by rwa [Fintype.card_fin] at hα
      rw [hMB, hMBi]
      exact (V.runB.route C0 X hα' cf wA).1
    pF := fun {α} _ _ X hα => by
      have hα' : a.p.h < Fintype.card α := by rwa [Fintype.card_fin] at hα
      exact V.runF.route_id C0 X hα' wB (noAdds_sound _ V.fin) }

/-- the ranks and the copy dimensions of the gauged label half (for the price) -/
theorem ValidE.glab_ranksZ {c : Raw} {a : Aux} {gl : List Nat} (V : ValidE c a gl) {C : Type} [Fintype C]
    [DecidableEq C] (S : OAI.PowerSaving.GX.GScalZ (Fin c.v) (Fin c.R) C)
    {emb : CR.L3 (Fin c.v) (Fin c.R) → Nat} {un : Nat → CR.L3 (Fin c.v) (Fin c.R)}
    (Ro : Roles c.v c.R emb un) (cf : Co → Coef) (hMA) (hMAi) (hMB) (hMBi)
    (reg : C → Nat) (hCc) (hreg) (hg) :
    (V.glabZ S Ro cf hMA hMAi hMB hMBi reg hCc hreg hg).rA = a.rA ∧
    (V.glabZ S Ro cf hMA hMAi hMB hMBi reg hCc hreg hg).rB = a.rB ∧
    (V.glabZ S Ro cf hMA hMAi hMB hMBi reg hCc hreg hg).rF = a.rF ∧
    ∀ k, ((V.glabZ S Ro cf hMA hMAi hMB hMBi reg hCc hreg hg).cen k).dim
      = cdim (a.p.tab.get (a.SA.get (reg k))) :=
  ⟨rfl, rfl, rfl, fun k => lab_dim (tabK_good _ _ V.tab _)⟩

end GLab

#print axioms GLab.ValidE.glabZ
#print axioms GLab.ValidE.glab_ranksZ
