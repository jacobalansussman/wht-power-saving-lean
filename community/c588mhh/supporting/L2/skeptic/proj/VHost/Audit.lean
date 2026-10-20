import XHost.Slot

/-! v-host-1 skeptic audit of h-L2: restate the plan's L2 in my own words and compose the pieces end to end. -/

set_option linter.unusedSectionVars false
namespace OAI.PowerSaving.GF.XHost
open Binary Matrix Finset RAM Module
noncomputable section

/-- L2 (a): a permutation conjugates dir into dir along the mapped direction. -/
theorem audit_dir {α : Type} [Fintype α] [DecidableEq α] (G : APerm α) (z : Space α) :
    permMat G * dir z * permMat (apInv G) = dir ((apInv G).π z) := permMat_conj_dir G z

/-- L2 (b): a FREE address permutation sending the slot directions z onto r of the host's n pending
directions v (chosen by an injection e), and carrying the block along z onto the block along v ∘ e. -/
theorem audit_transport {α : Type} [Fintype α] [DecidableEq α] {r n : ℕ} (z : Fin r → Space α)
    (v : Fin n → Space α) (e : Fin r → Fin n) (he : Function.Injective e)
    (hz : LinearIndependent F z) (hv : LinearIndependent F v) :
    ∃ σ : APerm α, Free (permMat σ) ∧ Free (permMat (apInv σ)) ∧ (∀ i, σ.π (z i) = v (e i)) ∧
      permMat (apInv σ) * blockMat z * permMat σ = blockMat (fun i => v (e i)) := by
  obtain ⟨σ, hσ⟩ := exists_transport_into z v e he hz hv
  exact ⟨σ, Free.perm σ, Free.perm _, hσ, transport_blockMat hσ⟩

/-- End to end, one slot: a host on its idle climb ⟨G,s⟩ → ⟨G,t⟩ takes the residual of one gauged slot
(GInvE hypotheses on another window map Y), then pays ONE block of rank |t \ s| - (card H' - dim E). -/
theorem audit_one {H H' α ρ : Type} [Fintype H] [DecidableEq H] [Fintype H'] [DecidableEq H']
    [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ] (l : ρ)
    (X : BXMap H α) (G : APerm H) {s t : Finset H} (hst : s ⊆ t)
    (Y : BXMap H' α) (E : SLbl H') {N Wn : CMat α} (hN : Y.Φ E * N = 1)
    (hW : Y.Φ ⟨E.G, univ⟩ = Wn * Y.Φ E)
    (hc1 : Fintype.card H' - E.dim < (t \ s).card)
    (hc2 : (t \ s).card - (Fintype.card H' - E.dim) < Fintype.card α) :
    ∃ P, Free P ∧ GStep l (X.Φ ⟨G, s⟩) (P * X.Φ ⟨G, s⟩) (fun _ => 0) ∧
      GStep l (Wn * (P * X.Φ ⟨G, s⟩)) (X.Φ ⟨G, t⟩)
        (fun φ => φ ((t \ s).card - (Fintype.card H' - E.dim))) := by
  obtain ⟨h0, htgt, hD, hcard⟩ := idle_host X G hst
  obtain ⟨P, hP, h1⟩ := hostOK_slot Y E hN hW h0 (by rw [hcard]; omega)
  refine ⟨P, hP, gstep_free hP l _, ?_⟩
  have h := hostOK_pad l h1 hD (by rw [hcard]; omega) (by rw [hcard]; omega)
  rw [← htgt, hcard, Nat.zero_add] at h
  exact h

/-- End to end, two slots on one host (stacking), fully hosted: the remainder is free. -/
theorem audit_two_full {H H1 H2 α ρ : Type} [Fintype H] [DecidableEq H] [Fintype H1] [DecidableEq H1]
    [Fintype H2] [DecidableEq H2] [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ] (l : ρ)
    (X : BXMap H α) (G : APerm H) {s t : Finset H} (hst : s ⊆ t)
    (Y1 : BXMap H1 α) (E1 : SLbl H1) {N1 W1 : CMat α} (hN1 : Y1.Φ E1 * N1 = 1)
    (hW1 : Y1.Φ ⟨E1.G, univ⟩ = W1 * Y1.Φ E1)
    (Y2 : BXMap H2 α) (E2 : SLbl H2) {N2 W2 : CMat α} (hN2 : Y2.Φ E2 * N2 = 1)
    (hW2 : Y2.Φ ⟨E2.G, univ⟩ = W2 * Y2.Φ E2)
    (hc : (Fintype.card H1 - E1.dim) + (Fintype.card H2 - E2.dim) = (t \ s).card) :
    ∃ P1 P2, Free P1 ∧ Free P2 ∧
      GStep l (W2 * (P2 * (W1 * (P1 * X.Φ ⟨G, s⟩)))) (X.Φ ⟨G, t⟩) (fun _ => 0) := by
  obtain ⟨h0, htgt, hD, hcard⟩ := idle_host X G hst
  obtain ⟨P1, hP1, h1⟩ := hostOK_slot Y1 E1 hN1 hW1 h0 (by rw [hcard]; omega)
  obtain ⟨P2, hP2, h2⟩ := hostOK_slot Y2 E2 hN2 hW2 h1 (by rw [hcard]; omega)
  have e : 0 + (Fintype.card H1 - E1.dim) + (Fintype.card H2 - E2.dim)
      = (X.liftS t \ X.liftS s).card := by rw [hcard]; omega
  rw [e] at h2
  refine ⟨P1, P2, hP1, hP2, ?_⟩
  have h := hostOK_full l h2 hD
  rw [← htgt] at h
  exact h

/-- The bank-form route does not need C' * C = 1: hostOK_use with ANY free C, C'. -/
theorem audit_noinv {α : Type} [Fintype α] [DecidableEq α] {S : Finset α} {B Z C C' : CMat α} {n : ℕ}
    (h : HostOK S B n Z) (hC : Free C) (hC' : Free C') (R : Finset α) (hc : n + R.card ≤ S.card) :
    ∃ P, Free P ∧ HostOK S B (n + R.card) ((C * pk R * C') * (P * Z)) :=
  hostOK_use h hC hC' R rfl hc

end
end OAI.PowerSaving.GF.XHost

#print axioms OAI.PowerSaving.GF.XHost.audit_dir
#print axioms OAI.PowerSaving.GF.XHost.audit_transport
#print axioms OAI.PowerSaving.GF.XHost.audit_one
#print axioms OAI.PowerSaving.GF.XHost.audit_two_full
#print axioms OAI.PowerSaving.GF.XHost.audit_noinv
