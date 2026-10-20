import XHost.HostBank
import Work.Bank401.Adapt

/-!
# x-hosting, L1 → L2 seam: a slot of a gauged invocation can be hosted

`B401.gauged_Wn` (Ryan Shea's complex588 tree, `Work/Bank401/Adapt.lean`, vendored byte-identical, sha256 d1e8804a…;
entrance gauges Dugongue's, PR #401) shows that the step `Wn` of a gauged invocation package (`GInvE.fwd`'s hypotheses
`X.Φ E * N = 1`, `X.Φ E.top = Wn * X.Φ E`, `E.top = ⟨E.G, univ⟩`) is `A * pk R * A'` with `A`, `A'` free and
`|R| = card H - dim E`.  With `hostOK_use` this is the hosting step of PROOF.md §2 for one slot:

* `hostOK_slot`  a host in state `HostOK S B n Z` with room (`n + (card H - dim E) ≤ |S|`) takes the slot's residual:
                 a free `P` on the host, the invocation (slot input `P * Z`, output `Wn * (P * Z)`), and the host is in
                 state `HostOK S B (n + (card H - dim E))`.

No `sorry`.  Written for the wht43 x-hosting study, tier L2 (13 Oct 2026), prepared with assistance from Claude.
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.GF
namespace XHost
open Binary Matrix Finset
noncomputable section

variable {H α : Type} [Fintype H] [DecidableEq H] [Fintype α] [DecidableEq α]

/-- **Hosting one slot of a gauged invocation.** -/
theorem hostOK_slot (X : BXMap H α) (E : SLbl H) {N Wn : CMat α} (hN : X.Φ E * N = 1)
    (hW : X.Φ ⟨E.G, univ⟩ = Wn * X.Φ E) {S : Finset α} {B Z : CMat α} {n : ℕ} (h : HostOK S B n Z)
    (hc : n + (Fintype.card H - E.dim) ≤ S.card) :
    ∃ P, Free P ∧ HostOK S B (n + (Fintype.card H - E.dim)) (Wn * (P * Z)) := by
  obtain ⟨A, A', R, hA, hA', hR, hWn⟩ := B401.gauged_Wn X E hN hW
  rw [← hR] at hc ⊢
  exact hostOK_use h hA hA' R hWn hc

end
end XHost
end PowerSaving.GF
end OAI

#print axioms OAI.PowerSaving.GF.XHost.hostOK_slot
