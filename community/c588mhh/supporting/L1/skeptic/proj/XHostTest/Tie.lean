import XHost.L1
import Work.GCert.Chain.WinB

/-! v-host-0 (skeptic) test: L1 for the network window `gwinB` by name, with z = unit vectors of `inl k`, k missing
from the entrance label.  Also re-states the plan form against h-L1's theorem with `#check` so its type is printed. -/

set_option linter.unusedSectionVars false
namespace OAI
namespace PowerSaving
namespace XHostTest
open Binary Matrix Finset RAM SS CB BR CR GF GX RF BG XHost
noncomputable section

theorem gwinB_L1 {B H T Sl C : Type} [Fintype B] [DecidableEq B] [Fintype H] [DecidableEq H]
    [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
    (I : GInvE H T Sl C) (hB : 1 ≤ Fintype.card B) (hH : 1 ≤ Fintype.card H) (d : RF.Orth (LB B H))
    (qq : Finset (B × H)) (q : Sl) :
    ∃ (Cm Cm' : CMat (LB B H)) (z : Fin (Fintype.card H - (I.ent q).dim) → Space (LB B H)),
      Free Cm ∧ Free Cm' ∧ Cm' * Cm = 1 ∧ Cm * Cm' = 1 ∧ LinearIndependent F z ∧
      Function.Injective z ∧
      (∀ i, ∃ k : H, k ∉ (I.ent q).s ∧ z i = eu (Sum.inl k)) ∧
      (gwinB hB hH d qq).X.Φ (I.ent q).top = Cm * blockMat z * Cm' * (gwinB hB hH d qq).X.Φ (I.ent q) := by
  obtain ⟨Cm, Cm', z, h1, h2, h3, h4, h5, h6, h7⟩ := ginvE_slot_residual I (gwinB hB hH d qq) q
  refine ⟨Cm, Cm', z, h1, h2, h3, h4, h5, h5.injective, fun i => ?_, h7⟩
  obtain ⟨k, hk, e⟩ := h6 i
  have hm : k ∈ missing (BXMap.ofConj d.1 d.2 (Orth.orth' d) qq) (I.ent q) := hk
  rw [ofConj_missing] at hm
  obtain ⟨j, hj, rfl⟩ := Finset.mem_map.mp hm
  exact ⟨j, (Finset.mem_sdiff.mp hj).2, e⟩

end
end XHostTest
end PowerSaving
end OAI

#check @OAI.PowerSaving.XHost.ginvE_slot_residual
#print axioms OAI.PowerSaving.XHostTest.gwinB_L1
