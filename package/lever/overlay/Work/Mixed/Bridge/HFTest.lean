import Work.Mixed.Bridge.HeadsFirst
import Work.XHostT.HF.GeomInst

/-!
# p-bridge: interface test against p-hfirst's REAL modules (v1, 19 Oct 12:07) and p-a3's geometric network

`engine3945_geom`: p-hfirst's `xnetGeom_heads_first` (the heads-first certificate of p-a3's `xnetGeom` for any fitting
host schedule) fed into p-bridge's `engine3945_of_cut` with p-hfirst's own cut data (`cutA`, `cutS`, `bph_cut`).
The remaining hypotheses are exactly the data facts integration must supply (see p-bridge/NEEDS.md): the unit profile
`hU`, the profile identity `hcut` (`hcut3945` when `U = unit3945`), the head count `hLlen`, the role count `hσ`, and
`card (LC (Fin 8) H l) = 73`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace PB
namespace HFTest
open Binary Matrix Finset RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT

noncomputable section

universe U

variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] {nR : ℕ}

/-- The 24 cut lines of every head of the geometric network (p-hfirst's `cutS`). -/
lemma cutS_card (l : H) (hH : Fintype.card H = 9) (r : BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 nR) :
    (cutS l hH r).card = 24 := by
  rcases r with ⟨_, d⟩ | ⟨_, b⟩
  · simp [cutS, cutLines]
  · exact b.elim0

/-- **Interface test**: the geometric network's heads-first certificate gives the #3945 engine theorem. -/
theorem engine3945_geom (I : GInvE H T Sl C) (Ic : Fin nR → Fin 9 → GInvE H T Sl C)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : T → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (S : HostSchedT nR Sl T) (hfit : S.fits (resIc Ic))
    (hb : Finset (Fin nR × Bool × T))
    (hunh : ∀ ρ l' q ρh j, (ρh, j) ∈ hb → S.host ρ l' q ≠ (ρh, 4, j))
    (h48 : ∀ ρ j, (ρ, j) ∈ hb → S.hrem ρ 4 j = 48)
    (U : List (ℕ × ℕ)) (hU : ∀ φ, unitPriceT Ic S φ = listCost φ U)
    (hcut : ∀ r, (if r = 24 then 57312 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn MX.cut3945 r + (if r = 24 + 24 then 57312 else 0))
    (hLlen : (headsL (Γ := Orth (BC.LC (Fin 8) H l)) (nb := 0) hb).length
      = Fintype.card (Orth (BC.LC (Fin 8) H l)) * 57312)
    (hσ : Fintype.card (BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 nR)
      = Fintype.card (Orth (BC.LC (Fin 8) H l)) * 322560)
    (hα : Fintype.card (BC.LC (Fin 8) H l) = 73)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1750110901:ℝ)/(10:ℝ)^12) (k i)) :=
  engine3945_of_cut Sum.inl Sum.inl_injective (headsL hb) (headsL_nodup hb) (zGeom base l hH)
    (zGeom_indep base l hH) (cutA base l) (cutS l hH) (bph_cut base l hH hb) (fun r _ => cutS_card l hH r)
    (Fintype.card (Orth (BC.LC (Fin 8) H l))) (unitPriceT Ic S) U hU
    (xnetGeom_heads_first I Ic hIc base l hbase hH S hfit hb hunh h48) hα hcut hLlen hσ Orth.card_pos m cl k v

end
end HFTest
end PB
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.PB.HFTest.engine3945_geom
