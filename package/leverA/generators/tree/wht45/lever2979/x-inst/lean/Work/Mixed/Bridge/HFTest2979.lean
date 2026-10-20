import Work.Mixed.Bridge.HeadsFirst2979
import Work.XHostT.HF.GeomInst
import Work.XHostT.HF.Inst2979

/-!
# x-inst: interface test of the #2979 instance against p-hfirst's REAL modules and p-a3's geometric network

p-bridge's `HFTest.lean` with the #2979 instance (round wht45/lever2979):
* `engine2979_geom`: p-hfirst's `xnetGeom_heads_first` (the heads-first certificate of p-a3's `xnetGeom` for any
  fitting host schedule) fed into `engine2979_of_cut` with p-hfirst's own cut data (`cutA`, `cutS`, `bph_cut`); the
  remaining hypotheses are the data facts (`hU`, `hcut`, `hLlen`, `hσ`, `hα`);
* `engine2979_geom_auto`: the same for the #2979 shape (74 replicas, classes `Fin 112`) with the CANONICAL head set
  `headBins S` and every count hypothesis discharged: the head count from the remainder tree (`cntKind rt 48 4 74 =
  16086`, one `decide +kernel` on x-tables' data), the live-role count (`card_live2979`), the label space
  (`card_LC401`), and `hcut2979` for the expected unit profile `unit2979`.  What is left are the schedule `S` with
  `hfit`, positive residuals, `hrem` (`rfl` for p-a3's `schedOfTabOK`), the count, and the unit profile `hU`.
Written by x-inst (wht45/lever2979, 19 Oct 2026) on p-bridge's file; prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace PB
namespace HFTest2979
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

/-- **Interface test**: the geometric network's heads-first certificate gives the #2979 engine theorem. -/
theorem engine2979_geom (I : GInvE H T Sl C) (Ic : Fin nR → Fin 9 → GInvE H T Sl C)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : T → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (S : HostSchedT nR Sl T) (hfit : S.fits (resIc Ic))
    (hb : Finset (Fin nR × Bool × T))
    (hunh : ∀ ρ l' q ρh j, (ρh, j) ∈ hb → S.host ρ l' q ≠ (ρh, 4, j))
    (h48 : ∀ ρ j, (ρ, j) ∈ hb → S.hrem ρ 4 j = 48)
    (U : List (ℕ × ℕ)) (hU : ∀ φ, unitPriceT Ic S φ = listCost φ U)
    (hcut : ∀ r, (if r = 24 then 16086 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn MX.cut2979 r + (if r = 24 + 24 then 16086 else 0))
    (hLlen : (headsL (Γ := Orth (BC.LC (Fin 8) H l)) (nb := 0) hb).length
      = Fintype.card (Orth (BC.LC (Fin 8) H l)) * 16086)
    (hσ : Fintype.card (BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 nR)
      = Fintype.card (Orth (BC.LC (Fin 8) H l)) * 66304)
    (hα : Fintype.card (BC.LC (Fin 8) H l) = 73)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1755647301:ℝ)/(10:ℝ)^12) (k i)) :=
  engine2979_of_cut Sum.inl Sum.inl_injective (headsL hb) (headsL_nodup hb) (zGeom base l hH)
    (zGeom_indep base l hH) (cutA base l) (cutS l hH) (bph_cut base l hH hb) (fun r _ => cutS_card l hH r)
    (Fintype.card (Orth (BC.LC (Fin 8) H l))) (unitPriceT Ic S) U hU
    (xnetGeom_heads_first I Ic hIc base l hbase hH S hfit hb hunh h48) hα hcut hLlen hσ Orth.card_pos m cl k v

/-- **The #2979 engine theorem from the network and the table data** (74 replicas, classes `Fin 112`, canonical
head set `headBins S`): only data facts left. -/
theorem engine2979_geom_auto (I : GInvE H (Fin 112) Sl C) (Ic : Fin 74 → Fin 9 → GInvE H (Fin 112) Sl C)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : Fin 112 → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (S : HostSchedT 74 Sl (Fin 112)) (hfit : S.fits (resIc Ic))
    (hpos : ∀ ρ l' q, 0 < resIc Ic ρ l' q) (rt : K.CT)
    (hrem : ∀ ρh c j, S.hrem ρh c j
      = K.remAt rt (2016 * ρh.val + 224 * c.val + (if j.1 then 112 else 0) + j.2.val))
    (hcnt : I2979.cntKind rt 48 4 74 = 16086)
    (U : List (ℕ × ℕ)) (hU : ∀ φ, unitPriceT Ic S φ = listCost φ U)
    (hcut : ∀ r, (if r = 24 then 16086 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn MX.cut2979 r + (if r = 24 + 24 then 16086 else 0))
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1755647301:ℝ)/(10:ℝ)^12) (k i)) :=
  engine2979_geom I Ic hIc base l hbase hH S hfit (headBins S) (headBins_unh S) (headBins_48 S) U hU hcut
    (I2979.hLlen S (resIc Ic) hfit hpos rt hrem hcnt) (I2979.card_live2979 _) (card_LC401 l hH) m cl k v

/-- The same with the expected unit profile `unit2979` (`hcut2979` discharges the profile identity). -/
theorem engine2979_geom_unit (I : GInvE H (Fin 112) Sl C) (Ic : Fin 74 → Fin 9 → GInvE H (Fin 112) Sl C)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : Fin 112 → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (S : HostSchedT 74 Sl (Fin 112)) (hfit : S.fits (resIc Ic))
    (hpos : ∀ ρ l' q, 0 < resIc Ic ρ l' q) (rt : K.CT)
    (hrem : ∀ ρh c j, S.hrem ρh c j
      = K.remAt rt (2016 * ρh.val + 224 * c.val + (if j.1 then 112 else 0) + j.2.val))
    (hcnt : I2979.cntKind rt 48 4 74 = 16086)
    (hU : ∀ φ, unitPriceT Ic S φ = listCost φ unit2979)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1755647301:ℝ)/(10:ℝ)^12) (k i)) :=
  engine2979_geom_auto I Ic hIc base l hbase hH S hfit hpos rt hrem hcnt unit2979 hU hcut2979 m cl k v

end
end HFTest2979
end PB
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.PB.HFTest2979.engine2979_geom
#print axioms OAI.PowerSaving.PB.HFTest2979.engine2979_geom_auto
#print axioms OAI.PowerSaving.PB.HFTest2979.engine2979_geom_unit
