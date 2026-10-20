import Work.Mixed.Bridge.HFTest2979
import Work.XHostT.HFW2979.GeomInstW

/-!
# x-inst: the #2979 engine theorem on the WARP network (the network of the real #2979 provider words)

The #2979 provider words are windowed packages (`KLW.GInvEW .. okW`, p-wgen / x-words), so the real network is p-a3's
`xnetGeomW`.  Here x-inst's warp heads-first certificate (`XHostTW.xnetGeomW_heads_first`, Work/XHostT/HFW2979/)
is fed into `engine2979_of_cut`, exactly as `HFTest2979.engine2979_geom` does for the plain network:
* `engine2979_geomW`: any head set `hb` (never hosted, remaining climb 48), data facts `hU` (unit price
  `unitPriceTW`), `hcut`, `hLlen`, `hσ`, `hα`;
* `engine2979_geomW_auto`: 74 replicas, classes `Fin 112`, the canonical head set `headBins S`, the head count in
  the Finset form of x-tables' `sched2979_heads`, positive residuals, any unit profile `U` with `hcut`;
* `engine2979_geomW_unit`: the same with `U = unit2979` (`hcut2979`).
Written by x-inst (wht45/lever2979, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace PB
namespace HFTest2979W
open Binary Matrix Finset RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT XHostTW

noncomputable section

universe U

variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] {nR : ℕ}
variable {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}

/-- **The #2979 engine theorem from the geometric WARP network.** -/
theorem engine2979_geomW (I : GInvW H T Sl C Ok) (Ic : Fin nR → Fin 9 → GInvW H T Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : T → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (hOk : ∀ (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT nR Sl T) (hfit : S.fits (resIcW Ic))
    (hb : Finset (Fin nR × Bool × T))
    (hunh : ∀ ρ l' q ρh j, (ρh, j) ∈ hb → S.host ρ l' q ≠ (ρh, 4, j))
    (h48 : ∀ ρ j, (ρ, j) ∈ hb → S.hrem ρ 4 j = 48)
    (U : List (ℕ × ℕ)) (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ U)
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
    (zGeom_indep base l hH) (cutA base l) (cutS l hH) (bph_cut base l hH hb)
    (fun r _ => HFTest2979.cutS_card l hH r)
    (Fintype.card (Orth (BC.LC (Fin 8) H l))) (unitPriceTW Ic S) U hU
    (xnetGeomW_heads_first I Ic hIc base l hbase hH hOk S hfit hb hunh h48) hα hcut hLlen hσ Orth.card_pos m cl k v

/-- **The #2979 engine theorem from the warp network and the table data**: 74 replicas, classes `Fin 112`, the
canonical head set, the head count in x-tables' form (`sched2979_heads`), positive residuals (`resK2979_pos` after
`resIcW Ic = resK ct resT 74 408`). -/
theorem engine2979_geomW_auto (I : GInvW H (Fin 112) Sl C Ok) (Ic : Fin 74 → Fin 9 → GInvW H (Fin 112) Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : Fin 112 → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (hOk : ∀ (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT 74 Sl (Fin 112)) (hfit : S.fits (resIcW Ic))
    (hpos : ∀ ρ l' q, 0 < resIcW Ic ρ l' q)
    (hheads : (Finset.univ.filter (fun p : Fin 74 × Bool × Fin 112 => S.hrem p.1 4 p.2 = 48)).card = 16086)
    (U : List (ℕ × ℕ)) (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ U)
    (hcut : ∀ r, (if r = 24 then 16086 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn MX.cut2979 r + (if r = 24 + 24 then 16086 else 0))
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1755647301:ℝ)/(10:ℝ)^12) (k i)) :=
  engine2979_geomW I Ic hIc base l hbase hH hOk S hfit (headBins S) (headBins_unh S) (headBins_48 S) U hU hcut
    (I2979.hLlen_card S (resIcW Ic) hfit hpos hheads) (I2979.card_live2979 _) (card_LC401 l hH) m cl k v

/-- The same with the expected unit profile `unit2979` (`hcut2979`). -/
theorem engine2979_geomW_unit (I : GInvW H (Fin 112) Sl C Ok) (Ic : Fin 74 → Fin 9 → GInvW H (Fin 112) Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : Fin 112 → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (hOk : ∀ (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT 74 Sl (Fin 112)) (hfit : S.fits (resIcW Ic))
    (hpos : ∀ ρ l' q, 0 < resIcW Ic ρ l' q)
    (hheads : (Finset.univ.filter (fun p : Fin 74 × Bool × Fin 112 => S.hrem p.1 4 p.2 = 48)).card = 16086)
    (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ unit2979)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1755647301:ℝ)/(10:ℝ)^12) (k i)) :=
  engine2979_geomW_auto I Ic hIc base l hbase hH hOk S hfit hpos hheads unit2979 hU hcut2979 m cl k v

end
end HFTest2979W
end PB
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.PB.HFTest2979W.engine2979_geomW
#print axioms OAI.PowerSaving.PB.HFTest2979W.engine2979_geomW_auto
#print axioms OAI.PowerSaving.PB.HFTest2979W.engine2979_geomW_unit
