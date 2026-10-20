import Work.Mixed.Bridge.HFTest
import Work.Mixed.Bridge.HeadsFirst3945H
import Work.XHostT.HFW2979.GeomInstW

/-!
# bB-data: the design-B (#3945H) engine theorem on the WARP network (the network of the real #3945 provider words)

i-final's `HFTest3945W.lean` with the numbers of design B (#3945's words, order and dispatch with the all-heads fill:
80,640 literal heads, figure 1.7581520e-3).  The provider words are windowed packages (`KLW.GInvEW .. okW`, p-wgen),
so the real network is p-a3's `xnetGeomW`; x-inst's generic warp heads-first certificate
(`XHostTW.xnetGeomW_heads_first`) is fed into `engine3945H_of_cut` (Work/Mixed/Bridge/HeadsFirst3945H.lean):
* `engine3945H_geomW`: any head set `hb` (never hosted, remaining climb 48), data facts `hU` (unit price
  `unitPriceTW`), `hcut`, `hLlen`, `hσ`, `hα`;
* `engine3945H_geomW_auto`: 360 replicas, classes `Fin 112`, the canonical head set `headBins S`, the head count in
  Finset form (80,640), positive residuals, any unit profile `U` with `hcut`;
* `engine3945H_geomW_unit`: the same with `U = unit3945H` (`hcut3945H`).
bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace PB
namespace HFTest3945HW
open Binary Matrix Finset RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT XHostTW

noncomputable section

universe U

variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] {nR : ℕ}
variable {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}

/-- the head-role count from the head count in Finset form (x-inst's `I2979.hLlen_card`, generic). -/
theorem hLlen_card {nR : ℕ} {Sl : Type} [Fintype Sl] [DecidableEq Sl] {Γ : Type} [Fintype Γ] [DecidableEq Γ]
    (S : HostSchedT nR Sl (Fin 112))
    (res : Fin nR → Fin 9 → Sl → ℕ) (hfit : S.fits res) (hpos : ∀ ρ l q, 0 < res ρ l q)
    {H0 : ℕ} (hc : (Finset.univ.filter (fun p : Fin nR × Bool × Fin 112 => S.hrem p.1 4 p.2 = 48)).card = H0) :
    (headsL (Γ := Γ) (nb := 0) (headBins S)).length = Fintype.card Γ * H0 := by
  rw [headsL_length, headBins_eq S res hfit hpos, hc, Nat.mul_comm]

/-- the live roles of a 360-replica network: `360 x 896` per gauge copy. -/
lemma card_live360 (Γ : Type) [Fintype Γ] :
    Fintype.card (BLiveT Γ (Fin 112) 0 360) = Fintype.card Γ * 322560 := by
  simp only [Fintype.card_sum, Fintype.card_prod, Fintype.card_fin]
  ring

/-- **The design-B (#3945H) engine theorem from the geometric WARP network.** -/
theorem engine3945H_geomW (I : GInvW H T Sl C Ok) (Ic : Fin nR → Fin 9 → GInvW H T Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : T → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (hOk : ∀ (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT nR Sl T) (hfit : S.fits (resIcW Ic))
    (hb : Finset (Fin nR × Bool × T))
    (hunh : ∀ ρ l' q ρh j, (ρh, j) ∈ hb → S.host ρ l' q ≠ (ρh, 4, j))
    (h48 : ∀ ρ j, (ρ, j) ∈ hb → S.hrem ρ 4 j = 48)
    (U : List (ℕ × ℕ)) (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ U)
    (hcut : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn MX.cut3945H r + (if r = 24 + 24 then 80640 else 0))
    (hLlen : (headsL (Γ := Orth (BC.LC (Fin 8) H l)) (nb := 0) hb).length
      = Fintype.card (Orth (BC.LC (Fin 8) H l)) * 80640)
    (hσ : Fintype.card (BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 nR)
      = Fintype.card (Orth (BC.LC (Fin 8) H l)) * 322560)
    (hα : Fintype.card (BC.LC (Fin 8) H l) = 73)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1758152001:ℝ)/(10:ℝ)^12) (k i)) :=
  engine3945H_of_cut Sum.inl Sum.inl_injective (headsL hb) (headsL_nodup hb) (zGeom base l hH)
    (zGeom_indep base l hH) (cutA base l) (cutS l hH) (bph_cut base l hH hb)
    (fun r _ => HFTest.cutS_card l hH r)
    (Fintype.card (Orth (BC.LC (Fin 8) H l))) (unitPriceTW Ic S) U hU
    (xnetGeomW_heads_first I Ic hIc base l hbase hH hOk S hfit hb hunh h48) hα hcut hLlen hσ Orth.card_pos m cl k v

/-- **The design-B (#3945H) engine theorem from the warp network and the table data**: 360 replicas, classes `Fin 112`, the
canonical head set, the head count in Finset form, positive residuals. -/
theorem engine3945H_geomW_auto (I : GInvW H (Fin 112) Sl C Ok) (Ic : Fin 360 → Fin 9 → GInvW H (Fin 112) Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : Fin 112 → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (hOk : ∀ (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT 360 Sl (Fin 112)) (hfit : S.fits (resIcW Ic))
    (hpos : ∀ ρ l' q, 0 < resIcW Ic ρ l' q)
    (hheads : (Finset.univ.filter (fun p : Fin 360 × Bool × Fin 112 => S.hrem p.1 4 p.2 = 48)).card = 80640)
    (U : List (ℕ × ℕ)) (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ U)
    (hcut : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn U r
      = BlockAccounting.histFn MX.cut3945H r + (if r = 24 + 24 then 80640 else 0))
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1758152001:ℝ)/(10:ℝ)^12) (k i)) :=
  engine3945H_geomW I Ic hIc base l hbase hH hOk S hfit (headBins S) (headBins_unh S) (headBins_48 S) U hU hcut
    (hLlen_card S (resIcW Ic) hfit hpos hheads) (card_live360 _) (card_LC401 l hH) m cl k v

/-- The same with the expected unit profile `unit3945H` (`hcut3945H`). -/
theorem engine3945H_geomW_unit (I : GInvW H (Fin 112) Sl C Ok) (Ic : Fin 360 → Fin 9 → GInvW H (Fin 112) Sl C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : Fin 112 → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (hOk : ∀ (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (S : HostSchedT 360 Sl (Fin 112)) (hfit : S.fits (resIcW Ic))
    (hpos : ∀ ρ l' q, 0 < resIcW Ic ρ l' q)
    (hheads : (Finset.univ.filter (fun p : Fin 360 × Bool × Fin 112 => S.hrem p.1 4 p.2 = 48)).card = 80640)
    (hU : ∀ φ, unitPriceTW Ic S φ = listCost φ unit3945H)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1758152001:ℝ)/(10:ℝ)^12) (k i)) :=
  engine3945H_geomW_auto I Ic hIc base l hbase hH hOk S hfit hpos hheads unit3945H hU hcut3945H m cl k v

end
end HFTest3945HW
end PB
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.PB.HFTest3945HW.engine3945H_geomW
#print axioms OAI.PowerSaving.PB.HFTest3945HW.engine3945H_geomW_auto
#print axioms OAI.PowerSaving.PB.HFTest3945HW.engine3945H_geomW_unit
