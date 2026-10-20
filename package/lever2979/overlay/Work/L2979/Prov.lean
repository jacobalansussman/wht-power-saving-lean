import Work.PWGen.S2979Pair
import Work.XHostT.HFW2979.Disp2979W
import Work.XHostT.Inst.S2979Facts
import Work.Mixed.Bridge.Inst2979

/-!
# x-final: the data of the #2979 warp network (providers, dispatch, window facts, residuals, unit price)

* `I0W`, `I1W`: x-words' two #2979 provider packages (`PWG2979.ginvE` = extra, provider 0; `PWG2979.ginvB` = gauge44,
  provider 1; windowed `KLW.GInvEW .. GXD.S2979B.okW`) as p-a3's warp providers `XHostTW.GInvW` (the same seven
  fields);
* `Ic2979W`: the per-call providers by x-tables' dispatch table `D2979.disp` (x-inst's `D2979W.IcW`), with
  `hIc2979W` (one target-vector family, x-words' `PWG2979.tv_eq`);
* `hOk2979`: every window of the network (Shea's `win401N`, conjugation windows `gwinC`) has a warp meaning
  (x-words' `GXD.S2979B.okC`, by p-warpsem's `WLift.ofConj`);
* `resIcW2979`: the entrance residual of every slot of every call is the host table's (x-words' `hent` against
  x-tables' residual table `DH2979.resT`, two kernel checks `resE_ok`, `resB_ok`, and x-tables' `dp_disp`); hence
  `hfit2979W` (x-tables' `sched2979_fits`) and `hpos2979W`;
* `hU2979W`: the unit price of the network is `listCost φ unit2979` (666 calls of the local profile
  {1: 2445, 2: 706, 3: 232, 4: 108, 5: 18, 8: 8} plus the remaining climbs {1: 17310, 3: 4956, 5: 1728, 48: 28190}).
i-final's `Work/L3945/Prov.lean` with the #2979 instance.  x-final (wht45/lever2979, 19 Oct 2026); prepared with
assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option maxRecDepth 100000
namespace OAI
namespace PowerSaving
namespace MX
namespace L2979
open Finset Binary Matrix RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT XHostTW

noncomputable section

/-- provider 0 of #2979's dispatch (the extra word) as a warp provider of p-a3's network. -/
def I0W : GInvW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S2979B.okW :=
  ⟨PWG2979.ginvE.tv, PWG2979.ginvE.hunit, PWG2979.ginvE.hgap, PWG2979.ginvE.cost, PWG2979.ginvE.ent,
   PWG2979.ginvE.fwd, PWG2979.ginvE.bwd⟩

/-- provider 1 of #2979's dispatch (the gauge44 word) as a warp provider of p-a3's network. -/
def I1W : GInvW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S2979B.okW :=
  ⟨PWG2979.ginvB.tv, PWG2979.ginvB.hunit, PWG2979.ginvB.hgap, PWG2979.ginvB.cost, PWG2979.ginvB.ent,
   PWG2979.ginvB.fwd, PWG2979.ginvB.bwd⟩

/-- the provider of every call of #2979 (x-tables' table: 416 extra, 250 gauge44). -/
def Ic2979W : Fin 74 → Fin 9 → GInvW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S2979B.okW :=
  D2979W.IcW D2979.disp I0W I1W

theorem hIc2979W : ∀ ρ l, (Ic2979W ρ l).tv = I0W.tv :=
  D2979W.hIcW D2979.disp I0W I1W PWG2979.tv_eq

/-- **every window of the network has a warp meaning** (Shea's windows are conjugation windows `gwinC`). -/
theorem hOk2979 : ∀ (l : Fin 9) (n : ℕ) (d : Orth (BC.LC (Fin 8) (Fin 9) l)),
    GXD.S2979B.okW (win401N l (Fintype.card_fin 9) n d) := by
  intro l n d
  unfold win401N
  split <;> exact GXD.S2979B.okC _ _ _ _

/-- the residual list `L` agrees with provider `p` of a residual table on the slots `q < m`. -/
def resEqL (resT p : Nat) (L : List Nat) : Nat → Bool
  | 0 => true
  | m + 1 => resEqL resT p L m && Nat.beq (L.getD m 0) (K.resAt resT p m)

theorem resEqL_sound (resT p : Nat) (L : List Nat) :
    ∀ m, resEqL resT p L m = true → ∀ q < m, L.getD q 0 = K.resAt resT p q
  | 0, _, q, hq => absurd hq (Nat.not_lt_zero q)
  | m + 1, h, q, hq => by
    simp only [resEqL, Bool.and_eq_true] at h
    rcases Nat.lt_succ_iff_lt_or_eq.mp hq with hq | rfl
    · exact resEqL_sound resT p L m h.1 q hq
    · exact Nat.eq_of_beq_eq_true h.2

/-- the extra word's slot residuals (x-words) are provider 0's of x-tables' table. -/
theorem resE_ok : resEqL K.DH2979.resT 0 GXD.S2979E.schedRes 408 = true := by decide +kernel
/-- the gauge44 word's slot residuals (x-words) are provider 1's of x-tables' table. -/
theorem resB_ok : resEqL K.DH2979.resT 1 GXD.S2979B.schedRes 408 = true := by decide +kernel

/-- **the residual of every slot of every call is the host table's.** -/
theorem resIcW2979 (ρ : Fin 74) (l : Fin 9) (q : Fin 408) :
    resIcW Ic2979W ρ l q = resK K.DH2979.ct K.DH2979.resT 74 408 ρ l q := by
  rw [Inst.resK2979_disp]
  have hp := (D2979.disp (ρ, l)).isLt
  unfold resIcW Ic2979W D2979W.IcW
  by_cases h0 : (D2979.disp (ρ, l)).val = 0
  · rw [if_pos h0, h0]
    show Fintype.card (Fin 9) - (GXD.S2979E.ginvW.ent q).dim = K.resAt K.DH2979.resT 0 q.val
    rw [GXD.S2979E.hent q]
    exact resEqL_sound _ _ _ 408 resE_ok q.val q.isLt
  · have h1 : (D2979.disp (ρ, l)).val = 1 := by omega
    rw [if_neg h0, h1]
    show Fintype.card (Fin 9) - (GXD.S2979B.ginvW.ent q).dim = K.resAt K.DH2979.resT 1 q.val
    rw [GXD.S2979B.hent q]
    exact resEqL_sound _ _ _ 408 resB_ok q.val q.isLt

/-- **the schedule fits the providers' residuals** (x-tables' `sched2979_fits`). -/
theorem hfit2979W : Inst.sched2979.fits (resIcW Ic2979W) := by
  have e : resIcW Ic2979W = resK K.DH2979.ct K.DH2979.resT 74 408 := by
    funext ρ l q
    exact resIcW2979 ρ l q
  rw [e]
  exact Inst.sched2979_fits

theorem hpos2979W : ∀ ρ l q, 0 < resIcW Ic2979W ρ l q := fun ρ l q => by
  rw [resIcW2979]
  exact Inst.resK2979_pos ρ l q

/-- **the unit price of the #2979 warp network is its profile** `unit2979`. -/
theorem hU2979W (φ : ℕ → ℝ) : unitPriceTW Ic2979W Inst.sched2979 φ = listCost φ PB.unit2979 := by
  rw [D2979W.unitPriceTW_const Ic2979W Inst.sched2979 I0W.cost
        (D2979W.IcW_cost D2979.disp I0W I1W (fun φ => (PWG2979.same_price φ).symm)) φ,
      Inst.sched2979_hcost]
  show ((9 * 74 : ℕ) : ℝ) * GXD.S2979E.ginvW.cost φ + (17310 * φ 1 + 4956 * φ 3 + 1728 * φ 5 + 28190 * φ 48)
    = listCost φ PB.unit2979
  rw [GXD.S2979E.ginv_cost]
  unfold GXD.S2979E.invPrice
  simp only [listCost, PB.unit2979, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  push_cast
  ring

end
end L2979
end MX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.MX.L2979.hIc2979W
#print axioms OAI.PowerSaving.MX.L2979.hOk2979
#print axioms OAI.PowerSaving.MX.L2979.resIcW2979
#print axioms OAI.PowerSaving.MX.L2979.hfit2979W
#print axioms OAI.PowerSaving.MX.L2979.hpos2979W
#print axioms OAI.PowerSaving.MX.L2979.hU2979W
