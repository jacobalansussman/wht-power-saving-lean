import Work.PWGen.S5Pair
import Work.XHostT.HFW2979.Disp2979W
import Work.XHostT.Inst.S5HFacts
import Work.Mixed.Bridge.Inst3945H

/-!
# bB-data: the data of the design-B (#3945H) warp network (providers, dispatch, window facts, residuals, unit price)

i-final's `Work/L3945/Prov.lean` with the schedule and host table of design B (`Inst.sched3945H`, `K.DS5H`): the
providers, dispatch and window facts do not depend on the fill (same words, order and dispatch); the residual check,
`hfit`, `hpos` and the unit price are restated on the S5H table.  Original description:

* `I0W`, `I1W`: p-wgen's two #3945 provider packages (`PWG.ginvE` = extra, provider 0; `PWG.ginvB` = gauge44, provider
  1; windowed `KLW.GInvEW .. GXD.S5B.okW`) as p-a3's warp providers `XHostTW.GInvW` (the same seven fields);
* `Ic3945W`: the per-call providers by p-disp's dispatch table `D3945.disp` (x-inst's generic `D2979W.IcW`), with
  `hIc3945W` (one target-vector family, p-wgen's `PWG.tv_eq`);
* `hOk3945`: every window of the network (Shea's `win401N`, conjugation windows `gwinC`) has a warp meaning
  (p-wgen's `GXD.S5B.okC`, by p-warpsem's `WLift.ofConj`);
* `resIcW3945`: the entrance residual of every slot of every call is the host table's (p-wgen's `hent` against
  p-hostck's residual table `DS5H.resT`, two kernel checks `resE_ok`, `resB_ok`, and p-disp's `dp_disp`); hence
  `hfit3945W` (p-a3's `sched3945H_fits`) and `hpos3945W`;
* `hU3945W`: the unit price of the network is `listCost φ unit3945H` (3,240 calls of the local profile
  {1: 2482, 2: 712, 3: 218, 4: 109, 5: 20, 8: 8} plus the remaining climbs {3: 23328, 5: 4320, 48: 137952}).
i-final (wht45/lever, 19 Oct 2026); instance 3945H by bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option maxRecDepth 100000
namespace OAI
namespace PowerSaving
namespace MX
namespace L3945H
open Finset Binary Matrix RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT XHostTW

noncomputable section

/-- provider 0 of #3945's dispatch (the extra word) as a warp provider of p-a3's network. -/
def I0W : GInvW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S5B.okW :=
  ⟨PWG.ginvE.tv, PWG.ginvE.hunit, PWG.ginvE.hgap, PWG.ginvE.cost, PWG.ginvE.ent, PWG.ginvE.fwd, PWG.ginvE.bwd⟩

/-- provider 1 of #3945's dispatch (the gauge44 word) as a warp provider of p-a3's network. -/
def I1W : GInvW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S5B.okW :=
  ⟨PWG.ginvB.tv, PWG.ginvB.hunit, PWG.ginvB.hgap, PWG.ginvB.cost, PWG.ginvB.ent, PWG.ginvB.fwd, PWG.ginvB.bwd⟩

/-- the provider of every call of #3945 (p-disp's table: 1,672 extra, 1,568 gauge44). -/
def Ic3945W : Fin 360 → Fin 9 → GInvW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S5B.okW :=
  D2979W.IcW D3945.disp I0W I1W

theorem hIc3945W : ∀ ρ l, (Ic3945W ρ l).tv = I0W.tv :=
  D2979W.hIcW D3945.disp I0W I1W PWG.tv_eq

/-- **every window of the network has a warp meaning** (Shea's windows are conjugation windows `gwinC`). -/
theorem hOk3945 : ∀ (l : Fin 9) (n : ℕ) (d : Orth (BC.LC (Fin 8) (Fin 9) l)),
    GXD.S5B.okW (win401N l (Fintype.card_fin 9) n d) := by
  intro l n d
  unfold win401N
  split <;> exact GXD.S5B.okC _ _ _ _

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

/-- the extra word's slot residuals (p-wgen) are provider 0's of p-hostck's table. -/
theorem resE_ok : resEqL K.DS5H.resT 0 GXD.S5E.schedRes 408 = true := by decide +kernel
/-- the gauge44 word's slot residuals (p-wgen) are provider 1's of p-hostck's table. -/
theorem resB_ok : resEqL K.DS5H.resT 1 GXD.S5B.schedRes 408 = true := by decide +kernel

/-- **the residual of every slot of every call is the host table's.** -/
theorem resIcW3945 (ρ : Fin 360) (l : Fin 9) (q : Fin 408) :
    resIcW Ic3945W ρ l q = resK K.DS5H.ct K.DS5H.resT 360 408 ρ l q := by
  rw [Inst.resK3945H_disp]
  have hp := (D3945.disp (ρ, l)).isLt
  unfold resIcW Ic3945W D2979W.IcW
  by_cases h0 : (D3945.disp (ρ, l)).val = 0
  · rw [if_pos h0, h0]
    show Fintype.card (Fin 9) - (GXD.S5E.ginvW.ent q).dim = K.resAt K.DS5H.resT 0 q.val
    rw [GXD.S5E.hent q]
    exact resEqL_sound _ _ _ 408 resE_ok q.val q.isLt
  · have h1 : (D3945.disp (ρ, l)).val = 1 := by omega
    rw [if_neg h0, h1]
    show Fintype.card (Fin 9) - (GXD.S5B.ginvW.ent q).dim = K.resAt K.DS5H.resT 1 q.val
    rw [GXD.S5B.hent q]
    exact resEqL_sound _ _ _ 408 resB_ok q.val q.isLt

/-- **the schedule fits the providers' residuals** (p-a3's `sched3945H_fits`). -/
theorem hfit3945W : Inst.sched3945H.fits (resIcW Ic3945W) := by
  have e : resIcW Ic3945W = resK K.DS5H.ct K.DS5H.resT 360 408 := by
    funext ρ l q
    exact resIcW3945 ρ l q
  rw [e]
  exact Inst.sched3945H_fits

theorem hpos3945W : ∀ ρ l q, 0 < resIcW Ic3945W ρ l q := fun ρ l q => by
  rw [resIcW3945]
  exact Inst.resK3945H_pos ρ l q

/-- **the unit price of the design-B warp network is its profile** `unit3945H` (= #3945's). -/
theorem hU3945W (φ : ℕ → ℝ) : unitPriceTW Ic3945W Inst.sched3945H φ = listCost φ PB.unit3945H := by
  rw [D2979W.unitPriceTW_const Ic3945W Inst.sched3945H I0W.cost
        (D2979W.IcW_cost D3945.disp I0W I1W (fun φ => (PWG.same_price φ).symm)) φ,
      Inst.sched3945H_hcost]
  show ((9 * 360 : ℕ) : ℝ) * GXD.S5E.ginvW.cost φ + (23328 * φ 3 + 4320 * φ 5 + 137952 * φ 48)
    = listCost φ PB.unit3945H
  rw [GXD.S5E.ginv_cost]
  unfold GXD.S5E.invPrice
  simp only [listCost, PB.unit3945H, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  push_cast
  ring

end
end L3945H
end MX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.MX.L3945H.hIc3945W
#print axioms OAI.PowerSaving.MX.L3945H.hOk3945
#print axioms OAI.PowerSaving.MX.L3945H.resIcW3945
#print axioms OAI.PowerSaving.MX.L3945H.hfit3945W
#print axioms OAI.PowerSaving.MX.L3945H.hpos3945W
#print axioms OAI.PowerSaving.MX.L3945H.hU3945W
