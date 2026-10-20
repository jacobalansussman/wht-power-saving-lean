import Work.XHostT.Inst.S2979Facts
import Work.XHostT.D2979Geom

/-!
# (key: bridge-host-cross) #2979's unit profile: the two providers' calls plus the remaining climbs = PROFILE.json

* `lp2979`: the paid profile of ONE call of either #2979 provider word (`PROVIDERS.json` actual_paid_profile, the same
  for `extra` and `gauge44`; our lowering's block histogram, x-gate / x-words);
* `unit2979`: #2979's whole unit profile (`allocation/PROFILE.json`, sha256 94a81974...);
* `unit2979_arith`: `416 c0 + 250 c1 + (idle climbs) = listCost φ unit2979` when each provider costs `listCost φ lp2979`
  (for either network variant, plain or warp);
* **`hU2979`**: for the plain network dispatched by #2979's table on `sched2979`, `unitPriceT .. sched2979 φ =
  listCost φ unit2979` (the bridge's `hU`), from the providers' cost lemmas.
x-tables (wht45/lever2979, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
namespace OAI.PowerSaving.XHostT.Inst
open Finset OAI.PowerSaving.XHostT
open Binary Matrix RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD

/-- the paid profile of one call of either #2979 provider (PROVIDERS.json `actual_paid_profile`). -/
def lp2979 : List (ℕ × ℕ) := [(1, 2445), (2, 706), (3, 232), (4, 108), (5, 18), (8, 8)]

/-- #2979's unit profile (`allocation/PROFILE.json`). -/
def unit2979 : List (ℕ × ℕ) :=
  [(1, 1645680), (2, 470196), (3, 159468), (4, 71928), (5, 13716), (8, 5328), (48, 28190)]

/-- the arithmetic of #2979's profile: 416 + 250 calls and the idle climbs. -/
theorem unit2979_arith (φ : ℕ → ℝ) (c0 c1 : ℝ) (h0 : c0 = listCost φ lp2979) (h1 : c1 = listCost φ lp2979) :
    416 * c0 + 250 * c1 + (17310 * φ 1 + 4956 * φ 3 + 1728 * φ 5 + 28190 * φ 48) = listCost φ unit2979 := by
  subst h0 h1
  simp only [listCost, lp2979, unit2979, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  push_cast
  ring

section Plain
variable {H C : Type} [Fintype H] [DecidableEq H] [Fintype C] [DecidableEq C]

/-- **#2979's unit price is its PROFILE.json** (plain providers dispatched by the real table, on `sched2979`). -/
theorem hU2979 (I0 I1 : GInvE H (Fin 112) (Fin 408) C) (h : I1.tv = I0.tv)
    (hc0 : ∀ φ, I0.cost φ = listCost φ lp2979) (hc1 : ∀ φ, I1.cost φ = listCost φ lp2979) (φ : ℕ → ℝ) :
    unitPriceT (dispIc (Prov2.mk' I0 I1 h) D2979.disp) sched2979 φ = listCost φ unit2979 := by
  rw [unitPriceT_disp2, D2979.nExtra, D2979.nGauge, sched2979_hcost]
  push_cast
  exact unit2979_arith φ _ _ (hc0 φ) (hc1 φ)

end Plain
end OAI.PowerSaving.XHostT.Inst

#print axioms OAI.PowerSaving.XHostT.Inst.unit2979_arith
#print axioms OAI.PowerSaving.XHostT.Inst.hU2979
