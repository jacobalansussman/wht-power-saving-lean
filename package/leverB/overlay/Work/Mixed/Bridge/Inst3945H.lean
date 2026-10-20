import Work.Mixed.Bridge.Core
import Work.Mixed.Bridge.Inst
import Work.Mixed.Rate3945H

/-!
# bB-data: the design-B instance (#3945H) of p-bridge's bridge (round wht45/allheads)

x-inst's `Inst2979.lean` (itself p-bridge's `Inst.lean` with other numbers) with the numbers of design B: hub PR #3945
(Dugongue's staged544 T360 terminal-gap fill, 360 replicas, 3,240 calls) with the all-heads fill of round wht45/reopt.
Target rule: `card α = 73`, heads of rank `48 = 24 + 24` cut by 24 lines (target 49), rank-24 remnants, every paid
rank at most 48 (`48 < 49`), `W0 = 322560` arrays per gauge copy (360 replicas x 896) of which `Hd = 80640` are cut,
the cut profile `cut3945H` of R5 (`Work.Mixed.Rate3945H`).

* `unit3945H`: the expected UNIT profile of the design-B network (= #3945's PROFILE.json, unchanged by the new fill;
  the 80,640 literal heads still rank 48; `unit3945H_eq : unit3945H = unit3945`), and `hcut3945H`: unit profile ->
  cut profile (checked by `decide`);
* `HFNet.bundle3945H`: the hypothesis bundle of `Eng.engine_seam_final` with kmin 49 and R5's exponent `zRate3945H`;
* `HFNet.engine3945H`: the engine theorem at the seam exponent `alphaZ = 1 - (100 N + 1)/10^12`, N = 17581520
  (e-r45m1's convention), i.e. the statement `hillsZ_program` takes.
Written by x-inst (wht45/lever2979) on p-bridge's Inst.lean; instance 3945H by bB-data (wht45/allheads, 19 Oct 2026, design B: #3945 with every child-8/10 bin a literal head); prepared with assistance from Claude.
No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.PB
open Binary Matrix Finset GF MX RAM

noncomputable section

universe U

/-- Expected unit profile of the design-B network (rank, count): #3945's profile; the 80,640 literal heads are rank 48. -/
def unit3945H : List (ℕ × ℕ) :=
  [(1, 8041680), (2, 2306880), (3, 729648), (4, 353160), (5, 69120), (8, 25920), (48, 137952)]

/-- the unit profile of design B is #3945's (only the fill changed, not the profile). -/
theorem unit3945H_eq : unit3945H = unit3945 := rfl

/-- Unit profile -> R5's cut profile for design B: 80,640 heads of rank 48 become rank-24 remnants. -/
theorem hcut3945H : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn unit3945H r
    = BlockAccounting.histFn cut3945H r + (if r = 24 + 24 then 80640 else 0) :=
  hcut_of_check unit3945H cut3945H 24 (24 + 24) 80640 [1, 2, 3, 4, 5, 8, 24, 48]
    (by decide) (by decide) (by decide) (by decide) (by decide)

/-- Negative control: the identity fails if the heads are NOT replaced (the unit profile is not the cut profile). -/
theorem unit3945H_ne_cut : BlockAccounting.histFn unit3945H 24 ≠ BlockAccounting.histFn cut3945H 24 := by decide

/-- Negative control: one head more or less breaks the profile identity (rank 48). -/
theorem hcut3945H_sharp : (if (48:ℕ) = 24 then 80641 else 0) + BlockAccounting.histFn unit3945H 48
    ≠ BlockAccounting.histFn cut3945H 48 + (if (48:ℕ) = 24 + 24 then 80641 else 0) := by decide

section Inst
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]

namespace HFNet
variable (N : HFNet α ρ σ)

/-- **(f) for design B (#3945H)**: the hypothesis bundle of `Eng.engine_seam_final` (kmin 49, rate fact at R5's `zRate3945H`). -/
theorem bundle3945H (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3945H r + (if r = 24 + 24 then 80640 else 0))
    (hLlen : N.L.length = N.G * 80640) (hσ : Fintype.card σ = N.G * 322560) (hG : 1 ≤ N.G) :
    ∃ (E : Eng α ρ σ) (R : MX.Rule) (typ : ℕ → Fin R.K),
      E.wd = N.cutWord ∧ E.e = N.e ∧ (∀ l, E.tr l = (cutTgt N.L N.A N.S l).1) ∧ E.mx = 73 ∧
      E.BridgeHyp R typ 49 ∧ R.RateQ zRate3945H 1 := by
  have hcut' : ∀ r, (if r = N.b then 80640 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3945H r + (if r = N.a + N.b then 80640 else 0) := by
    rw [ha, hb]; exact hcut
  have hmom : BlockAccounting.moment (Fintype.card α) zRate3945H cut3945H
      ≤ ((322560 - 80640 : ℕ) : ℝ) + BlockAccounting.T (Fintype.card α) zRate3945H (Fintype.card α - N.a) 80640 := by
    rw [hα, ha]
    exact rate3945H_moment
  obtain ⟨E, R, typ, h1, h2, h3, h4, hB, hR⟩ := N.cut_bundle cut3945H hcut' hLlen hσ hG (by norm_num) (by norm_num)
    48 (fun q hq => (cut3945H_proper q hq).2) (by rw [hα, ha]; norm_num) (by rw [ha]; norm_num)
    (by rw [hb]; norm_num) (by rw [hα, ha, hb]; norm_num) (by rw [hα, ha]; norm_num) hmom
  refine ⟨E, R, typ, h1, h2, h3, by rw [h4, hα], ?_, hR⟩
  have : Fintype.card α - N.a = 49 := by rw [hα, ha]
  rw [this] at hB
  exact hB

/-- **THE ENGINE THEOREM FOR DESIGN B, #3945H (half cut)** from the network: the `hillsZ_program` statement at the seam
exponent `1 - 1758152001/10^12` (DFT figure 1.7581520e-3 with e-r45m1's chain convention). -/
theorem engine3945H (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 80640 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3945H r + (if r = 24 + 24 then 80640 else 0))
    (hLlen : N.L.length = N.G * 80640) (hσ : Fintype.card σ = N.G * 322560) (hG : 1 ≤ N.G)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1758152001:ℝ)/(10:ℝ)^12) (k i)) := by
  obtain ⟨E, R, typ, -, -, -, -, hB, hR⟩ := N.bundle3945H hα ha hb hcut hLlen hσ hG
  exact E.engine_seam_final hB E.mx le_rfl hR (by norm_num [zRate3945H]) (by norm_num [zRate3945H])
    (by norm_num) m cl k v

end HFNet
end Inst

end
end PowerSaving.PB
end OAI

#print axioms OAI.PowerSaving.PB.hcut3945H
#print axioms OAI.PowerSaving.PB.unit3945H_ne_cut
#print axioms OAI.PowerSaving.PB.hcut3945H_sharp
#print axioms OAI.PowerSaving.PB.HFNet.bundle3945H
#print axioms OAI.PowerSaving.PB.HFNet.engine3945H
