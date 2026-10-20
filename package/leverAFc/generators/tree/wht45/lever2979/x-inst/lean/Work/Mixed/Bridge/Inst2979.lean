import Work.Mixed.Bridge.Core
import Work.Mixed.Rate2979

/-!
# x-inst: the #2979 instance of p-bridge's bridge (round wht45/lever2979)

p-bridge's `Inst.lean` (#3945) with the numbers of hub PR #2979 (mdumitrean, on Jeremy Kothe's #2957 and Mason Hale's
#2917 words; T74, 74 replicas, 666 calls).  Target rule: `card α = 73`, heads of rank `48 = 24 + 24` cut by 24 lines
(target 49), rank-24 remnants, every paid rank at most 48 (`48 < 49`), `W0 = 66304` arrays per gauge copy
(74 replicas x 896) of which `Hd = 16086` are cut, the cut profile `cut2979` of R5 (`Work.Mixed.Rate2979`).

* `unit2979`: the expected UNIT profile of the #2979 network (p-alt's measured full profile = the hub's
  `allocation/PROFILE.json`; the 16,086 literal heads still rank 48), and `hcut2979`: unit profile -> cut profile
  (checked by `decide`);
* `HFNet.bundle2979`: the hypothesis bundle of `Eng.engine_seam_final` with kmin 49 and R5's exponent `zRate2979`;
* `HFNet.engine2979`: the engine theorem at the seam exponent `alphaZ = 1 - (100 N + 1)/10^12`, N = 17556473
  (e-r45m1's convention), i.e. the statement `hillsZ_program` takes.
Written by x-inst (wht45/lever2979, 19 Oct 2026) on p-bridge's Inst.lean; prepared with assistance from Claude.
No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.PB
open Binary Matrix Finset GF MX RAM

noncomputable section

universe U

/-- Expected unit profile of the #2979 network (rank, count): the 16,086 literal heads are rank 48. -/
def unit2979 : List (ℕ × ℕ) :=
  [(1, 1645680), (2, 470196), (3, 159468), (4, 71928), (5, 13716), (8, 5328), (48, 28190)]

/-- Unit profile -> R5's cut profile for #2979: 16,086 heads of rank 48 become rank-24 remnants. -/
theorem hcut2979 : ∀ r, (if r = 24 then 16086 else 0) + BlockAccounting.histFn unit2979 r
    = BlockAccounting.histFn cut2979 r + (if r = 24 + 24 then 16086 else 0) :=
  hcut_of_check unit2979 cut2979 24 (24 + 24) 16086 [1, 2, 3, 4, 5, 8, 24, 48]
    (by decide) (by decide) (by decide) (by decide) (by decide)

/-- Negative control: the identity fails if the heads are NOT replaced (the unit profile is not the cut profile). -/
theorem unit2979_ne_cut : BlockAccounting.histFn unit2979 24 ≠ BlockAccounting.histFn cut2979 24 := by decide

/-- Negative control: one head more or less breaks the profile identity (rank 48). -/
theorem hcut2979_sharp : (if (48:ℕ) = 24 then 16087 else 0) + BlockAccounting.histFn unit2979 48
    ≠ BlockAccounting.histFn cut2979 48 + (if (48:ℕ) = 24 + 24 then 16087 else 0) := by decide

section Inst
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]

namespace HFNet
variable (N : HFNet α ρ σ)

/-- **(f) for #2979**: the hypothesis bundle of `Eng.engine_seam_final` (kmin 49, rate fact at R5's `zRate2979`). -/
theorem bundle2979 (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 16086 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut2979 r + (if r = 24 + 24 then 16086 else 0))
    (hLlen : N.L.length = N.G * 16086) (hσ : Fintype.card σ = N.G * 66304) (hG : 1 ≤ N.G) :
    ∃ (E : Eng α ρ σ) (R : MX.Rule) (typ : ℕ → Fin R.K),
      E.wd = N.cutWord ∧ E.e = N.e ∧ (∀ l, E.tr l = (cutTgt N.L N.A N.S l).1) ∧ E.mx = 73 ∧
      E.BridgeHyp R typ 49 ∧ R.RateQ zRate2979 1 := by
  have hcut' : ∀ r, (if r = N.b then 16086 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut2979 r + (if r = N.a + N.b then 16086 else 0) := by
    rw [ha, hb]; exact hcut
  have hmom : BlockAccounting.moment (Fintype.card α) zRate2979 cut2979
      ≤ ((66304 - 16086 : ℕ) : ℝ) + BlockAccounting.T (Fintype.card α) zRate2979 (Fintype.card α - N.a) 16086 := by
    rw [hα, ha]
    exact rate2979_moment
  obtain ⟨E, R, typ, h1, h2, h3, h4, hB, hR⟩ := N.cut_bundle cut2979 hcut' hLlen hσ hG (by norm_num) (by norm_num)
    48 (fun q hq => (cut2979_proper q hq).2) (by rw [hα, ha]; norm_num) (by rw [ha]; norm_num)
    (by rw [hb]; norm_num) (by rw [hα, ha, hb]; norm_num) (by rw [hα, ha]; norm_num) hmom
  refine ⟨E, R, typ, h1, h2, h3, by rw [h4, hα], ?_, hR⟩
  have : Fintype.card α - N.a = 49 := by rw [hα, ha]
  rw [this] at hB
  exact hB

/-- **THE ENGINE THEOREM FOR #2979 (half cut)** from the network: the `hillsZ_program` statement at the seam
exponent `1 - 1755647301/10^12` (DFT figure 1.7556473e-3 with e-r45m1's chain convention). -/
theorem engine2979 (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 16086 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut2979 r + (if r = 24 + 24 then 16086 else 0))
    (hLlen : N.L.length = N.G * 16086) (hσ : Fintype.card σ = N.G * 66304) (hG : 1 ≤ N.G)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1755647301:ℝ)/(10:ℝ)^12) (k i)) := by
  obtain ⟨E, R, typ, -, -, -, -, hB, hR⟩ := N.bundle2979 hα ha hb hcut hLlen hσ hG
  exact E.engine_seam_final hB E.mx le_rfl hR (by norm_num [zRate2979]) (by norm_num [zRate2979])
    (by norm_num) m cl k v

end HFNet
end Inst

end
end PowerSaving.PB
end OAI

#print axioms OAI.PowerSaving.PB.hcut2979
#print axioms OAI.PowerSaving.PB.unit2979_ne_cut
#print axioms OAI.PowerSaving.PB.hcut2979_sharp
#print axioms OAI.PowerSaving.PB.HFNet.bundle2979
#print axioms OAI.PowerSaving.PB.HFNet.engine2979
