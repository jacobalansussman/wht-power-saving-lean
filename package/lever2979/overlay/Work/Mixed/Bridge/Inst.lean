import Work.Mixed.Bridge.Core
import Work.Mixed.Rate3945
import Work.Mixed.Rate3815

/-!
# p-bridge: the #3945 (and #3815) instance of the bridge (round wht45/lever)

Target rule (g-plan PLAN.md section 6): `card α = 73`, heads of rank `48 = 24 + 24` cut by 24 lines (target 49),
rank-24 remnants, every paid rank at most 48 (`48 < 49`), `W0 = 322560` arrays per gauge copy of which
`Hd = 57312` are cut, the cut profile `cut3945` of R5 (e-r45m1's `Work.Mixed.Rate3945`).

* `unit3945` / `unit3815`: the expected UNIT profile of the network (lever-reprice: 3,240 calls + idle blocks;
  the 57,312 heads still rank 48), and `hcut3945` / `hcut3815`: unit profile -> cut profile (checked by `decide`);
* `HFNet.bundle3945`: the hypothesis bundle of `Eng.engine_seam_final` with kmin 49 and R5's exponent `zRate3945`;
* `HFNet.engine3945`: the engine theorem at the seam exponent `alphaZ = 1 - (100 N + 1)/10^12`, N = 17501109
  (e-r45m1's convention), i.e. the statement `hillsZ_program` takes;
* the same for #3815 (N = 17500242).
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.PB
open Binary Matrix Finset GF MX RAM

noncomputable section

universe U

/-- Expected unit profile of the #3945 network (rank, count): the 57,312 literal heads are rank 48. -/
def unit3945 : List (ℕ × ℕ) :=
  [(1, 8041680), (2, 2306880), (3, 729648), (4, 353160), (5, 69120), (8, 25920), (48, 137952)]

/-- Expected unit profile of the #3815 network. -/
def unit3815 : List (ℕ × ℕ) :=
  [(1, 7989840), (2, 2349000), (3, 749088), (4, 330480), (5, 69120), (8, 25920), (48, 137952)]

/-- Unit profile -> R5's cut profile for #3945: 57,312 heads of rank 48 become rank-24 remnants. -/
theorem hcut3945 : ∀ r, (if r = 24 then 57312 else 0) + BlockAccounting.histFn unit3945 r
    = BlockAccounting.histFn cut3945 r + (if r = 24 + 24 then 57312 else 0) :=
  hcut_of_check unit3945 cut3945 24 (24 + 24) 57312 [1, 2, 3, 4, 5, 8, 24, 48]
    (by decide) (by decide) (by decide) (by decide) (by decide)

theorem hcut3815 : ∀ r, (if r = 24 then 57312 else 0) + BlockAccounting.histFn unit3815 r
    = BlockAccounting.histFn cut3815 r + (if r = 24 + 24 then 57312 else 0) :=
  hcut_of_check unit3815 cut3815 24 (24 + 24) 57312 [1, 2, 3, 4, 5, 8, 24, 48]
    (by decide) (by decide) (by decide) (by decide) (by decide)

/-- Negative control: the identity fails if the heads are NOT replaced (the unit profile is not the cut profile). -/
theorem unit3945_ne_cut : BlockAccounting.histFn unit3945 24 ≠ BlockAccounting.histFn cut3945 24 := by decide

section Inst
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]

namespace HFNet
variable (N : HFNet α ρ σ)

/-- **(f) for #3945**: the hypothesis bundle of `Eng.engine_seam_final` (kmin 49, rate fact at R5's `zRate3945`). -/
theorem bundle3945 (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 57312 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3945 r + (if r = 24 + 24 then 57312 else 0))
    (hLlen : N.L.length = N.G * 57312) (hσ : Fintype.card σ = N.G * 322560) (hG : 1 ≤ N.G) :
    ∃ (E : Eng α ρ σ) (R : MX.Rule) (typ : ℕ → Fin R.K),
      E.wd = N.cutWord ∧ E.e = N.e ∧ (∀ l, E.tr l = (cutTgt N.L N.A N.S l).1) ∧ E.mx = 73 ∧
      E.BridgeHyp R typ 49 ∧ R.RateQ zRate3945 1 := by
  have hcut' : ∀ r, (if r = N.b then 57312 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3945 r + (if r = N.a + N.b then 57312 else 0) := by
    rw [ha, hb]; exact hcut
  have hmom : BlockAccounting.moment (Fintype.card α) zRate3945 cut3945
      ≤ ((322560 - 57312 : ℕ) : ℝ) + BlockAccounting.T (Fintype.card α) zRate3945 (Fintype.card α - N.a) 57312 := by
    rw [hα, ha]
    exact rate3945_moment
  obtain ⟨E, R, typ, h1, h2, h3, h4, hB, hR⟩ := N.cut_bundle cut3945 hcut' hLlen hσ hG (by norm_num) (by norm_num)
    48 (fun q hq => (cut3945_proper q hq).2) (by rw [hα, ha]; norm_num) (by rw [ha]; norm_num)
    (by rw [hb]; norm_num) (by rw [hα, ha, hb]; norm_num) (by rw [hα, ha]; norm_num) hmom
  refine ⟨E, R, typ, h1, h2, h3, by rw [h4, hα], ?_, hR⟩
  have : Fintype.card α - N.a = 49 := by rw [hα, ha]
  rw [this] at hB
  exact hB

/-- **THE ENGINE THEOREM FOR #3945 (half cut)** from the network: the `hillsZ_program` statement at the seam
exponent `1 - 1750110901/10^12` (DFT figure 1.7501109e-3 with e-r45m1's chain convention). -/
theorem engine3945 (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 57312 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3945 r + (if r = 24 + 24 then 57312 else 0))
    (hLlen : N.L.length = N.G * 57312) (hσ : Fintype.card σ = N.G * 322560) (hG : 1 ≤ N.G)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1750110901:ℝ)/(10:ℝ)^12) (k i)) := by
  obtain ⟨E, R, typ, -, -, -, -, hB, hR⟩ := N.bundle3945 hα ha hb hcut hLlen hσ hG
  exact E.engine_seam_final hB E.mx le_rfl hR (by norm_num [zRate3945]) (by norm_num [zRate3945])
    (by norm_num) m cl k v

/-- (f) for #3815 (fallback): kmin 49, rate fact at R5's `zRate3815`. -/
theorem bundle3815 (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 57312 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3815 r + (if r = 24 + 24 then 57312 else 0))
    (hLlen : N.L.length = N.G * 57312) (hσ : Fintype.card σ = N.G * 322560) (hG : 1 ≤ N.G) :
    ∃ (E : Eng α ρ σ) (R : MX.Rule) (typ : ℕ → Fin R.K),
      E.wd = N.cutWord ∧ E.e = N.e ∧ (∀ l, E.tr l = (cutTgt N.L N.A N.S l).1) ∧ E.mx = 73 ∧
      E.BridgeHyp R typ 49 ∧ R.RateQ zRate3815 1 := by
  have hcut' : ∀ r, (if r = N.b then 57312 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3815 r + (if r = N.a + N.b then 57312 else 0) := by
    rw [ha, hb]; exact hcut
  have hmom : BlockAccounting.moment (Fintype.card α) zRate3815 cut3815
      ≤ ((322560 - 57312 : ℕ) : ℝ) + BlockAccounting.T (Fintype.card α) zRate3815 (Fintype.card α - N.a) 57312 := by
    rw [hα, ha]
    exact rate3815_moment
  obtain ⟨E, R, typ, h1, h2, h3, h4, hB, hR⟩ := N.cut_bundle cut3815 hcut' hLlen hσ hG (by norm_num) (by norm_num)
    48 (fun q hq => (cut3815_proper q hq).2) (by rw [hα, ha]; norm_num) (by rw [ha]; norm_num)
    (by rw [hb]; norm_num) (by rw [hα, ha, hb]; norm_num) (by rw [hα, ha]; norm_num) hmom
  refine ⟨E, R, typ, h1, h2, h3, by rw [h4, hα], ?_, hR⟩
  have : Fintype.card α - N.a = 49 := by rw [hα, ha]
  rw [this] at hB
  exact hB

/-- The engine theorem for #3815 at `1 - 1750024201/10^12` (DFT figure 1.7500242e-3). -/
theorem engine3815 (hα : Fintype.card α = 73) (ha : N.a = 24) (hb : N.b = 24)
    (hcut : ∀ r, (if r = 24 then 57312 else 0) + BlockAccounting.histFn N.U r
      = BlockAccounting.histFn cut3815 r + (if r = 24 + 24 then 57312 else 0))
    (hLlen : N.L.length = N.G * 57312) (hσ : Fintype.card σ = N.G * 322560) (hG : 1 ≤ N.G)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1750024201:ℝ)/(10:ℝ)^12) (k i)) := by
  obtain ⟨E, R, typ, -, -, -, -, hB, hR⟩ := N.bundle3815 hα ha hb hcut hLlen hσ hG
  exact E.engine_seam_final hB E.mx le_rfl hR (by norm_num [zRate3815]) (by norm_num [zRate3815])
    (by norm_num) m cl k v

end HFNet
end Inst

end
end PowerSaving.PB
end OAI

#print axioms OAI.PowerSaving.PB.hcut3945
#print axioms OAI.PowerSaving.PB.hcut3815
#print axioms OAI.PowerSaving.PB.unit3945_ne_cut
#print axioms OAI.PowerSaving.PB.HFNet.bundle3945
#print axioms OAI.PowerSaving.PB.HFNet.engine3945
#print axioms OAI.PowerSaving.PB.HFNet.bundle3815
#print axioms OAI.PowerSaving.PB.HFNet.engine3815
