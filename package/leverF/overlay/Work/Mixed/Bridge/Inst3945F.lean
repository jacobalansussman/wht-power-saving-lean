import Work.Mixed.Bridge.CoreF
import Work.Mixed.Bridge.Inst3945H
import Work.Mixed.Rate3945F

/-!
# fc-engine: the design-B FULL-cut instance (#3945F) of the bridge (round wht45/fullcut)

Design B (hub PR #3945, Dugongue's staged544 T360 terminal-gap fill, with the every-head-kept fill of round
wht45/reopt): the unit profile is #3945's, `unit3945H` (bB-data), with 80,640 literal heads of rank 48 per gauge copy.
The FULL cut: every head is removed entirely (its 48 lines; target `73 - 48 = 25` on the head array), and every other
rank-48 block (57,312 per copy) is split `24 + 24`.  Cut profile `cut3945F` (fc-rate's R5 file `Work.Mixed.Rate3945F`):
{1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 24: 114624}; every paid rank at most 24 < 25
(separated); `W0 = 322560` arrays per copy, `Hd = 80640` of them cut.

* `hcutF3945F`: unit profile -> full-cut profile (checked by `decide`), with negatives `hcutF3945F_sharp` (one head
  more breaks it) and `unit3945H_ne_cutF` (the unit profile is not the cut profile);
* `engine3945F_of_cut`: the engine theorem from p-hfirst's heads-first conclusion with the full-cut data given
  directly (frames `A`, the 48 head lines `S`, `bph (z l) = bph ((A l).fam (S l))`), at the seam exponent
  `alphaZ = 1 - (100 N + 1)/10^12`, N = 17806935 (e-r45m1's convention), i.e. the statement `hillsZ_program` takes;
* `bundle3945F_of_cut`: the hypothesis bundle of `Eng.engine_seam_final` (kmin 25, rate fact at `zRate3945F`).
fc-engine (wht45/fullcut, 19-20 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.PB
open Binary Matrix Finset GF MX RAM

noncomputable section

universe U

/-- Unit profile -> full-cut profile for design B: the 80,640 heads of rank 48 disappear, the other 57,312 rank-48
blocks are split into 2 x 57,312 = 114,624 blocks of rank 24. -/
theorem hcutF3945F : ∀ r, BlockAccounting.histFn cut3945F r
    = (if r = 24 + 24 then 0 else BlockAccounting.histFn unit3945H r)
      + (if r = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 80640 else 0)
      + (if r = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 80640 else 0) :=
  hcutF_of_check unit3945H cut3945F 24 24 80640 [1, 2, 3, 4, 5, 8, 24, 48]
    (by decide) (by decide) (by decide) (by decide) (by decide)

/-- Negative control: the identity fails with one head more (80,641). -/
theorem hcutF3945F_sharp : BlockAccounting.histFn cut3945F 24
    ≠ (if (24:ℕ) = 24 + 24 then 0 else BlockAccounting.histFn unit3945H 24)
      + (if (24:ℕ) = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 80641 else 0)
      + (if (24:ℕ) = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 80641 else 0) := by decide

/-- Negative control: the unit profile is not the full-cut profile. -/
theorem unit3945H_ne_cutF : BlockAccounting.histFn unit3945H 48 ≠ BlockAccounting.histFn cut3945F 48 := by decide

section Inst
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]

/-- **(f) for design B, full cut (#3945F)**: the hypothesis bundle of `Eng.engine_seam_final` (kmin 25, rate fact at
R5's `zRate3945F`), from p-hfirst's heads-first conclusion and the full-cut data. -/
theorem bundle3945F_of_cut (e : σ → ρ) (he : Function.Injective e) (L : List σ) (hL : L.Nodup)
    (z : σ → Fin 48 → Space α) (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (S : σ → Finset α)
    (hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x)
    (hS : ∀ l ∈ L, (S l).card = 48)
    (G : ℕ) (c0 : (ℕ → ℝ) → ℝ) (U : List (ℕ × ℕ)) (hU : ∀ φ, c0 φ = listCost φ U)
    (hnet : ∃ body : GWord α ρ, body.Proper (Fintype.card α) ∧
      (∀ φ, body.costR φ = (G : ℝ) * c0 φ - (L.length : ℝ) * φ 48) ∧
      GLiveKernel e (headWord e L (fun _ => 48) z ind ++ body))
    (hα : Fintype.card α = 73)
    (hcutF : ∀ r, BlockAccounting.histFn cut3945F r
      = (if r = 24 + 24 then 0 else BlockAccounting.histFn U r)
        + (if r = 24 then BlockAccounting.histFn U (24 + 24) - 80640 else 0)
        + (if r = 24 then BlockAccounting.histFn U (24 + 24) - 80640 else 0))
    (hLlen : L.length = G * 80640) (hσ : Fintype.card σ = G * 322560) (hG : 1 ≤ G) :
    ∃ (E : Eng α ρ σ) (R : MX.Rule) (typ : ℕ → Fin R.K),
      E.e = e ∧ E.mx = 73 ∧ E.BridgeHyp R typ 25 ∧ R.RateQ zRate3945F 1 := by
  obtain ⟨body, hP, hcost, hk⟩ := hnet
  have hcost' : ∀ φ, body.costR φ = (G : ℝ) * listCost φ U - (L.length : ℝ) * φ 48 := fun φ => by
    rw [hcost φ, hU φ]
  have hmom : BlockAccounting.moment (Fintype.card α) zRate3945F cut3945F
      ≤ ((322560 - 80640 : ℕ) : ℝ) + BlockAccounting.T (Fintype.card α) zRate3945F (Fintype.card α - 48) 80640 := by
    rw [hα]
    exact rate3945F_moment
  obtain ⟨E, R, typ, -, h2, -, h4, hB, hR⟩ :=
    full_bundle e he L hL (fun _ => 48) z ind A S hz hS body hP hk 24 24 (by norm_num) (by norm_num)
      (by rw [hα]; norm_num) (by rw [hα]; norm_num) cut3945F hcost' hcutF hLlen hσ hG (by norm_num) (by norm_num)
      24 (fun q hq => (cut3945F_proper q hq).2) (by rw [hα]; norm_num) (by norm_num) (by rw [hα]; norm_num)
      (by rw [hα]; norm_num) hmom
  refine ⟨E, R, typ, h2, by rw [h4, hα], ?_, hR⟩
  have : Fintype.card α - 48 = 25 := by rw [hα]
  rw [this] at hB
  exact hB

/-- **THE ENGINE THEOREM FOR DESIGN B, FULL CUT (#3945F)** from p-hfirst's heads-first conclusion with the full-cut
data given directly: the `hillsZ_program` statement at the seam exponent `1 - 1780693501/10^12` (DFT figure
1.7806935e-3 with e-r45m1's chain convention). -/
theorem engine3945F_of_cut (e : σ → ρ) (he : Function.Injective e) (L : List σ) (hL : L.Nodup)
    (z : σ → Fin 48 → Space α) (ind : ∀ l, LinearIndependent F (z l)) (A : σ → OBase α) (S : σ → Finset α)
    (hz : ∀ l ∈ L, ∀ x, bph (z l) x = bph ((A l).fam (S l)) x)
    (hS : ∀ l ∈ L, (S l).card = 48)
    (G : ℕ) (c0 : (ℕ → ℝ) → ℝ) (U : List (ℕ × ℕ)) (hU : ∀ φ, c0 φ = listCost φ U)
    (hnet : ∃ body : GWord α ρ, body.Proper (Fintype.card α) ∧
      (∀ φ, body.costR φ = (G : ℝ) * c0 φ - (L.length : ℝ) * φ 48) ∧
      GLiveKernel e (headWord e L (fun _ => 48) z ind ++ body))
    (hα : Fintype.card α = 73)
    (hcutF : ∀ r, BlockAccounting.histFn cut3945F r
      = (if r = 24 + 24 then 0 else BlockAccounting.histFn U r)
        + (if r = 24 then BlockAccounting.histFn U (24 + 24) - 80640 else 0)
        + (if r = 24 then BlockAccounting.histFn U (24 + 24) - 80640 else 0))
    (hLlen : L.length = G * 80640) (hσ : Fintype.card σ = G * 322560) (hG : 1 ≤ G)
    (m : Bool) (cl : Paint) {ι : Type U} (k : ι → ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2 ^ k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2 ^ k i * envelope (1 - (1780693501:ℝ)/(10:ℝ)^12) (k i)) := by
  obtain ⟨E, R, typ, -, -, hB, hR⟩ :=
    bundle3945F_of_cut e he L hL z ind A S hz hS G c0 U hU hnet hα hcutF hLlen hσ hG
  exact E.engine_seam_final hB E.mx le_rfl hR (by norm_num [zRate3945F]) (by norm_num [zRate3945F])
    (by norm_num) m cl k v

end Inst

end
end PowerSaving.PB
end OAI

#print axioms OAI.PowerSaving.PB.hcutF3945F
#print axioms OAI.PowerSaving.PB.hcutF3945F_sharp
#print axioms OAI.PowerSaving.PB.unit3945H_ne_cutF
#print axioms OAI.PowerSaving.PB.bundle3945F_of_cut
#print axioms OAI.PowerSaving.PB.engine3945F_of_cut
