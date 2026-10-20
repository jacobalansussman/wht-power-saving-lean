import Work.GFrame.Engine.Main
import Work.GFrame.Engine.Cert

/-!
# ExLegit: the engine accepts a unit that is only known to EXIST

Research file of 14 October 2026; it is not part of any record build.  It imports this repository's engine
files `Work/GFrame/Engine/{Main,Cert}.lean` (with `Letters.lean` and `Work/Block/Word.lean` behind them).  It was
first elaborated against Shea's Lean tree, whose copies of those four files are byte-identical (sha256) to this
repository's; see `README.md` beside this file.

* `exlegit_engine_of_exists`   the engine's kernel program from `∃ w, GLiveKernel e w ∧ Proper ∧ price`,
                               by `obtain`: no word is ever written down.
* `exlegit_engine_of_choose`   the same with the unit named by `Classical.choose` (a noncomputable word).
* `exlegit_costR_localM`       the price of a mixing letter does not depend on its rational matrix.
* `exlegit_engine_of_coeff_exists`   a fixed SHAPE with unknown coefficients: price proved once for the
                               shape, correctness only for some (non-explicit) coefficient vector.
* `exlegit_engine_of_gcert`    the same from the engine's own existential certificate `GCert`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI.PowerSaving
open RAM Cluster

namespace RAM
open Binary Matrix Ty Finset GF
universe U

/-- **(A)** A unit given only by an existence statement.  The conclusion is the engine's kernel program
(`engine_program_g_costR`, Work/GFrame/Engine/Main.lean), i.e. the shape of upstream `hills_program`. -/
theorem exlegit_engine_of_exists {α ρ σ : Type} [Fintype α] [DecidableEq α] [Fintype ρ]
    [DecidableEq ρ] [Fintype σ] [DecidableEq σ]
    (u a : ℕ) (hα : Fintype.card α = u) (hσ : Fintype.card σ = 2^a) (hu : 3 ≤ u)
    (e : σ → ρ) (he : Function.Injective e) (z : ℝ) (hz : 0 ≤ z)
    (hex : ∃ w : GWord α ρ, GLiveKernel e w ∧ w.Proper u ∧
        w.costR (fun r => ((r:ℝ)/(u:ℝ))^z) < (2:ℝ)^a)
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) := by
  obtain ⟨w, hw, hP, hmom⟩ := hex
  exact engine_program_g_costR u a hα hσ hu e he w hw hP z hz hmom cl m k v

/-- **(B)** The same, with the unit NAMED by `Classical.choose` (a noncomputable term). -/
theorem exlegit_engine_of_choose {α ρ σ : Type} [Fintype α] [DecidableEq α] [Fintype ρ]
    [DecidableEq ρ] [Fintype σ] [DecidableEq σ]
    (u a : ℕ) (hα : Fintype.card α = u) (hσ : Fintype.card σ = 2^a) (hu : 3 ≤ u)
    (e : σ → ρ) (he : Function.Injective e) (z : ℝ) (hz : 0 ≤ z)
    (hex : ∃ w : GWord α ρ, GLiveKernel e w ∧ w.Proper u ∧
        w.costR (fun r => ((r:ℝ)/(u:ℝ))^z) < (2:ℝ)^a)
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  engine_program_g_costR u a hα hσ hu e he (Classical.choose hex) (Classical.choose_spec hex).1
    (Classical.choose_spec hex).2.1 z hz (Classical.choose_spec hex).2.2 cl m k v

/-- **(C)** The price of a free mixing letter does not depend on its (rational) matrix: `rfl`. -/
theorem exlegit_costR_localM {α ρ : Type} (φ : ℕ → ℝ) (g g' : Matrix ρ ρ ℚ) :
    (GMove.old (BMove.localM g) : GMove α ρ).costR φ =
      (GMove.old (BMove.localM g') : GMove α ρ).costR φ := rfl

/-- **(D)** A fixed SHAPE `shape : κ → GWord α ρ` (κ = the unknown coefficients).  The price is proved
once, for the shape; correctness (`GLiveKernel`) is only known for SOME coefficient vector. -/
theorem exlegit_engine_of_coeff_exists {α ρ σ κ : Type} [Fintype α] [DecidableEq α] [Fintype ρ]
    [DecidableEq ρ] [Fintype σ] [DecidableEq σ]
    (u a : ℕ) (hα : Fintype.card α = u) (hσ : Fintype.card σ = 2^a) (hu : 3 ≤ u)
    (e : σ → ρ) (he : Function.Injective e) (z : ℝ) (hz : 0 ≤ z)
    (shape : κ → GWord α ρ) (c0 : κ)
    (hprice : ∀ c, (shape c).costR (fun r => ((r:ℝ)/(u:ℝ))^z) =
      (shape c0).costR (fun r => ((r:ℝ)/(u:ℝ))^z))
    (hP : ∀ c, (shape c).Proper u)
    (hmom : (shape c0).costR (fun r => ((r:ℝ)/(u:ℝ))^z) < (2:ℝ)^a)
    (hex : ∃ c, GLiveKernel e (shape c))
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) := by
  obtain ⟨c, hc⟩ := hex
  exact engine_program_g_costR u a hα hσ hu e he (shape c) hc (hP c) z hz
    (by rw [hprice c]; exact hmom) cl m k v

/-- **(E)** From the engine's own existential certificate `GCert` (Work/GFrame/Engine/Cert.lean:104),
whose definition is already `∃ w, …`. -/
theorem exlegit_engine_of_gcert {α ρ σ : Type} [Fintype α] [DecidableEq α] [Fintype ρ]
    [DecidableEq ρ] [Fintype σ] [DecidableEq σ]
    (u a : ℕ) (hα : Fintype.card α = u) (hσ : Fintype.card σ = 2^a) (hu : 3 ≤ u)
    (e : σ → ρ) (he : Function.Injective e) (z : ℝ) (hz : 0 ≤ z)
    (c : (ℕ → ℝ) → ℝ) (hc : GCert (α := α) e c)
    (hmom : c (fun r => ((r:ℝ)/(u:ℝ))^z) < (2:ℝ)^a)
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) := by
  obtain ⟨w, hP, hcost, hw⟩ := hc
  rw [hα] at hP
  exact engine_program_g_costR u a hα hσ hu e he w hw hP z hz (by rw [hcost]; exact hmom) cl m k v

end RAM
end OAI.PowerSaving

#print axioms OAI.PowerSaving.RAM.exlegit_engine_of_exists
#print axioms OAI.PowerSaving.RAM.exlegit_engine_of_choose
#print axioms OAI.PowerSaving.RAM.exlegit_costR_localM
#print axioms OAI.PowerSaving.RAM.exlegit_engine_of_coeff_exists
#print axioms OAI.PowerSaving.RAM.exlegit_engine_of_gcert
