import Work.GCert.Data.C588MGlue401
import Work.GCert.Data.InstC588M
import Work.Bridge401.End

/-!
# THE SEAM for the Fourier transfer: our kernel program in OpenAI's `hills_program` form

New file, not a copy of an OpenAI file: Jacob Sussman's Work/FourierE8/Seam.lean (Apache-2.0) for the
entrance-gauged nine-stage layout of PR #401 (Dugongue, CrocSwap/integer-mult-bounds; on jacklightChen's PR #377 layout;
Lean: gauged halves DFT Bridge3, network DFT New Sums, end theorem Work/Bridge401/End.lean), certificate `GXD.C588M`,
suffix `C588M`; written by tools/mkseam401.py (DFT Lean). Sections 2 and 3 of Jacob's file are replaced by one theorem
that calls `GX.engine_of_ginv_401`.  It contains lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of
lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math
(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hillsC588M`, `alphaC588M`.

Twin of Work/Fourier/Seam.lean for the certificate `GXD.C588M` (label width 9; Work/GCert/Data/C588M.lean) in the
place of `GXD.P193`: the same text with the names and numbers of that certificate, and with every declared name
changed (`C588M` in the place of `Z`, `C588M` appended to the other names, namespace `SeamC588M` for `Seam`), so that the
two files have no declaration in common.  This file imports the C588M instance only, not the published seam; that
is why the general lemmas of sections 1 and 2 are stated here again, under new names, with the same proofs.

OpenAI's Fourier proof uses the tensor network through ONE theorem, `RAM.hills_program`
(OAI/Computability/FourierTransform/TensorProgram.lean:95), with ONE call site (`Bench.gleamT`,
SectorAlgorithm.lean:45), plus four envelope facts (`hills_pos`, `hills_mono`, `hills_le`, `hills_real`,
TensorSaving.lean:116-133) and the number facts `alpha_pos`, `alpha_lt`, `Skies.alpha_theta`,
`Skies.decimal_range`.  This file gives all of them at OUR exponent

  `alphaC588M = 1 - 13566303/10^10`      (OpenAI: `alpha = 1 - 2/10^11`),

under OpenAI's names with `C588M` after `hills` / `alpha`, in OpenAI's own namespaces.

* `GX.engine_of_ginvC588M`, `GX.engine_of_halvesC588M`   twins of `GX.wht_of_ginv`, `GX.wht_of_halves`
  (Work/GCert/Chain/End.lean) ending in the kernel program (`RAM.engine_program_ggroup_list`) instead of the
  Walsh-Hadamard corollary: every colour `cl`, both modes `m`, every batch of inputs.
* `SeamC588M.engine_C588Mx_of_scal`, `_of_ids`   twins of `WHT.wht_C588Mx_of_scal`, `_of_ids`
  (Work/GCert/Data/C588M.lean): same certificate `GXD.C588M`, same unit list, same rate fact
  `BlockAccounting.fold_C588Mx`, fill parameter `s = 40`.
* `RAM.hillsC588M_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.C588M.scalOK`.
* `RAM.envelope_program_geC588M`   the same at every exponent `z ≥ alphaC588M`.

No `sorry`, no `native_decide`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
noncomputable section

/-! ## 1. The exponent and the envelope facts -/

/-- Our exponent, in the place of OpenAI's `alpha` (TensorSaving.lean:81). -/
def alphaC588M : ℝ := 1 - 13566303/(10:ℝ)^10

theorem alphaC588M_pos : 0 < alphaC588M := by norm_num [alphaC588M]
theorem alphaC588M_lt : alphaC588M < 1 := by norm_num [alphaC588M]

lemma envelope_monoC588M (z : ℝ) (hz : 0 ≤ z) {a b : ℕ} (h : a ≤ b) : envelope z a ≤ envelope z b := by
  unfold envelope
  gcongr

-- `envelope_le` (for `z ≤ 1`) and `envelope_real` (for `0 ≤ z`) already exist: Work/Block/WHT.lean:26, :31.

lemma envelope_mono_expC588M {z z' : ℝ} (h : z ≤ z') (k : ℕ) : envelope z k ≤ envelope z' k := by
  unfold envelope
  exact Nat.ceil_le_ceil (Real.rpow_le_rpow_of_exponent_le (by simp) h)

/-- Natural envelope at our exponent, in the place of OpenAI's `hills` (TensorSaving.lean:114). -/
def hillsC588M (k : ℕ) := Nat.ceil (((k:ℝ)+1)^alphaC588M)

lemma hillsC588M_eq_envelope : hillsC588M = envelope alphaC588M := rfl

lemma hillsC588M_pos (k : ℕ) : 1≤hillsC588M k := envelope_pos alphaC588M k

lemma hillsC588M_mono {a b : ℕ} (h : a≤b) : hillsC588M a ≤ hillsC588M b := envelope_monoC588M alphaC588M alphaC588M_pos.le h

lemma hillsC588M_le (k : ℕ) : hillsC588M k ≤ k+1 := envelope_le alphaC588M alphaC588M_lt.le k

lemma hillsC588M_real (k : ℕ) : (hillsC588M k:ℝ) ≤ 2*((k:ℝ)+1)^alphaC588M := envelope_real alphaC588M alphaC588M_pos.le k

/-! ## 3. The certificate `GXD.C588M` with entrance gauges (PR #401 layout) -/

namespace SeamC588M
open Binary Matrix Finset RAM RAM.Ty
universe U

/-- **Kernel program at `1 - 13566303/10^10`** for the certificate `GXD.C588M`: `GX.engine_of_ginv_401` with its gauged
invocation package `GXD.C588M.ginv` (`GXD.ginvOf` of its kernel checks), its bank schedule `GXD.C588M.sched` (`Work/GCert/Data/C588MGlue401.lean`), the
unit list `BlockAccounting.blocksC588Mx` and the rate fact `BlockAccounting.fold_C588Mx` (fill parameter `s = 40`). -/
theorem engine_C588Mx (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 13566303/(10:ℝ)^10) (k i)) :=
  GX.engine_of_ginv_401 GXD.C588M.ginv GXD.C588M.h9 GXD.C588M.sched GXD.C588M.fitsS
    73 11173 rfl GXD.C588M.width (by norm_num) BlockAccounting.blocksC588Mx
    (fun φ => by rw [GXD.C588M.price φ]; exact BridgeRate.C588M.unit_cost_x φ)
    40 (1 - 13566303/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_C588Mx cl m k v

end SeamC588M

/-! ## 4. THE SEAM -/

namespace RAM
open Binary Matrix Ty Finset
universe U

/-- **THE SEAM.**  OpenAI's `hills_program` (TensorProgram.lean:95) word for word, with `hillsC588M`
(exponent `1 - 13566303/10^10`) in the place of `hills` (exponent `1 - 2/10^11`): the kernel program of
the certificate `GXD.C588M` (the one behind `WHT.wht_main_block_C588Mx`), for every colour,
both modes and every batch of inputs, in `2^k * hillsC588M k` operations. -/
theorem hillsC588M_program (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*hillsC588M (k i)) :=
  SeamC588M.engine_C588Mx cl m k v

/-- **General exponent**: the same at every exponent `z ≥ alphaC588M`, envelope `⌈(k+1)^z⌉`. -/
theorem envelope_program_geC588M (z : ℝ) (hz : alphaC588M ≤ z) (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  (hillsC588M_program cl m k v).mono fun i => Nat.mul_le_mul_left _ (envelope_mono_expC588M hz (k i))

end RAM

/-! ## 5. Number facts for the copy of `Main.lean` -/

namespace SeamC588M

lemma alphaC588M_lt_alpha : alphaC588M < alpha := by norm_num [alphaC588M, alpha]

/-- for `Skies.alpha_theta` (Main.lean:39): below the paper's exponent. -/
lemma alphaC588M_theta : alphaC588M < theta := alphaC588M_lt_alpha.trans Skies.alpha_theta

/-- for `Skies.decimal_range` (Main.lean:59). -/
lemma alphaC588M_decimal : alphaC588M < decimalExponent ∧ decimalExponent < 1 := by
  norm_num [decimalExponent, alphaC588M]

/-- Target exponent of the final Fourier statement: strictly above `alphaC588M`, so that the
`(log log n)^2` factor of the reduction is absorbed. -/
def targetExponentC588M : ℝ := 1 - 13566302/(10:ℝ)^10
/-- Fallback exponent. -/
def fallbackExponentC588M : ℝ := 1 - 102210/(10:ℝ)^8

lemma alphaC588M_lt_target : alphaC588M < targetExponentC588M := by norm_num [alphaC588M, targetExponentC588M]
lemma alphaC588M_lt_fallback : alphaC588M < fallbackExponentC588M := by norm_num [alphaC588M, fallbackExponentC588M]
lemma target_lt_fallbackC588M : targetExponentC588M < fallbackExponentC588M := by
  norm_num [targetExponentC588M, fallbackExponentC588M]
lemma target_posC588M : 0 < targetExponentC588M := by norm_num [targetExponentC588M]
lemma target_lt_oneC588M : targetExponentC588M < 1 := by norm_num [targetExponentC588M]
lemma fallback_lt_oneC588M : fallbackExponentC588M < 1 := by norm_num [fallbackExponentC588M]
lemma target_lt_decimalC588M : targetExponentC588M < decimalExponent := by
  norm_num [targetExponentC588M, decimalExponent]
lemma fallback_lt_decimalC588M : fallbackExponentC588M < decimalExponent := by
  norm_num [fallbackExponentC588M, decimalExponent]
lemma target_lt_thetaC588M : targetExponentC588M < theta :=
  (show targetExponentC588M < alpha by norm_num [targetExponentC588M, alpha]).trans Skies.alpha_theta
lemma fallback_lt_thetaC588M : fallbackExponentC588M < theta :=
  (show fallbackExponentC588M < alpha by norm_num [fallbackExponentC588M, alpha]).trans Skies.alpha_theta

end SeamC588M
end
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.RAM.hillsC588M_program
#print axioms OAI.PowerSaving.RAM.envelope_program_geC588M
#print axioms OAI.PowerSaving.SeamC588M.alphaC588M_theta
#print axioms OAI.PowerSaving.hillsC588M_real
