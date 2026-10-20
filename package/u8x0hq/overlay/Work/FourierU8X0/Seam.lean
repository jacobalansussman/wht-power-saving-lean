import Work.GCert.Data.U8X0Glue401
import Work.GCert.Data.InstU8X0
import Work.Bridge401.End

/-!
# THE SEAM for the Fourier transfer: our kernel program in OpenAI's `hills_program` form

New file, not a copy of an OpenAI file: Jacob Sussman's Work/FourierE8/Seam.lean (Apache-2.0) for the
entrance-gauged nine-stage layout of PR #401 (Dugongue, CrocSwap/integer-mult-bounds; on jacklightChen's PR #377 layout;
Lean: gauged halves DFT Bridge3, network DFT New Sums, end theorem Work/Bridge401/End.lean), certificate `GXD.U8X0`,
suffix `U8X0`; written by tools/mkseam401.py (DFT Lean). Sections 2 and 3 of Jacob's file are replaced by one theorem
that calls `GX.engine_of_ginv_401`.  It contains lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of
lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math
(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hillsU8X0`, `alphaU8X0`.

Twin of Work/Fourier/Seam.lean for the certificate `GXD.U8X0` (label width 9; Work/GCert/Data/U8X0.lean) in the
place of `GXD.P193`: the same text with the names and numbers of that certificate, and with every declared name
changed (`U8X0` in the place of `Z`, `U8X0` appended to the other names, namespace `SeamU8X0` for `Seam`), so that the
two files have no declaration in common.  This file imports the U8X0 instance only, not the published seam; that
is why the general lemmas of sections 1 and 2 are stated here again, under new names, with the same proofs.

OpenAI's Fourier proof uses the tensor network through ONE theorem, `RAM.hills_program`
(OAI/Computability/FourierTransform/TensorProgram.lean:95), with ONE call site (`Bench.gleamT`,
SectorAlgorithm.lean:45), plus four envelope facts (`hills_pos`, `hills_mono`, `hills_le`, `hills_real`,
TensorSaving.lean:116-133) and the number facts `alpha_pos`, `alpha_lt`, `Skies.alpha_theta`,
`Skies.decimal_range`.  This file gives all of them at OUR exponent

  `alphaU8X0 = 1 - 13928619/10^10`      (OpenAI: `alpha = 1 - 2/10^11`),

under OpenAI's names with `U8X0` after `hills` / `alpha`, in OpenAI's own namespaces.

* `GX.engine_of_ginvU8X0`, `GX.engine_of_halvesU8X0`   twins of `GX.wht_of_ginv`, `GX.wht_of_halves`
  (Work/GCert/Chain/End.lean) ending in the kernel program (`RAM.engine_program_ggroup_list`) instead of the
  Walsh-Hadamard corollary: every colour `cl`, both modes `m`, every batch of inputs.
* `SeamU8X0.engine_U8X0x_of_scal`, `_of_ids`   twins of `WHT.wht_U8X0x_of_scal`, `_of_ids`
  (Work/GCert/Data/U8X0.lean): same certificate `GXD.U8X0`, same unit list, same rate fact
  `BlockAccounting.fold_U8X0x`, fill parameter `s = 40`.
* `RAM.hillsU8X0_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.U8X0.scalOK`.
* `RAM.envelope_program_geU8X0`   the same at every exponent `z ≥ alphaU8X0`.

No `sorry`, no `native_decide`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
noncomputable section

/-! ## 1. The exponent and the envelope facts -/

/-- Our exponent, in the place of OpenAI's `alpha` (TensorSaving.lean:81). -/
def alphaU8X0 : ℝ := 1 - 13928619/(10:ℝ)^10

theorem alphaU8X0_pos : 0 < alphaU8X0 := by norm_num [alphaU8X0]
theorem alphaU8X0_lt : alphaU8X0 < 1 := by norm_num [alphaU8X0]

lemma envelope_monoU8X0 (z : ℝ) (hz : 0 ≤ z) {a b : ℕ} (h : a ≤ b) : envelope z a ≤ envelope z b := by
  unfold envelope
  gcongr

-- `envelope_le` (for `z ≤ 1`) and `envelope_real` (for `0 ≤ z`) already exist: Work/Block/WHT.lean:26, :31.

lemma envelope_mono_expU8X0 {z z' : ℝ} (h : z ≤ z') (k : ℕ) : envelope z k ≤ envelope z' k := by
  unfold envelope
  exact Nat.ceil_le_ceil (Real.rpow_le_rpow_of_exponent_le (by simp) h)

/-- Natural envelope at our exponent, in the place of OpenAI's `hills` (TensorSaving.lean:114). -/
def hillsU8X0 (k : ℕ) := Nat.ceil (((k:ℝ)+1)^alphaU8X0)

lemma hillsU8X0_eq_envelope : hillsU8X0 = envelope alphaU8X0 := rfl

lemma hillsU8X0_pos (k : ℕ) : 1≤hillsU8X0 k := envelope_pos alphaU8X0 k

lemma hillsU8X0_mono {a b : ℕ} (h : a≤b) : hillsU8X0 a ≤ hillsU8X0 b := envelope_monoU8X0 alphaU8X0 alphaU8X0_pos.le h

lemma hillsU8X0_le (k : ℕ) : hillsU8X0 k ≤ k+1 := envelope_le alphaU8X0 alphaU8X0_lt.le k

lemma hillsU8X0_real (k : ℕ) : (hillsU8X0 k:ℝ) ≤ 2*((k:ℝ)+1)^alphaU8X0 := envelope_real alphaU8X0 alphaU8X0_pos.le k

/-! ## 3. The certificate `GXD.U8X0` with entrance gauges (PR #401 layout) -/

namespace SeamU8X0
open Binary Matrix Finset RAM RAM.Ty
universe U

/-- **Kernel program at `1 - 13928619/10^10`** for the certificate `GXD.U8X0`: `GX.engine_of_ginv_401` with its gauged
invocation package `GXD.U8X0.ginv` (`GXD.ginvOf` of its kernel checks), its bank schedule `GXD.U8X0.sched` (`Work/GCert/Data/U8X0Glue401.lean`), the
unit list `BlockAccounting.blocksU8X0x` and the rate fact `BlockAccounting.fold_U8X0x` (fill parameter `s = 40`). -/
theorem engine_U8X0x (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 13928619/(10:ℝ)^10) (k i)) :=
  GX.engine_of_ginv_401 GXD.U8X0.ginv GXD.U8X0.h9 GXD.U8X0.sched GXD.U8X0.fitsS
    73 10895 rfl GXD.U8X0.width (by norm_num) BlockAccounting.blocksU8X0x
    (fun φ => by rw [GXD.U8X0.price φ]; exact BridgeRate.U8X0.unit_cost_x φ)
    40 (1 - 13928619/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_U8X0x cl m k v

end SeamU8X0

/-! ## 4. THE SEAM -/

namespace RAM
open Binary Matrix Ty Finset
universe U

/-- **THE SEAM.**  OpenAI's `hills_program` (TensorProgram.lean:95) word for word, with `hillsU8X0`
(exponent `1 - 13928619/10^10`) in the place of `hills` (exponent `1 - 2/10^11`): the kernel program of
the certificate `GXD.U8X0` (the one behind `WHT.wht_main_block_U8X0x`), for every colour,
both modes and every batch of inputs, in `2^k * hillsU8X0 k` operations. -/
theorem hillsU8X0_program (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*hillsU8X0 (k i)) :=
  SeamU8X0.engine_U8X0x cl m k v

/-- **General exponent**: the same at every exponent `z ≥ alphaU8X0`, envelope `⌈(k+1)^z⌉`. -/
theorem envelope_program_geU8X0 (z : ℝ) (hz : alphaU8X0 ≤ z) (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  (hillsU8X0_program cl m k v).mono fun i => Nat.mul_le_mul_left _ (envelope_mono_expU8X0 hz (k i))

end RAM

/-! ## 5. Number facts for the copy of `Main.lean` -/

namespace SeamU8X0

lemma alphaU8X0_lt_alpha : alphaU8X0 < alpha := by norm_num [alphaU8X0, alpha]

/-- for `Skies.alpha_theta` (Main.lean:39): below the paper's exponent. -/
lemma alphaU8X0_theta : alphaU8X0 < theta := alphaU8X0_lt_alpha.trans Skies.alpha_theta

/-- for `Skies.decimal_range` (Main.lean:59). -/
lemma alphaU8X0_decimal : alphaU8X0 < decimalExponent ∧ decimalExponent < 1 := by
  norm_num [decimalExponent, alphaU8X0]

/-- Target exponent of the final Fourier statement: strictly above `alphaU8X0`, so that the
`(log log n)^2` factor of the reduction is absorbed. -/
def targetExponentU8X0 : ℝ := 1 - 13928618/(10:ℝ)^10
/-- Fallback exponent. -/
def fallbackExponentU8X0 : ℝ := 1 - 102210/(10:ℝ)^8

lemma alphaU8X0_lt_target : alphaU8X0 < targetExponentU8X0 := by norm_num [alphaU8X0, targetExponentU8X0]
lemma alphaU8X0_lt_fallback : alphaU8X0 < fallbackExponentU8X0 := by norm_num [alphaU8X0, fallbackExponentU8X0]
lemma target_lt_fallbackU8X0 : targetExponentU8X0 < fallbackExponentU8X0 := by
  norm_num [targetExponentU8X0, fallbackExponentU8X0]
lemma target_posU8X0 : 0 < targetExponentU8X0 := by norm_num [targetExponentU8X0]
lemma target_lt_oneU8X0 : targetExponentU8X0 < 1 := by norm_num [targetExponentU8X0]
lemma fallback_lt_oneU8X0 : fallbackExponentU8X0 < 1 := by norm_num [fallbackExponentU8X0]
lemma target_lt_decimalU8X0 : targetExponentU8X0 < decimalExponent := by
  norm_num [targetExponentU8X0, decimalExponent]
lemma fallback_lt_decimalU8X0 : fallbackExponentU8X0 < decimalExponent := by
  norm_num [fallbackExponentU8X0, decimalExponent]
lemma target_lt_thetaU8X0 : targetExponentU8X0 < theta :=
  (show targetExponentU8X0 < alpha by norm_num [targetExponentU8X0, alpha]).trans Skies.alpha_theta
lemma fallback_lt_thetaU8X0 : fallbackExponentU8X0 < theta :=
  (show fallbackExponentU8X0 < alpha by norm_num [fallbackExponentU8X0, alpha]).trans Skies.alpha_theta

end SeamU8X0
end
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.RAM.hillsU8X0_program
#print axioms OAI.PowerSaving.RAM.envelope_program_geU8X0
#print axioms OAI.PowerSaving.SeamU8X0.alphaU8X0_theta
#print axioms OAI.PowerSaving.hillsU8X0_real
