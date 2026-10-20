import Work.XHost.C588MHH

/-!
# THE SEAM for the Fourier transfer: our kernel program in OpenAI's `hills_program` form

New file, not a copy of an OpenAI file: Jacob Sussman's Work/FourierE8/Seam.lean (Apache-2.0) for the
entrance-gauged nine-stage layout of PR #401 (Dugongue, CrocSwap/integer-mult-bounds; on jacklightChen's PR #377 layout;
Lean: gauged halves DFT Bridge3, network DFT New Sums, end theorem Work/Bridge401/End.lean), certificate `GXD.C588MHH`,
suffix `C588MHH`; written by tools/mkseam401.py (DFT Lean). Sections 2 and 3 of Jacob's file are replaced by one theorem
that calls `XHost.engine_of_ginv_host` (h-L3-build: the HOSTED network).  It contains lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of
lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math
(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hillsC588MHH`, `alphaC588MHH`.

Twin of Work/Fourier/Seam.lean for the certificate `GXD.C588MHH` (label width 9; Work/GCert/Data/C588MHH.lean) in the
place of `GXD.P193`: the same text with the names and numbers of that certificate, and with every declared name
changed (`C588MHH` in the place of `Z`, `C588MHH` appended to the other names, namespace `SeamC588MHH` for `Seam`), so that the
two files have no declaration in common.  This file imports the C588MHH instance only, not the published seam; that
is why the general lemmas of sections 1 and 2 are stated here again, under new names, with the same proofs.

OpenAI's Fourier proof uses the tensor network through ONE theorem, `RAM.hills_program`
(OAI/Computability/FourierTransform/TensorProgram.lean:95), with ONE call site (`Bench.gleamT`,
SectorAlgorithm.lean:45), plus four envelope facts (`hills_pos`, `hills_mono`, `hills_le`, `hills_real`,
TensorSaving.lean:116-133) and the number facts `alpha_pos`, `alpha_lt`, `Skies.alpha_theta`,
`Skies.decimal_range`.  This file gives all of them at OUR exponent

  `alphaC588MHH = 1 - 15699306/10^10`      (OpenAI: `alpha = 1 - 2/10^11`),

under OpenAI's names with `C588MHH` after `hills` / `alpha`, in OpenAI's own namespaces.

* `GX.engine_of_ginvC588MHH`, `GX.engine_of_halvesC588MHH`   twins of `GX.wht_of_ginv`, `GX.wht_of_halves`
  (Work/GCert/Chain/End.lean) ending in the kernel program (`RAM.engine_program_ggroup_list`) instead of the
  Walsh-Hadamard corollary: every colour `cl`, both modes `m`, every batch of inputs.
* `SeamC588MHH.engine_C588MHHx_of_scal`, `_of_ids`   twins of `WHT.wht_C588MHHx_of_scal`, `_of_ids`
  (Work/GCert/Data/C588MHH.lean): same certificate `GXD.C588MHH`, same unit list, same rate fact
  `BlockAccounting.fold_C588Mx`, fill parameter `s = 40`.
* `RAM.hillsC588MHH_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.C588M.scalOK`.
* `RAM.envelope_program_geC588MHH`   the same at every exponent `z ≥ alphaC588MHH`.

No `sorry`, no `native_decide`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
noncomputable section

/-! ## 1. The exponent and the envelope facts -/

/-- Our exponent, in the place of OpenAI's `alpha` (TensorSaving.lean:81). -/
def alphaC588MHH : ℝ := 1 - 15699306/(10:ℝ)^10

theorem alphaC588MHH_pos : 0 < alphaC588MHH := by norm_num [alphaC588MHH]
theorem alphaC588MHH_lt : alphaC588MHH < 1 := by norm_num [alphaC588MHH]

lemma envelope_monoC588MHH (z : ℝ) (hz : 0 ≤ z) {a b : ℕ} (h : a ≤ b) : envelope z a ≤ envelope z b := by
  unfold envelope
  gcongr

-- `envelope_le` (for `z ≤ 1`) and `envelope_real` (for `0 ≤ z`) already exist: Work/Block/WHT.lean:26, :31.

lemma envelope_mono_expC588MHH {z z' : ℝ} (h : z ≤ z') (k : ℕ) : envelope z k ≤ envelope z' k := by
  unfold envelope
  exact Nat.ceil_le_ceil (Real.rpow_le_rpow_of_exponent_le (by simp) h)

/-- Natural envelope at our exponent, in the place of OpenAI's `hills` (TensorSaving.lean:114). -/
def hillsC588MHH (k : ℕ) := Nat.ceil (((k:ℝ)+1)^alphaC588MHH)

lemma hillsC588MHH_eq_envelope : hillsC588MHH = envelope alphaC588MHH := rfl

lemma hillsC588MHH_pos (k : ℕ) : 1≤hillsC588MHH k := envelope_pos alphaC588MHH k

lemma hillsC588MHH_mono {a b : ℕ} (h : a≤b) : hillsC588MHH a ≤ hillsC588MHH b := envelope_monoC588MHH alphaC588MHH alphaC588MHH_pos.le h

lemma hillsC588MHH_le (k : ℕ) : hillsC588MHH k ≤ k+1 := envelope_le alphaC588MHH alphaC588MHH_lt.le k

lemma hillsC588MHH_real (k : ℕ) : (hillsC588MHH k:ℝ) ≤ 2*((k:ℝ)+1)^alphaC588MHH := envelope_real alphaC588MHH alphaC588MHH_pos.le k

/-! ## 3. The certificate `GXD.C588MHH` with entrance gauges (PR #401 layout) -/

namespace SeamC588MHH
open Binary Matrix Finset RAM RAM.Ty
universe U

/-- **Kernel program at `1 - 15699306/10^10`** for the certificate `GXD.C588MHH`: `GX.engine_of_ginv_401` with its gauged
invocation package `GXD.C588M.ginv` (`GXD.ginvOf` of its kernel checks), its bank schedule `GXD.C588M.sched` (`Work/GCert/Data/C588MHHGlue401.lean`), the
unit list `BlockAccounting.blocksC588Mx` and the rate fact `BlockAccounting.fold_C588Mx` (fill parameter `s = 40`). -/
theorem engine_C588MHHx (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 15699306/(10:ℝ)^10) (k i)) :=
  XHost.engine_of_ginv_host GXD.C588M.ginv GXD.C588M.h9 GXD.C588MHH.hsched GXD.C588MHH.hfits
    73 7279 rfl GXD.C588MHH.hwidth (by norm_num) BlockAccounting.blocksC588MHHx
    (fun φ => by rw [GXD.C588MHH.hprice φ]; exact BridgeRate.C588MHH.unit_cost_x φ)
    40 (1 - 15699306/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_C588MHHx cl m k v

end SeamC588MHH

/-! ## 4. THE SEAM -/

namespace RAM
open Binary Matrix Ty Finset
universe U

/-- **THE SEAM.**  OpenAI's `hills_program` (TensorProgram.lean:95) word for word, with `hillsC588MHH`
(exponent `1 - 15699306/10^10`) in the place of `hills` (exponent `1 - 2/10^11`): the kernel program of
the certificate `GXD.C588MHH` (the one behind `WHT.wht_main_block_C588MHHx`), for every colour,
both modes and every batch of inputs, in `2^k * hillsC588MHH k` operations. -/
theorem hillsC588MHH_program (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*hillsC588MHH (k i)) :=
  SeamC588MHH.engine_C588MHHx cl m k v

/-- **General exponent**: the same at every exponent `z ≥ alphaC588MHH`, envelope `⌈(k+1)^z⌉`. -/
theorem envelope_program_geC588MHH (z : ℝ) (hz : alphaC588MHH ≤ z) (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  (hillsC588MHH_program cl m k v).mono fun i => Nat.mul_le_mul_left _ (envelope_mono_expC588MHH hz (k i))

end RAM

/-! ## 5. Number facts for the copy of `Main.lean` -/

namespace SeamC588MHH

lemma alphaC588MHH_lt_alpha : alphaC588MHH < alpha := by norm_num [alphaC588MHH, alpha]

/-- for `Skies.alpha_theta` (Main.lean:39): below the paper's exponent. -/
lemma alphaC588MHH_theta : alphaC588MHH < theta := alphaC588MHH_lt_alpha.trans Skies.alpha_theta

/-- for `Skies.decimal_range` (Main.lean:59). -/
lemma alphaC588MHH_decimal : alphaC588MHH < decimalExponent ∧ decimalExponent < 1 := by
  norm_num [decimalExponent, alphaC588MHH]

/-- Target exponent of the final Fourier statement: strictly above `alphaC588MHH`, so that the
`(log log n)^2` factor of the reduction is absorbed. -/
def targetExponentC588MHH : ℝ := 1 - 15699305/(10:ℝ)^10
/-- Fallback exponent. -/
def fallbackExponentC588MHH : ℝ := 1 - 102210/(10:ℝ)^8

lemma alphaC588MHH_lt_target : alphaC588MHH < targetExponentC588MHH := by norm_num [alphaC588MHH, targetExponentC588MHH]
lemma alphaC588MHH_lt_fallback : alphaC588MHH < fallbackExponentC588MHH := by norm_num [alphaC588MHH, fallbackExponentC588MHH]
lemma target_lt_fallbackC588MHH : targetExponentC588MHH < fallbackExponentC588MHH := by
  norm_num [targetExponentC588MHH, fallbackExponentC588MHH]
lemma target_posC588MHH : 0 < targetExponentC588MHH := by norm_num [targetExponentC588MHH]
lemma target_lt_oneC588MHH : targetExponentC588MHH < 1 := by norm_num [targetExponentC588MHH]
lemma fallback_lt_oneC588MHH : fallbackExponentC588MHH < 1 := by norm_num [fallbackExponentC588MHH]
lemma target_lt_decimalC588MHH : targetExponentC588MHH < decimalExponent := by
  norm_num [targetExponentC588MHH, decimalExponent]
lemma fallback_lt_decimalC588MHH : fallbackExponentC588MHH < decimalExponent := by
  norm_num [fallbackExponentC588MHH, decimalExponent]
lemma target_lt_thetaC588MHH : targetExponentC588MHH < theta :=
  (show targetExponentC588MHH < alpha by norm_num [targetExponentC588MHH, alpha]).trans Skies.alpha_theta
lemma fallback_lt_thetaC588MHH : fallbackExponentC588MHH < theta :=
  (show fallbackExponentC588MHH < alpha by norm_num [fallbackExponentC588MHH, alpha]).trans Skies.alpha_theta

end SeamC588MHH
end
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.RAM.hillsC588MHH_program
#print axioms OAI.PowerSaving.RAM.envelope_program_geC588MHH
#print axioms OAI.PowerSaving.SeamC588MHH.alphaC588MHH_theta
#print axioms OAI.PowerSaving.hillsC588MHH_real
