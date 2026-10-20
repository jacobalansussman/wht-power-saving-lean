import Work.XHost.C588MH

/-!
# THE SEAM for the Fourier transfer: our kernel program in OpenAI's `hills_program` form

New file, not a copy of an OpenAI file: Jacob Sussman's Work/FourierE8/Seam.lean (Apache-2.0) for the
entrance-gauged nine-stage layout of PR #401 (Dugongue, CrocSwap/integer-mult-bounds; on jacklightChen's PR #377 layout;
Lean: gauged halves DFT Bridge3, network DFT New Sums, end theorem Work/Bridge401/End.lean), certificate `GXD.C588MH`,
suffix `C588MH`; written by tools/mkseam401.py (DFT Lean). Sections 2 and 3 of Jacob's file are replaced by one theorem
that calls `XHost.engine_of_ginv_host` (h-L3-build: the HOSTED network).  It contains lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of
lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math
(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hillsC588MH`, `alphaC588MH`.

Twin of Work/Fourier/Seam.lean for the certificate `GXD.C588MH` (label width 9; Work/GCert/Data/C588MH.lean) in the
place of `GXD.P193`: the same text with the names and numbers of that certificate, and with every declared name
changed (`C588MH` in the place of `Z`, `C588MH` appended to the other names, namespace `SeamC588MH` for `Seam`), so that the
two files have no declaration in common.  This file imports the C588MH instance only, not the published seam; that
is why the general lemmas of sections 1 and 2 are stated here again, under new names, with the same proofs.

OpenAI's Fourier proof uses the tensor network through ONE theorem, `RAM.hills_program`
(OAI/Computability/FourierTransform/TensorProgram.lean:95), with ONE call site (`Bench.gleamT`,
SectorAlgorithm.lean:45), plus four envelope facts (`hills_pos`, `hills_mono`, `hills_le`, `hills_real`,
TensorSaving.lean:116-133) and the number facts `alpha_pos`, `alpha_lt`, `Skies.alpha_theta`,
`Skies.decimal_range`.  This file gives all of them at OUR exponent

  `alphaC588MH = 1 - 13552607/10^10`      (OpenAI: `alpha = 1 - 2/10^11`),

under OpenAI's names with `C588MH` after `hills` / `alpha`, in OpenAI's own namespaces.

* `GX.engine_of_ginvC588MH`, `GX.engine_of_halvesC588MH`   twins of `GX.wht_of_ginv`, `GX.wht_of_halves`
  (Work/GCert/Chain/End.lean) ending in the kernel program (`RAM.engine_program_ggroup_list`) instead of the
  Walsh-Hadamard corollary: every colour `cl`, both modes `m`, every batch of inputs.
* `SeamC588MH.engine_C588MHx_of_scal`, `_of_ids`   twins of `WHT.wht_C588MHx_of_scal`, `_of_ids`
  (Work/GCert/Data/C588MH.lean): same certificate `GXD.C588MH`, same unit list, same rate fact
  `BlockAccounting.fold_C588Mx`, fill parameter `s = 40`.
* `RAM.hillsC588MH_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.C588M.scalOK`.
* `RAM.envelope_program_geC588MH`   the same at every exponent `z ≥ alphaC588MH`.

No `sorry`, no `native_decide`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
noncomputable section

/-! ## 1. The exponent and the envelope facts -/

/-- Our exponent, in the place of OpenAI's `alpha` (TensorSaving.lean:81). -/
def alphaC588MH : ℝ := 1 - 13552607/(10:ℝ)^10

theorem alphaC588MH_pos : 0 < alphaC588MH := by norm_num [alphaC588MH]
theorem alphaC588MH_lt : alphaC588MH < 1 := by norm_num [alphaC588MH]

lemma envelope_monoC588MH (z : ℝ) (hz : 0 ≤ z) {a b : ℕ} (h : a ≤ b) : envelope z a ≤ envelope z b := by
  unfold envelope
  gcongr

-- `envelope_le` (for `z ≤ 1`) and `envelope_real` (for `0 ≤ z`) already exist: Work/Block/WHT.lean:26, :31.

lemma envelope_mono_expC588MH {z z' : ℝ} (h : z ≤ z') (k : ℕ) : envelope z k ≤ envelope z' k := by
  unfold envelope
  exact Nat.ceil_le_ceil (Real.rpow_le_rpow_of_exponent_le (by simp) h)

/-- Natural envelope at our exponent, in the place of OpenAI's `hills` (TensorSaving.lean:114). -/
def hillsC588MH (k : ℕ) := Nat.ceil (((k:ℝ)+1)^alphaC588MH)

lemma hillsC588MH_eq_envelope : hillsC588MH = envelope alphaC588MH := rfl

lemma hillsC588MH_pos (k : ℕ) : 1≤hillsC588MH k := envelope_pos alphaC588MH k

lemma hillsC588MH_mono {a b : ℕ} (h : a≤b) : hillsC588MH a ≤ hillsC588MH b := envelope_monoC588MH alphaC588MH alphaC588MH_pos.le h

lemma hillsC588MH_le (k : ℕ) : hillsC588MH k ≤ k+1 := envelope_le alphaC588MH alphaC588MH_lt.le k

lemma hillsC588MH_real (k : ℕ) : (hillsC588MH k:ℝ) ≤ 2*((k:ℝ)+1)^alphaC588MH := envelope_real alphaC588MH alphaC588MH_pos.le k

/-! ## 3. The certificate `GXD.C588MH` with entrance gauges (PR #401 layout) -/

namespace SeamC588MH
open Binary Matrix Finset RAM RAM.Ty
universe U

/-- **Kernel program at `1 - 13552607/10^10`** for the certificate `GXD.C588MH`: `GX.engine_of_ginv_401` with its gauged
invocation package `GXD.C588M.ginv` (`GXD.ginvOf` of its kernel checks), its bank schedule `GXD.C588M.sched` (`Work/GCert/Data/C588MHGlue401.lean`), the
unit list `BlockAccounting.blocksC588Mx` and the rate fact `BlockAccounting.fold_C588Mx` (fill parameter `s = 40`). -/
theorem engine_C588MHx (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 13552607/(10:ℝ)^10) (k i)) :=
  XHost.engine_of_ginv_host GXD.C588M.ginv GXD.C588M.h9 GXD.C588MH.hsched GXD.C588MH.hfits
    73 11200 rfl GXD.C588MH.hwidth (by norm_num) BlockAccounting.blocksC588Mx
    (fun φ => by rw [GXD.C588MH.hprice φ]; exact BridgeRate.C588M.unit_cost_x φ)
    40 (1 - 13552607/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_C588Mx cl m k v

end SeamC588MH

/-! ## 4. THE SEAM -/

namespace RAM
open Binary Matrix Ty Finset
universe U

/-- **THE SEAM.**  OpenAI's `hills_program` (TensorProgram.lean:95) word for word, with `hillsC588MH`
(exponent `1 - 13552607/10^10`) in the place of `hills` (exponent `1 - 2/10^11`): the kernel program of
the certificate `GXD.C588MH` (the one behind `WHT.wht_main_block_C588MHx`), for every colour,
both modes and every batch of inputs, in `2^k * hillsC588MH k` operations. -/
theorem hillsC588MH_program (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*hillsC588MH (k i)) :=
  SeamC588MH.engine_C588MHx cl m k v

/-- **General exponent**: the same at every exponent `z ≥ alphaC588MH`, envelope `⌈(k+1)^z⌉`. -/
theorem envelope_program_geC588MH (z : ℝ) (hz : alphaC588MH ≤ z) (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  (hillsC588MH_program cl m k v).mono fun i => Nat.mul_le_mul_left _ (envelope_mono_expC588MH hz (k i))

end RAM

/-! ## 5. Number facts for the copy of `Main.lean` -/

namespace SeamC588MH

lemma alphaC588MH_lt_alpha : alphaC588MH < alpha := by norm_num [alphaC588MH, alpha]

/-- for `Skies.alpha_theta` (Main.lean:39): below the paper's exponent. -/
lemma alphaC588MH_theta : alphaC588MH < theta := alphaC588MH_lt_alpha.trans Skies.alpha_theta

/-- for `Skies.decimal_range` (Main.lean:59). -/
lemma alphaC588MH_decimal : alphaC588MH < decimalExponent ∧ decimalExponent < 1 := by
  norm_num [decimalExponent, alphaC588MH]

/-- Target exponent of the final Fourier statement: strictly above `alphaC588MH`, so that the
`(log log n)^2` factor of the reduction is absorbed. -/
def targetExponentC588MH : ℝ := 1 - 13552606/(10:ℝ)^10
/-- Fallback exponent. -/
def fallbackExponentC588MH : ℝ := 1 - 102210/(10:ℝ)^8

lemma alphaC588MH_lt_target : alphaC588MH < targetExponentC588MH := by norm_num [alphaC588MH, targetExponentC588MH]
lemma alphaC588MH_lt_fallback : alphaC588MH < fallbackExponentC588MH := by norm_num [alphaC588MH, fallbackExponentC588MH]
lemma target_lt_fallbackC588MH : targetExponentC588MH < fallbackExponentC588MH := by
  norm_num [targetExponentC588MH, fallbackExponentC588MH]
lemma target_posC588MH : 0 < targetExponentC588MH := by norm_num [targetExponentC588MH]
lemma target_lt_oneC588MH : targetExponentC588MH < 1 := by norm_num [targetExponentC588MH]
lemma fallback_lt_oneC588MH : fallbackExponentC588MH < 1 := by norm_num [fallbackExponentC588MH]
lemma target_lt_decimalC588MH : targetExponentC588MH < decimalExponent := by
  norm_num [targetExponentC588MH, decimalExponent]
lemma fallback_lt_decimalC588MH : fallbackExponentC588MH < decimalExponent := by
  norm_num [fallbackExponentC588MH, decimalExponent]
lemma target_lt_thetaC588MH : targetExponentC588MH < theta :=
  (show targetExponentC588MH < alpha by norm_num [targetExponentC588MH, alpha]).trans Skies.alpha_theta
lemma fallback_lt_thetaC588MH : fallbackExponentC588MH < theta :=
  (show fallbackExponentC588MH < alpha by norm_num [fallbackExponentC588MH, alpha]).trans Skies.alpha_theta

end SeamC588MH
end
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.RAM.hillsC588MH_program
#print axioms OAI.PowerSaving.RAM.envelope_program_geC588MH
#print axioms OAI.PowerSaving.SeamC588MH.alphaC588MH_theta
#print axioms OAI.PowerSaving.hillsC588MH_real
