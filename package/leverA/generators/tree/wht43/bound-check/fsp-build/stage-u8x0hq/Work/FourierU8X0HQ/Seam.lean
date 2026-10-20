import Work.XHost.U8X0HQ

/-!
# THE SEAM for the Fourier transfer: our kernel program in OpenAI's `hills_program` form

New file, not a copy of an OpenAI file: Jacob Sussman's Work/FourierE8/Seam.lean (Apache-2.0) for the
entrance-gauged nine-stage layout of PR #401 (Dugongue, CrocSwap/integer-mult-bounds; on jacklightChen's PR #377 layout;
Lean: gauged halves DFT Bridge3, network DFT New Sums, end theorem Work/Bridge401/End.lean), certificate `GXD.U8X0HQ`,
suffix `U8X0HQ`; written by tools/mkseam401.py (DFT Lean). Sections 2 and 3 of Jacob's file are replaced by one theorem
that calls `XHost.engine_of_ginv_host` (h-L3-build: the HOSTED network).  It contains lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of
lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math
(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hillsU8X0HQ`, `alphaU8X0HQ`.

Twin of Work/Fourier/Seam.lean for the certificate `GXD.U8X0HQ` (label width 9; Work/GCert/Data/U8X0HQ.lean) in the
place of `GXD.P193`: the same text with the names and numbers of that certificate, and with every declared name
changed (`U8X0HQ` in the place of `Z`, `U8X0HQ` appended to the other names, namespace `SeamU8X0HQ` for `Seam`), so that the
two files have no declaration in common.  This file imports the U8X0HQ instance only, not the published seam; that
is why the general lemmas of sections 1 and 2 are stated here again, under new names, with the same proofs.

OpenAI's Fourier proof uses the tensor network through ONE theorem, `RAM.hills_program`
(OAI/Computability/FourierTransform/TensorProgram.lean:95), with ONE call site (`Bench.gleamT`,
SectorAlgorithm.lean:45), plus four envelope facts (`hills_pos`, `hills_mono`, `hills_le`, `hills_real`,
TensorSaving.lean:116-133) and the number facts `alpha_pos`, `alpha_lt`, `Skies.alpha_theta`,
`Skies.decimal_range`.  This file gives all of them at OUR exponent

  `alphaU8X0HQ = 1 - 16256630/10^10`      (OpenAI: `alpha = 1 - 2/10^11`),

under OpenAI's names with `U8X0HQ` after `hills` / `alpha`, in OpenAI's own namespaces.

* `GX.engine_of_ginvU8X0HQ`, `GX.engine_of_halvesU8X0HQ`   twins of `GX.wht_of_ginv`, `GX.wht_of_halves`
  (Work/GCert/Chain/End.lean) ending in the kernel program (`RAM.engine_program_ggroup_list`) instead of the
  Walsh-Hadamard corollary: every colour `cl`, both modes `m`, every batch of inputs.
* `SeamU8X0HQ.engine_U8X0HQx_of_scal`, `_of_ids`   twins of `WHT.wht_U8X0HQx_of_scal`, `_of_ids`
  (Work/GCert/Data/U8X0HQ.lean): same certificate `GXD.U8X0HQ`, same unit list, same rate fact
  `BlockAccounting.fold_U8X0x`, fill parameter `s = 40`.
* `RAM.hillsU8X0HQ_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.U8X0.scalOK`.
* `RAM.envelope_program_geU8X0HQ`   the same at every exponent `z ≥ alphaU8X0HQ`.

No `sorry`, no `native_decide`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
noncomputable section

/-! ## 1. The exponent and the envelope facts -/

/-- Our exponent, in the place of OpenAI's `alpha` (TensorSaving.lean:81). -/
def alphaU8X0HQ : ℝ := 1 - 16256630/(10:ℝ)^10

theorem alphaU8X0HQ_pos : 0 < alphaU8X0HQ := by norm_num [alphaU8X0HQ]
theorem alphaU8X0HQ_lt : alphaU8X0HQ < 1 := by norm_num [alphaU8X0HQ]

lemma envelope_monoU8X0HQ (z : ℝ) (hz : 0 ≤ z) {a b : ℕ} (h : a ≤ b) : envelope z a ≤ envelope z b := by
  unfold envelope
  gcongr

-- `envelope_le` (for `z ≤ 1`) and `envelope_real` (for `0 ≤ z`) already exist: Work/Block/WHT.lean:26, :31.

lemma envelope_mono_expU8X0HQ {z z' : ℝ} (h : z ≤ z') (k : ℕ) : envelope z k ≤ envelope z' k := by
  unfold envelope
  exact Nat.ceil_le_ceil (Real.rpow_le_rpow_of_exponent_le (by simp) h)

/-- Natural envelope at our exponent, in the place of OpenAI's `hills` (TensorSaving.lean:114). -/
def hillsU8X0HQ (k : ℕ) := Nat.ceil (((k:ℝ)+1)^alphaU8X0HQ)

lemma hillsU8X0HQ_eq_envelope : hillsU8X0HQ = envelope alphaU8X0HQ := rfl

lemma hillsU8X0HQ_pos (k : ℕ) : 1≤hillsU8X0HQ k := envelope_pos alphaU8X0HQ k

lemma hillsU8X0HQ_mono {a b : ℕ} (h : a≤b) : hillsU8X0HQ a ≤ hillsU8X0HQ b := envelope_monoU8X0HQ alphaU8X0HQ alphaU8X0HQ_pos.le h

lemma hillsU8X0HQ_le (k : ℕ) : hillsU8X0HQ k ≤ k+1 := envelope_le alphaU8X0HQ alphaU8X0HQ_lt.le k

lemma hillsU8X0HQ_real (k : ℕ) : (hillsU8X0HQ k:ℝ) ≤ 2*((k:ℝ)+1)^alphaU8X0HQ := envelope_real alphaU8X0HQ alphaU8X0HQ_pos.le k

/-! ## 3. The certificate `GXD.U8X0HQ` with entrance gauges (PR #401 layout) -/

namespace SeamU8X0HQ
open Binary Matrix Finset RAM RAM.Ty
universe U

/-- **Kernel program at `1 - 16256630/10^10`** for the certificate `GXD.U8X0HQ`: `GX.engine_of_ginv_401` with its gauged
invocation package `GXD.U8X0.ginv` (`GXD.ginvOf` of its kernel checks), its bank schedule `GXD.U8X0.sched` (`Work/GCert/Data/U8X0HQGlue401.lean`), the
unit list `BlockAccounting.blocksU8X0x` and the rate fact `BlockAccounting.fold_U8X0x` (fill parameter `s = 40`). -/
theorem engine_U8X0HQx (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 16256630/(10:ℝ)^10) (k i)) :=
  XHost.engine_of_ginv_host GXD.U8X0.ginv GXD.U8X0.h9 GXD.U8X0HQ.hsched GXD.U8X0HQ.hfits
    73 7491 rfl GXD.U8X0HQ.hwidth (by norm_num) BlockAccounting.blocksU8X0HQx
    (fun φ => by rw [GXD.U8X0HQ.hprice φ]; exact BridgeRate.U8X0HQ.unit_cost_x φ)
    40 (1 - 16256630/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_U8X0HQx cl m k v

end SeamU8X0HQ

/-! ## 4. THE SEAM -/

namespace RAM
open Binary Matrix Ty Finset
universe U

/-- **THE SEAM.**  OpenAI's `hills_program` (TensorProgram.lean:95) word for word, with `hillsU8X0HQ`
(exponent `1 - 16256630/10^10`) in the place of `hills` (exponent `1 - 2/10^11`): the kernel program of
the certificate `GXD.U8X0HQ` (the one behind `WHT.wht_main_block_U8X0HQx`), for every colour,
both modes and every batch of inputs, in `2^k * hillsU8X0HQ k` operations. -/
theorem hillsU8X0HQ_program (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*hillsU8X0HQ (k i)) :=
  SeamU8X0HQ.engine_U8X0HQx cl m k v

/-- **General exponent**: the same at every exponent `z ≥ alphaU8X0HQ`, envelope `⌈(k+1)^z⌉`. -/
theorem envelope_program_geU8X0HQ (z : ℝ) (hz : alphaU8X0HQ ≤ z) (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  (hillsU8X0HQ_program cl m k v).mono fun i => Nat.mul_le_mul_left _ (envelope_mono_expU8X0HQ hz (k i))

end RAM

/-! ## 5. Number facts for the copy of `Main.lean` -/

namespace SeamU8X0HQ

lemma alphaU8X0HQ_lt_alpha : alphaU8X0HQ < alpha := by norm_num [alphaU8X0HQ, alpha]

/-- for `Skies.alpha_theta` (Main.lean:39): below the paper's exponent. -/
lemma alphaU8X0HQ_theta : alphaU8X0HQ < theta := alphaU8X0HQ_lt_alpha.trans Skies.alpha_theta

/-- for `Skies.decimal_range` (Main.lean:59). -/
lemma alphaU8X0HQ_decimal : alphaU8X0HQ < decimalExponent ∧ decimalExponent < 1 := by
  norm_num [decimalExponent, alphaU8X0HQ]

/-- Target exponent of the final Fourier statement: strictly above `alphaU8X0HQ`, so that the
`(log log n)^2` factor of the reduction is absorbed. -/
def targetExponentU8X0HQ : ℝ := 1 - 16256629/(10:ℝ)^10
/-- Fallback exponent. -/
def fallbackExponentU8X0HQ : ℝ := 1 - 102210/(10:ℝ)^8

lemma alphaU8X0HQ_lt_target : alphaU8X0HQ < targetExponentU8X0HQ := by norm_num [alphaU8X0HQ, targetExponentU8X0HQ]
lemma alphaU8X0HQ_lt_fallback : alphaU8X0HQ < fallbackExponentU8X0HQ := by norm_num [alphaU8X0HQ, fallbackExponentU8X0HQ]
lemma target_lt_fallbackU8X0HQ : targetExponentU8X0HQ < fallbackExponentU8X0HQ := by
  norm_num [targetExponentU8X0HQ, fallbackExponentU8X0HQ]
lemma target_posU8X0HQ : 0 < targetExponentU8X0HQ := by norm_num [targetExponentU8X0HQ]
lemma target_lt_oneU8X0HQ : targetExponentU8X0HQ < 1 := by norm_num [targetExponentU8X0HQ]
lemma fallback_lt_oneU8X0HQ : fallbackExponentU8X0HQ < 1 := by norm_num [fallbackExponentU8X0HQ]
lemma target_lt_decimalU8X0HQ : targetExponentU8X0HQ < decimalExponent := by
  norm_num [targetExponentU8X0HQ, decimalExponent]
lemma fallback_lt_decimalU8X0HQ : fallbackExponentU8X0HQ < decimalExponent := by
  norm_num [fallbackExponentU8X0HQ, decimalExponent]
lemma target_lt_thetaU8X0HQ : targetExponentU8X0HQ < theta :=
  (show targetExponentU8X0HQ < alpha by norm_num [targetExponentU8X0HQ, alpha]).trans Skies.alpha_theta
lemma fallback_lt_thetaU8X0HQ : fallbackExponentU8X0HQ < theta :=
  (show fallbackExponentU8X0HQ < alpha by norm_num [fallbackExponentU8X0HQ, alpha]).trans Skies.alpha_theta

end SeamU8X0HQ
end
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.RAM.hillsU8X0HQ_program
#print axioms OAI.PowerSaving.RAM.envelope_program_geU8X0HQ
#print axioms OAI.PowerSaving.SeamU8X0HQ.alphaU8X0HQ_theta
#print axioms OAI.PowerSaving.hillsU8X0HQ_real
