import Work.XHost.U8XMHH

/-!
# THE SEAM for the Fourier transfer: our kernel program in OpenAI's `hills_program` form

New file, not a copy of an OpenAI file: Jacob Sussman's Work/FourierE8/Seam.lean (Apache-2.0) for the
entrance-gauged nine-stage layout of PR #401 (Dugongue, CrocSwap/integer-mult-bounds; on jacklightChen's PR #377 layout;
Lean: gauged halves DFT Bridge3, network DFT New Sums, end theorem Work/Bridge401/End.lean), certificate `GXD.U8XMHH`,
suffix `U8XMHH`; written by tools/mkseam401.py (DFT Lean). Sections 2 and 3 of Jacob's file are replaced by one theorem
that calls `XHost.engine_of_ginv_host` (h-L3-build: the HOSTED network).  It contains lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of
lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math
(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hillsU8XMHH`, `alphaU8XMHH`.

Twin of Work/Fourier/Seam.lean for the certificate `GXD.U8XMHH` (label width 9; Work/GCert/Data/U8XMHH.lean) in the
place of `GXD.P193`: the same text with the names and numbers of that certificate, and with every declared name
changed (`U8XMHH` in the place of `Z`, `U8XMHH` appended to the other names, namespace `SeamU8XMHH` for `Seam`), so that the
two files have no declaration in common.  This file imports the U8XMHH instance only, not the published seam; that
is why the general lemmas of sections 1 and 2 are stated here again, under new names, with the same proofs.

OpenAI's Fourier proof uses the tensor network through ONE theorem, `RAM.hills_program`
(OAI/Computability/FourierTransform/TensorProgram.lean:95), with ONE call site (`Bench.gleamT`,
SectorAlgorithm.lean:45), plus four envelope facts (`hills_pos`, `hills_mono`, `hills_le`, `hills_real`,
TensorSaving.lean:116-133) and the number facts `alpha_pos`, `alpha_lt`, `Skies.alpha_theta`,
`Skies.decimal_range`.  This file gives all of them at OUR exponent

  `alphaU8XMHH = 1 - 16197706/10^10`      (OpenAI: `alpha = 1 - 2/10^11`),

under OpenAI's names with `U8XMHH` after `hills` / `alpha`, in OpenAI's own namespaces.

* `GX.engine_of_ginvU8XMHH`, `GX.engine_of_halvesU8XMHH`   twins of `GX.wht_of_ginv`, `GX.wht_of_halves`
  (Work/GCert/Chain/End.lean) ending in the kernel program (`RAM.engine_program_ggroup_list`) instead of the
  Walsh-Hadamard corollary: every colour `cl`, both modes `m`, every batch of inputs.
* `SeamU8XMHH.engine_U8XMHHx_of_scal`, `_of_ids`   twins of `WHT.wht_U8XMHHx_of_scal`, `_of_ids`
  (Work/GCert/Data/U8XMHH.lean): same certificate `GXD.U8XMHH`, same unit list, same rate fact
  `BlockAccounting.fold_U8XMx`, fill parameter `s = 40`.
* `RAM.hillsU8XMHH_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.U8XM.scalOK`.
* `RAM.envelope_program_geU8XMHH`   the same at every exponent `z ≥ alphaU8XMHH`.

No `sorry`, no `native_decide`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
noncomputable section

/-! ## 1. The exponent and the envelope facts -/

/-- Our exponent, in the place of OpenAI's `alpha` (TensorSaving.lean:81). -/
def alphaU8XMHH : ℝ := 1 - 16197706/(10:ℝ)^10

theorem alphaU8XMHH_pos : 0 < alphaU8XMHH := by norm_num [alphaU8XMHH]
theorem alphaU8XMHH_lt : alphaU8XMHH < 1 := by norm_num [alphaU8XMHH]

lemma envelope_monoU8XMHH (z : ℝ) (hz : 0 ≤ z) {a b : ℕ} (h : a ≤ b) : envelope z a ≤ envelope z b := by
  unfold envelope
  gcongr

-- `envelope_le` (for `z ≤ 1`) and `envelope_real` (for `0 ≤ z`) already exist: Work/Block/WHT.lean:26, :31.

lemma envelope_mono_expU8XMHH {z z' : ℝ} (h : z ≤ z') (k : ℕ) : envelope z k ≤ envelope z' k := by
  unfold envelope
  exact Nat.ceil_le_ceil (Real.rpow_le_rpow_of_exponent_le (by simp) h)

/-- Natural envelope at our exponent, in the place of OpenAI's `hills` (TensorSaving.lean:114). -/
def hillsU8XMHH (k : ℕ) := Nat.ceil (((k:ℝ)+1)^alphaU8XMHH)

lemma hillsU8XMHH_eq_envelope : hillsU8XMHH = envelope alphaU8XMHH := rfl

lemma hillsU8XMHH_pos (k : ℕ) : 1≤hillsU8XMHH k := envelope_pos alphaU8XMHH k

lemma hillsU8XMHH_mono {a b : ℕ} (h : a≤b) : hillsU8XMHH a ≤ hillsU8XMHH b := envelope_monoU8XMHH alphaU8XMHH alphaU8XMHH_pos.le h

lemma hillsU8XMHH_le (k : ℕ) : hillsU8XMHH k ≤ k+1 := envelope_le alphaU8XMHH alphaU8XMHH_lt.le k

lemma hillsU8XMHH_real (k : ℕ) : (hillsU8XMHH k:ℝ) ≤ 2*((k:ℝ)+1)^alphaU8XMHH := envelope_real alphaU8XMHH alphaU8XMHH_pos.le k

/-! ## 3. The certificate `GXD.U8XMHH` with entrance gauges (PR #401 layout) -/

namespace SeamU8XMHH
open Binary Matrix Finset RAM RAM.Ty
universe U

/-- **Kernel program at `1 - 16197706/10^10`** for the certificate `GXD.U8XMHH`: `GX.engine_of_ginv_401` with its gauged
invocation package `GXD.U8XM.ginv` (`GXD.ginvOf` of its kernel checks), its bank schedule `GXD.U8XM.sched` (`Work/GCert/Data/U8XMHHGlue401.lean`), the
unit list `BlockAccounting.blocksU8XMx` and the rate fact `BlockAccounting.fold_U8XMx` (fill parameter `s = 40`). -/
theorem engine_U8XMHHx (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 16197706/(10:ℝ)^10) (k i)) :=
  XHost.engine_of_ginv_host GXD.U8XM.ginv GXD.U8XM.h9 GXD.U8XMHH.hsched GXD.U8XMHH.hfits
    73 7478 rfl GXD.U8XMHH.hwidth (by norm_num) BlockAccounting.blocksU8XMHHx
    (fun φ => by rw [GXD.U8XMHH.hprice φ]; exact BridgeRate.U8XMHH.unit_cost_x φ)
    40 (1 - 16197706/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_U8XMHHx cl m k v

end SeamU8XMHH

/-! ## 4. THE SEAM -/

namespace RAM
open Binary Matrix Ty Finset
universe U

/-- **THE SEAM.**  OpenAI's `hills_program` (TensorProgram.lean:95) word for word, with `hillsU8XMHH`
(exponent `1 - 16197706/10^10`) in the place of `hills` (exponent `1 - 2/10^11`): the kernel program of
the certificate `GXD.U8XMHH` (the one behind `WHT.wht_main_block_U8XMHHx`), for every colour,
both modes and every batch of inputs, in `2^k * hillsU8XMHH k` operations. -/
theorem hillsU8XMHH_program (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*hillsU8XMHH (k i)) :=
  SeamU8XMHH.engine_U8XMHHx cl m k v

/-- **General exponent**: the same at every exponent `z ≥ alphaU8XMHH`, envelope `⌈(k+1)^z⌉`. -/
theorem envelope_program_geU8XMHH (z : ℝ) (hz : alphaU8XMHH ≤ z) (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  (hillsU8XMHH_program cl m k v).mono fun i => Nat.mul_le_mul_left _ (envelope_mono_expU8XMHH hz (k i))

end RAM

/-! ## 5. Number facts for the copy of `Main.lean` -/

namespace SeamU8XMHH

lemma alphaU8XMHH_lt_alpha : alphaU8XMHH < alpha := by norm_num [alphaU8XMHH, alpha]

/-- for `Skies.alpha_theta` (Main.lean:39): below the paper's exponent. -/
lemma alphaU8XMHH_theta : alphaU8XMHH < theta := alphaU8XMHH_lt_alpha.trans Skies.alpha_theta

/-- for `Skies.decimal_range` (Main.lean:59). -/
lemma alphaU8XMHH_decimal : alphaU8XMHH < decimalExponent ∧ decimalExponent < 1 := by
  norm_num [decimalExponent, alphaU8XMHH]

/-- Target exponent of the final Fourier statement: strictly above `alphaU8XMHH`, so that the
`(log log n)^2` factor of the reduction is absorbed. -/
def targetExponentU8XMHH : ℝ := 1 - 16197705/(10:ℝ)^10
/-- Fallback exponent. -/
def fallbackExponentU8XMHH : ℝ := 1 - 102210/(10:ℝ)^8

lemma alphaU8XMHH_lt_target : alphaU8XMHH < targetExponentU8XMHH := by norm_num [alphaU8XMHH, targetExponentU8XMHH]
lemma alphaU8XMHH_lt_fallback : alphaU8XMHH < fallbackExponentU8XMHH := by norm_num [alphaU8XMHH, fallbackExponentU8XMHH]
lemma target_lt_fallbackU8XMHH : targetExponentU8XMHH < fallbackExponentU8XMHH := by
  norm_num [targetExponentU8XMHH, fallbackExponentU8XMHH]
lemma target_posU8XMHH : 0 < targetExponentU8XMHH := by norm_num [targetExponentU8XMHH]
lemma target_lt_oneU8XMHH : targetExponentU8XMHH < 1 := by norm_num [targetExponentU8XMHH]
lemma fallback_lt_oneU8XMHH : fallbackExponentU8XMHH < 1 := by norm_num [fallbackExponentU8XMHH]
lemma target_lt_decimalU8XMHH : targetExponentU8XMHH < decimalExponent := by
  norm_num [targetExponentU8XMHH, decimalExponent]
lemma fallback_lt_decimalU8XMHH : fallbackExponentU8XMHH < decimalExponent := by
  norm_num [fallbackExponentU8XMHH, decimalExponent]
lemma target_lt_thetaU8XMHH : targetExponentU8XMHH < theta :=
  (show targetExponentU8XMHH < alpha by norm_num [targetExponentU8XMHH, alpha]).trans Skies.alpha_theta
lemma fallback_lt_thetaU8XMHH : fallbackExponentU8XMHH < theta :=
  (show fallbackExponentU8XMHH < alpha by norm_num [fallbackExponentU8XMHH, alpha]).trans Skies.alpha_theta

end SeamU8XMHH
end
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.RAM.hillsU8XMHH_program
#print axioms OAI.PowerSaving.RAM.envelope_program_geU8XMHH
#print axioms OAI.PowerSaving.SeamU8XMHH.alphaU8XMHH_theta
#print axioms OAI.PowerSaving.hillsU8XMHH_real
