import Work.GCert.Data.B2Ge8Main

/-!
# THE SEAM for the Fourier transfer: our kernel program in OpenAI's `hills_program` form

New file of this repository, not a copy of an OpenAI file.  It contains lines 95-99 of
lean/OAI/Computability/FourierTransform/TensorProgram.lean (the statement of `hills_program`) and line 114 of
lean/OAI/Computability/FourierTransform/TensorSaving.lean (the definition of `hills`) of github.com/openai/math
(commit fd4aeeb, Apache-2.0), with `hills`, `alpha` replaced by `hillsE8`, `alphaE8`.

Twin of Work/Fourier/Seam.lean for the certificate `GXD.E8` (label width 9; Work/GCert/Data/B2Ge8.lean) in the
place of `GXD.P193`: the same text with the names and numbers of that certificate, and with every declared name
changed (`E8` in the place of `Z`, `E8` appended to the other names, namespace `SeamE8` for `Seam`), so that the
two files have no declaration in common.  This file imports the E8 instance only, not the published seam; that
is why the general lemmas of sections 1 and 2 are stated here again, under new names, with the same proofs.

OpenAI's Fourier proof uses the tensor network through ONE theorem, `RAM.hills_program`
(OAI/Computability/FourierTransform/TensorProgram.lean:95), with ONE call site (`Bench.gleamT`,
SectorAlgorithm.lean:45), plus four envelope facts (`hills_pos`, `hills_mono`, `hills_le`, `hills_real`,
TensorSaving.lean:116-133) and the number facts `alpha_pos`, `alpha_lt`, `Skies.alpha_theta`,
`Skies.decimal_range`.  This file gives all of them at OUR exponent

  `alphaE8 = 1 - 8762479/10^10`      (OpenAI: `alpha = 1 - 2/10^11`),

under OpenAI's names with `E8` after `hills` / `alpha`, in OpenAI's own namespaces.

* `GX.engine_of_ginvE8`, `GX.engine_of_halvesE8`   twins of `GX.wht_of_ginv`, `GX.wht_of_halves`
  (Work/GCert/Chain/End.lean) ending in the kernel program (`RAM.engine_program_ggroup_list`) instead of the
  Walsh-Hadamard corollary: every colour `cl`, both modes `m`, every batch of inputs.
* `SeamE8.engine_B2Ge8x_of_scal`, `_of_ids`   twins of `WHT.wht_B2Ge8x_of_scal`, `_of_ids`
  (Work/GCert/Data/B2Ge8.lean): same certificate `GXD.E8`, same unit list, same rate fact
  `BlockAccounting.fold_B2Ge8x`, fill parameter `s = 40`.
* `RAM.hillsE8_program`   THE SEAM, closed with the kernel-checked scalar half `GXD.E8.scalOK`.
* `RAM.envelope_program_geE8`   the same at every exponent `z ≥ alphaE8`.

No `sorry`, no `native_decide`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
noncomputable section

/-! ## 1. The exponent and the envelope facts -/

/-- Our exponent, in the place of OpenAI's `alpha` (TensorSaving.lean:81). -/
def alphaE8 : ℝ := 1 - 8762479/(10:ℝ)^10

theorem alphaE8_pos : 0 < alphaE8 := by norm_num [alphaE8]
theorem alphaE8_lt : alphaE8 < 1 := by norm_num [alphaE8]

lemma envelope_monoE8 (z : ℝ) (hz : 0 ≤ z) {a b : ℕ} (h : a ≤ b) : envelope z a ≤ envelope z b := by
  unfold envelope
  gcongr

-- `envelope_le` (for `z ≤ 1`) and `envelope_real` (for `0 ≤ z`) already exist: Work/Block/WHT.lean:26, :31.

lemma envelope_mono_expE8 {z z' : ℝ} (h : z ≤ z') (k : ℕ) : envelope z k ≤ envelope z' k := by
  unfold envelope
  exact Nat.ceil_le_ceil (Real.rpow_le_rpow_of_exponent_le (by simp) h)

/-- Natural envelope at our exponent, in the place of OpenAI's `hills` (TensorSaving.lean:114). -/
def hillsE8 (k : ℕ) := Nat.ceil (((k:ℝ)+1)^alphaE8)

lemma hillsE8_eq_envelope : hillsE8 = envelope alphaE8 := rfl

lemma hillsE8_pos (k : ℕ) : 1≤hillsE8 k := envelope_pos alphaE8 k

lemma hillsE8_mono {a b : ℕ} (h : a≤b) : hillsE8 a ≤ hillsE8 b := envelope_monoE8 alphaE8 alphaE8_pos.le h

lemma hillsE8_le (k : ℕ) : hillsE8 k ≤ k+1 := envelope_le alphaE8 alphaE8_lt.le k

lemma hillsE8_real (k : ℕ) : (hillsE8 k:ℝ) ≤ 2*((k:ℝ)+1)^alphaE8 := envelope_real alphaE8 alphaE8_pos.le k

/-! ## 2. Two halves + rate fact => the kernel program (twins of `Work/GCert/Chain/End.lean`) -/

namespace GX
open Binary Matrix Finset RAM RAM.Ty SS CB RF BR BG GF FoldRate
universe U

section End
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

/-- **Kernel program from a generalised invocation package and a rate fact** (twin of `wht_of_ginv`). -/
theorem engine_of_ginvE8 (I : GInv H T Sl C) (hH : 1 ≤ Fintype.card H)
    (u W : ℕ) (hu : 5 * Fintype.card H = u)
    (hW : 4 * Fintype.card T + Fintype.card Sl = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ,
      5 * I.cost φ + (Fintype.card T : ℝ) * bankPrice (Fintype.card H) φ = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ))
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  engine_program_ggroup_list u (card_L5.trans hu) (by omega)
    (Sum.inl : BLive (Orth (L5 H)) T Sl → BRole (Orth (L5 H)) T Sl C) Sum.inl_injective
    (Fintype.card (Orth (L5 H))) W Orth.card_pos hW1 (by rw [card_BLive, hW])
    _ (bridge2_gcert I hH) U (fun φ => by rw [hU φ]) s z hz hz1 hnum cl m k v

/-- **Scalar half + label half + rate fact => the kernel program** (twin of `wht_of_halves`). -/
theorem engine_of_halvesE8 (S : Scal T Sl C) (L : Lab H S) (hH : 1 ≤ Fintype.card H)
    (u W : ℕ) (hu : 5 * Fintype.card H = u)
    (hW : 4 * Fintype.card T + Fintype.card Sl = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ,
      5 * L.cost φ + (Fintype.card T : ℝ) * bankPrice (Fintype.card H) φ = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ))
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  engine_of_ginvE8 L.inv hH u W hu hW hW1 U hU s z hz hz1 hnum cl m k v

end End
end GX

/-! ## 3. The certificate `GXD.E8` (twins of `Work/GCert/Data/B2Ge8.lean`) -/

namespace SeamE8
open Binary Matrix Finset RAM RAM.Ty
universe U

variable {C : Type} [Fintype C] [DecidableEq C] (S : GX.Scal (Fin GXD.E8.raw.v) (Fin GXD.E8.raw.R) C)
  {emb : CR.L3 (Fin GXD.E8.raw.v) (Fin GXD.E8.raw.R) → Nat} {un : Nat → CR.L3 (Fin GXD.E8.raw.v) (Fin GXD.E8.raw.R)}
  (cf : GXD.Co → SSC.Coef)

/-- **Kernel program at `1 - 8762479/10^10`** over any scalar half of `GXD.E8.raw`
(twin of `WHT.wht_B2Ge8x_of_scal`: same arguments, `GX.engine_of_halvesE8` for `GX.wht_of_halves`). -/
theorem engine_B2Ge8x_of_scal (Ro : GLab.Roles GXD.E8.raw.v GXD.E8.raw.R emb un)
    (hMA : S.MA = SSC.matP un (GLab.mic cf (GXD.E8.raw.gatesA.flatMap GXD.Gate.adds)))
    (hMAi : S.MAi = SSC.matPinv un (GLab.mic cf (GXD.E8.raw.gatesA.flatMap GXD.Gate.adds)))
    (hMB : S.MB = SSC.matP un (GLab.mic cf (GXD.E8.raw.gatesB.flatMap GXD.Gate.adds)))
    (hMBi : S.MBi = SSC.matPinv un (GLab.mic cf (GXD.E8.raw.gatesB.flatMap GXD.Gate.adds)))
    (reg : C → Nat) (hCc : ∀ k (q : Fin GXD.E8.raw.R), S.Cc k q ≠ 0 → reg k = 2 * GXD.E8.raw.v + q.val)
    (hreg : ∀ k, reg k ∈ GXD.E8.raw.ret.map Prod.fst) (hC : Fintype.card C = 8)
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 8762479/(10:ℝ)^10) (k i)) :=
  GX.engine_of_halvesE8 S (GXD.E8.valid.lab S Ro cf hMA hMAi hMB hMBi reg hCc hreg)
    (by rw [Fintype.card_fin]; decide) 45 1263 (by rw [Fintype.card_fin]; rfl)
    (by rw [Fintype.card_fin, Fintype.card_fin]; rfl) (by norm_num) BlockAccounting.blocksB2Ge8x
    (fun φ => by
      rw [GXD.E8.lab_price S Ro cf hMA hMAi hMB hMBi reg hCc hreg hC φ, Fintype.card_fin, Fintype.card_fin]
      exact BridgeRate.B2Ge8.unit_cost_x φ)
    40 (1 - 8762479/(10:ℝ)^10) (by norm_num) (by norm_num) BlockAccounting.fold_B2Ge8x cl m k v

/-- the same from the FIVE scalar identities of the matrices `GS.MA GS.MAi GS.MB GS.MBi GS.Jrm GS.Ccm` of
`GXD.E8.raw` (twin of `WHT.wht_B2Ge8x_of_ids`). -/
theorem engine_B2Ge8x_of_ids (hv : 0 < GXD.E8.raw.v)
    (hAi : GS.MAi GXD.E8.raw hv * GS.MA GXD.E8.raw hv = 1) (hBi : GS.MBi GXD.E8.raw hv * GS.MB GXD.E8.raw hv = 1)
    (hx : ∀ (t : Fin GXD.E8.raw.v) (j : GS.Rl GXD.E8.raw), (GS.MB GXD.E8.raw hv * CR.scatM (GS.Jrm GXD.E8.raw * GS.Ccm GXD.E8.raw)
      * GS.MA GXD.E8.raw hv) (CR.lX t) j = if CR.lX t = j then 1 else 0)
    (hy : ∀ (i : GS.Rl GXD.E8.raw) (t : Fin GXD.E8.raw.v), (GS.MB GXD.E8.raw hv * CR.scatM (GS.Jrm GXD.E8.raw * GS.Ccm GXD.E8.raw)
      * GS.MA GXD.E8.raw hv) i (CR.lY t) = if i = CR.lY t then 1 else 0)
    (hid : ∀ S t : Fin GXD.E8.raw.v, (GS.MB GXD.E8.raw hv * CR.scatM (GS.Jrm GXD.E8.raw * GS.Ccm GXD.E8.raw)
      * GS.MA GXD.E8.raw hv) (CR.lY S) (CR.lX t) = if S = t then 1 else 0)
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope (1 - 8762479/(10:ℝ)^10) (k i)) :=
  engine_B2Ge8x_of_scal (GXD.scalOf GXD.E8.raw hv hAi hBi hx hy hid) GS.toCoef (GXD.roles GXD.E8.raw hv)
    (GXD.scalOf_MA _ hv hAi hBi hx hy hid) (GXD.scalOf_MAi _ hv hAi hBi hx hy hid)
    (GXD.scalOf_MB _ hv hAi hBi hx hy hid) (GXD.scalOf_MBi _ hv hAi hBi hx hy hid)
    (GXD.reg GXD.E8.raw) (GXD.reg_slot GXD.E8.raw hv hAi hBi hx hy hid) (GXD.reg_mem GXD.E8.raw)
    (by rw [Fintype.card_fin]; decide) cl m k v

end SeamE8

/-! ## 4. THE SEAM -/

namespace RAM
open Binary Matrix Ty Finset
universe U

/-- **THE SEAM.**  OpenAI's `hills_program` (TensorProgram.lean:95) word for word, with `hillsE8`
(exponent `1 - 8762479/10^10`) in the place of `hills` (exponent `1 - 2/10^11`): the kernel program of
the certificate `GXD.E8` (the one behind `WHT.wht_main_block_B2Ge8x`), for every colour,
both modes and every batch of inputs, in `2^k * hillsE8 k` operations. -/
theorem hillsE8_program (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*hillsE8 (k i)) :=
  SeamE8.engine_B2Ge8x_of_ids GXD.E8.scalOK.hv GXD.E8.scalOK.hAi GXD.E8.scalOK.hBi
    GXD.E8.scalOK.hx GXD.E8.scalOK.hy GXD.E8.scalOK.hid cl m k v

/-- **General exponent**: the same at every exponent `z ≥ alphaE8`, envelope `⌈(k+1)^z⌉`. -/
theorem envelope_program_geE8 (z : ℝ) (hz : alphaE8 ≤ z) (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ)
    (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  (hillsE8_program cl m k v).mono fun i => Nat.mul_le_mul_left _ (envelope_mono_expE8 hz (k i))

end RAM

/-! ## 5. Number facts for the copy of `Main.lean` -/

namespace SeamE8

lemma alphaE8_lt_alpha : alphaE8 < alpha := by norm_num [alphaE8, alpha]

/-- for `Skies.alpha_theta` (Main.lean:39): below the paper's exponent. -/
lemma alphaE8_theta : alphaE8 < theta := alphaE8_lt_alpha.trans Skies.alpha_theta

/-- for `Skies.decimal_range` (Main.lean:59). -/
lemma alphaE8_decimal : alphaE8 < decimalExponent ∧ decimalExponent < 1 := by
  norm_num [decimalExponent, alphaE8]

/-- Target exponent of the final Fourier statement: strictly above `alphaE8`, so that the
`(log log n)^2` factor of the reduction is absorbed. -/
def targetExponentE8 : ℝ := 1 - 8762478/(10:ℝ)^10
/-- Fallback exponent. -/
def fallbackExponentE8 : ℝ := 1 - 87624/(10:ℝ)^8

lemma alphaE8_lt_target : alphaE8 < targetExponentE8 := by norm_num [alphaE8, targetExponentE8]
lemma alphaE8_lt_fallback : alphaE8 < fallbackExponentE8 := by norm_num [alphaE8, fallbackExponentE8]
lemma target_lt_fallbackE8 : targetExponentE8 < fallbackExponentE8 := by
  norm_num [targetExponentE8, fallbackExponentE8]
lemma target_posE8 : 0 < targetExponentE8 := by norm_num [targetExponentE8]
lemma target_lt_oneE8 : targetExponentE8 < 1 := by norm_num [targetExponentE8]
lemma fallback_lt_oneE8 : fallbackExponentE8 < 1 := by norm_num [fallbackExponentE8]
lemma target_lt_decimalE8 : targetExponentE8 < decimalExponent := by
  norm_num [targetExponentE8, decimalExponent]
lemma fallback_lt_decimalE8 : fallbackExponentE8 < decimalExponent := by
  norm_num [fallbackExponentE8, decimalExponent]
lemma target_lt_thetaE8 : targetExponentE8 < theta :=
  (show targetExponentE8 < alpha by norm_num [targetExponentE8, alpha]).trans Skies.alpha_theta
lemma fallback_lt_thetaE8 : fallbackExponentE8 < theta :=
  (show fallbackExponentE8 < alpha by norm_num [fallbackExponentE8, alpha]).trans Skies.alpha_theta

end SeamE8
end
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.RAM.hillsE8_program
#print axioms OAI.PowerSaving.RAM.envelope_program_geE8
#print axioms OAI.PowerSaving.SeamE8.alphaE8_theta
#print axioms OAI.PowerSaving.hillsE8_real
