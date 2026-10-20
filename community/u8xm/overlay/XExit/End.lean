import Work.Bridge401.End
import XExit.Halves

/-!
# (key: x-exit) THE END THEOREMS for exit-frame helpers (PR #401 nine-stage layout, Shea's bridge unchanged)

`engine_of_ginvx_401` / `wht_of_ginvx_401`: an exit-frame package `I : GInvX` (slot `q` from `ent q` to `ex q`, same
base) + a bank schedule whose loads are the EXIT residuals `dim ex q - dim ent q` (`fitsX`) + the rate fact
=> the kernel program / Walsh-Hadamard program, with unit price `72 * I.cost + card T * bankPrice377 + pads`.
Proof: XExit/Pack's reduction `GInvX.toE'` (same price, stand-in entrance labels of dimension
`9 - (dim ex - dim ent)`, free adapters) and Shea's `GX.engine_of_ginv_401` / `GX.wht_of_ginv_401`.
`engine_of_halvesX_401`: the same from a scalar half `GScal` and an exit-labelled label half `GLabX` (XExit/Halves).

Credits: Shea's / Mamano's `engine_of_ginv_401` (Work/Bridge401/End, shea256/fourier-transform-below-nlogn
lean/complex588) on Dugongue's #401, Lydakis's #275, jacklightChen's #377, DaysSky's #363; the engine is Jacob
Sussman's (wht-power-saving-lean, Apache-2.0); the exit frames are the hub's (#593/#599/#604, Kothe's #799 lineage).
No `sorry`.  NOT compiled on the Mac (needs Shea's 20 bridge modules; see x-exit-lean/COLAB.sh).
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace XEX
open Binary Matrix Finset RAM RAM.Ty SS CB RF BR BG GF FoldRate GX
noncomputable section

universe U

section End
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

/-- **Walsh-Hadamard transform from an exit-frame package and a rate fact (PR #401 layout).** -/
theorem wht_of_ginvx_401 (I : GInvX H T Sl C) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : fitsX S I)
    (u W : ℕ) (hu : 73 = u)
    (hW : 64 * Fintype.card T + S.nb = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ,
      72 * I.cost φ + (Fintype.card T : ℝ) * bankPrice377 φ + ∑ b, bcost φ (S.pad b)
        = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z < 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ)) :
    ∃ solve W', WHT.WHTProgram solve W' ∧ BlockWHT.WHTTimeBoundsAt z W' :=
  GX.wht_of_ginv_401 I.toE' hH S (fits_toE' S I hres) u W hu hW hW1 U hU s z hz hz1 hnum

/-- **Kernel program from an exit-frame package and a rate fact (PR #401 layout)**: the form the Fourier seam uses. -/
theorem engine_of_ginvx_401 (I : GInvX H T Sl C) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : fitsX S I)
    (u W : ℕ) (hu : 73 = u)
    (hW : 64 * Fintype.card T + S.nb = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ,
      72 * I.cost φ + (Fintype.card T : ℝ) * bankPrice377 φ + ∑ b, bcost φ (S.pad b)
        = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ))
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  GX.engine_of_ginv_401 I.toE' hH S (fits_toE' S I hres) u W hu hW hW1 U hU s z hz hz1 hnum cl m k v

/-- the same from the two halves, exit labels (twin of Shea's `engine_of_halves_401`). -/
theorem engine_of_halvesX_401 (G : GScal T Sl C) (L : GLabX H G) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : fitsX S L.inv)
    (u W : ℕ) (hu : 73 = u)
    (hW : 64 * Fintype.card T + S.nb = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ,
      72 * L.cost φ + (Fintype.card T : ℝ) * bankPrice377 φ + ∑ b, bcost φ (S.pad b)
        = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ))
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) :=
  engine_of_ginvx_401 L.inv hH S hres u W hu hW hW1 U hU s z hz hz1 hnum cl m k v

/-- the Walsh-Hadamard twin from the two halves, exit labels. -/
theorem wht_of_halvesX_401 (G : GScal T Sl C) (L : GLabX H G) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : fitsX S L.inv)
    (u W : ℕ) (hu : 73 = u)
    (hW : 64 * Fintype.card T + S.nb = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ,
      72 * L.cost φ + (Fintype.card T : ℝ) * bankPrice377 φ + ∑ b, bcost φ (S.pad b)
        = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z < 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ)) :
    ∃ solve W', WHT.WHTProgram solve W' ∧ BlockWHT.WHTTimeBoundsAt z W' :=
  wht_of_ginvx_401 L.inv hH S hres u W hu hW hW1 U hU s z hz hz1 hnum

end End
end
end XEX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XEX.wht_of_ginvx_401
#print axioms OAI.PowerSaving.XEX.engine_of_ginvx_401
#print axioms OAI.PowerSaving.XEX.engine_of_halvesX_401
#print axioms OAI.PowerSaving.XEX.wht_of_halvesX_401
