import Work.Bridge401.GGeom
import Work.PWGen.B401W.GGeomW

import Work.Gauge.Iface
import Work.GCert.Chain.End
/-! (key: p-wgen) WINDOWED copy, made by p-wgen's tools/mk_b401w.py from Work/Bridge401 (round wht45/lever, 19 Oct
2026): the package is p-warpchk's `KLW.GInvEW H T Sl C Ok` (routes only at windows with `Ok`), every declared name
carries a W, and the window facts `Ok` are threaded to the two places that call the package's routes.  The original
text follows. -/

/-!
# (key: bridge401) THE END THEOREMS for the entrance-gauged nine-stage layout: gauged halves + bank schedule + rate fact
=> WHT / kernel program

Credits: the entrance gauges and the residual banks are Dugongue's (CrocSwap/integer-mult-bounds PR #401); the bank
packing bound is Lydakis's (PR #275); the compressed nine-stage layout, its banks and padding are jacklightChen's (PR #377);
the bridged word B_K is DaysSky's (PR #363). The engine (`BlockWHT.wht_main_of_ggroup_list`,
`RAM.engine_program_ggroup_list`), the invocation package and the stage geometry are Jacob Sussman's
(wht-power-saving-lean, Apache-2.0); this file follows Work/Bridge377/End.lean. The gauged halves (`GScal`, `GLab`,
`GLab.inv`, a definition with `GLab.inv_cost` by `rfl`) are DFT Bridge3's (L1), the network theorem `bridge401W_gcert` is L3's (DFT New Sums). Written by DFT Lean
(#401 L4).

One unit: label space `BC.LC (Fin 8) H l` of dimension `u = 73` (h = 9, K = 4), `W = 64 * card T + S.nb` live roles
(`S.nb` banks over all 9 stages), unit price `unitPrice401W I S.pad φ = 72 * I.cost φ + card T * bankPrice377 φ
+ ∑ b, bcost φ (S.pad b)` (one final pad block per bank). The pivot `l` is any element of `H`. No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace GX
open Binary Matrix Finset RAM RAM.Ty SS CB RF BR BG GF FoldRate
noncomputable section

universe U

section EndW
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
variable {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}

/-- **Walsh-Hadamard transform from a generalised invocation package and a rate fact (PR #401 layout).** -/
theorem wht_of_ginv_401W (I : KLW.GInvEW H T Sl C Ok) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : S.fits I.ent)
    (hOk : ∀ (l : H) (hH2 : 2 ≤ Fintype.card H) (d : Orth (BC.LC (Fin 8) H l)) (q : Finset (Fin 8 × BC.Hp H l)),
      Ok (gwinC hB401 hH2 d q))
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
    ∃ solve W', WHT.WHTProgram solve W' ∧ BlockWHT.WHTTimeBoundsAt z W' := by
  obtain ⟨l⟩ : Nonempty H := Fintype.card_pos_iff.mp (by omega)
  exact BlockWHT.wht_main_of_ggroup_list u (by rw [BC.card_LC, hH, Fintype.card_fin]; omega) (by omega)
    (Sum.inl : BLive401 (Orth (BC.LC (Fin 8) H l)) T S.nb → BRole401 (Orth (BC.LC (Fin 8) H l)) T S.nb C)
    Sum.inl_injective (Fintype.card (Orth (BC.LC (Fin 8) H l))) W Orth.card_pos hW1 (by rw [card_BLive401, hW])
    _ (bridge401W_gcert I hH l S hres (hOk l (by omega))) U (fun φ => by rw [unitPrice401W, hU φ]) s z hz hz1 hnum

/-- **Kernel program from a generalised invocation package and a rate fact (PR #401 layout)**: the form that
OpenAI's Fourier proof uses (`hills_program`); for the seam of the Fourier chain. -/
theorem engine_of_ginv_401W (I : KLW.GInvEW H T Sl C Ok) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : S.fits I.ent)
    (hOk : ∀ (l : H) (hH2 : 2 ≤ Fintype.card H) (d : Orth (BC.LC (Fin 8) H l)) (q : Finset (Fin 8 × BC.Hp H l)),
      Ok (gwinC hB401 hH2 d q))
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
      (fun i => 2^k i*envelope z (k i)) := by
  obtain ⟨l⟩ : Nonempty H := Fintype.card_pos_iff.mp (by omega)
  exact engine_program_ggroup_list u (by rw [BC.card_LC, hH, Fintype.card_fin]; omega) (by omega)
    (Sum.inl : BLive401 (Orth (BC.LC (Fin 8) H l)) T S.nb → BRole401 (Orth (BC.LC (Fin 8) H l)) T S.nb C)
    Sum.inl_injective (Fintype.card (Orth (BC.LC (Fin 8) H l))) W Orth.card_pos hW1 (by rw [card_BLive401, hW])
    _ (bridge401W_gcert I hH l S hres (hOk l (by omega))) U (fun φ => by rw [unitPrice401W, hU φ]) s z hz hz1 hnum cl m k v

end EndW
end
end GX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.GX.wht_of_ginv_401W
#print axioms OAI.PowerSaving.GX.engine_of_ginv_401W
