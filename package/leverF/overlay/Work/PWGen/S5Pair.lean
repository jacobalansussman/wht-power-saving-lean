import Work.GCert.Data.S5BGlue401
import Work.GCert.Data.S5EGlue401

/-!
# (key: p-wgen) the two #3945 provider packages side by side, for the two-provider dispatch (p-disp)

p-wgen, round wht45/lever, 19 Oct 2026.  `GXD.S5B.ginvW` (gauge44) and `GXD.S5E.ginvW` (extra) are windowed packages
`KLW.GInvEW (Fin 9) (Fin 112) (Fin 408) (Fin 99) okW` with the same window condition (each word's `okW` is the same
lambda) and the same target vectors (the two words have the same 112 ports, checked by the kernel).  Here both are
stated over `GXD.S5B.okW`, with `tv_eq` in the shape of `XHostD.Prov2.mk'`'s hypothesis (`I1.tv = I0.tv` with
`I0 := extra`, `I1 := gauge44`), and their prices.  No `sorry`.
-/

namespace PWG
open OAI.PowerSaving OAI.PowerSaving.CB

/-- provider 0 of p-disp: the extra word -/
noncomputable def ginvE : KLW.GInvEW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S5B.okW := GXD.S5E.ginvW
/-- provider 1 of p-disp: the gauge44 word -/
noncomputable def ginvB : KLW.GInvEW (Fin 9) (Fin 112) (Fin 408) (Fin 99) GXD.S5B.okW := GXD.S5B.ginvW

theorem okW_eq : @GXD.S5E.okW = @GXD.S5B.okW := rfl

theorem ports_eq : GXD.S5B.wc.ports = GXD.S5E.wc.ports := by decide +kernel

/-- the same target vectors (`I1.tv = I0.tv` with `I0 := ginvE`, `I1 := ginvB`) -/
theorem tv_eq : ginvB.tv = ginvE.tv := by
  funext t
  show GXD.S5B.wc.tv t = GXD.S5E.wc.tv t
  unfold KLW.WCert.tv
  rw [ports_eq]
  try rfl

theorem ginvE_cost (φ : ℕ → ℝ) : ginvE.cost φ = GXD.S5E.invPrice φ := GXD.S5E.ginv_cost φ
theorem ginvB_cost (φ : ℕ → ℝ) : ginvB.cost φ = GXD.S5B.invPrice φ := GXD.S5B.ginv_cost φ

/-- the two words cost the same at every price list (same local profile) -/
theorem same_price (φ : ℕ → ℝ) : ginvE.cost φ = ginvB.cost φ := by
  rw [ginvE_cost, ginvB_cost]
  unfold GXD.S5E.invPrice GXD.S5B.invPrice
  try rfl

end PWG

#print axioms PWG.tv_eq
#print axioms PWG.same_price
