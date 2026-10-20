import Work.XHostT.TDispW
import Work.XHostT.Data.D3945.Disp

/-!
# (key: bridge-host-cross) #3945's two-provider network: 1,672 calls of `extra`, 1,568 of `gauge44`

The network theorem of `TDisp` / `TDispW` specialised to #3945's real dispatch table (`D3945.disp`, 360 replicas,
generated from the hub's `DISPATCH.json` / `ORDER.json` by p-disp/tools/mk_disp3945.py): for ANY cross-replica
network `N` (geometry, global order, host table) whose call `(ρ, l)` runs provider `D3945.disp (ρ, l)`
(`0` = `extra`, `1` = `gauge44`), the certificate price is `|Γ| * (1672 * cost_extra + 1568 * cost_gauge44 + restU)`;
with private scratch per call (`e : C ↪ C'`) the same.  The providers are generic here: i-final instantiates
`I0 := extra` (p-wgen's `GXD.S5E`), `I1 := gauge44` (`GXD.S5B`).
Written by p-disp (wht45/lever, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHostT
namespace D3945
open Binary Matrix Finset RAM SS CB RF BR GF GX XHost XHostD
noncomputable section

section Net
variable {H T Sl C C' α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype C'] [DecidableEq C'] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {nb : ℕ}

/-- **the price of #3945's 3,240 calls**: 1,672 of provider 0 (`extra`), 1,568 of provider 1 (`gauge44`). -/
theorem calls (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (φ : ℕ → ℝ) :
    ∑ ρ : Fin 360, ∑ l : Fin 9, (dispIc (Prov2.mk' I0 I1 h) disp ρ l).cost φ = 1672 * I0.cost φ + 1568 * I1.cost φ := by
  rw [sum_calls_disp2, nExtra, nGauge]
  norm_num

/-- the unit price of a network dispatched by #3945's table. -/
theorem unit (N : XNet H T Sl C α Γ nb 360) (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv)
    (hN : ∀ ρ l, N.Ic ρ l = dispIc (Prov2.mk' I0 I1 h) disp ρ l) (φ : ℕ → ℝ) :
    N.unit φ = 1672 * I0.cost φ + 1568 * I1.cost φ + N.restU φ := by
  rw [N.unit_disp2 I0 I1 h disp hN φ, nExtra, nGauge]
  norm_num

/-- **#3945's cross-replica network with its two provider words** (any geometry, order and host table). -/
theorem certificate (N : XNet H T Sl C α Γ nb 360) (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv)
    (hN : ∀ ρ l, N.Ic ρ l = dispIc (Prov2.mk' I0 I1 h) disp ρ l) :
    GCert α (Sum.inl : BLiveT Γ T nb 360 → BRoleT Γ T nb C 360)
      (fun φ => (Fintype.card Γ : ℝ) * (1672 * I0.cost φ + 1568 * I1.cost φ + N.restU φ)) :=
  (xnet_certificate N).cast (fun φ => by rw [unit N I0 I1 h hN φ])

/-- **the same with private scratch in every call** (copies re-indexed by `e`; e.g. the 91 zero buffers). -/
theorem certificate_scratch (N : XNet H T Sl C α Γ nb 360) (e : C ↪ C') (I0 I1 : GInvE H T Sl C)
    (h : I1.tv = I0.tv) (hN : ∀ ρ l, N.Ic ρ l = dispIc (Prov2.mk' I0 I1 h) disp ρ l) :
    GCert α (Sum.inl : BLiveT Γ T nb 360 → BRoleT Γ T nb C' 360)
      (fun φ => (Fintype.card Γ : ℝ) * (1672 * I0.cost φ + 1568 * I1.cost φ + N.restU φ)) :=
  (xnet_reC_certificate N e).cast (fun φ => by rw [unit N I0 I1 h hN φ])

end Net
end
end D3945
end XHostT
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHostT.D3945.calls
#print axioms OAI.PowerSaving.XHostT.D3945.unit
#print axioms OAI.PowerSaving.XHostT.D3945.certificate
#print axioms OAI.PowerSaving.XHostT.D3945.certificate_scratch
