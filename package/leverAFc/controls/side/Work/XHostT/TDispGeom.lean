import Work.XHostT.TGeom
import Work.XHostT.D3945Net

/-!
# (key: bridge-host-cross) The cross-replica geometry instance with TWO PROVIDERS and PRIVATE SCRATCH

p-a3's `bridgeHostT_gcert` (TGeom: the network `xnetGeom` on Shea's label space `LC (Fin 8) H l` from a cross-replica
host schedule `HostSchedT` with every slot hosted) for the per-call invocations `dispIc (Prov2.mk' I0 I1 h) disp`,
every call with its copies re-indexed by `e : C ↪ C'` (private scratch, e.g. the 91 zero buffers):

* `resIc_disp`: the residuals the host check reads are the providers' (`D.res (disp (ρ, l)) q`, `rfl`);
  `resIc_reC`: the private scratch does not change them (`rfl`);
* `unitPriceT_disp2`, `unitPriceT_reC`: the price `n0 * cost0 + n1 * cost1 + Σ hcostT`, unchanged by the scratch;
* **`bridgeHostT_disp_scratch_gcert`**, `bridgeHostT_disp_gcert` (no scratch) and `D3945.bridge_gcert` (#3945's
  real table: `1672 * cost_extra + 1568 * cost_gauge44 + Σ hcostT`).
Written by p-disp (wht45/lever, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHostT
open Binary Matrix Finset RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD
noncomputable section

section Geo
variable {H T Sl C C' : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype C'] [DecidableEq C'] {nR : ℕ}

/-- the residuals a host check reads: those of the call's provider. -/
lemma resIc_disp {ι : Type} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) (ρ : Fin nR) (l : Fin 9) (q : Sl) :
    resIc (dispIc D disp) ρ l q = D.res (disp (ρ, l)) q := rfl

/-- the private scratch does not change the residuals (the host check is the same). -/
lemma resIc_reC (Ic : Fin nR → Fin 9 → GInvE H T Sl C) (e : C ↪ C') :
    resIc (fun ρ l => (Ic ρ l).reC e) = resIc Ic := rfl

/-- the private scratch is not paid. -/
theorem unitPriceT_reC (Ic : Fin nR → Fin 9 → GInvE H T Sl C) (e : C ↪ C') (S : HostSchedT nR Sl T) (φ : ℕ → ℝ) :
    unitPriceT (fun ρ l => (Ic ρ l).reC e) S φ = unitPriceT Ic S φ := rfl

/-- the price with two providers. -/
theorem unitPriceT_disp2 (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (disp : Fin nR × Fin 9 → Fin 2)
    (S : HostSchedT nR Sl T) (φ : ℕ → ℝ) :
    unitPriceT (dispIc (Prov2.mk' I0 I1 h) disp) S φ
      = (nCalls disp 0 : ℝ) * I0.cost φ + (nCalls disp 1 : ℝ) * I1.cost φ + ∑ ρ : Fin nR, ∑ c : Fin 9, hcostT S ρ c φ := by
  unfold unitPriceT
  rw [sum_calls_disp2]

/-- **THE GEOMETRY INSTANCE WITH TWO PROVIDERS AND PRIVATE SCRATCH**: calls dispatched by `disp` to `I0`/`I1`, every
call with its copies re-indexed by `e` (private scratch), every slot hosted by the checked cross-replica schedule `S`
(residuals: the providers'). -/
theorem bridgeHostT_disp_scratch_gcert (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (disp : Fin nR × Fin 9 → Fin 2)
    (e : C ↪ C') (hH : Fintype.card H = 9) (l : H) (S : HostSchedT nR Sl T)
    (hfit : S.fits (fun ρ l q => (Prov2.mk' I0 I1 h).res (disp (ρ, l)) q)) :
    GCert (BC.LC (Fin 8) H l)
      (Sum.inl : BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 nR → BRoleT (Orth (BC.LC (Fin 8) H l)) T 0 C' nR)
      (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * ((nCalls disp 0 : ℝ) * I0.cost φ
        + (nCalls disp 1 : ℝ) * I1.cost φ + ∑ ρ : Fin nR, ∑ c : Fin 9, hcostT S ρ c φ)) :=
  (bridgeHostT_gcert (I0.reC e) (fun ρ l => (dispIc (Prov2.mk' I0 I1 h) disp ρ l).reC e) (fun _ _ => rfl) hH l S
    hfit).cast (fun φ => by
      rw [unitPriceT_reC (dispIc (Prov2.mk' I0 I1 h) disp) e S φ, unitPriceT_disp2])

/-- the same without private scratch. -/
theorem bridgeHostT_disp_gcert (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (disp : Fin nR × Fin 9 → Fin 2)
    (hH : Fintype.card H = 9) (l : H) (S : HostSchedT nR Sl T)
    (hfit : S.fits (fun ρ l q => (Prov2.mk' I0 I1 h).res (disp (ρ, l)) q)) :
    GCert (BC.LC (Fin 8) H l)
      (Sum.inl : BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 nR → BRoleT (Orth (BC.LC (Fin 8) H l)) T 0 C nR)
      (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * ((nCalls disp 0 : ℝ) * I0.cost φ
        + (nCalls disp 1 : ℝ) * I1.cost φ + ∑ ρ : Fin nR, ∑ c : Fin 9, hcostT S ρ c φ)) :=
  (bridgeHostT_gcert I0 (dispIc (Prov2.mk' I0 I1 h) disp) (fun _ _ => rfl) hH l S hfit).cast
    (fun φ => by rw [unitPriceT_disp2])

/-- **#3945's real table on the geometry instance**: 1,672 calls of `I0` (extra), 1,568 of `I1` (gauge44), private
scratch per call. -/
theorem D3945.bridge_gcert (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (e : C ↪ C') (hH : Fintype.card H = 9) (l : H)
    (S : HostSchedT 360 Sl T) (hfit : S.fits (fun ρ l q => (Prov2.mk' I0 I1 h).res (D3945.disp (ρ, l)) q)) :
    GCert (BC.LC (Fin 8) H l)
      (Sum.inl : BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 360 → BRoleT (Orth (BC.LC (Fin 8) H l)) T 0 C' 360)
      (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * (1672 * I0.cost φ + 1568 * I1.cost φ
        + ∑ ρ : Fin 360, ∑ c : Fin 9, hcostT S ρ c φ)) :=
  (bridgeHostT_disp_scratch_gcert I0 I1 h D3945.disp e hH l S hfit).cast
    (fun φ => by rw [D3945.nExtra, D3945.nGauge]; norm_num)

end Geo
end
end XHostT
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHostT.unitPriceT_disp2
#print axioms OAI.PowerSaving.XHostT.bridgeHostT_disp_scratch_gcert
#print axioms OAI.PowerSaving.XHostT.bridgeHostT_disp_gcert
#print axioms OAI.PowerSaving.XHostT.D3945.bridge_gcert
