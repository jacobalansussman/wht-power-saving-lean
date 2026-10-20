import Work.XHostT.TDispGeom
import Work.XHostT.D2979Net

/-!
# (key: bridge-host-cross) #2979's geometry instance with its two providers and private scratch

p-disp's `D3945.bridge_gcert` (TDispGeom.lean) with the instance changed to #2979's real dispatch table
(`D2979.disp`, 74 replicas, 416 calls of `extra`, 250 of `gauge44`): price
`|Orth| * (416 * cost_extra + 250 * cost_gauge44 + Σ hcostT)` for ANY cross-replica host schedule `S` that fits the
dispatched providers' residuals (e.g. `Inst.sched2979`, with `hfit` from `sched2979_fits` and the residual link).
x-tables (wht45/lever2979, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
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

/-- **#2979's geometry instance**: two providers by the real table, private scratch per call. -/
theorem D2979.bridge_gcert (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (e : C ↪ C') (hH : Fintype.card H = 9) (l : H)
    (S : HostSchedT 74 Sl T) (hfit : S.fits (fun ρ l q => (Prov2.mk' I0 I1 h).res (D2979.disp (ρ, l)) q)) :
    GCert (BC.LC (Fin 8) H l)
      (Sum.inl : BLiveT (Orth (BC.LC (Fin 8) H l)) T 0 74 → BRoleT (Orth (BC.LC (Fin 8) H l)) T 0 C' 74)
      (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * (416 * I0.cost φ + 250 * I1.cost φ
        + ∑ ρ : Fin 74, ∑ c : Fin 9, hcostT S ρ c φ)) :=
  (bridgeHostT_disp_scratch_gcert I0 I1 h D2979.disp e hH l S hfit).cast
    (fun φ => by rw [D2979.nExtra, D2979.nGauge]; norm_num)

end Geo
end
end XHostT
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHostT.D2979.bridge_gcert
