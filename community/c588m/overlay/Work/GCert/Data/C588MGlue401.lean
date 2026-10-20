import Work.Gauge.Glue
import Work.XMut.Glue
import Work.Bridge401.GGeom
import Work.GCert.Data.C588MCert
import Work.GCert.Data.C588MScalOK
import Work.GCert.Data.Gen.C588MEChk
import Work.GCert.Data.InstC588M
import Work.GCert.Data.Gen.C588MSched

/-!
# (key: gx-data) the gauged invocation package and the bank schedule of `GXD.C588M` (PR #401 layout)
Adapted from the G5 glue by DFT Lean in Nil Mamano's record submission as in Ryan Shea's generate_global_lean.py; mutable sources: m-gen's generate_global_lean_z.py.  `ginv` is x-mutable's `GXD.ginvOfZ` (Work/XMut/Glue.lean, on DFT Bridge3's `GXD.ginvOf`) on the
kernel checks of this certificate: scalar half `scalOKZ` (sources mutable), E7 `echk<j>_ok` (three batches, SPar ⟨10, 8, 4, 8, 8⟩, sources mutable), label half `validE`.
-/

set_option maxRecDepth 100000
namespace GXD.C588M
open GXD OAI.PowerSaving OAI.PowerSaving.CB
noncomputable section

/-- the E7 batches: parameters and gauge registers -/
def es : List (GS.SPar × List Nat) := [(⟨10, 8, 4, 8, 8⟩, glb0), (⟨10, 8, 4, 8, 8⟩, glb1), (⟨10, 8, 4, 8, 8⟩, glb2)]

theorem es_gl : es.flatMap Prod.snd = gl := by first | rfl | decide

theorem es_ok : ∀ e ∈ es, GS.eCheckZ raw e.1 e.2 = true := by
  intro e he
  unfold es at he
  rcases List.mem_cons.mp he with h | he
  · subst h; exact echk0_ok
  rcases List.mem_cons.mp he with h | he
  · subst h; exact echk1_ok
  rcases List.mem_cons.mp he with h | he
  · subst h; exact echk2_ok
  exact absurd he List.not_mem_nil

/-- **the gauged invocation package** of the certificate -/
def ginv : GX.GInvE (Fin aux.p.h) (Fin raw.v) (Fin raw.R) (Fin raw.ret.length) :=
  GXD.ginvOfZ raw scalOKZ es es_ok aux (es_gl.symm ▸ validE)

theorem h9 : Fintype.card (Fin aux.p.h) = 9 := by rw [Fintype.card_fin]; rfl

/-- its price is the invocation price `invPrice` (kernel counts) -/
theorem ginv_cost (φ : ℕ → ℝ) : ginv.cost φ = invPrice φ := by
  have h5 : ∀ k : Fin raw.ret.length,
      bcost φ (GLab.cdim (aux.p.tab.get (aux.SA.get (GXD.reg raw k)))) = bcost φ 8 :=
    fun k => congrArg (bcost φ) (cen_dim _ (GXD.reg_mem raw k))
  rw [ginv, GXD.ginvOfZ_cost, Finset.sum_congr rfl (fun k _ => h5 k), Finset.sum_const, Finset.card_univ,
    nsmul_eq_mul, Fintype.card_fin, ← inv_cost φ, show raw.ret.length = 8 from by decide]
  norm_num [bcost]
  try rfl

/-- the bank schedule of the certificate (tables `Gen.C588MSched`, checker `GX.schedK` of DFT New Sums) -/
def sched : GX.BankSched (Fin raw.R) := GX.BankSched.ofTab 4032 raw.R schedRows schedPad schedRes sched_ok

theorem res_ok : ∀ q : Fin 504, 9 - GLab.cdim (par.tab.get (st0.get (2 * 112 + q.val))) = schedRes.getD q.val 0 := by
  decide +kernel

theorem hent : ∀ q : Fin raw.R, Fintype.card (Fin aux.p.h) - (ginv.ent q).dim = schedRes.getD q.val 0 := fun q => by
  rw [ginv, GXD.ginvOfZ_ent_dim, Fintype.card_fin]
  exact res_ok q

theorem fitsS : sched.fits ginv.ent := GX.BankSched.ofTab_fits _ _ _ _ _ sched_ok ginv.ent hent

theorem width : 64 * Fintype.card (Fin raw.v) + sched.nb = 11200 := by
  rw [Fintype.card_fin, sched, GX.BankSched.ofTab_nb, raw_v]

theorem pad_price (φ : ℕ → ℝ) : ∑ b, bcost φ (sched.pad b) = BridgeRate.C588M.padPrice φ := by
  rw [sched, GX.BankSched.ofTab_padPrice, GXD.rcost_wsumK φ schedPad 28 rankBoundP]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, rankP0_ok, rankP1_ok, rankP2_ok, rankP3_ok, rankP4_ok, rankP5_ok, rankP6_ok, rankP7_ok, rankP8_ok, rankP9_ok, rankP10_ok, rankP11_ok, rankP12_ok, rankP13_ok, rankP14_ok, rankP15_ok, rankP16_ok, rankP17_ok, rankP18_ok, rankP19_ok, rankP20_ok, rankP21_ok, rankP22_ok, rankP23_ok, rankP24_ok, rankP25_ok, rankP26_ok, rankP27_ok, rankP28_ok]
  unfold BridgeRate.C588M.padPrice
  norm_num [bcost]

/-- the unit price identity of the PR #401 unit (`GX.unitPrice401` shape) -/
theorem price (φ : ℕ → ℝ) : 72 * ginv.cost φ + (Fintype.card (Fin raw.v) : ℝ) * GX.bankPrice377 φ
    + ∑ b, bcost φ (sched.pad b) = BridgeRate.C588M.unitPrice φ := by
  rw [ginv_cost, pad_price, Fintype.card_fin, raw_v]
  unfold BridgeRate.C588M.unitPrice BridgeRate.C588M.invPrice invPrice
  push_cast
  ring

end
end GXD.C588M
#print axioms GXD.C588M.ginv_cost
#print axioms GXD.C588M.price
