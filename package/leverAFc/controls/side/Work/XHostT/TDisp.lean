import Work.XHostT.TCert

/-!
# (key: bridge-host-cross) TWO PROVIDER WORDS on the cross-replica hosted network: the dispatch table as data

The owner layer of #3945 (Dugongue, on dataisfire's #3734 staged544 owner interface) runs each of its 3,240 calls
(360 replicas, 9 stages) with ONE of two provider words, `extra` (1,672 calls) or `gauge44` (1,568 calls).  The hub
gives the dispatch as a table (`owner/inputs/DISPATCH.json`, in the global call order of `ORDER.json`).  On p-a3's
network `XNet` (any replica count, any global order, one invocation `Ic ρ l` per call; `Work/XHostT/TDef.lean`) this
file fixes the per-call invocations to a family of providers and a table:

* `dispIc D disp`: call `(ρ, l)` runs provider `disp (ρ, l)` of the family `D` (lg-a2a4's `XHostD.ProvF`: all
  providers have the same target vectors); `dispOfG tm htm dG`: the same table given by GLOBAL call (`dG τ` is the
  provider of the call at time `τ`, read through the order `tm`);
* the price is summed per call: `sum_calls_disp` (`Σ_i nCalls disp i * cost_i`), `sum_calls_disp2` (two providers),
  `sum_calls_global` / `sum_calls_dispG` (by global call); `XNet.unit_split` (calls + `restU`: the remaining climbs
  of the hosts, the final pair-4 blocks and the pads), `XNet.unit_disp`, `XNet.unit_disp2`, `XNet.unit_dispG`;
* **`xnet_disp_certificate`**: the network's `GCert` with the two-provider price;
* `XNet.withIc` (generated part, `TDispW.lean`): the same network with other per-call invocations of the same
  entrance labels; used by the controls.

Written by p-disp (wht45/lever, 19 Oct 2026) after lg-a2a4's `Work/XHostD/Prov.lean` (wht43, 14 Oct); network by
p-a3 (wht45/lever) after fsp-perrun / h-L3-build on Ryan Shea's complex588 tree (Nil Mamano's gen401); hosting by
DaysSky (PR #442); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHostT
open Binary Matrix Finset RAM SS CB RF BR GF GX XHost XHostD
noncomputable section

section Disp
variable {ι H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

/-- **the per-call invocations of a dispatched network**: call `(ρ, l)` runs provider `disp (ρ, l)` of `D`. -/
def dispIc {nR : ℕ} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) : Fin nR → Fin 9 → GInvE H T Sl C :=
  fun ρ l => D.call disp (ρ, l)

lemma dispIc_tv {nR : ℕ} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) (ρ : Fin nR) (l : Fin 9) :
    (dispIc D disp ρ l).tv = (D.P D.i0).tv := rfl

lemma dispIc_ent {nR : ℕ} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) (ρ : Fin nR) (l : Fin 9) :
    (dispIc D disp ρ l).ent = (D.P (disp (ρ, l))).ent := rfl

lemma dispIc_cost {nR : ℕ} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) (ρ : Fin nR) (l : Fin 9) :
    (dispIc D disp ρ l).cost = (D.P (disp (ρ, l))).cost := rfl

/-- the residual of slot `q` in call `(ρ, l)` is that of the call's provider (what a host check reads). -/
lemma dispIc_res {nR : ℕ} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) (ρ : Fin nR) (l : Fin 9) (q : Sl) :
    Fintype.card H - ((dispIc D disp ρ l).ent q).dim = D.res (disp (ρ, l)) q := rfl

/-- the reference invocation of a dispatched network (its target vectors are every provider's). -/
lemma dispIc_hIc {nR : ℕ} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) :
    ∀ ρ l, (dispIc D disp ρ l).tv = (D.P D.i0).tv := fun _ _ => rfl

/-- a dispatch table indexed by GLOBAL call (`dG τ`: the provider of the call at time `τ`), read through `tm`. -/
def dispOfG {nR : ℕ} (tm : Fin nR → Fin 9 → ℕ) (htm : ∀ ρ l, tm ρ l < 9 * nR) (dG : Fin (9 * nR) → ι) :
    Fin nR × Fin 9 → ι := fun c => dG ⟨tm c.1 c.2, htm c.1 c.2⟩

end Disp

section Price
variable {ι H T Sl C : Type} [Fintype ι] [DecidableEq ι] [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

/-- **the price of the calls under a dispatch table**: every call pays its own provider. -/
theorem sum_calls_disp {nR : ℕ} (D : ProvF ι H T Sl C) (disp : Fin nR × Fin 9 → ι) (φ : ℕ → ℝ) :
    ∑ ρ : Fin nR, ∑ l : Fin 9, (dispIc D disp ρ l).cost φ = ∑ i, (nCalls disp i : ℝ) * (D.P i).cost φ :=
  (sum_runs_stages (fun c : Fin nR × Fin 9 => (D.P (disp c)).cost φ)).trans
    (sum_dispatch (fun i => (D.P i).cost φ) disp)

/-- two providers: `nCalls disp 0` calls of word 0 and `nCalls disp 1` of word 1. -/
theorem sum_calls_disp2 {nR : ℕ} (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (disp : Fin nR × Fin 9 → Fin 2)
    (φ : ℕ → ℝ) :
    ∑ ρ : Fin nR, ∑ l : Fin 9, (dispIc (Prov2.mk' I0 I1 h) disp ρ l).cost φ
      = (nCalls disp 0 : ℝ) * I0.cost φ + (nCalls disp 1 : ℝ) * I1.cost φ := by
  rw [sum_calls_disp, Fin.sum_univ_two]
  rfl

/-- a sum over the calls `(ρ, l)` of a function of their global time is the sum over the global times. -/
theorem sum_calls_global {M : Type*} [AddCommMonoid M] {nR : ℕ} (tm : Fin nR → Fin 9 → ℕ)
    (htm : ∀ ρ l, tm ρ l < 9 * nR) (cl : Fin (9 * nR) → Fin nR × Fin 9) (hbij : Function.Bijective cl)
    (hcl : ∀ τ, tm (cl τ).1 (cl τ).2 = τ.val) (f : Fin (9 * nR) → M) :
    ∑ ρ : Fin nR, ∑ l : Fin 9, f ⟨tm ρ l, htm ρ l⟩ = ∑ τ, f τ := by
  have e2 := Fintype.sum_bijective cl hbij f (fun c : Fin nR × Fin 9 => f ⟨tm c.1 c.2, htm c.1 c.2⟩)
    (fun τ => congrArg f (Fin.ext (hcl τ).symm))
  exact (Fintype.sum_prod_type' (fun ρ l => f ⟨tm ρ l, htm ρ l⟩)).symm.trans e2.symm

/-- **the price of the calls under a table by GLOBAL call**. -/
theorem sum_calls_dispG {nR : ℕ} (D : ProvF ι H T Sl C) (tm : Fin nR → Fin 9 → ℕ)
    (htm : ∀ ρ l, tm ρ l < 9 * nR) (cl : Fin (9 * nR) → Fin nR × Fin 9) (hbij : Function.Bijective cl)
    (hcl : ∀ τ, tm (cl τ).1 (cl τ).2 = τ.val) (dG : Fin (9 * nR) → ι) (φ : ℕ → ℝ) :
    ∑ ρ : Fin nR, ∑ l : Fin 9, (dispIc D (dispOfG tm htm dG) ρ l).cost φ
      = ∑ i, (nCalls dG i : ℝ) * (D.P i).cost φ :=
  (sum_calls_global tm htm cl hbij hcl (fun τ => (D.P (dG τ)).cost φ)).trans
    (sum_dispatch (fun i => (D.P i).cost φ) dG)

/-- the call counts of a table by global call, read through the order, are the table's counts. -/
theorem nCalls_dispOfG {nR : ℕ} (tm : Fin nR → Fin 9 → ℕ) (htm : ∀ ρ l, tm ρ l < 9 * nR)
    (cl : Fin (9 * nR) → Fin nR × Fin 9) (hbij : Function.Bijective cl) (hcl : ∀ τ, tm (cl τ).1 (cl τ).2 = τ.val)
    (dG : Fin (9 * nR) → ι) (i : ι) : nCalls (dispOfG tm htm dG) i = nCalls dG i := by
  unfold nCalls
  rw [Finset.card_filter, Finset.card_filter]
  refine (Fintype.sum_bijective cl hbij (fun τ => if dG τ = i then 1 else 0)
    (fun c => if dispOfG tm htm dG c = i then 1 else 0) (fun τ => ?_)).symm
  have e : (⟨tm (cl τ).1 (cl τ).2, htm (cl τ).1 (cl τ).2⟩ : Fin (9 * nR)) = τ := Fin.ext (hcl τ)
  show (if dG τ = i then 1 else 0) = (if dG ⟨tm (cl τ).1 (cl τ).2, htm (cl τ).1 (cl τ).2⟩ = i then 1 else 0)
  rw [e]

end Price

section Net
variable {H T Sl C α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {nb nR : ℕ} (N : XNet H T Sl C α Γ nb nR)

/-- **everything in the unit price but the calls**: the remaining climbs of the hosts (`tailU`, with the final
pair-4 blocks at stage 8) and one pad per bank.  It does not mention the invocations. -/
def XNet.restU (φ : ℕ → ℝ) : ℝ := (∑ ρ : Fin nR, ∑ l : Fin 9, N.tailU ρ l φ) + ∑ b, bcost φ (N.pad b)

/-- the unit price is the calls plus the rest. -/
theorem XNet.unit_split (φ : ℕ → ℝ) :
    N.unit φ = (∑ ρ : Fin nR, ∑ l : Fin 9, (N.Ic ρ l).cost φ) + N.restU φ := by
  unfold XNet.unit XNet.callU XNet.restU
  simp only [Finset.sum_add_distrib]
  ring

/-- **the unit price of a dispatched network**. -/
theorem XNet.unit_disp {ι : Type} [Fintype ι] [DecidableEq ι] (D : ProvF ι H T Sl C)
    (disp : Fin nR × Fin 9 → ι) (hN : ∀ ρ l, N.Ic ρ l = dispIc D disp ρ l) (φ : ℕ → ℝ) :
    N.unit φ = (∑ i, (nCalls disp i : ℝ) * (D.P i).cost φ) + N.restU φ := by
  rw [N.unit_split φ, ← sum_calls_disp D disp φ]
  simp only [hN]

/-- **the unit price with two providers**: `n0 * cost0 + n1 * cost1 + restU`. -/
theorem XNet.unit_disp2 (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (disp : Fin nR × Fin 9 → Fin 2)
    (hN : ∀ ρ l, N.Ic ρ l = dispIc (Prov2.mk' I0 I1 h) disp ρ l) (φ : ℕ → ℝ) :
    N.unit φ = (nCalls disp 0 : ℝ) * I0.cost φ + (nCalls disp 1 : ℝ) * I1.cost φ + N.restU φ := by
  rw [N.unit_split φ, ← sum_calls_disp2 I0 I1 h disp φ]
  simp only [hN]

/-- the unit price with a table by GLOBAL call (`dG τ`: the provider of the call at time `τ` of `N`'s order). -/
theorem XNet.unit_dispG {ι : Type} [Fintype ι] [DecidableEq ι] (D : ProvF ι H T Sl C) (dG : Fin (9 * nR) → ι)
    (hN : ∀ ρ l, N.Ic ρ l = dispIc D (dispOfG N.tm N.htm_lt dG) ρ l) (φ : ℕ → ℝ) :
    N.unit φ = (∑ i, (nCalls dG i : ℝ) * (D.P i).cost φ) + N.restU φ := by
  rw [N.unit_split φ, ← sum_calls_dispG D N.tm N.htm_lt N.cl N.cl_bij N.hcl dG φ]
  simp only [hN]

/-- **THE NETWORK THEOREM WITH TWO PROVIDERS**: the cross-replica hosted word whose call `(ρ, l)` runs provider
`disp (ρ, l)` (`0` = `I0`, `1` = `I1`) is a scratch certificate, in whole blocks, with the price
`|Γ| * (n0 * cost0 + n1 * cost1 + restU)` for every price list. -/
theorem xnet_disp_certificate (I0 I1 : GInvE H T Sl C) (h : I1.tv = I0.tv) (disp : Fin nR × Fin 9 → Fin 2)
    (hN : ∀ ρ l, N.Ic ρ l = dispIc (Prov2.mk' I0 I1 h) disp ρ l) :
    GCert α (Sum.inl : BLiveT Γ T nb nR → BRoleT Γ T nb C nR)
      (fun φ => (Fintype.card Γ : ℝ) * ((nCalls disp 0 : ℝ) * I0.cost φ + (nCalls disp 1 : ℝ) * I1.cost φ
        + N.restU φ)) :=
  (xnet_certificate N).cast (fun φ => by rw [N.unit_disp2 I0 I1 h disp hN φ])

/-- the same with a family of any size. -/
theorem xnet_dispF_certificate {ι : Type} [Fintype ι] [DecidableEq ι] (D : ProvF ι H T Sl C)
    (disp : Fin nR × Fin 9 → ι) (hN : ∀ ρ l, N.Ic ρ l = dispIc D disp ρ l) :
    GCert α (Sum.inl : BLiveT Γ T nb nR → BRoleT Γ T nb C nR)
      (fun φ => (Fintype.card Γ : ℝ) * ((∑ i, (nCalls disp i : ℝ) * (D.P i).cost φ) + N.restU φ)) :=
  (xnet_certificate N).cast (fun φ => by rw [N.unit_disp D disp hN φ])

end Net
end
end XHostT
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHostT.sum_calls_disp
#print axioms OAI.PowerSaving.XHostT.sum_calls_disp2
#print axioms OAI.PowerSaving.XHostT.sum_calls_dispG
#print axioms OAI.PowerSaving.XHostT.nCalls_dispOfG
#print axioms OAI.PowerSaving.XHostT.XNet.unit_split
#print axioms OAI.PowerSaving.XHostT.XNet.unit_disp2
#print axioms OAI.PowerSaving.XHostT.XNet.unit_dispG
#print axioms OAI.PowerSaving.XHostT.xnet_disp_certificate
#print axioms OAI.PowerSaving.XHostT.xnet_dispF_certificate
