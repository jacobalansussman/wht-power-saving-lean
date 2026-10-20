import Work.XHost.HEnd

/-!
# (key: bridge-host) A bank schedule is a host schedule with no host

`HostSched.ofBank T S`: Shea's / Mamano's bank schedule `S` (Work/Bridge401/Sched.lean) with the empty host table; every
idle climb keeps its full rank (`hrem c j = rankOf c`).  Its residual conditions follow from `S.fits`, and the hosted unit
price is Shea's `unitPrice401` (`8 * Σ_c hcost = card T * bankPrice377`).  This is the end-to-end control of the hosted
network: through `engine_of_ginv_host` it must give exactly the figure of the banked certificate.
Written by h-L3-build (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX B401
noncomputable section

/-- a bank schedule as a host schedule with no host. -/
def HostSched.ofBank {Sl : Type} [Fintype Sl] (T : Type) (S : BankSched Sl) : HostSched Sl T where
  nb := S.nb
  hostAt := fun _ _ => none
  hostInj := fun l q q' j h _ => absurd h (by simp)
  hostWin := fun l q c j h => absurd h (by simp)
  bk := S.use
  hinj := fun ρ l q q' _ _ e => S.hinj ρ l e
  pad := S.pad
  hrem := fun c _ => rankOf c.val

section Fits
variable {Sl T H : Type} [Fintype Sl] [DecidableEq Sl] [DecidableEq T] [Fintype H]

lemma bUse_none_LH (S : BankSched Sl) (r : Sl → ℕ) (b : Fin S.nb) :
    LH (bUse (fun (_ : Fin 9) (_ : Sl) => (none : Option (Fin 9 × (Bool × T)))) S.use) r 72 b = L S.useK r 72 b := by
  unfold LH rrH L rr
  refine Finset.sum_congr rfl (fun k _ => ?_)
  refine Finset.sum_congr ?_ (fun _ _ => rfl)
  ext q
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  show bUse (fun (_ : Fin 9) (_ : Sl) => (none : Option (Fin 9 × (Bool × T)))) S.use k q = some b ↔ S.useK k q = b
  unfold bUse BankSched.useK
  simp

lemma hUse_none_LH (c : Fin 9) (k0 K : ℕ) (r : Sl → ℕ) (j : Bool × T) :
    LH (hUse (fun (_ : Fin 9) (_ : Sl) => (none : Option (Fin 9 × (Bool × T)))) c k0) r K j = 0 := by
  unfold LH rrH
  refine Finset.sum_eq_zero (fun k _ => Finset.sum_eq_zero (fun q hq => ?_))
  have h := (Finset.mem_filter.mp hq).2
  exact absurd h (by unfold hUse; split_ifs <;> simp)

/-- **the residual conditions of a bank schedule, as a host schedule.** -/
theorem HostSched.ofBank_fits (S : BankSched Sl) (ent : Sl → SLbl H) (hres : S.fits ent) :
    (HostSched.ofBank T S).fits ent := by
  unfold HostSched.fits HostSched.fitsR
  refine ⟨fun (b : Fin S.nb) => ?_, hres.2, fun c j => ?_⟩
  · show LH (bUse (fun (_ : Fin 9) (_ : Sl) => (none : Option (Fin 9 × (Bool × T)))) S.use)
      (fun q => Fintype.card H - (ent q).dim) 72 b + S.pad b = 73
    rw [bUse_none_LH, S.L_eq]
    exact hres.1 b
  · show LH (hUse (fun (_ : Fin 9) (_ : Sl) => (none : Option (Fin 9 × (Bool × T)))) c (kz c))
      (fun q => Fintype.card H - (ent q).dim) (kK c) j + rankOf c.val = rankOf c.val
    rw [hUse_none_LH, Nat.zero_add]

end Fits

section Price
variable {Sl T : Type} [Fintype Sl] [Fintype T]

lemma hcostS_ofBank (S : BankSched Sl) (c : Fin 9) (φ : ℕ → ℝ) :
    hcostS (HostSched.ofBank T S) c φ = (Fintype.card T : ℝ) * (2 * bcost φ (rankOf c.val)) := by
  unfold hcostS
  simp only [HostSched.ofBank, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  ring

/-- **the remaining climbs of the empty host table cost what Shea's idle climbs cost.** -/
theorem hcost_ofBank (S : BankSched Sl) (φ : ℕ → ℝ) :
    8 * (hcostS (HostSched.ofBank T S) 0 φ + hcostS (HostSched.ofBank T S) 1 φ + hcostS (HostSched.ofBank T S) 2 φ
      + hcostS (HostSched.ofBank T S) 3 φ + hcostS (HostSched.ofBank T S) 4 φ + hcostS (HostSched.ofBank T S) 5 φ
      + hcostS (HostSched.ofBank T S) 6 φ + hcostS (HostSched.ofBank T S) 7 φ + hcostS (HostSched.ofBank T S) 8 φ)
      = (Fintype.card T : ℝ) * bankPrice377 φ := by
  simp only [hcostS_ofBank]
  have r0 : rankOf (0 : Fin 9).val = 16 := rfl
  have r1 : rankOf (1 : Fin 9).val = 8 := rfl
  have r2 : rankOf (2 : Fin 9).val = 32 := rfl
  have r3 : rankOf (3 : Fin 9).val = 8 := rfl
  have r4 : rankOf (4 : Fin 9).val = 48 := rfl
  have r5 : rankOf (5 : Fin 9).val = 8 := rfl
  have r6 : rankOf (6 : Fin 9).val = 48 := rfl
  have r7 : rankOf (7 : Fin 9).val = 32 := rfl
  have r8 : rankOf (8 : Fin 9).val = 16 := rfl
  rw [r0, r1, r2, r3, r4, r5, r6, r7, r8]
  unfold bankPrice377
  ring

end Price
end
end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.HostSched.ofBank_fits
#print axioms OAI.PowerSaving.XHost.hcost_ofBank
