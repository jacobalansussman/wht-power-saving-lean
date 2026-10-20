import Work.XHostT.Inst.S5
import Work.XHostT.TEnd
import Work.XHostT.KHead
import Work.XHostT.Data.S5.Glob
import Work.XHostT.Data.S5.Heads
import Work.XHostT.Data.S5.DispLink

/-!
# (key: bridge-host-cross) The data facts of #3945's schedule that the bridge needs

For p-a3's `sched3945` (Inst/S5.lean, the cross-replica host schedule of hub PR #3945 from p-hostck's checked table):
* **`sched3945_heads`**: 57,312 hosts of kind A4 have remainder 48 (the literal heads: with positive residuals,
  p-hfirst's `headBins_eq` makes them exactly `headBins sched3945`);
* **`sched3945_hcost`**: the remaining climbs of all host chains cost `23328 φ3 + 4320 φ5 + 137952 φ48` (the idle part
  of #3945's profile; remainder-0 hosts are free), from p-hostck's `rem_hist`;
* **`resK3945_pos`** (every residual is positive, `headBins_eq`'s `hpos`), `resK3945_ge5`, `resK3945_eq0`,
  `resK3945_disp` (the residual read is the dispatched provider's, p-disp's `dp_disp`);
* **`card_live3945`**: `322,560 = 360 x 896` live roles per gauge copy (`hσ`).
The proofs are x-tables' (wht45/lever2979, Inst/S2979Facts.lean) with the instance changed.
i-final (wht45/lever, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.XHostT.Inst
open Finset OAI.PowerSaving.XHostT OAI.PowerSaving.XHostT.K OAI.PowerSaving.XHostT.K.DS5
open Binary Matrix RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD

/-- **the head count of #3945**: 57,312 hosts of kind A4 (both sides) have remainder 48. -/
theorem sched3945_heads :
    #{p ∈ (univ : Finset (Fin 360 × Bool × Fin 112)) | sched3945.hrem p.1 4 p.2 = 48} = 57312 := by
  show #{p ∈ (univ : Finset (Fin 360 × Bool × Fin 112)) |
    remAt rt (2016 * p.1.val + 224 * (4 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 57312
  rw [← cntK_card]
  exact heads_cnt

/-- the hosts of kind F1 with remainder 48 (unhosted, not heads). -/
theorem sched3945_f1 :
    #{p ∈ (univ : Finset (Fin 360 × Bool × Fin 112)) | sched3945.hrem p.1 6 p.2 = 48} = 80640 := by
  show #{p ∈ (univ : Finset (Fin 360 × Bool × Fin 112)) |
    remAt rt (2016 * p.1.val + 224 * (6 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 80640
  rw [← cntK_card]
  exact f1_cnt

/-- **the price of the remaining climbs of #3945's host chains.** -/
theorem sched3945_hcost (φ : ℕ → ℝ) :
    ∑ ρ : Fin 360, ∑ c : Fin 9, hcostT sched3945 ρ c φ
      = 23328 * φ 3 + 4320 * φ 5 + 137952 * φ 48 := by
  have h := rem_hist (M := ℝ) (bcost φ)
  rw [show 224 * NK = 224 * (9 * 360) from rfl, sum_bins] at h
  have e : ∀ ρ : Fin 360, ∀ c : Fin 9, hcostT sched3945 ρ c φ = ∑ t : Fin 112,
      (bcost φ (remAt rt (2016 * ρ.val + 224 * c.val + t.val))
        + bcost φ (remAt rt (2016 * ρ.val + 224 * c.val + 112 + t.val))) := by
    intro ρ c
    unfold hcostT
    refine Fintype.sum_congr _ _ (fun t => ?_)
    rw [sched3945_hrem, sched3945_hrem]
    simp only [Bool.false_eq_true, ↓reduceIte, Nat.add_zero]
  simp only [e]
  rw [h]
  simp only [bcost, nsmul_eq_mul]
  norm_num

/-- the residuals of the table's two providers. -/
lemma resAt3945_ge (p : Fin 2) (q : Fin 408) : 5 ≤ resAt resT p.val q.val := by
  rcases p with ⟨p, hp⟩
  interval_cases p
  · exact resGE_sound resT 0 5 408 res_ge0 q.val q.isLt
  · exact resGE_sound resT 1 5 408 res_ge1 q.val q.isLt

/-- the residual of slot `q` of call `(ρ, l)` is the dispatched provider's. -/
theorem resK3945_disp (ρ : Fin 360) (l : Fin 9) (q : Fin 408) :
    resK ct resT 360 408 ρ l q = resAt resT (D3945.disp (ρ, l)).val q.val := by
  unfold resK
  rw [dp_disp]

theorem resK3945_ge5 (ρ : Fin 360) (l : Fin 9) (q : Fin 408) : 5 ≤ resK ct resT 360 408 ρ l q := by
  rw [resK3945_disp]
  exact resAt3945_ge _ q

/-- **every residual of the table is positive** (`headBins_eq`'s `hpos`). -/
theorem resK3945_pos (ρ : Fin 360) (l : Fin 9) (q : Fin 408) : 0 < resK ct resT 360 408 ρ l q :=
  lt_of_lt_of_le (by norm_num) (resK3945_ge5 ρ l q)

/-- the residual of a slot does not depend on the provider (the two residual vectors of #3945 are equal). -/
theorem resK3945_eq0 (ρ : Fin 360) (l : Fin 9) (q : Fin 408) : resK ct resT 360 408 ρ l q = resAt resT 0 q.val := by
  rw [resK3945_disp]
  rcases hd : D3945.disp (ρ, l) with ⟨p, hp⟩
  interval_cases p
  · rfl
  · exact (resEqP_sound resT 408 res_eq01 q.val q.isLt).symm

/-- **the live roles of #3945's network**: `360 x 8 x 112 = 322,560` per gauge copy. -/
theorem card_live3945 (Γ : Type) [Fintype Γ] :
    Fintype.card (BLiveT Γ (Fin 112) 0 360) = Fintype.card Γ * 322560 := by
  rw [card_BLiveT, Fintype.card_fin]

end OAI.PowerSaving.XHostT.Inst

#print axioms OAI.PowerSaving.XHostT.Inst.sched3945_heads
#print axioms OAI.PowerSaving.XHostT.Inst.sched3945_f1
#print axioms OAI.PowerSaving.XHostT.Inst.sched3945_hcost
#print axioms OAI.PowerSaving.XHostT.Inst.resK3945_pos
#print axioms OAI.PowerSaving.XHostT.Inst.resK3945_eq0
#print axioms OAI.PowerSaving.XHostT.Inst.card_live3945
