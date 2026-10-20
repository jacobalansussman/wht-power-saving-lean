import Work.XHostT.Inst.S2979
import Work.XHostT.TEnd
import Work.XHostT.KHead
import Work.XHostT.Data.H2979.Glob
import Work.XHostT.Data.H2979.Heads
import Work.XHostT.Data.H2979.DispLink

/-!
# (key: bridge-host-cross) The data facts of #2979's schedule that the bridge needs

For `sched2979` (Inst/S2979.lean, the cross-replica host schedule of hub PR #2979):
* **`sched2979_heads`**: 16,086 hosts of kind A4 have remainder 48 (the literal heads: with positive residuals, p-hfirst's
  `headBins_eq` makes them exactly `headBins sched2979`);
* **`sched2979_hcost`**: the remaining climbs of all host chains cost `17310 φ1 + 4956 φ3 + 1728 φ5 + 28190 φ48`
  (the idle part of #2979's PROFILE.json; remainder-0 hosts are free);
* **`resK2979_pos`** (every residual is positive, `headBins_eq`'s `hpos`), `resK2979_ge4`, `resK2979_eq0` (the residual
  of a slot does not depend on the call's provider), `resK2979_disp` (the residual read is the dispatched provider's);
* **`card_live2979`**: `66,304 = 74 x 896` live roles per gauge copy (`hσ`).
x-tables (wht45/lever2979, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI.PowerSaving.XHostT.Inst
open Finset OAI.PowerSaving.XHostT OAI.PowerSaving.XHostT.K OAI.PowerSaving.XHostT.K.DH2979
open Binary Matrix RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD

/-- **the head count of #2979**: 16,086 hosts of kind A4 (both sides) have remainder 48. -/
theorem sched2979_heads :
    #{p ∈ (univ : Finset (Fin 74 × Bool × Fin 112)) | sched2979.hrem p.1 4 p.2 = 48} = 16086 := by
  show #{p ∈ (univ : Finset (Fin 74 × Bool × Fin 112)) |
    remAt rt (2016 * p.1.val + 224 * (4 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 16086
  rw [← cntK_card]
  exact heads_cnt

/-- the hosts of kind F1 with remainder 48 (unhosted, not heads). -/
theorem sched2979_f1 :
    #{p ∈ (univ : Finset (Fin 74 × Bool × Fin 112)) | sched2979.hrem p.1 6 p.2 = 48} = 12104 := by
  show #{p ∈ (univ : Finset (Fin 74 × Bool × Fin 112)) |
    remAt rt (2016 * p.1.val + 224 * (6 : Fin 9).val + (if p.2.1 then 112 else 0) + p.2.2.val) = 48} = 12104
  rw [← cntK_card]
  exact f1_cnt

/-- **the price of the remaining climbs of #2979's host chains.** -/
theorem sched2979_hcost (φ : ℕ → ℝ) :
    ∑ ρ : Fin 74, ∑ c : Fin 9, hcostT sched2979 ρ c φ
      = 17310 * φ 1 + 4956 * φ 3 + 1728 * φ 5 + 28190 * φ 48 := by
  have h := rem_hist (M := ℝ) (bcost φ)
  rw [show 224 * NK = 224 * (9 * 74) from rfl, sum_bins] at h
  have e : ∀ ρ : Fin 74, ∀ c : Fin 9, hcostT sched2979 ρ c φ = ∑ t : Fin 112,
      (bcost φ (remAt rt (2016 * ρ.val + 224 * c.val + t.val))
        + bcost φ (remAt rt (2016 * ρ.val + 224 * c.val + 112 + t.val))) := by
    intro ρ c
    unfold hcostT
    refine Fintype.sum_congr _ _ (fun t => ?_)
    rw [sched2979_hrem, sched2979_hrem]
    simp only [Bool.false_eq_true, ↓reduceIte, Nat.add_zero]
  simp only [e]
  rw [h]
  simp only [bcost, nsmul_eq_mul]
  norm_num

/-- the residuals of the table's two providers. -/
lemma resAt2979_ge (p : Fin 2) (q : Fin 408) : 4 ≤ resAt resT p.val q.val := by
  rcases p with ⟨p, hp⟩
  interval_cases p
  · exact resGE_sound resT 0 4 408 res_ge0 q.val q.isLt
  · exact resGE_sound resT 1 4 408 res_ge1 q.val q.isLt

/-- the residual of slot `q` of call `(ρ, l)` is the dispatched provider's. -/
theorem resK2979_disp (ρ : Fin 74) (l : Fin 9) (q : Fin 408) :
    resK ct resT 74 408 ρ l q = resAt resT (D2979.disp (ρ, l)).val q.val := by
  unfold resK
  rw [dp_disp]

theorem resK2979_ge4 (ρ : Fin 74) (l : Fin 9) (q : Fin 408) : 4 ≤ resK ct resT 74 408 ρ l q := by
  rw [resK2979_disp]
  exact resAt2979_ge _ q

/-- **every residual of the table is positive** (`headBins_eq`'s `hpos`). -/
theorem resK2979_pos (ρ : Fin 74) (l : Fin 9) (q : Fin 408) : 0 < resK ct resT 74 408 ρ l q :=
  lt_of_lt_of_le (by norm_num) (resK2979_ge4 ρ l q)

/-- the residual of a slot does not depend on the provider (the two residual vectors of #2979 are equal). -/
theorem resK2979_eq0 (ρ : Fin 74) (l : Fin 9) (q : Fin 408) : resK ct resT 74 408 ρ l q = resAt resT 0 q.val := by
  rw [resK2979_disp]
  rcases hd : D2979.disp (ρ, l) with ⟨p, hp⟩
  interval_cases p
  · rfl
  · exact (resEqP_sound resT 408 res_eq01 q.val q.isLt).symm

/-- **the live roles of #2979's network**: `74 x 8 x 112 = 66,304` per gauge copy. -/
theorem card_live2979 (Γ : Type) [Fintype Γ] :
    Fintype.card (BLiveT Γ (Fin 112) 0 74) = Fintype.card Γ * 66304 := by
  rw [card_BLiveT, Fintype.card_fin]

end OAI.PowerSaving.XHostT.Inst

#print axioms OAI.PowerSaving.XHostT.Inst.sched2979_heads
#print axioms OAI.PowerSaving.XHostT.Inst.sched2979_f1
#print axioms OAI.PowerSaving.XHostT.Inst.sched2979_hcost
#print axioms OAI.PowerSaving.XHostT.Inst.resK2979_pos
#print axioms OAI.PowerSaving.XHostT.Inst.resK2979_eq0
#print axioms OAI.PowerSaving.XHostT.Inst.card_live2979
