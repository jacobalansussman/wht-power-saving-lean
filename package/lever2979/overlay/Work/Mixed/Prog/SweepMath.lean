import Work.Mixed.Prog.Acc

/-!
# Mixed engine, R3 program part 13: facts about the routing (agent e-r3prog)

* `slotR_spec`    every parent unit `p` feeds exactly the child slot `slotR p` (separation is used here);
* `ownSrc_len`, `ownSrc_getD`   the source slots in order;
* `cjob_src`      a source slot carries the supply vector at the running offset plus its rank;
* `stepV_block`   a block letter on the class state is `ferryM` of its splitting's kernel.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)

lemma tq_ge (a : ℕ) (ha : a < E.W) (r : ℕ) (hr : ∀ l, r < E.tr l) : r < E.tq a := by
  unfold tq; rw [dif_pos ha]; exact hr _

/-- **Every parent unit feeds exactly one child slot.** -/
theorem slotR_spec (n h μ r : ℕ) (hr : ∀ l, r < E.tr l) (p : Fin (E.NUc n (h+1) μ)) :
    E.slotR n h μ r p < E.NUc n h (μ * r) * E.W ∧ E.flR n h μ r (E.slotR n h μ r p) ≠ 0 ∧
      E.ppR n h μ r (E.slotR n h μ r p) = p ∧
      E.slotSize n h (μ * r) (E.slotR n h μ r p) = r * E.cUc n (h+1) μ p := by
  obtain ⟨hreg, hlt, hpc⟩ := MXR.U_spec p.2 (W:=E.W) (tq:=E.tq) (mx:=E.mx) (n:=n) (h:=h+1) (ρ:=μ)
  set u := E.Uc n (h+1) μ p with hu
  have hW := E.W_pos
  set τ := u / E.W
  set a0 := u % E.W
  have ha0 : a0 < E.W := Nat.mod_lt _ hW
  have hu' : a0 + E.W * τ = u := by rw [Nat.add_comm]; exact Nat.div_add_mod u E.W
  have hτ : τ < E.W^h := by rw [Nat.div_lt_iff_lt_mul hW, ← pow_succ]; exact hlt
  -- the child unit τ is in the child's region (separation)
  have hcol : colS E.W E.tq (h+1) u = E.tq a0 * colS E.W E.tq h τ := by
    rw [← hu']; exact colS_cons E.W E.tq h a0 τ ha0
  have hregc : inRS E.W E.tq E.mx n h (μ * r) τ := by
    refine ⟨hτ, ?_⟩
    have h1 := hreg.2
    rw [hcol] at h1
    have h2 : r ≤ E.tq a0 := (E.tq_ge a0 ha0 r hr).le
    calc E.mx * (μ * r * colS E.W E.tq h τ) = E.mx * (μ * (r * colS E.W E.tq h τ)) := by ring
      _ ≤ E.mx * (μ * (E.tq a0 * colS E.W E.tq h τ)) := by gcongr
      _ ≤ n := h1
  have hpc' := PC_lt hregc
  have hUPC : E.Uc n h (μ * r) (E.PCc n h (μ * r) τ) = τ := U_PC hregc
  have hslot : E.slotR n h μ r p = E.PCc n h (μ * r) τ * E.W + a0 := rfl
  have e1 : E.slotR n h μ r p % E.W = a0 := by
    rw [hslot, Nat.add_comm, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt ha0]
  have e2 : E.slotR n h μ r p / E.W = E.PCc n h (μ * r) τ := by
    rw [hslot, Nat.add_comm, Nat.add_mul_div_right _ _ hW, Nat.div_eq_of_lt ha0, zero_add]
  have hlt' : E.slotR n h μ r p < E.NUc n h (μ * r) * E.W := by
    rw [hslot]
    calc E.PCc n h (μ * r) τ * E.W + a0 < E.PCc n h (μ * r) τ * E.W + E.W := by omega
      _ = (E.PCc n h (μ * r) τ + 1) * E.W := by ring
      _ ≤ E.NUc n h (μ * r) * E.W := Nat.mul_le_mul_right _ hpc'
  have hpar : E.parCode n h μ r (E.slotR n h μ r p) = u := by
    unfold parCode; rw [e1, e2]
    show a0 + E.W * E.Uc n h (μ * r) (E.PCc n h (μ * r) τ) = u
    rw [hUPC, hu']
  refine ⟨hlt', ?_, ?_, ?_⟩
  · unfold flR; rw [if_pos (show E.fedP n h μ r _ from ⟨hlt', by rw [hpar]; exact hreg⟩)]; exact one_ne_zero
  · unfold ppR; rw [hpar]; exact hpc
  · unfold slotSize cUc
    rw [e1, e2]
    show E.tq a0 * ((μ * r) * colS E.W E.tq h (E.Uc n h (μ * r) (E.PCc n h (μ * r) τ))) =
      r * (μ * colS E.W E.tq (h+1) u)
    rw [hUPC, hcol]; ring

lemma ownSrc_eq (n h μ r : ℕ) :
    E.ownSrc n h μ r = (List.range (Nat.count (E.srcP n h μ r) (E.NUc n h (μ * r) * E.W))).map
      (fun k => E.slotSize n h (μ * r) (Nat.nth (E.srcP n h μ r) k)) := by
  unfold ownSrc
  have : (List.range (E.NUc n h (μ * r) * E.W)).filterMap
      (fun t => if E.srcP n h μ r t then some (E.slotSize n h (μ * r) t) else none) =
      (List.range (Nat.count (E.srcP n h μ r) (E.NUc n h (μ * r) * E.W))).map
      (fun k => E.slotSize n h (μ * r) (Nat.nth (E.srcP n h μ r) k)) :=
    filterMap_range_eq (E.srcP n h μ r) (E.slotSize n h (μ * r)) _
  exact this

lemma ownSrc_len (n h μ r : ℕ) :
    (E.ownSrc n h μ r).length = Nat.count (E.srcP n h μ r) (E.NUc n h (μ * r) * E.W) := by
  rw [ownSrc_eq]; simp

/-- Source slots are exactly the slots that are not fed. -/
lemma fl_zero_iff (n h μ r t : ℕ) (ht : t < E.NUc n h (μ * r) * E.W) :
    E.flR n h μ r t = 0 ↔ E.srcP n h μ r t := by
  unfold flR
  constructor
  · intro h0
    refine ⟨ht, fun hin => ?_⟩
    rw [if_pos (show E.fedP n h μ r t from ⟨ht, hin⟩)] at h0
    exact one_ne_zero h0
  · intro hs
    rw [if_neg (show ¬ E.fedP n h μ r t from fun hh => hs.2 hh.2)]

/-- **A block letter on the class state is `ferryM` of the kernel on its splitting.** -/
lemma stepV_block {n h μ : ℕ} (V : E.CState n h μ) (l : ρ) (r : ℕ) (z : Fin r → Space α)
    (ind : LinearIndependent F z) :
    E.stepV V (.old (.block l r z ind)) =
      fun (p : Fin (E.NUc n h μ)) q => ferryM 0 (E.cUc n h μ p) l ((bS z ind).lift (kernel (Fin r))) (V p q) := by
  funext p q x
  have hM : blockMat z = (bS z ind).lift (kernel (Fin r)) :=
    blockMat_eq_lift (bS z ind) z (Classical.choose_spec (exists_split z ind))
  unfold stepV gferry ferryM
  simp only [GMove.act, BMove.act]
  rw [hM]

end Eng

end
end PowerSaving.RAM
end OAI
