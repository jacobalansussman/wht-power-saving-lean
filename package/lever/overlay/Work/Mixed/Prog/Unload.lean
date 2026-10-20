import Work.Mixed.Prog.RouteTabs

/-!
# Mixed engine, R3 program part 10: loading a class from its jobs, and answering them (agent e-r3prog)

* `Eng.unl`         own slot `t = p*W + a` read from unit `p`'s live role `a` through the target placement;
* `Eng.unl_run`     **correctness of a class**: after the whole word, `unl` is `hlow` of the slot's job on its
                    `tq a * cU p` H-bits (the mixed certificate + `hlow_placeS`);
* `CCtx.load`, `CCtx.unload`   the two RAM passes (linear work).
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

universe U

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)

/-- The live role at position `t % W`. -/
def roleAt (t : ℕ) : σ := E.Q.atFin ⟨t % E.W, Nat.mod_lt _ E.W_pos⟩

lemma loc_roleAt (t : ℕ) : E.Q.loc (E.roleAt t) = t % E.W := by
  unfold roleAt; rw [Layout.loc_atFin]

lemma tq_roleAt (t : ℕ) : E.tq (t % E.W) = E.tr (E.roleAt t) := by
  unfold tq roleAt; rw [dif_pos (Nat.mod_lt _ E.W_pos)]

/-- Own slot `t` read from the class state. -/
def unl (n h μ : ℕ) (V : E.CState n h μ) (t : ℕ) (e : Bits n) : ℂ :=
  if hp : t / E.W < E.NUc n h μ then
    V ⟨t / E.W, hp⟩ ((E.placeT n h μ ⟨t / E.W, hp⟩ (E.roleAt t)).symm e).1
      (E.e (E.roleAt t), (0, ((E.placeT n h μ ⟨t / E.W, hp⟩ (E.roleAt t)).symm e).2))
  else 0

lemma slot_decomp (t : ℕ) : t / E.W * E.W + E.Q.loc (E.roleAt t) = t := by
  rw [loc_roleAt, Nat.mul_comm]; exact Nat.div_add_mod t E.W

/-- **Correctness of a class**: after the word, every own slot holds `hlow` of its job. -/
theorem unl_run (n h μ : ℕ) (J : ℕ → Bits n → ℂ) (t : ℕ) (ht : t < E.NUc n h μ * E.W) :
    E.unl n h μ (E.runV (E.V0 n h μ J) E.wd) t = hlow n (E.slotSize n h μ t) (J t) := by
  funext e
  have hp : t / E.W < E.NUc n h μ := (Nat.div_lt_iff_lt_mul E.W_pos).2 ht
  unfold unl
  rw [dif_pos hp]
  obtain ⟨⟨q, y⟩, hqy⟩ : ∃ z, E.placeT n h μ ⟨t / E.W, hp⟩ (E.roleAt t) z = e := ⟨_, Equiv.apply_symm_apply _ e⟩
  subst hqy
  rw [Equiv.symm_apply_apply]
  -- the word on one instance
  have hv : ∀ s, (∀ l', E.e l' ≠ s) →
      (fun i j => E.V0 n h μ J ⟨t / E.W, hp⟩ q (i, (0 : Bits 0), j) : Data α ρ (E.cUc n h μ (t / E.W))) s = 0 := by
    intro s hs; funext j
    unfold V0
    dsimp only
    rw [E.liveOf_none s hs]
    rfl
  have key := E.hw _ _ hv (E.roleAt t)
  change gwalk (E.cUc n h μ (t / E.W)) E.wd (fun i j => E.V0 n h μ J ⟨t / E.W, hp⟩ q (i, (0 : Bits 0), j))
    (E.e (E.roleAt t)) y = _
  rw [key]
  have hT : blockMat (E.tz (E.roleAt t)) = (E.St (E.roleAt t)).lift (kernel (Fin (E.tr (E.roleAt t)))) :=
    blockMat_eq_lift (E.St (E.roleAt t)) (E.tz (E.roleAt t)) (E.St_spec (E.roleAt t))
  rw [hT]
  have hJ : (fun j => E.V0 n h μ J ⟨t / E.W, hp⟩ q (E.e (E.roleAt t), (0 : Bits 0), j)) =
      fun j => J t (E.placeT n h μ ⟨t / E.W, hp⟩ (E.roleAt t) (q, j)) := by
    funext j
    unfold V0
    dsimp only
    rw [E.liveOf_e (E.roleAt t)]
    dsimp only
    rw [E.slot_decomp t]
  rw [hJ]
  have hsz : E.slotSize n h μ t = E.tr (E.roleAt t) * E.cUc n h μ (t / E.W) := by
    unfold slotSize; rw [E.tq_roleAt]
  rw [hsz]
  exact (hlow_placeS (E.St (E.roleAt t)) (E.Xt (E.roleAt t)) _ _ n _ (J t) q y).symm

end Eng

end
end PowerSaving.RAM
end OAI
