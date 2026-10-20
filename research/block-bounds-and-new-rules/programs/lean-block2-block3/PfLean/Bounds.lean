import Work.PfLean.Dem2

/-!
# Lower bounds for units of the engine game (pf-lean, 13 October 2026)

`w.cost id` = total block rank of the word (a block of rank `r` counts `r` moves).  The live
roles are `e : σ → ρ` (injective, `W = card σ`); any number of scratch roles; free letters:
`localM` (any rational mixing), `shift`, `perm` (address permutations), `phase`.

* `hook_bound`      : for every block size and every `j`:
                      `W * raise(hook, kernel) ≤ (2^m / 2) * C(2^m - 2, j) * total rank`.
* `block2_bound`    : block size 2: `total rank ≥ 2 W`  (fully proved, no hypothesis).
* `block3_bound_of_dem` : block size 3: `5 * total rank ≥ 14 W`, GIVEN the demand fact
                      `224 ≤ raise (hookN 3) (kernel)` at `f = 5` (stated as a hypothesis).
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.PfLean
open Binary GF Module

/-- Theorem R with the hook certificate, any block size `|α|` and any `j` (`f = j + 2`). -/
theorem hook_bound {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
    [Fintype σ] (j : ℕ) {e : σ → ρ} (he : Function.Injective e) (w : GWord α ρ)
    (hw : GLiveKernel e w) :
    Fintype.card σ * raise (hookN (α := α) j) (matL (j+2) (kernel α)) ≤
      ((Fintype.card (Space α) / 2) * (Fintype.card (Space α) - 2).choose j) * w.cost id :=
  engine_lower_bound (hookN j) he (hook_kept j) _ (fun z hz => hook_cap j z hz) w hw

lemma card_space (m : ℕ) : Fintype.card (Space (Fin m)) = 2 ^ m := by
  simp [Space, Fintype.card_fun, ZMod.card]

/-- **Block size 2: every unit needs at least `2 W` moves** (total block rank). -/
theorem block2_bound {ρ σ : Type*} [Fintype ρ] [DecidableEq ρ] [Fintype σ]
    {e : σ → ρ} (he : Function.Injective e) (w : GWord (Fin 2) ρ) (hw : GLiveKernel e w) :
    2 * Fintype.card σ ≤ w.cost id := by
  have h := hook_bound 1 he w hw
  rw [card_space] at h
  have hc : (2 ^ 2 / 2) * (2 ^ 2 - 2).choose 1 = 4 := by decide
  rw [hc] at h
  have h2 := Nat.mul_le_mul_left (Fintype.card σ)
    (show 8 ≤ raise (hookN (α := Fin 2) 1) (matL (1+2) (kernel (Fin 2))) from dem2)
  omega

/-- **Block size 3, conditional on the demand fact**: `5 * (total rank) ≥ 14 W`. -/
theorem block3_bound_of_dem
    (hdem : 224 ≤ raise (hookN (α := Fin 3) 3) (matL 5 (kernel (Fin 3))))
    {ρ σ : Type*} [Fintype ρ] [DecidableEq ρ] [Fintype σ]
    {e : σ → ρ} (he : Function.Injective e) (w : GWord (Fin 3) ρ) (hw : GLiveKernel e w) :
    14 * Fintype.card σ ≤ 5 * w.cost id := by
  have h := hook_bound 3 he w hw
  rw [card_space] at h
  have hc : (2 ^ 3 / 2) * (2 ^ 3 - 2).choose 3 = 80 := by decide
  rw [hc] at h
  have h2 := Nat.mul_le_mul_left (Fintype.card σ)
    (show 224 ≤ raise (hookN (α := Fin 3) 3) (matL (3+2) (kernel (Fin 3))) from hdem)
  omega

/-- The cap at block size 3 (`f = 5`): every direction raises the hook by at most `80`. -/
theorem block3_cap (z : Space (Fin 3)) (hz : z ≠ 0) :
    raise (hookN (α := Fin 3) 3) (matL 5 (dir z)) ≤ 80 := by
  have h := hook_cap (α := Fin 3) 3 z hz
  rw [card_space] at h
  exact h

end PowerSaving.PfLean
end OAI

#print axioms OAI.PowerSaving.PfLean.hook_bound
#print axioms OAI.PowerSaving.PfLean.block2_bound
#print axioms OAI.PowerSaving.PfLean.block3_bound_of_dem
#print axioms OAI.PowerSaving.PfLean.block3_cap
