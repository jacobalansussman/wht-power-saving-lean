import Work.PfLean.Dem3
import Work.PfLean.Dem3Check

/-!
# Block size 3: `14 W / 5`, unconditionally (px-lean, 14 October 2026)

`dem3`: the demand fact `224 ≤ raise (hookN 3) kernel` at `f = 5`, from `Dem3.dem3_of_check` and the
eight kernel-evaluated block checks `Dem3Comp.blkOK`.  `block3_bound`: with pf-lean's
`block3_bound_of_dem` (Theorem R + hook kept + cap 80), every engine word that computes the block-3
kernel on `W` live roles has total block rank at least `14 W / 5`.
-/

namespace OAI
namespace PowerSaving.PfLean
open Binary GF

theorem dem3 : 224 ≤ raise (hookN (α := Fin 3) 3) (matL 5 (kernel (Fin 3))) :=
  Dem3.dem3_of_check Dem3Comp.blkOK

/-- **Block size 3: every unit of the engine game needs total rank at least `14 W / 5`.** -/
theorem block3_bound {ρ σ : Type*} [Fintype ρ] [DecidableEq ρ] [Fintype σ]
    {e : σ → ρ} (he : Function.Injective e) (w : GWord (Fin 3) ρ) (hw : GLiveKernel e w) :
    14 * Fintype.card σ ≤ 5 * w.cost id :=
  block3_bound_of_dem dem3 he w hw

end PowerSaving.PfLean
end OAI

#print axioms OAI.PowerSaving.PfLean.dem3
#print axioms OAI.PowerSaving.PfLean.block3_bound
