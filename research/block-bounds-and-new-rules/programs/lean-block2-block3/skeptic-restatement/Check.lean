import Work.PfLean.Dem3Final
-- v-px-lean-logic (skeptic): restate the claimed theorem in my own words, with fully qualified names,
-- and close it with px-lean's term. If the statement differed from this, elaboration would fail.
open OAI.PowerSaving in
theorem skeptic_block3 :
    ∀ {ρ : Type} {σ : Type} [Fintype ρ] [DecidableEq ρ] [Fintype σ] (e : σ → ρ),
      Function.Injective e →
      ∀ w : GF.GWord (Fin 3) ρ, GF.GLiveKernel e w →
        14 * Fintype.card σ ≤ 5 * GF.GWord.cost id w :=
  fun e he w hw => OAI.PowerSaving.PfLean.block3_bound he w hw

open OAI.PowerSaving in
theorem skeptic_dem3 :
    224 ≤ PfLean.raise (PfLean.hookN (α := Fin 3) 3) (PfLean.matL 5 (Binary.kernel (Fin 3))) :=
  PfLean.dem3

-- the definitions the statement rests on, printed from the environment
#print OAI.PowerSaving.GF.GLiveKernel
#print OAI.PowerSaving.GF.gwalk
#print OAI.PowerSaving.GF.GMove.act
#print OAI.PowerSaving.Binary.BMove.act
#print OAI.PowerSaving.Binary.Data
#print OAI.PowerSaving.Binary.roleAct
#print OAI.PowerSaving.Binary.kernel
#print OAI.PowerSaving.Binary.littleC
#print OAI.PowerSaving.Binary.ca
#print OAI.PowerSaving.Binary.cb
#print OAI.PowerSaving.Binary.matAct
#print OAI.PowerSaving.PfLean.raise
#print OAI.PowerSaving.PfLean.matL
#print OAI.PowerSaving.PfLean.hookN
#print OAI.PowerSaving.PfLean.Dsub
#print OAI.PowerSaving.PfLean.Alt
#print OAI.PowerSaving.PfLean.Xf
#print axioms skeptic_block3
#print axioms skeptic_dem3
#print axioms OAI.PowerSaving.PfLean.Dem3Comp.blk0
#print axioms OAI.PowerSaving.PfLean.Dem3Comp.blk7
#print axioms OAI.PowerSaving.PfLean.Dem3.block_det
#print axioms OAI.PowerSaving.PfLean.Dem3.data_ok
