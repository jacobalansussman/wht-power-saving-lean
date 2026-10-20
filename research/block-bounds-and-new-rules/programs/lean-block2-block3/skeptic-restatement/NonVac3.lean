import Work.PfLean.Dem3Final
open OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.GF
-- v-px-lean-logic: non-vacuity at block 3. The one-role word with one rank-3 block (the three unit
-- directions) satisfies GLiveKernel, costs 3, and block3_bound gives 14 <= 15 for it: the hypothesis is
-- satisfiable and the bound is not contradicted by the trivial unit.
lemma eu3_li : LinearIndependent (ZMod 2) (fun i : Fin 3 => (eu i : Space (Fin 3))) := by
  have h : (fun i : Fin 3 => (eu i : Space (Fin 3))) = fun i => Pi.basisFun (ZMod 2) (Fin 3) i := by
    funext i j; simp [eu, Pi.basisFun_apply, Pi.single_apply, eq_comm]
  rw [h]; exact (Pi.basisFun (ZMod 2) (Fin 3)).linearIndependent

noncomputable def w3 : GWord (Fin 3) Unit := [GMove.old (BMove.block () 3 (fun i => eu i) eu3_li)]

theorem w3_unit : GLiveKernel (id : Unit → Unit) w3 := by
  intro f v hv l
  cases l
  simp [w3, gwalk, GMove.act, BMove.act, roleAct, kernel_eq_blockMat]

theorem w3_cost : w3.cost id = 3 := by
  simp [w3, GWord.cost, GWord.shadow, BWord.cost, BMove.cost, BMove.rk?, GMove.shadow]

example : 14 * Fintype.card Unit ≤ 5 * w3.cost id :=
  PfLean.block3_bound (e := (id : Unit → Unit)) Function.injective_id w3 w3_unit

-- and the bound would be violated by a hypothetical cost-2 one-role unit:
example : ¬ (14 * Fintype.card Unit ≤ 5 * 2) := by decide

#print axioms w3_unit
#print axioms w3_cost
