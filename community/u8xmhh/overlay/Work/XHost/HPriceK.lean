import Work.XHost.HSchedK
import Work.Combine.XG

/-!
# (key: bridge-host) The price of the remaining climbs of a checked host table

`remTot`: the 9 × 112 × 2 remaining climbs of a run are the 2,016 entries of the list `rem` (host index
`encH c (s, t) = 224 c + 112 s + t` is a bijection onto `range 2016`), so their price is `rcost φ rem`.
Written by h-L3-build (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI
namespace PowerSaving
namespace XHost
open Finset CB

/-- the host indices as an equivalence. -/
def encE : Fin 9 × (Bool × Fin 112) ≃ Fin 2016 where
  toFun x := ⟨encH x.1 x.2, encH_lt x.1 x.2⟩
  invFun i := decH i.val
  left_inv x := by
    obtain ⟨c, j⟩ := x
    exact decH_encH c j
  right_inv i := Fin.ext (encH_decH i.isLt)

theorem remTot {M : Type} [AddCommMonoid M] (f : ℕ → M) :
    ∑ c : Fin 9, ∑ t : Fin 112, (f (encH c (false, t)) + f (encH c (true, t))) = ∑ i ∈ range 2016, f i := by
  calc ∑ c : Fin 9, ∑ t : Fin 112, (f (encH c (false, t)) + f (encH c (true, t)))
      = ∑ x : Fin 9 × (Bool × Fin 112), f (encH x.1 x.2) := by
        rw [Fintype.sum_prod_type]
        refine Finset.sum_congr rfl (fun c _ => ?_)
        rw [Fintype.sum_prod_type, Fintype.sum_bool, Finset.sum_add_distrib, add_comm]
    _ = ∑ i : Fin 2016, f i.val := Fintype.sum_equiv encE _ _ (fun x => rfl)
    _ = ∑ i ∈ range 2016, f i := Fin.sum_univ_eq_sum_range f 2016

theorem rcost_range (φ : ℕ → ℝ) (l : List ℕ) : ∑ i ∈ range l.length, bcost φ (l.getD i 0) = rcost φ l := by
  unfold rcost
  exact (map_sum_eqH (bcost φ) 0 l).symm

theorem rcost_fin (φ : ℕ → ℝ) (l : List ℕ) (n : ℕ) (hl : l.length = n) :
    ∑ b : Fin n, bcost φ (l.getD b.val 0) = rcost φ l := by
  rw [← rcost_range, hl]
  exact Fin.sum_univ_eq_sum_range (fun b => bcost φ (l.getD b 0)) n

theorem sum_fin9 {M : Type} [AddCommMonoid M] (g : Fin 9 → M) :
    ∑ c : Fin 9, g c = g 0 + g 1 + g 2 + g 3 + g 4 + g 5 + g 6 + g 7 + g 8 := by
  rw [Fin.sum_univ_succ, Fin.sum_univ_succ, Fin.sum_univ_succ, Fin.sum_univ_succ, Fin.sum_univ_succ,
    Fin.sum_univ_succ, Fin.sum_univ_succ, Fin.sum_univ_succ, Fin.sum_univ_succ, Fin.sum_univ_zero]
  simp only [add_assoc, add_zero]
  rfl

end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.remTot
#print axioms OAI.PowerSaving.XHost.sum_fin9
