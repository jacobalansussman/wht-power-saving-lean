import Work.XHostT.Data.S5.Calls
import Work.XHostT.Data.D3945.Disp

/-!
# The provider field of p-hostck's host table is p-disp's dispatch table

p-hostck's call tree `DS5.ct` (from p-own's #3945 host table, `Work/XHostT/Data/S5/Calls.lean`) carries, in the leaf of
call `k = 9 ρ + l`, the provider of the call (`(ct.get k).2.2`; 0 = extra, 1 = gauge44).  p-disp's table
`D3945.disp` (from the hub's DISPATCH.json / ORDER.json) is bit `k` of `D3945.dispKB`.  `dpChk_ok` checks the two agree
on all 3,240 calls (`decide +kernel`, one bit read and one tree lookup per call); `dp_disp` states it for every call
`(ρ, l)`, so the residuals of p-hostck's check (`resAt resT (ct.get k).2.2 q`) are those of the dispatched providers.
Written by p-disp (wht45/lever, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5
open OAI.PowerSaving.XHostT

/-- the provider of every call `k < n` of the tree is bit `k` of the dispatch table. -/
noncomputable def dpChk : Nat → Bool
  | 0 => true
  | n + 1 => dpChk n && Nat.beq (ct.get n).2.2 ((D3945.dispKB >>> n) &&& 1)

theorem dpChk_ok : dpChk 3240 = true := by decide +kernel

theorem dpChk_sound : ∀ n, dpChk n = true → ∀ k < n, (ct.get k).2.2 = (D3945.dispKB >>> k) &&& 1 := by
  intro n
  induction n with
  | zero => intro _ k hk; exact absurd hk (Nat.not_lt_zero k)
  | succ n ih =>
    intro h k hk
    have h' : dpChk n = true ∧ Nat.beq (ct.get n).2.2 ((D3945.dispKB >>> n) &&& 1) = true := by
      simpa [dpChk] using h
    rcases Nat.lt_succ_iff_lt_or_eq.mp hk with hk | rfl
    · exact ih h'.1 k hk
    · exact Nat.eq_of_beq_eq_true h'.2

/-- **the provider of call `(ρ, l)` in p-hostck's table is `D3945.disp (ρ, l)`**. -/
theorem dp_disp (ρ : Fin 360) (l : Fin 9) : (ct.get (9 * ρ.val + l.val)).2.2 = (D3945.disp (ρ, l)).val := by
  rw [dpChk_sound 3240 dpChk_ok _ (by omega)]
  show (D3945.dispKB >>> (9 * ρ.val + l.val)) &&& 1 = (D3945.bitP D3945.dispKB (9 * ρ.val + l.val)).val
  unfold D3945.bitP
  generalize D3945.dispKB >>> (9 * ρ.val + l.val) = x
  rw [Nat.and_one_is_mod]
  rcases Nat.mod_two_eq_zero_or_one x with h | h <;> simp [h]

end OAI.PowerSaving.XHostT.K.DS5
#print axioms OAI.PowerSaving.XHostT.K.DS5.dpChk_ok
#print axioms OAI.PowerSaving.XHostT.K.DS5.dp_disp
