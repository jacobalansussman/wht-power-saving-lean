import Work.PfLean.Dem3Defs

/-! Negative control (px-lean): one certificate entry changed by 1; the kernel must refute the check. -/

namespace OAI.PowerSaving.PfLean.Dem3Comp

def BfBad (s i k : ℕ) : ℤ := if s = 0 ∧ i = 0 ∧ k = 0 then Bf s i k + 1 else Bf s i k

theorem neg_control : ¬ (∀ i < 28, ∀ j < 28,
    (dotN (BfBad 0 i) (fun k => Gf 0 k j) 28 - if i = j then 1 else 0) % 65537 = 0) := by
  decide +kernel

end OAI.PowerSaving.PfLean.Dem3Comp

#print axioms OAI.PowerSaving.PfLean.Dem3Comp.neg_control
