import Work.PfLean.Dem3Defs

/-! The eight block checks, evaluated by the kernel (`decide +kernel`; no `native_decide`). -/

namespace OAI.PowerSaving.PfLean.Dem3Comp

set_option profiler true

theorem blk0 : BlockCheck 0 := by decide +kernel
theorem blk1 : BlockCheck 1 := by decide +kernel
theorem blk2 : BlockCheck 2 := by decide +kernel
theorem blk3 : BlockCheck 3 := by decide +kernel
theorem blk4 : BlockCheck 4 := by decide +kernel
theorem blk5 : BlockCheck 5 := by decide +kernel
theorem blk6 : BlockCheck 6 := by decide +kernel
theorem blk7 : BlockCheck 7 := by decide +kernel

theorem blkOK (s : ℕ) (hs : s < 8) : BlockCheck s := by
  interval_cases s
  exacts [blk0, blk1, blk2, blk3, blk4, blk5, blk6, blk7]

end OAI.PowerSaving.PfLean.Dem3Comp

#print axioms OAI.PowerSaving.PfLean.Dem3Comp.blkOK
