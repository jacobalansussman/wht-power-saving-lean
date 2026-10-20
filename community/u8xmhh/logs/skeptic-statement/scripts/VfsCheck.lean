import Work.FourierU8XMHH.Main
open Lean Elab Command in
elab "#vfs_mod " n:ident : command => do
  let nm := n.getId
  let r ← findModuleOf? nm
  logInfo m!"{nm} is declared in module {r}"
#print axioms OAI.PowerSaving.transform_mainU8XMHH
#print axioms OAI.PowerSaving.convolution_mainU8XMHH
#check (OAI.PowerSaving.transform_mainU8XMHH : OAI.PowerSaving.DFTGoalU8XMHH)
#check (OAI.PowerSaving.convolution_mainU8XMHH : OAI.PowerSaving.ConvGoalU8XMHH)
#print OAI.PowerSaving.decimalExponentU8XMHH
#print OAI.PowerSaving.DFTGoalU8XMHH
#print OAI.PowerSaving.TimeBoundsU8XMHH
example : OAI.PowerSaving.decimalExponentU8XMHH = 1 - 16197705/(10^(10:ℕ)) := rfl
example : OAI.PowerSaving.decimalExponentU8XMHH = 1 - (0.0016197705 : ℝ) := by
  norm_num [OAI.PowerSaving.decimalExponentU8XMHH]
#vfs_mod OAI.PowerSaving.DFTGoalU8XMHH
#vfs_mod OAI.PowerSaving.transform_mainU8XMHH
#vfs_mod OAI.PowerSaving.DFTProgram
#vfs_mod OAI.PowerSaving.ConvProgram
#vfs_mod OAI.PowerSaving.paperTime
#vfs_mod OAI.PowerSaving.nlogn
#vfs_mod OAI.PowerSaving.XHost.engine_of_ginv_host
#vfs_mod OAI.PowerSaving.RAM.hillsU8XMHH_program
