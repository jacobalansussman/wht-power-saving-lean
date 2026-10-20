import Work.WHTC588Q.Main

/-!
Comparator solution `Work.WHTC588Q.Solution` (not part of github.com/openai/math).

The whole proof is in `Work.WHTC588Q.Main` (and its imports), which ends with

    theorem WHTC588Q.wht_main_C588Qx :
        ∃ solve W, WHTProgram solve W ∧ BlockWHT.WHTTimeBoundsAt (1 - 12383292/(10:ℝ)^10) W

This module adds nothing mathematical: it defines `WHT.WHTTimeBoundsAt` with the text of the
challenge (the same definition as `BlockWHT.WHTTimeBoundsAt` of `Work.Block.WHT`) and restates
the theorem under the name the challenge uses.  Its proof is that single term.
`WHTProgram`, `wht`, `whtInKind`, `whtOutKind` come from `WHTCheck.Solution`.
Same text as `Work/GCert/Data/SolutionB2Ge8x.lean` of wht-power-saving-lean (9c94857), with the
exponent, the theorem name and the import changed.
-/

namespace OAI
namespace PowerSaving
namespace WHT
open Filter Asymptotics

/-- `W = O(2^k (k+1)^β)` and `W = o(2^k k)` (`WHTTimeBounds` is the case `β = 1 - 2/10^11`). -/
def WHTTimeBoundsAt (β : ℝ) (W : ℕ → ℕ) : Prop :=
  let W' := fun k : ℕ => (W k : ℝ)
  W' =O[atTop] (fun k : ℕ => (2:ℝ)^k * ((k:ℝ)+1)^β) ∧
  W' =o[atTop] (fun k : ℕ => (2:ℝ)^k * (k:ℝ))

/-- **Walsh-Hadamard transform with the improved exponent**: the statement `WHTGoal` of
`WHTCheck/Challenge.lean` with `1 - 2/10^11` replaced by `1 - 12383292/10^10`. -/
theorem wht_main_block_C588Qx :
    ∃ solve W, WHTProgram solve W ∧ WHTTimeBoundsAt (1 - 12383292/(10:ℝ)^10) W :=
  WHTC588Q.wht_main_C588Qx

end WHT
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.WHT.wht_main_block_C588Qx
