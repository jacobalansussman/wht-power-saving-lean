import Work.FourierC588.Seam
import Work.Block.WHT

/-!
# Walsh-Hadamard transform from Ryan Shea's C588 kernel program

New file (13 Oct 2026, agent c-wht, for Jacob Sussman); not part of github.com/openai/math and not part of
Ryan Shea's tree (github.com/shea256/fourier-transform-below-nlogn, lean/complex588 at 5239484), which it
is laid over without changing any of his files.

The kernel program is Shea's `SeamC588.engine_C588x` (Work/FourierC588/Seam.lean): the certificate
`GXD.C588` of Kothe's hub #588 circuit (on the community line Dugongue #401, jacklightChen #377, Mamano's
gauged #401 Lean extension), checked in Lean through the engine and checkers of Jacob Sussman's
wht-power-saving-lean, at the exponent `1 - 12331686/10^10`.  The Walsh-Hadamard corollary is
`BlockWHT.wht_main_of_engine` (Work/Block/WHT.lean, Jacob Sussman's engine file, byte-identical in Shea's
tree and in wht-power-saving-lean 9c94857): pre-scale, kernel engine, post-scale.

This file adds one term and nothing else.  No `sorry`, no `native_decide`, no axiom.
-/

namespace OAI
namespace PowerSaving
namespace WHTC588
open RAM

/-- **Walsh-Hadamard transform at `1 - 12331686/10^10`**: `BlockWHT.wht_main_of_engine` applied to Shea's
kernel program `SeamC588.engine_C588x` (colour `Paint.left`, transform mode `false`). -/
theorem wht_main_C588x :
    ∃ solve W, WHT.WHTProgram solve W ∧ BlockWHT.WHTTimeBoundsAt (1 - 12331686/(10:ℝ)^10) W :=
  BlockWHT.wht_main_of_engine (1 - 12331686/(10:ℝ)^10) (by norm_num) (by norm_num)
    (fun k v => SeamC588.engine_C588x Paint.left false k v)

end WHTC588
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.WHTC588.wht_main_C588x
