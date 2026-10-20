import Work.FourierC588M.Seam
import Work.Block.WHT

/-!
# Walsh-Hadamard transform from the C588M kernel program

New file (13 Oct 2026); not part of github.com/openai/math and not part of Ryan Shea's tree
(github.com/shea256/fourier-transform-below-nlogn, lean/complex588 at 5239484), which it is laid over
without changing any of his files.

The kernel program is `SeamC588M.engine_C588Mx` (Work/FourierC588M/Seam.lean of the C588M package): the
certificate `GXD.C588M` of Jeremy Kothe's #588 circuit (CrocSwap/integer-mult-bounds hub) as Ryan Shea
converted it to Lean (his conversion pipeline, generators and `complex588` instantiation), with Shea's 112
source copies removed: the sources are written in place and restored, which the mutable-sources engine
`Work/XMut/*` (Jacob Sussman's Lean proof, after Zhihao Chen's (jacklightChen) written scalar-interface
argument on hub #425) allows.  Bank geometry: Dugongue (#401) and jacklightChen (#377); gauged #401 Lean
extension `GX.engine_of_ginv_401` and Fourier chain templates: Nil Mamano.  Everything is checked in Lean
through the engine and checkers of Jacob Sussman's wht-power-saving-lean, at the exponent
`1 - 13552607/10^10` (annealed bank table).

The Walsh-Hadamard corollary is `BlockWHT.wht_main_of_engine` (Work/Block/WHT.lean, Jacob Sussman's engine
file, byte-identical in Shea's tree and in wht-power-saving-lean 9c94857): pre-scale, kernel engine,
post-scale.

This file adds one term and nothing else.  No `sorry`, no `native_decide`, no axiom.
-/

namespace OAI
namespace PowerSaving
namespace WHTC588M
open RAM

/-- **Walsh-Hadamard transform at `1 - 13552607/10^10`**: `BlockWHT.wht_main_of_engine` applied to the
kernel program `SeamC588M.engine_C588Mx` (colour `Paint.left`, transform mode `false`). -/
theorem wht_main_C588Mx :
    ∃ solve W, WHT.WHTProgram solve W ∧ BlockWHT.WHTTimeBoundsAt (1 - 13552607/(10:ℝ)^10) W :=
  BlockWHT.wht_main_of_engine (1 - 13552607/(10:ℝ)^10) (by norm_num) (by norm_num)
    (fun k v => SeamC588M.engine_C588Mx Paint.left false k v)

end WHTC588M
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.WHTC588M.wht_main_C588Mx
