import Work.FourierC588MHH.Seam
import Work.Block.WHT

/-!
# Walsh-Hadamard transform from the C588MHH kernel program

New file (13 Oct 2026); not part of github.com/openai/math and not part of Ryan Shea's tree
(github.com/shea256/fourier-transform-below-nlogn, lean/complex588 at 5239484), which it is laid over
without changing any of his files.

The kernel program is `SeamC588MHH.engine_C588MHHx` (Work/FourierC588MHH/Seam.lean of the C588MHH package): the
certificate `GXD.C588M` of Jeremy Kothe's #588 circuit (CrocSwap/integer-mult-bounds hub) as Ryan Shea
converted it to Lean (his conversion pipeline, generators and `complex588` instantiation), with Shea's 112
source copies removed: the sources are written in place and restored, which the mutable-sources engine
`Work/XMut/*` (Jacob Sussman's Lean proof, after Zhihao Chen's (jacklightChen) written scalar-interface
argument on hub #425) allows; and with same-run hosting (DaysSky, hub #442): bank entries are held in data
slots that are idle in the same run, through the hosted network `Work/XHost/*` (Jacob Sussman's Lean proof,
ending in `XHost.engine_of_ginv_host`; host and bank tables checked by `XHost.hostSchedK`, whose soundness
is `HostSched.ofTab_fitsR`), with 111 banks and unit width 7279.  The hosted network restates lemmas that
were already in Shea's tree (`B401.gauged_step`, `B401.gauged_Wn`, `GX.entW_spec`, `GX.entW_form`,
`GX.gauge_step`), built on Dugongue's #401 entrance gauges, Nil Mamano's gauged #401 extension and Jacob
Sussman's `core_nested`.  Bank geometry: Dugongue (#401) and jacklightChen (#377); Fourier chain templates:
Nil Mamano.  Everything is checked in Lean through the engine and checkers of Jacob Sussman's
wht-power-saving-lean, at the exponent `1 - 15699306/10^10`.

The Walsh-Hadamard corollary is `BlockWHT.wht_main_of_engine` (Work/Block/WHT.lean, Jacob Sussman's engine
file, byte-identical in Shea's tree and in wht-power-saving-lean 9c94857): pre-scale, kernel engine,
post-scale.

This file adds one term and nothing else.  No `sorry`, no `native_decide`, no axiom.
-/

namespace OAI
namespace PowerSaving
namespace WHTC588MHH
open RAM

/-- **Walsh-Hadamard transform at `1 - 15699306/10^10`**: `BlockWHT.wht_main_of_engine` applied to the
kernel program `SeamC588MHH.engine_C588MHHx` (colour `Paint.left`, transform mode `false`). -/
theorem wht_main_C588MHHx :
    ∃ solve W, WHT.WHTProgram solve W ∧ BlockWHT.WHTTimeBoundsAt (1 - 15699306/(10:ℝ)^10) W :=
  BlockWHT.wht_main_of_engine (1 - 15699306/(10:ℝ)^10) (by norm_num) (by norm_num)
    (fun k v => SeamC588MHH.engine_C588MHHx Paint.left false k v)

end WHTC588MHH
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.WHTC588MHH.wht_main_C588MHHx
