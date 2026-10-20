import Work.KMX.GlueXZ
import Work.GCert.Data.Gen.U8XMScal0
import Work.GCert.Data.Gen.U8XMEnd
set_option maxRecDepth 100000
/-! # (key: k-combine) Lean negative controls of the new kernel checks of `GXD.U8XM` (all must evaluate to `false`;
the exact Python mirrors predicted `false` for each, k-combine/logs/controls-U8XM.json)
* the mutable-source zero check sees the content of source x_0 in another row (x_0 is added into x_7 and into
  helpers: sources are written, so an x row can read another register);
* Sussman's CLASSIC scalar check refuses the certificate (a phase-A gate writes a source: kind test `7 += 237`);
* k-lower's exit-aware end check with an EMPTY exit list rejects it (the exit helpers do not end at full frames). -/
namespace GXD.U8XM
open GXD
theorem zcheck_control : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ [0] (GXD.rowsXS (GXD.nyRaw raw)) = false := by
  decide +kernel
theorem classic_rejected : GS.scalCheck raw ⟨10, 22, 8, 16, 8⟩ 0 28 = false := by
  decide +kernel
theorem endkx_control : GLab.endKX par 112 ports (raw.ret.map Prod.fst) gl [] st0 st3 st6 = false := by
  decide +kernel
end GXD.U8XM
#print axioms GXD.U8XM.zcheck_control
#print axioms GXD.U8XM.classic_rejected
#print axioms GXD.U8XM.endkx_control
