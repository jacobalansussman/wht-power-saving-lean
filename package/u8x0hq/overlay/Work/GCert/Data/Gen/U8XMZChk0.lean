import Work.KLow.NoY
import Work.KMX.GlueXZ
import Work.GCert.Data.Gen.U8XMZL
set_option maxRecDepth 100000
namespace GXD.U8XM
open GXD
theorem zb0f : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 0 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb0t : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 0 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb1f : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 1 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb1t : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 1 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb2f : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 2 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb2t : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 2 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb3f : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 3 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb3t : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 3 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
end GXD.U8XM
#print axioms GXD.U8XM.zb0f
#print axioms GXD.U8XM.zb0t
#print axioms GXD.U8XM.zb1f
#print axioms GXD.U8XM.zb1t
#print axioms GXD.U8XM.zb2f
#print axioms GXD.U8XM.zb2t
#print axioms GXD.U8XM.zb3f
#print axioms GXD.U8XM.zb3t
