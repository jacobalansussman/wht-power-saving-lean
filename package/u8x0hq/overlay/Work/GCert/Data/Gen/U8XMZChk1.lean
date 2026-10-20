import Work.KLow.NoY
import Work.KMX.GlueXZ
import Work.GCert.Data.Gen.U8XMZL
set_option maxRecDepth 100000
namespace GXD.U8XM
open GXD
theorem zb4f : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112, 113] ++ GXD.bitB lel 4 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb4t : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 4 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb5f : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 5 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb5t : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 5 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb6f : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 6 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb6t : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 6 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zh0 : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ [587]) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zc0 : GS.zCheckZ (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk0 el = true := by decide +kernel
end GXD.U8XM
#print axioms GXD.U8XM.zb4f
#print axioms GXD.U8XM.zb4t
#print axioms GXD.U8XM.zb5f
#print axioms GXD.U8XM.zb5t
#print axioms GXD.U8XM.zb6f
#print axioms GXD.U8XM.zb6t
#print axioms GXD.U8XM.zh0
#print axioms GXD.U8XM.zc0
