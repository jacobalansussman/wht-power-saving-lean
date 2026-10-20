import Work.KLow.NoY
import CSX.GlueS
import CSX.TwinX
import Work.GCert.Data.Gen.U8X0ZL
set_option maxRecDepth 100000
namespace GXD.U8X0
open GXD
theorem zb4f : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112, 113] ++ GXD.bitB lel 4 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb4t : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 4 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb5f : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 5 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb5t : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 5 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb6f : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 6 false) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb6t : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 6 true) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zh0 : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ [587]) (GXD.rowsXS (GXD.nyRaw raw)) = true := by decide +kernel
theorem zc0 : GS.zCheckZX (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk0 el = true := by decide +kernel
end GXD.U8X0
#print axioms GXD.U8X0.zb4f
#print axioms GXD.U8X0.zb4t
#print axioms GXD.U8X0.zb5f
#print axioms GXD.U8X0.zb5t
#print axioms GXD.U8X0.zb6f
#print axioms GXD.U8X0.zb6t
#print axioms GXD.U8X0.zh0
#print axioms GXD.U8X0.zc0
