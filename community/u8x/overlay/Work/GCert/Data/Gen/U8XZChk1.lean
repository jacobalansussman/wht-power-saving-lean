import Work.KLow.NoY
import Work.KLow.GlueX
import Work.GCert.Data.Gen.U8XZL
set_option maxRecDepth 100000
namespace GXD.U8X
open GXD
theorem zb4f : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112, 113] ++ GXD.bitB lel 4 false) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb4t : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 4 true) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb5f : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 5 false) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb5t : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 5 true) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb6f : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 6 false) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb6t : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 6 true) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zh0 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ [587]) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zc0 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk0 el = true := by decide +kernel
end GXD.U8X
#print axioms GXD.U8X.zb4f
#print axioms GXD.U8X.zb4t
#print axioms GXD.U8X.zb5f
#print axioms GXD.U8X.zb5t
#print axioms GXD.U8X.zb6f
#print axioms GXD.U8X.zb6t
#print axioms GXD.U8X.zh0
#print axioms GXD.U8X.zc0
