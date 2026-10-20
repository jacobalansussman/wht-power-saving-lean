import Work.KLow.NoY
import Work.KLow.GlueX
import Work.GCert.Data.Gen.U8XZL
set_option maxRecDepth 100000
namespace GXD.U8X
open GXD
theorem zb0f : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 0 false) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb0t : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 0 true) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb1f : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 1 false) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb1t : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 1 true) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb2f : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 2 false) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb2t : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 2 true) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb3f : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([112] ++ GXD.bitB lel 3 false) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
theorem zb3t : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ ([] ++ GXD.bitB lel 3 true) (GXD.slotsL (GXD.nyRaw raw)) = true := by decide +kernel
end GXD.U8X
#print axioms GXD.U8X.zb0f
#print axioms GXD.U8X.zb0t
#print axioms GXD.U8X.zb1f
#print axioms GXD.U8X.zb1t
#print axioms GXD.U8X.zb2f
#print axioms GXD.U8X.zb2t
#print axioms GXD.U8X.zb3f
#print axioms GXD.U8X.zb3t
