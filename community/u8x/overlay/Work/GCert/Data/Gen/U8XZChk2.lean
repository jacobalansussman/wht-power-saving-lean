import Work.KLow.NoY
import Work.KLow.GlueX
import Work.GCert.Data.Gen.U8XZL
set_option maxRecDepth 100000
namespace GXD.U8X
open GXD
theorem zc1 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk1 el = true := by decide +kernel
theorem zc2 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk2 el = true := by decide +kernel
theorem zc3 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk3 el = true := by decide +kernel
theorem zc4 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk4 el = true := by decide +kernel
theorem zc5 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk5 el = true := by decide +kernel
theorem zc6 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk6 el = true := by decide +kernel
theorem zc7 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk7 el = true := by decide +kernel
theorem zc8 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk8 el = true := by decide +kernel
end GXD.U8X
#print axioms GXD.U8X.zc1
#print axioms GXD.U8X.zc2
#print axioms GXD.U8X.zc3
#print axioms GXD.U8X.zc4
#print axioms GXD.U8X.zc5
#print axioms GXD.U8X.zc6
#print axioms GXD.U8X.zc7
#print axioms GXD.U8X.zc8
