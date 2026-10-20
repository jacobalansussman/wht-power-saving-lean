import Work.KLow.NoY
import Work.KLow.GlueX
import Work.GCert.Data.Gen.U8XZL
set_option maxRecDepth 100000
namespace GXD.U8X
open GXD
theorem zc25 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk25 el = true := by decide +kernel
theorem zc26 : GS.zCheck (GXD.nyRaw raw) ⟨10, 36, 32, 32, 1⟩ chunk26 el = true := by decide +kernel
end GXD.U8X
#print axioms GXD.U8X.zc25
#print axioms GXD.U8X.zc26
