import Work.Gauge.ECheck
import Work.GCert.Data.Gen.U8XGl
set_option maxRecDepth 100000
namespace GXD.U8X
open GXD
theorem echk0_ok : GS.eCheck raw ⟨10, 8, 1, 16, 8⟩ glb0 = true := by decide +kernel
theorem echk1_ok : GS.eCheck raw ⟨10, 8, 1, 16, 8⟩ glb1 = true := by decide +kernel
theorem echk2_ok : GS.eCheck raw ⟨10, 8, 1, 16, 8⟩ glb2 = true := by decide +kernel
end GXD.U8X
#print axioms GXD.U8X.echk0_ok
#print axioms GXD.U8X.echk1_ok
#print axioms GXD.U8X.echk2_ok
