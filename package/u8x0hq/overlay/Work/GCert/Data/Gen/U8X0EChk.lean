import CSX.ECheckZX
import Work.GCert.Data.Gen.U8X0Gl
set_option maxRecDepth 100000
namespace GXD.U8X0
open GXD
theorem echk0_ok : GS.eCheckZX raw ⟨10, 8, 8, 16, 8⟩ glb0 = true := by decide +kernel
theorem echk1_ok : GS.eCheckZX raw ⟨10, 8, 8, 16, 8⟩ glb1 = true := by decide +kernel
theorem echk2_ok : GS.eCheckZX raw ⟨10, 8, 8, 16, 8⟩ glb2 = true := by decide +kernel
end GXD.U8X0
#print axioms GXD.U8X0.echk0_ok
#print axioms GXD.U8X0.echk1_ok
#print axioms GXD.U8X0.echk2_ok
