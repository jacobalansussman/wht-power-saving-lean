import Work.XMut.ECheckZ
import Work.GCert.Data.Gen.U8XMGl
set_option maxRecDepth 100000
namespace GXD.U8XM
open GXD
theorem echk0_ok : GS.eCheckZ raw ⟨10, 8, 8, 16, 8⟩ glb0 = true := by decide +kernel
theorem echk1_ok : GS.eCheckZ raw ⟨10, 8, 8, 16, 8⟩ glb1 = true := by decide +kernel
theorem echk2_ok : GS.eCheckZ raw ⟨10, 8, 8, 16, 8⟩ glb2 = true := by decide +kernel
end GXD.U8XM
#print axioms GXD.U8XM.echk0_ok
#print axioms GXD.U8XM.echk1_ok
#print axioms GXD.U8XM.echk2_ok
