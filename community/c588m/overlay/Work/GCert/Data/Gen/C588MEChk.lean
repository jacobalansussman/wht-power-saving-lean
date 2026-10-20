import Work.XMut.ECheckZ
import Work.GCert.Data.Gen.C588MGl
set_option maxRecDepth 100000
namespace GXD.C588M
open GXD
theorem echk0_ok : GS.eCheckZ raw ⟨10, 8, 4, 8, 8⟩ glb0 = true := by decide +kernel
theorem echk1_ok : GS.eCheckZ raw ⟨10, 8, 4, 8, 8⟩ glb1 = true := by decide +kernel
theorem echk2_ok : GS.eCheckZ raw ⟨10, 8, 4, 8, 8⟩ glb2 = true := by decide +kernel
end GXD.C588M
#print axioms GXD.C588M.echk0_ok
#print axioms GXD.C588M.echk1_ok
#print axioms GXD.C588M.echk2_ok
