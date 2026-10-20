import Work.GCert.Data.Gen.C588MScal0
import Work.GCert.Data.Gen.C588MSched
set_option maxRecDepth 100000
namespace GXD.C588MControls
open GXD
/-- the sign of the first coefficient of the first scatter row flipped -/
def flipRows : List (List (Nat × Co)) → List (List (Nat × Co))
  | ((k, c) :: r) :: rows => ((k, {c with neg := !c.neg}) :: r) :: rows
  | rows => rows
def flipScat : Scat → Scat
  | .table rows => .table (flipRows rows)
  | s => s
noncomputable def badRaw : Raw := {C588M.raw with scat := flipScat C588M.raw.scat}
/-- the mutable-source scalar check rejects a corrupted scatter coefficient (y_0 no longer receives exactly x_0) -/
theorem bad_scalar_rejected : GS.scalCheckZ badRaw ⟨10, 21, 4, 8, 8⟩ 0 28 = false := by decide +kernel
/-- Sussman's classic scalar check refuses the copy-free certificate itself: it writes sources -/
theorem classic_rejected : GS.scalCheck C588M.raw ⟨10, 21, 4, 8, 8⟩ 0 28 = false := by decide +kernel
noncomputable def badRows : List (List Nat) := match C588M.schedRows with
  | (x::_y::xs)::rs => (x::x::xs)::rs
  | other => other
theorem duplicate_owner_rejected : OAI.PowerSaving.GX.schedK 4005 C588M.raw.R badRows C588M.schedPad C588M.schedRes = false := by decide +kernel
end GXD.C588MControls
#print axioms GXD.C588MControls.bad_scalar_rejected
#print axioms GXD.C588MControls.classic_rejected
#print axioms GXD.C588MControls.duplicate_owner_rejected
