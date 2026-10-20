import Work.Mixed.Bridge.Inst3945F

/-!
# fc-engine: negative controls of the design-B full-cut instance (#3945F); never imported by the chain

* `hcutF3945F_old`: the full-cut profile identity is FALSE with #3945's own head count (57,312 heads, the published fill):
  the profile `cut3945F` is design B's (80,640 heads), not #3945's;
* `hcutF3945F_more`: it is false with one head more (80,641), for all ranks at once (the per-rank form is
  `hcutF3945F_sharp` in `Inst3945F`);
* `hcutF3945F_unsplit`: the profile with the 57,312 non-head idle blocks left at rank 48 (no split) is not `cut3945F`:
  without the split a rank-48 block would remain, above the cut target 25 (separation needs every paid rank below 25).
fc-engine (wht45/fullcut, 20 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI
namespace PowerSaving.PB
open MX

theorem hcutF3945F_old : ¬ (∀ r, BlockAccounting.histFn cut3945F r
    = (if r = 24 + 24 then 0 else BlockAccounting.histFn unit3945H r)
      + (if r = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 57312 else 0)
      + (if r = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 57312 else 0)) := by
  intro h
  have h24 := h 24
  revert h24
  decide

theorem hcutF3945F_more : ¬ (∀ r, BlockAccounting.histFn cut3945F r
    = (if r = 24 + 24 then 0 else BlockAccounting.histFn unit3945H r)
      + (if r = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 80641 else 0)
      + (if r = 24 then BlockAccounting.histFn unit3945H (24 + 24) - 80641 else 0)) := by
  intro h
  have h24 := h 24
  revert h24
  decide

theorem hcutF3945F_unsplit : BlockAccounting.histFn cut3945F 48 ≠ BlockAccounting.histFn unit3945H 48 - 80640 := by
  decide

end PowerSaving.PB
end OAI

#print axioms OAI.PowerSaving.PB.hcutF3945F_old
#print axioms OAI.PowerSaving.PB.hcutF3945F_more
#print axioms OAI.PowerSaving.PB.hcutF3945F_unsplit
