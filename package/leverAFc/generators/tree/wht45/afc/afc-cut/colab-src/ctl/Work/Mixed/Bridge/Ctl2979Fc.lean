import Work.Mixed.Bridge.Inst2979Fc
import Work.Mixed.Bridge.Inst2979H

/-!
# afc-cut: negative controls of the A-fc full-cut instance (#2979Fc); never imported by the chain

* `hcutF2979Fc_fewerAll`: the full-cut profile identity is FALSE with one head fewer (16,575), for all ranks at once
  (a rank-48 block would be left and split; the per-rank form is `hcutF2979Fc_fewer` in `Inst2979Fc`);
* `hcutF2979Fc_designA`: it is false for design A's unit profile `PB.unit2979H` (the every-head-kept fill 764395fc of
  the same words, 28,190 rank-48 blocks = 16,576 heads + 11,614 others): `cut2979Fc` is A-fc's profile, not design A's;
* `hcutF2979Fc_noSplitA`: design A's profile with the heads removed still has 11,614 rank-48 blocks, above the cut
  target 25, while `cut2979Fc` has none (on A-fc nothing is left to split; on design A there would be);
* `unit2979Fc_ne_designA`: A-fc's unit profile differs from design A's (rank 23: 25,905 vs 0).
afc-cut (wht45/afc, 20 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

namespace OAI
namespace PowerSaving.PB

theorem hcutF2979Fc_fewerAll : ¬ (∀ r, BlockAccounting.histFn MX.cut2979Fc r
    = (if r = 24 + 24 then 0 else BlockAccounting.histFn MX.unit2979Fc r)
      + (if r = 24 then BlockAccounting.histFn MX.unit2979Fc (24 + 24) - 16575 else 0)
      + (if r = 24 then BlockAccounting.histFn MX.unit2979Fc (24 + 24) - 16575 else 0)) := by
  intro h
  have h24 := h 24
  revert h24
  decide

theorem hcutF2979Fc_designA : ¬ (∀ r, BlockAccounting.histFn MX.cut2979Fc r
    = (if r = 24 + 24 then 0 else BlockAccounting.histFn unit2979H r)
      + (if r = 24 then BlockAccounting.histFn unit2979H (24 + 24) - 16576 else 0)
      + (if r = 24 then BlockAccounting.histFn unit2979H (24 + 24) - 16576 else 0)) := by
  intro h
  have h1 := h 1
  revert h1
  decide

theorem hcutF2979Fc_noSplitA : BlockAccounting.histFn unit2979H 48 - 16576 = 11614 ∧
    BlockAccounting.histFn MX.cut2979Fc 48 = 0 := by
  decide

theorem unit2979Fc_ne_designA : BlockAccounting.histFn MX.unit2979Fc 23 ≠ BlockAccounting.histFn unit2979H 23 := by
  decide

end PowerSaving.PB
end OAI

#print axioms OAI.PowerSaving.PB.hcutF2979Fc_fewerAll
#print axioms OAI.PowerSaving.PB.hcutF2979Fc_designA
#print axioms OAI.PowerSaving.PB.hcutF2979Fc_noSplitA
#print axioms OAI.PowerSaving.PB.unit2979Fc_ne_designA
