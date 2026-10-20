import Work.XHost.HBank
import Work.GCert.Data.C588MGlue401
import Work.GCert.Data.InstC588M

/-!
# (key: bridge-host) the control instance `C588MH`: m-gen's copy-free #588 certificate `GXD.C588M` through the HOSTED
network with the empty host table

Same gauged invocation package (`GXD.C588M.ginv`: Kothe's #588 circuit, x-mutable's mutable sources, Mamano's gauges),
same bank table (`GXD.C588M.sched`, 4,032 banks), as a host schedule with no host (`HostSched.ofBank`).  The hosted unit
price equals Shea's `unitPrice401`, so the rate certificate `BlockAccounting.fold_C588Mx` (s = 40) applies unchanged.
Written by h-L3-build (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

namespace GXD.C588MH
open GXD OAI.PowerSaving OAI.PowerSaving.CB OAI.PowerSaving.XHost
noncomputable section

/-- the bank table of `GXD.C588M` as a host schedule with no host -/
def hsched : HostSched (Fin GXD.C588M.raw.R) (Fin GXD.C588M.raw.v) :=
  HostSched.ofBank (Fin GXD.C588M.raw.v) GXD.C588M.sched

theorem hfits : hsched.fits GXD.C588M.ginv.ent :=
  HostSched.ofBank_fits GXD.C588M.sched GXD.C588M.ginv.ent GXD.C588M.fitsS

theorem hwidth : 64 * Fintype.card (Fin GXD.C588M.raw.v) + hsched.nb = 11200 := GXD.C588M.width

/-- the hosted unit price of the control is the banked unit price of `GXD.C588M` -/
theorem hprice (φ : ℕ → ℝ) : unitPriceHost GXD.C588M.ginv hsched φ = BridgeRate.C588M.unitPrice φ := by
  rw [← GXD.C588M.price φ]
  unfold unitPriceHost hsched
  rw [hcost_ofBank]
  rfl

end
end GXD.C588MH

#print axioms GXD.C588MH.hfits
#print axioms GXD.C588MH.hprice
