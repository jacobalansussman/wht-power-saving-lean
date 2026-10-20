import Work.XHost.HBank
import Work.XHost.HPriceK
import Work.XHost.Data.HU8XMSched
import Work.XHost.Data.InstU8XMHH
import Work.GCert.Data.U8XMGlue401

/-!
# (key: bridge-host) the instance `U8XMHH`: k-combine's U8XM certificate `GXD.U8XM` through the HOSTED network

Same gauged invocation package (`GXD.U8XM.ginv = ginvX.toE'`: the U8 word of Kothe #869/#870, Boukhalfa #890, Capalbo #942,
hub exit frames, x-mutable's mutable sources with 5 centre-source copies, Mamano's gauges);
host schedule: h-L3-design's same-run host table (`host_price.py --samerun 1`, after x-hosting / DaysSky's hosting_core
layout) with the 310 banks of s-stack's re-pack, checked by `XHost.hostSchedK`.  Unit: W = 7478, kernel numerator 16197706.
Written by fs-price/tools/gen_hu8xm.py (after h-L3-build's gen_h588m.py; wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option maxRecDepth 100000
namespace GXD.U8XMHH
open GXD OAI.PowerSaving OAI.PowerSaving.CB OAI.PowerSaving.XHost
noncomputable section

/-- the host schedule (checked table) -/
def hsched : HostSched (Fin GXD.U8XM.raw.R) (Fin GXD.U8XM.raw.v) :=
  HostSched.ofTab 509 310 hrowsL browsL hpadL hremL GXD.U8XM.schedRes hs_ok

theorem hfits : hsched.fits GXD.U8XM.ginv.ent := by
  unfold HostSched.fits
  have e : (fun q : Fin GXD.U8XM.raw.R => Fintype.card (Fin GXD.U8XM.aux.p.h) - (GXD.U8XM.ginv.ent q).dim)
      = (fun q => GXD.U8XM.schedRes.getD q.val 0) := funext GXD.U8XM.hent
  rw [e]
  exact HostSched.ofTab_fitsR 509 310 hrowsL browsL hpadL hremL GXD.U8XM.schedRes hs_ok

theorem hwidth : 64 * Fintype.card (Fin GXD.U8XM.raw.v) + hsched.nb = 7478 := by
  rw [Fintype.card_fin]
  rfl

theorem hrem_price (φ : ℕ → ℝ) : hcostS hsched 0 φ + hcostS hsched 1 φ + hcostS hsched 2 φ + hcostS hsched 3 φ
    + hcostS hsched 4 φ + hcostS hsched 5 φ + hcostS hsched 6 φ + hcostS hsched 7 φ + hcostS hsched 8 φ
      = rcost φ hremL := by
  rw [← sum_fin9 (fun c => hcostS hsched c φ)]
  unfold hcostS
  refine (remTot (fun i => bcost φ (hremL.getD i 0))).trans ?_
  rw [← rem_len]
  exact rcost_range φ hremL

theorem hpad_price (φ : ℕ → ℝ) : ∑ b, bcost φ (hsched.pad b) = rcost φ hpadL :=
  rcost_fin φ hpadL 310 pad_len

/-- the hosted unit price of `U8XMHH` -/
theorem hprice (φ : ℕ → ℝ) : unitPriceHost GXD.U8XM.ginv hsched φ = BridgeRate.U8XMHH.unitPrice φ := by
  unfold unitPriceHost
  rw [GXD.U8XM.ginv_cost, hrem_price, hpad_price, GXD.rcost_wsumK φ hremL 48 remBound,
    GXD.rcost_wsumK φ hpadL 61 padBound]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, rankR0_ok, rankR1_ok, rankR2_ok, rankR3_ok, rankR4_ok, rankR5_ok, rankR6_ok, rankR7_ok, rankR8_ok, rankR9_ok, rankR10_ok, rankR11_ok, rankR12_ok, rankR13_ok, rankR14_ok, rankR15_ok, rankR16_ok, rankR17_ok, rankR18_ok, rankR19_ok, rankR20_ok, rankR21_ok, rankR22_ok, rankR23_ok, rankR24_ok, rankR25_ok, rankR26_ok, rankR27_ok, rankR28_ok, rankR29_ok, rankR30_ok, rankR31_ok, rankR32_ok, rankR33_ok, rankR34_ok, rankR35_ok, rankR36_ok, rankR37_ok, rankR38_ok, rankR39_ok, rankR40_ok, rankR41_ok, rankR42_ok, rankR43_ok, rankR44_ok, rankR45_ok, rankR46_ok, rankR47_ok, rankR48_ok, rankP0_ok, rankP1_ok, rankP2_ok, rankP3_ok, rankP4_ok, rankP5_ok, rankP6_ok, rankP7_ok, rankP8_ok, rankP9_ok, rankP10_ok, rankP11_ok, rankP12_ok, rankP13_ok, rankP14_ok, rankP15_ok, rankP16_ok, rankP17_ok, rankP18_ok, rankP19_ok, rankP20_ok, rankP21_ok, rankP22_ok, rankP23_ok, rankP24_ok, rankP25_ok, rankP26_ok, rankP27_ok, rankP28_ok, rankP29_ok, rankP30_ok, rankP31_ok, rankP32_ok, rankP33_ok, rankP34_ok, rankP35_ok, rankP36_ok, rankP37_ok, rankP38_ok, rankP39_ok, rankP40_ok, rankP41_ok, rankP42_ok, rankP43_ok, rankP44_ok, rankP45_ok, rankP46_ok, rankP47_ok, rankP48_ok, rankP49_ok, rankP50_ok, rankP51_ok, rankP52_ok, rankP53_ok, rankP54_ok, rankP55_ok, rankP56_ok, rankP57_ok, rankP58_ok, rankP59_ok, rankP60_ok, rankP61_ok]
  unfold BridgeRate.U8XMHH.unitPrice BridgeRate.U8XMHH.remPrice BridgeRate.U8XMHH.padPrice BridgeRate.U8XM.invPrice
    GXD.U8XM.invPrice
  norm_num [bcost]
  try ring

end
end GXD.U8XMHH

#print axioms GXD.U8XMHH.hfits
#print axioms GXD.U8XMHH.hprice
