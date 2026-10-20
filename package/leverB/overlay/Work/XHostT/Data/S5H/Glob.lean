import Work.XHostT.KHist
import Work.XHostT.Data.S5H.GlobC0
import Work.XHostT.Data.S5H.GlobC1
import Work.XHostT.Data.S5H.GlobC2
import Work.XHostT.Data.S5H.GlobC3
import Work.XHostT.Data.S5H.GlobC4
import Work.XHostT.Data.S5H.GlobC5
import Work.XHostT.Data.S5H.GlobC6
import Work.XHostT.Data.S5H.GlobC7
import Work.XHostT.Data.S5H.GlobC8
import Work.XHostT.Data.S5H.GlobC9
import Work.XHostT.Data.S5H.GlobC10
import Work.XHostT.Data.S5H.GlobC11
import Work.XHostT.Data.S5H.GlobC12
import Work.XHostT.Data.S5H.GlobC13
import Work.XHostT.Data.S5H.GlobC14
import Work.XHostT.Data.S5H.GlobC15
import Work.XHostT.Data.S5H.GlobC16
import Work.XHostT.Data.S5H.GlobC17
import Work.XHostT.Data.S5H.GlobC18
import Work.XHostT.Data.S5H.GlobC19
import Work.XHostT.Data.S5H.GlobC20
import Work.XHostT.Data.S5H.GlobC21
import Work.XHostT.Data.S5H.GlobC22
import Work.XHostT.Data.S5H.GlobC23
import Work.XHostT.Data.S5H.GlobC24
import Work.XHostT.Data.S5H.GlobC25
import Work.XHostT.Data.S5H.GlobC26
/-! p-hostck mk_parts.py: global facts of S5H: remainder histogram {0: 560160, 3: 23328, 5: 4320, 48: 137952}, residual sums {0: 3304, 1: 3304} (27 chunks of 120 kinds). -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DS5H
open OAI.PowerSaving.XHostT.K
/-- every remainder is one of [0, 3, 5, 48]. -/
theorem rem_in : inR rt [0, 3, 5, 48] 0 NK = true :=
  inR_list rt [0, 3, 5, 48] 120 27 0 (fun j hj => by
    interval_cases j
    exacts [gin0, gin1, gin2, gin3, gin4, gin5, gin6, gin7, gin8, gin9, gin10, gin11, gin12, gin13, gin14, gin15, gin16, gin17, gin18, gin19, gin20, gin21, gin22, gin23, gin24, gin25, gin26])
theorem rem_c0 : cntR rt 0 0 NK = 560160 := by
  have h := cntR_list rt 0 120 [20384, 20832, 20788, 20380, 20832, 20832, 20384, 20832, 20832, 21056, 20368, 20832, 21056, 20384, 20832, 21056, 20384, 20832, 21056, 20832, 20384, 21056, 20832, 20384, 21056, 20832, 20832] 0 27 rfl (fun j hj => by
    interval_cases j
    exacts [gc0_0, gc0_1, gc0_2, gc0_3, gc0_4, gc0_5, gc0_6, gc0_7, gc0_8, gc0_9, gc0_10, gc0_11, gc0_12, gc0_13, gc0_14, gc0_15, gc0_16, gc0_17, gc0_18, gc0_19, gc0_20, gc0_21, gc0_22, gc0_23, gc0_24, gc0_25, gc0_26])
  rw [show ([20384, 20832, 20788, 20380, 20832, 20832, 20384, 20832, 20832, 21056, 20368, 20832, 21056, 20384, 20832, 21056, 20384, 20832, 21056, 20832, 20384, 21056, 20832, 20384, 21056, 20832, 20832] : List ℕ).sum = 560160 by decide] at h
  exact h
theorem rem_c3 : cntR rt 3 0 NK = 23328 := by
  have h := cntR_list rt 3 120 [672, 1792, 224, 672, 2240, 224, 672, 1120, 1344, 224, 224, 2464, 0, 448, 2240, 224, 672, 896, 1120, 448, 0, 2240, 0, 224, 1824, 448, 672] 0 27 rfl (fun j hj => by
    interval_cases j
    exacts [gc3_0, gc3_1, gc3_2, gc3_3, gc3_4, gc3_5, gc3_6, gc3_7, gc3_8, gc3_9, gc3_10, gc3_11, gc3_12, gc3_13, gc3_14, gc3_15, gc3_16, gc3_17, gc3_18, gc3_19, gc3_20, gc3_21, gc3_22, gc3_23, gc3_24, gc3_25, gc3_26])
  rw [show ([672, 1792, 224, 672, 2240, 224, 672, 1120, 1344, 224, 224, 2464, 0, 448, 2240, 224, 672, 896, 1120, 448, 0, 2240, 0, 224, 1824, 448, 672] : List ℕ).sum = 23328 by decide] at h
  exact h
theorem rem_c5 : cntR rt 5 0 NK = 4320 := by
  have h := cntR_list rt 5 120 [672, 0, 44, 676, 0, 0, 672, 0, 0, 0, 464, 0, 0, 448, 0, 0, 448, 0, 0, 0, 448, 0, 0, 448, 0, 0, 0] 0 27 rfl (fun j hj => by
    interval_cases j
    exacts [gc5_0, gc5_1, gc5_2, gc5_3, gc5_4, gc5_5, gc5_6, gc5_7, gc5_8, gc5_9, gc5_10, gc5_11, gc5_12, gc5_13, gc5_14, gc5_15, gc5_16, gc5_17, gc5_18, gc5_19, gc5_20, gc5_21, gc5_22, gc5_23, gc5_24, gc5_25, gc5_26])
  rw [show ([672, 0, 44, 676, 0, 0, 672, 0, 0, 0, 464, 0, 0, 448, 0, 0, 448, 0, 0, 0, 448, 0, 0, 448, 0, 0, 0] : List ℕ).sum = 4320 by decide] at h
  exact h
theorem rem_c48 : cntR rt 48 0 NK = 137952 := by
  have h := cntR_list rt 48 120 [5152, 4256, 5824, 5152, 3808, 5824, 5152, 4928, 4704, 5600, 5824, 3584, 5824, 5600, 3808, 5600, 5376, 5152, 4704, 5600, 6048, 3584, 6048, 5824, 4000, 5600, 5376] 0 27 rfl (fun j hj => by
    interval_cases j
    exacts [gc48_0, gc48_1, gc48_2, gc48_3, gc48_4, gc48_5, gc48_6, gc48_7, gc48_8, gc48_9, gc48_10, gc48_11, gc48_12, gc48_13, gc48_14, gc48_15, gc48_16, gc48_17, gc48_18, gc48_19, gc48_20, gc48_21, gc48_22, gc48_23, gc48_24, gc48_25, gc48_26])
  rw [show ([5152, 4256, 5824, 5152, 3808, 5824, 5152, 4928, 4704, 5600, 5824, 3584, 5824, 5600, 3808, 5600, 5376, 5152, 4704, 5600, 6048, 3584, 6048, 5824, 4000, 5600, 5376] : List ℕ).sum = 137952 by decide] at h
  exact h
/-- **the remainder histogram** of S5H: any function of the remainders summed over all host indices. -/
theorem rem_hist {M : Type*} [AddCommMonoid M] (f : ℕ → M) :
    ∑ i ∈ Finset.range (224 * NK), f (remAt rt i) = (560160 : ℕ) • f 0 + (23328 : ℕ) • f 3 + (4320 : ℕ) • f 5 + (137952 : ℕ) • f 48 := by
  rw [rem_sum rt NK [0, 3, 5, 48] (by decide) rem_in f]
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rem_c0, rem_c3, rem_c5, rem_c48, add_zero, add_assoc]
/-- the residual sum of provider 0. -/
theorem res_sum0 : ∑ q ∈ Finset.range R, resAt resT 0 q = 3304 := by
  rw [← resS_eq]; decide +kernel
/-- the residual sum of provider 1. -/
theorem res_sum1 : ∑ q ∈ Finset.range R, resAt resT 1 q = 3304 := by
  rw [← resS_eq]; decide +kernel
end OAI.PowerSaving.XHostT.K.DS5H
#print axioms OAI.PowerSaving.XHostT.K.DS5H.rem_hist
#print axioms OAI.PowerSaving.XHostT.K.DS5H.res_sum0
#print axioms OAI.PowerSaving.XHostT.K.DS5H.res_sum1
