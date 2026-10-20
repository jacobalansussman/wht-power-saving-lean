import Work.XHostT.KHist
import Work.XHostT.Data.H2979.GlobC0
import Work.XHostT.Data.H2979.GlobC1
import Work.XHostT.Data.H2979.GlobC2
import Work.XHostT.Data.H2979.GlobC3
import Work.XHostT.Data.H2979.GlobC4
import Work.XHostT.Data.H2979.GlobC5
/-! p-hostck mk_parts.py: global facts of H2979: remainder histogram {0: 97000, 1: 17310, 3: 4956, 5: 1728, 48: 28190}, residual sums {0: 3283, 1: 3283} (6 chunks of 111 kinds). -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979
open OAI.PowerSaving.XHostT.K
/-- every remainder is one of [0, 1, 3, 5, 48]. -/
theorem rem_in : inR rt [0, 1, 3, 5, 48] 0 NK = true :=
  inR_list rt [0, 1, 3, 5, 48] 111 6 0 (fun j hj => by
    interval_cases j
    exacts [gin0, gin1, gin2, gin3, gin4, gin5])
theorem rem_c0 : cntR rt 0 0 NK = 97000 := by
  have h := cntR_list rt 0 111 [15158, 16647, 16695, 15158, 16647, 16695] 0 6 rfl (fun j hj => by
    interval_cases j
    exacts [gc0_0, gc0_1, gc0_2, gc0_3, gc0_4, gc0_5])
  rw [show ([15158, 16647, 16695, 15158, 16647, 16695] : List ℕ).sum = 97000 by decide] at h
  exact h
theorem rem_c1 : cntR rt 1 0 NK = 17310 := by
  have h := cntR_list rt 1 111 [3504, 2581, 2570, 3504, 2581, 2570] 0 6 rfl (fun j hj => by
    interval_cases j
    exacts [gc1_0, gc1_1, gc1_2, gc1_3, gc1_4, gc1_5])
  rw [show ([3504, 2581, 2570, 3504, 2581, 2570] : List ℕ).sum = 17310 by decide] at h
  exact h
theorem rem_c3 : cntR rt 3 0 NK = 4956 := by
  have h := cntR_list rt 3 111 [384, 533, 1561, 384, 533, 1561] 0 6 rfl (fun j hj => by
    interval_cases j
    exacts [gc3_0, gc3_1, gc3_2, gc3_3, gc3_4, gc3_5])
  rw [show ([384, 533, 1561, 384, 533, 1561] : List ℕ).sum = 4956 by decide] at h
  exact h
theorem rem_c5 : cntR rt 5 0 NK = 1728 := by
  have h := cntR_list rt 5 111 [826, 37, 1, 826, 37, 1] 0 6 rfl (fun j hj => by
    interval_cases j
    exacts [gc5_0, gc5_1, gc5_2, gc5_3, gc5_4, gc5_5])
  rw [show ([826, 37, 1, 826, 37, 1] : List ℕ).sum = 1728 by decide] at h
  exact h
theorem rem_c48 : cntR rt 48 0 NK = 28190 := by
  have h := cntR_list rt 48 111 [4992, 5066, 4037, 4992, 5066, 4037] 0 6 rfl (fun j hj => by
    interval_cases j
    exacts [gc48_0, gc48_1, gc48_2, gc48_3, gc48_4, gc48_5])
  rw [show ([4992, 5066, 4037, 4992, 5066, 4037] : List ℕ).sum = 28190 by decide] at h
  exact h
/-- **the remainder histogram** of H2979: any function of the remainders summed over all host indices. -/
theorem rem_hist {M : Type*} [AddCommMonoid M] (f : ℕ → M) :
    ∑ i ∈ Finset.range (224 * NK), f (remAt rt i) = (97000 : ℕ) • f 0 + (17310 : ℕ) • f 1 + (4956 : ℕ) • f 3 + (1728 : ℕ) • f 5 + (28190 : ℕ) • f 48 := by
  rw [rem_sum rt NK [0, 1, 3, 5, 48] (by decide) rem_in f]
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, rem_c0, rem_c1, rem_c3, rem_c5, rem_c48, add_zero, add_assoc]
/-- the residual sum of provider 0. -/
theorem res_sum0 : ∑ q ∈ Finset.range R, resAt resT 0 q = 3283 := by
  rw [← resS_eq]; decide +kernel
/-- the residual sum of provider 1. -/
theorem res_sum1 : ∑ q ∈ Finset.range R, resAt resT 1 q = 3283 := by
  rw [← resS_eq]; decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979
#print axioms OAI.PowerSaving.XHostT.K.DH2979.rem_hist
#print axioms OAI.PowerSaving.XHostT.K.DH2979.res_sum0
#print axioms OAI.PowerSaving.XHostT.K.DH2979.res_sum1
