import Work.Mixed.Prog.LoadRAM

/-!
# Mixed engine, R3 program part 12: the source accumulator (agent e-r3prog)

A class returns, besides its own answers, the answered supply vectors its subtree consumed, in
consumption order, with their sizes.  After each child call it appends the child's own source
slots (in slot order) and the child's consumed list.

* `filterMap_range_eq`   the source sizes in slot order are `slotSize (nth srcP k)`;
* `concat3`              three tapes concatenated (pointer copies, linear in the number of entries).
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

universe U

lemma filterMap_range_eq (P : ℕ → Prop) [DecidablePred P] (f : ℕ → ℕ) (N : ℕ) :
    (List.range N).filterMap (fun t => if P t then some (f t) else none) =
      (List.range (Nat.count P N)).map (fun k => f (Nat.nth P k)) := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [List.range_succ, List.filterMap_append, ih, Nat.count_succ]
    by_cases hp : P N
    · have e1 : List.filterMap (fun t => if P t then some (f t) else none) [N] = [f N] := by simp [hp]
      rw [e1, if_pos hp, List.range_succ, List.map_append]
      simp [Nat.nth_count hp]
    · simp [hp]

lemma getD_map_range (g : ℕ → ℕ) (N k : ℕ) (hk : k < N) : ((List.range N).map g).getD k 0 = g k := by
  rw [List.getD_eq_getElem _ _ (by simpa using hk)]
  simp

section Concat
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s t : Ty} {x : ι → T s}

/-- Three tapes concatenated. -/
def cat3 {μ : Type*} (L1 L2 L3 : ℕ) (f1 f2 f3 : ℕ → μ) (k : ℕ) : μ :=
  if k < L1 then f1 k else if k < L1 + L2 then f2 (k - L1) else f3 (k - L1 - L2)

lemma concat3 {L1 L2 L3 : ι → ℕ} {f1 f2 f3 : ι → ℕ → T t}
    (h1 : Has m B s (a t) x (fun i => (Layout.std (L1 i)).tape (fun k => f1 i k)))
    (h2 : Has m B s (a t) x (fun i => (Layout.std (L2 i)).tape (fun k => f2 i k)))
    (h3 : Has m B s (a t) x (fun i => (Layout.std (L3 i)).tape (fun k => f3 i k)))
    (k1 : Knows m B s x L1) (k2 : Knows m B s x L2) (k3 : Knows m B s x L3)
    (s1 : Small B L1) (s2 : Small B L2) (s3 : Small B L3) :
    Can m B s (a t) x (fun i => (Layout.std (L1 i + L2 i + L3 i)).tape (fun k => cat3 (L1 i) (L2 i) (L3 i) (f1 i) (f2 i) (f3 i) k))
      (fun i => L1 i + L2 i + L3 i) := by
  let N (i : ι) := L1 i + L2 i + L3 i
  have hN : Knows m B s x N := (k1.plus k2 s1 s2).plus k3 (s1.add s2) s3
  have hNs : Small B N := (s1.add s2).add s3
  let UN := UnitF N
  let yU (u : UN) : T (p s w) := (x u.1, u.2.val)
  have hk : Knows m (fun u : UN => B u.1) (p s w) yU (fun u => u.2.val) := Can.second
  have hks : Small (fun u : UN => B u.1) (fun u => u.2.val) := (hNs.reindex (fun u : UN => u.1)).mono (fun u => u.2.2.le)
  have K1 : Knows m (fun u : UN => B u.1) (p s w) yU (fun u => L1 u.1) := Can.first.then_do (k1.reindex _)
  have K2 : Knows m (fun u : UN => B u.1) (p s w) yU (fun u => L2 u.1) := Can.first.then_do (k2.reindex _)
  have S1 : Small (fun u : UN => B u.1) (fun u => L1 u.1) := s1.reindex _
  have S2 : Small (fun u : UN => B u.1) (fun u => L2 u.1) := s2.reindex _
  have T1U : Has m (fun u : UN => B u.1) (p s w) (a t) yU (fun u => (Layout.std (L1 u.1)).tape (fun k => f1 u.1 k)) :=
    Can.first.then_do (h1.reindex (fun u : UN => u.1))
  have T2U : Has m (fun u : UN => B u.1) (p s w) (a t) yU (fun u => (Layout.std (L2 u.1)).tape (fun k => f2 u.1 k)) :=
    Can.first.then_do (h2.reindex (fun u : UN => u.1))
  have T3U : Has m (fun u : UN => B u.1) (p s w) (a t) yU (fun u => (Layout.std (L3 u.1)).tape (fun k => f3 u.1 k)) :=
    Can.first.then_do (h3.reindex (fun u : UN => u.1))
  have hent : Has m (fun u : UN => B u.1) (p s w) t yU
      (fun u => cat3 (L1 u.1) (L2 u.1) (L3 u.1) (f1 u.1) (f2 u.1) (f3 u.1) u.2.val) := by
    apply Can.testz (hk.wlt K1 |> fun hh => ((Can.wc 1).wsub hh ((Small.const _ 1).mono (fun _ => Nat.sub_le _ _))))
    · -- k < L1
      let K := {u : UN // 1 - (if u.2.val < L1 u.1 then 1 else 0) = 0}
      have hlt : ∀ κ : K, κ.1.2.val < L1 κ.1.1 := fun κ => by
        have := κ.2; by_contra hc; simp [hc] at this
      refine (Can.fetch (t:=t) (L:=fun κ : K => Layout.std (L1 κ.1.1)) (v:=fun κ k => f1 κ.1.1 k)
        (i:=fun κ => ⟨κ.1.2.val, hlt κ⟩) (T1U.reindex (fun κ : K => κ.1))
        (hk.reindex (fun κ : K => κ.1))).cong (fun _ => rfl) (fun κ => ?_)
      simp [cat3, hlt κ]
    · apply Can.testz (((hk.wlt (K1.plus K2 S1 S2)) |> fun hh => ((Can.wc 1).wsub hh
        ((Small.const _ 1).mono (fun _ => Nat.sub_le _ _)))).reindex Subtype.val)
      · let K := {κ : {u : UN // 1 - (if u.2.val < L1 u.1 then 1 else 0) ≠ 0} //
          1 - (if κ.1.2.val < L1 κ.1.1 + L2 κ.1.1 then 1 else 0) = 0}
        have hge : ∀ κ : K, L1 κ.1.1.1 ≤ κ.1.1.2.val := fun κ => by
          have := κ.1.2; by_contra hc; simp [Nat.lt_of_not_le hc] at this
        have hlt : ∀ κ : K, κ.1.1.2.val - L1 κ.1.1.1 < L2 κ.1.1.1 := fun κ => by
          have := κ.2
          have h2 : κ.1.1.2.val < L1 κ.1.1.1 + L2 κ.1.1.1 := by by_contra hc; simp [hc] at this
          have := hge κ; omega
        refine (Can.fetch (t:=t) (L:=fun κ : K => Layout.std (L2 κ.1.1.1)) (v:=fun κ k => f2 κ.1.1.1 k)
          (i:=fun κ => ⟨κ.1.1.2.val - L1 κ.1.1.1, hlt κ⟩)
          (T2U.reindex (fun κ : K => κ.1.1))
          ((hk.wsub K1 hks.sub).reindex (fun κ : K => κ.1.1))).cong (fun _ => rfl) (fun κ => ?_)
        have h1 : ¬ κ.1.1.2.val < L1 κ.1.1.1 := Nat.not_lt.mpr (hge κ)
        have h2 : κ.1.1.2.val < L1 κ.1.1.1 + L2 κ.1.1.1 := by have := hlt κ; omega
        simp [cat3, h1, h2]
      · let K := {κ : {u : UN // 1 - (if u.2.val < L1 u.1 then 1 else 0) ≠ 0} //
          1 - (if κ.1.2.val < L1 κ.1.1 + L2 κ.1.1 then 1 else 0) ≠ 0}
        have hge : ∀ κ : K, L1 κ.1.1.1 + L2 κ.1.1.1 ≤ κ.1.1.2.val := fun κ => by
          have := κ.2; by_contra hc; simp [Nat.lt_of_not_le hc] at this
        have hlt : ∀ κ : K, κ.1.1.2.val - L1 κ.1.1.1 - L2 κ.1.1.1 < L3 κ.1.1.1 := fun κ => by
          have := κ.1.1.2.2; have := hge κ; simp only [N] at *; omega
        refine (Can.fetch (t:=t) (L:=fun κ : K => Layout.std (L3 κ.1.1.1)) (v:=fun κ k => f3 κ.1.1.1 k)
          (i:=fun κ => ⟨κ.1.1.2.val - L1 κ.1.1.1 - L2 κ.1.1.1, hlt κ⟩)
          (T3U.reindex (fun κ : K => κ.1.1))
          (((hk.wsub K1 hks.sub).wsub K2 hks.sub.sub).reindex (fun κ : K => κ.1.1))).cong
          (fun _ => rfl) (fun κ => ?_)
        have h1 : ¬ κ.1.1.2.val < L1 κ.1.1.1 := by have := hge κ; omega
        have h2 : ¬ κ.1.1.2.val < L1 κ.1.1.1 + L2 κ.1.1.1 := Nat.not_lt.mpr (hge κ)
        simp [cat3, h1, h2]
  have H := Can.layout (t:=t) (L:=fun i => Layout.std (N i)) (v:=fun i k => cat3 (L1 i) (L2 i) (L3 i) (f1 i) (f2 i) (f3 i) k)
    (P:=fun _ => 0) (W:=fun _ => 0) hN (hent.mono (fun _ => Nat.zero_le _)) hNs
  refine H.mono (fun i => ?_)
  change 0 + N i + N i * 0 ≤ N i
  omega

end Concat

end
end PowerSaving.RAM
end OAI
