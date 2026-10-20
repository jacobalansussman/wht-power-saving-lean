import XHost.Transport
import Work.GFrame.Clifford.Core
import Work.GFrame.Labels.BXMap

/-!
# x-hosting, L2 in the coordinate (bank) form: a hosted idle climb is a bank

The bank layer of the #401 layout (`Work/Bank401/Adapt.lean`, key bank-401: `bank_adapt` K6, `bank_pad` K7; entrance
gauges and residual banks Dugongue's, CrocSwap/integer-mult-bounds PR #401; bank patterns Lydakis's, #275) keeps a bank's
frame as `A * pk U * B` (`A` free, `U` the coordinates already climbed, `B` the start frame) and receives a gauged
residual `C * pk R * C'` (`C`, `C'` free; `gauged_Wn` there).  Idle-axis hosting (DaysSky, PR #442, PROOF.md §2) puts
such a residual on a HOST: an idle data role whose climb `M ↦ D * pk S * D' * M` (`D = kg G`, `D' = kg G⁻¹`, `S` the
pending coordinates; `idle_climb_factor`) has not run yet.  The L3 shortcut treats the host as a bank whose start frame
is the host's frame `M` before its climb and whose end frame is the frame after it; the remainder is the pad:

* `exists_pk_transport`  **the transport σ, coordinate form**: for `|R| = |T|` there is an additive bijection `σ` with
                         `permMat σ⁻¹ * pk R * permMat σ = pk T` (from `exists_transport`, `transport_blockMat`);
* `host_start`           `M = D * pk ∅ * (D' * M)` when `D * D' = 1`;
* `host_adapt`           **K6 for a host**: if `U ⊆ S` and `|U| + |R| ≤ |S|` a FREE `P` before the invocation carries the
                         residual coordinates `R` onto pending coordinates `T ⊆ S \ U`:
                         `C * pk R * C' * (P * (A * pk U * B)) = A2 * pk U2 * B`, `A2` free, `U ⊆ U2 ⊆ S`,
                         `|U2| = |U| + |R|`;
* `host_pad`             **K7 for a host, the remainder climb**: from `A * pk U * B` (`U ⊆ S`) to `D * pk S * B` is ONE
                         block of rank `|S| - |U|`;
* `host_full`            fully hosted (`U = S`): the remainder is free;
* `HostOK`, `hostOK_start`, `hostOK_use`, `hostOK_pad`, `hostOK_full`   the same as an invariant for a generator: the
                         host's state after `n` hosted coordinates; one lemma per hosted use; the remainder at the end;
* `idle_host`            a host at `X.Φ ⟨G, s⟩` climbing to `X.Φ ⟨G, t⟩` starts in `HostOK (liftS t \ liftS s) _ 0`;
* `idle_climb_factor`    the host's idle climb in one base: `X.Φ ⟨G, t⟩ = kg * pk (liftS t \ liftS s) * (kg⁻¹ * X.Φ ⟨G, s⟩)`
                         with `|liftS t \ liftS s| = |t \ s|` and `kg * kg⁻¹ = 1`.

The engine (`pk`, `kg`, `core`, `Free`, `free_inv`, `gstep_of_factor`, `core_nested`) is Jacob Sussman's; `host_adapt`,
`host_pad`, `host_full` follow the proofs of `bank_adapt`, `bank_pad`, `bank_full` (bank-401) with the target set `S`
of the host in place of `univ`, and a linear transport in place of the coordinate permutation.
No `sorry`.  Written for the wht43 x-hosting study, tier L2 (13 Oct 2026), prepared with assistance from Claude.
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.GF
namespace XHost
open Binary Matrix Finset RAM Module
noncomputable section

section PkTransport
variable {α : Type*} [Fintype α] [DecidableEq α]

/-- **The transport σ, coordinate form.**  The kernel on any `|R|` coordinates is conjugated by a free address
permutation onto the kernel on any other `|R|` coordinates. -/
theorem exists_pk_transport (R T : Finset α) (h : R.card = T.card) :
    ∃ σ : APerm α, permMat (apInv σ) * pk R * permMat σ = pk T := by
  let z : Fin T.card → Space α := fun i => (OBase.canonical α).fam R (finCongr h.symm i)
  have hz : LinearIndependent F z :=
    ((OBase.canonical α).fam_indep R).comp _ (finCongr h.symm).injective
  obtain ⟨σ, hσ⟩ := exists_transport z ((OBase.canonical α).fam T) hz ((OBase.canonical α).fam_indep T)
  refine ⟨σ, ?_⟩
  have hR : pk R = blockMat z := by
    rw [pk_eq_blockMat]
    exact (blockMat_comp_equiv _ (finCongr h.symm)).symm
  rw [hR, transport_blockMat hσ, ← pk_eq_blockMat]

/-- the transported kernel commutes past the transport: `pk R * permMat σ = permMat σ * pk T`. -/
lemma pk_transport_comm {R T : Finset α} {σ : APerm α} (hσ : permMat (apInv σ) * pk R * permMat σ = pk T) :
    pk R * permMat σ = permMat σ * pk T := by
  rw [← hσ]
  simp only [← Matrix.mul_assoc, permMat_mul_inv, Matrix.one_mul]

/-- **The start of a hosted climb**: the host's frame `M` is `D * pk ∅ * (D' * M)`. -/
theorem host_start {D D' : CMat α} (hDD : D * D' = 1) (M : CMat α) : M = D * pk (∅ : Finset α) * (D' * M) := by
  rw [pk_empty, Matrix.mul_one, ← Matrix.mul_assoc, hDD, Matrix.one_mul]

/-- **K6 for a host, the hosted use.**  A free `P` before the invocation carries the residual coordinates `R` of the
gauged step onto pending coordinates of the host (inside `S`, off the coordinates `U` already hosted). -/
theorem host_adapt {A B C C' : CMat α} (hA : Free A) (hC : Free C) (hC' : Free C')
    (S U R : Finset α) (hUS : U ⊆ S) (hc : U.card + R.card ≤ S.card) :
    ∃ P, Free P ∧ ∃ (A2 : CMat α) (U2 : Finset α), Free A2 ∧ U ⊆ U2 ∧ U2 ⊆ S ∧
      U2.card = U.card + R.card ∧ C * pk R * C' * (P * (A * pk U * B)) = A2 * pk U2 * B := by
  obtain ⟨Ai, hAi, hAiA, -⟩ := free_inv hA
  obtain ⟨Ci, hCi, -, hC'Ci⟩ := free_inv hC'
  have hUc : R.card ≤ (S \ U).card := by
    rw [Finset.card_sdiff_of_subset hUS]
    omega
  obtain ⟨T, hTU, hTc⟩ := Finset.exists_subset_card_eq hUc
  obtain ⟨σ, hσ⟩ := exists_pk_transport R T hTc.symm
  have hPi := pk_transport_comm hσ
  have hdisj : Disjoint T U := by
    rw [Finset.disjoint_left]
    intro a ha hu
    exact (Finset.mem_sdiff.mp (hTU ha)).2 hu
  have hTS : T ⊆ S := fun a ha => (Finset.mem_sdiff.mp (hTU ha)).1
  refine ⟨Ci * permMat σ * Ai, Free.mul (Free.mul hCi (Free.perm σ)) hAi,
    C * permMat σ, T ∪ U, Free.mul hC (Free.perm σ), Finset.subset_union_right,
    Finset.union_subset hTS hUS, ?_, ?_⟩
  · rw [Finset.card_union_of_disjoint hdisj, hTc, Nat.add_comm]
  · have e : C * pk R * C' * (Ci * permMat σ * Ai * (A * pk U * B))
        = C * (pk R * permMat σ) * (Ai * A) * pk U * B := by
      simp only [Matrix.mul_assoc]
      rw [← Matrix.mul_assoc C' Ci, hC'Ci, Matrix.one_mul]
    rw [e, hAiA, Matrix.mul_one, hPi, pk_union T U hdisj]
    simp only [Matrix.mul_assoc]

end PkTransport

section Pad
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- **K7 for a host, the remainder climb**: from `A * pk U * B` (hosted coordinates `U ⊆ S`) to the host's target
`D * pk S * B` is ONE block of rank `|S| - |U|`. -/
theorem host_pad (l : ρ) {A D : CMat α} (hA : Free A) (hD : Free D) (B : CMat α) (S U : Finset α)
    (hUS : U ⊆ S) (h1 : U.card < S.card) (h2 : S.card - U.card < Fintype.card α) :
    GStep l (A * pk U * B) (D * pk S * B) (fun φ => φ (S.card - U.card)) := by
  obtain ⟨Ai, hAi, hAiA, -⟩ := free_inv hA
  have hk : (S \ U).card = S.card - U.card := Finset.card_sdiff_of_subset hUS
  have h := gstep_of_factor l (A * pk U * B) ((OBase.canonical α).fam (S \ U))
    ((OBase.canonical α).fam_indep _) hD hAi (by omega) (by omega)
  have e : D * blockMat ((OBase.canonical α).fam (S \ U)) * Ai * (A * pk U * B) = D * pk S * B := by
    rw [← pk_eq_blockMat, pk_sdiff hUS]
    have e2 : D * pk (S \ U) * Ai * (A * pk U * B) = D * pk (S \ U) * (Ai * A) * pk U * B := by
      simp only [Matrix.mul_assoc]
    rw [e2, hAiA, Matrix.mul_one]
    simp only [Matrix.mul_assoc]
  rw [e, hk] at h
  exact h

/-- **A fully hosted climb**: from `A * pk S * B` to `D * pk S * B` is free. -/
theorem host_full (l : ρ) {A D : CMat α} (hA : Free A) (hD : Free D) (B : CMat α) (S : Finset α) :
    GStep l (A * pk S * B) (D * pk S * B) (fun _ => 0) := by
  obtain ⟨Ai, hAi, hAiA, -⟩ := free_inv hA
  have h := gstep_free (Free.mul hD hAi) l (A * pk S * B)
  have e : D * Ai * (A * pk S * B) = D * pk S * B := by
    calc D * Ai * (A * pk S * B) = D * (Ai * A) * pk S * B := by simp only [Matrix.mul_assoc]
      _ = D * pk S * B := by rw [hAiA, Matrix.mul_one]
  rw [e] at h
  exact h

end Pad

section Chain
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- **The state of a host** with pending coordinates `S` and fixed tail `B`: its frame is `A * pk U * B` with `A`
free, `U ⊆ S` the coordinates already hosted, `|U| = n` (the bank invariant of bank-401 with `S` for `univ`). -/
def HostOK (S : Finset α) (B : CMat α) (n : ℕ) (Z : CMat α) : Prop :=
  ∃ (A : CMat α) (U : Finset α), Free A ∧ U ⊆ S ∧ U.card = n ∧ Z = A * pk U * B

/-- the host before its first hosted use (its frame `M`, nothing hosted). -/
theorem hostOK_start {D D' : CMat α} (hDD : D * D' = 1) (hD : Free D) (S : Finset α) (M : CMat α) :
    HostOK S (D' * M) 0 M :=
  ⟨D, ∅, hD, Finset.empty_subset S, Finset.card_empty, host_start hDD M⟩

/-- **One hosted use**: a free `P` on the host, then the invocation's gauged residual `W = C * pk R * C'` on the slot
input `P * Z`; the host stays in its state with `|R|` more coordinates hosted. -/
theorem hostOK_use {S : Finset α} {B Z W C C' : CMat α} {n : ℕ} (h : HostOK S B n Z) (hC : Free C)
    (hC' : Free C') (R : Finset α) (hW : W = C * pk R * C') (hc : n + R.card ≤ S.card) :
    ∃ P, Free P ∧ HostOK S B (n + R.card) (W * (P * Z)) := by
  obtain ⟨A, U, hA, hUS, hU, rfl⟩ := h
  subst hW
  obtain ⟨P, hP, A2, U2, hA2, -, hU2S, hU2, he⟩ :=
    host_adapt (B := B) hA hC hC' S U R hUS (by rw [hU]; exact hc)
  exact ⟨P, hP, A2, U2, hA2, hU2S, by rw [hU2, hU], he⟩

/-- **The remainder of a partly hosted climb**: ONE block of rank `|S| - n` to the host's target `D * pk S * B`. -/
theorem hostOK_pad (l : ρ) {S : Finset α} {B Z D : CMat α} {n : ℕ} (h : HostOK S B n Z) (hD : Free D)
    (h1 : n < S.card) (h2 : S.card - n < Fintype.card α) :
    GStep l Z (D * pk S * B) (fun φ => φ (S.card - n)) := by
  obtain ⟨A, U, hA, hUS, hU, rfl⟩ := h
  subst hU
  exact host_pad l hA hD B S U hUS h1 h2

/-- **A fully hosted climb**: nothing left to pay. -/
theorem hostOK_full (l : ρ) {S : Finset α} {B Z D : CMat α} (h : HostOK S B S.card Z) (hD : Free D) :
    GStep l Z (D * pk S * B) (fun _ => 0) := by
  obtain ⟨A, U, hA, hUS, hU, rfl⟩ := h
  have e : U = S := Finset.eq_of_subset_of_card_le hUS (by rw [hU])
  rw [e]
  exact host_full l hA hD B S

end Chain

section Idle
variable {H α : Type} [Fintype H] [DecidableEq H] [Fintype α] [DecidableEq α]

/-- **The host's idle climb in one base** (`core_nested` lifted through a `BXMap`): with `D = kg (liftG G)`,
`D' = kg (liftG G)⁻¹`, `S = liftS t \ liftS s`: `X.Φ ⟨G, t⟩ = D * pk S * (D' * X.Φ ⟨G, s⟩)`, `|S| = |t \ s|`, `D * D' = 1`. -/
theorem idle_climb_factor (X : BXMap H α) (G : APerm H) {s t : Finset H} (h : s ⊆ t) :
    X.Φ ⟨G, t⟩ = kg (X.liftG G) * pk (X.liftS t \ X.liftS s) * (kg (apInv (X.liftG G)) * X.Φ ⟨G, s⟩) ∧
      (X.liftS t \ X.liftS s).card = (t \ s).card ∧ kg (X.liftG G) * kg (apInv (X.liftG G)) = 1 := by
  refine ⟨?_, X.card h, kg_mul_inv _⟩
  show core (X.liftG G) (X.liftS t) = _ * (_ * core (X.liftG G) (X.liftS s))
  rw [core_nested _ (X.mono h), ← pk_eq_blockMat]
  simp only [Matrix.mul_assoc]

/-- **A host on its idle climb, ready for hosted uses**: the host's frame `X.Φ ⟨G, s⟩` is in state `HostOK S B 0`
with `S = liftS t \ liftS s` (`|S| = |t \ s|`) and its climb target is `kg * pk S * B`. -/
theorem idle_host (X : BXMap H α) (G : APerm H) {s t : Finset H} (h : s ⊆ t) :
    HostOK (X.liftS t \ X.liftS s) (kg (apInv (X.liftG G)) * X.Φ ⟨G, s⟩) 0 (X.Φ ⟨G, s⟩) ∧
      X.Φ ⟨G, t⟩ = kg (X.liftG G) * pk (X.liftS t \ X.liftS s) * (kg (apInv (X.liftG G)) * X.Φ ⟨G, s⟩) ∧
      Free (kg (X.liftG G)) ∧ (X.liftS t \ X.liftS s).card = (t \ s).card := by
  obtain ⟨e, hc, hDD⟩ := idle_climb_factor X G h
  exact ⟨hostOK_start hDD (free_kg _) _ _, e, free_kg _, hc⟩

end Idle
end
end XHost
end PowerSaving.GF
end OAI

#print axioms OAI.PowerSaving.GF.XHost.exists_pk_transport
#print axioms OAI.PowerSaving.GF.XHost.host_start
#print axioms OAI.PowerSaving.GF.XHost.host_adapt
#print axioms OAI.PowerSaving.GF.XHost.host_pad
#print axioms OAI.PowerSaving.GF.XHost.host_full
#print axioms OAI.PowerSaving.GF.XHost.idle_climb_factor
#print axioms OAI.PowerSaving.GF.XHost.hostOK_use
#print axioms OAI.PowerSaving.GF.XHost.hostOK_pad
#print axioms OAI.PowerSaving.GF.XHost.hostOK_full
#print axioms OAI.PowerSaving.GF.XHost.idle_host
