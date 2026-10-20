import XHost.Core
import Work.GFrame.Clifford.Basic

/-!
# x-hosting, L2: the transport σ

Idle-axis hosting (DaysSky, CrocSwap/integer-mult-bounds PR #442, research/idle-axis-hosting/PROOF.md §2, item 2) lets
a helper residual of an invocation be owned by an idle data stream "(F, d·σ⁻¹)", where σ is a coordinate isometry
carrying the local residual block onto a set S' of w PENDING coordinates of the host's idle climb ("composed with
U_q ⊕ I" when a gauge is involved).  In our engine (Jacob Sussman's GFrame engine, wht-power-saving-lean, Apache-2.0,
as used unchanged by Ryan Shea's complex588) this is:

* `permMat_conj_shift`, `permMat_conj_dir`, `permMat_conj_blockMat`, `permMat_inv_conj_blockMat`
      **an address permutation conjugates `dir` into `dir` along the mapped direction**:
      `permMat G⁻¹ * blockMat z * permMat G = blockMat (G ∘ z)`;
* `exists_transport`, `exists_transport_into`
      **the transport σ**: for ANY linearly independent `z` (the images of the missing label coordinates of a slot,
      L1) and ANY `r` linearly independent targets `u` (e.g. `r` of the host's pending directions `v ∘ e`, `e`
      injective) there is a FREE address permutation `σ` (an additive bijection) with `σ (z i) = u i`;
* `transport_twisted`
      a gauged residual `C * blockMat z * C'` (free `C`, `C'`, `C' * C = 1`, the L1 form) is carried by free
      `T = D * permMat σ⁻¹ * C'`, `T' = C * permMat σ * D'` (`T' * T = T * T' = 1`) onto the host's gauged factor
      `D * blockMat u * D'` (free `D`, `D'`, `D' * D = 1`; e.g. `D = kg G`, `D' = kg G⁻¹` of `core_nested`);
* `host_stack`     several hosted residuals on one host stack: `blockMat (Fin.append a u)`;
* `blockMat_comp_equiv`, `blockMat_split`   a block does not depend on the order of its directions;
* `hosted_climb`, `hosted_climb_full`, `hosted_climb_of_split`
      the host's route around an invocation: free `T'` before, the invocation's residual on the slot, free `T`
      after, then ONE remaining block of rank `s` (or nothing when `s = 0`) reaches the host's old target
      `D * blockMat v * D' * M`.

No `sorry`.  Written for the wht43 x-hosting study, tier L2 (13 Oct 2026), prepared with assistance from Claude.
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.GF
namespace XHost
open Binary Matrix Finset RAM Module
noncomputable section

section Conj
variable {α : Type*} [Fintype α] [DecidableEq α]

lemma permMat_apInv_apInv (G : APerm α) : permMat (apInv (apInv G)) = permMat G := by
  ext x y
  by_cases h : G.π x = y <;> simp [permMat, apInv, h]

lemma apInv_apInv_apply (G : APerm α) (x : Space α) : (apInv (apInv G)).π x = G.π x := by
  simp [apInv]

/-- **An address permutation conjugates a translation into a translation** along the mapped vector. -/
theorem permMat_conj_shift (G : APerm α) (z : Space α) :
    permMat G * shift z * permMat (apInv G) = shift ((apInv G).π z) := by
  ext x y
  rw [mul_permMat_apply, permMat_mul_apply]
  have hs : (apInv G).π.symm y = G.π y := by simp [apInv]
  rw [hs]
  simp only [shift]
  by_cases h : x + (apInv G).π z = y
  · have h' : G.π x + z = G.π y := by rw [← h, G.add, apply_apInv]
    rw [ite_eq_left h', ite_eq_left h]
  · have h' : ¬ G.π x + z = G.π y := by
      intro h'
      apply h
      have h'' := congrArg (apInv G).π h'
      rwa [(apInv G).add, apInv_apply, apInv_apply] at h''
    rw [ite_eq_right h', ite_eq_right h]

/-- **An address permutation conjugates `dir` into `dir`** along the mapped direction. -/
theorem permMat_conj_dir (G : APerm α) (z : Space α) :
    permMat G * dir z * permMat (apInv G) = dir ((apInv G).π z) := by
  unfold dir
  rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.mul_smul, Matrix.smul_mul, Matrix.smul_mul,
    Matrix.mul_one, permMat_mul_inv, permMat_conj_shift]

/-- conjugation by `P` (with `Q * P = 1 = P * Q`) is multiplicative on list products. -/
lemma conj_list_prod {P Q : CMat α} (hQP : Q * P = 1) (hPQ : P * Q = 1) (L : List (CMat α)) :
    P * L.prod * Q = (L.map fun M => P * M * Q).prod := by
  induction L with
  | nil => simp [hPQ]
  | cons a L ih =>
    rw [List.prod_cons, List.map_cons, List.prod_cons, ← ih]
    calc P * (a * L.prod) * Q = P * a * 1 * L.prod * Q := by simp only [mul_one, mul_assoc]
      _ = P * a * (Q * P) * L.prod * Q := by rw [hQP]
      _ = P * a * Q * (P * L.prod * Q) := by simp only [mul_assoc]

/-- **An address permutation conjugates a block into a block** along the mapped directions. -/
theorem permMat_conj_blockMat (G : APerm α) {r : ℕ} (z : Fin r → Space α) :
    permMat G * blockMat z * permMat (apInv G) = blockMat (fun i => (apInv G).π (z i)) := by
  unfold blockMat
  rw [conj_list_prod (permMat_inv_mul G) (permMat_mul_inv G), List.map_ofFn]
  congr 1
  apply congrArg List.ofFn
  funext i
  exact permMat_conj_dir G (z i)

/-- The same, read the other way: `permMat G⁻¹ * blockMat z * permMat G = blockMat (G ∘ z)`. -/
theorem permMat_inv_conj_blockMat (G : APerm α) {r : ℕ} (z : Fin r → Space α) :
    permMat (apInv G) * blockMat z * permMat G = blockMat (fun i => G.π (z i)) := by
  have h := permMat_conj_blockMat (apInv G) z
  rw [permMat_apInv_apInv] at h
  rw [h]
  congr 1

end Conj

section Transport
variable {α : Type*} [Fintype α] [DecidableEq α]

/-- The additive bijection of a linear automorphism of the label space. -/
def apOfLin (L : Space α ≃ₗ[F] Space α) : APerm α := ⟨L.toEquiv, fun x y => map_add L x y⟩

@[simp] lemma apOfLin_apply (L : Space α ≃ₗ[F] Space α) (x : Space α) : (apOfLin L).π x = L x := rfl

lemma sumExtend_inl {ι : Type*} {z : ι → Space α} (hz : LinearIndependent F z) (i : ι) :
    Basis.sumExtend hz (Sum.inl i) = z i := by
  simp only [Basis.sumExtend, Basis.reindex_apply, Basis.coe_extend]
  rfl

/-- **The transport σ.**  Any `r` linearly independent directions `z` are carried onto any `r` linearly independent
directions `u` by an additive bijection (a FREE address permutation). -/
theorem exists_transport {r : ℕ} (z u : Fin r → Space α) (hz : LinearIndependent F z)
    (hu : LinearIndependent F u) : ∃ σ : APerm α, ∀ i, σ.π (z i) = u i := by
  classical
  let bz := Basis.sumExtend hz
  let bu := Basis.sumExtend hu
  have cz := Module.finrank_eq_card_basis bz
  have cu := Module.finrank_eq_card_basis bu
  rw [Fintype.card_sum, Fintype.card_fin] at cz cu
  have hcard : Fintype.card (Basis.sumExtendIndex hz) = Fintype.card (Basis.sumExtendIndex hu) := by omega
  let e := Fintype.equivOfCardEq hcard
  refine ⟨apOfLin (bz.equiv bu (Equiv.sumCongr (Equiv.refl _) e)), fun i => ?_⟩
  rw [apOfLin_apply, ← sumExtend_inl hz i]
  show (bz.equiv bu (Equiv.sumCongr (Equiv.refl _) e)) (bz (Sum.inl i)) = u i
  rw [Basis.equiv_apply]
  simp only [Equiv.sumCongr_apply, Sum.map_inl, Equiv.refl_apply]
  exact sumExtend_inl hu i

/-- **The transport σ onto pending directions.**  The images `z` of the missing label coordinates of a slot are
carried onto `r` of the host's pending directions `v` (chosen by an injection `e`, the hosted set `S'`). -/
theorem exists_transport_into {r n : ℕ} (z : Fin r → Space α) (v : Fin n → Space α) (e : Fin r → Fin n)
    (he : Function.Injective e) (hz : LinearIndependent F z) (hv : LinearIndependent F v) :
    ∃ σ : APerm α, ∀ i, σ.π (z i) = v (e i) :=
  exists_transport z (fun i => v (e i)) hz (hv.comp e he)

/-- The transported block: `permMat σ⁻¹ * blockMat z * permMat σ = blockMat u` when `σ ∘ z = u`. -/
theorem transport_blockMat {r : ℕ} {z u : Fin r → Space α} {σ : APerm α} (hσ : ∀ i, σ.π (z i) = u i) :
    permMat (apInv σ) * blockMat z * permMat σ = blockMat u := by
  rw [permMat_inv_conj_blockMat, show (fun i => σ.π (z i)) = u from funext hσ]

/-- **Transport of a gauged residual.**  The residual `C * blockMat z * C'` of a slot (L1 form) is conjugated by free
adapters into the host's gauged factor `D * blockMat u * D'`. -/
theorem transport_twisted {r : ℕ} (z u : Fin r → Space α) (hz : LinearIndependent F z)
    (hu : LinearIndependent F u) {C C' D D' : CMat α} (hC : Free C) (hC' : Free C') (hCC : C' * C = 1)
    (hD : Free D) (hD' : Free D') (hDD : D' * D = 1) :
    ∃ T T' : CMat α, Free T ∧ Free T' ∧ T' * T = 1 ∧ T * T' = 1 ∧
      T * (C * blockMat z * C') * T' = D * blockMat u * D' := by
  obtain ⟨σ, hσ⟩ := exists_transport z u hz hu
  have hCC' : C * C' = 1 := mul_eq_one_comm.mp hCC
  have hDD' : D * D' = 1 := mul_eq_one_comm.mp hDD
  refine ⟨D * permMat (apInv σ) * C', C * permMat σ * D', (hD.mul (Free.perm _)).mul hC',
    (hC.mul (Free.perm _)).mul hD', ?_, ?_, ?_⟩
  · calc C * permMat σ * D' * (D * permMat (apInv σ) * C')
        = C * permMat σ * (D' * D) * permMat (apInv σ) * C' := by simp only [mul_assoc]
      _ = C * (permMat σ * permMat (apInv σ)) * C' := by rw [hDD]; simp only [mul_one, mul_assoc]
      _ = 1 := by rw [permMat_mul_inv, mul_one, hCC']
  · calc D * permMat (apInv σ) * C' * (C * permMat σ * D')
        = D * permMat (apInv σ) * (C' * C) * permMat σ * D' := by simp only [mul_assoc]
      _ = D * (permMat (apInv σ) * permMat σ) * D' := by rw [hCC]; simp only [mul_one, mul_assoc]
      _ = 1 := by rw [permMat_inv_mul, mul_one, hDD']
  · calc D * permMat (apInv σ) * C' * (C * blockMat z * C') * (C * permMat σ * D')
        = D * (permMat (apInv σ) * (C' * C) * blockMat z * (C' * C) * permMat σ) * D' := by
          simp only [mul_assoc]
      _ = D * (permMat (apInv σ) * blockMat z * permMat σ) * D' := by rw [hCC]; simp only [mul_one, mul_assoc]
      _ = D * blockMat u * D' := by rw [transport_blockMat hσ]

end Transport

section Hosted
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- **Hosted residuals stack.**  A host already carrying hosted factors `a` in its gauge receives a transported
residual on further directions `u`: the factors concatenate. -/
theorem host_stack {r k : ℕ} (a : Fin k → Space α) (u : Fin r → Space α) {D D' R : CMat α}
    (hDD : D' * D = 1) (hR : R = D * blockMat u * D') (M : CMat α) :
    R * (D * blockMat a * D' * M) = D * blockMat (Fin.append a u) * D' * M := by
  subst hR
  calc D * blockMat u * D' * (D * blockMat a * D' * M)
      = D * blockMat u * (D' * D) * blockMat a * D' * M := by simp only [mul_assoc]
    _ = D * (blockMat u * blockMat a) * D' * M := by rw [hDD]; simp only [mul_one, mul_assoc]
    _ = D * blockMat (Fin.append a u) * D' * M := by rw [blockMat_comm u a, ← blockMat_append]

/-- **The hosted climb.**  A host at `M` whose idle climb is ONE gauged block `D * blockMat (Fin.append u y) * D'`
hosts the gauged residual `C * blockMat z * C'` of a slot (with `card z = card u`): free `T'` before the
invocation, the invocation's residual on the slot input `T' * M`, then (free `T` and) ONE block of rank `s`. -/
theorem hosted_climb (l : ρ) (M : CMat α) {r s : ℕ} (z u : Fin r → Space α) (y : Fin s → Space α)
    (hz : LinearIndependent F z) (hv : LinearIndependent F (Fin.append u y))
    (h1 : 1 ≤ s) (h2 : s < Fintype.card α)
    {C C' D D' : CMat α} (hC : Free C) (hC' : Free C') (hCC : C' * C = 1)
    (hD : Free D) (hD' : Free D') (hDD : D' * D = 1) :
    ∃ T T' : CMat α, Free T ∧ Free T' ∧ T' * T = 1 ∧
      GStep l M (T' * M) (fun _ => 0) ∧
      GStep l (C * blockMat z * C' * (T' * M)) (D * blockMat (Fin.append u y) * D' * M) (fun φ => φ s) := by
  have hu : LinearIndependent F u := by
    have h := hv.comp (Fin.castAdd s) (Fin.castAdd_injective r s)
    have e : (Fin.append u y ∘ Fin.castAdd s) = u := by funext i; simp
    rwa [e] at h
  have hy : LinearIndependent F y := by
    have h := hv.comp (Fin.natAdd r) (Fin.natAdd_injective s r)
    have e : (Fin.append u y ∘ Fin.natAdd r) = y := by funext i; simp
    rwa [e] at h
  obtain ⟨T, T', hT, hT', hTT, -, hconj⟩ := transport_twisted z u hz hu hC hC' hCC hD hD' hDD
  refine ⟨T, T', hT, hT', hTT, gstep_free hT' l M, ?_⟩
  have a := gstep_free hT l (C * blockMat z * C' * (T' * M))
  have b := host_remaining_climb_twisted l M hD hD' hDD u y hy h1 h2
  have e : T * (C * blockMat z * C' * (T' * M)) = D * blockMat u * D' * M := by
    rw [← hconj]; simp only [mul_assoc]
  rw [e] at a
  exact (a.trans b).cast (fun φ => by simp)

/-- **The hosted climb, fully hosted** (`s = 0`): the transported residual IS the host's idle climb, which then costs
nothing beyond the invocation. -/
theorem hosted_climb_full (l : ρ) (M : CMat α) {r : ℕ} (z u : Fin r → Space α)
    (hz : LinearIndependent F z) (hu : LinearIndependent F u)
    {C C' D D' : CMat α} (hC : Free C) (hC' : Free C') (hCC : C' * C = 1)
    (hD : Free D) (hD' : Free D') (hDD : D' * D = 1) :
    ∃ T T' : CMat α, Free T ∧ Free T' ∧ T' * T = 1 ∧
      GStep l M (T' * M) (fun _ => 0) ∧
      GStep l (C * blockMat z * C' * (T' * M)) (D * blockMat u * D' * M) (fun _ => 0) := by
  obtain ⟨T, T', hT, hT', hTT, -, hconj⟩ := transport_twisted z u hz hu hC hC' hCC hD hD' hDD
  refine ⟨T, T', hT, hT', hTT, gstep_free hT' l M, ?_⟩
  have a := gstep_free hT l (C * blockMat z * C' * (T' * M))
  have e : T * (C * blockMat z * C' * (T' * M)) = D * blockMat u * D' * M := by
    rw [← hconj]; simp only [mul_assoc]
  rw [e] at a
  exact a

/-- A block does not depend on the order of its directions. -/
theorem blockMat_comp_equiv {m n : ℕ} (v : Fin n → Space α) (e : Fin m ≃ Fin n) :
    blockMat (fun i => v (e i)) = blockMat v := by
  rw [blockMat_wrap, blockMat_wrap]
  congr 1
  funext x
  rw [List.map_ofFn, List.map_ofFn, List.prod_ofFn, List.prod_ofFn]
  exact e.prod_comp (fun j => tint (dot (v j) x))

lemma append_split {r s n : ℕ} (v : Fin n → Space α) (e : Fin (r + s) ≃ Fin n) :
    Fin.append (fun i => v (e (Fin.castAdd s i))) (fun j => v (e (Fin.natAdd r j))) = fun k => v (e k) := by
  funext k
  refine Fin.addCases (fun i => ?_) (fun j => ?_) k <;> simp

/-- A block along `v` splits along any partition `e` of its directions into a hosted part and a remainder. -/
theorem blockMat_split {r s n : ℕ} (v : Fin n → Space α) (e : Fin (r + s) ≃ Fin n) :
    blockMat v =
      blockMat (Fin.append (fun i => v (e (Fin.castAdd s i))) (fun j => v (e (Fin.natAdd r j)))) := by
  rw [append_split, blockMat_comp_equiv]

/-- **The hosted climb on the host's own directions.**  The host's idle climb is ONE gauged block along its pending
directions `v`; a partition `e` puts `r` of them (the hosted set `S'`) first.  The slot residual is transported onto
them and the remainder is ONE block of rank `s`. -/
theorem hosted_climb_of_split (l : ρ) (M : CMat α) {r s n : ℕ} (z : Fin r → Space α) (v : Fin n → Space α)
    (e : Fin (r + s) ≃ Fin n) (hz : LinearIndependent F z) (hv : LinearIndependent F v)
    (h1 : 1 ≤ s) (h2 : s < Fintype.card α)
    {C C' D D' : CMat α} (hC : Free C) (hC' : Free C') (hCC : C' * C = 1)
    (hD : Free D) (hD' : Free D') (hDD : D' * D = 1) :
    ∃ T T' : CMat α, Free T ∧ Free T' ∧ T' * T = 1 ∧
      GStep l M (T' * M) (fun _ => 0) ∧
      GStep l (C * blockMat z * C' * (T' * M)) (D * blockMat v * D' * M) (fun φ => φ s) := by
  have hv' : LinearIndependent F
      (Fin.append (fun i => v (e (Fin.castAdd s i))) (fun j => v (e (Fin.natAdd r j)))) := by
    rw [append_split]
    exact hv.comp e e.injective
  rw [blockMat_split v e]
  exact hosted_climb l M z _ _ hz hv' h1 h2 hC hC' hCC hD hD' hDD

end Hosted
end
end XHost
end PowerSaving.GF
end OAI

#print axioms OAI.PowerSaving.GF.XHost.permMat_conj_dir
#print axioms OAI.PowerSaving.GF.XHost.permMat_conj_blockMat
#print axioms OAI.PowerSaving.GF.XHost.exists_transport
#print axioms OAI.PowerSaving.GF.XHost.exists_transport_into
#print axioms OAI.PowerSaving.GF.XHost.transport_twisted
#print axioms OAI.PowerSaving.GF.XHost.host_stack
#print axioms OAI.PowerSaving.GF.XHost.hosted_climb
#print axioms OAI.PowerSaving.GF.XHost.hosted_climb_full
#print axioms OAI.PowerSaving.GF.XHost.hosted_climb_of_split
