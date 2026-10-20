import Work.KMX.HalvesZ

/-!
# fsp-cc (14 Oct 2026): split centres on MUTABLE SOURCE registers: the scatter also reads x roles

Sussman's scatter (`Work/Carrier/{Scat,Glue}.lean`) reads the retained totals from SLOTS only:
copies `c += Cc s` (`Cc : Matrix C Sl`), the copies come down, `y += Jr c`; on the live roles this is
`scatM (Jr * Cc)` (`y += (Jr * Cc) s`).  With x-mutable's mutable sources a split centre may end phase A in a SOURCE
register; here a second copy matrix `Cx : Matrix C T` lets copy `k` also read x roles:

    c += Cc s ;  c += Cx x ;  (copies down) ;  y += Jr c                       (forward)
    c += Jrᵀ y ;  (copies up) ;  s -= Ccᵀ c ;  x -= Cxᵀ c                      (backward)

* `embYX K`, `scatX K = 1 + embYX K` (`y += K x`), `scatXi K = 1 - embYX K`, `scatX_inv`, `act_scatX`, `act_scatXT`;
* `scat2_glue` / `scat2_glueT`: on zero copies the gate sequences above act on the live roles as
  `scatM (Jr * Cc) ∘ scatX (Jr * Cx)` and `(scatMi (Jr * Cc))ᵀ ∘ (scatXi (Jr * Cx))ᵀ`;
* `fwd_scalarX` / `bwd_scalarX`: for x-mutable's scalar half `K : ScalZ` whose phase-A product is
  `K.MA = scatX (K.Jr * Cx) * MA0` (and `K.MAi = MA0i * scatXi (K.Jr * Cx)`), the real word (phase A = `MA0`, the
  scatter reading slots AND sources) is `y += x` (resp. `x -= y`) with sources and slots restored, i.e. exactly
  `ScalZ.fwd_scalar` / `ScalZ.bwd_scalar`.  So every scalar fact of `K` (hy, hid, E7, the exit conditions) is a fact
  about the true total matrix `MB * scatM (Jr Cc) * scatX (Jr Cx) * MA0`, and nothing downstream changes.

Credits: Jacob Sussman (engine, scatter, dirt gates; Apache-2.0); x-mutable (mutable sources, `ScalZ`, `g8z`);
Zhihao Chen (#425, the mutable-source interface).  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace CSX
open Binary Matrix Finset RAM SS CB BR CR GF RF GX B401 XEX
noncomputable section

section Scat
variable {T Sl : Type} [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl]

/-- the block `y ← x` of a square matrix `K : T × T`, as a matrix on the live roles. -/
def embYX (K : Matrix T T ℚ) : Matrix (L3 T Sl) (L3 T Sl) ℚ := fun i j =>
  Sum.elim (Sum.elim (fun _ => 0)
    (fun S => Sum.elim (Sum.elim (fun t => K S t) (fun _ => 0)) (fun _ => 0) j)) (fun _ => 0) i

/-- the source part of the scatter: `y += K x`. -/
def scatX (K : Matrix T T ℚ) : Matrix (L3 T Sl) (L3 T Sl) ℚ := 1 + embYX K
/-- its inverse: `y -= K x`. -/
def scatXi (K : Matrix T T ℚ) : Matrix (L3 T Sl) (L3 T Sl) ℚ := 1 - embYX K

lemma embYX_mul (A B : Matrix T T ℚ) : (embYX A : Matrix (L3 T Sl) (L3 T Sl) ℚ) * embYX B = 0 := by
  ext i j
  rw [Matrix.mul_apply, Matrix.zero_apply]
  apply Finset.sum_eq_zero
  intro k _
  rcases k with ((t|S)|q)
  · simp [embYX]
  · rcases i with ((t'|S')|q') <;> simp [embYX]
  · simp [embYX]

lemma scatX_inv (K : Matrix T T ℚ) : (scatXi K : Matrix (L3 T Sl) (L3 T Sl) ℚ) * scatX K = 1 := by
  unfold scatXi scatX
  rw [sub_mul, one_mul, mul_add, mul_one, embYX_mul, add_zero]
  abel

/-- scalar map of the source scatter: `y_S += ∑_t K S t * x_t`. -/
lemma act_scatX (K : Matrix T T ℚ) (u : L3 T Sl → ℂ) :
    actPoint (scatX K) u = fun i => u i +
      Sum.elim (Sum.elim (fun _ => 0) (fun S => ap K (fun t => u (Sum.inl (Sum.inl t))) S)) (fun _ => 0) i := by
  have h : actPoint (scatX K) u = u + ap (embYX K) u := by
    have e := ap_madd (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) (embYX K) u
    rw [ap_one] at e
    exact e
  funext i
  have e2 : actPoint (scatX K) u i = u i + ap (embYX K) u i := congrFun h i
  rw [e2]
  congr 1
  rcases i with ((t|S)|q)
  · show ∑ j, ((embYX K (Sum.inl (Sum.inl t)) j : ℚ) : ℂ) * u j = 0
    apply Finset.sum_eq_zero
    intro j _
    simp [embYX]
  · show ∑ j, ((embYX K (Sum.inl (Sum.inr S)) j : ℚ) : ℂ) * u j
      = ∑ t, ((K S t : ℚ) : ℂ) * u (Sum.inl (Sum.inl t))
    rw [Fintype.sum_sum_type, Fintype.sum_sum_type]
    simp [embYX]
  · show ∑ j, ((embYX K (Sum.inr q) j : ℚ) : ℂ) * u j = 0
    apply Finset.sum_eq_zero
    intro j _
    simp [embYX]

/-- scalar map of the transposed inverse source scatter: `x_t -= ∑_S K S t * y_S`. -/
lemma act_scatXT (K : Matrix T T ℚ) (u : L3 T Sl → ℂ) :
    actPoint ((scatXi K)ᵀ) u = fun i => u i -
      Sum.elim (Sum.elim (fun t => ap Kᵀ (fun S => u (Sum.inl (Sum.inr S))) t) (fun _ => 0)) (fun _ => 0) i := by
  have e' : (scatXi K : Matrix (L3 T Sl) (L3 T Sl) ℚ)ᵀ = 1 + -(embYX K)ᵀ := by
    unfold scatXi
    rw [Matrix.transpose_sub, Matrix.transpose_one, sub_eq_add_neg]
  have h : actPoint (1 + -(embYX K : Matrix (L3 T Sl) (L3 T Sl) ℚ)ᵀ) u = u + -(ap (embYX K)ᵀ u) := by
    have e := ap_madd (1 : Matrix (L3 T Sl) (L3 T Sl) ℚ) (-(embYX K)ᵀ) u
    rw [ap_one, ap_mneg] at e
    exact e
  funext i
  have e2 : actPoint (1 + -(embYX K : Matrix (L3 T Sl) (L3 T Sl) ℚ)ᵀ) u i = u i + -(ap (embYX K)ᵀ u i) :=
    congrFun h i
  rw [e', e2, ← sub_eq_add_neg]
  congr 1
  rcases i with ((t|S)|q)
  · show ∑ j, (((embYX K)ᵀ (Sum.inl (Sum.inl t)) j : ℚ) : ℂ) * u j
      = ∑ S, ((Kᵀ t S : ℚ) : ℂ) * u (Sum.inl (Sum.inr S))
    rw [Fintype.sum_sum_type, Fintype.sum_sum_type]
    simp [embYX, Matrix.transpose_apply]
  · show ∑ j, (((embYX K)ᵀ (Sum.inl (Sum.inr S)) j : ℚ) : ℂ) * u j = 0
    apply Finset.sum_eq_zero
    intro j _
    rcases j with ((t'|S')|q') <;> simp [embYX, Matrix.transpose_apply]
  · show ∑ j, (((embYX K)ᵀ (Sum.inr q) j : ℚ) : ℂ) * u j = 0
    apply Finset.sum_eq_zero
    intro j _
    rcases j with ((t'|S')|q') <;> simp [embYX, Matrix.transpose_apply]

/-- product with the source scatter: `(scatX K * X) i j = X i j + [i = y_S] ∑_t K S t * X (x_t) j`. -/
lemma scatX_mul_apply {m : Type} (K : Matrix T T ℚ) (X : Matrix (L3 T Sl) m ℚ) (i : L3 T Sl) (j : m) :
    (scatX (Sl := Sl) K * X) i j = X i j +
      Sum.elim (Sum.elim (fun _ => 0) (fun S => ∑ t, K S t * X (Sum.inl (Sum.inl t)) j)) (fun _ => 0) i := by
  unfold scatX
  rw [Matrix.add_mul, Matrix.one_mul, Matrix.add_apply, Matrix.mul_apply, Fintype.sum_sum_type,
    Fintype.sum_sum_type]
  rcases i with (a | S) | q
  · simp [embYX]
  · simp [embYX]
  · simp [embYX]

end Scat

section Glue
variable {T Sl C : Type} [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl]
  [Fintype C] [DecidableEq C]

lemma liveAct_glueX (φ : (L3 T Sl → ℂ) → (L3 T Sl → ℂ)) (u : L3 T Sl → ℂ) (c : C → ℂ) :
    liveAct φ (glue u c) = glue (φ u) c :=
  congrArg (fun v => glue (φ v) c) (glue_up u c)

/-- **on zero copies, `c += Cc s ; c += Cx x ; y += Jr c` is `scatM (Jr * Cc) ∘ scatX (Jr * Cx)` on the live
roles.** -/
lemma scat2_glue (Jr : Matrix T C ℚ) (Cc : Matrix C Sl ℚ) (Cx : Matrix C T ℚ) (z : Box T Sl C → ℂ)
    (hc : z ∘ bC = 0) :
    gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) ∘ (upE : L3 T Sl ↪ Box T Sl C)
      = actPoint (scatM (Jr * Cc)) (actPoint (scatX (Jr * Cx)) (z ∘ upE)) := by
  have hX : gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) ∘ bX = z ∘ bX := by
    simp only [gateS, aY_X, aC_X]
  have hY : gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) ∘ bY
      = z ∘ bY + ap Jr (ap Cc (z ∘ bS) + ap Cx (z ∘ bX)) := by
    simp only [gateS, aY_Y, aC_Y, aC_C, aC_S, aC_X, hc, zero_add]
  have hS : gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) ∘ bS = z ∘ bS := by
    simp only [gateS, aY_S, aC_S]
  rw [act_scatX]
  unfold scatM
  rw [act_scat]
  funext r
  rcases r with ((t|S)|q)
  · show gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) (bX t) = (z (bX t) + 0) + 0
    rw [add_zero, add_zero]
    exact congrFun hX t
  · have h : gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) (bY S)
        = z (bY S) + ap Jr (ap Cc (z ∘ bS) + ap Cx (z ∘ bX)) S := congrFun hY S
    show gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) (bY S)
      = (z (bY S) + ap (Jr * Cx) (fun t => z (bX t)) S) + ap (Jr * Cc) (fun q => z (bS q) + 0) S
    simp only [add_zero]
    rw [h, ap_add, Pi.add_apply]
    show z (bY S) + (ap Jr (ap Cc (z ∘ bS)) S + ap Jr (ap Cx (z ∘ bX)) S)
      = (z (bY S) + ap (Jr * Cx) (z ∘ bX) S) + ap (Jr * Cc) (z ∘ bS) S
    rw [ap_mul, ap_mul]
    abel
  · show gateS bY bC Jr (gateS bC bX Cx (gateS bC bS Cc z)) (bS q) = (z (bS q) + 0) + 0
    rw [add_zero, add_zero]
    exact congrFun hS q

/-- **on zero copies, `c += Jrᵀ y ; s -= Ccᵀ c ; x -= Cxᵀ c` is `(scatMi (Jr Cc))ᵀ ∘ (scatXi (Jr Cx))ᵀ` on the live
roles.** -/
lemma scat2_glueT (Jr : Matrix T C ℚ) (Cc : Matrix C Sl ℚ) (Cx : Matrix C T ℚ) (z : Box T Sl C → ℂ)
    (hc : z ∘ bC = 0) :
    gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) ∘ (upE : L3 T Sl ↪ Box T Sl C)
      = actPoint ((scatMi (Jr * Cc))ᵀ) (actPoint ((scatXi (Jr * Cx))ᵀ) (z ∘ upE)) := by
  have hX : gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) ∘ bX
      = z ∘ bX + -(ap (Jr * Cx)ᵀ (z ∘ bY)) := by
    simp only [gateS, aX_X, aS_X, aC_X, aX_C, aS_C, aC_C, aC_Y, hc, zero_add, ap_mneg, Matrix.transpose_mul,
      ap_mul]
  have hY : gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) ∘ bY = z ∘ bY := by
    simp only [gateS, aX_Y, aS_Y, aC_Y]
  have hS : gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) ∘ bS
      = z ∘ bS + -(ap (Jr * Cc)ᵀ (z ∘ bY)) := by
    simp only [gateS, aX_S, aS_S, aC_S, aC_C, aC_Y, hc, zero_add, ap_mneg, Matrix.transpose_mul, ap_mul]
  rw [act_scatXT]
  unfold scatMi
  rw [act_scatT]
  funext r
  rcases r with ((t|S)|q)
  · have h : gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) (bX t)
        = z (bX t) + -(ap (Jr * Cx)ᵀ (z ∘ bY) t) := congrFun hX t
    show gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) (bX t)
      = (z (bX t) - ap (Jr * Cx)ᵀ (fun S => z (bY S)) t) - 0
    rw [sub_zero, h, sub_eq_add_neg]
    try rfl
  · show gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) (bY S) = (z (bY S) - 0) - 0
    rw [sub_zero, sub_zero]
    exact congrFun hY S
  · have h : gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) (bS q)
        = z (bS q) + -(ap (Jr * Cc)ᵀ (z ∘ bY) q) := congrFun hS q
    show gateS bX bC (-Cxᵀ) (gateS bS bC (-Ccᵀ) (gateS bC bY Jrᵀ z)) (bS q)
      = (z (bS q) - 0) - ap (Jr * Cc)ᵀ (fun S => z (bY S) - 0) q
    simp only [sub_zero]
    rw [h, sub_eq_add_neg]
    try rfl

end Glue

section Word
variable {T Sl C : Type} [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl]
  [Fintype C] [DecidableEq C]
variable (K : ScalZ T Sl C) (Cx : Matrix C T ℚ) (MA0 MA0i : Matrix (L3 T Sl) (L3 T Sl) ℚ)

/-- the real live map before the scatter: `y -= K s`, then the TRUE phase A (`MA0`). -/
def f1x (u : L3 T Sl → ℂ) : L3 T Sl → ℂ := actPoint MA0 (actPoint (g1 K.Mt) u)
/-- the same for the backward word. -/
def b1x (u : L3 T Sl → ℂ) : L3 T Sl → ℂ := actPoint MA0iᵀ (actPoint (g1i K.Mt)ᵀ u)

/-- **Scalar identity of the forward word, centres on sources**: phase A, then `c += Cc s ; c += Cx x`, copies down,
`y += Jr c`, phase B, x and slot rows of the inverse: `y += x`, sources and slots restored. -/
theorem fwd_scalarX (hMA : K.MA = scatX (K.Jr * Cx) * MA0) (z : Box T Sl C → ℂ) :
    eraseC (liveAct K.f2 (gateS bY bC K.Jr (gateS bC bX Cx (gateS bC bS K.Cc
      (liveAct (f1x K MA0) (eraseC z)))))) = shearF z := by
  rw [← K.fwd_scalar z]
  have e1 : liveAct (f1x K MA0) (eraseC z) = glue (f1x K MA0 (z ∘ upE)) (fun _ => 0) := by
    rw [erase_glue z]; exact liveAct_glueX _ _ _
  have e2 : liveAct K.f1 (eraseC z) = glue (K.f1 (z ∘ upE)) (fun _ => 0) := by
    rw [erase_glue z]; exact liveAct_glueX _ _ _
  have hw : gateS bY bC K.Jr (gateS bC bX Cx (gateS bC bS K.Cc (glue (f1x K MA0 (z ∘ upE)) (fun _ => 0))))
        ∘ (upE : L3 T Sl ↪ Box T Sl C)
      = gateS bY bC K.Jr (gateS bC bS K.Cc (glue (K.f1 (z ∘ upE)) (fun _ => 0))) ∘ upE := by
    have h1 := scat2_glue K.Jr K.Cc Cx (glue (f1x K MA0 (z ∘ upE)) (fun _ => 0)) rfl
    have h2 := scat_glue K.Jr K.Cc (glue (K.f1 (z ∘ upE)) (fun _ => 0)) rfl
    rw [glue_up] at h1 h2
    rw [h1, h2]
    congr 1
    show actPoint (scatX (K.Jr * Cx)) (actPoint MA0 (actPoint (g1 K.Mt) (z ∘ upE)))
      = actPoint K.MA (actPoint (g1 K.Mt) (z ∘ upE))
    rw [hMA, CR.actPoint_mul]
  rw [e1, e2]
  unfold liveAct
  rw [erase_of_glue, erase_of_glue, hw]

/-- **Scalar identity of the backward word, centres on sources**: `x -= y`, sources and slots restored. -/
theorem bwd_scalarX (hMAi : K.MAi = MA0i * scatXi (K.Jr * Cx)) (z : Box T Sl C → ℂ) :
    eraseC (liveAct K.b2 (gateS bX bC (-Cxᵀ) (gateS bS bC (-K.Ccᵀ) (gateS bC bY K.Jrᵀ
      (liveAct (b1x K MA0i) (eraseC z)))))) = shearB z := by
  rw [← K.bwd_scalar z]
  have e1 : liveAct (b1x K MA0i) (eraseC z) = glue (b1x K MA0i (z ∘ upE)) (fun _ => 0) := by
    rw [erase_glue z]; exact liveAct_glueX _ _ _
  have e2 : liveAct K.b1 (eraseC z) = glue (K.b1 (z ∘ upE)) (fun _ => 0) := by
    rw [erase_glue z]; exact liveAct_glueX _ _ _
  have hw : gateS bX bC (-Cxᵀ) (gateS bS bC (-K.Ccᵀ) (gateS bC bY K.Jrᵀ
        (glue (b1x K MA0i (z ∘ upE)) (fun _ => 0)))) ∘ (upE : L3 T Sl ↪ Box T Sl C)
      = gateS bS bC (-K.Ccᵀ) (gateS bC bY K.Jrᵀ (glue (K.b1 (z ∘ upE)) (fun _ => 0))) ∘ upE := by
    have h1 := scat2_glueT K.Jr K.Cc Cx (glue (b1x K MA0i (z ∘ upE)) (fun _ => 0)) rfl
    have h2 := scat_glueT K.Jr K.Cc (glue (K.b1 (z ∘ upE)) (fun _ => 0)) rfl
    rw [glue_up] at h1 h2
    rw [h1, h2]
    congr 1
    show actPoint (scatXi (K.Jr * Cx))ᵀ (actPoint MA0iᵀ (actPoint (g1i K.Mt)ᵀ (z ∘ upE)))
      = actPoint K.MAiᵀ (actPoint (g1i K.Mt)ᵀ (z ∘ upE))
    rw [hMAi, Matrix.transpose_mul, CR.actPoint_mul]
  rw [e1, e2]
  unfold liveAct
  rw [erase_of_glue, erase_of_glue, hw]

end Word
end
end CSX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.CSX.scatX_inv
#print axioms OAI.PowerSaving.CSX.scat2_glue
#print axioms OAI.PowerSaving.CSX.scat2_glueT
#print axioms OAI.PowerSaving.CSX.fwd_scalarX
#print axioms OAI.PowerSaving.CSX.bwd_scalarX
