import Work.KMX.GlueXZ
import CSX.ScalarZX

/-!
# fsp-cc (14 Oct 2026): the y-free twin route with split centres on SOURCES

k-lower's y-free twin (`GXD.nyRaw`, Work/KLow/NoY.lean) keeps the non-y rows of the total matrix (`nyRaw_Mt`); the
zero checks run on the twin.  With centres on sources the total matrix is `GS.MtX` (Work/CSX/ScalarZX.lean).  Here:
* `MtX_nonY`: the non-y rows of `MtX c` are those of `Mt c` (the source part of the scatter, like the slot part,
  writes y roles only, and no add of phase B with a non-y target reads a y role);
* `exit_of_twinZX`: the three exit conditions on `MtX c` from those on `MtX (nyRaw c)` (zero checks `zCheckZX` on the
  twin), through `MtX_nonY` (twice) and k-lower's `nyRaw_Mt`.
No `sorry`.  Credits: k-lower (nyRaw, nyRaw_Mt, diag_matP), k-combine (exit_of_twinZ), Jacob Sussman (engine).
-/

set_option linter.unusedVariables false
set_option linter.unusedSectionVars false
set_option linter.deprecated false

namespace GXD
open SSC OAI.PowerSaving OAI.PowerSaving.CR OAI.PowerSaving.CB Matrix Finset

/-- **the non-y rows of the total matrix ignore the source part of the scatter** -/
theorem MtX_nonY (c : Raw) (hv : 0 < c.v)
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (i j : GS.Rl c) (hi : iY i = 0) :
    GS.MtX c hv i j = GS.Mt c hv i j := by
  have hd : ∀ r, isYN c.v r = true → dN (GS.unemb c hv r) = 0 := fun r => (isYN_unemb c hv r).1
  have hd1 : ∀ r, isYN c.v r = false → dN (GS.unemb c hv r) = 1 := fun r => (isYN_unemb c hv r).2
  have hBy : ∀ x ∈ GS.addsB c, isYN c.v x.2.1 = true → isYN c.v x.1 = true := by
    intro x hx h
    simp only [isYN, Bool.and_eq_true, Nat.ble_eq, Nat.blt_eq] at h ⊢
    exact hB x hx h.1 h.2
  have key := diag_matP dN (GS.unemb c hv) (isYN c.v) hd hd1 (GS.addsB c) hBy
  have hSy : ∀ K : Matrix (Fin c.v) (Fin c.R) ℚ,
      Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * scatM K
        = Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) := by
    intro K
    unfold scatM
    rw [Matrix.mul_add, Matrix.mul_one, add_eq_left]
    ext a b
    rw [Matrix.diagonal_mul]
    rcases a with ((t|t)|q)
    · simp [embYS]
    · simp [dN, iY]
    · simp [embYS]
  have hSX : ∀ K : Matrix (Fin c.v) (Fin c.v) ℚ,
      Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * CSX.scatX (Sl := Fin c.R) K
        = Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) := by
    intro K
    unfold CSX.scatX
    rw [Matrix.mul_add, Matrix.mul_one, add_eq_left]
    ext a b
    rw [Matrix.diagonal_mul]
    rcases a with ((t|t)|q)
    · simp [CSX.embYX]
    · simp [dN, iY]
    · simp [CSX.embYX]
  have eX : Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * GS.MtX c hv
      = matP (GS.unemb c hv) (GS.mic ((GS.addsB c).filter (fun x => !isYN c.v x.1)))
        * (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * GS.MA c hv) := by
    show Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * (matP (GS.unemb c hv) (GS.mic (GS.addsB c))
      * scatM (GS.Jrm c * GS.Ccm c) * (CSX.scatX (Sl := Fin c.R) (GS.Jrm c * GS.Cxm c) * GS.MA c hv)) = _
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R))) (matP _ _), key]
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R))) (scatM _), hSy,
      ← Matrix.mul_assoc (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R))) (CSX.scatX _), hSX]
  have e0 : Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * GS.Mt c hv
      = matP (GS.unemb c hv) (GS.mic ((GS.addsB c).filter (fun x => !isYN c.v x.1)))
        * (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * GS.MA c hv) := by
    show Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * (matP (GS.unemb c hv) (GS.mic (GS.addsB c))
      * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv) = _
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R))) (matP _ _), key]
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R))) (scatM _), hSy]
  have h1 := congrFun (congrFun (eX.trans e0.symm) i) j
  simp only [Matrix.diagonal_mul] at h1
  have hdi : (dN (T := Fin c.v) (Sl := Fin c.R)) i = 1 := by simp [dN, hi]
  rw [hdi, one_mul, one_mul] at h1
  exact h1

/-- the non-y rows of `MtX` are those of the y-free twin's `MtX` -/
theorem MtX_twin (c : Raw) (hv : 0 < c.v)
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (i j : GS.Rl c) (hi : iY i = 0) : GS.MtX c hv i j = GS.MtX (nyRaw c) hv i j := by
  have hBn : ∀ x ∈ GS.addsB (nyRaw c), c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v := by
    intro x hx
    rw [addsB_ny] at hx
    exact hB x (List.mem_filter.mp hx).1
  rw [MtX_nonY c hv hB i j hi, nyRaw_Mt c hv hB i j hi]
  exact (MtX_nonY (nyRaw c) hv hBn i j hi).symm

/-- the y columns of the twin's `MtX` at the slot rows (from `hy` of `MtX c`) -/
theorem twin_hyX (c : Raw) (hv : 0 < c.v)
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (hy : ∀ (i : GS.Rl c) (t : Fin c.v), GS.MtX c hv i (lY t) = if i = lY t then 1 else 0) :
    ∀ (q : Fin c.R) (t : Fin c.v), GS.MtX (nyRaw c) hv (lS q) (lY t) = 0 := by
  intro q t
  have h0 := hy (lS q) t
  rw [if_neg (by simp [lS, lY])] at h0
  exact (MtX_twin c hv hB (lS q) (lY t) (by simp [iY, lS])).symm.trans h0

/-- **the three exit conditions on `MtX c` from those on the y-free twin** (zero checks `zCheckZX` on `nyRaw c`) -/
theorem exit_of_twinZX (c : Raw) (hv : 0 < c.v) (el : List Nat)
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (h : XEX.ExitRows (GLab.exitSet (nyRaw c) el) (GS.MtX (nyRaw c) hv) ∧
      XEX.ExitCols (GLab.exitSet (nyRaw c) el) (GS.MtX (nyRaw c) hv) ∧
      XEX.ExitColsX (GLab.exitSet (nyRaw c) el) (GS.MtX (nyRaw c) hv)) :
    XEX.ExitRows (GLab.exitSet c el) (GS.MtX c hv) ∧ XEX.ExitCols (GLab.exitSet c el) (GS.MtX c hv) ∧
      XEX.ExitColsX (GLab.exitSet c el) (GS.MtX c hv) := by
  have hBn : ∀ x ∈ GS.addsB (nyRaw c), c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v := by
    intro x hx
    rw [addsB_ny] at hx
    exact hB x (List.mem_filter.mp hx).1
  have tr : ∀ (i j : GS.Rl c), iY i = 0 → GS.MtX c hv i j = GS.MtX (nyRaw c) hv i j := by
    intro i j hi
    rw [MtX_nonY c hv hB i j hi, nyRaw_Mt c hv hB i j hi]
    exact (MtX_nonY (nyRaw c) hv hBn i j hi).symm
  refine ⟨fun e he j hj => ?_, fun e he q hq => ?_, fun e he t => ?_⟩
  · rw [tr (lS e) j (by simp [iY, lS])]
    exact h.1 e he j hj
  · rw [tr (lS q) (lS e) (by simp [iY, lS])]
    exact h.2.1 e he q hq
  · rw [tr (lX t) (lS e) (by simp [iY, lX])]
    exact h.2.2 e he t

end GXD

#print axioms GXD.MtX_nonY
#print axioms GXD.exit_of_twinZX
#print axioms GXD.twin_hyX
