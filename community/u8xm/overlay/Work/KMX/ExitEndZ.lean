import Work.KLow.ExitEnd
import Work.XMut.GInv

/-!
# (key: k-combine) With MUTABLE sources, x-mutable's end gate `g8z Mti` never touches an exit helper

k-combine, 13 Oct 2026.  k-lower's `exit_hGr` / `exit_hHr` (Work/KLow/ExitEnd.lean) for x-mutable's derived end
gates `XMut.g8z Mti` (x AND slot rows of the inverse) and `(XMut.h8z Mt)ᵀ`.  With sources mutable the x rows of `Mt`
are no longer unit rows (no `hx`), so a THIRD finite condition replaces `hx`:

* `ExitColsX E Mt`  every exit column of `Mt` is zero at every x row (no source reads an exit helper; this is the
  matrix form of x-mutable's condition (d): with `ExitRows` and `ExitCols` it is equivalent to "no x row of `Mti`
  reads an exit helper", the form k-word checked on the U8 word).

Then the rest-gate fields `hGr` / `hHr` of `GLabXZ` hold for `Gr = g8z Mti`, `Hr = (h8z Mt)ᵀ` (`Gt = Ht = 1`).
Credits: x-mutable (`g8z`, `h8z`, `offdiag_z`), k-lower (`ExitRows`, `ExitCols`, `exit_row_inv`), x-exit-lean
(`restOk`), Jacob Sussman (the dirt rule).  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace XEX
open Binary Matrix Finset RAM SS CB BR CR GF RF GX B401
noncomputable section

section ExitEndZ
variable {T Sl : Type} [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl]

/-- every exit column of `M` is zero at every x row -/
def ExitColsX (E : Finset Sl) (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) : Prop :=
  ∀ e ∈ E, ∀ t : T, M (lX t) (lS e) = 0

/-- off the diagonal, `g8z` / `h8z` keep the non-y rows of their argument -/
lemma g8z_off (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) {i j : L3 T Sl} (e : i ≠ j) :
    XMut.g8z M i j = (1 - iY i) * M i j := by
  simp only [XMut.g8z, PY, Matrix.add_apply, Matrix.sub_mul, Matrix.one_mul, Matrix.sub_apply,
    Matrix.diagonal_mul, Matrix.diagonal_apply_ne _ e]
  ring

lemma h8z_off (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) {i j : L3 T Sl} (e : i ≠ j) :
    XMut.h8z M i j = (1 - iY i) * M i j := by
  simp only [XMut.h8z, PY, Matrix.add_apply, Matrix.sub_mul, Matrix.one_mul, Matrix.sub_apply,
    Matrix.diagonal_mul, Matrix.diagonal_apply_ne _ e]
  ring

/-- unit y columns: a non-y row has no y entry. -/
lemma col_y_zero' {M : Matrix (L3 T Sl) (L3 T Sl) ℚ} (h : M * PY = PY) {i : L3 T Sl} (hi : iY i = 0) (t : T) :
    M i (lY t) = 0 := by
  rcases i with ((t'|t')|q)
  · have h' := congrFun (congrFun h (lX t')) (lY t)
    simpa [PY, Matrix.mul_diagonal, iY, Matrix.diagonal_apply] using h'
  · simp [iY] at hi
  · exact col_y_zero h q t

variable {E : Finset Sl} {M Mi : Matrix (L3 T Sl) (L3 T Sl) ℚ}

/-- an exit column of the inverse is zero at every non-y row other than the exit itself. -/
lemma exit_col_invZ (him : Mi * M = 1) (hmi : M * Mi = 1) (hR : ExitRows E M) (hC : ExitCols E M)
    (hX : ExitColsX E M) (hiy : Mi * PY = PY)
    {e : Sl} (he : e ∈ E) {i : L3 T Sl} (hi : iY i = 0) (hie : i ≠ lS e) : Mi i (lS e) = 0 := by
  have hne := (exit_row_inv hmi hR he).1
  have key : (Mi * M) i (lS e) = Mi i (lS e) * M (lS e) (lS e) := by
    rw [Matrix.mul_apply, Finset.sum_eq_single (lS e)]
    · intro k _ hk
      rcases k with ((t|t)|q')
      · rw [hX e he t, mul_zero]
      · rw [col_y_zero' hiy hi t, zero_mul]
      · have hq' : q' ≠ e := fun h => hk (by rw [h])
        rw [hC e he q' hq', mul_zero]
    · intro h; exact absurd (Finset.mem_univ _) h
  have h0 : (Mi * M) i (lS e) = 0 := by
    rw [him, Matrix.one_apply_ne hie]
  rw [key] at h0
  rcases mul_eq_zero.mp h0 with h | h
  · exact h
  · exact absurd h hne

/-- **The forward rest gate `g8z Mti` touches neither y roles nor exit slots.** -/
theorem exit_hGrZ (him : Mi * M = 1) (hmi : M * Mi = 1) (hR : ExitRows E M) (hC : ExitCols E M)
    (hX : ExitColsX E M) (hiy : Mi * PY = PY) :
    ∀ i j, i ≠ j → XMut.g8z Mi i j ≠ 0 → restOk E i ∧ restOk E j := by
  intro i j e h
  have f := XMut.g8z_facts hiy
  obtain ⟨hiy', hjy⟩ := XMut.offdiag_z f.1 f.2 h e
  rw [g8z_off _ e, hiy', sub_zero, one_mul] at h
  refine ⟨⟨hiy', ?_⟩, ⟨hjy, ?_⟩⟩
  · rcases i with ((t|t)|q)
    · simp [iE]
    · simp [iY] at hiy'
    · rw [iE_lS]
      by_cases hq : q ∈ E
      · exact absurd ((exit_row_inv hmi hR hq).2 j (Ne.symm e)) h
      · rw [if_neg hq]
  · rcases j with ((t|t)|e')
    · simp [iE]
    · simp [iY] at hjy
    · rw [iE_lS]
      by_cases he' : e' ∈ E
      · exact absurd (exit_col_invZ him hmi hR hC hX hiy he' hiy' e) h
      · rw [if_neg he']

/-- **The backward rest gate `(h8z Mt)ᵀ` touches neither y roles nor exit slots.** -/
theorem exit_hHrZ (hmy : M * PY = PY) (hR : ExitRows E M) (hC : ExitCols E M) (hX : ExitColsX E M) :
    ∀ i j, i ≠ j → (XMut.h8z M)ᵀ i j ≠ 0 → restOk E i ∧ restOk E j := by
  intro i j e h
  have h' : XMut.h8z M j i ≠ 0 := h
  have f := XMut.h8z_facts hmy
  obtain ⟨hjy, hiy⟩ := XMut.offdiag_z f.1 f.2 h' (Ne.symm e)
  rw [h8z_off _ (Ne.symm e), hjy, sub_zero, one_mul] at h'
  refine ⟨⟨hiy, ?_⟩, ⟨hjy, ?_⟩⟩
  · rcases i with ((t|t)|e')
    · simp [iE]
    · simp [iY] at hiy
    · rw [iE_lS]
      by_cases he' : e' ∈ E
      · rcases j with ((t|t)|q)
        · exact absurd (hX e' he' t) h'
        · simp [iY] at hjy
        · have hqe : q ≠ e' := fun h2 => e (by rw [h2])
          exact absurd (hC e' he' q hqe) h'
      · rw [if_neg he']
  · rcases j with ((t|t)|q)
    · simp [iE]
    · simp [iY] at hjy
    · rw [iE_lS]
      by_cases hq : q ∈ E
      · exact absurd (hR q hq i e) h'
      · rw [if_neg hq]

/-- the two fields for a gauged mutable-source scalar half. -/
theorem exit_hGrZ_scal {C : Type} [Fintype C] [DecidableEq C] (S : GScalZ T Sl C)
    (hR : ExitRows E S.toScalZ.Mt) (hC : ExitCols E S.toScalZ.Mt) (hX : ExitColsX E S.toScalZ.Mt) :
    ∀ i j, i ≠ j → XMut.g8z S.toScalZ.Mti i j ≠ 0 → restOk E i ∧ restOk E j :=
  exit_hGrZ S.toScalZ.im S.toScalZ.mi hR hC hX S.toScalZ.hiy

theorem exit_hHrZ_scal {C : Type} [Fintype C] [DecidableEq C] (S : GScalZ T Sl C)
    (hR : ExitRows E S.toScalZ.Mt) (hC : ExitCols E S.toScalZ.Mt) (hX : ExitColsX E S.toScalZ.Mt) :
    ∀ i j, i ≠ j → (XMut.h8z S.toScalZ.Mt)ᵀ i j ≠ 0 → restOk E i ∧ restOk E j :=
  exit_hHrZ S.toScalZ.hmy hR hC hX

end ExitEndZ
end
end XEX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XEX.exit_hGrZ
#print axioms OAI.PowerSaving.XEX.exit_hHrZ
#print axioms OAI.PowerSaving.XEX.exit_hGrZ_scal
