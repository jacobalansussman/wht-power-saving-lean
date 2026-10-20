import Work.Carrier.Scat

/-!
# x-mutable (13 Oct 2026): the dirt rule WITHOUT read-only sources

Jacob Sussman's `CR.dirt_fwd` / `CR.dirt_bwd` (`Work/Carrier/Ring.lean`) assume `px * m = px`: no gate of the
main phase writes a source (x) role, so the x rows of the product `m` are unit rows.  Here that hypothesis is
DROPPED.  The derived undo gate at the end restores the x rows together with the slot rows:

    g8z mi = py + (1 - py) * mi        (the x AND slot rows of `m⁻¹`),
    h8z m  = py + (1 - py) * m         (its inverse, for the transposed backward word),

and the forward / backward words are still the shears `1 + py m px` / `1 - py m px`.  Only `m * py = py`
(no target feeds a source or a helper; the y columns are unit columns) and `mi * m = 1` are used.

Both derived gates have unit y rows and unit y columns, so every off-diagonal entry joins two NON-target roles
(`offdiag_z`): at the end of an invocation those are the sources and the slots, which all stand at the old frame
of the full label, so the gates are legal there (as Sussman's `g8` is).

No `sorry`.  Credits: engine and the original dirt rule: Jacob Sussman (Apache-2.0).
-/

namespace XMut
open OAI.PowerSaving OAI.PowerSaving.CR Matrix

section Ring
variable {R : Type*} [Ring R] {px py m mi : R}

/-- **The dirt rule, forward, with mutable sources.** -/
theorem dirt_fwdZ (yy : py * py = py) (xy : px * py = 0) (yx : py * px = 0)
    (my : m * py = py) (im : mi * m = 1) :
    (py + (1 - py) * mi) * m * (1 - py * m * (1 - px - py)) = 1 + py * m * px := by
  have yy' : ∀ z, py * (py * z) = py * z := fun z => by rw [← mul_assoc, yy]
  have xy' : ∀ z, px * (py * z) = 0 := fun z => by rw [← mul_assoc, xy, zero_mul]
  have yx' : ∀ z, py * (px * z) = 0 := fun z => by rw [← mul_assoc, yx, zero_mul]
  have my' : ∀ z, m * (py * z) = py * z := fun z => by rw [← mul_assoc, my]
  have im' : ∀ z, mi * (m * z) = z := fun z => by rw [← mul_assoc, im, one_mul]
  simp only [mul_add, add_mul, mul_sub, sub_mul, mul_one, one_mul, mul_zero, zero_mul,
    mul_assoc, yy, xy, yx, my, im, yy', xy', yx', my', im', sub_zero, add_zero, zero_add]
  abel

/-- **The dirt rule, backward, with mutable sources** (the inverse of the forward product). -/
theorem dirt_bwdZ (yy : py * py = py) (xy : px * py = 0) (yx : py * px = 0)
    (my : m * py = py) (im : mi * m = 1) :
    (1 + py * m * (1 - px - py)) * mi * (py + (1 - py) * m) = 1 - py * m * px := by
  have iy : mi * py = py := dirt_iy my im
  have yy' : ∀ z, py * (py * z) = py * z := fun z => by rw [← mul_assoc, yy]
  have xy' : ∀ z, px * (py * z) = 0 := fun z => by rw [← mul_assoc, xy, zero_mul]
  have yx' : ∀ z, py * (px * z) = 0 := fun z => by rw [← mul_assoc, yx, zero_mul]
  have my' : ∀ z, m * (py * z) = py * z := fun z => by rw [← mul_assoc, my]
  have im' : ∀ z, mi * (m * z) = z := fun z => by rw [← mul_assoc, im, one_mul]
  have iy' : ∀ z, mi * (py * z) = py * z := fun z => by rw [← mul_assoc, iy]
  simp only [mul_add, add_mul, mul_sub, sub_mul, mul_one, one_mul, mul_zero, zero_mul,
    mul_assoc, yy, xy, yx, my, im, iy, yy', xy', yx', my', im', iy', sub_zero, add_zero, zero_add]
  abel

/-- the y rows of the undo gates are unit rows -/
lemma sideZ_y (yy : py * py = py) (w : R) : py * (py + (1 - py) * w) = py := by
  rw [mul_add, yy, ← mul_assoc, mul_sub, mul_one, yy, sub_self, zero_mul, add_zero]

/-- the y columns of the undo gates are unit columns when `w` has unit y columns -/
lemma sideZ_ycol (yy : py * py = py) {w : R} (hw : w * py = py) : (py + (1 - py) * w) * py = py := by
  rw [add_mul, yy, mul_assoc, hw, sub_mul, one_mul, yy, sub_self, add_zero]

end Ring

section Mat
variable {T Sl : Type} [Fintype T] [DecidableEq T] [Fintype Sl] [DecidableEq Sl]

/-- the derived undo gate of the forward word: the x and slot rows of `Mi`. -/
def g8z (Mi : Matrix (L3 T Sl) (L3 T Sl) ℚ) : Matrix (L3 T Sl) (L3 T Sl) ℚ := PY + (1 - PY) * Mi
/-- its inverse, for the transposed backward word. -/
def h8z (M : Matrix (L3 T Sl) (L3 T Sl) ℚ) : Matrix (L3 T Sl) (L3 T Sl) ℚ := PY + (1 - PY) * M

lemma g8z_facts {Mi : Matrix (L3 T Sl) (L3 T Sl) ℚ} (hiy : Mi * PY = PY) :
    PY * g8z Mi = PY ∧ g8z Mi * PY = PY :=
  ⟨sideZ_y PY_sq Mi, sideZ_ycol PY_sq hiy⟩

lemma h8z_facts {M : Matrix (L3 T Sl) (L3 T Sl) ℚ} (hmy : M * PY = PY) :
    PY * h8z M = PY ∧ h8z M * PY = PY :=
  ⟨sideZ_y PY_sq M, sideZ_ycol PY_sq hmy⟩

/-- a gate with unit y rows and unit y columns: a nonzero entry off the diagonal joins two non-target roles. -/
lemma offdiag_z {G : Matrix (L3 T Sl) (L3 T Sl) ℚ} (hy : PY * G = PY) (hyc : G * PY = PY)
    {i j : L3 T Sl} (hij : G i j ≠ 0) (e : i ≠ j) : iY i = 0 ∧ iY j = 0 := by
  refine ⟨?_, ?_⟩
  · by_contra h
    exact hij (row_zero_of hy ((iY_cases i).resolve_left h) e)
  · by_contra h
    exact hij (col_zero_of hyc ((iY_cases j).resolve_left h) e)

/-- the forward dirt rule for matrices, with mutable sources -/
lemma hFz {M Mi : Matrix (L3 T Sl) (L3 T Sl) ℚ} (my : M * PY = PY) (im : Mi * M = 1) :
    g8z Mi * M * g1 M = 1 + PY * M * PX := by
  unfold g8z g1
  exact dirt_fwdZ PY_sq PX_PY PY_PX my im

/-- the backward dirt rule for matrices, with mutable sources -/
lemma hBz {M Mi : Matrix (L3 T Sl) (L3 T Sl) ℚ} (my : M * PY = PY) (im : Mi * M = 1) :
    g1i M * Mi * h8z M = 1 - PY * M * PX := by
  unfold h8z g1i
  exact dirt_bwdZ PY_sq PX_PY PY_PX my im

end Mat
end XMut

#print axioms XMut.dirt_fwdZ
#print axioms XMut.dirt_bwdZ
#print axioms XMut.hFz
#print axioms XMut.hBz
#print axioms XMut.offdiag_z
