import Work.KLow.ZCheck
import Work.KLow.LabelsX

/-!
# (key: k-lower) The y-free twin of a certificate: same non-y rows of the total matrix

k-lower, 13 Oct 2026.  `nyRaw c` = `c` with every add whose TARGET is a y register removed from phase B and an empty
scatter (v empty rows).  Since no add of B with a non-y target reads a y register (Sussman's kind `okBK`), and the
scatter writes only y roles, the non-y rows of the total matrix are unchanged (`nyRaw_Mt`).  The zero checks
(Work/KLow/ZCheck) run on the twin: its y registers stay empty, so the packed tags stay small.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option linter.deprecated false

namespace GXD
open SSC OAI.PowerSaving OAI.PowerSaving.CR OAI.PowerSaving.CB Matrix Finset

/-- `v ≤ r < 2 v` -/
def isYN (v r : Nat) : Bool := Nat.ble v r && Nat.blt r (2 * v)

/-- a gate without its adds into y registers -/
def nyGate (v : Nat) : Gate → Gate
  | .out f s ts ex => .out f s (ts.filter (fun p => !isYN v p.1)) ex
  | .inn f t ss => .inn f t (if isYN v t then [] else ss)

/-- **the y-free twin** -/
def nyRaw (c : Raw) : Raw :=
  { c with B := c.B.map (fun gs => gs.map (nyGate c.v)), scat := .table (List.replicate c.v []) }

theorem nyGate_adds (v : Nat) (g : Gate) :
    (nyGate v g).adds = g.adds.filter (fun x => !isYN v x.1) := by
  cases g with
  | out f s ts ex =>
    simp only [nyGate, Gate.adds, List.filter_map]
    rfl
  | inn f t ss =>
    by_cases h : isYN v t = true
    · simp [nyGate, Gate.adds, h, List.filter_map]
    · have h' : isYN v t = false := by simpa using h
      have hf : List.filter ((fun x => !isYN v x.1) ∘ fun p => (t, p)) ss = ss :=
        List.filter_eq_self.mpr (fun p _ => by simp [Function.comp, h'])
      simp only [nyGate, Gate.adds, h', Bool.false_eq_true, if_false, List.filter_map, hf]

theorem addsB_ny (c : Raw) : GS.addsB (nyRaw c) = (GS.addsB c).filter (fun x => !isYN c.v x.1) := by
  unfold GS.addsB Raw.gatesB
  show ((c.B.map (fun gs => gs.map (nyGate c.v))).flatten.flatMap Gate.adds) = _
  rw [← List.map_flatten, List.flatMap_map, List.filter_flatMap]
  congr 1
  funext g
  exact nyGate_adds c.v g

section Proj
variable {ρ : Type} [Fintype ρ] [DecidableEq ρ]

lemma diag_addMat_y (d : ρ → ℚ) (a b : ρ) (cc : ℚ) (ha : d a = 0) :
    Matrix.diagonal d * addMat a b cc = Matrix.diagonal d := by
  ext i j
  rw [Matrix.diagonal_mul]
  simp only [addMat, Matrix.diagonal_apply]
  by_cases h : i = a ∧ j = b
  · rw [h.1, ha]; simp
  · rw [if_neg h, add_zero]
    by_cases e : i = j <;> simp [e]

lemma diag_addMat_comm (d : ρ → ℚ) (a b : ρ) (cc : ℚ) (hab : d a = d b) :
    Matrix.diagonal d * addMat a b cc = addMat a b cc * Matrix.diagonal d := by
  ext i j
  rw [Matrix.diagonal_mul, Matrix.mul_diagonal]
  simp only [addMat]
  by_cases h : i = a ∧ j = b
  · rw [if_pos h, h.1, h.2, hab]
    by_cases e : a = b
    · subst e; ring
    · rw [if_neg e]; ring
  · rw [if_neg h, add_zero]
    by_cases e : i = j
    · subst e; simp
    · simp [e]

/-- the non-y rows of a gate product: drop the adds into y roles -/
theorem diag_matP (d : ρ → ℚ) (un : Nat → ρ) (isy : Nat → Bool) (hd : ∀ r, isy r = true → d (un r) = 0)
    (hd1 : ∀ r, isy r = false → d (un r) = 1) :
    ∀ l : List (Nat × Nat × Co), (∀ x ∈ l, isy x.2.1 = true → isy x.1 = true) →
      Matrix.diagonal d * matP un (GS.mic l) = matP un (GS.mic (l.filter (fun x => !isy x.1))) * Matrix.diagonal d := by
  intro l
  induction l with
  | nil => intro _; simp [GS.mic, matP]
  | cons x l ih =>
    intro hl
    have hl' : ∀ y ∈ l, isy y.2.1 = true → isy y.1 = true := fun y hy => hl y (List.mem_cons_of_mem _ hy)
    show Matrix.diagonal d * (matP un (GS.mic l) * addMat (un x.1) (un x.2.1) (GS.toCoef x.2.2).val) = _
    rw [← Matrix.mul_assoc, ih hl', Matrix.mul_assoc]
    cases hx : isy x.1 with
    | true =>
      rw [diag_addMat_y d _ _ _ (hd _ hx)]
      have e : (x :: l).filter (fun x => !isy x.1) = l.filter (fun x => !isy x.1) := by
        rw [List.filter_cons]; simp [hx]
      rw [e]
    | false =>
      have hs : isy x.2.1 = false := by
        cases h2 : isy x.2.1 with
        | false => rfl
        | true => have := hl x List.mem_cons_self h2; rw [hx] at this; exact absurd this Bool.false_ne_true
      rw [diag_addMat_comm d _ _ _ (by rw [hd1 _ hx, hd1 _ hs]), ← Matrix.mul_assoc]
      have e : (x :: l).filter (fun x => !isy x.1) = x :: l.filter (fun x => !isy x.1) := by
        rw [List.filter_cons]; simp [hx]
      rw [e]
      rfl

end Proj

/-- the non-y projector of the live roles -/
def dN {T Sl : Type} : L3 T Sl → ℚ := fun i => 1 - iY i

theorem isYN_unemb (c : Raw) (hv : 0 < c.v) (r : Nat) :
    (isYN c.v r = true → dN (GS.unemb c hv r) = 0) ∧ (isYN c.v r = false → dN (GS.unemb c hv r) = 1) := by
  unfold isYN
  constructor
  · intro h
    simp only [Bool.and_eq_true, Nat.ble_eq, Nat.blt_eq] at h
    rw [GS.unemb_y hv r h.1 h.2]
    simp [dN, iY, lY]
  · intro h
    by_cases h1 : r < c.v
    · rw [GS.unemb_x hv r h1]; simp [dN, iY, lX]
    · by_cases h2 : r < 2 * c.v
      · exfalso
        have : (Nat.ble c.v r && Nat.blt r (2 * c.v)) = true := by
          simp only [Bool.and_eq_true, Nat.ble_eq, Nat.blt_eq]; omega
        rw [h] at this; exact Bool.false_ne_true this
      · unfold GS.unemb
        rw [if_neg h1, if_neg h2]
        split <;> simp [dN, iY]

/-- **the twin keeps the non-y rows of the total matrix** -/
theorem nyRaw_Mt (c : Raw) (hv : 0 < c.v)
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (i j : GS.Rl c) (hi : iY i = 0) :
    GS.Mt c hv i j = GS.Mt (nyRaw c) hv i j := by
  have hd : ∀ r, isYN c.v r = true → dN (GS.unemb c hv r) = 0 := fun r => (isYN_unemb c hv r).1
  have hd1 : ∀ r, isYN c.v r = false → dN (GS.unemb c hv r) = 1 := fun r => (isYN_unemb c hv r).2
  have hBy : ∀ x ∈ GS.addsB c, isYN c.v x.2.1 = true → isYN c.v x.1 = true := by
    intro x hx h
    simp only [isYN, Bool.and_eq_true, Nat.ble_eq, Nat.blt_eq] at h ⊢
    exact hB x hx h.1 h.2
  have key := diag_matP dN (GS.unemb c hv) (isYN c.v) hd hd1 (GS.addsB c) hBy
  have hBy' : ∀ x ∈ (GS.addsB c).filter (fun x => !isYN c.v x.1), isYN c.v x.2.1 = true → isYN c.v x.1 = true := by
    intro x hx h
    exact hBy x (List.mem_filter.mp hx).1 h
  have key2 := diag_matP dN (GS.unemb c hv) (isYN c.v) hd hd1 _ hBy'
  rw [List.filter_filter] at key2
  simp only [Bool.and_self] at key2
  -- the scatter of the twin is the identity; the non-y rows ignore the scatter of `c`
  have hJ : GS.Jrm (nyRaw c) = 0 := by
    funext S k
    show rowOf ((((List.replicate c.v ([] : List (Nat × Co))) : List (List (Nat × Co))).getD S.val []).map
      fun x => (x.1, GS.toCoef x.2)) k = 0
    rw [List.getD_eq_getElem _ _ (by rw [List.length_replicate]; exact S.isLt), List.getElem_replicate]
    rfl
  have hS1 : @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (scatM (GS.Jrm (nyRaw c) * GS.Ccm (nyRaw c))) = 1 := by
    rw [hJ, Matrix.zero_mul]
    unfold scatM
    dsimp only [id]
    ext a b
    rw [Matrix.add_apply]
    rcases a with ((t|t)|q) <;> rcases b with ((t'|t')|q') <;> simp [embYS, Matrix.one_apply] <;> rfl
  have hSy : ∀ K : Matrix (Fin c.v) (Fin c.R) ℚ, Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * scatM K = Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) := by
    intro K
    unfold scatM
    rw [Matrix.mul_add, Matrix.mul_one, add_eq_left]
    ext a b
    rw [Matrix.diagonal_mul]
    rcases a with ((t|t)|q)
    · simp [embYS]
    · simp [dN, iY]
    · simp [embYS]
  have eMB : GS.MB c hv = matP (GS.unemb c hv) (GS.mic (GS.addsB c)) := rfl
  have eB : @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (GS.MB (nyRaw c) hv) = matP (GS.unemb c hv) (GS.mic ((GS.addsB c).filter (fun x => !isYN c.v x.1))) := by
    show matP (GS.unemb (nyRaw c) hv) (GS.mic (GS.addsB (nyRaw c))) = _
    rw [addsB_ny]
    rfl
  have e2 : @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (GS.Mt (nyRaw c) hv)
      = matP (GS.unemb c hv) (GS.mic ((GS.addsB c).filter (fun x => !isYN c.v x.1))) * GS.MA c hv := by
    have e : @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (GS.Mt (nyRaw c) hv) = @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (GS.MB (nyRaw c) hv)
        * @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (scatM (GS.Jrm (nyRaw c) * GS.Ccm (nyRaw c))) * GS.MA c hv := rfl
    rw [e, hS1, Matrix.mul_one, eB]
  have hMt : Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * GS.Mt c hv = Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (GS.Mt (nyRaw c) hv) := by
    have e1 : GS.Mt c hv = GS.MB c hv * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv := rfl
    rw [e1, e2, eMB]
    calc Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * (matP (GS.unemb c hv) (GS.mic (GS.addsB c)) * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv)
        = (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * matP (GS.unemb c hv) (GS.mic (GS.addsB c))) * scatM (GS.Jrm c * GS.Ccm c)
            * GS.MA c hv := by simp only [Matrix.mul_assoc]
      _ = matP (GS.unemb c hv) (GS.mic ((GS.addsB c).filter (fun x => !isYN c.v x.1)))
            * (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * scatM (GS.Jrm c * GS.Ccm c)) * GS.MA c hv := by
          rw [key]; simp only [Matrix.mul_assoc]
      _ = (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * matP (GS.unemb c hv) (GS.mic ((GS.addsB c).filter (fun x => !isYN c.v x.1))))
            * GS.MA c hv := by rw [hSy, key2]
      _ = _ := by rw [Matrix.mul_assoc]
  have h1 : (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * GS.Mt c hv) i j = (Matrix.diagonal (dN (T := Fin c.v) (Sl := Fin c.R)) * @id (Matrix (GS.Rl c) (GS.Rl c) ℚ) (GS.Mt (nyRaw c) hv)) i j := by
    rw [hMt]
  simp only [Matrix.diagonal_mul] at h1
  have hdi : dN i = 1 := by simp [dN, hi]
  rw [hdi, one_mul, one_mul] at h1
  exact h1

/-- the exit conditions of `c` from those of its twin (exit rows and slot rows are not y rows) -/
theorem exit_of_twin (c : Raw) (hv : 0 < c.v)
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (hy : ∀ (i : GS.Rl c) (t : Fin c.v), (GS.MB c hv * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv) i (lY t)
      = if i = lY t then 1 else 0)
    (el : List Nat) (h : OAI.PowerSaving.XEX.ExitRows (GLab.exitSet (nyRaw c) el) (GS.Mt (nyRaw c) hv) ∧
      OAI.PowerSaving.XEX.ExitCols (GLab.exitSet (nyRaw c) el) (GS.Mt (nyRaw c) hv)) :
    OAI.PowerSaving.XEX.ExitRows (GLab.exitSet c el) (GS.Mt c hv) ∧
      OAI.PowerSaving.XEX.ExitCols (GLab.exitSet c el) (GS.Mt c hv) := by
  refine ⟨fun e he j hj => ?_, fun e he q hq => ?_⟩
  · rw [nyRaw_Mt c hv hB (lS e) j (by simp [iY, lS])]
    exact h.1 e he j hj
  · rw [nyRaw_Mt c hv hB (lS q) (lS e) (by simp [iY, lS])]
    exact h.2 e he q hq

/-- the y columns of the twin at the slot rows -/
theorem twin_hy (c : Raw) (hv : 0 < c.v)
    (hB : ∀ x ∈ GS.addsB c, c.v ≤ x.2.1 → x.2.1 < 2 * c.v → c.v ≤ x.1 ∧ x.1 < 2 * c.v)
    (hy : ∀ (i : GS.Rl c) (t : Fin c.v), (GS.MB c hv * scatM (GS.Jrm c * GS.Ccm c) * GS.MA c hv) i (lY t)
      = if i = lY t then 1 else 0) :
    ∀ (q : Fin c.R) (t : Fin c.v), GS.Mt (nyRaw c) hv (lS q) (lY t) = 0 := by
  intro q t
  have h0 := hy (lS q) t
  rw [if_neg (by simp [lS, lY])] at h0
  exact (nyRaw_Mt c hv hB (lS q) (lY t) (by simp [iY, lS])).symm.trans h0

end GXD

#print axioms GXD.nyRaw_Mt
