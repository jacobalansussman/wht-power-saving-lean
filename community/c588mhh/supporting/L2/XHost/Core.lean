import Work.GFrame.Engine.Free

/-!
# x-hosting, core lemmas: idle-axis hosting in our engine

Idle-axis hosting (DaysSky, CrocSwap/integer-mult-bounds PR #442, research/idle-axis-hosting/PROOF.md) puts a dirty-helper
residual of an invocation on the PENDING coordinates of an idle data role instead of a bank role.  In our engine
(Jacob Sussman's GFrame engine, wht-power-saving-lean, Apache-2.0; as used unchanged by Ryan Shea's complex588 and Nil
Mamano's gauged extension) the invocation package already admits ANY slot frame (`GInvE.fwd`: `∀ Z`, slot ends in
`Wn * Z`, scalar map `shearF` keeps slot scalars).  What hosting needs on top of that is:

* `GWord.Avoids`, `gwalk_roleAct_comm`  the WINDOW: a word that neither acts on role `d` nor reads it commutes with any
                                        matrix applied to `d` (hub PROOF.md §2, "nothing touches F between ι and τ");
* `blockMat_comm`, `blockMat_append`    blocks along independent directions commute and concatenate (hub §1(b));
* `host_remaining_climb`                a host whose idle climb is ONE block along `Fin.append z y` and which received the
                                        residual `blockMat z` from an invocation finishes with ONE block of rank `card y`;
* `host_remaining_climb_twisted`        the same for a gauged residual `C * blockMat z * C'` with free `C`, `C'`,
                                        `C' * C = 1` (hub §2: "σ composed with U_q ⊕ I").

No `sorry`.  Written for the wht43 x-hosting study (13 Oct 2026), prepared with assistance from Claude.
-/

set_option linter.unusedSectionVars false

namespace OAI
namespace PowerSaving.GF
open Binary Matrix Finset RAM
noncomputable section

section Window
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- A letter AVOIDS role `d`: it does not act on `d`, and a gate neither reads `d` into another role nor changes `d`. -/
def GMove.Avoids (d : ρ) : GMove α ρ → Prop
  | .old (.localM g) => (∀ j, g d j = if j = d then 1 else 0) ∧ (∀ i, i ≠ d → g i d = 0)
  | .old (.shift l _) => l ≠ d
  | .old (.block l _ _ _) => l ≠ d
  | .perm l _ => l ≠ d
  | .phase l _ _ => l ≠ d

/-- A word avoids `d` when every letter does. -/
def GWord.Avoids (d : ρ) (w : GWord α ρ) : Prop := ∀ mv ∈ w, GMove.Avoids d mv

lemma roleAct_comm_of_ne (f : ℕ) {l d : ρ} (h : l ≠ d) (M A : CMat α) (x : Data α ρ f) :
    roleAct f l M (roleAct f d A x) = roleAct f d A (roleAct f l M x) := by
  funext j
  by_cases hl : j = l
  · subst hl
    simp [roleAct, h]
  · by_cases hd : j = d
    · subst hd
      simp [roleAct, hl]
    · simp [roleAct, hl, hd]

lemma point_actPoint_comm (f : ℕ) (d : ρ) (g : Matrix ρ ρ ℚ)
    (hrow : ∀ j, g d j = if j = d then 1 else 0) (hcol : ∀ i, i ≠ d → g i d = 0)
    (A : CMat α) (x : Data α ρ f) :
    point f (actPoint g) (roleAct f d A x) = roleAct f d A (point f (actPoint g) x) := by
  have key : ∀ (y : Data α ρ f) (i : Sky α f), actPoint g (fun r => y r i) d = y d i := by
    intro y i
    unfold actPoint
    rw [Finset.sum_eq_single d]
    · simp [hrow]
    · intro j _ hj
      simp [hrow, hj]
    · intro h
      exact absurd (Finset.mem_univ d) h
  have hd : point f (actPoint g) x d = x d := by
    funext i
    exact key x i
  funext r i
  by_cases hr : r = d
  · subst hr
    have e1 : point f (actPoint g) (roleAct f r A x) r i = roleAct f r A x r i := key (roleAct f r A x) i
    rw [e1]
    simp [roleAct, hd]
  · have e : roleAct f d A (point f (actPoint g) x) r = point f (actPoint g) x r := by
      simp [roleAct, hr]
    rw [e]
    show actPoint g (fun r' => roleAct f d A x r' i) r = actPoint g (fun r' => x r' i) r
    unfold actPoint
    refine Finset.sum_congr rfl (fun j _ => ?_)
    by_cases hj : j = d
    · subst hj
      simp [hcol r hr]
    · simp [roleAct, hj]

/-- **One letter that avoids `d` commutes with any matrix on `d`.** -/
lemma GMove.act_comm_of_avoids (f : ℕ) (d : ρ) (mv : GMove α ρ) (h : GMove.Avoids d mv) (A : CMat α)
    (x : Data α ρ f) : mv.act f (roleAct f d A x) = roleAct f d A (mv.act f x) := by
  rcases mv with (mv | ⟨l, G⟩ | ⟨l, z, c⟩)
  · rcases mv with (g | ⟨l, z⟩ | ⟨l, rk, z, ind⟩)
    · obtain ⟨h1, h2⟩ := h
      exact point_actPoint_comm f d g h1 h2 A x
    · exact roleAct_comm_of_ne f h _ A x
    · exact roleAct_comm_of_ne f h _ A x
  · exact roleAct_comm_of_ne f h _ A x
  · exact roleAct_comm_of_ne f h _ A x

/-- **The window.**  A word that avoids role `d` commutes with any matrix applied to `d`: a residual applied to an idle
host at the hosting invocation may be moved to the host's idle climb, whatever runs in between. -/
theorem gwalk_roleAct_comm (f : ℕ) (d : ρ) (w : GWord α ρ) (h : GWord.Avoids d w) (A : CMat α)
    (x : Data α ρ f) : gwalk f w (roleAct f d A x) = roleAct f d A (gwalk f w x) := by
  induction w generalizing x with
  | nil => rfl
  | cons mv w ih =>
    rw [gwalk_cons, gwalk_cons, GMove.act_comm_of_avoids f d mv (h mv List.mem_cons_self) A x]
    exact ih (fun m hm => h m (List.mem_cons_of_mem _ hm)) _

end Window

section Blocks
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

lemma blockMat_comm {r s : ℕ} (z : Fin r → Space α) (y : Fin s → Space α) :
    blockMat z * blockMat y = blockMat y * blockMat z := by
  rw [blockMat_wrap, blockMat_wrap]
  apply wrap_comm

lemma blockMat_append {r s : ℕ} (z : Fin r → Space α) (y : Fin s → Space α) :
    blockMat (Fin.append z y) = blockMat z * blockMat y := by
  unfold blockMat
  rw [List.ofFn_add, List.prod_append]
  simp [Fin.append_left, Fin.append_right]

/-- **Hosting, the remaining climb.**  A host role whose idle climb is ONE block along `Fin.append z y` (rank
`r + s`), and which already received the residual `blockMat z` (from an invocation that used it as a dirty helper),
reaches the same frame with ONE block of rank `s`.  Price `φ s` instead of `φ (r + s)`. -/
theorem host_remaining_climb (l : ρ) (Y : CMat α) {r s : ℕ} (z : Fin r → Space α) (y : Fin s → Space α)
    (ind : LinearIndependent F y) (h1 : 1 ≤ s) (h2 : s < Fintype.card α) :
    GStep l (blockMat z * Y) (blockMat (Fin.append z y) * Y) (fun φ => φ s) := by
  have h := gstep_block (ρ := ρ) l (blockMat z * Y) y ind h1 h2
  have e : blockMat y * (blockMat z * Y) = blockMat (Fin.append z y) * Y := by
    rw [blockMat_append, ← mul_assoc, blockMat_comm y z]
  rw [e] at h
  exact h

/-- **Hosting a gauged residual.**  The residual is `C * blockMat z * C'` with free adapters (`C' * C = 1`); the host's
idle climb is `C * blockMat (Fin.append z y) * C'`; the remaining climb is ONE block of rank `s`. -/
theorem host_remaining_climb_twisted (l : ρ) (Y : CMat α) {C C' : CMat α} (hC : Free C) (hC' : Free C')
    (hCC : C' * C = 1) {r s : ℕ} (z : Fin r → Space α) (y : Fin s → Space α)
    (ind : LinearIndependent F y) (h1 : 1 ≤ s) (h2 : s < Fintype.card α) :
    GStep l (C * blockMat z * C' * Y) (C * blockMat (Fin.append z y) * C' * Y) (fun φ => φ s) := by
  have h := gstep_of_factor (ρ := ρ) l (C * blockMat z * C' * Y) y ind hC hC' h1 h2
  have e : C * blockMat y * C' * (C * blockMat z * C' * Y) = C * blockMat (Fin.append z y) * C' * Y := by
    calc C * blockMat y * C' * (C * blockMat z * C' * Y)
        = C * blockMat y * (C' * C) * blockMat z * C' * Y := by simp only [mul_assoc]
      _ = C * (blockMat y * blockMat z) * C' * Y := by rw [hCC]; simp only [mul_one, mul_assoc]
      _ = C * blockMat (Fin.append z y) * C' * Y := by rw [blockMat_comm y z, ← blockMat_append]
  rw [e] at h
  exact h

/-- **Hosting, whole idle climb.**  With the residual applied, the host's full climb splits as the residual (paid by the
invocation, rank `r`) followed by the remaining block (rank `s`): composing with any route that took the host from `Y`
to `blockMat z * Y` at price `c` gives the old target at price `c + φ s`. -/
theorem host_route_split (l : ρ) (Y : CMat α) {r s : ℕ} (z : Fin r → Space α) (y : Fin s → Space α)
    (ind : LinearIndependent F y) (h1 : 1 ≤ s) (h2 : s < Fintype.card α) {c : (ℕ → ℝ) → ℝ}
    (hres : GStep l Y (blockMat z * Y) c) :
    GStep l Y (blockMat (Fin.append z y) * Y) (fun φ => c φ + φ s) :=
  hres.trans (host_remaining_climb l Y z y ind h1 h2)

end Blocks
end
end PowerSaving.GF
end OAI

#print axioms OAI.PowerSaving.GF.gwalk_roleAct_comm
#print axioms OAI.PowerSaving.GF.host_remaining_climb
#print axioms OAI.PowerSaving.GF.host_remaining_climb_twisted
#print axioms OAI.PowerSaving.GF.host_route_split
