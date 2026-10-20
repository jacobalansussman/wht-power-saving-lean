import Work.GFrame.Engine.Level
import Work.Mixed.Prog.Sem
import Work.Mixed.Prog.Seq

/-!
# Mixed engine, R3 program part 2: the class table (agent e-r3prog, 19 Oct 2026)

A class of the lockstep tree holds `NU` units; unit `p` is a word instance at column `c p` with
`2^(ib p)` instances.  Its data is the CLASS TABLE (a tape of units, each a tape of instances,
each the role table of one instance, no leftover bits):

* `ctab X R NU c ib V`  the class table of the element data `V p q : Grain α ρ 0 (c p)`;
* `UnitF NU`, `ElemF NU ib`  the families of units and of elements (instances) of a family of classes;
* `CCtx`              the per-unit tables a class keeps in its environment (columns, powers,
                      XOR tables, `lum` tables), from which every element is a `GRig`;
* `CCtx.erig`         the element `GRig` (column `c p`, no leftover bits);
* `ctab_build`        a class table is rebuilt from per-element role tables in linear work;
* `CCtx.free_letter`  a free letter (`localM`, `shift`, `perm`, `phase`) is one pass over the table.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF
noncomputable section

universe U

section Tab
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]

/-- The class table of element data `V`. -/
def ctab (X : Layout α) (R : Layout ρ) (NU : ℕ) (c ib : ℕ → ℕ)
    (V : ∀ p : Fin NU, Bits (ib p) → Grain α ρ 0 (c p)) : Tape (Tape (Tape ℂ)) :=
  (Layout.std NU).tape (fun p => (bits (ib p)).tape (fun q => (fleet X R 0 (c p)).tape (V p q)))

end Tab

section Fam
variable {ι : Type U}

/-- Units of a family of classes. -/
abbrev UnitF (NU : ι → ℕ) := Σ i, Fin (NU i)
/-- Elements (instances of units) of a family of classes. -/
abbrev ElemF (NU : ι → ℕ) (ib : ι → ℕ → ℕ) := Σ u : UnitF NU, Bits (ib u.1 u.2)

end Fam

/-- Tables a class keeps for its units: column, `2^column`, instance bits, `2^instance bits`,
XOR table and `lum` table of the column. -/
structure CCtx {ι : Type U} (m : Bool) (B : ι → ℕ) (s : Ty) (x : ι → T s) (mx : ℕ)
    (n NU : ι → ℕ) (col ib : ι → ℕ → ℕ) : Prop where
  nK : Knows m B s x n
  p2n : Knows m B s x (fun i => 2^(n i))
  nuK : Knows m B s x NU
  imag : Has m B s sc x (fun _ => Complex.I)
  tc : Has m B s (a w) x (fun i => (Layout.std (NU i)).tape (fun p => col i p))
  tpc : Has m B s (a w) x (fun i => (Layout.std (NU i)).tape (fun p => 2^(col i p)))
  tib : Has m B s (a w) x (fun i => (Layout.std (NU i)).tape (fun p => ib i p))
  tpib : Has m B s (a w) x (fun i => (Layout.std (NU i)).tape (fun p => 2^(ib i p)))
  txr : Has m B s (a (a w)) x (fun i => (Layout.std (NU i)).tape (fun p => xorTab (col i p)))
  tlum : Has m B s (a (a sc)) x
    (fun i => (Layout.std (NU i)).tape (fun p => (bits (col i p)).tape (fun y : Bits (col i p) => lum y)))
  big : Small B (fun i => 2^(n i))
  nuS : Small B NU
  fit : ∀ i (p : Fin (NU i)), mx * col i p + ib i p = n i
  mx1 : 1 ≤ mx

namespace CCtx
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s t : Ty} {x : ι → T s} {mx : ℕ}
  {n NU : ι → ℕ} {col ib : ι → ℕ → ℕ}
variable (h : CCtx m B s x mx n NU col ib)
include h

/-- The tables survive extending the environment. -/
lemma ext {t : Ty} (y : ι → T t) :
    CCtx m B (p s t) (fun i => (x i, y i)) mx n NU col ib := by
  have F : ∀ {u : Ty} {z : ι → T u}, Has m B s u x z → Has m B (p s t) u (fun i => (x i, y i)) z :=
    fun g => Can.first.then_do g
  exact ⟨F h.nK, F h.p2n, F h.nuK, F h.imag, F h.tc, F h.tpc, F h.tib, F h.tpib, F h.txr, F h.tlum,
    h.big, h.nuS, h.fit, h.mx1⟩

omit h in
/-- A table entry of unit `p`, in the environment of a unit (`(x, p)`). -/
lemma _root_.OAI.PowerSaving.RAM.unit_fetch {t : Ty} (g : ∀ i, Fin (NU i) → T t)
    (ht : Has m B s (a t) x (fun i => (Layout.std (NU i)).tape (g i))) :
    Has m (fun u : UnitF NU => B u.1) (p s w) t (fun u => (x u.1, u.2.val)) (fun u => g u.1 u.2) :=
  Can.fetch (t:=t) (L:=fun u : UnitF NU => Layout.std (NU u.1)) (v:=fun u => g u.1) (i:=fun u => u.2)
    (Can.first.then_do (ht.reindex (fun u : UnitF NU => u.1))) Can.second

lemma two_c_le (i : ι) (p : Fin (NU i)) : 2^(col i p) ≤ 2^(n i) := by
  apply Nat.pow_le_pow_right (by omega)
  have h1 := h.fit i p
  have h2 : col i p ≤ mx * col i p := Nat.le_mul_of_pos_left _ h.mx1
  omega

lemma two_ib_le (i : ι) (p : Fin (NU i)) : 2^(ib i p) ≤ 2^(n i) := by
  apply Nat.pow_le_pow_right (by omega)
  have h1 := h.fit i p
  omega

lemma uBig : Small (fun u : UnitF NU => B u.1) (fun u => 2^(col u.1 u.2)) :=
  (h.big.reindex (fun u : UnitF NU => u.1)).mono (fun u => h.two_c_le u.1 u.2)

lemma uBigIb : Small (fun u : UnitF NU => B u.1) (fun u => 2^(ib u.1 u.2)) :=
  (h.big.reindex (fun u : UnitF NU => u.1)).mono (fun u => h.two_ib_le u.1 u.2)

/-- The environment of an element `((x, p), loc q)` contains the environment of its unit. -/
lemma elem_unit {t : Ty} {y : UnitF NU → T t}
    (g : Has m (fun u : UnitF NU => B u.1) (p s w) t (fun u => (x u.1, u.2.val)) y) :
    Has m (fun e : ElemF NU ib => B e.1.1) (p (p s w) w) t
      (fun e => ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2)) (fun e => y e.1) :=
  Can.first.then_do (g.reindex (fun e : ElemF NU ib => e.1))

section Elem
variable {α ρ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
variable {cl : Paint} {X : Layout α} {R : Layout ρ}

/-- **Every element of a class table is a `GRig`** (column `c p`, no leftover bits). -/
lemma erig {V : ∀ i (p : Fin (NU i)), Bits (ib i p) → Grain α ρ 0 (col i p)}
    (ct : Has m B s (a (a (a (c cl)))) x (fun i => ctab X R (NU i) (col i) (ib i) (V i))) :
    GRig m (fun e : ElemF NU ib => B e.1.1) (fun e => col e.1.1 e.1.2) (p (p s w) w)
      (fun e => ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2)) cl X R 0
      (fun e => V e.1.1 e.1.2 e.2) := by
  have hpc := h.elem_unit (ib:=ib) (unit_fetch _ h.tpc)
  have hc := h.elem_unit (ib:=ib) (unit_fetch _ h.tc)
  have hxr := h.elem_unit (ib:=ib) (unit_fetch _ h.txr)
  have hlm := h.elem_unit (ib:=ib) (unit_fetch _ h.tlum)
  have hbd : Small (fun e : ElemF NU ib => B e.1.1) (fun e => 2^(col e.1.1 e.1.2)) :=
    h.uBig.reindex (fun e : ElemF NU ib => e.1)
  -- the unit's tape of instances, then the instance
  have hu : Has m (fun u : UnitF NU => B u.1) (p s w) (a (a (c cl))) (fun u => (x u.1, u.2.val))
      (fun u => (bits (ib u.1 u.2)).tape (fun q => (fleet X R 0 (col u.1 u.2)).tape (V u.1 u.2 q))) :=
    unit_fetch (m:=m) (B:=B) (s:=s) (x:=x) (NU:=NU) (t:=a (a (c cl)))
      (fun i (p : Fin (NU i)) => (bits (ib i p)).tape (fun q => (fleet X R 0 (col i p)).tape (V i p q))) ct
  have hit : Has m (fun e : ElemF NU ib => B e.1.1) (p (p s w) w) (a (c cl))
      (fun e => ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2))
      (fun e => (fleet X R 0 (col e.1.1 e.1.2)).tape (V e.1.1 e.1.2 e.2)) :=
    Can.fetch (t:=a (c cl)) (L:=fun e : ElemF NU ib => bits (ib e.1.1 e.1.2))
      (v:=fun e q => (fleet X R 0 (col e.1.1 e.1.2)).tape (V e.1.1 e.1.2 q)) (i:=fun e => e.2)
      (h.elem_unit hu) Can.second
  have him : Has m (fun e : ElemF NU ib => B e.1.1) (p (p s w) w) sc
      (fun e => ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2)) (fun _ => Complex.I) :=
    (Can.first.then_do Can.first).then_do (h.imag.reindex (fun e : ElemF NU ib => e.1.1))
  exact ⟨⟨⟨hpc, hc, hxr, hbd⟩, him, hit⟩, hlm⟩

/-- **Rebuilding a class table** from per-element role tables: linear work. -/
lemma ctab_build {V' : ∀ i (p : Fin (NU i)), Bits (ib i p) → Grain α ρ 0 (col i p)}
    (hel : Can m (fun e : ElemF NU ib => B e.1.1) (p (p s w) w) (a (c cl))
      (fun e => ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2))
      (fun e => (fleet X R 0 (col e.1.1 e.1.2)).tape (V' e.1.1 e.1.2 e.2))
      (fun e => R.n * 2^(X.n * col e.1.1 e.1.2))) (hX : X.n = mx) :
    Can m B s (a (a (a (c cl)))) x (fun i => ctab X R (NU i) (col i) (ib i) (V' i))
      (fun i => NU i * 2^(n i)) := by
  -- inner layout: one unit's tape of instances
  have hlen : Knows m (fun u : UnitF NU => B u.1) (p s w) (fun u => (x u.1, u.2.val))
      (fun u => (bits (ib u.1 u.2)).n) := unit_fetch _ h.tpib
  have HI := Can.layout (t:=a (c cl)) (L:=fun u : UnitF NU => bits (ib u.1 u.2))
    (v:=fun u q => (fleet X R 0 (col u.1 u.2)).tape (V' u.1 u.2 q))
    (P:=fun u => R.n * 2^(X.n * col u.1 u.2)) (W:=fun _ => 0) hlen hel h.uBigIb
  have HI' : Can m (fun u : UnitF NU => B u.1) (p s w) (a (a (c cl))) (fun u => (x u.1, u.2.val))
      (fun u => (bits (ib u.1 u.2)).tape (fun q => (fleet X R 0 (col u.1 u.2)).tape (V' u.1 u.2 q)))
      (fun u => (R.n + 1) * 2^(n u.1)) := by
    refine HI.mono (fun u => ?_)
    have hf := h.fit u.1 u.2
    have e1 : 2^(ib u.1 u.2) * (R.n * 2^(X.n * col u.1 u.2)) = R.n * 2^(n u.1) := by
      rw [← hf, hX, pow_add]; ring
    have e2 := h.two_ib_le u.1 u.2
    change 0 + 2^(ib u.1 u.2) + 2^(ib u.1 u.2) * (R.n * 2^(X.n * col u.1 u.2)) ≤ _
    rw [e1]; nlinarith
  have HO := Can.layout (t:=a (a (c cl))) (L:=fun i => Layout.std (NU i))
    (v:=fun i p => (bits (ib i p)).tape (fun q => (fleet X R 0 (col i p)).tape (V' i p q)))
    (P:=fun i => (R.n + 1) * 2^(n i)) (W:=fun _ => 0) h.nuK HI' h.nuS
  refine HO.weaken ⟨R.n + 3, fun i => ?_⟩
  change 0 + NU i + NU i * ((R.n + 1) * 2^(n i)) + 1 ≤ (R.n + 3) * (NU i * 2^(n i) + 1)
  have h1 : 1 ≤ 2^(n i) := Nat.one_le_two_pow
  have h2 : NU i ≤ NU i * 2^(n i) := Nat.le_mul_of_pos_right _ h1
  nlinarith

/-- **A free letter is one pass over the class table.** -/
lemma free_letter {V : ∀ i (p : Fin (NU i)), Bits (ib i p) → Grain α ρ 0 (col i p)}
    (ct : Has m B s (a (a (a (c cl)))) x (fun i => ctab X R (NU i) (col i) (ib i) (V i)))
    (hX : X.n = mx) (mv : GMove α ρ) (hfree : mv.rk? = none) :
    Can m B s (a (a (a (c cl)))) x
      (fun i => ctab X R (NU i) (col i) (ib i) (fun p q => gferry 0 (col i p) mv (V i p q)))
      (fun i => NU i * 2^(n i)) := by
  have G := h.erig ct
  apply h.ctab_build _ hX
  have hn : ∀ e : ElemF NU ib, R.n * 2^(X.n * col e.1.1 e.1.2) = (fleet X R 0 (col e.1.1 e.1.2)).n := by
    intro e; simp [fleet]
  have key : Can m (fun e : ElemF NU ib => B e.1.1) (p (p s w) w) (a (c cl))
      (fun e => ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2))
      (fun e => (fleet X R 0 (col e.1.1 e.1.2)).tape (gferry 0 (col e.1.1 e.1.2) mv (V e.1.1 e.1.2 e.2)))
      (fun e => (fleet X R 0 (col e.1.1 e.1.2)).n) := by
    cases mv with
    | old mv =>
      cases mv with
      | localM M => exact G.rig.point_cost M
      | shift l z => exact G.rig.shift_cost l z
      | block l rk z ind => simp [GMove.rk?, GMove.shadow, BMove.rk?] at hfree
    | perm l Gp => exact G.rig.perm_cost l Gp
    | phase l z c' => exact G.phase_cost l z c'
  exact key.mono (fun e => le_of_eq (hn e).symm)

end Elem
end CCtx

end
end PowerSaving.RAM
end OAI
