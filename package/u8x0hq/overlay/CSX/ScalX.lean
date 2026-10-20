import Work.XMut.ScalarZ
import CSX.Scat

/-!
# fsp-cc (14 Oct 2026): the kernel scalar replay with split centres on SOURCES: the scatter step

Sussman's kernel scalar check (`Work/GCert/Scalar/*`, reused by x-mutable's `ScalarZ`) replays the scatter with the
kind `okSK` (a y role receives from a SLOT) and requires every retained register to be a slot (`retOKK`, the field
`PGood.hret`).  Here, for centres on mutable sources:
* `okSXK`: a y role receives from a SOURCE or a slot; `retOKXK` / `gPreX` / `PGoodX`: every retained register is a
  source or a slot (the parameter tests of `gPre` otherwise unchanged);
* `scatRowX_run`, `scatFromX_run`: Sussman's `scatRow_run` / `scatFrom_run` with the reference state agreeing on the
  sources AND the slots (the scatter writes y roles only);
* `Cxm` (copies += sources) next to Sussman's `Ccm` (copies += slots; zero for a centre on a source), and
  **`cm_scatX`**: the scatter adds act on the content matrix as `scatM (Jrm Ccm) * scatX (Jrm Cxm)`, the source part
  first (`CSX.scatX`), exactly the factorisation `S.MA = scatX (Jr Cx) * MA0` that `CSX.GLabXZS` consumes.
For a certificate with every centre on a slot, `Cxm = 0` and `cm_scatX` is Sussman's `cm_scat`.
No `sorry`.  Credits: Jacob Sussman (kernel scalar check, replay lemmas; Apache-2.0); x-mutable (ScalarZ).
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option linter.unusedSimpArgs false
set_option linter.deprecated false

namespace GS
open SSC GXD OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.RAM OAI.PowerSaving.SS
  OAI.PowerSaving.CB OAI.PowerSaving.CR Finset Matrix

/-- kinds of the scatter, centres on sources: a y role receives from a source or a slot -/
def okSXK (v v2 nr t s : Nat) : Bool :=
  Nat.ble v t && (Nat.ble (Nat.succ t) v2 &&
    (Nat.ble (Nat.succ s) v || (Nat.ble v2 s && Nat.ble (Nat.succ s) nr)))

theorem okSX_spec {v v2 nr t s : Nat} (h : okSXK v v2 nr t s = true) :
    v ≤ t ∧ t < v2 ∧ (s < v ∨ (v2 ≤ s ∧ s < nr)) := by
  simp only [okSXK, Bool.and_eq_true, Bool.or_eq_true, Nat.ble_eq] at h
  omega

abbrev okSX (c : Raw) : Nat → Nat → Bool := okSXK c.v (2 * c.v) (2 * c.v + c.R)

/-- every retained register is a source or a slot -/
noncomputable def retOKXK (v v2 nr : Nat) : List (Nat × Nat) → Bool :=
  List.rec (motive := fun _ => Bool) true
    (fun x _ ih => guard (Nat.ble (Nat.succ x.1) v || (Nat.ble v2 x.1 && Nat.ble (Nat.succ x.1) nr)) ih)

theorem retOKXK_sound (v v2 nr : Nat) : ∀ (l : List (Nat × Nat)), retOKXK v v2 nr l = true →
    ∀ x ∈ l, x.1 < v ∨ (v2 ≤ x.1 ∧ x.1 < nr) := by
  intro l
  induction l with
  | nil => intro _ x hx; exact absurd hx List.not_mem_nil
  | cons y l ih =>
    intro h x hx
    have h' : guard (Nat.ble (Nat.succ y.1) v || (Nat.ble v2 y.1 && Nat.ble (Nat.succ y.1) nr))
        (retOKXK v v2 nr l) = true := h
    obtain ⟨h1, h2⟩ := guard_true h'
    rcases List.mem_cons.mp hx with e | e
    · rw [e]
      simp only [Bool.or_eq_true, Bool.and_eq_true, Nat.ble_eq] at h1
      omega
    · exact ih h2 x e

/-- Sussman's `gPre` with the retained-register test `retOKXK` (sources or slots) -/
noncomputable def gPreX (c : Raw) (p : SPar) (lo n : Nat)
    (k : Nat → Nat → Nat → Nat → Nat → Nat → Bool) : Bool :=
  forceN (Nat.mul 2 c.v) fun v2 => forceN (Nat.add v2 c.R) fun nr => forceN (Nat.add lo n) fun hi =>
    forceN (Nat.pow 2 (Nat.sub p.sw 1)) fun half => forceN (Nat.add 64 (Nat.mul p.sw n)) fun K =>
      guard (Nat.ble 1 p.ux) (guard (Nat.ble 1 p.us) (guard (Nat.ble 1 p.uy)
        (guard (Nat.ble 1 p.sw) (guard (Nat.ble nr (Nat.pow 2 p.d)) (guard (Nat.ble hi c.v)
          (guard (Nat.ble (Nat.succ p.ux) half) (guard (Nat.ble (Nat.succ p.uy) half)
            (guard (Nat.ble (Nat.succ p.ux) 2147483648) (guard (Nat.ble (Nat.succ p.uy) 2147483648)
              (guard (retOKXK c.v v2 nr c.ret)
                (forceN (Nat.lor (Nat.shiftLeft (Nat.shiftLeft (rep p.sw n) (Nat.sub p.sw 1)) 64)
                    (Nat.pow 2 63)) fun tm =>
                  forceN (Nat.pow 2 K) fun m2 => k v2 nr hi tm K m2)))))))))))

theorem gPreX_sound (c : Raw) (p : SPar) (lo n : Nat) (k : Nat → Nat → Nat → Nat → Nat → Nat → Bool)
    (h : gPreX c p lo n k = true) :
    (1 ≤ p.ux ∧ 1 ≤ p.us ∧ 1 ≤ p.uy ∧ 1 ≤ p.sw) ∧ 2 * c.v + c.R ≤ 2^p.d ∧ lo + n ≤ c.v ∧
      (p.ux < 2^(p.sw - 1) ∧ p.uy < 2^(p.sw - 1) ∧ p.ux < 2^31 ∧ p.uy < 2^31) ∧
      (∀ x ∈ c.ret, x.1 < c.v ∨ (2 * c.v ≤ x.1 ∧ x.1 < 2 * c.v + c.R)) ∧
      k (2 * c.v) (2 * c.v + c.R) (lo + n) (maskW p.sw n) (64 + p.sw * n) (2^(64 + p.sw * n))
        = true := by
  unfold gPreX at h
  rw [forceN_eq, forceN_eq, forceN_eq, forceN_eq, forceN_eq] at h
  obtain ⟨h1, h⟩ := guard_true h
  obtain ⟨h2, h⟩ := guard_true h
  obtain ⟨h3, h⟩ := guard_true h
  obtain ⟨h4, h⟩ := guard_true h
  obtain ⟨h5, h⟩ := guard_true h
  obtain ⟨h6, h⟩ := guard_true h
  obtain ⟨h7, h⟩ := guard_true h
  obtain ⟨h8, h⟩ := guard_true h
  obtain ⟨h9, h⟩ := guard_true h
  obtain ⟨h10, h⟩ := guard_true h
  obtain ⟨h11, h⟩ := guard_true h
  rw [forceN_eq, forceN_eq] at h
  exact ⟨⟨ble_true h1, ble_true h2, ble_true h3, ble_true h4⟩, ble_true h5, ble_true h6,
    ⟨ble_true h7, ble_true h8, ble_true h9, ble_true h10⟩, retOKXK_sound _ _ _ _ h11, h⟩

/-- the tests on the parameters, retained registers sources or slots -/
structure PGoodX (c : Raw) (p : SPar) (lo n : Nat) : Prop where
  hx : 1 ≤ p.ux
  hs : 1 ≤ p.us
  hy : 1 ≤ p.uy
  hsw : 1 ≤ p.sw
  hd : 2 * c.v + c.R ≤ 2^p.d
  hn : lo + n ≤ c.v
  hxs : p.ux < 2^(p.sw - 1)
  hys : p.uy < 2^(p.sw - 1)
  hx2 : p.ux < 2^31
  hy2 : p.uy < 2^31
  hret : ∀ x ∈ c.ret, x.1 < c.v ∨ (2 * c.v ≤ x.1 ∧ x.1 < 2 * c.v + c.R)

theorem pgoodX_of {c : Raw} {p : SPar} {lo n : Nat} {X : Prop}
    (h : (1 ≤ p.ux ∧ 1 ≤ p.us ∧ 1 ≤ p.uy ∧ 1 ≤ p.sw) ∧ 2 * c.v + c.R ≤ 2^p.d ∧ lo + n ≤ c.v ∧
      (p.ux < 2^(p.sw - 1) ∧ p.uy < 2^(p.sw - 1) ∧ p.ux < 2^31 ∧ p.uy < 2^31) ∧
      (∀ x ∈ c.ret, x.1 < c.v ∨ (2 * c.v ≤ x.1 ∧ x.1 < 2 * c.v + c.R)) ∧ X) :
    PGoodX c p lo n ∧ X :=
  ⟨⟨h.1.1, h.1.2.1, h.1.2.2.1, h.1.2.2.2, h.2.1, h.2.2.1, h.2.2.2.1.1, h.2.2.2.1.2.1,
    h.2.2.2.1.2.2.1, h.2.2.2.1.2.2.2, h.2.2.2.2.1⟩, h.2.2.2.2.2⟩

/-! ## the scatter rows, sources allowed -/

section
variable (sw n : Nat) (u : Nat → Nat) (okf : Nat → Nat → Bool) (v1 v2 : Nat)
  (hok : ∀ t s, okf t s = true → (s < v1 ∨ v2 ≤ s) ∧ v1 ≤ t ∧ t < v2) (ret : List (Nat × Nat))
include hok

/-- one scatter row (Sussman's `scatRow_run`): the reference state `a0` agrees with the current one on the sources
and the slots (the registers outside `[v1, v2)`). -/
theorem scatRowX_run (t : Nat) (l : List (Nat × Co)) : ∀ (a b a0 : SSt),
    (∀ r, (r < v1 ∨ v2 ≤ r) → a.P r = a0.P r ∧ a.N r = a0.N r) →
    GRun sw n u okf (fun _ => true) a (scatRowAdds ret t l) b →
    (∀ r, r ≠ t → b.P r = a.P r ∧ b.N r = a.N r) ∧
    ∀ T, T < n → gcont sw u b t T = gcont sw u a t T
      + (l.map fun x => (toCoef x.2).val * gcont sw u a0 (ret.getD x.1 (0, 0)).1 T).sum := by
  induction l with
  | nil =>
    intro a b a0 _ h
    have e : b = a := h
    rw [e]
    exact ⟨fun r _ => ⟨rfl, rfl⟩, fun T _ => by simp⟩
  | cons x l ih =>
    intro a b a0 h0 h
    have h' : if (fun _ => true) (ret.getD x.1 (0, 0)).1 = true then
        (okf t (ret.getD x.1 (0, 0)).1 = true ∧ ∃ m, GStep sw n u t (ret.getD x.1 (0, 0)).1
          (toCoef x.2) a m ∧ GRun sw n u okf (fun _ => true) m (scatRowAdds ret t l) b)
        else GRun sw n u okf (fun _ => true) a (scatRowAdds ret t l) b := h
    rw [if_pos rfl] at h'
    obtain ⟨ho, m, hstep, hrun⟩ := h'
    obtain ⟨hs2, ht1, ht2⟩ := hok _ _ ho
    obtain ⟨_, hoth, hcont⟩ := hstep
    have hm0 : ∀ r, (r < v1 ∨ v2 ≤ r) → m.P r = a0.P r ∧ m.N r = a0.N r := by
      intro r hr
      obtain ⟨p1, p2⟩ := hoth r (by omega)
      obtain ⟨q1, q2⟩ := h0 r hr
      exact ⟨p1.trans q1, p2.trans q2⟩
    obtain ⟨i1, i2⟩ := ih m b a0 hm0 hrun
    refine ⟨fun r hr => ?_, fun T hT => ?_⟩
    · obtain ⟨p1, p2⟩ := hoth r hr
      obtain ⟨q1, q2⟩ := i1 r hr
      exact ⟨q1.trans p1, q2.trans p2⟩
    · obtain ⟨q1, q2⟩ := h0 _ hs2
      rw [i2 T hT, hcont T hT, gcont_congr sw u a0 a _ T q1 q2, List.map_cons, List.sum_cons]
      ring

/-- all scatter rows (Sussman's `scatFrom_run`), sources allowed. -/
theorem scatFromX_run (rows : List (List (Nat × Co))) :
    ∀ (t0 : Nat) (a b a0 : SSt), v1 ≤ t0 → t0 + rows.length ≤ v2 →
    (∀ r, (r < v1 ∨ v2 ≤ r) → a.P r = a0.P r ∧ a.N r = a0.N r) →
    GRun sw n u okf (fun _ => true) a (scatAddsFrom ret t0 rows) b →
    (∀ r, (r < t0 ∨ t0 + rows.length ≤ r) → b.P r = a.P r ∧ b.N r = a.N r) ∧
    ∀ j, j < rows.length → ∀ T, T < n → gcont sw u b (t0 + j) T = gcont sw u a (t0 + j) T
      + ((rows.getD j []).map fun x =>
          (toCoef x.2).val * gcont sw u a0 (ret.getD x.1 (0, 0)).1 T).sum := by
  induction rows with
  | nil =>
    intro t0 a b a0 _ _ _ h
    have e : b = a := h
    rw [e]
    exact ⟨fun r _ => ⟨rfl, rfl⟩, fun j hj => absurd hj (by simp)⟩
  | cons l rows ih =>
    intro t0 a b a0 hv1 ht0 h0 h
    have hlen : (l :: rows).length = rows.length + 1 := rfl
    rw [hlen] at ht0
    obtain ⟨m, h1, h2⟩ := GRun_split sw n u okf _ (scatRowAdds ret t0 l)
      (scatAddsFrom ret (t0+1) rows) a b h
    obtain ⟨r1, r2⟩ := scatRowX_run sw n u okf v1 v2 hok ret t0 l a m a0 h0 h1
    have hm0 : ∀ r, (r < v1 ∨ v2 ≤ r) → m.P r = a0.P r ∧ m.N r = a0.N r := by
      intro r hr
      obtain ⟨p1, p2⟩ := r1 r (by omega)
      obtain ⟨q1, q2⟩ := h0 r hr
      exact ⟨p1.trans q1, p2.trans q2⟩
    obtain ⟨i1, i2⟩ := ih (t0+1) m b a0 (by omega) (by omega) hm0 h2
    refine ⟨fun r hr => ?_, fun j hj T hT => ?_⟩
    · rw [hlen] at hr
      obtain ⟨p1, p2⟩ := r1 r (by omega)
      obtain ⟨q1, q2⟩ := i1 r (by omega)
      exact ⟨q1.trans p1, q2.trans p2⟩
    · cases j with
      | zero =>
        obtain ⟨q1, q2⟩ := i1 t0 (Or.inl (by omega))
        show gcont sw u b (t0 + 0) T = gcont sw u a (t0 + 0) T
          + (l.map fun x => (toCoef x.2).val * gcont sw u a0 (ret.getD x.1 (0, 0)).1 T).sum
        rw [Nat.add_zero, gcont_congr sw u m b t0 T q1 q2]
        exact r2 T hT
      | succ j =>
        have hj' : j < rows.length := by rw [hlen] at hj; omega
        have e : t0 + (j + 1) = t0 + 1 + j := by omega
        obtain ⟨p1, p2⟩ := r1 (t0 + 1 + j) (by omega)
        show gcont sw u b (t0 + (j + 1)) T = gcont sw u a (t0 + (j + 1)) T
          + ((rows.getD j []).map fun x =>
              (toCoef x.2).val * gcont sw u a0 (ret.getD x.1 (0, 0)).1 T).sum
        rw [e, i2 j hj' T hT, gcont_congr sw u a m (t0 + 1 + j) T p1 p2]

end

/-! ## the scatter on the content matrix -/

/-- copies += sources: copy `k` reads the SOURCE register of the retained total `k` (zero row for a slot) -/
def Cxm (c : Raw) : Matrix (Fin c.ret.length) (Fin c.v) ℚ :=
  fun k t => if t.val = (c.ret.getD k.val (0, 0)).1 then 1 else 0

variable {c : Raw}

/-- **the scatter adds act on the content matrix as `scatM (Jrm Ccm) * scatX (Jrm Cxm)`** (centres on sources or
slots; Sussman's `cm_scat` when every centre is on a slot). -/
theorem cm_scatX (p : SPar) (lo n : Nat) (hp : PGoodX c p lo n) (a1 a2 : SSt)
    (h : GRun p.sw n (U c p) (okSX c) (fun _ => true) a1 (scatAdds c) a2) (hsh : ScatShape c) :
    Cm p.sw n (U c p) (emb c) a2
      = scatM (Jrm c * Ccm c) * (CSX.scatX (Sl := Fin c.R) (Jrm c * Cxm c) * Cm p.sw n (U c p) (emb c) a1) := by
  obtain ⟨hlen, hent⟩ := hsh
  have hl : (scatRows c).length = c.v := by omega
  obtain ⟨u1, u2⟩ := scatFromX_run p.sw n (U c p) (okSX c) c.v (2 * c.v)
    (fun t s ho => by
      obtain ⟨h1, h2, h3⟩ := okSX_spec ho
      exact ⟨by omega, h1, h2⟩) c.ret (scatRows c) c.v a1 a2 a1
    le_rfl (by omega) (fun r _ => ⟨rfl, rfl⟩) h
  ext i T
  have eY : ∀ j : Rl c, (CSX.scatX (Sl := Fin c.R) (Jrm c * Cxm c) * Cm p.sw n (U c p) (emb c) a1) j T
      = Cm p.sw n (U c p) (emb c) a1 j T + Sum.elim (Sum.elim (fun _ => 0)
          (fun S => ∑ t, (Jrm c * Cxm c) S t * Cm p.sw n (U c p) (emb c) a1 (Sum.inl (Sum.inl t)) T))
          (fun _ => 0) j :=
    fun j => CSX.scatX_mul_apply _ _ j T
  unfold scatM
  rw [Matrix.add_mul, Matrix.one_mul, Matrix.add_apply, embYS_mul_apply]
  rcases i with (t' | S') | q
  · obtain ⟨p1, p2⟩ := u1 t'.val (Or.inl t'.isLt)
    simp only [eY, Sum.elim_inl, Sum.elim_inr, add_zero]
    show gcont p.sw (U c p) a2 t'.val T.val = gcont p.sw (U c p) a1 t'.val T.val
    exact gcont_congr p.sw (U c p) a1 a2 _ _ p1 p2
  · have hS : S'.val < (scatRows c).length := by rw [hl]; exact S'.isLt
    have hrow : ∀ x ∈ (scatRows c).getD S'.val [], x.1 < c.ret.length :=
      hent _ (getD_mem_lt (scatRows c) [] S'.val hS)
    -- the register of a retained total, read through the two copy matrices
    have hg : ∀ k : Fin c.ret.length,
        (∑ q, Ccm c k q * Cm p.sw n (U c p) (emb c) a1 (Sum.inr q) T)
          + (∑ t, Cxm c k t * Cm p.sw n (U c p) (emb c) a1 (Sum.inl (Sum.inl t)) T)
          = gcont p.sw (U c p) a1 (c.ret.getD k.val (0, 0)).1 T.val := by
      intro k
      rcases hp.hret _ (getD_mem_lt c.ret (0, 0) k.val k.isLt) with hx | ⟨r1, r2⟩
      · have h0 : ∑ q, Ccm c k q * Cm p.sw n (U c p) (emb c) a1 (Sum.inr q) T = 0 := by
          apply Finset.sum_eq_zero
          intro q _
          have hne : ¬ 2 * c.v + q.val = (c.ret.getD k.val (0, 0)).1 := by omega
          show (if 2 * c.v + q.val = (c.ret.getD k.val (0, 0)).1 then (1:ℚ) else 0) * _ = 0
          rw [if_neg hne, zero_mul]
        rw [h0, zero_add, Finset.sum_eq_single (⟨(c.ret.getD k.val (0, 0)).1, hx⟩ : Fin c.v)]
        · show (if (c.ret.getD k.val (0, 0)).1 = (c.ret.getD k.val (0, 0)).1 then (1:ℚ) else 0)
            * gcont p.sw (U c p) a1 (c.ret.getD k.val (0, 0)).1 T.val = _
          rw [if_pos rfl, one_mul]
        · intro t _ ht
          have hne : ¬ t.val = (c.ret.getD k.val (0, 0)).1 := fun e => ht (Fin.ext e)
          show (if t.val = (c.ret.getD k.val (0, 0)).1 then (1:ℚ) else 0) * _ = 0
          rw [if_neg hne, zero_mul]
        · intro hb
          exact absurd (Finset.mem_univ _) hb
      · have h0 : ∑ t, Cxm c k t * Cm p.sw n (U c p) (emb c) a1 (Sum.inl (Sum.inl t)) T = 0 := by
          apply Finset.sum_eq_zero
          intro t _
          have hne : ¬ t.val = (c.ret.getD k.val (0, 0)).1 := by have := t.isLt; omega
          show (if t.val = (c.ret.getD k.val (0, 0)).1 then (1:ℚ) else 0) * _ = 0
          rw [if_neg hne, zero_mul]
        rw [h0, add_zero, Finset.sum_eq_single (⟨(c.ret.getD k.val (0, 0)).1 - 2 * c.v, by omega⟩ : Fin c.R)]
        · have e : 2 * c.v + ((c.ret.getD k.val (0, 0)).1 - 2 * c.v) = (c.ret.getD k.val (0, 0)).1 := by
            omega
          show (if 2 * c.v + ((c.ret.getD k.val (0, 0)).1 - 2 * c.v) = (c.ret.getD k.val (0, 0)).1
            then (1:ℚ) else 0) * gcont p.sw (U c p) a1
              (2 * c.v + ((c.ret.getD k.val (0, 0)).1 - 2 * c.v)) T.val = _
          rw [if_pos e, one_mul, e]
        · intro q _ hq
          have hne : ¬ 2 * c.v + q.val = (c.ret.getD k.val (0, 0)).1 := fun e =>
            hq (Fin.ext (by show q.val = (c.ret.getD k.val (0, 0)).1 - 2 * c.v; omega))
          show (if 2 * c.v + q.val = (c.ret.getD k.val (0, 0)).1 then (1:ℚ) else 0) * _ = 0
          rw [if_neg hne, zero_mul]
        · intro hb
          exact absurd (Finset.mem_univ _) hb
    have step : (∑ q, (Jrm c * Ccm c) S' q * Cm p.sw n (U c p) (emb c) a1 (Sum.inr q) T)
        + (∑ t, (Jrm c * Cxm c) S' t * Cm p.sw n (U c p) (emb c) a1 (Sum.inl (Sum.inl t)) T)
        = ∑ k, Jrm c S' k * gcont p.sw (U c p) a1 (c.ret.getD k.val (0, 0)).1 T.val := by
      have e1 : ∑ q, (Jrm c * Ccm c) S' q * Cm p.sw n (U c p) (emb c) a1 (Sum.inr q) T
          = ∑ k, Jrm c S' k * ∑ q, Ccm c k q * Cm p.sw n (U c p) (emb c) a1 (Sum.inr q) T := by
        simp only [Matrix.mul_apply, Finset.sum_mul, Finset.mul_sum, mul_assoc]
        exact Finset.sum_comm
      have e2 : ∑ t, (Jrm c * Cxm c) S' t * Cm p.sw n (U c p) (emb c) a1 (Sum.inl (Sum.inl t)) T
          = ∑ k, Jrm c S' k * ∑ t, Cxm c k t * Cm p.sw n (U c p) (emb c) a1 (Sum.inl (Sum.inl t)) T := by
        simp only [Matrix.mul_apply, Finset.sum_mul, Finset.mul_sum, mul_assoc]
        exact Finset.sum_comm
      rw [e1, e2, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl (fun k _ => ?_)
      rw [← mul_add, hg k]
    have hrowsum : (((scatRows c).getD S'.val []).map fun x =>
          (toCoef x.2).val * gcont p.sw (U c p) a1 (c.ret.getD x.1 (0, 0)).1 T.val).sum
        = ∑ k, Jrm c S' k * gcont p.sw (U c p) a1 (c.ret.getD k.val (0, 0)).1 T.val := by
      have hl' : ∀ x ∈ ((scatRows c).getD S'.val []).map (fun x => (x.1, toCoef x.2)),
          x.1 < c.ret.length := by
        intro x hx
        obtain ⟨y, hy, rfl⟩ := List.mem_map.mp hx
        exact hrow y hy
      have := rowOf_sum (((scatRows c).getD S'.val []).map fun x => (x.1, toCoef x.2)) hl'
        (fun k => gcont p.sw (U c p) a1 (c.ret.getD k.val (0, 0)).1 T.val)
        (fun k' => gcont p.sw (U c p) a1 (c.ret.getD k' (0, 0)).1 T.val) (fun _ => rfl)
      rw [List.map_map] at this
      exact this.symm
    have e2 : Cm p.sw n (U c p) (emb c) a2 (Sum.inl (Sum.inr S')) T
        = Cm p.sw n (U c p) (emb c) a1 (Sum.inl (Sum.inr S')) T
          + ∑ k, Jrm c S' k * gcont p.sw (U c p) a1 (c.ret.getD k.val (0, 0)).1 T.val := by
      show gcont p.sw (U c p) a2 (c.v + S'.val) T.val = gcont p.sw (U c p) a1 (c.v + S'.val) T.val + _
      rw [u2 S'.val hS T.val T.isLt, hrowsum]
    rw [e2, ← step]
    simp only [eY, Sum.elim_inl, Sum.elim_inr, add_zero]
    abel
  · obtain ⟨p1, p2⟩ := u1 (2 * c.v + q.val) (Or.inr (by omega))
    simp only [eY, Sum.elim_inl, Sum.elim_inr, add_zero]
    show gcont p.sw (U c p) a2 (2 * c.v + q.val) T.val = gcont p.sw (U c p) a1 (2 * c.v + q.val) T.val
    exact gcont_congr p.sw (U c p) a1 a2 _ _ p1 p2

end GS

#print axioms GS.gPreX_sound
#print axioms GS.scatFromX_run
#print axioms GS.cm_scatX
