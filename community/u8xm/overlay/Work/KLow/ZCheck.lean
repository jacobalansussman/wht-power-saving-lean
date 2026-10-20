import Work.Gauge.ECheck

/-!
# (key: k-lower) The zero check: entries of the total matrix that vanish, by one packed replay

k-lower, 13 Oct 2026.  A twin of DFT Bridge3's E7 check `GS.eCheck` (Work/Gauge/ECheck.lean, Ryan Shea's
lean/complex588), itself a mirror of Jacob Sussman's `GS.scalCheck` (Work/GCert/Scalar, Apache-2.0).

`zCheck c p gl rows`: the sources are the registers of the list `gl` (strictly increasing, any registers: a y register
is a harmless filler that shifts the digits of the others); every other register starts empty; the adds of `A`, the scatter and the adds of `B` are replayed with
the digit rule of `gAdd`; at the end every register of `rows` that is NOT a source holds nothing.
Soundness `zCheck_sound`: for every role `j` whose register is in `gl` and every role `i` whose register is in
`rows \ gl`, `Mt i j = 0` (`Mt = MB * scatM (Jrm * Ccm) * MA`).  The replay, its invariant and the content matrix
are Sussman's / Bridge3's unchanged (`gInitL`, `gItems`, `gItems_sound`, `init_okL`, `cm_scat`, `GRun_Cm`); only the
final test (`zFin`) is new.  A source of unit `U r` holding `us` gives the column
`(us / U r) e_j`, a nonzero multiple of the unit column: enough for a ZERO entry.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option linter.unusedSimpArgs false
set_option linter.deprecated false

namespace GS
open SSC GXD OAI.PowerSaving OAI.PowerSaving.Binary OAI.PowerSaving.RAM OAI.PowerSaving.SS
  OAI.PowerSaving.CB OAI.PowerSaving.CR Finset Matrix

/-- `r ∈ gl`, as a kernel term -/
noncomputable def memZ (r : Nat) (gl : List Nat) : Bool :=
  List.rec (motive := fun _ => Bool) false
    (fun x _ ih => Bool.rec (motive := fun _ => Bool) ih true (Nat.beq x r)) gl

theorem memZ_sound (r : Nat) : ∀ gl : List Nat, memZ r gl = true → r ∈ gl := by
  intro gl
  induction gl with
  | nil => intro h; exact absurd h Bool.false_ne_true
  | cons x gl ih =>
    intro h
    have h' : Bool.rec (motive := fun _ => Bool) (memZ r gl) true (Nat.beq x r) = true := h
    cases e : Nat.beq x r with
    | true => rw [← beq_true e]; exact List.mem_cons_self
    | false => rw [e] at h'; exact List.mem_cons_of_mem _ (ih h')

/-- no source is a y register: `r < v` or `2 v ≤ r` -/
noncomputable def noYK (v v2 : Nat) : List Nat → Bool :=
  List.rec (motive := fun _ => Bool) true
    (fun r _ ih => guard (Bool.rec (motive := fun _ => Bool) (Nat.ble v2 r) true (Nat.ble (Nat.succ r) v)) ih)

theorem noYK_sound (v v2 : Nat) : ∀ gl : List Nat, noYK v v2 gl = true → ∀ r ∈ gl, r < v ∨ v2 ≤ r := by
  intro gl
  induction gl with
  | nil => intro _ r hr; exact absurd hr List.not_mem_nil
  | cons x gl ih =>
    intro h r hr
    have h' : guard (Bool.rec (motive := fun _ => Bool) (Nat.ble v2 x) true (Nat.ble (Nat.succ x) v))
      (noYK v v2 gl) = true := h
    obtain ⟨h1, h2⟩ := guard_true h'
    rcases List.mem_cons.mp hr with e | e
    · subst e
      cases e2 : Nat.ble (Nat.succ r) v with
      | true => exact Or.inl (ble_true e2)
      | false => rw [e2] at h1; exact Or.inr (ble_true h1)
    · exact ih h2 r e

/-- at the end every register of `rows` that is not a source holds nothing (tags dropped) -/
noncomputable def zFin (K m2 : Nat) (gl : List Nat) (tr : Trie) : List Nat → Bool :=
  List.rec (motive := fun _ => Bool) true
    (fun r _ ih => guard (Bool.rec (motive := fun _ => Bool)
        (tr.getK r fun x => Nat.beq (Nat.shiftRight (Nat.mod x m2) 64) (Nat.shiftRight (Nat.shiftRight x K) 64))
        true (memZ r gl)) ih)

theorem zFin_sound (d K : Nat) (gl : List Nat) (tr : Trie) (hwf : tr.wf d) :
    ∀ rows : List Nat, zFin K (2^K) gl tr rows = true → ∀ r ∈ rows, r ∉ gl →
      (absG d K tr).P r = (absG d K tr).N r := by
  intro rows
  induction rows with
  | nil => intro _ r hr; exact absurd hr List.not_mem_nil
  | cons r0 rows ih =>
    intro h r hr hng
    have h' : guard (Bool.rec (motive := fun _ => Bool)
        (tr.getK r0 fun x => Nat.beq (Nat.shiftRight (Nat.mod x (2^K)) 64) (Nat.shiftRight (Nat.shiftRight x K) 64))
        true (memZ r0 gl)) (zFin K (2^K) gl tr rows) = true := h
    obtain ⟨h1, h2⟩ := guard_true h'
    rcases List.mem_cons.mp hr with e | e
    · subst e
      cases hm : memZ r gl with
      | true => exact absurd (memZ_sound r gl hm) hng
      | false =>
        rw [hm] at h1
        have h1' : (tr.getK r fun x => Nat.beq (Nat.shiftRight (Nat.mod x (2^K)) 64)
            (Nat.shiftRight (Nat.shiftRight x K) 64)) = true := h1
        rw [Trie.getK_eq d tr hwf] at h1'
        simp only [Bool.and_eq_true, decide_eq_true_eq] at h1'
        obtain ⟨hr', hb⟩ := h1'
        have hv := beq_true hb
        show (absS (2^d) tr r % 2^K) >>> 64 = (absS (2^d) tr r >>> K) >>> 64
        rw [absS_apply d tr _ hr']
        exact hv
    · exact ih h2 r e hng

/-- **The zero check** (one evaluation). -/
noncomputable def zCheck (c : Raw) (p : SPar) (gl rows : List Nat) : Bool :=
  gPre c p 0 gl.length fun v2 nr _ tm K m2 =>
    guard (incK nr gl 0) (guard (Nat.ble (Nat.succ p.us) (Nat.pow 2 (Nat.sub p.sw 1)))
      (guard (Nat.ble (Nat.succ p.us) 2147483648)
        (forceN (List.length c.ret) fun rl =>
          gInitL p.sw p.us gl 0 (Trie.mk p.d) fun tr0 =>
            gItems c p v2 nr tm K m2 rl (prog c) tr0 fun tr3 =>
              zFin K m2 gl tr3 rows)))

/-- a scaled unit column -/
theorem mul_scaledCol {ρ m : Type} [Fintype ρ] [DecidableEq ρ] (M : Matrix ρ ρ ℚ) (X : Matrix ρ m ℚ) (j0 : ρ)
    (lam : ℚ) (T : m) (hX : ∀ j, X j T = if j = j0 then lam else 0) (i : ρ) : (M * X) i T = M i j0 * lam := by
  rw [Matrix.mul_apply, Finset.sum_eq_single j0]
  · rw [hX j0, if_pos rfl]
  · intro k _ hk
    rw [hX k, if_neg hk, mul_zero]
  · intro hb
    exact absurd (Finset.mem_univ _) hb

/-- **Soundness of the zero check.** -/
theorem zCheck_sound {c : Raw} {p : SPar} {gl rows : List Nat} (hv : 0 < c.v)
    (h : zCheck c p gl rows = true) (i j : Rl c) (hj : emb c j ∈ gl) (hi : emb c i ∈ rows)
    (hni : emb c i ∉ gl) : Mt c hv i j = 0 := by
  unfold zCheck at h
  obtain ⟨hp, h1⟩ := pgood_of (gPre_sound c p 0 gl.length _ h)
  obtain ⟨hinc, h3⟩ := guard_true h1
  obtain ⟨hus, h4⟩ := guard_true h3
  obtain ⟨hus2, h5⟩ := guard_true h4
  rw [forceN_eq] at h5
  obtain ⟨hgl, hnd⟩ := incK_sound _ gl _ hinc
  have hu1 := ble_true hus
  have hu2 := ble_true hus2
  have hlt : ∀ r ∈ gl, r < 2^p.d := fun r hr => by have := (hgl r hr).2; have := hp.hd; omega
  obtain ⟨t0, hwf, hk, hval⟩ := gInitL_sound p.d p.sw p.us gl 0 (Trie.mk p.d) _ (Trie.wf_mk p.d) hlt hnd h5
  have hval' : ∀ r, absS (2^p.d) t0 r = if r ∈ gl then leaf0 p.sw p.us (0 + gl.idxOf r)
      else absS (2^p.d) (Trie.mk p.d) r := hval
  have husw : p.us < 2^(p.sw - 1) := by
    have e : Nat.pow 2 (Nat.sub p.sw 1) = 2^(p.sw - 1) := rfl
    omega
  obtain ⟨hI0, hinit⟩ := init_okL p.d p.sw p.us gl hp.hsw husw (by omega) t0 hwf hval'
  obtain ⟨tr3, hI3, hk3, hr⟩ := gItems_sound c p 0 gl.length hp (prog c) t0 _ hI0 hk
  have hfin := zFin_sound p.d (64 + p.sw * gl.length) gl tr3 hI3.wf rows hk3 (emb c i) hi hni
  obtain ⟨a1, r1, r2⟩ := ItemsRun_split c p gl.length _ _ _ _ hr
  obtain ⟨a2, ⟨rS, hsh⟩, r3⟩ := r2
  have rA := ItemsRun_cA c p gl.length _ _ a1 r1
  have rB := ItemsRun_cB c p gl.length _ a2 _ r3
  have hokA : ∀ t s, okA c t s = true → emb c (unemb c hv t) = t ∧ emb c (unemb c hv s) = s :=
    fun t s ho => ⟨emb_unemb hv t (okA_spec ho).1, emb_unemb hv s (okA_spec ho).2.1⟩
  have hokB : ∀ t s, okB c t s = true → emb c (unemb c hv t) = t ∧ emb c (unemb c hv s) = s :=
    fun t s ho => ⟨emb_unemb hv t (okB_spec ho).1, emb_unemb hv s (okB_spec ho).2.1⟩
  have hCA : Cm p.sw gl.length (U c p) (emb c) a1
      = MA c hv * Cm p.sw gl.length (U c p) (emb c) (absG p.d (64 + p.sw * gl.length) t0) := by
    have := GRun_Cm p.sw gl.length (U c p) (okA c) (fun _ => true) (emb c) (unemb c hv) (unemb_emb hv) hokA
      (addsA c) _ a1 rA
    rw [micF_true] at this
    exact this
  have hCB : Cm p.sw gl.length (U c p) (emb c) (absG p.d (64 + p.sw * gl.length) tr3)
      = MB c hv * Cm p.sw gl.length (U c p) (emb c) a2 := by
    have := GRun_Cm p.sw gl.length (U c p) (okB c) (fun _ => true) (emb c) (unemb c hv) (unemb_emb hv) hokB
      (addsB c) a2 _ rB
    rw [micF_true] at this
    exact this
  have hCS := cm_scat p 0 gl.length hp a1 a2 rS hsh
  have hus' : p.us < 2^p.sw := lt_of_lt_of_le husw (Nat.pow_le_pow_right (by norm_num) (by omega))
  have hT : gl.idxOf (emb c j) < gl.length := List.idxOf_lt_length_of_mem hj
  have hUj : (1 : ℚ) ≤ (U c p (emb c j) : ℚ) := by
    have := unitK_pos p.ux p.us p.uy c.v (2 * c.v) (emb c j) hp.hx hp.hs hp.hy
    exact_mod_cast this
  have hus1 : (1 : ℚ) ≤ (p.us : ℚ) := by exact_mod_cast hp.hs
  have hlam : (p.us : ℚ) / (U c p (emb c j) : ℚ) ≠ 0 := by
    apply div_ne_zero <;> linarith
  have hC0 : ∀ j' : Rl c, Cm p.sw gl.length (U c p) (emb c) (absG p.d (64 + p.sw * gl.length) t0) j'
      ⟨gl.idxOf (emb c j), hT⟩ = if j' = j then (p.us : ℚ) / (U c p (emb c j) : ℚ) else 0 := by
    intro j'
    show gcont p.sw (U c p) (absG p.d (64 + p.sw * gl.length) t0) (emb c j') (gl.idxOf (emb c j)) = _
    unfold gcont
    rw [(hinit (emb c j')).1, (hinit (emb c j')).2, dg_zero_left]
    by_cases hj' : emb c j' ∈ gl
    · rw [if_pos hj', dg_mul_pow p.sw p.us hus']
      by_cases e : j' = j
      · subst e
        rw [if_pos rfl, if_pos rfl]
        simp
      · have e' : ¬ gl.idxOf (emb c j) = gl.idxOf (emb c j') := by
          intro he
          apply e
          have := (List.idxOf_inj hj).mp he
          rw [← unemb_emb hv j', ← this]
          exact unemb_emb hv j
        rw [if_neg e', if_neg e]
        simp
    · have e : ¬ j' = j := fun e => hj' (by rw [e]; exact hj)
      rw [if_neg hj', if_neg e, dg_zero_left]
      simp
  have key : Mt c hv i j * ((p.us : ℚ) / (U c p (emb c j) : ℚ)) = Cm p.sw gl.length (U c p) (emb c)
      (absG p.d (64 + p.sw * gl.length) tr3) i ⟨gl.idxOf (emb c j), hT⟩ := by
    have e1 := mul_scaledCol (Mt c hv) (Cm p.sw gl.length (U c p) (emb c) (absG p.d (64 + p.sw * gl.length) t0))
      j _ _ hC0 i
    rw [← e1]
    unfold Mt
    rw [Matrix.mul_assoc, Matrix.mul_assoc, ← hCA, ← hCS, ← hCB]
  have hz : Cm p.sw gl.length (U c p) (emb c) (absG p.d (64 + p.sw * gl.length) tr3) i
      ⟨gl.idxOf (emb c j), hT⟩ = 0 := by
    show gcont p.sw (U c p) (absG p.d (64 + p.sw * gl.length) tr3) (emb c i) (gl.idxOf (emb c j)) = 0
    unfold gcont
    rw [hfin]
    simp
  rw [hz] at key
  rcases mul_eq_zero.mp key with h0 | h0
  · exact h0
  · exact absurd h0 hlam

end GS

#print axioms GS.zCheck_sound
