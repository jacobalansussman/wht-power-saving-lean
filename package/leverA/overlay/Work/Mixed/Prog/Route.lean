import Work.Mixed.Prog.CTab

/-!
# Mixed engine, R3 program part 3: routing a block letter to the child class (agent e-r3prog)

A block letter of rank `k` on role `l` turns every unit of the class into ONE job of `2^n`
entries (`blockJob`: the role read through `placeS`, block coordinates low).  The child class
of the lockstep tree receives, for each of its input slots, either such a job of a parent unit
(the slot is FED) or a fresh top vector (the slot is a SOURCE).  This file builds that input
(`gather`) and reads the answers back into the class table (`scatter`).

The routing is abstract here: tables `fl` (fed flag), `pp` (parent unit of a fed slot), `si`
(supply index of a source slot) and `slot` (child slot of a parent unit).  Their values for the
lockstep tree are computed elsewhere.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF
noncomputable section

universe U

/-! ## 1. Arithmetic of `placeS` addresses -/

lemma split3 {F Rr Qq A Bb : ℕ} (hF : F < 2^A) (hR : Rr < 2^Bb) :
    (F + 2^A * Rr + 2^(A+Bb) * Qq) % 2^A = F ∧
    ((F + 2^A * Rr + 2^(A+Bb) * Qq) / 2^A) % 2^Bb = Rr ∧
    (F + 2^A * Rr + 2^(A+Bb) * Qq) / 2^(A+Bb) = Qq := by
  have hA : 0 < 2^A := Nat.two_pow_pos A
  have hB : 0 < 2^Bb := Nat.two_pow_pos Bb
  have e1 : F + 2^A * Rr + 2^(A+Bb) * Qq = F + 2^A * (Rr + 2^Bb * Qq) := by
    rw [pow_add]; ring
  have e2 : (F + 2^A * (Rr + 2^Bb * Qq)) / 2^A = Rr + 2^Bb * Qq := by
    rw [Nat.add_mul_div_left _ _ hA, Nat.div_eq_of_lt hF, zero_add]
  refine ⟨?_, ?_, ?_⟩
  · rw [e1, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hF]
  · rw [e1, e2, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hR]
  · rw [e1, pow_add, ← Nat.div_div_eq_div_mul, e2, Nat.add_mul_div_left _ _ hB,
      Nat.div_eq_of_lt hR, zero_add]

section Place
variable {α β ρ : Type*} [Fintype α] [DecidableEq α] [Fintype β] [DecidableEq β]
  [Fintype ρ] [DecidableEq ρ]
variable {k : ℕ}

/-- The three parts of the address of a job entry. -/
lemma placeS_parts (S : Split α β (Fin k)) (X' : Layout β) (f r n : ℕ)
    (h : k*f + X'.n*f + r = n) (e : Bits n) :
    let qy := (placeS S X' f r n h).symm e
    (bits n).loc e % 2^(k*f) = (bits (k*f)).loc (S.peel 0 f (0, qy.2)).2 ∧
    ((bits n).loc e / 2^(k*f)) % 2^(X'.n*f) = (book X' 0 f).loc (S.peel 0 f (0, qy.2)).1 ∧
    (bits n).loc e / 2^(k*f + X'.n*f) = (bits r).loc qy.1 := by
  intro qy
  have he : e = placeS S X' f r n h (qy.1, qy.2) := by simp [qy]
  have hl := placeS_loc S X' f r n h qy.1 qy.2
  rw [← he] at hl
  rw [hl]
  have hF := (bits (k*f)).bd (S.peel 0 f (0, qy.2)).2
  have hR := (book X' 0 f).bd (S.peel 0 f (0, qy.2)).1
  rw [bits_len] at hF
  rw [book_len, Nat.add_zero] at hR
  exact split3 hF hR

/-- The job of one unit for a block of rank `k` on role `l`: the role read through `placeS`. -/
def blockJob (S : Split α β (Fin k)) (X' : Layout β) (f r n : ℕ)
    (h : k*f + X'.n*f + r = n) (l : ρ) (v : Bits r → Grain α ρ 0 f) : Bits n → ℂ :=
  fun e => v ((placeS S X' f r n h).symm e).1 (l, (0, ((placeS S X' f r n h).symm e).2))

lemma blockJob_placeS (S : Split α β (Fin k)) (X' : Layout β) (f r n : ℕ)
    (h : k*f + X'.n*f + r = n) (l : ρ) (v : Bits r → Grain α ρ 0 f) (q : Bits r) (y : Sky α f) :
    blockJob S X' f r n h l v (placeS S X' f r n h (q, y)) = v q (l, (0, y)) := by
  simp [blockJob]

/-- **The child's answer for a fed slot is the block letter on the parent unit.** -/
lemma blockJob_hlow (S : Split α β (Fin k)) (X' : Layout β) (f r n : ℕ)
    (h : k*f + X'.n*f + r = n) (l : ρ) (v : Bits r → Grain α ρ 0 f) (q : Bits r) (b0 : Bits 0)
    (y : Sky α f) :
    hlow n (k*f) (blockJob S X' f r n h l v) (placeS S X' f r n h (q, y)) =
      ferryM 0 f l (S.lift (kernel (Fin k))) (v q) (l, (b0, y)) := by
  rw [hlow_placeS]
  have hb : b0 = 0 := Subsingleton.elim _ _
  subst hb
  unfold ferryM roleAct
  simp only [if_true]
  congr 1
  funext y'
  rw [blockJob_placeS]

lemma fitk {X : Layout α} (X' : Layout β) (hX' : X'.n + k = X.n) {f r n : ℕ}
    (hf : X.n * f + r = n) : k*f + X'.n*f + r = n := by
  rw [← hf, ← hX', Nat.add_mul]; ring

end Place

/-! ## 2. The child's input -/

section Gather
variable {α β ρ : Type*} [Fintype α] [DecidableEq α] [Fintype β] [DecidableEq β]
  [Fintype ρ] [DecidableEq ρ]
variable {k : ℕ}

/-- Input job of child slot `t`: a parent unit's block job if fed, else a supply vector. -/
def cjob (X : Layout α) (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n) (l : ρ)
    {NU n : ℕ} {col ib : ℕ → ℕ} (hfit : ∀ p : Fin NU, X.n * col p + ib p = n)
    (V : ∀ p : Fin NU, Bits (ib p) → Grain α ρ 0 (col p))
    (fl pp si : ℕ → ℕ) (hpp : ∀ t, fl t ≠ 0 → pp t < NU) (sup : ℕ → Bits n → ℂ) (t : ℕ) :
    Bits n → ℂ :=
  if hf : fl t ≠ 0 then
    blockJob S X' (col (pp t)) (ib (pp t)) n (fitk X' hX' (hfit ⟨pp t, hpp t hf⟩)) l (V ⟨pp t, hpp t hf⟩)
  else sup (si t)

end Gather

section GatherRAM
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s}
  {n NU : ι → ℕ} {col ib : ι → ℕ → ℕ}
variable {α β ρ : Type*} [Fintype α] [DecidableEq α] [Fintype β] [DecidableEq β]
  [Fintype ρ] [DecidableEq ρ]
variable {cl : Paint} {X : Layout α} {R : Layout ρ} {k : ℕ}

/-- Entries of the child's input: slot `t`, entry `e`. -/
abbrev EntF (NJ : ι → ℕ) (n : ι → ℕ) := Σ u : UnitF NJ, Bits (n u.1)

/-- **The child's input** (fed slots from the parent's units, source slots from the supply):
linear work. -/
theorem CCtx.gather (h : CCtx m B s x X.n n NU col ib)
    {V : ∀ i (p : Fin (NU i)), Bits (ib i p) → Grain α ρ 0 (col i p)}
    (ct : Has m B s (a (a (a (c cl)))) x (fun i => ctab X R (NU i) (col i) (ib i) (V i)))
    (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n) (l : ρ)
    {NJ BS : ι → ℕ} {fl pp si : ι → ℕ → ℕ} {sup : ∀ i, ℕ → Bits (n i) → ℂ}
    (hpp : ∀ i t, fl i t ≠ 0 → pp i t < NU i)
    (hsi : ∀ i t, t < NJ i → fl i t = 0 → si i t < BS i)
    (hNJ : Knows m B s x NJ) (hNJs : Small B NJ)
    (tfl : Has m B s (a w) x (fun i => (Layout.std (NJ i)).tape (fun t => fl i t)))
    (tpp : Has m B s (a w) x (fun i => (Layout.std (NJ i)).tape (fun t => pp i t)))
    (tsi : Has m B s (a w) x (fun i => (Layout.std (NJ i)).tape (fun t => si i t)))
    (tsup : Has m B s (a (a (c cl))) x
      (fun i => (Layout.std (BS i)).tape (fun t => (bits (n i)).tape (sup i t)))) :
    Can m B s (a (a (c cl))) x
      (fun i => (Layout.std (NJ i)).tape (fun t => (bits (n i)).tape
        (cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) t)))
      (fun i => NJ i * 2^(n i)) := by
  let J (i : ι) (t : ℕ) : Bits (n i) → ℂ :=
    cjob X S X' hX' l (h.fit i) (V i) (fl i) (pp i) (si i) (hpp i) (sup i) t
  let yE (en : EntF NJ n) : T (p (p s w) w) := ((x en.1.1, en.1.2.val), (bits (n en.1.1)).loc en.2)
  let BE (en : EntF NJ n) := B en.1.1
  -- the environment of an entry contains the slot environment and the class environment
  have hslot : ∀ {t : Ty} {y : UnitF NJ → T t},
      Has m (fun u : UnitF NJ => B u.1) (p s w) t (fun u => (x u.1, u.2.val)) y →
      Has m BE (p (p s w) w) t yE (fun en => y en.1) :=
    fun g => Can.first.then_do (g.reindex (fun en : EntF NJ n => en.1))
  have hcls : ∀ {t : Ty} {y : ι → T t}, Has m B s t x y → Has m BE (p (p s w) w) t yE (fun en => y en.1.1) :=
    fun g => (Can.first.then_do Can.first).then_do (g.reindex (fun en : EntF NJ n => en.1.1))
  have hE : Knows m BE (p (p s w) w) yE (fun en => (bits (n en.1.1)).loc en.2) := Can.second
  have hbig : Small BE (fun en => 2^(n en.1.1)) := h.big.reindex (fun en : EntF NJ n => en.1.1)
  have hEs : Small BE (fun en => (bits (n en.1.1)).loc en.2) :=
    hbig.mono (fun en => ((bits (n en.1.1)).bd en.2).le)
  have hfl : Knows m BE (p (p s w) w) yE (fun en => fl en.1.1 en.1.2) := hslot (unit_fetch _ tfl)
  -- one entry
  have hent : Has m BE (p (p s w) w) (c cl) yE (fun en => J en.1.1 en.1.2 en.2) := by
    apply Can.testz hfl
    · -- source slot
      let K0 := {en : EntF NJ n // fl en.1.1 en.1.2 = 0}
      have hsiK : Knows m (fun κ : K0 => BE κ.1) (p (p s w) w) (fun κ => yE κ.1)
          (fun κ => si κ.1.1.1 κ.1.1.2) :=
        (hslot (unit_fetch _ tsi)).reindex (fun κ : K0 => κ.1)
      let ix (κ : K0) : Fin (BS κ.1.1.1) := ⟨si κ.1.1.1 κ.1.1.2, hsi _ _ κ.1.1.2.2 κ.2⟩
      have hv := Can.fetch (t:=a (c cl)) (L:=fun κ : K0 => Layout.std (BS κ.1.1.1))
        (v:=fun κ j => (bits (n κ.1.1.1)).tape (sup κ.1.1.1 j)) (i:=ix)
        ((hcls tsup).reindex (fun κ : K0 => κ.1)) hsiK
      have hw := Can.fetch (t:=c cl) (L:=fun κ : K0 => bits (n κ.1.1.1))
        (v:=fun κ => sup κ.1.1.1 (ix κ).val) (i:=fun κ => κ.1.2) hv (hE.reindex (fun κ : K0 => κ.1))
      refine hw.cong (fun _ => rfl) (fun κ => ?_)
      change sup κ.1.1.1 (si κ.1.1.1 κ.1.1.2) κ.1.2 = J κ.1.1.1 κ.1.1.2 κ.1.2
      simp only [J, cjob, κ.2, ne_eq, not_true_eq_false, dif_neg, not_false_eq_true]
    · -- fed slot: read the parent unit's role `l` through `placeS`
      let K1 := {en : EntF NJ n // fl en.1.1 en.1.2 ≠ 0}
      let BK (κ : K1) := BE κ.1
      let yK (κ : K1) := yE κ.1
      let i_ (κ : K1) := κ.1.1.1
      let p_ (κ : K1) : Fin (NU (i_ κ)) := ⟨pp (i_ κ) κ.1.1.2, hpp _ _ κ.2⟩
      let fit_ (κ : K1) := fitk X' hX' (h.fit (i_ κ) (p_ κ))
      let PL (κ : K1) := placeS S X' (col (i_ κ) (p_ κ)) (ib (i_ κ) (p_ κ)) (n (i_ κ)) (fit_ κ)
      let qy (κ : K1) := (PL κ).symm κ.1.2
      let j (κ : K1) : ElemF NU ib := ⟨⟨i_ κ, p_ κ⟩, (qy κ).1⟩
      have hppK : Knows m BK (p (p s w) w) yK (fun κ => (p_ κ).val) :=
        (hslot (unit_fetch _ tpp)).reindex (fun κ : K1 => κ.1)
      have hpcK : Knows m BK (p (p s w) w) yK (fun κ => 2^(col (i_ κ) (p_ κ))) :=
        Can.fetch (t:=w) (L:=fun κ : K1 => Layout.std (NU (i_ κ)))
          (v:=fun κ q => 2^(col (i_ κ) q)) (i:=p_)
          ((hcls h.tpc).reindex (fun κ : K1 => κ.1)) hppK
      have hcbig : Small BK (fun κ => 2^(col (i_ κ) (p_ κ))) :=
        (hbig.reindex (fun κ : K1 => κ.1)).mono (fun κ => h.two_c_le (i_ κ) (p_ κ))
      have hXc : Knows m BK (p (p s w) w) yK (fun κ => 2^(k * col (i_ κ) (p_ κ) +
          X'.n * col (i_ κ) (p_ κ))) :=
        (hpcK.wpow_fixed X.n hcbig).cong (fun _ => rfl) (fun κ => by
          rw [← pow_mul, ← hX']; ring_nf)
      have hq : Knows m BK (p (p s w) w) yK (fun κ => (bits (ib (i_ κ) (p_ κ))).loc (qy κ).1) :=
        ((hE.reindex (fun κ : K1 => κ.1)).dvdW hXc (hEs.reindex (fun κ : K1 => κ.1))).cong
          (fun _ => rfl) (fun κ => (placeS_parts S X' _ _ _ (fit_ κ) κ.1.2).2.2)
      have hx0 : Has m BE (p (p s w) w) s yE (fun en => x en.1.1) :=
        Can.first.then_do Can.first
      have hxK : Has m BK (p (p s w) w) s yK (fun κ => x (i_ κ)) :=
        hx0.reindex (fun κ : K1 => κ.1)
      have g : Has m (fun κ => B (j κ).1.1) (p (p s w) w) (p (p s w) w) yK
          (fun κ => ((x (j κ).1.1, (j κ).1.2.val), (bits (ib (j κ).1.1 (j κ).1.2)).loc (j κ).2)) :=
        (hxK.pair hppK).pair hq
      have G1 := (h.erig ct).gather j yK g
      have hk1 : Knows m BK (p (p s w) w) yK (fun κ => 2^(k * col (i_ κ) (p_ κ))) :=
        G1.rig.gear.powerB k
      have hk2 : Knows m BK (p (p s w) w) yK (fun κ => 2^(X'.n * col (i_ κ) (p_ κ))) :=
        G1.rig.gear.powerB X'.n
      have hfib : Knows m BK (p (p s w) w) yK
          (fun κ => (bits (k * col (i_ κ) (p_ κ))).loc (S.peel 0 (col (i_ κ) (p_ κ)) (0, (qy κ).2)).2) :=
        ((hE.reindex (fun κ : K1 => κ.1)).modW hk1 (hEs.reindex (fun κ : K1 => κ.1))).cong
          (fun _ => rfl) (fun κ => (placeS_parts S X' _ _ _ (fit_ κ) κ.1.2).1)
      have hrest : Knows m BK (p (p s w) w) yK
          (fun κ => (book X' 0 (col (i_ κ) (p_ κ))).loc (S.peel 0 (col (i_ κ) (p_ κ)) (0, (qy κ).2)).1) :=
        (((hE.reindex (fun κ : K1 => κ.1)).dvdW hk1 (hEs.reindex (fun κ : K1 => κ.1))).modW hk2
          ((hEs.reindex (fun κ : K1 => κ.1)).div)).cong
          (fun _ => rfl) (fun κ => (placeS_parts S X' _ _ _ (fit_ κ) κ.1.2).2.1)
      have hgrid := G1.rig.gear.to_gridB X 0 S X'
        (fun κ => S.peel 0 (col (i_ κ) (p_ κ)) (0, (qy κ).2)) hrest hfib
      have hb : Knows m BK (p (p s w) w) yK
          (fun κ => (book X 0 (col (i_ κ) (p_ κ))).loc ((0 : Bits 0), (qy κ).2)) :=
        hgrid.cong (fun _ => rfl) (fun κ => by rw [Equiv.symm_apply_apply])
      have hr := G1.rig.read (fun κ => (l, ((0 : Bits 0), (qy κ).2))) (Can.wc (R.loc l)) hb
      refine hr.cong (fun _ => rfl) (fun κ => ?_)
      change V (i_ κ) (p_ κ) (qy κ).1 (l, (0, (qy κ).2)) = J κ.1.1.1 κ.1.1.2 κ.1.2
      simp only [J, cjob, dif_pos κ.2]
      rfl
  -- inner layout over the entries of one slot, then outer layout over the slots
  have hlen : Knows m (fun u : UnitF NJ => B u.1) (p s w) (fun u => (x u.1, u.2.val))
      (fun u => (bits (n u.1)).n) :=
    Can.first.then_do (h.p2n.reindex (fun u : UnitF NJ => u.1))
  have HI := Can.layout (t:=c cl) (L:=fun u : UnitF NJ => bits (n u.1))
    (v:=fun u e => J u.1 u.2 e) (P:=fun _ => 0) (W:=fun _ => 0) hlen
    (hent.mono (fun _ => Nat.zero_le _)) (h.big.reindex (fun u : UnitF NJ => u.1))
  have HO := Can.layout (t:=a (c cl)) (L:=fun i => Layout.std (NJ i))
    (v:=fun i t => (bits (n i)).tape (J i t)) (P:=fun i => 2^(n i)) (W:=fun _ => 0) hNJ
    (HI.mono (fun u => by simp)) hNJs
  refine HO.weaken ⟨3, fun i => ?_⟩
  change 0 + NJ i + NJ i * 2^(n i) + 1 ≤ 3 * (NJ i * 2^(n i) + 1)
  have h1 : 1 ≤ 2^(n i) := Nat.one_le_two_pow
  have h2 : NJ i ≤ NJ i * 2^(n i) := Nat.le_mul_of_pos_right _ h1
  omega

end GatherRAM

section ScatterRAM
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s}
  {n NU : ι → ℕ} {col ib : ι → ℕ → ℕ}
variable {α β ρ : Type*} [Fintype α] [DecidableEq α] [Fintype β] [DecidableEq β]
  [Fintype ρ] [DecidableEq ρ]
variable {cl : Paint} {X : Layout α} {R : Layout ρ} {k : ℕ}

/-- **Reading the child's answers back** into the class table: role `l` of unit `p` comes from
child slot `slot p` (read through `placeS`), the other roles are kept.  Linear work. -/
theorem CCtx.scatter (h : CCtx m B s x X.n n NU col ib)
    {V : ∀ i (p : Fin (NU i)), Bits (ib i p) → Grain α ρ 0 (col i p)}
    (ct : Has m B s (a (a (a (c cl)))) x (fun i => ctab X R (NU i) (col i) (ib i) (V i)))
    (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n) (l : ρ)
    {NJ : ι → ℕ} {slot : ι → ℕ → ℕ} {Y : ∀ i, ℕ → Bits (n i) → ℂ}
    (hslot : ∀ i (p : Fin (NU i)), slot i p < NJ i)
    (tslot : Has m B s (a w) x (fun i => (Layout.std (NU i)).tape (fun p => slot i p)))
    (tout : Has m B s (a (a (c cl))) x
      (fun i => (Layout.std (NJ i)).tape (fun t => (bits (n i)).tape (Y i t))))
    (hY : ∀ i (p : Fin (NU i)) (q : Bits (ib i p)) (y : Sky α (col i p)),
      Y i (slot i p) (placeS S X' (col i p) (ib i p) (n i) (fitk X' hX' (h.fit i p)) (q, y)) =
        ferryM 0 (col i p) l (S.lift (kernel (Fin k))) (V i p q) (l, (0, y))) :
    Can m B s (a (a (a (c cl)))) x
      (fun i => ctab X R (NU i) (col i) (ib i)
        (fun p q => ferryM 0 (col i p) l (S.lift (kernel (Fin k))) (V i p q)))
      (fun i => NU i * 2^(n i)) := by
  have G := h.erig ct
  apply h.ctab_build _ rfl
  let yE (e : ElemF NU ib) : T (p (p s w) w) :=
    ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2)
  let BE (e : ElemF NU ib) := B e.1.1
  have hel : Can m BE (p (p s w) w) (a (c cl)) yE
      (fun e => (fleet X R 0 (col e.1.1 e.1.2)).tape
        (ferryM 0 (col e.1.1 e.1.2) l (S.lift (kernel (Fin k))) (V e.1.1 e.1.2 e.2)))
      (fun e => (fleet X R 0 (col e.1.1 e.1.2)).n) := by
    refine G.rig.chart _ ?_
    -- one beach of one element
    let D := Σ e : ElemF NU ib, Beach α ρ 0 (col e.1.1 e.1.2)
    let yD (d : D) : T (p (p (p s w) w) w) := (yE d.1, (fleet X R 0 (col d.1.1.1 d.1.1.2)).loc d.2)
    let BD (d : D) := BE d.1
    have GD := G.rig.gather (fun d : D => d.1) yD (t:=p (p (p s w) w) w) Can.first
    obtain ⟨hl, hk⟩ := GD.addrs (fun d => d.2) Can.second
    have hbD : Has m BD (p (p (p s w) w) w) (c cl) yD (fun d => V d.1.1.1 d.1.1.2 d.1.2 d.2) :=
      GD.read (fun d => d.2) hl hk
    apply Can.cases_code (q:=fun d : D => d.2.1) R hl
    intro r'
    by_cases hr : r' = l
    · -- role `l`: from the child's answer
      let K := {d : D // d.2.1 = r'}
      let BK (κ : K) := BD κ.1
      have GK := GD.gather (fun κ : K => κ.1) (fun κ => yD κ.1) (t:=p (p (p s w) w) w) Can.id
      obtain ⟨hrest, hfib⟩ := GK.gear.from_gridB X 0 S X' (fun κ => κ.1.2.2)
        (hk.reindex (fun κ : K => κ.1))
      have hsU : Knows m (fun u : UnitF NU => B u.1) (p s w) (fun u => (x u.1, u.2.val))
          (fun u => slot u.1 u.2) := unit_fetch _ tslot
      have hsE : Knows m BE (p (p s w) w) yE (fun e => slot e.1.1 e.1.2) :=
        Can.first.then_do (hsU.reindex (fun e : ElemF NU ib => e.1))
      have hsD : Knows m BD (p (p (p s w) w) w) yD (fun d => slot d.1.1.1 d.1.1.2) :=
        Can.first.then_do (hsE.reindex (fun d : D => d.1))
      have hslotK : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => slot κ.1.1.1.1 κ.1.1.1.2) := hsD.reindex (fun κ : K => κ.1)
      have hqD : Knows m BD (p (p (p s w) w) w) yD (fun d => (bits (ib d.1.1.1 d.1.1.2)).loc d.1.2) :=
        Can.first.then_do Can.second
      have hq : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => (bits (ib κ.1.1.1.1 κ.1.1.1.2)).loc κ.1.1.2) := hqD.reindex (fun κ : K => κ.1)
      have hk1 : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => 2^(k * col κ.1.1.1.1 κ.1.1.1.2)) := GK.gear.powerB k
      have hk3 : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => 2^(k * col κ.1.1.1.1 κ.1.1.1.2 + X'.n * col κ.1.1.1.1 κ.1.1.1.2)) :=
        (GK.gear.powerB X.n).cong (fun _ => rfl) (fun κ => by rw [← hX']; ring_nf)
      have hbig : Small BK (fun κ => 2^(n κ.1.1.1.1)) :=
        h.big.reindex (fun κ : K => κ.1.1.1.1)
      have hcb : Small BK (fun κ => 2^(col κ.1.1.1.1 κ.1.1.1.2)) :=
        hbig.mono (fun κ => h.two_c_le _ _)
      have hs1 : Small BK (fun κ => 2^(k * col κ.1.1.1.1 κ.1.1.1.2)) :=
        GK.gear.boundB k
      have hs3 : Small BK (fun κ => 2^(k * col κ.1.1.1.1 κ.1.1.1.2 + X'.n * col κ.1.1.1.1 κ.1.1.1.2)) :=
        (GK.gear.boundB X.n).mono (fun κ => le_of_eq (by rw [← hX']; ring_nf))
      have hsr : Small BK (fun κ => (book X' 0 (col κ.1.1.1.1 κ.1.1.1.2)).loc
          (S.peel 0 (col κ.1.1.1.1 κ.1.1.1.2) κ.1.2.2).1) :=
        (GK.gear.grid_small X' 0).mono (fun κ => ((book X' 0 _).bd _).le)
      have hsq : Small BK (fun κ => (bits (ib κ.1.1.1.1 κ.1.1.1.2)).loc κ.1.1.2) :=
        hbig.mono (fun κ => le_trans ((bits _).bd _).le (h.two_ib_le _ _))
      have hsf : Small BK (fun κ => (bits (k * col κ.1.1.1.1 κ.1.1.1.2)).loc
          (S.peel 0 (col κ.1.1.1.1 κ.1.1.1.2) κ.1.2.2).2) :=
        hs1.mono (fun κ => ((bits _).bd _).le)
      -- the address inside the child's job
      have haddr : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => (bits (k * col κ.1.1.1.1 κ.1.1.1.2)).loc
              (S.peel 0 (col κ.1.1.1.1 κ.1.1.1.2) κ.1.2.2).2 +
            2^(k * col κ.1.1.1.1 κ.1.1.1.2) * (book X' 0 (col κ.1.1.1.1 κ.1.1.1.2)).loc
              (S.peel 0 (col κ.1.1.1.1 κ.1.1.1.2) κ.1.2.2).1 +
            2^(k * col κ.1.1.1.1 κ.1.1.1.2 + X'.n * col κ.1.1.1.1 κ.1.1.1.2) *
              (bits (ib κ.1.1.1.1 κ.1.1.1.2)).loc κ.1.1.2) :=
        (hfib.plus (hk1.times hrest hs1 hsr) hsf (hs1.mul hsr)).plus (hk3.times hq hs3 hsq)
          (hsf.add (hs1.mul hsr)) (hs3.mul hsq)
      let jq (κ : K) := κ.1.1.2
      let jy (κ : K) := κ.1.2.2.2
      let PL (κ : K) := placeS S X' (col κ.1.1.1.1 κ.1.1.1.2) (ib κ.1.1.1.1 κ.1.1.1.2)
        (n κ.1.1.1.1) (fitk X' hX' (h.fit κ.1.1.1.1 κ.1.1.1.2))
      have haddr' : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => (bits (n κ.1.1.1.1)).loc (PL κ (jq κ, jy κ))) := by
        refine haddr.cong (fun _ => rfl) (fun κ => ?_)
        rw [placeS_loc]
        have hb0 : κ.1.2.2 = ((0 : Bits 0), jy κ) := Prod.ext (Subsingleton.elim _ _) rfl
        rw [hb0]
      have hoU : Has m (fun u : UnitF NU => B u.1) (p s w) (a (a (c cl))) (fun u => (x u.1, u.2.val))
          (fun u => (Layout.std (NJ u.1)).tape (fun t => (bits (n u.1)).tape (Y u.1 t))) :=
        Can.first.then_do (tout.reindex (fun u : UnitF NU => u.1))
      have hoE : Has m BE (p (p s w) w) (a (a (c cl))) yE
          (fun e => (Layout.std (NJ e.1.1)).tape (fun t => (bits (n e.1.1)).tape (Y e.1.1 t))) :=
        Can.first.then_do (hoU.reindex (fun e : ElemF NU ib => e.1))
      have hoD : Has m BD (p (p (p s w) w) w) (a (a (c cl))) yD
          (fun d => (Layout.std (NJ d.1.1.1)).tape (fun t => (bits (n d.1.1.1)).tape (Y d.1.1.1 t))) :=
        Can.first.then_do (hoE.reindex (fun d : D => d.1))
      have hI := hoD.reindex (fun κ : K => κ.1)
      have hJ := Can.fetch (t:=a (c cl)) (L:=fun κ : K => Layout.std (NJ κ.1.1.1.1))
        (v:=fun κ t => (bits (n κ.1.1.1.1)).tape (Y κ.1.1.1.1 t))
        (i:=fun κ => (⟨slot κ.1.1.1.1 κ.1.1.1.2, hslot _ _⟩ : Fin (NJ κ.1.1.1.1))) hI hslotK
      have hv := Can.fetch (t:=c cl) (L:=fun κ : K => bits (n κ.1.1.1.1))
        (v:=fun κ => Y κ.1.1.1.1 (slot κ.1.1.1.1 κ.1.1.1.2)) (i:=fun κ => PL κ (jq κ, jy κ)) hJ haddr'
      refine hv.cong (fun _ => rfl) (fun κ => ?_)
      change Y κ.1.1.1.1 (slot κ.1.1.1.1 κ.1.1.1.2) (PL κ (jq κ, jy κ)) =
        ferryM 0 _ l (S.lift (kernel (Fin k))) (V κ.1.1.1.1 κ.1.1.1.2 κ.1.1.2) κ.1.2
      rw [hY]
      have hd : κ.1.2 = (l, ((0 : Bits 0), jy κ)) :=
        Prod.ext (κ.2.trans hr) (Prod.ext (Subsingleton.elim _ _) rfl)
      rw [hd]
    · -- other roles: unchanged
      refine (hbD.reindex (fun κ : {d : D // d.2.1 = r'} => κ.1)).cong (fun _ => rfl) (fun κ => ?_)
      have H : κ.1.2.1 ≠ l := fun hh => hr (κ.2.symm.trans hh)
      simp only [ferryM, roleAct, H, ite_false, Prod.mk.eta]
  refine hel.mono (fun e => le_of_eq ?_)
  simp [fleet]

end ScatterRAM

end
end PowerSaving.RAM
end OAI
