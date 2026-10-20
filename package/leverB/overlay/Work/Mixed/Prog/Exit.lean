import Work.Mixed.Prog.Block

/-!
# Mixed engine, R3 program part 5: the exits of a leaf class (agent e-r3prog)

In a leaf class every unit has the same column `c0`.  A block letter of rank `k` makes ONE job per
unit (`2^n` entries, `k*c0` H-bits low); the `NU` jobs are cut into `NU * 2^(n - k*c0)` vectors of
`2^(k*c0)` entries and handed to ONE fresh top call `Do(k*c0, NU * 2^(n - k*c0))`
(rl-prove Lemma D).

* `hlow_catB`        H on the low `k` bits of `catB (x, y)` is `ripple k` in `x` for each `y`;
* `toVec`, `fromVec`  reshaping jobs into vectors and back (`can_toVec`, `can_fromVec`);
* `CCtx.exit_letter` gather (all slots fed), reshape, ONE top call, reshape, read back.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF
noncomputable section

universe U

/-- **H on the low bits of a split vector.** -/
lemma hlow_catB (k b : ℕ) (v : Bits (k+b) → ℂ) (x : Bits k) (y : Bits b) :
    hlow (k+b) k v (catB k b (x, y)) = ripple k (fun x' => v (catB k b (x', y))) x := by
  unfold hlow ripple
  simp only [Matrix.mulVec, dotProduct]
  rw [← Equiv.sum_comp (catB k b), Fintype.sum_prod_type]
  simp_rw [maskK_catB]
  apply Finset.sum_congr rfl
  intro x' _
  rw [Finset.sum_eq_single y]
  · simp
  · intro y' _ hy; rw [if_neg (Ne.symm hy)]; simp
  · intro hh; exact absurd (Finset.mem_univ _) hh

lemma hlow_castB_catB (n k : ℕ) (hk : k ≤ n) (v : Bits n → ℂ) (x : Bits k) (y : Bits (n-k)) :
    hlow n k v (castB (Nat.add_sub_cancel' hk) (catB k (n-k) (x, y))) =
      ripple k (fun x' => v (castB (Nat.add_sub_cancel' hk) (catB k (n-k) (x', y)))) x := by
  have H := hlow_catB k (n-k) (fun z => v (castB (Nat.add_sub_cancel' hk) z)) x y
  have e : ∀ (m : ℕ) (hm : k + (n-k) = m) (u : Bits m → ℂ) (z : Bits (k + (n-k))),
      hlow m k u (castB hm z) = hlow (k + (n-k)) k (fun z' => u (castB hm z')) z := by
    intro m hm u z; subst hm; rfl
  rw [e n _ v]
  exact H

/-- An element of `Bits b` from an address (taken modulo `2^b`). -/
def bitsOf (b j : ℕ) : Bits b := (bits b).atFin ⟨j % 2^b, Nat.mod_lt _ (Nat.two_pow_pos b)⟩

lemma bitsOf_loc (b j : ℕ) : (bits b).loc (bitsOf b j) = j % 2^b := by
  unfold bitsOf; rw [Layout.loc_atFin]

lemma bitsOf_of_loc (b : ℕ) (y : Bits b) : bitsOf b ((bits b).loc y) = y := by
  apply (bits b).inj
  rw [bitsOf_loc, Nat.mod_eq_of_lt (by have := (bits b).bd y; rwa [bits_len] at this)]

/-- Jobs (`2^n` entries, H on the low `k` bits) cut into vectors of `2^k` entries:
vector `t = p * 2^(n-k) + loc y` of job `p`. -/
def toVec (n k : ℕ) (hk : k ≤ n) (J : ℕ → Bits n → ℂ) (t : ℕ) : Bits k → ℂ :=
  fun x => J (t / 2^(n-k)) (castB (Nat.add_sub_cancel' hk) (catB k (n-k) (x, bitsOf (n-k) t)))

/-- Vectors put back into jobs. -/
def fromVec (n k : ℕ) (Y : ℕ → Bits k → ℂ) (p : ℕ) : Bits n → ℂ :=
  fun e => Y (p * 2^(n-k) + (bits n).loc e / 2^k) (bitsOf k ((bits n).loc e))

lemma loc_castB_catB (n k : ℕ) (hk : k ≤ n) (x : Bits k) (y : Bits (n-k)) :
    (bits n).loc (castB (Nat.add_sub_cancel' hk) (catB k (n-k) (x, y))) =
      (bits k).loc x + 2^k * (bits (n-k)).loc y := by
  rw [castB_loc, catB_loc]

/-- The exits answered and put back are `hlow` of the jobs. -/
lemma fromVec_ripple (n k : ℕ) (hk : k ≤ n) (J : ℕ → Bits n → ℂ) (p : ℕ) (e : Bits n) :
    fromVec n k (fun t => ripple k (toVec n k hk J t)) p e = hlow n k (J p) e := by
  -- write `e` as a split vector
  obtain ⟨⟨x, y⟩, hxy⟩ : ∃ z : Bits k × Bits (n-k),
      castB (Nat.add_sub_cancel' hk) (catB k (n-k) z) = e :=
    ⟨(catB k (n-k)).symm ((castB (Nat.add_sub_cancel' hk)).symm e), by simp⟩
  subst hxy
  rw [hlow_castB_catB n k hk (J p) x y]
  unfold fromVec toVec
  have hx : (bits k).loc x < 2^k := by have := (bits k).bd x; rwa [bits_len] at this
  have hy : (bits (n-k)).loc y < 2^(n-k) := by have := (bits (n-k)).bd y; rwa [bits_len] at this
  have h2 : 0 < 2^k := Nat.two_pow_pos k
  have h3 : 0 < 2^(n-k) := Nat.two_pow_pos (n-k)
  rw [loc_castB_catB n k hk]
  have e1 : ((bits k).loc x + 2^k * (bits (n-k)).loc y) / 2^k = (bits (n-k)).loc y := by
    rw [Nat.add_mul_div_left _ _ h2, Nat.div_eq_of_lt hx, zero_add]
  have e2 : bitsOf k ((bits k).loc x + 2^k * (bits (n-k)).loc y) = x := by
    apply (bits k).inj
    rw [bitsOf_loc, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt hx]
  rw [e1, e2]
  have e3 : (p * 2^(n-k) + (bits (n-k)).loc y) / 2^(n-k) = p := by
    rw [Nat.add_comm, Nat.add_mul_div_right _ _ h3, Nat.div_eq_of_lt hy, zero_add]
  have e4 : bitsOf (n-k) (p * 2^(n-k) + (bits (n-k)).loc y) = y := by
    apply (bits (n-k)).inj
    rw [bitsOf_loc, Nat.add_comm, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt hy]
  simp only [e3, e4]

section Reshape
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s} {cl : Paint}
variable {n k NJ : ι → ℕ}

/-- Cutting jobs into vectors: linear work. -/
lemma can_toVec (hk : ∀ i, k i ≤ n i) {J : ∀ i, ℕ → Bits (n i) → ℂ}
    (hJ : Has m B s (jobsT cl) x (fun i => jobsTape (n i) (NJ i) (J i)))
    (hNJ : Knows m B s x NJ) (hp2 : Knows m B s x (fun i => 2^(n i - k i)))
    (hp1 : Knows m B s x (fun i => 2^(k i)))
    (hbig : Small B (fun i => 2^(n i))) (hNJs : Small B NJ) :
    Can m B s (jobsT cl) x
      (fun i => jobsTape (k i) (NJ i * 2^(n i - k i)) (toVec (n i) (k i) (hk i) (J i)))
      (fun i => NJ i * 2^(n i)) := by
  have hle2 : ∀ i, 2^(n i - k i) ≤ 2^(n i) := fun i => Nat.pow_le_pow_right (by omega) (Nat.sub_le _ _)
  have hle1 : ∀ i, 2^(k i) ≤ 2^(n i) := fun i => Nat.pow_le_pow_right (by omega) (hk i)
  have hs2 : Small B (fun i => 2^(n i - k i)) := hbig.mono hle2
  have hs1 : Small B (fun i => 2^(k i)) := hbig.mono hle1
  have hNV : Knows m B s x (fun i => NJ i * 2^(n i - k i)) := hNJ.times hp2 hNJs hs2
  have hNVs : Small B (fun i => NJ i * 2^(n i - k i)) := hNJs.mul hs2
  -- one entry: vector t, entry x
  let Dv := Σ u : UnitF (fun i => NJ i * 2^(n i - k i)), Bits (k u.1)
  let yv (d : Dv) : T (p (p s w) w) := ((x d.1.1, d.1.2.val), (bits (k d.1.1)).loc d.2)
  let Bv (d : Dv) := B d.1.1
  have hcls : ∀ {t : Ty} {y : ι → T t}, Has m B s t x y → Has m Bv (p (p s w) w) t yv (fun d => y d.1.1) :=
    fun g => (Can.first.then_do Can.first).then_do (g.reindex (fun d : Dv => d.1.1))
  have ht : Knows m Bv (p (p s w) w) yv (fun d => d.1.2.val) := Can.first.then_do Can.second
  have hxl : Knows m Bv (p (p s w) w) yv (fun d => (bits (k d.1.1)).loc d.2) := Can.second
  have hts : Small Bv (fun d => d.1.2.val) := (hNVs.reindex (fun d : Dv => d.1.1)).mono (fun d => d.1.2.2.le)
  have hjob : Knows m Bv (p (p s w) w) yv (fun d => d.1.2.val / 2^(n d.1.1 - k d.1.1)) :=
    ht.dvdW (hcls hp2) hts
  have hjlt : ∀ d : Dv, d.1.2.val / 2^(n d.1.1 - k d.1.1) < NJ d.1.1 := fun d =>
    (Nat.div_lt_iff_lt_mul (Nat.two_pow_pos _)).2 d.1.2.2
  have hJt := Can.fetch (t:=a (c cl)) (L:=fun d : Dv => Layout.std (NJ d.1.1))
    (v:=fun d j => (bits (n d.1.1)).tape (J d.1.1 j))
    (i:=fun d => (⟨d.1.2.val / 2^(n d.1.1 - k d.1.1), hjlt d⟩ : Fin (NJ d.1.1))) (hcls hJ) hjob
  have hmod : Knows m Bv (p (p s w) w) yv (fun d => d.1.2.val % 2^(n d.1.1 - k d.1.1)) :=
    ht.modW (hcls hp2) hts
  have hxs : Small Bv (fun d => (bits (k d.1.1)).loc d.2) :=
    (hs1.reindex (fun d : Dv => d.1.1)).mono (fun d => ((bits _).bd _).le)
  have hms : Small Bv (fun d => d.1.2.val % 2^(n d.1.1 - k d.1.1)) := hts.mod
  have haddr : Knows m Bv (p (p s w) w) yv (fun d => (bits (n d.1.1)).loc
      (castB (Nat.add_sub_cancel' (hk d.1.1)) (catB (k d.1.1) (n d.1.1 - k d.1.1)
        (d.2, bitsOf (n d.1.1 - k d.1.1) d.1.2.val)))) :=
    (hxl.plus ((hcls hp1).times hmod (hs1.reindex _) hms) hxs ((hs1.reindex _).mul hms)).cong
      (fun _ => rfl) (fun d => by rw [loc_castB_catB _ _ (hk d.1.1), bitsOf_loc])
  have hv := Can.fetch (t:=c cl) (L:=fun d : Dv => bits (n d.1.1))
    (v:=fun d => J d.1.1 (d.1.2.val / 2^(n d.1.1 - k d.1.1)))
    (i:=fun d => castB (Nat.add_sub_cancel' (hk d.1.1)) (catB (k d.1.1) (n d.1.1 - k d.1.1)
        (d.2, bitsOf (n d.1.1 - k d.1.1) d.1.2.val))) hJt haddr
  have hlen : Knows m (fun u : UnitF (fun i => NJ i * 2^(n i - k i)) => B u.1) (p s w)
      (fun u => (x u.1, u.2.val)) (fun u => (bits (k u.1)).n) :=
    Can.first.then_do (hp1.reindex (fun u : UnitF (fun i => NJ i * 2^(n i - k i)) => u.1))
  have HI := Can.layout (t:=c cl) (L:=fun u : UnitF (fun i => NJ i * 2^(n i - k i)) => bits (k u.1))
    (v:=fun u e => toVec (n u.1) (k u.1) (hk u.1) (J u.1) u.2 e) (P:=fun _ => 0) (W:=fun _ => 0) hlen
    (hv.mono (fun _ => Nat.zero_le _)) (hs1.reindex _)
  have HO := Can.layout (t:=a (c cl)) (L:=fun i => Layout.std (NJ i * 2^(n i - k i)))
    (v:=fun i t => (bits (k i)).tape (toVec (n i) (k i) (hk i) (J i) t)) (P:=fun i => 2^(k i))
    (W:=fun _ => 0) hNV (HI.mono (fun u => by simp)) hNVs
  refine HO.weaken ⟨3, fun i => ?_⟩
  change 0 + NJ i * 2^(n i - k i) + NJ i * 2^(n i - k i) * 2^(k i) + 1 ≤ 3 * (NJ i * 2^(n i) + 1)
  have e1 : NJ i * 2^(n i - k i) * 2^(k i) = NJ i * 2^(n i) := by
    rw [mul_assoc, ← pow_add, Nat.sub_add_cancel (hk i)]
  have e2 : NJ i * 2^(n i - k i) ≤ NJ i * 2^(n i) := Nat.mul_le_mul_left _ (hle2 i)
  omega

/-- Putting vectors back into jobs: linear work. -/
lemma can_fromVec (hk : ∀ i, k i ≤ n i) {Y : ∀ i, ℕ → Bits (k i) → ℂ}
    (hY : Has m B s (jobsT cl) x (fun i => jobsTape (k i) (NJ i * 2^(n i - k i)) (Y i)))
    (hNJ : Knows m B s x NJ) (hp2 : Knows m B s x (fun i => 2^(n i - k i)))
    (hp1 : Knows m B s x (fun i => 2^(k i))) (hpn : Knows m B s x (fun i => 2^(n i)))
    (hbig : Small B (fun i => 2^(n i))) (hNJs : Small B NJ) :
    Can m B s (jobsT cl) x
      (fun i => jobsTape (n i) (NJ i) (fromVec (n i) (k i) (Y i)))
      (fun i => NJ i * 2^(n i)) := by
  have hle2 : ∀ i, 2^(n i - k i) ≤ 2^(n i) := fun i => Nat.pow_le_pow_right (by omega) (Nat.sub_le _ _)
  have hle1 : ∀ i, 2^(k i) ≤ 2^(n i) := fun i => Nat.pow_le_pow_right (by omega) (hk i)
  have hs2 : Small B (fun i => 2^(n i - k i)) := hbig.mono hle2
  have hs1 : Small B (fun i => 2^(k i)) := hbig.mono hle1
  let De := Σ u : UnitF NJ, Bits (n u.1)
  let ye (d : De) : T (p (p s w) w) := ((x d.1.1, d.1.2.val), (bits (n d.1.1)).loc d.2)
  let Be (d : De) := B d.1.1
  have hcls : ∀ {t : Ty} {y : ι → T t}, Has m B s t x y → Has m Be (p (p s w) w) t ye (fun d => y d.1.1) :=
    fun g => (Can.first.then_do Can.first).then_do (g.reindex (fun d : De => d.1.1))
  have hp : Knows m Be (p (p s w) w) ye (fun d => d.1.2.val) := Can.first.then_do Can.second
  have hE : Knows m Be (p (p s w) w) ye (fun d => (bits (n d.1.1)).loc d.2) := Can.second
  have hEs : Small Be (fun d => (bits (n d.1.1)).loc d.2) :=
    (hbig.reindex (fun d : De => d.1.1)).mono (fun d => ((bits _).bd _).le)
  have hps : Small Be (fun d => d.1.2.val) := (hNJs.reindex (fun d : De => d.1.1)).mono (fun d => d.1.2.2.le)
  have hvi : Knows m Be (p (p s w) w) ye
      (fun d => d.1.2.val * 2^(n d.1.1 - k d.1.1) + (bits (n d.1.1)).loc d.2 / 2^(k d.1.1)) :=
    (hp.times (hcls hp2) hps (hs2.reindex _)).plus (hE.dvdW (hcls hp1) hEs)
      (hps.mul (hs2.reindex _)) hEs.div
  have hvlt : ∀ d : De, d.1.2.val * 2^(n d.1.1 - k d.1.1) + (bits (n d.1.1)).loc d.2 / 2^(k d.1.1) <
      NJ d.1.1 * 2^(n d.1.1 - k d.1.1) := by
    intro d
    have h1 : (bits (n d.1.1)).loc d.2 < 2^(k d.1.1) * 2^(n d.1.1 - k d.1.1) := by
      have := (bits (n d.1.1)).bd d.2
      rw [bits_len, ← pow_add, Nat.add_sub_cancel' (hk d.1.1)] at *
      exact this
    have h2 : (bits (n d.1.1)).loc d.2 / 2^(k d.1.1) < 2^(n d.1.1 - k d.1.1) :=
      (Nat.div_lt_iff_lt_mul (Nat.two_pow_pos _)).2 (by rw [Nat.mul_comm]; exact h1)
    have h3 : d.1.2.val + 1 ≤ NJ d.1.1 := d.1.2.2
    calc d.1.2.val * 2^(n d.1.1 - k d.1.1) + (bits (n d.1.1)).loc d.2 / 2^(k d.1.1)
        < d.1.2.val * 2^(n d.1.1 - k d.1.1) + 2^(n d.1.1 - k d.1.1) := by omega
      _ = (d.1.2.val + 1) * 2^(n d.1.1 - k d.1.1) := by ring
      _ ≤ NJ d.1.1 * 2^(n d.1.1 - k d.1.1) := Nat.mul_le_mul_right _ h3
  have hVt := Can.fetch (t:=a (c cl)) (L:=fun d : De => Layout.std (NJ d.1.1 * 2^(n d.1.1 - k d.1.1)))
    (v:=fun d j => (bits (k d.1.1)).tape (Y d.1.1 j))
    (i:=fun d => (⟨_, hvlt d⟩ : Fin (NJ d.1.1 * 2^(n d.1.1 - k d.1.1)))) (hcls hY) hvi
  have hxa : Knows m Be (p (p s w) w) ye
      (fun d => (bits (k d.1.1)).loc (bitsOf (k d.1.1) ((bits (n d.1.1)).loc d.2))) :=
    (hE.modW (hcls hp1) hEs).cong (fun _ => rfl) (fun d => by rw [bitsOf_loc])
  have hv := Can.fetch (t:=c cl) (L:=fun d : De => bits (k d.1.1))
    (v:=fun d => Y d.1.1 (d.1.2.val * 2^(n d.1.1 - k d.1.1) + (bits (n d.1.1)).loc d.2 / 2^(k d.1.1)))
    (i:=fun d => bitsOf (k d.1.1) ((bits (n d.1.1)).loc d.2)) hVt hxa
  have hlen : Knows m (fun u : UnitF NJ => B u.1) (p s w) (fun u => (x u.1, u.2.val))
      (fun u => (bits (n u.1)).n) := Can.first.then_do (hpn.reindex (fun u : UnitF NJ => u.1))
  have HI := Can.layout (t:=c cl) (L:=fun u : UnitF NJ => bits (n u.1))
    (v:=fun u e => fromVec (n u.1) (k u.1) (Y u.1) u.2 e) (P:=fun _ => 0) (W:=fun _ => 0) hlen
    (hv.mono (fun _ => Nat.zero_le _)) (hbig.reindex _)
  have HO := Can.layout (t:=a (c cl)) (L:=fun i => Layout.std (NJ i))
    (v:=fun i t => (bits (n i)).tape (fromVec (n i) (k i) (Y i) t)) (P:=fun i => 2^(n i))
    (W:=fun _ => 0) hNJ (HI.mono (fun u => by simp)) hNJs
  refine HO.weaken ⟨3, fun i => ?_⟩
  change 0 + NJ i + NJ i * 2^(n i) + 1 ≤ 3 * (NJ i * 2^(n i) + 1)
  have h2 : NJ i ≤ NJ i * 2^(n i) := Nat.le_mul_of_pos_right _ (Nat.two_pow_pos _)
  omega

end Reshape

/-- Some closed program of any type (only used on empty families). -/
def anyProg (m : Bool) : (s t : Ty) → Prog m s t
  | _, .w => .atom (.lit 0)
  | _, .c e => .atom (.cz e)
  | s, .p u v => .fork (anyProg m s u) (anyProg m s v)
  | s, .a v => .tab (.atom (.lit 0)) (anyProg m (p s w) v)

lemma can_of_empty {ι : Type U} {m : Bool} {B : ι → ℕ} {s t : Ty} {x : ι → T s} {y : ι → T t}
    {W : ι → ℕ} (hE : ∀ i : ι, False) : Can m B s t x y W :=
  ⟨anyProg m s t, 0, 0, fun i => (hE i).elim⟩

/-- A table of values computable at constant cost from the index: linear work. -/
lemma can_tabIx {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s} {t : Ty}
    {N : ι → ℕ} (g : ∀ i, Fin (N i) → T t) (hN : Knows m B s x N) (hNs : Small B N)
    (hg : Has m (fun u : UnitF N => B u.1) (p s w) t (fun u => (x u.1, u.2.val)) (fun u => g u.1 u.2)) :
    Can m B s (a t) x (fun i => (Layout.std (N i)).tape (g i)) N := by
  have H := Can.layout (t:=t) (L:=fun i => Layout.std (N i)) (v:=g) (P:=fun _ => 0) (W:=fun _ => 0)
    hN (hg.mono (fun _ => Nat.zero_le _)) hNs
  refine H.mono (fun i => ?_)
  change 0 + N i + N i * 0 ≤ N i
  omega

/-- The empty tape of jobs. -/
lemma has_noJobs {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s} {cl : Paint}
    (n : ι → ℕ) : Has m B s (jobsT cl) x (fun i => jobsTape (n i) 0 (fun _ _ => 0)) := by
  have H := can_tabIx (t:=a (c cl)) (N:=fun _ => 0) (x:=x) (B:=B) (m:=m)
    (fun i (t : Fin 0) => (bits (n i)).tape (fun _ => (0:ℂ))) (Can.wc 0) (Small.const _ 0)
    (can_of_empty (fun u => Nat.not_lt_zero _ u.2.2))
  exact H

section ExitRAM
variable {ι : Type U} {m : Bool} {cl : Paint} {b : Bank (portM cl) ι} {s : Ty} {x : ι → T s}
  {n NU : ι → ℕ} {col ib : ι → ℕ → ℕ}
variable {α β ρ : Type*} [Fintype α] [DecidableEq α] [Fintype β] [DecidableEq β]
  [Fintype ρ] [DecidableEq ρ]
variable {X : Layout α} {R : Layout ρ} {k : ℕ}

/-- The fed flag of the exit routing: every unit is fed (to itself). -/
def exFl (N t : ℕ) : ℕ := if t < N then 1 else 0

lemma exFl_pp (N t : ℕ) (h : exFl N t ≠ 0) : t < N := by
  unfold exFl at h; split_ifs at h with h1
  · exact h1
  · exact absurd rfl h

/-- The jobs of a leaf class for a block letter (one per unit). -/
def exJobs (X : Layout α) (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n) (l : ρ)
    {NU n : ℕ} {col ib : ℕ → ℕ} (hfit : ∀ p : Fin NU, X.n * col p + ib p = n)
    (V : ∀ p : Fin NU, Bits (ib p) → Grain α ρ 0 (col p)) : ℕ → Bits n → ℂ :=
  cjob X S X' hX' l hfit V (exFl NU) id (fun _ => 0) (exFl_pp NU) (fun _ _ => 0)

/-- **One block letter of a leaf class: ONE top call on all its units' jobs (an exit).** -/
theorem CCtx.exit_letter (h : CCtx m b.B s x X.n n NU col ib)
    {V : ∀ i (p : Fin (NU i)), Bits (ib i p) → Grain α ρ 0 (col i p)}
    (ct : Has m b.B s (a (a (a (c cl)))) x (fun i => ctab X R (NU i) (col i) (ib i) (V i)))
    (S : Split α β (Fin k)) (X' : Layout β) (hX' : X'.n + k = X.n) (l : ρ)
    (c0 : ι → ℕ) (hc0 : ∀ i (p : Fin (NU i)), col i p = c0 i) (hkn : ∀ i, k * c0 i ≤ n i)
    (hp1 : Knows m b.B s x (fun i => 2^(k * c0 i)))
    (hp2 : Knows m b.B s x (fun i => 2^(n i - k * c0 i)))
    (hke : Knows m b.B s x (fun i => k * c0 i))
    {ACC : ι → Tape (Tape ℂ)} {SZ : ι → Tape ℕ} {τ : ι → ℕ}
    (hk : ∀ i, (b.use i (mkIn cl (k * c0 i) 0 0 0
        (jobsTape (k * c0 i) (NU i * 2^(n i - k * c0 i))
          (toVec (n i) (k * c0 i) (hkn i) (exJobs X S X' hX' l (h.fit i) (V i))))
        (jobsTape (k * c0 i) 0 (fun _ _ => 0)))).OK
      (jobsTape (k * c0 i) (NU i * 2^(n i - k * c0 i))
        (fun t => ripple (k * c0 i) (toVec (n i) (k * c0 i) (hkn i) (exJobs X S X' hX' l (h.fit i) (V i)) t)),
        ACC i, SZ i) (τ i) (b.ext i)) :
    Able m b s (p (a (a (a (c cl)))) (portOut cl)) x
      (fun i => (ctab X R (NU i) (col i) (ib i)
          (fun p q => ferryM 0 (col i p) l (S.lift (kernel (Fin k))) (V i p q)),
        (jobsTape (k * c0 i) (NU i * 2^(n i - k * c0 i))
          (fun t => ripple (k * c0 i) (toVec (n i) (k * c0 i) (hkn i) (exJobs X S X' hX' l (h.fit i) (V i)) t)),
          ACC i, SZ i)))
      (fun i => NU i * 2^(n i)) τ := by
  -- trivial routing tables: every unit feeds itself
  let TT : Ty := p (a w) (p (a w) (a w))
  let tabs (i : ι) : T TT := ((Layout.std (NU i)).tape (fun t => exFl (NU i) t),
    (Layout.std (NU i)).tape (fun t => id t.val), (Layout.std (NU i)).tape (fun _ => 0))
  have Htab : Can m b.B s TT x tabs NU := by
    have h1 := can_tabIx (t:=w) (fun i (t : Fin (NU i)) => exFl (NU i) t) h.nuK h.nuS
      ((Can.wc (m:=m) (s:=p s w) (B:=fun u : UnitF NU => b.B u.1) (x:=fun u => (x u.1, u.2.val)) 1).cong
        (fun _ => rfl) (fun u => by simp [exFl, u.2.2]))
    have h2 := can_tabIx (t:=w) (fun i (t : Fin (NU i)) => id t.val) h.nuK h.nuS Can.second
    have h3 := can_tabIx (t:=w) (fun i (_ : Fin (NU i)) => 0) h.nuK h.nuS (Can.wc 0)
    exact h1.pair (h2.pair h3)
  let x1 (i : ι) : T (p s TT) := (x i, tabs i)
  have h1 := h.ext tabs
  have F1 : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z → Has m b.B (p s TT) u x1 z :=
    fun g => Can.first.then_do g
  have tfl : Has m b.B (p s TT) (a w) x1 (fun i => (Layout.std (NU i)).tape (fun t => exFl (NU i) t)) :=
    Can.fst Can.second
  have tpp : Has m b.B (p s TT) (a w) x1 (fun i => (Layout.std (NU i)).tape (fun t => id t.val)) :=
    Can.fst (Can.snd Can.second)
  have tsi : Has m b.B (p s TT) (a w) x1 (fun i => (Layout.std (NU i)).tape (fun _ => 0)) :=
    Can.snd (Can.snd Can.second)
  have tsup : Has m b.B (p s TT) (jobsT cl) x1 (fun i => jobsTape (n i) 0 (fun _ _ => 0)) :=
    has_noJobs n
  -- the jobs (one per unit)
  have HG := h1.gather (F1 ct) S X' hX' l (fl:=fun i => exFl (NU i)) (pp:=fun _ => id)
    (si:=fun _ _ => 0) (BS:=fun _ => 0) (sup:=fun _ _ _ => 0) (NJ:=NU)
    (fun i => exFl_pp (NU i)) (fun i t ht hf => by simp [exFl, ht] at hf)
    (F1 h.nuK) h.nuS tfl tpp tsi tsup
  -- `HG` is the tape of `exJobs`
  let JB (i : ι) := jobsTape (n i) (NU i) (exJobs X S X' hX' l (h.fit i) (V i))
  have HG' : Can m b.B (p s TT) (jobsT cl) x1 JB (fun i => NU i * 2^(n i)) := HG
  let x2 (i : ι) : T (p (p s TT) (jobsT cl)) := (x1 i, JB i)
  have F2 : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z → Has m b.B (p (p s TT) (jobsT cl)) u x2 z :=
    fun g => Can.first.then_do (F1 g)
  have hJ2 : Has m b.B (p (p s TT) (jobsT cl)) (jobsT cl) x2 JB := Can.second
  have HV := can_toVec (cl:=cl) (k:=fun i => k * c0 i) (NJ:=NU) hkn hJ2 (F2 h.nuK) (F2 hp2) (F2 hp1) h.big h.nuS
  have Hin : Can m b.B (p (p s TT) (jobsT cl)) (portIn cl) x2
      (fun i => mkIn cl (k * c0 i) 0 0 0
        (jobsTape (k * c0 i) (NU i * 2^(n i - k * c0 i))
          (toVec (n i) (k * c0 i) (hkn i) (exJobs X S X' hX' l (h.fit i) (V i))))
        (jobsTape (k * c0 i) 0 (fun _ _ => 0))) (fun i => NU i * 2^(n i)) :=
    ((F2 hke).mono (fun _ => Nat.zero_le _)).pair (((F2 h.imag).mono (fun _ => Nat.zero_le _)).pair
      ((Can.wc 0).pair ((Can.wc 0).pair ((Can.wc 0).pair (HV.pair
        ((has_noJobs (fun i => k * c0 i)).mono (fun _ => Nat.zero_le _)))))))
  have Hcall := (Able.of_can (b:=b) Hin).comp (Able.use_call (m:=m) (b:=b) hk)
  -- back into jobs, then into the class table
  let outp (i : ι) : T (portOut cl) :=
    (jobsTape (k * c0 i) (NU i * 2^(n i - k * c0 i))
      (fun t => ripple (k * c0 i) (toVec (n i) (k * c0 i) (hkn i) (exJobs X S X' hX' l (h.fit i) (V i)) t)),
      ACC i, SZ i)
  let x3 (i : ι) : T (p (p (p s TT) (jobsT cl)) (portOut cl)) := (x2 i, outp i)
  have F3 : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z → Has m b.B (p (p (p s TT) (jobsT cl)) (portOut cl)) u x3 z :=
    fun g => Can.first.then_do (F2 g)
  have hO3 : Has m b.B (p (p (p s TT) (jobsT cl)) (portOut cl)) (jobsT cl) x3
      (fun i => jobsTape (k * c0 i) (NU i * 2^(n i - k * c0 i))
        (fun t => ripple (k * c0 i) (toVec (n i) (k * c0 i) (hkn i) (exJobs X S X' hX' l (h.fit i) (V i)) t))) :=
    Can.fst Can.second
  have HB := can_fromVec (cl:=cl) (k:=fun i => k * c0 i) (NJ:=NU) hkn hO3 (F3 h.nuK) (F3 hp2) (F3 hp1)
    (F3 h.p2n) h.big h.nuS
  let YB (i : ι) := fromVec (n i) (k * c0 i)
    (fun t => ripple (k * c0 i) (toVec (n i) (k * c0 i) (hkn i) (exJobs X S X' hX' l (h.fit i) (V i)) t))
  let x4 (i : ι) := (x3 i, jobsTape (n i) (NU i) (YB i))
  have h4 := ((h1.ext (t:=jobsT cl) JB).ext (t:=portOut cl) outp).ext (t:=jobsT cl)
    (fun i => jobsTape (n i) (NU i) (YB i))
  have F4 : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z →
      Has m b.B (p (p (p (p s TT) (jobsT cl)) (portOut cl)) (jobsT cl)) u x4 z :=
    fun g => Can.first.then_do (F3 g)
  have tslot : Has m b.B (p (p (p (p s TT) (jobsT cl)) (portOut cl)) (jobsT cl)) (a w) x4
      (fun i => (Layout.std (NU i)).tape (fun p => id p.val)) :=
    Can.first.then_do (Can.first.then_do (Can.first.then_do tpp))
  have hY : ∀ i (p : Fin (NU i)) (q : Bits (ib i p)) (y : Sky α (col i p)),
      YB i (id p.val) (placeS S X' (col i p) (ib i p) (n i) (fitk X' hX' (h4.fit i p)) (q, y)) =
      ferryM 0 (col i p) l (S.lift (kernel (Fin k))) (V i p q) (l, (0, y)) := by
    intro i p q y
    simp only [YB, id]
    rw [fromVec_ripple]
    have hf : exFl (NU i) p.val ≠ 0 := by simp [exFl, p.2]
    have hc := hc0 i p
    have e1 : exJobs X S X' hX' l (h.fit i) (V i) p.val =
        blockJob S X' (col i p) (ib i p) (n i) (fitk X' hX' (h.fit i p)) l (V i p) :=
      cjob_fed X S X' hX' l (h.fit i) (V i) _ id _ _ _ p.val p hf rfl
    rw [e1, ← hc]
    exact blockJob_hlow S X' _ _ _ _ l (V i p) q 0 y
  have HS := h4.scatter (Y:=YB) (F4 ct) S X' hX' l (slot:=fun _ t => id t) (NJ:=NU)
    (fun i p => p.2) tslot Can.second hY
  -- assemble
  have Hpost : Can m b.B (p (p (p s TT) (jobsT cl)) (portOut cl)) (p (a (a (a (c cl)))) (portOut cl)) x3
      (fun i => (ctab X R (NU i) (col i) (ib i)
          (fun p q => ferryM 0 (col i p) l (S.lift (kernel (Fin k))) (V i p q)), outp i))
      (fun i => NU i * 2^(n i)) :=
    HB.bind ((HS.mono (fun _ => le_rfl)).pair (Can.first.then_do Can.second |>.mono (fun _ => Nat.zero_le _)))
  have Hall := ((Able.of_can (b:=b) Htab).save ((Able.of_can (b:=b) HG').save
    (Hcall.save (Able.of_can (b:=b) Hpost))))
  refine Hall.weaken ⟨5, fun i => ?_⟩ (fun i => by simp)
  have h1' : NU i ≤ NU i * 2^(n i) := Nat.le_mul_of_pos_right _ (Nat.two_pow_pos _)
  simp only
  nlinarith

end ExitRAM

end
end PowerSaving.RAM
end OAI
