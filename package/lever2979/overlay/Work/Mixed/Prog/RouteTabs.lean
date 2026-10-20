import Work.Mixed.Prog.ClassTabs

/-!
# Mixed engine, R3 program part 9: the routing tables of a block letter (agent e-r3prog)

For a block letter of rank `r` of a class `(n, h+1, μ)` the class computes (brute force,
polynomial bookkeeping) the index part of its child `(n, h, μ*r)` and the tables the block letter
reads (`Route.lean`): fed flags, parent units, supply indices, the child slot of each parent unit,
and, for the source accumulator, the source slots in order and their sizes.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

universe U

/-- A table of `N` entries, each computed at cost `P`: work `N*(P+1)`. -/
lemma tabC {ι : Type U} {m : Bool} {B : ι → ℕ} {s t : Ty} {x : ι → T s} {N P : ι → ℕ}
    (g : ∀ i, Fin (N i) → T t) (hN : Knows m B s x N) (hNs : Small B N)
    (hg : Can m (fun u : UnitF N => B u.1) (p s w) t (fun u => (x u.1, u.2.val)) (fun u => g u.1 u.2)
      (fun u => P u.1)) :
    Can m B s (a t) x (fun i => (Layout.std (N i)).tape (g i)) (fun i => N i * (P i + 1)) := by
  have H := Can.layout (t:=t) (L:=fun i => Layout.std (N i)) (v:=g) (P:=P) (W:=fun _ => 0) hN hg hNs
  refine H.mono (fun i => ?_)
  change 0 + N i + N i * P i ≤ N i * (P i + 1)
  nlinarith

/-- The child's knowledge from the parent's. -/
lemma CIn.child {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s} {n h μ : ι → ℕ}
    (c : CIn m B s x n (fun i => h i + 1) μ) (r : ℕ) :
    CIn m B s x n h (fun i => μ i * r) :=
  ⟨c.nK, (c.hK.wsub (Can.wc 1) (c.hS.sub)).cong (fun _ => rfl) (fun i => by simp), c.μK.times (Can.wc r)
    c.μS (Small.const _ r), c.p2n, c.imag, c.big, fun i => by have := c.hle i; omega,
    c.μS.mul (Small.const _ r)⟩

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s} {n h μ : ι → ℕ}

lemma parCode_lt (n h μ r t : ℕ) (ht : t < E.NUc n h (μ * r) * E.W) : E.parCode n h μ r t < E.W^(h+1) := by
  unfold parCode
  have hp : t / E.W < E.NUc n h (μ * r) := (Nat.div_lt_iff_lt_mul E.W_pos).2 ht
  exact lt_of_cons E.W h _ _ (Nat.mod_lt _ E.W_pos) (MXR.U_lt hp)

/-- The routing tables of a block letter of rank `r` (child `(n, h, μ*r)` of `(n, h+1, μ)`). -/
abbrev RtT : Ty := p IdxT (p (a w) (p (a w) (p (a w) (p (a w) (p (a w) (a w))))))
def rtVal (n h μ r o : ℕ) : T RtT :=
  (E.idxVal n h (μ * r),
   (Layout.std (E.NUc n h (μ * r) * E.W)).tape (fun t => E.flR n h μ r t),
   (Layout.std (E.NUc n h (μ * r) * E.W)).tape (fun t => E.ppR n h μ r t),
   (Layout.std (E.NUc n h (μ * r) * E.W)).tape (fun t => E.siR n h μ r o t),
   (Layout.std (E.NUc n (h+1) μ)).tape (fun p => E.slotR n h μ r p),
   (Layout.std (Nat.count (E.srcP n h μ r) (E.NUc n h (μ * r) * E.W))).tape
      (fun k => Nat.nth (E.srcP n h μ r) k),
   (Layout.std (Nat.count (E.srcP n h μ r) (E.NUc n h (μ * r) * E.W))).tape
      (fun k => E.slotSize n h (μ * r) (Nat.nth (E.srcP n h μ r) k)))

/-- Cost of the routing tables. -/
def rtCost (n h μ r : ℕ) : ℕ :=
  E.idxCost n h (μ * r) + 8 * (E.W^(h+1) + 1) * (E.W^(h+1) + 1) * (h + 3) + 8

/-- **The routing tables of a block letter.** -/
lemma can_route (c : CIn m B s x n (fun i => h i + 1) μ) (r : ℕ) {o : ι → ℕ} (ho : Knows m B s x o)
    (hos : Small B o)
    (hWp : Knows m B s x (fun i => E.W^(h i + 1)))
    (hNUp : Knows m B s x (fun i => E.NUc (n i) (h i + 1) (μ i)))
    (hPCp : Has m B s (a w) x (fun i => (Layout.std (E.W^(h i + 1))).tape (fun t => E.PCc (n i) (h i + 1) (μ i) t)))
    (hUp : Has m B s (a w) x (fun i => (Layout.std (E.NUc (n i) (h i + 1) (μ i))).tape
      (fun p => E.Uc (n i) (h i + 1) (μ i) p))) :
    Can m B s RtT x (fun i => E.rtVal (n i) (h i) (μ i) r (o i)) (fun i => E.rtCost (n i) (h i) (μ i) r) := by
  have c' := c.child r
  have Hidx := E.can_index c'
  -- environment x1 = (x, child index tables)
  let x1 (i : ι) : T (p s IdxT) := (x i, E.idxVal (n i) (h i) (μ i * r))
  let NU' (i : ι) := E.NUc (n i) (h i) (μ i * r)
  let NJ (i : ι) := NU' i * E.W
  have c1 := c.ext (t:=IdxT) (fun i => E.idxVal (n i) (h i) (μ i * r))
  have F1 : ∀ {t : Ty} {y : ι → T t}, Has m B s t x y → Has m B (p s IdxT) t x1 y :=
    fun g => Can.first.then_do g
  have hNU' : Knows m B (p s IdxT) x1 NU' := Can.fst (Can.snd Can.second)
  have hPC' : Has m B (p s IdxT) (a w) x1
      (fun i => (Layout.std (E.W^(h i))).tape (fun t => E.PCc (n i) (h i) (μ i * r) t)) :=
    Can.fst (Can.snd (Can.snd Can.second))
  have hU' : Has m B (p s IdxT) (a w) x1
      (fun i => (Layout.std (NU' i)).tape (fun p => E.Uc (n i) (h i) (μ i * r) p)) :=
    Can.snd (Can.snd (Can.snd Can.second))
  have hNU's : Small B NU' := (c'.powS E.W).mono (fun i => E.NU_le' _ _ _)
  have hNJ : Knows m B (p s IdxT) x1 NJ := hNU'.times (Can.wc E.W) hNU's (Small.const _ _)
  have hNJs : Small B NJ := hNU's.mul (Small.const _ _)
  have hWps : Small B (fun i => E.W^(h i + 1)) := c.powS E.W
  -- the parent code of a child slot, at any index below NJ
  have pcAt : ∀ {κ : Type U} {sy : Ty} (jx : κ → ι) (y : κ → T sy) (jj : κ → ℕ),
      Has m (fun k => B (jx k)) sy (p s IdxT) y (fun k => x1 (jx k)) →
      Knows m (fun k => B (jx k)) sy y jj → (∀ k, jj k < NJ (jx k)) →
      Knows m (fun k => B (jx k)) sy y (fun k => E.parCode (n (jx k)) (h (jx k)) (μ (jx k)) r (jj k)) := by
    intro κ sy jx y jj hy hj hlt
    have hjs : Small (fun k => B (jx k)) jj := (hNJs.reindex jx).mono (fun k => (hlt k).le)
    have hq : Knows m (fun k => B (jx k)) sy y (fun k => jj k / E.W) := hj.dvdW (Can.wc E.W) hjs
    have hqlt : ∀ k, jj k / E.W < NU' (jx k) := fun k => (Nat.div_lt_iff_lt_mul E.W_pos).2 (hlt k)
    have hUq := Can.fetch (t:=w) (L:=fun k => Layout.std (NU' (jx k)))
      (v:=fun k p => E.Uc (n (jx k)) (h (jx k)) (μ (jx k) * r) p) (i:=fun k => ⟨jj k / E.W, hqlt k⟩)
      (hy.then_do (hU'.reindex jx)) hq
    have hUs : Small (fun k => B (jx k)) (fun k => E.Uc (n (jx k)) (h (jx k)) (μ (jx k) * r) (jj k / E.W)) :=
      ((c'.powS E.W).reindex jx).mono (fun k => (MXR.U_lt (hqlt k)).le)
    exact (hj.modW (Can.wc E.W) hjs).plus ((Can.wc E.W).times hUq (Small.const _ _) hUs) hjs.mod
      ((Small.const _ _).mul hUs)
  have pcAt_s : ∀ {κ : Type U} (jx : κ → ι) (jj : κ → ℕ), (∀ k, jj k < NJ (jx k)) →
      Small (fun k => B (jx k)) (fun k => E.parCode (n (jx k)) (h (jx k)) (μ (jx k)) r (jj k)) := by
    intro κ jx jj hlt
    exact ((hWps.reindex jx)).mono (fun k => (E.parCode_lt _ _ _ _ _ (hlt k)).le)
  -- the region test of the parent at a parent code
  have inAt : ∀ {κ : Type U} {sy : Ty} (jx : κ → ι) (y : κ → T sy) (jj : κ → ℕ),
      Has m (fun k => B (jx k)) sy (p s IdxT) y (fun k => x1 (jx k)) →
      Knows m (fun k => B (jx k)) sy y jj → (∀ k, jj k < NJ (jx k)) →
      Can m (fun k => B (jx k)) sy w y
        (fun k => if inRS E.W E.tq E.mx (n (jx k)) (h (jx k) + 1) (μ (jx k))
          (E.parCode (n (jx k)) (h (jx k)) (μ (jx k)) r (jj k)) then 1 else 0) (fun k => h (jx k) + 1 + 1) := by
    intro κ sy jx y jj hy hj hlt
    have cK := c1.reindex jx y hy
    exact E.knows_inR cK (pcAt jx y jj hy hj hlt) (pcAt_s jx jj hlt) (hy.then_do ((F1 hWp).reindex jx))
  -- slots
  let SL := UnitF NJ
  let yS (u : SL) : T (p (p s IdxT) w) := (x1 u.1, u.2.val)
  have hfl : Can m (fun u : SL => B u.1) (p (p s IdxT) w) w yS
      (fun u => E.flR (n u.1) (h u.1) (μ u.1) r u.2.val) (fun u => h u.1 + 2) := by
    refine (inAt (fun u : SL => u.1) yS (fun u => u.2.val) Can.first Can.second (fun u => u.2.2)).cong
      (fun _ => rfl) (fun u => ?_)
    have h1 : u.2.val < NJ u.1 := u.2.2
    unfold flR
    by_cases h2 : inRS E.W E.tq E.mx (n u.1) (h u.1 + 1) (μ u.1) (E.parCode (n u.1) (h u.1) (μ u.1) r u.2.val)
    · rw [if_pos h2, if_pos (show E.fedP (n u.1) (h u.1) (μ u.1) r u.2.val from ⟨h1, h2⟩)]
    · rw [if_neg h2, if_neg (show ¬ E.fedP (n u.1) (h u.1) (μ u.1) r u.2.val from fun hh => h2 hh.2)]
  have hpp : Knows m (fun u : SL => B u.1) (p (p s IdxT) w) yS (fun u => E.ppR (n u.1) (h u.1) (μ u.1) r u.2.val) :=
    Can.fetch (t:=w) (L:=fun u : SL => Layout.std (E.W^(h u.1 + 1)))
      (v:=fun u t => E.PCc (n u.1) (h u.1 + 1) (μ u.1) t)
      (i:=fun u => ⟨E.parCode (n u.1) (h u.1) (μ u.1) r u.2.val, E.parCode_lt _ _ _ _ _ u.2.2⟩)
      (Can.first.then_do ((F1 hPCp).reindex (fun u : SL => u.1)))
      (pcAt (fun u : SL => u.1) yS (fun u => u.2.val) Can.first Can.second (fun u => u.2.2))
  -- supply index: the running offset plus the source slots before
  have hsrc : ∀ {κ : Type U} {sy : Ty} (jx : κ → ι) (y : κ → T sy) (jj : κ → ℕ),
      Has m (fun k => B (jx k)) sy (p s IdxT) y (fun k => x1 (jx k)) →
      Knows m (fun k => B (jx k)) sy y jj → (∀ k, jj k < NJ (jx k)) →
      Can m (fun k => B (jx k)) sy w y
        (fun k => if E.srcP (n (jx k)) (h (jx k)) (μ (jx k)) r (jj k) then 1 else 0) (fun k => h (jx k) + 2) := by
    intro κ sy jx y jj hy hj hlt
    refine ((Can.wc 1).wsub (inAt jx y jj hy hj hlt) ((Small.const _ 1).mono (fun _ => Nat.sub_le _ _))).cong
      (fun _ => rfl) (fun k => ?_)
    by_cases h2 : inRS E.W E.tq E.mx (n (jx k)) (h (jx k) + 1) (μ (jx k)) (E.parCode (n (jx k)) (h (jx k)) (μ (jx k)) r (jj k))
    · rw [if_pos h2, if_neg (show ¬ E.srcP (n (jx k)) (h (jx k)) (μ (jx k)) r (jj k) from fun hh => hh.2 h2)]
    · rw [if_neg h2, if_pos (show E.srcP (n (jx k)) (h (jx k)) (μ (jx k)) r (jj k) from ⟨hlt k, h2⟩)]
  have hsi : Can m (fun u : SL => B u.1) (p (p s IdxT) w) w yS
      (fun u => E.siR (n u.1) (h u.1) (μ u.1) r (o u.1) u.2.val) (fun u => NJ u.1 * (h u.1 + 3) + 2) := by
    let DD := During (fun u : SL => u.2.val)
    have hyD : Has m (fun d : DD => B d.i.1) (p (p (p s IdxT) w) w) (p s IdxT) (fun d => (yS d.i, d.j))
        (fun d => x1 d.i.1) := Can.first.then_do Can.first
    have hjD : Knows m (fun d : DD => B d.i.1) (p (p (p s IdxT) w) w) (fun d => (yS d.i, d.j)) (fun d => d.j) :=
      Can.second
    have hcnt := can_count (P:=fun u : SL => E.srcP (n u.1) (h u.1) (μ u.1) r) (Wt:=fun u => h u.1 + 2)
      (T:=fun u : SL => u.2.val) (x:=yS) Can.second ((hNJs.reindex (fun u : SL => u.1)).mono (fun u => u.2.2.le))
      (hsrc (sy:=p (p (p s IdxT) w) w) (fun d : DD => d.i.1) (fun d => (yS d.i, d.j)) (fun d => d.j)
        hyD hjD (fun d => lt_trans d.hj d.i.2.2))
    refine ((Can.first.then_do ((F1 ho).reindex (fun u : SL => u.1))).mono (fun _ => Nat.zero_le _) |>.plus
      hcnt ((hos.reindex (fun u : SL => u.1))) (((hNJs.reindex (fun u : SL => u.1))).mono
        (fun u => le_trans (MXR.count_le_self _) u.2.2.le))).mono (fun u => ?_)
    have := u.2.2
    change u.2.val * (h u.1 + 2 + 1) + 1 ≤ NJ u.1 * (h u.1 + 3) + 2
    nlinarith
  -- child slot of each parent unit
  let NUp (i : ι) := E.NUc (n i) (h i + 1) (μ i)
  let UP := UnitF NUp
  let yP (u : UP) : T (p (p s IdxT) w) := (x1 u.1, u.2.val)
  have hUpp : Knows m (fun u : UP => B u.1) (p (p s IdxT) w) yP (fun u => E.Uc (n u.1) (h u.1 + 1) (μ u.1) u.2) :=
    unit_fetch (m:=m) (B:=B) (s:=p s IdxT) (t:=w) (NU:=NUp) (x:=x1)
      (fun i (p : Fin (NUp i)) => E.Uc (n i) (h i + 1) (μ i) p) (F1 hUp)
  have hUps : Small (fun u : UP => B u.1) (fun u => E.Uc (n u.1) (h u.1 + 1) (μ u.1) u.2) :=
    ((hWps.reindex (fun u : UP => u.1))).mono (fun u => (MXR.U_lt u.2.2).le)
  have hτlt : ∀ u : UP, E.Uc (n u.1) (h u.1 + 1) (μ u.1) u.2 / E.W < E.W^(h u.1) := fun u => by
    rw [Nat.div_lt_iff_lt_mul E.W_pos, ← pow_succ]; exact MXR.U_lt u.2.2
  have hPCτ := Can.fetch (t:=w) (L:=fun u : UP => Layout.std (E.W^(h u.1)))
    (v:=fun u t => E.PCc (n u.1) (h u.1) (μ u.1 * r) t)
    (i:=fun u => ⟨E.Uc (n u.1) (h u.1 + 1) (μ u.1) u.2 / E.W, hτlt u⟩)
    (Can.first.then_do (hPC'.reindex (fun u : UP => u.1))) (hUpp.dvdW (Can.wc E.W) hUps)
  have hPCτs : Small (fun u : UP => B u.1) (fun u => E.PCc (n u.1) (h u.1) (μ u.1 * r)
      (E.Uc (n u.1) (h u.1 + 1) (μ u.1) u.2 / E.W)) :=
    ((c'.powS E.W).reindex (fun u : UP => u.1)).mono (fun u => le_trans (MXR.count_le_self _) (hτlt u).le)
  have hslot : Knows m (fun u : UP => B u.1) (p (p s IdxT) w) yP (fun u => E.slotR (n u.1) (h u.1) (μ u.1) r u.2) :=
    (hPCτ.times (Can.wc E.W) hPCτs (Small.const _ _)).plus (hUpp.modW (Can.wc E.W) hUps)
      (hPCτs.mul (Small.const _ _)) hUps.mod
  have hNUps : Small B NUp := hWps.mono (fun i => E.NU_le' _ _ _)
  have Tslot := tabC (t:=w) (P:=fun _ => 0) (fun i (p : Fin (NUp i)) => E.slotR (n i) (h i) (μ i) r p)
    (F1 hNUp) hNUps hslot
  -- tables over slots
  have Tfl := tabC (t:=w) (P:=fun i => h i + 2) (fun i (t : Fin (NJ i)) => E.flR (n i) (h i) (μ i) r t)
    hNJ hNJs hfl
  have Tpp := tabC (t:=w) (P:=fun _ => 0) (fun i (t : Fin (NJ i)) => E.ppR (n i) (h i) (μ i) r t) hNJ hNJs hpp
  have Tsi := tabC (t:=w) (P:=fun i => NJ i * (h i + 3) + 2)
    (fun i (t : Fin (NJ i)) => E.siR (n i) (h i) (μ i) r (o i) t) hNJ hNJs hsi
  -- number of source slots
  let ns (i : ι) := Nat.count (E.srcP (n i) (h i) (μ i) r) (NJ i)
  have hns : Can m B (p s IdxT) w x1 ns (fun i => NJ i * (h i + 2 + 1) + 1) := by
    let DD := During NJ
    have hyD : Has m (fun d : DD => B d.i) (p (p s IdxT) w) (p s IdxT) (fun d => (x1 d.i, d.j))
        (fun d => x1 d.i) := Can.first
    have hjD : Knows m (fun d : DD => B d.i) (p (p s IdxT) w) (fun d => (x1 d.i, d.j)) (fun d => d.j) := Can.second
    exact can_count (P:=fun i => E.srcP (n i) (h i) (μ i) r) (Wt:=fun i => h i + 2) hNJ hNJs
      (hsrc (sy:=p (p s IdxT) w) (fun d : DD => d.i) (fun d => (x1 d.i, d.j)) (fun d => d.j) hyD hjD
        (fun d => d.hj))
  have hnsle : ∀ i, ns i ≤ NJ i := fun i => MXR.count_le_self _
  have hnss : Small B ns := hNJs.mono hnsle
  -- environment x2 = (x1, ns)
  let x2 (i : ι) : T (p (p s IdxT) w) := (x1 i, ns i)
  let SK := UnitF ns
  let yK (u : SK) : T (p (p (p s IdxT) w) w) := (x2 u.1, u.2.val)
  have hnth : Can m (fun u : SK => B u.1) (p (p (p s IdxT) w) w) w yK
      (fun u => Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2.val) (fun u => NJ u.1 * (h u.1 + 2 + 1) + 1) := by
    let DD := During (fun u : SK => NJ u.1)
    have hyD : Has m (fun d : DD => B d.i.1) (p (p (p (p s IdxT) w) w) w) (p s IdxT) (fun d => (yK d.i, d.j))
        (fun d => x1 d.i.1) := Can.first.then_do (Can.first.then_do Can.first)
    have hjD : Knows m (fun d : DD => B d.i.1) (p (p (p (p s IdxT) w) w) w) (fun d => (yK d.i, d.j))
        (fun d => d.j) := Can.second
    have H := can_nthc (P:=fun u : SK => E.srcP (n u.1) (h u.1) (μ u.1) r) (Wt:=fun u => h u.1 + 2)
      (T:=fun u : SK => NJ u.1) (q:=fun u => u.2.val) (x:=yK)
      (Can.first.then_do (Can.first.then_do (hNJ.reindex (fun u : SK => u.1))))
      (hNJs.reindex (fun u : SK => u.1)) Can.second
      (hsrc (sy:=p (p (p (p s IdxT) w) w) w) (fun d : DD => d.i.1) (fun d => (yK d.i, d.j)) (fun d => d.j)
        hyD hjD (fun d => d.hj))
    refine H.cong (fun _ => rfl) (fun u => ?_)
    have : u.2.val < Nat.count (E.srcP (n u.1) (h u.1) (μ u.1) r) (NJ u.1) := u.2.2
    rw [if_pos this]
  have hnthlt : ∀ u : SK, Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2.val < NJ u.1 := fun u =>
    (MXR.nth_spec u.2.2).2.1
  have hnths : Small (fun u : SK => B u.1) (fun u => Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2.val) :=
    (hNJs.reindex (fun u : SK => u.1)).mono (fun u => (hnthlt u).le)
  -- the size of the k-th source slot
  have hsz : Can m (fun u : SK => B u.1) (p (p (p s IdxT) w) w) w yK
      (fun u => E.slotSize (n u.1) (h u.1) (μ u.1 * r) (Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2.val))
      (fun u => NJ u.1 * (h u.1 + 3) + h u.1 + 2) := by
    let tt (u : SK) := Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2.val
    have hsave := hnth.save (t:=w) (u:=w) (z:=fun u => E.slotSize (n u.1) (h u.1) (μ u.1 * r) (tt u))
      (W':=fun u => h u.1 + 1) (by
        let yT (u : SK) : T (p (p (p (p s IdxT) w) w) w) := (yK u, tt u)
        have ht : Knows m (fun u : SK => B u.1) (p (p (p (p s IdxT) w) w) w) yT tt := Can.second
        have hq : Knows m (fun u : SK => B u.1) (p (p (p (p s IdxT) w) w) w) yT (fun u => tt u / E.W) :=
          ht.dvdW (Can.wc E.W) hnths
        have hqlt : ∀ u : SK, tt u / E.W < NU' u.1 := fun u => (Nat.div_lt_iff_lt_mul E.W_pos).2 (hnthlt u)
        have hUq := Can.fetch (t:=w) (L:=fun u : SK => Layout.std (NU' u.1))
          (v:=fun u p => E.Uc (n u.1) (h u.1) (μ u.1 * r) p) (i:=fun u => ⟨tt u / E.W, hqlt u⟩)
          (Can.first.then_do (Can.first.then_do (Can.first.then_do (hU'.reindex (fun u : SK => u.1))))) hq
        have hUs : Small (fun u : SK => B u.1) (fun u => E.Uc (n u.1) (h u.1) (μ u.1 * r) (tt u / E.W)) :=
          ((c'.powS E.W).reindex (fun u : SK => u.1)).mono (fun u => (MXR.U_lt (hqlt u)).le)
        have cK := c'.reindex (fun u : SK => u.1) (t:=p (p (p (p s IdxT) w) w) w) yT
          (Can.first.then_do (Can.first.then_do (Can.first.then_do Can.first)))
        have hcol := E.can_colS cK.hK cK.hS hUq hUs (cK.powS E.mx)
        have hcs : Small (fun u : SK => B u.1) (fun u => colS E.W E.tq (h u.1) (E.Uc (n u.1) (h u.1) (μ u.1 * r) (tt u / E.W))) :=
          (cK.powS E.mx).mono (fun u => (E.stepC_iter_le _ _).1)
        have hcU := (cK.μK.mono (fun _ => Nat.zero_le _)).times hcol cK.μS hcs
        have htq := E.knows_tq (ht.modW (Can.wc E.W) hnths) (fun u => Nat.mod_lt _ E.W_pos)
        refine (htq.mono (fun _ => Nat.zero_le _) |>.times hcU ((Small.const _ E.mx).mono (fun u => E.tq_le _))
          (cK.μS.mul hcs)).cong (fun _ => rfl) (fun u => ?_)
        rfl)
    exact hsave.mono (fun u => by have := u.2.2; nlinarith)
  have hnsK : Knows m B (p (p s IdxT) w) x2 ns := Can.second
  have TN := tabC (t:=w) (P:=fun i => NJ i * (h i + 2 + 1) + 1)
    (fun i (k : Fin (ns i)) => Nat.nth (E.srcP (n i) (h i) (μ i) r) k) hnsK hnss hnth
  have TS := tabC (t:=w) (P:=fun i => NJ i * (h i + 3) + h i + 2)
    (fun i (k : Fin (ns i)) => E.slotSize (n i) (h i) (μ i * r) (Nat.nth (E.srcP (n i) (h i) (μ i) r) k))
    hnsK hnss hsz
  have Tns := hns.save (TN.fork TS)
  have hidx1 : Has m B (p s IdxT) IdxT x1 (fun i => E.idxVal (n i) (h i) (μ i * r)) := Can.second
  have Hrest := hidx1.fork (Tfl.fork (Tpp.fork (Tsi.fork (Tslot.fork Tns))))
  have H := Hidx.save Hrest
  refine H.weaken ⟨20, fun i => ?_⟩
  have e1 : NJ i ≤ E.W^(h i + 1) := by
    have := E.NU_le' (n i) (h i) (μ i * r)
    calc NJ i = NU' i * E.W := rfl
      _ ≤ E.W^(h i) * E.W := Nat.mul_le_mul_right _ this
      _ = E.W^(h i + 1) := (pow_succ _ _).symm
  have e2 : NUp i ≤ E.W^(h i + 1) := E.NU_le' _ _ _
  have e3 : ns i ≤ NJ i := hnsle i
  simp only [rtCost]
  set X := E.W^(h i + 1) with hX
  set Y := E.idxCost (n i) (h i) (μ i * r)
  change Y + (0 + (NJ i * (h i + 2 + 1) + (NJ i * (0 + 1) + (NJ i * (NJ i * (h i + 3) + 2 + 1) +
    (NUp i * (0 + 1) + (NJ i * (h i + 2 + 1) + 1 + (ns i * (NJ i * (h i + 2 + 1) + 1 + 1) +
      ns i * (NJ i * (h i + 3) + h i + 2 + 1)))))))) + 1 ≤ 20 * (Y + 8 * (X + 1) * (X + 1) * (h i + 3) + 8 + 1)
  have f1 : NJ i * (NJ i * (h i + 3) + 2 + 1) ≤ X * X * (h i + 3) + 3 * X := by
    have := Nat.mul_le_mul e1 (Nat.add_le_add_right (Nat.mul_le_mul_right (h i + 3) e1) 3)
    calc _ = NJ i * (NJ i * (h i + 3) + 3) := by ring
      _ ≤ X * (X * (h i + 3) + 3) := this
      _ = _ := by ring
  have f2 : ns i * (NJ i * (h i + 2 + 1) + 1 + 1) ≤ X * X * (h i + 3) + 2 * X := by
    have := Nat.mul_le_mul (le_trans e3 e1) (Nat.add_le_add_right (Nat.mul_le_mul_right (h i + 3) e1) 2)
    calc _ = ns i * (NJ i * (h i + 3) + 2) := by ring
      _ ≤ X * (X * (h i + 3) + 2) := this
      _ = _ := by ring
  have f3 : ns i * (NJ i * (h i + 3) + h i + 2 + 1) ≤ X * X * (h i + 3) + X * (h i + 3) := by
    have := Nat.mul_le_mul (le_trans e3 e1) (Nat.add_le_add_right (Nat.mul_le_mul_right (h i + 3) e1) (h i + 3))
    calc _ = ns i * (NJ i * (h i + 3) + (h i + 3)) := by ring
      _ ≤ X * (X * (h i + 3) + (h i + 3)) := this
      _ = _ := by ring
  have f4 : NJ i * (h i + 2 + 1) ≤ X * (h i + 3) := Nat.mul_le_mul_right _ e1
  have f5 : NJ i * (0 + 1) ≤ X := by simpa using e1
  have f6 : NUp i * (0 + 1) ≤ X := by simpa using e2
  have f7 : 8 * (X + 1) * (X + 1) * (h i + 3) = 8 * (X * X * (h i + 3)) + 16 * (X * (h i + 3)) + 8 * (h i + 3) := by
    ring
  have f8 : X ≤ X * (h i + 3) := Nat.le_mul_of_pos_right _ (by omega)
  rw [f7]
  omega


end Eng

end
end PowerSaving.RAM
end OAI
