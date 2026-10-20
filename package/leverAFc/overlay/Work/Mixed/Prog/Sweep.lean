import Work.Mixed.Prog.SweepMath

/-!
# Mixed engine, R3 program part 14: one letter of a class with its source accumulator (agent e-r3prog)

The sweep state of a class: `(class table, consumed answers, consumed sizes, running offset)`.

* `Eng.block_step`  a block letter of a class of height `h+1`: routing tables, ONE child call, read back,
                    accumulate the child's own sources and its consumed list.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

universe U

/-- The sweep state of a class. -/
abbrev StT (cl : Paint) : Ty := p (a (a (a (c cl)))) (p (jobsT cl) (p (a w) w))

/-- Answered consumed supply vectors (sizes `A`, from offset `o`). -/
def accT (n o : ℕ) (sup : ℕ → Bits n → ℂ) (A : List ℕ) : Tape (Tape ℂ) :=
  jobsTape n A.length (fun k => hlow n (A.getD k 0) (sup (o + k)))
/-- Their sizes. -/
def szT (A : List ℕ) : Tape ℕ := (Layout.std A.length).tape (fun k => A.getD k 0)

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)
variable {ι : Type U} {m : Bool} {cl : Paint} {b : Bank (portM cl) ι} {s : Ty} {x : ι → T s}
variable {n h μ o BS : ι → ℕ} {sup : ∀ i, ℕ → Bits (n i) → ℂ}

/-- The class table of a class family. -/
abbrev ctabE (n h μ : ℕ) (V : E.CState n h μ) : Tape (Tape (Tape ℂ)) :=
  ctab E.X E.R (E.NUc n h μ) (fun p => E.cUc n h μ p) (fun p => E.iUc n h μ p) V

/-- The accumulated part of one block letter of rank `r`. -/
def accR (n h μ r : ℕ) : List ℕ := E.ownSrc n h μ r ++ E.srcSizes n h (μ * r)

lemma tape_std_ext {μ' : Type*} {N N' : ℕ} (hN : N = N') (f : Fin N → μ') (g : Fin N' → μ')
    (hfg : ∀ k (h1 : k < N) (h2 : k < N'), f ⟨k, h1⟩ = g ⟨k, h2⟩) :
    (Layout.std N).tape f = (Layout.std N').tape g := by
  subst hN
  unfold Layout.tape
  congr 1
  funext j
  have hj : (Layout.std N).atFin j = j := by
    apply Fin.ext; exact Layout.loc_atFin (Layout.std N) j
  rw [hj]
  exact hfg j.val j.2 j.2

/-- **The accumulator after one child.** -/
lemma acc_cat (n o : ℕ) (sup : ℕ → Bits n → ℂ) (A own Ac : List ℕ) (Y' : ℕ → Bits n → ℂ) (nthS : ℕ → ℕ)
    (hY : ∀ k, k < own.length → Y' (nthS k) = hlow n (own.getD k 0) (sup (o + A.length + k))) :
    (Layout.std (A.length + own.length + Ac.length)).tape (fun k => cat3 A.length own.length Ac.length
      (fun k => (bits n).tape (hlow n (A.getD k 0) (sup (o + k))))
      (fun k => (bits n).tape (Y' (nthS k)))
      (fun k => (bits n).tape (hlow n (Ac.getD k 0) (sup (o + A.length + own.length + k)))) k) =
    accT n o sup (A ++ own ++ Ac) := by
  unfold accT jobsTape
  apply tape_std_ext (by simp [Nat.add_assoc])
  intro k h1 h2
  simp only [cat3]
  by_cases ha : k < A.length
  · rw [if_pos ha, List.append_assoc, List.getD_append _ _ _ _ ha]
  · rw [if_neg ha]
    by_cases hb : k < A.length + own.length
    · rw [if_pos hb, hY (k - A.length) (by omega), List.append_assoc,
        List.getD_append_right _ _ _ _ (by omega), List.getD_append _ _ _ _ (by omega)]
      congr 3; omega
    · rw [if_neg hb, List.append_assoc, List.getD_append_right _ _ _ _ (by omega),
        List.getD_append_right _ _ _ _ (by omega)]
      congr 3
      omega

lemma sz_cat (A own Ac : List ℕ) :
    (Layout.std (A.length + own.length + Ac.length)).tape (fun k => cat3 A.length own.length Ac.length
      (fun k => A.getD k 0) (fun k => own.getD k 0) (fun k => Ac.getD k 0) k) = szT (A ++ own ++ Ac) := by
  unfold szT
  apply tape_std_ext (by simp [Nat.add_assoc])
  intro k h1 h2
  simp only [cat3]
  by_cases ha : k < A.length
  · rw [if_pos ha, List.append_assoc, List.getD_append _ _ _ _ ha]
  · rw [if_neg ha]
    by_cases hb : k < A.length + own.length
    · rw [if_pos hb, List.append_assoc, List.getD_append_right _ _ _ _ (by omega),
        List.getD_append _ _ _ _ (by omega)]
    · rw [if_neg hb, List.append_assoc, List.getD_append_right _ _ _ _ (by omega),
        List.getD_append_right _ _ _ _ (by omega)]

/-- Cost of one block letter of a class of height `h+1`. -/
def blkCost (n h μ r L : ℕ) : ℕ :=
  E.rtCost n h μ r + (E.NUc n h (μ * r) * E.W + E.NUc n (h+1) μ) * 2^n + 2 * (L + (E.accR n h μ r).length) + 20

/-- **One block letter of a class of height `h+1`**: routing tables, ONE child call, read-back,
accumulation of the child's sources. -/
theorem block_step
    (cin : CIn m b.B s x n (fun i => h i + 1) μ)
    (cc : CCtx m b.B s x E.X.n n (fun i => E.NUc (n i) (h i + 1) (μ i)) (fun i p => E.cUc (n i) (h i + 1) (μ i) p)
      (fun i p => E.iUc (n i) (h i + 1) (μ i) p))
    (hWp : Knows m b.B s x (fun i => E.W^(h i + 1)))
    (hPCp : Has m b.B s (a w) x (fun i => (Layout.std (E.W^(h i + 1))).tape (fun t => E.PCc (n i) (h i + 1) (μ i) t)))
    (hUp : Has m b.B s (a w) x (fun i => (Layout.std (E.NUc (n i) (h i + 1) (μ i))).tape
      (fun p => E.Uc (n i) (h i + 1) (μ i) p)))
    (tsup : Has m b.B s (jobsT cl) x (fun i => jobsTape (n i) (BS i) (sup i))) (hBSs : Small b.B BS)
    {V : ∀ i, E.CState (n i) (h i + 1) (μ i)} {A : ι → List ℕ}
    (ct : Has m b.B s (a (a (a (c cl)))) x (fun i => E.ctabE (n i) (h i + 1) (μ i) (V i)))
    (tacc : Has m b.B s (jobsT cl) x (fun i => accT (n i) (o i) (sup i) (A i)))
    (tsz : Has m b.B s (a w) x (fun i => szT (A i)))
    (horun : Knows m b.B s x (fun i => o i + (A i).length))
    (l : ρ) (r : ℕ) (z : Fin r → Space α) (ind : LinearIndependent F z) (hr : r ≤ E.X.n)
    (hsep : ∀ l', r < E.tr l')
    (hfit : ∀ i, o i + (A i).length + (E.accR (n i) (h i) (μ i) r).length ≤ BS i)
    {τ : ι → ℕ}
    (hk : ∀ i, (b.use i (E.enc cl (E.childC (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) l r z ind hr))).OK
      (E.out cl (E.childC (n i) (h i) (μ i) (BS i) (sup i) (V i) (o i + (A i).length) l r z ind hr)) (τ i) (b.ext i)) :
    Able m b s (StT cl) x
      (fun i => (E.ctabE (n i) (h i + 1) (μ i) (E.stepV (V i) (.old (.block l r z ind))),
        accT (n i) (o i) (sup i) (A i ++ E.accR (n i) (h i) (μ i) r), szT (A i ++ E.accR (n i) (h i) (μ i) r),
        o i + (A i ++ E.accR (n i) (h i) (μ i) r).length))
      (fun i => E.blkCost (n i) (h i) (μ i) r (A i).length) τ := by
  let orun (i : ι) := o i + (A i).length
  have hos : Small b.B orun := hBSs.mono (fun i => by have := hfit i; simp only [orun]; omega)
  -- 1. routing tables
  have HR := E.can_route cin r horun hos hWp cc.nuK hPCp hUp
  let rt (i : ι) := E.rtVal (n i) (h i) (μ i) r (orun i)
  let x1 (i : ι) : T (p s RtT) := (x i, rt i)
  have F1 : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z → Has m b.B (p s RtT) u x1 z :=
    fun g => Can.first.then_do g
  have cc1 := cc.ext rt
  have hR : Has m b.B (p s RtT) RtT x1 rt := Can.second
  let NU' (i : ι) := E.NUc (n i) (h i) (μ i * r)
  let NJ (i : ι) := NU' i * E.W
  have hNU' : Knows m b.B (p s RtT) x1 NU' := Can.fst (Can.snd (Can.fst hR))
  have hNU's : Small b.B NU' := (cin.child r).powS E.W |>.mono (fun i => E.NU_le' _ _ _)
  have hNJ : Knows m b.B (p s RtT) x1 NJ := hNU'.times (Can.wc E.W) hNU's (Small.const _ _)
  have hNJs : Small b.B NJ := hNU's.mul (Small.const _ _)
  have tfl : Has m b.B (p s RtT) (a w) x1
      (fun i => (Layout.std (NJ i)).tape (fun t => E.flR (n i) (h i) (μ i) r t)) := Can.fst (Can.snd hR)
  have tpp : Has m b.B (p s RtT) (a w) x1
      (fun i => (Layout.std (NJ i)).tape (fun t => E.ppR (n i) (h i) (μ i) r t)) := Can.fst (Can.snd (Can.snd hR))
  have tsi : Has m b.B (p s RtT) (a w) x1
      (fun i => (Layout.std (NJ i)).tape (fun t => E.siR (n i) (h i) (μ i) r (orun i) t)) :=
    Can.fst (Can.snd (Can.snd (Can.snd hR)))
  have tslot : Has m b.B (p s RtT) (a w) x1
      (fun i => (Layout.std (E.NUc (n i) (h i + 1) (μ i))).tape (fun p => E.slotR (n i) (h i) (μ i) r p)) :=
    Can.fst (Can.snd (Can.snd (Can.snd (Can.snd hR))))
  let ns (i : ι) := Nat.count (E.srcP (n i) (h i) (μ i) r) (NJ i)
  have tnth : Has m b.B (p s RtT) (a w) x1
      (fun i => (Layout.std (ns i)).tape (fun k => Nat.nth (E.srcP (n i) (h i) (μ i) r) k)) :=
    Can.fst (Can.snd (Can.snd (Can.snd (Can.snd (Can.snd hR)))))
  have tssz : Has m b.B (p s RtT) (a w) x1
      (fun i => (Layout.std (ns i)).tape (fun k => E.slotSize (n i) (h i) (μ i * r) (Nat.nth (E.srcP (n i) (h i) (μ i) r) k))) :=
    Can.snd (Can.snd (Can.snd (Can.snd (Can.snd (Can.snd hR)))))
  have hnsle : ∀ i, ns i = (E.ownSrc (n i) (h i) (μ i) r).length := fun i => (E.ownSrc_len _ _ _ _).symm
  have hnss : Small b.B ns := hNJs.mono (fun i => MXR.count_le_self _)
  have hns : Knows m b.B (p s RtT) x1 ns := tnth.len (hnss.mono (fun i => le_of_eq (by simp [Layout.tape, Layout.std])))
  -- 2. the block letter (ONE child call)
  have hpp : ∀ i t, E.flR (n i) (h i) (μ i) r t ≠ 0 → E.ppR (n i) (h i) (μ i) r t < E.NUc (n i) (h i + 1) (μ i) :=
    fun i t => E.ppR_lt _ _ _ _ t
  have hsi : ∀ i t, t < NJ i → E.flR (n i) (h i) (μ i) r t = 0 → E.siR (n i) (h i) (μ i) r (orun i) t < BS i := by
    intro i t ht hf
    have hsrc := (E.fl_zero_iff _ _ _ _ t ht).1 hf
    have h1 : Nat.count (E.srcP (n i) (h i) (μ i) r) t < ns i := count_lt_of_mem hsrc ht
    have h2 := hfit i
    have h3 : (E.ownSrc (n i) (h i) (μ i) r).length ≤ (E.accR (n i) (h i) (μ i) r).length := by
      simp [accR]
    unfold siR
    rw [hnsle] at h1
    simp only [orun]; omega
  have hslot := fun i (p : Fin (E.NUc (n i) (h i + 1) (μ i))) => E.slotR_spec (n i) (h i) (μ i) r hsep p
  let S := bS z ind
  let X' := bX (α:=α) r
  have hX' : X'.n + r = E.X.n := E.bX_n hr
  have hkd : Knows m b.B (p s RtT) x1 (fun i => h i + 1) := F1 cin.hK
  have hμ' : Knows m b.B (p s RtT) x1 (fun i => μ i * r) := F1 (cin.μK.times (Can.wc r) cin.μS (Small.const _ _))
  have ho' : Knows m b.B (p s RtT) x1 (fun i => orun i + ns i) := (F1 horun).plus hns hos hnss
  have hk' : ∀ i, (b.use i (mkIn cl (n i) (h i + 1) (μ i * r) (orun i + ns i)
      (jobsTape (n i) (NJ i) (cjob E.X S X' hX' l (cc1.fit i) (V i) (E.flR (n i) (h i) (μ i) r)
        (E.ppR (n i) (h i) (μ i) r) (E.siR (n i) (h i) (μ i) r (orun i)) (hpp i) (sup i)))
      (jobsTape (n i) (BS i) (sup i)))).OK
      (jobsTape (n i) (NJ i) (fun t => hlow (n i) (E.slotSize (n i) (h i) (μ i * r) t)
          (cjob E.X S X' hX' l (cc1.fit i) (V i) (E.flR (n i) (h i) (μ i) r) (E.ppR (n i) (h i) (μ i) r)
            (E.siR (n i) (h i) (μ i) r (orun i)) (hpp i) (sup i) t)),
        jobsTape (n i) (E.cons (n i) (h i) (μ i * r))
          (fun k => hlow (n i) (E.csz (n i) (h i) (μ i * r) k) (sup i (orun i + ns i + k))),
        (Layout.std (E.cons (n i) (h i) (μ i * r))).tape (fun k => E.csz (n i) (h i) (μ i * r) k))
      (τ i) (b.ext i) := by
    intro i
    have := hk i
    rw [hnsle i]
    exact this
  have HB := cc1.block_letter (F1 ct) S X' hX' l (NJ:=NJ) (BS:=BS) (kd:=fun i => h i + 1)
    (μ:=fun i => μ i * r) (o:=fun i => orun i + ns i) (fl:=fun i => E.flR (n i) (h i) (μ i) r)
    (pp:=fun i => E.ppR (n i) (h i) (μ i) r) (si:=fun i => E.siR (n i) (h i) (μ i) r (orun i))
    (slot:=fun i => E.slotR (n i) (h i) (μ i) r) (sz:=fun i => E.slotSize (n i) (h i) (μ i * r))
    hpp hsi hslot hNJ hNJs tfl tpp tsi tslot (F1 tsup) hkd hμ' ho' hk'
  -- 3. the accumulator, in the environment holding the call's answer
  let CJ (i : ι) := cjob E.X S X' hX' l (cc1.fit i) (V i) (E.flR (n i) (h i) (μ i) r) (E.ppR (n i) (h i) (μ i) r)
    (E.siR (n i) (h i) (μ i) r (orun i)) (hpp i) (sup i)
  let Y' (i : ι) (t : ℕ) := hlow (n i) (E.slotSize (n i) (h i) (μ i * r) t) (CJ i t)
  let ACC' (i : ι) := jobsTape (n i) (E.cons (n i) (h i) (μ i * r))
    (fun k => hlow (n i) (E.csz (n i) (h i) (μ i * r) k) (sup i (orun i + ns i + k)))
  let SZ' (i : ι) := (Layout.std (E.cons (n i) (h i) (μ i * r))).tape (fun k => E.csz (n i) (h i) (μ i * r) k)
  let ct' (i : ι) := ctab E.X E.R (E.NUc (n i) (h i + 1) (μ i)) (fun p => E.cUc (n i) (h i + 1) (μ i) p)
    (fun p => E.iUc (n i) (h i + 1) (μ i) p)
    (fun p q => ferryM 0 (E.cUc (n i) (h i + 1) (μ i) p) l (S.lift (kernel (Fin r))) (V i p q))
  let R2 (i : ι) : T (p (a (a (a (c cl)))) (portOut cl)) := (ct' i, jobsTape (n i) (NJ i) (Y' i), ACC' i, SZ' i)
  let x2 (i : ι) : T (p (p s RtT) (p (a (a (a (c cl)))) (portOut cl))) := (x1 i, R2 i)
  have F2 : ∀ {u : Ty} {z : ι → T u}, Has m b.B (p s RtT) u x1 z →
      Has m b.B (p (p s RtT) (p (a (a (a (c cl)))) (portOut cl))) u x2 z := fun g => Can.first.then_do g
  have hR2 : Has m b.B (p (p s RtT) (p (a (a (a (c cl)))) (portOut cl))) (p (a (a (a (c cl)))) (portOut cl)) x2 R2 :=
    Can.second
  -- the child's own source slots, answered, in slot order
  let SK := UnitF ns
  let yK (u : SK) : T (p (p (p s RtT) (p (a (a (a (c cl)))) (portOut cl))) w) := (x2 u.1, u.2.val)
  have hnthK : Knows m (fun u : SK => b.B u.1) (p (p (p s RtT) (p (a (a (a (c cl)))) (portOut cl))) w) yK
      (fun u => Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2) := unit_fetch _ (F2 tnth)
  have hnthlt : ∀ u : SK, Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2 < NJ u.1 := fun u =>
    (MXR.nth_spec u.2.2).2.1
  have hOS := Can.fetch (t:=a (c cl)) (L:=fun u : SK => Layout.std (NJ u.1))
    (v:=fun u t => (bits (n u.1)).tape (Y' u.1 t))
    (i:=fun u => (⟨Nat.nth (E.srcP (n u.1) (h u.1) (μ u.1) r) u.2, hnthlt u⟩ : Fin (NJ u.1)))
    (Can.first.then_do ((Can.fst (Can.snd hR2)).reindex (fun u : SK => u.1))) hnthK
  have TOS := tabC (t:=a (c cl)) (P:=fun _ => 0)
    (fun i (k : Fin (ns i)) => (bits (n i)).tape (Y' i (Nat.nth (E.srcP (n i) (h i) (μ i) r) k)))
    (F2 hns) hnss hOS
  let OS (i : ι) := (Layout.std (ns i)).tape (fun k : Fin (ns i) =>
    (bits (n i)).tape (Y' i (Nat.nth (E.srcP (n i) (h i) (μ i) r) k)))
  let x3 (i : ι) := (x2 i, OS i)
  let s3 : Ty := p (p (p s RtT) (p (a (a (a (c cl)))) (portOut cl))) (jobsT cl)
  have F3 : ∀ {u : Ty} {z : ι → T u}, Has m b.B s u x z → Has m b.B s3 u x3 z :=
    fun g => Can.first.then_do (F2 (F1 g))
  have hR3 : Has m b.B s3 (p (a (a (a (c cl)))) (portOut cl)) x3 R2 := Can.first.then_do hR2
  -- lengths
  have hAs : Small b.B (fun i => (A i).length) := hBSs.mono (fun i => by have := hfit i; omega)
  have hconsle : ∀ i, E.cons (n i) (h i) (μ i * r) ≤ BS i := fun i => by
    have := hfit i; simp [accR, cons] at *; omega
  have hcs : Small b.B (fun i => E.cons (n i) (h i) (μ i * r)) := hBSs.mono hconsle
  have kA : Knows m b.B s3 x3 (fun i => (A i).length) := (F3 tacc).len (hAs.mono (fun i => le_of_eq rfl))
  have kns : Knows m b.B s3 x3 ns := Can.first.then_do (F2 hns)
  have kc : Knows m b.B s3 x3 (fun i => E.cons (n i) (h i) (μ i * r)) :=
    (Can.fst (Can.snd (Can.snd hR3))).len (hcs.mono (fun i => le_of_eq rfl))
  have HACC := concat3 (t:=a (c cl)) (L1:=fun i => (A i).length) (L2:=ns) (L3:=fun i => E.cons (n i) (h i) (μ i * r))
    (f1:=fun i k => (bits (n i)).tape (hlow (n i) ((A i).getD k 0) (sup i (o i + k))))
    (f2:=fun i k => (bits (n i)).tape (Y' i (Nat.nth (E.srcP (n i) (h i) (μ i) r) k)))
    (f3:=fun i k => (bits (n i)).tape (hlow (n i) (E.csz (n i) (h i) (μ i * r) k) (sup i (orun i + ns i + k))))
    (F3 tacc) Can.second (Can.fst (Can.snd (Can.snd hR3))) kA kns kc hAs hnss hcs
  have HSZ := concat3 (t:=w) (L1:=fun i => (A i).length) (L2:=ns) (L3:=fun i => E.cons (n i) (h i) (μ i * r))
    (f1:=fun i k => (A i).getD k 0)
    (f2:=fun i k => E.slotSize (n i) (h i) (μ i * r) (Nat.nth (E.srcP (n i) (h i) (μ i) r) k))
    (f3:=fun i k => E.csz (n i) (h i) (μ i * r) k)
    (F3 tsz) (Can.first.then_do (F2 tssz)) (Can.snd (Can.snd (Can.snd hR3))) kA kns kc hAs hnss hcs
  have hor' : Knows m b.B s3 x3 (fun i => orun i + ns i + E.cons (n i) (h i) (μ i * r)) :=
    ((F3 horun).plus kns hos hnss).plus kc (hos.add hnss) hcs
  have HST := (Can.fst hR3).fork (HACC.fork (HSZ.fork hor'))
  -- 4. the values
  have eY : ∀ i k, k < (E.ownSrc (n i) (h i) (μ i) r).length →
      Y' i (Nat.nth (E.srcP (n i) (h i) (μ i) r) k) =
        hlow (n i) ((E.ownSrc (n i) (h i) (μ i) r).getD k 0) (sup i (o i + (A i).length + k)) := by
    intro i k hk
    rw [← hnsle i] at hk
    obtain ⟨hP, hlt, hcnt⟩ := MXR.nth_spec hk
    have hf0 : E.flR (n i) (h i) (μ i) r (Nat.nth (E.srcP (n i) (h i) (μ i) r) k) = 0 :=
      (E.fl_zero_iff _ _ _ _ _ hlt).2 hP
    simp only [Y', CJ]
    unfold cjob
    rw [dif_neg (fun hh => hh hf0)]
    rw [E.ownSrc_eq, getD_map_range _ _ _ hk]
    unfold siR
    rw [hcnt]
  have eACC : ∀ i, (Layout.std ((A i).length + ns i + E.cons (n i) (h i) (μ i * r))).tape (fun k =>
      cat3 (A i).length (ns i) (E.cons (n i) (h i) (μ i * r))
        (fun k => (bits (n i)).tape (hlow (n i) ((A i).getD k 0) (sup i (o i + k))))
        (fun k => (bits (n i)).tape (Y' i (Nat.nth (E.srcP (n i) (h i) (μ i) r) k)))
        (fun k => (bits (n i)).tape (hlow (n i) (E.csz (n i) (h i) (μ i * r) k) (sup i (orun i + ns i + k)))) k) =
      accT (n i) (o i) (sup i) (A i ++ E.accR (n i) (h i) (μ i) r) := by
    intro i
    have := acc_cat (n i) (o i) (sup i) (A i) (E.ownSrc (n i) (h i) (μ i) r) (E.srcSizes (n i) (h i) (μ i * r))
      (Y' i) (Nat.nth (E.srcP (n i) (h i) (μ i) r)) (eY i)
    show _ = accT (n i) (o i) (sup i) (A i ++ (E.ownSrc (n i) (h i) (μ i) r ++ E.srcSizes (n i) (h i) (μ i * r)))
    rw [← List.append_assoc, ← this, hnsle i]
    rfl
  have eSZ : ∀ i, (Layout.std ((A i).length + ns i + E.cons (n i) (h i) (μ i * r))).tape (fun k =>
      cat3 (A i).length (ns i) (E.cons (n i) (h i) (μ i * r)) (fun k => (A i).getD k 0)
        (fun k => E.slotSize (n i) (h i) (μ i * r) (Nat.nth (E.srcP (n i) (h i) (μ i) r) k))
        (fun k => E.csz (n i) (h i) (μ i * r) k) k) = szT (A i ++ E.accR (n i) (h i) (μ i) r) := by
    intro i
    show _ = szT (A i ++ (E.ownSrc (n i) (h i) (μ i) r ++ E.srcSizes (n i) (h i) (μ i * r)))
    rw [← List.append_assoc, ← sz_cat]
    apply tape_std_ext (by rw [hnsle i]; rfl)
    intro k h1 h2
    simp only [cat3]
    by_cases ha : k < (A i).length
    · rw [if_pos ha, if_pos ha]
    · rw [if_neg ha, if_neg ha]
      by_cases hb : k < (A i).length + ns i
      · have hns_def : ns i = Nat.count (E.srcP (n i) (h i) (μ i) r) (E.NUc (n i) (h i) (μ i * r) * E.W) := rfl
        rw [if_pos hb, if_pos (by rw [← hnsle i]; exact hb), E.ownSrc_eq, getD_map_range _ _ _ (by omega)]
      · rw [if_neg hb, if_neg (by rw [← hnsle i]; exact hb), hnsle i]
        rfl
  have ect : ∀ i, ct' i = E.ctabE (n i) (h i + 1) (μ i) (E.stepV (V i) (.old (.block l r z ind))) := by
    intro i; rw [E.stepV_block]
  have eor : ∀ i, orun i + ns i + E.cons (n i) (h i) (μ i * r) = o i + (A i ++ E.accR (n i) (h i) (μ i) r).length := by
    intro i
    show o i + (A i).length + ns i + (E.srcSizes (n i) (h i) (μ i * r)).length =
      o i + (A i ++ (E.ownSrc (n i) (h i) (μ i) r ++ E.srcSizes (n i) (h i) (μ i * r))).length
    rw [List.length_append, List.length_append, hnsle i]; omega
  have Hfin : Can m b.B s3 (StT cl) x3
      (fun i => (E.ctabE (n i) (h i + 1) (μ i) (E.stepV (V i) (.old (.block l r z ind))),
        accT (n i) (o i) (sup i) (A i ++ E.accR (n i) (h i) (μ i) r), szT (A i ++ E.accR (n i) (h i) (μ i) r),
        o i + (A i ++ E.accR (n i) (h i) (μ i) r).length))
      (fun i => 0 + ((A i).length + ns i + E.cons (n i) (h i) (μ i * r) +
        ((A i).length + ns i + E.cons (n i) (h i) (μ i * r) + 0))) :=
    HST.cong (fun _ => rfl) (fun i => Prod.ext (ect i) (Prod.ext (eACC i) (Prod.ext (eSZ i) (eor i))))
  have Hpost := TOS.save Hfin
  have H := (Able.of_can (b:=b) HR).save (HB.save (Able.of_can (b:=b) Hpost))
  refine H.weaken ⟨4, fun i => ?_⟩ (fun i => by simp)
  have hns' := hnsle i
  have hc' : E.cons (n i) (h i) (μ i * r) = (E.srcSizes (n i) (h i) (μ i * r)).length := rfl
  have hacc : (E.accR (n i) (h i) (μ i) r).length = ns i + E.cons (n i) (h i) (μ i * r) := by
    rw [hns', hc']; simp [accR]
  have hbk : E.blkCost (n i) (h i) (μ i) r (A i).length = E.rtCost (n i) (h i) (μ i) r +
      (NJ i + E.NUc (n i) (h i + 1) (μ i)) * 2^(n i) + 2 * ((A i).length + (ns i + E.cons (n i) (h i) (μ i * r))) + 20 := by
    unfold blkCost; rw [hacc]
  dsimp only
  rw [hbk, Nat.add_mul]
  generalize NJ i * 2^(n i) = X1
  generalize E.NUc (n i) (h i + 1) (μ i) * 2^(n i) = X2
  omega

/-- What the letters of a word add to the accumulator. -/
def accW (n h μ : ℕ) (ws : GWord α ρ) : List ℕ := (blockRanks ws).flatMap (fun r => E.accR n h μ r)

lemma accW_cons (n h μ : ℕ) (mv : GMove α ρ) (ws : GWord α ρ) :
    E.accW n h μ (mv :: ws) = (match mv.rk? with | some r => E.accR n h μ r | none => []) ++ E.accW n h μ ws := by
  unfold accW blockRanks
  rw [List.filterMap_cons]
  cases mv.rk? <;> simp

/-- Sweep cost of a word for a class of height `h+1` (accumulator length `L` before it). -/
def swCost (n h μ : ℕ) : ℕ → GWord α ρ → ℕ
  | _, [] => 0
  | L, mv :: ws => (match mv.rk? with
      | some r => E.blkCost n h μ r L
      | none => E.NUc n (h+1) μ * 2^n) +
      swCost n h μ (L + (match mv.rk? with | some r => (E.accR n h μ r).length | none => 0)) ws

lemma runV_cons {n h μ : ℕ} (V : E.CState n h μ) (mv : GMove α ρ) (ws : GWord α ρ) :
    E.runV V (mv :: ws) = E.runV (E.stepV V mv) ws := by
  funext p q; unfold runV stepV; exact gcoast_cons 0 _ mv ws (V p q)

lemma runV_nil {n h μ : ℕ} (V : E.CState n h μ) : E.runV V [] = V := by
  funext p q x; rfl

end Eng

end
end PowerSaving.RAM
end OAI
