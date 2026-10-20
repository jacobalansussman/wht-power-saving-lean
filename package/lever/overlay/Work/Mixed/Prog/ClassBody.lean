import Work.Mixed.Prog.SweepLeaf

/-!
# Mixed engine, R3 program part 17: the body of a class call (agent e-r3prog)

From the port input `(n, I, h+1, μ, o, own jobs, supply)`: `2^n`, the class tables, load, the sweep of
the whole word (child calls or exits), unload, answer `(own answers, consumed answers, sizes)`.

* `Eng.class_up`    a class of height `h+1` (children);
* `Eng.class_leaf`  a class of height `0` (exits).
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

universe U

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)
variable {ι : Type U} {m : Bool} {cl : Paint} {b : Bank (portM cl) ι}
variable {n h μ o BS : ι → ℕ} {J sup : ∀ i, ℕ → Bits (n i) → ℂ} {τ : ι → MDesc → ℕ}

/-- The empty tapes are all equal. -/
lemma jobsTape_zero (n : ℕ) (f g : ℕ → Bits n → ℂ) : jobsTape n 0 f = jobsTape n 0 g := by
  unfold jobsTape; exact tape_std_ext rfl _ _ (fun k h1 _ => absurd h1 (Nat.not_lt_zero _))

/-- The own answers after the word are `hlow` of the own jobs. -/
lemma own_answers (n h μ : ℕ) (J : ℕ → Bits n → ℂ) :
    jobsTape n (E.NUc n h μ * E.W) (E.unl n h μ (E.runV (E.V0 n h μ J) E.wd)) =
      jobsTape n (E.NUc n h μ * E.W) (fun t => hlow n (E.slotSize n h μ t) (J t)) := by
  unfold jobsTape
  apply tape_std_ext rfl
  intro k h1 _
  show (bits n).tape (E.unl n h μ (E.runV (E.V0 n h μ J) E.wd) k) = (bits n).tape (hlow n (E.slotSize n h μ k) (J k))
  rw [E.unl_run n h μ J k h1]

/-- The decoded port input of a class. -/
lemma decode {x : ι → T (portIn cl)}
    (hx : ∀ i, x i = mkIn cl (n i) (h i + 1) (μ i) (o i) (jobsTape (n i) (E.NUc (n i) (h i) (μ i) * E.W) (J i))
      (jobsTape (n i) (BS i) (sup i))) :
    Knows m b.B (portIn cl) x n ∧ Has m b.B (portIn cl) sc x (fun _ => Complex.I) ∧
    Knows m b.B (portIn cl) x (fun i => h i + 1) ∧ Knows m b.B (portIn cl) x μ ∧ Knows m b.B (portIn cl) x o ∧
    Has m b.B (portIn cl) (jobsT cl) x (fun i => jobsTape (n i) (E.NUc (n i) (h i) (μ i) * E.W) (J i)) ∧
    Has m b.B (portIn cl) (jobsT cl) x (fun i => jobsTape (n i) (BS i) (sup i)) := by
  refine ⟨Can.first.cong (fun _ => rfl) (fun i => by rw [hx i]; rfl),
    (Can.fst Can.second).cong (fun _ => rfl) (fun i => by rw [hx i]; rfl),
    (Can.fst (Can.snd Can.second)).cong (fun _ => rfl) (fun i => by rw [hx i]; rfl),
    (Can.fst (Can.snd (Can.snd Can.second))).cong (fun _ => rfl) (fun i => by rw [hx i]; rfl),
    (Can.fst (Can.snd (Can.snd (Can.snd Can.second)))).cong (fun _ => rfl) (fun i => by rw [hx i]; rfl),
    (Can.fst (Can.snd (Can.snd (Can.snd (Can.snd Can.second))))).cong (fun _ => rfl) (fun i => by rw [hx i]; rfl),
    (Can.snd (Can.snd (Can.snd (Can.snd (Can.snd Can.second))))).cong (fun _ => rfl) (fun i => by rw [hx i]; rfl)⟩

/-- Local work of a class of height `h+1`. -/
def locUp (n h μ : ℕ) : ℕ :=
  (n + 1) + E.idxCost n (h+1) μ + E.NUc n (h+1) μ * (4 * n + 4 * (h + 1) + 8) +
    E.NUc n (h+1) μ * (2 * 2^n + 2) + E.NUc n (h+1) μ * 2^n + E.swCost n h μ 0 E.wd +
    E.NUc n (h+1) μ * E.W * 2^n

lemma szT_nil : szT [] = (Layout.std 0).tape (fun _ : Fin 0 => (0 : ℕ)) := by
  unfold szT; exact tape_std_ext rfl _ _ (fun k h1 _ => absurd h1 (Nat.not_lt_zero _))

/-- **The body of a class of height `h+1`.** -/
theorem class_up (hB2 : Small b.B (fun i => 2^(n i))) (hBS : Small b.B BS)
    (hle : ∀ i, h i + 1 ≤ n i) (hμ : Small b.B μ)
    (hfit : ∀ i, o i + E.cons (n i) (h i + 1) (μ i) ≤ BS i)
    (hk : ∀ i c, c ∈ E.kidsC (MDesc.cls (n i) (h i + 1) (μ i) (o i) (BS i) (J i) (sup i)) →
      (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) :
    Able m b (portIn cl) (portOut cl)
      (fun i => E.enc cl (MDesc.cls (n i) (h i + 1) (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.out cl (MDesc.cls (n i) (h i + 1) (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.locUp (n i) (h i) (μ i))
      (fun i => ((E.kidsC (MDesc.cls (n i) (h i + 1) (μ i) (o i) (BS i) (J i) (sup i))).map (τ i)).sum) := by
  let x0 (i : ι) := E.enc cl (MDesc.cls (n i) (h i + 1) (μ i) (o i) (BS i) (J i) (sup i))
  obtain ⟨hn, hI, hkd, hμK, hoK, hown, hsup⟩ := E.decode (b:=b) (x:=x0) (n:=n) (h:=fun i => h i + 1) (μ:=μ)
    (o:=o) (BS:=BS) (J:=J) (sup:=sup) (fun i => rfl)
  -- 2^n
  have hnS : Small b.B n := hB2.mono (fun i => Nat.lt_two_pow_self.le)
  have H2n : Can m b.B (portIn cl) w x0 (fun i => 2^(n i)) (fun i => n i + 1) :=
    (hn.mono (fun _ => Nat.zero_le _)).wTwo hnS hB2 (Dom.of_le (fun i => by omega))
  let x1 (i : ι) : T (p (portIn cl) w) := (x0 i, 2^(n i))
  have c1 : CIn m b.B (p (portIn cl) w) x1 n (fun i => h i + 1) μ :=
    ⟨Can.first.then_do hn, Can.first.then_do ((hkd.wsub (Can.wc 1) (hnS.mono (fun i => by have := hle i; omega))).cong
      (fun _ => rfl) (fun i => by simp)), Can.first.then_do hμK, Can.second, Can.first.then_do hI, hB2,
      hle, hμ⟩
  have Hidx := E.can_index c1
  have Hun := E.can_units c1
  have Hge := E.can_gears c1
  have cc := E.cctx c1
  let idx (i : ι) := E.idxVal (n i) (h i + 1) (μ i)
  let un (i : ι) := E.unitVal (n i) (h i + 1) (μ i)
  let ge (i : ι) := E.gearVal (n i) (h i + 1) (μ i)
  let s4 : Ty := ClsT (p (portIn cl) w)
  let x4 (i : ι) : T s4 := (((x1 i, idx i), un i), ge i)
  have F4 : ∀ {u : Ty} {z : ι → T u}, Has m b.B (portIn cl) u x0 z → Has m b.B s4 u x4 z :=
    fun g => Can.first.then_do (Can.first.then_do (Can.first.then_do (Can.first.then_do g)))
  have hidx4 : Has m b.B s4 IdxT x4 idx := Can.first.then_do (Can.first.then_do Can.second)
  have c4 : CIn m b.B s4 x4 n (fun i => h i + 1) μ :=
    ((c1.ext (t:=IdxT) idx).ext (t:=UnT) un).ext (t:=GearT) ge
  -- load
  have Hld := E.load cc (J:=J) (F4 hown)
  let ct0 (i : ι) : T (a (a (a (c cl)))) := E.ctabE (n i) (h i + 1) (μ i) (E.V0 (n i) (h i + 1) (μ i) (J i))
  let s5 : Ty := p s4 (a (a (a (c cl))))
  let x5 (i : ι) : T s5 := (x4 i, ct0 i)
  have F5 : ∀ {u : Ty} {z : ι → T u}, Has m b.B s4 u x4 z → Has m b.B s5 u x5 z := fun g => Can.first.then_do g
  -- the sweep
  have tacc : Has m b.B s5 (jobsT cl) x5 (fun i => accT (n i) (o i) (sup i) ([] : List ℕ)) :=
    (has_noJobs (cl:=cl) n).cong (fun _ => rfl) (fun i => jobsTape_zero _ _ _)
  have tsz : Has m b.B s5 (a w) x5 (fun i => szT ([] : List ℕ)) := by
    have := can_tabIx (t:=w) (N:=fun _ => 0) (x:=x5) (B:=b.B) (m:=m) (fun i (t : Fin 0) => (0 : ℕ)) (Can.wc 0)
      (Small.const _ 0) (can_of_empty (fun u => Nat.not_lt_zero _ u.2.2))
    exact this.cong (fun _ => rfl) (fun i => szT_nil.symm)
  have horun : Knows m b.B s5 x5 (fun i => o i + ([] : List ℕ).length) := F5 (F4 hoK)
  have Hsw := E.sweepUp (b:=b) (n:=n) (h:=h) (μ:=μ) (o:=o) (BS:=BS) (sup:=sup) (τ:=τ) E.wd
    (s:=s5) (x:=x5) (V:=fun i => E.V0 (n i) (h i + 1) (μ i) (J i)) (A:=fun _ => [])
    (c4.ext (t:=a (a (a (c cl)))) ct0) (cc.ext (t:=a (a (a (c cl)))) ct0)
    (F5 (Can.fst hidx4)) (F5 (Can.fst (Can.snd (Can.snd hidx4)))) (F5 (Can.snd (Can.snd (Can.snd hidx4))))
    (F5 (F4 hsup)) hBS Can.second tacc tsz horun E.hP E.sep
    (fun i => by
      have h1 := hfit i
      have e : (E.accW (n i) (h i) (μ i) E.wd).length = E.cons (n i) (h i + 1) (μ i) := rfl
      simp only [List.length_nil, Nat.add_zero]; omega)
    (fun i c hc => hk i c hc)
  -- unload and answer
  let stF (i : ι) : T (StT cl) := (E.ctabE (n i) (h i + 1) (μ i) (E.runV (E.V0 (n i) (h i + 1) (μ i) (J i)) E.wd),
      accT (n i) (o i) (sup i) ([] ++ E.accW (n i) (h i) (μ i) E.wd),
      szT ([] ++ E.accW (n i) (h i) (μ i) E.wd), o i + ([] ++ E.accW (n i) (h i) (μ i) E.wd).length)
  let s6 : Ty := p s5 (StT cl)
  let x6 (i : ι) : T s6 := (x5 i, stF i)
  have cc6 := (cc.ext (t:=a (a (a (c cl)))) ct0).ext (t:=StT cl) stF
  have hst : Has m b.B s6 (StT cl) x6 stF := Can.second
  have HU := E.unload cc6 (V:=fun i => E.runV (E.V0 (n i) (h i + 1) (μ i) (J i)) E.wd) (Can.fst hst)
  have Hout : Can m b.B s6 (portOut cl) x6
      (fun i => E.out cl (MDesc.cls (n i) (h i + 1) (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.NUc (n i) (h i + 1) (μ i) * E.W * 2^(n i)) := by
    refine ((HU.fork ((Can.fst (Can.snd hst)).fork (Can.fst (Can.snd (Can.snd hst))))).cong (fun _ => rfl)
      (fun i => ?_)).mono (fun i => by simp)
    show (_, _, _) = (_, _, _)
    rw [E.own_answers]
    rfl
  have H := (Able.of_can (b:=b) H2n).save ((Able.of_can (b:=b) Hidx).save ((Able.of_can (b:=b) Hun).save
    ((Able.of_can (b:=b) Hge).save ((Able.of_can (b:=b) Hld).save (Hsw.save (Able.of_can (b:=b) Hout))))))
  refine H.mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
  · simp only [locUp, List.length_nil]; ring
  · simp only [Nat.zero_add, Nat.add_zero, List.length_nil]
    rfl

/-- Local work of a leaf class. -/
def locLeaf (n μ : ℕ) : ℕ :=
  (n + 1) + E.idxCost n 0 μ + E.NUc n 0 μ * (4 * n + 4 * 0 + 8) +
    E.NUc n 0 μ * (2 * 2^n + 2) + E.NUc n 0 μ * 2^n + E.lfCost n μ E.wd + E.NUc n 0 μ * E.W * 2^n

/-- **The body of a leaf class.** -/
theorem class_leaf (hB2 : Small b.B (fun i => 2^(n i))) (hμ : Small b.B μ) (hfit : ∀ i, E.X.n * μ i ≤ n i)
    (hk : ∀ i c, c ∈ E.kidsC (MDesc.cls (n i) 0 (μ i) (o i) (BS i) (J i) (sup i)) →
      (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) :
    Able m b (portIn cl) (portOut cl)
      (fun i => E.enc cl (MDesc.cls (n i) 0 (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.out cl (MDesc.cls (n i) 0 (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.locLeaf (n i) (μ i))
      (fun i => ((E.kidsC (MDesc.cls (n i) 0 (μ i) (o i) (BS i) (J i) (sup i))).map (τ i)).sum) := by
  let x0 (i : ι) := E.enc cl (MDesc.cls (n i) 0 (μ i) (o i) (BS i) (J i) (sup i))
  obtain ⟨hn, hI, hkd, hμK, hoK, hown, hsup⟩ := E.decode (b:=b) (x:=x0) (n:=n) (h:=fun _ => 0) (μ:=μ)
    (o:=o) (BS:=BS) (J:=J) (sup:=sup) (fun i => rfl)
  have hnS : Small b.B n := hB2.mono (fun i => Nat.lt_two_pow_self.le)
  have H2n : Can m b.B (portIn cl) w x0 (fun i => 2^(n i)) (fun i => n i + 1) :=
    (hn.mono (fun _ => Nat.zero_le _)).wTwo hnS hB2 (Dom.of_le (fun i => by omega))
  let x1 (i : ι) : T (p (portIn cl) w) := (x0 i, 2^(n i))
  have c1 : CIn m b.B (p (portIn cl) w) x1 n (fun _ => 0) μ :=
    ⟨Can.first.then_do hn, Can.wc 0, Can.first.then_do hμK, Can.second, Can.first.then_do hI, hB2,
      fun i => Nat.zero_le _, hμ⟩
  have Hidx := E.can_index c1
  have Hun := E.can_units c1
  have Hge := E.can_gears c1
  have cc := E.cctx c1
  let idx (i : ι) := E.idxVal (n i) 0 (μ i)
  let un (i : ι) := E.unitVal (n i) 0 (μ i)
  let ge (i : ι) := E.gearVal (n i) 0 (μ i)
  let s4 : Ty := ClsT (p (portIn cl) w)
  let x4 (i : ι) : T s4 := (((x1 i, idx i), un i), ge i)
  have F4 : ∀ {u : Ty} {z : ι → T u}, Has m b.B (portIn cl) u x0 z → Has m b.B s4 u x4 z :=
    fun g => Can.first.then_do (Can.first.then_do (Can.first.then_do (Can.first.then_do g)))
  have c4 : CIn m b.B s4 x4 n (fun _ => 0) μ :=
    ((c1.ext (t:=IdxT) idx).ext (t:=UnT) un).ext (t:=GearT) ge
  have Hld := E.load cc (J:=J) (F4 hown)
  let ct0 (i : ι) : T (a (a (a (c cl)))) := E.ctabE (n i) 0 (μ i) (E.V0 (n i) 0 (μ i) (J i))
  let s5 : Ty := p s4 (a (a (a (c cl))))
  let x5 (i : ι) : T s5 := (x4 i, ct0 i)
  have Hsw := E.sweepLeaf (b:=b) (n:=n) (μ:=μ) (τ:=τ) E.wd (s:=s5) (x:=x5)
    (V:=fun i => E.V0 (n i) 0 (μ i) (J i)) (c4.ext (t:=a (a (a (c cl)))) ct0) (cc.ext (t:=a (a (a (c cl)))) ct0)
    Can.second E.hP hfit (fun i c hc => hk i c hc)
  let ctF (i : ι) : T (a (a (a (c cl)))) := E.ctabE (n i) 0 (μ i) (E.runV (E.V0 (n i) 0 (μ i) (J i)) E.wd)
  let s6 : Ty := p s5 (a (a (a (c cl))))
  let x6 (i : ι) : T s6 := (x5 i, ctF i)
  have cc6 := (cc.ext (t:=a (a (a (c cl)))) ct0).ext (t:=a (a (a (c cl)))) ctF
  have HU := E.unload cc6 (V:=fun i => E.runV (E.V0 (n i) 0 (μ i) (J i)) E.wd) Can.second
  have hz : Has m b.B s6 (jobsT cl) x6 (fun i => jobsTape (n i) 0 (fun _ _ => 0)) := has_noJobs n
  have hz2 : Has m b.B s6 (a w) x6 (fun i => (Layout.std 0).tape (fun _ : Fin 0 => (0 : ℕ))) :=
    can_tabIx (t:=w) (N:=fun _ => 0) (x:=x6) (B:=b.B) (m:=m) (fun i (t : Fin 0) => (0 : ℕ)) (Can.wc 0)
      (Small.const _ 0) (can_of_empty (fun u => Nat.not_lt_zero _ u.2.2))
  have Hout : Can m b.B s6 (portOut cl) x6
      (fun i => E.out cl (MDesc.cls (n i) 0 (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.NUc (n i) 0 (μ i) * E.W * 2^(n i)) := by
    refine ((HU.fork (hz.fork hz2)).cong (fun _ => rfl) (fun i => ?_)).mono (fun i => by simp)
    show (_, _, _) = (_, _, _)
    rw [E.own_answers]
    refine Prod.ext rfl (Prod.ext ?_ ?_)
    · dsimp only
      unfold jobsTape
      exact tape_std_ext (show 0 = E.cons (n i) 0 (μ i) from rfl) _ _ (fun k h1 _ => absurd h1 (Nat.not_lt_zero _))
    · dsimp only
      exact tape_std_ext (show 0 = E.cons (n i) 0 (μ i) from rfl) _ _ (fun k h1 _ => absurd h1 (Nat.not_lt_zero _))
  have H := (Able.of_can (b:=b) H2n).save ((Able.of_can (b:=b) Hidx).save ((Able.of_can (b:=b) Hun).save
    ((Able.of_can (b:=b) Hge).save ((Able.of_can (b:=b) Hld).save (Hsw.save (Able.of_can (b:=b) Hout))))))
  refine H.mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
  · simp only [locLeaf]; ring
  · simp only [Nat.zero_add, Nat.add_zero]
    rfl

/-- Local work of a class. -/
def clsLoc (n h μ : ℕ) : ℕ := if h = 0 then E.locLeaf n μ else E.locUp n (h - 1) μ

/-- **The body of a class call** (any height): a leaf class makes exits, a higher class makes child calls. -/
theorem class_body (hB2 : Small b.B (fun i => 2^(n i))) (hBS : Small b.B BS)
    (hv1 : ∀ i, h i = 0 → E.X.n * μ i ≤ n i) (hv2 : ∀ i, h i ≤ n i) (hv3 : ∀ i, μ i ≤ n i)
    (hv4 : ∀ i, o i + E.cons (n i) (h i) (μ i) ≤ BS i)
    (hk : ∀ i c, c ∈ E.kidsC (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i)) →
      (b.use i (E.enc cl c)).OK (E.out cl c) (τ i c) (b.ext i)) :
    Able m b (portIn cl) (portOut cl)
      (fun i => E.enc cl (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.out cl (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i)))
      (fun i => E.clsLoc (n i) (h i) (μ i))
      (fun i => ((E.kidsC (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i))).map (τ i)).sum) := by
  let x0 (i : ι) := E.enc cl (MDesc.cls (n i) (h i) (μ i) (o i) (BS i) (J i) (sup i))
  have hkd : Knows m b.B (portIn cl) x0 (fun i => h i + 1) := Can.fst (Can.snd Can.second)
  have hnS : Small b.B n := hB2.mono (fun i => Nat.lt_two_pow_self.le)
  have hh : Knows m b.B (portIn cl) x0 h :=
    (hkd.wsub (Can.wc 1) (hnS.mono (fun i => by have := hv2 i; omega))).cong (fun _ => rfl) (fun i => by simp)
  have hμs : Small b.B μ := hnS.mono hv3
  apply Able.branch (hh.mono (W':=fun i => E.clsLoc (n i) (h i) (μ i)) (fun _ => Nat.zero_le _))
    (fun i => h i = 0) (fun i => Iff.rfl)
  · -- leaf classes
    let K := {i // h i = 0}
    have H := E.class_leaf (m:=m) (b:=b.comap (fun i : K => i.1)) (n:=fun i : K => n i.1) (μ:=fun i : K => μ i.1)
      (o:=fun i : K => o i.1) (BS:=fun i : K => BS i.1) (J:=fun i : K => J i.1) (sup:=fun i : K => sup i.1)
      (τ:=fun i : K => τ i.1) (hB2.reindex _) (hμs.reindex _) (fun i => hv1 i.1 i.2)
      (fun i c hc => hk i.1 c (by rw [i.2]; exact hc))
    refine (H.cong (fun i => ?_) (fun i => ?_)).mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
    · show E.enc cl (MDesc.cls (n i.1) 0 _ _ _ _ _) = E.enc cl (MDesc.cls (n i.1) (h i.1) _ _ _ _ _)
      rw [i.2]
    · show E.out cl (MDesc.cls (n i.1) 0 _ _ _ _ _) = E.out cl (MDesc.cls (n i.1) (h i.1) _ _ _ _ _)
      rw [i.2]
    · show E.locLeaf (n i.1) (μ i.1) = E.clsLoc (n i.1) (h i.1) (μ i.1)
      unfold clsLoc; rw [if_pos i.2]
    · show _ = ((E.kidsC (MDesc.cls (n i.1) (h i.1) (μ i.1) (o i.1) (BS i.1) (J i.1) (sup i.1))).map (τ i.1)).sum
      rw [i.2]
  · -- classes of height h'+1
    let K := {i // ¬ h i = 0}
    have hsucc : ∀ i : K, h i.1 - 1 + 1 = h i.1 := fun i => by have := i.2; omega
    have H := E.class_up (m:=m) (b:=b.comap (fun i : K => i.1)) (n:=fun i : K => n i.1) (h:=fun i : K => h i.1 - 1)
      (μ:=fun i : K => μ i.1) (o:=fun i : K => o i.1) (BS:=fun i : K => BS i.1) (J:=fun i : K => J i.1)
      (sup:=fun i : K => sup i.1) (τ:=fun i : K => τ i.1) (hB2.reindex _) (hBS.reindex _)
      (fun i => by rw [hsucc i]; exact hv2 i.1) (hμs.reindex _)
      (fun i => by rw [hsucc i]; exact hv4 i.1)
      (fun i c hc => hk i.1 c (by rw [← hsucc i]; exact hc))
    refine (H.cong (fun i => ?_) (fun i => ?_)).mono (fun i => le_of_eq ?_) (fun i => le_of_eq ?_)
    · show E.enc cl (MDesc.cls (n i.1) (h i.1 - 1 + 1) _ _ _ _ _) = E.enc cl (MDesc.cls (n i.1) (h i.1) _ _ _ _ _)
      rw [hsucc i]
    · show E.out cl (MDesc.cls (n i.1) (h i.1 - 1 + 1) _ _ _ _ _) = E.out cl (MDesc.cls (n i.1) (h i.1) _ _ _ _ _)
      rw [hsucc i]
    · show E.locUp (n i.1) (h i.1 - 1) (μ i.1) = E.clsLoc (n i.1) (h i.1) (μ i.1)
      unfold clsLoc; rw [if_neg i.2]
    · show _ = ((E.kidsC (MDesc.cls (n i.1) (h i.1) (μ i.1) (o i.1) (BS i.1) (J i.1) (sup i.1))).map (τ i.1)).sum
      rw [hsucc i]

end Eng

end
end PowerSaving.RAM
end OAI
