import Work.Mixed.Prog.Unload

/-!
# Mixed engine, R3 program part 11: the load and unload passes (agent e-r3prog)

* `CCtx.egear`    the gear of every element (from the per-unit tables; no data needed);
* `Eng.load`      the class table `V0` from the own jobs (live roles through the target placement, scratch 0);
* `Eng.unload`    the own answers `unl` from the class table.
Both are one pass, work `O(NU * 2^n)`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving.RAM
open Binary Matrix Ty Finset GF MXR
noncomputable section

universe U

namespace CCtx
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s} {mx : ℕ}
  {n NU : ι → ℕ} {col ib : ι → ℕ → ℕ}

/-- The gear of every element. -/
lemma egear (h : CCtx m B s x mx n NU col ib) :
    Gear m (fun e : ElemF NU ib => B e.1.1) (fun e => col e.1.1 e.1.2) (p (p s w) w)
      (fun e => ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2)) :=
  ⟨h.elem_unit (ib:=ib) (unit_fetch _ h.tpc), h.elem_unit (ib:=ib) (unit_fetch _ h.tc),
    h.elem_unit (ib:=ib) (unit_fetch _ h.txr), h.uBig.reindex (fun e : ElemF NU ib => e.1)⟩

end CCtx

namespace Eng
variable {α ρ σ : Type*} [Fintype α] [DecidableEq α] [Fintype ρ] [DecidableEq ρ]
  [Fintype σ] [DecidableEq σ]
variable (E : Eng α ρ σ)
variable {ι : Type U} {m : Bool} {B : ι → ℕ} {s : Ty} {x : ι → T s} {cl : Paint}
variable {n hh μ : ι → ℕ}

/-- **Loading the class table from the own jobs.** -/
lemma load
    (h : CCtx m B s x E.X.n n (fun i => E.NUc (n i) (hh i) (μ i)) (fun i p => E.cUc (n i) (hh i) (μ i) p)
      (fun i p => E.iUc (n i) (hh i) (μ i) p))
    {J : ∀ i, ℕ → Bits (n i) → ℂ}
    (tJ : Has m B s (jobsT cl) x (fun i => jobsTape (n i) (E.NUc (n i) (hh i) (μ i) * E.W) (J i))) :
    Can m B s (a (a (a (c cl)))) x
      (fun i => ctab E.X E.R (E.NUc (n i) (hh i) (μ i)) (fun p => E.cUc (n i) (hh i) (μ i) p)
        (fun p => E.iUc (n i) (hh i) (μ i) p) (E.V0 (n i) (hh i) (μ i) (J i)))
      (fun i => E.NUc (n i) (hh i) (μ i) * 2^(n i)) := by
  let NU (i : ι) := E.NUc (n i) (hh i) (μ i)
  let col (i : ι) (p : ℕ) := E.cUc (n i) (hh i) (μ i) p
  let ib (i : ι) (p : ℕ) := E.iUc (n i) (hh i) (μ i) p
  have g := h.egear
  apply h.ctab_build _ rfl
  let yE (e : ElemF NU ib) : T (p (p s w) w) := ((x e.1.1, e.1.2.val), (bits (ib e.1.1 e.1.2)).loc e.2)
  let BE (e : ElemF NU ib) := B e.1.1
  have hsize : Knows m BE (p (p s w) w) yE (fun e => (fleet E.X E.R 0 (col e.1.1 e.1.2)).n) :=
    (Can.wc E.R.n).times (g.grid_size E.X 0) (Small.const _ _) (g.grid_small E.X 0)
  have hlarge : Small BE (fun e => (fleet E.X E.R 0 (col e.1.1 e.1.2)).n) :=
    (Small.const _ E.R.n).mul (g.grid_small E.X 0)
  -- one beach
  let D := Σ e : ElemF NU ib, Beach α ρ 0 (col e.1.1 e.1.2)
  let yD (d : D) : T (p (p (p s w) w) w) := (yE d.1, (fleet E.X E.R 0 (col d.1.1.1 d.1.1.2)).loc d.2)
  let BD (d : D) := BE d.1
  have gD := g.descend (fun d : D => d.1) yD (t:=p (p (p s w) w) w) Can.first
  have hl : Knows m BD (p (p (p s w) w) w) yD (fun d => E.R.loc d.2.1) :=
    Can.addr_high (u:=fun d : D => d.2) (L:=fun _ => E.R) (R:=fun d => book E.X 0 (col d.1.1.1 d.1.1.2))
      Can.second (gD.grid_size E.X 0) (hlarge.reindex (fun d : D => d.1))
  have hk : Knows m BD (p (p (p s w) w) w) yD (fun d => (book E.X 0 (col d.1.1.1 d.1.1.2)).loc d.2.2) :=
    Can.addr_low (u:=fun d : D => d.2) (L:=fun _ => E.R) (R:=fun d => book E.X 0 (col d.1.1.1 d.1.1.2))
      Can.second (gD.grid_size E.X 0) (hlarge.reindex (fun d : D => d.1))
  have hv : Has m BD (p (p (p s w) w) w) (c cl) yD
      (fun d => E.V0 (n d.1.1.1) (hh d.1.1.1) (μ d.1.1.1) (J d.1.1.1) d.1.1.2 d.1.2 d.2) := by
    apply Can.cases_code (q:=fun d : D => d.2.1) E.R hl
    intro r'
    let K := {d : D // d.2.1 = r'}
    let BK (κ : K) := BD κ.1
    cases hL : E.liveOf r' with
    | none =>
      refine (Can.zero (cl:=cl)).cong (fun _ => rfl) (fun κ => ?_)
      show (0 : ℂ) = E.V0 _ _ _ _ κ.1.1.1.2 κ.1.1.2 κ.1.2
      unfold V0
      rw [κ.2, hL]
    | some l0 =>
      let gK := gD.descend (fun κ : K => κ.1) (fun κ => yD κ.1) (t:=p (p (p s w) w) w) Can.id
      obtain ⟨hrest, hfib⟩ := gK.from_gridB E.X 0 (E.St l0) (E.Xt l0) (fun κ => κ.1.2.2)
        (hk.reindex (fun κ : K => κ.1))
      let i_ (κ : K) := κ.1.1.1.1
      let p_ (κ : K) := κ.1.1.1.2
      let c_ (κ : K) := col (i_ κ) (p_ κ)
      have hbig : Small BK (fun κ => 2^(n (i_ κ))) := h.big.reindex (fun κ : K => i_ κ)
      have hk1 : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1) (fun κ => 2^(E.tr l0 * c_ κ)) :=
        gK.powerB (E.tr l0)
      have hk3 : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => 2^(E.tr l0 * c_ κ + (E.Xt l0).n * c_ κ)) :=
        (gK.powerB E.X.n).cong (fun _ => rfl) (fun κ => by simp only [c_, col, i_, p_]; rw [← E.Xt_n l0]; ring_nf)
      have hs1 : Small BK (fun κ => 2^(E.tr l0 * c_ κ)) := gK.boundB _
      have hs3 : Small BK (fun κ => 2^(E.tr l0 * c_ κ + (E.Xt l0).n * c_ κ)) :=
        (gK.boundB E.X.n).mono (fun κ => le_of_eq (by simp only [c_, col, i_, p_]; rw [← E.Xt_n l0]; ring_nf))
      have hsr : Small BK (fun κ => (book (E.Xt l0) 0 (c_ κ)).loc (E.St l0 |>.peel 0 (c_ κ) κ.1.2.2).1) :=
        (gK.grid_small (E.Xt l0) 0).mono (fun κ => ((book _ 0 _).bd _).le)
      have hsf : Small BK (fun κ => (bits (E.tr l0 * c_ κ)).loc (E.St l0 |>.peel 0 (c_ κ) κ.1.2.2).2) :=
        hs1.mono (fun κ => ((bits _).bd _).le)
      have hq : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => (bits (ib (i_ κ) (p_ κ))).loc κ.1.1.2) :=
        ((Can.first.then_do Can.second : Knows m BD (p (p (p s w) w) w) yD
          (fun d => (bits (ib d.1.1.1 d.1.1.2)).loc d.1.2))).reindex (fun κ : K => κ.1)
      have hsq : Small BK (fun κ => (bits (ib (i_ κ) (p_ κ))).loc κ.1.1.2) :=
        hbig.mono (fun κ => le_trans ((bits _).bd _).le (h.two_ib_le _ _))
      have haddr : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1)
          (fun κ => (bits (n (i_ κ))).loc (E.placeT (n (i_ κ)) (hh (i_ κ)) (μ (i_ κ)) (p_ κ) l0
            (κ.1.1.2, κ.1.2.2.2))) := by
        refine ((hfib.plus (hk1.times hrest hs1 hsr) hsf (hs1.mul hsr)).plus (hk3.times hq hs3 hsq)
          (hsf.add (hs1.mul hsr)) (hs3.mul hsq)).cong (fun _ => rfl) (fun κ => ?_)
        unfold placeT
        rw [placeS_loc]
        have hb0 : κ.1.2.2 = ((0 : Bits 0), κ.1.2.2.2) := Prod.ext (Subsingleton.elim _ _) rfl
        rw [hb0]
      -- the job of the slot `p*W + loc l0`
      have hpv : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1) (fun κ => (p_ κ).val) :=
        ((Can.first.then_do (Can.first.then_do Can.second) : Knows m BD (p (p (p s w) w) w) yD
          (fun d => d.1.1.2.val))).reindex (fun κ : K => κ.1)
      have hNUs : Small BK (fun κ => NU (i_ κ)) := h.nuS.reindex (fun κ : K => i_ κ)
      have hpvs : Small BK (fun κ => (p_ κ).val) := hNUs.mono (fun κ => (p_ κ).2.le)
      have hjix : Knows m BK (p (p (p s w) w) w) (fun κ => yD κ.1) (fun κ => (p_ κ).val * E.W + E.Q.loc l0) :=
        (hpv.times (Can.wc E.W) hpvs (Small.const _ _)).plus (Can.wc _) (hpvs.mul (Small.const _ _))
          (Small.const _ _)
      have hjlt : ∀ κ : K, (p_ κ).val * E.W + E.Q.loc l0 < NU (i_ κ) * E.W := fun κ => by
        have h1 := (p_ κ).2
        have h2 : E.Q.loc l0 < E.W := E.Q.bd l0
        calc (p_ κ).val * E.W + E.Q.loc l0 < (p_ κ).val * E.W + E.W := by omega
          _ = ((p_ κ).val + 1) * E.W := by ring
          _ ≤ NU (i_ κ) * E.W := Nat.mul_le_mul_right _ h1
      have tJE : Has m BE (p (p s w) w) (jobsT cl) yE
          (fun e => jobsTape (n e.1.1) (NU e.1.1 * E.W) (J e.1.1)) :=
        (Can.first.then_do Can.first).then_do (tJ.reindex (fun e : ElemF NU ib => e.1.1))
      have tJD : Has m BK (p (p (p s w) w) w) (jobsT cl) (fun κ => yD κ.1)
          (fun κ => jobsTape (n (i_ κ)) (NU (i_ κ) * E.W) (J (i_ κ))) :=
        ((Can.first.then_do (tJE.reindex (fun d : D => d.1)) : Has m BD (p (p (p s w) w) w) (jobsT cl) yD
          (fun d => jobsTape (n d.1.1.1) (NU d.1.1.1 * E.W) (J d.1.1.1)))).reindex (fun κ : K => κ.1)
      have hJt := Can.fetch (t:=a (c cl)) (L:=fun κ : K => Layout.std (NU (i_ κ) * E.W))
        (v:=fun κ t => (bits (n (i_ κ))).tape (J (i_ κ) t))
        (i:=fun κ => (⟨(p_ κ).val * E.W + E.Q.loc l0, hjlt κ⟩ : Fin (NU (i_ κ) * E.W))) tJD hjix
      have hval := Can.fetch (t:=c cl) (L:=fun κ : K => bits (n (i_ κ)))
        (v:=fun κ => J (i_ κ) ((p_ κ).val * E.W + E.Q.loc l0))
        (i:=fun κ => E.placeT (n (i_ κ)) (hh (i_ κ)) (μ (i_ κ)) (p_ κ) l0 (κ.1.1.2, κ.1.2.2.2)) hJt haddr
      refine hval.cong (fun _ => rfl) (fun κ => ?_)
      show _ = E.V0 _ _ _ _ κ.1.1.1.2 κ.1.1.2 κ.1.2
      unfold V0
      rw [κ.2, hL]
  have HL := Can.layout (t:=c cl) (L:=fun e : ElemF NU ib => fleet E.X E.R 0 (col e.1.1 e.1.2))
    (v:=fun e => E.V0 (n e.1.1) (hh e.1.1) (μ e.1.1) (J e.1.1) e.1.2 e.2) (P:=fun _ => 0) (W:=fun _ => 0)
    hsize hv hlarge
  refine HL.mono (fun e => ?_)
  simp [fleet, col]

lemma roleAt_of (t : ℕ) (σ0 : σ) (h : E.Q.loc σ0 = t % E.W) : E.roleAt t = σ0 := by
  apply E.Q.inj; rw [loc_roleAt, h]

/-- **Answering the own slots from the class table.** -/
lemma unload
    (h : CCtx m B s x E.X.n n (fun i => E.NUc (n i) (hh i) (μ i)) (fun i p => E.cUc (n i) (hh i) (μ i) p)
      (fun i p => E.iUc (n i) (hh i) (μ i) p))
    {V : ∀ i, E.CState (n i) (hh i) (μ i)}
    (ct : Has m B s (a (a (a (c cl)))) x (fun i => ctab E.X E.R (E.NUc (n i) (hh i) (μ i))
      (fun p => E.cUc (n i) (hh i) (μ i) p) (fun p => E.iUc (n i) (hh i) (μ i) p) (V i))) :
    Can m B s (jobsT cl) x
      (fun i => jobsTape (n i) (E.NUc (n i) (hh i) (μ i) * E.W) (E.unl (n i) (hh i) (μ i) (V i)))
      (fun i => E.NUc (n i) (hh i) (μ i) * E.W * 2^(n i)) := by
  let NU (i : ι) := E.NUc (n i) (hh i) (μ i)
  let NJ (i : ι) := NU i * E.W
  let col (i : ι) (p : ℕ) := E.cUc (n i) (hh i) (μ i) p
  let ib (i : ι) (p : ℕ) := E.iUc (n i) (hh i) (μ i) p
  have G := h.erig ct
  let EN := EntF NJ n
  let yE (en : EN) : T (p (p s w) w) := ((x en.1.1, en.1.2.val), (bits (n en.1.1)).loc en.2)
  let BE (en : EN) := B en.1.1
  have hE : Knows m BE (p (p s w) w) yE (fun en => (bits (n en.1.1)).loc en.2) := Can.second
  have hbig : Small BE (fun en => 2^(n en.1.1)) := h.big.reindex (fun en : EN => en.1.1)
  have hEs : Small BE (fun en => (bits (n en.1.1)).loc en.2) := hbig.mono (fun en => ((bits _).bd _).le)
  have hNJs : Small B NJ := h.nuS.mul (Small.const _ _)
  have ht : Knows m BE (p (p s w) w) yE (fun en => en.1.2.val) := Can.first.then_do Can.second
  have hts : Small BE (fun en => en.1.2.val) := (hNJs.reindex (fun en : EN => en.1.1)).mono (fun en => en.1.2.2.le)
  have hpt : Knows m BE (p (p s w) w) yE (fun en => en.1.2.val / E.W) := ht.dvdW (Can.wc E.W) hts
  have hat : Knows m BE (p (p s w) w) yE (fun en => en.1.2.val % E.W) := ht.modW (Can.wc E.W) hts
  have hplt : ∀ en : EN, en.1.2.val / E.W < NU en.1.1 := fun en => (Nat.div_lt_iff_lt_mul E.W_pos).2 en.1.2.2
  have hent : Has m BE (p (p s w) w) (c cl) yE
      (fun en => E.unl (n en.1.1) (hh en.1.1) (μ en.1.1) (V en.1.1) en.1.2.val en.2) := by
    apply Can.cases_code (q:=fun en : EN => E.roleAt en.1.2.val) E.Q
      (hat.cong (fun _ => rfl) (fun en => (E.loc_roleAt _).symm))
    intro σ0
    let K := {en : EN // E.roleAt en.1.2.val = σ0}
    let BK (κ : K) := BE κ.1
    let yK (κ : K) := yE κ.1
    let i_ (κ : K) := κ.1.1.1
    let p_ (κ : K) : Fin (NU (i_ κ)) := ⟨κ.1.1.2.val / E.W, hplt κ.1⟩
    let PL (κ : K) := E.placeT (n (i_ κ)) (hh (i_ κ)) (μ (i_ κ)) (p_ κ) σ0
    let qy (κ : K) := (PL κ).symm κ.1.2
    let j (κ : K) : ElemF NU ib := ⟨⟨i_ κ, p_ κ⟩, (qy κ).1⟩
    have hEK := hE.reindex (fun κ : K => κ.1)
    have hEsK := hEs.reindex (fun κ : K => κ.1)
    have hppK : Knows m BK (p (p s w) w) yK (fun κ => (p_ κ).val) := hpt.reindex (fun κ : K => κ.1)
    have hpcK : Knows m BK (p (p s w) w) yK (fun κ => 2^(col (i_ κ) (p_ κ))) :=
      Can.fetch (t:=w) (L:=fun κ : K => Layout.std (NU (i_ κ))) (v:=fun κ q => 2^(col (i_ κ) q)) (i:=p_)
        (((Can.first.then_do Can.first).then_do (h.tpc.reindex (fun en : EN => en.1.1)) : Has m BE (p (p s w) w)
          (a w) yE (fun en => (Layout.std (NU en.1.1)).tape (fun q => 2^(col en.1.1 q)))).reindex
          (fun κ : K => κ.1)) hppK
    have hcbig : Small BK (fun κ => 2^(col (i_ κ) (p_ κ))) :=
      (hbig.reindex (fun κ : K => κ.1)).mono (fun κ => h.two_c_le (i_ κ) (p_ κ))
    have hXc : Knows m BK (p (p s w) w) yK (fun κ => 2^(E.tr σ0 * col (i_ κ) (p_ κ) +
        (E.Xt σ0).n * col (i_ κ) (p_ κ))) :=
      (hpcK.wpow_fixed E.X.n hcbig).cong (fun _ => rfl) (fun κ => by
        rw [← pow_mul, ← E.Xt_n σ0]; ring_nf)
    have hq : Knows m BK (p (p s w) w) yK (fun κ => (bits (ib (i_ κ) (p_ κ))).loc (qy κ).1) :=
      (hEK.dvdW hXc hEsK).cong (fun _ => rfl) (fun κ => (placeS_parts (E.St σ0) (E.Xt σ0) _ _ _ _ κ.1.2).2.2)
    have hxK : Has m BK (p (p s w) w) s yK (fun κ => x (i_ κ)) :=
      ((Can.first.then_do Can.first : Has m BE (p (p s w) w) s yE (fun en => x en.1.1))).reindex (fun κ : K => κ.1)
    have g : Has m (fun κ => B (j κ).1.1) (p (p s w) w) (p (p s w) w) yK
        (fun κ => ((x (j κ).1.1, (j κ).1.2.val), (bits (ib (j κ).1.1 (j κ).1.2)).loc (j κ).2)) :=
      (hxK.pair hppK).pair hq
    have G1 := G.gather j yK g
    have hk1 : Knows m BK (p (p s w) w) yK (fun κ => 2^(E.tr σ0 * col (i_ κ) (p_ κ))) := G1.rig.gear.powerB _
    have hk2 : Knows m BK (p (p s w) w) yK (fun κ => 2^((E.Xt σ0).n * col (i_ κ) (p_ κ))) := G1.rig.gear.powerB _
    have hfib : Knows m BK (p (p s w) w) yK
        (fun κ => (bits (E.tr σ0 * col (i_ κ) (p_ κ))).loc (E.St σ0 |>.peel 0 (col (i_ κ) (p_ κ)) (0, (qy κ).2)).2) :=
      (hEK.modW hk1 hEsK).cong (fun _ => rfl) (fun κ => (placeS_parts (E.St σ0) (E.Xt σ0) _ _ _ _ κ.1.2).1)
    have hrest : Knows m BK (p (p s w) w) yK
        (fun κ => (book (E.Xt σ0) 0 (col (i_ κ) (p_ κ))).loc (E.St σ0 |>.peel 0 (col (i_ κ) (p_ κ)) (0, (qy κ).2)).1) :=
      ((hEK.dvdW hk1 hEsK).modW hk2 hEsK.div).cong (fun _ => rfl)
        (fun κ => (placeS_parts (E.St σ0) (E.Xt σ0) _ _ _ _ κ.1.2).2.1)
    have hgrid := G1.rig.gear.to_gridB E.X 0 (E.St σ0) (E.Xt σ0)
      (fun κ => (E.St σ0).peel 0 (col (i_ κ) (p_ κ)) (0, (qy κ).2)) hrest hfib
    have hb : Knows m BK (p (p s w) w) yK (fun κ => (book E.X 0 (col (i_ κ) (p_ κ))).loc ((0 : Bits 0), (qy κ).2)) :=
      hgrid.cong (fun _ => rfl) (fun κ => by rw [Equiv.symm_apply_apply])
    have hr := G1.rig.read (fun κ => (E.e σ0, ((0 : Bits 0), (qy κ).2))) (Can.wc (E.R.loc (E.e σ0))) hb
    refine hr.cong (fun _ => rfl) (fun κ => ?_)
    show V (i_ κ) (p_ κ) (qy κ).1 (E.e σ0, (0, (qy κ).2)) = E.unl _ _ _ _ κ.1.1.2.val κ.1.2
    unfold unl
    rw [dif_pos (hplt κ.1)]
    have hr0 : E.roleAt κ.1.1.2.val = σ0 := κ.2
    simp only [hr0]
    rfl
  have hlen : Knows m (fun u : UnitF NJ => B u.1) (p s w) (fun u => (x u.1, u.2.val)) (fun u => (bits (n u.1)).n) :=
    Can.first.then_do (h.p2n.reindex (fun u : UnitF NJ => u.1))
  have HI := Can.layout (t:=c cl) (L:=fun u : UnitF NJ => bits (n u.1))
    (v:=fun u e => E.unl (n u.1) (hh u.1) (μ u.1) (V u.1) u.2.val e) (P:=fun _ => 0) (W:=fun _ => 0) hlen
    (hent.mono (fun _ => Nat.zero_le _)) (h.big.reindex (fun u : UnitF NJ => u.1))
  have HO := Can.layout (t:=a (c cl)) (L:=fun i => Layout.std (NJ i))
    (v:=fun i t => (bits (n i)).tape (E.unl (n i) (hh i) (μ i) (V i) t)) (P:=fun i => 2^(n i))
    (W:=fun _ => 0) (h.nuK.times (Can.wc E.W) h.nuS (Small.const _ _)) (HI.mono (fun u => by simp)) hNJs
  refine HO.weaken ⟨3, fun i => ?_⟩
  change 0 + NJ i + NJ i * 2^(n i) + 1 ≤ 3 * (NJ i * 2^(n i) + 1)
  have h2 : NJ i ≤ NJ i * 2^(n i) := Nat.le_mul_of_pos_right _ (Nat.two_pow_pos _)
  omega

end Eng

end
end PowerSaving.RAM
end OAI
