import Work.Gauge.Inv
import Work.XMut.Scal

/-!
# x-mutable (13 Oct 2026): THE GAUGED INVOCATION THEOREM WITH MUTABLE SOURCES, `GLabZ.inv : GInvE`

The gauged invocation theorem of Nil Mamano's extension (`Work/Gauge/Inv.lean`, `GLab.inv`, on Jacob Sussman's
`Lab.inv`, with Dugongue's #401 entrance gauges), for a scalar half WITHOUT `hx`: the main phase may write the
source roles (helper-to-source and source-to-source additions at a common label, before and after the scatter).
Only change to the word: the end gate is `g8z` (x and slot rows of the inverse; legal because sources and slots
all stand at the old frame of the full label there) instead of `g8` (slot rows only).  The conclusion is the SAME
package `GX.GInvE` (same `cost`, same `ent`), so everything downstream (`GX.engine_of_ginv_401`, the bank
schedule, the rate, the Fourier seam) applies unchanged.

`GLabZ` is `GLab` (`Work/Gauge/Defs.lean`) over `GScalZ` (= `GScal` without `hx`), field for field.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace GX
open Binary Matrix Finset RAM SS CB BR CR GF RF
noncomputable section

/-- **Scalar half with entrance gauges, sources mutable**: a gauged slot never reaches a target (E7). -/
structure GScalZ (T Sl C : Type) [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] extends ScalZ T Sl C where
  gauged : Finset Sl
  hE : ∀ q ∈ gauged, ∀ t, (MB * scatM (Jr * Cc) * MA) (lY t) (lS q) = 0

/-- **Label half with entrance labels over a mutable-source scalar half** (`GLab`, field for field). -/
structure GLabZ (H : Type) [Fintype H] [DecidableEq H] {T Sl C : Type} [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] (S : GScalZ T Sl C) where
  tv : T → Space H
  hunit : ∀ t, dot (tv t) (tv t) = 1
  hgap : ∀ t, ∃ p, tv t p = 0
  ent : Sl → SLbl H
  lab0 : L3 T Sl → SLbl H
  labA : L3 T Sl → SLbl H
  labB : L3 T Sl → SLbl H
  labF : L3 T Sl → SLbl H
  cen : C → SLbl H
  z0 : SLbl H
  hz0 : z0.dim = 0
  h0x : ∀ t, (lab0 (lX t)).dim = 1 ∧ (lab0 (lX t)).Mem (tv t)
  h0y : ∀ t, (lab0 (lY t)).dim = 0
  h0s' : ∀ q, lab0 (lS q) = ent q
  hent : ∀ q, q ∉ S.gauged → (ent q).dim = 0
  hAy : ∀ t, labA (lY t) = z0
  hAc : ∀ k q, S.Cc k q ≠ 0 → cen k = labA (lS q)
  hcen : ∀ k, 0 < (cen k).dim
  hFx : ∀ t, (labF (lX t)).dim = Fintype.card H
  hFs : ∀ q, (labF (lS q)).dim = Fintype.card H
  hFy : ∀ t, (labF (lY t)).dim + 1 = Fintype.card H ∧
    ∀ x, (labF (lY t)).Mem x → dot (tv t) x = 0
  rA : List ℕ
  rB : List ℕ
  rF : List ℕ
  pA : ∀ {α : Type} [Fintype α] [DecidableEq α] (X : BXMap H α),
    Fintype.card H < Fintype.card α →
    GRoute (fun r => X.Φ (lab0 r)) (fun r => X.Φ (labA r)) (actPoint S.MA)
        (fun φ => rcost φ rA) ∧
    GRoute (fun r => X.Φ (lab0 r)) (fun r => X.Φ (labA r)) (actPoint S.MAiᵀ)
        (fun φ => rcost φ rA)
  pB : ∀ {α : Type} [Fintype α] [DecidableEq α] (X : BXMap H α),
    Fintype.card H < Fintype.card α →
    GRoute (fun r => X.Φ (labA r)) (fun r => X.Φ (labB r)) (actPoint S.MB)
        (fun φ => rcost φ rB) ∧
    GRoute (fun r => X.Φ (labA r)) (fun r => X.Φ (labB r)) (actPoint S.MBiᵀ)
        (fun φ => rcost φ rB)
  pF : ∀ {α : Type} [Fintype α] [DecidableEq α] (X : BXMap H α),
    Fintype.card H < Fintype.card α →
    GRoute (fun r => X.Φ (labB r)) (fun r => X.Φ (labF r)) id (fun φ => rcost φ rF)

/-- price of one invocation (as `GLab.cost`). -/
def GLabZ.cost {H : Type} [Fintype H] [DecidableEq H] {T Sl C : Type} [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] {S : GScalZ T Sl C} (L : GLabZ H S)
    (φ : ℕ → ℝ) : ℝ :=
  rcost φ L.rA + (∑ k, bcost φ (L.cen k).dim) + rcost φ L.rB + rcost φ L.rF

/-- at the END of an invocation the sources and the slots share the old frame of the full label. -/
lemma endZ {H T Sl α : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T] [Fintype Sl]
    [DecidableEq Sl] [Fintype α] [DecidableEq α] (Xm : XMap H α) (tv : T → Space H) (i j : L3 T Sl)
    (hi : iY i = 0) (hj : iY j = 0) : Xm.Φ (labE tv i) = Xm.Φ (labE tv j) := by
  rcases i with ((t|t)|q)
  · rcases j with ((t'|t')|q')
    · rfl
    · exact absurd (show (1:ℚ) = 0 from hj) one_ne_zero
    · rfl
  · exact absurd (show (1:ℚ) = 0 from hi) one_ne_zero
  · rcases j with ((t'|t')|q')
    · rfl
    · exact absurd (show (1:ℚ) = 0 from hj) one_ne_zero
    · rfl

section Thm
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
variable {α : Type} [Fintype α] [DecidableEq α]
variable {S : GScalZ T Sl C} (L : GLabZ H S) (W : GWin H α)

theorem GLabZ.pre : GRoute (gst0 W L.tv L.ent) (gst1 W L.tv L.ent S.gauged) id (fun _ => 0) := by
  refine (GRoute.reach_all _ _ (fun _ _ => 0) ?_).cast (fun φ => by simp)
  rintro ((t|t)|q)
  · exact GReach.refl _
  · exact GReach.refl _
  · show GReach (W.X.Φ (L.ent q)) (if q ∈ S.gauged then W.X.Φ (L.ent q) else W.Xm.Φ 0) _
    by_cases hq : q ∈ S.gauged
    · rw [if_pos hq]; exact GReach.refl _
    · rw [if_neg hq]; exact W.zeroB _ (L.hent q hq)

theorem GLabZ.ad0 :
    GRoute (gst1 W L.tv L.ent S.gauged) (fun r => W.X.Φ (L.lab0 r)) id (fun _ => 0) := by
  refine (GRoute.reach_all _ _ (fun _ _ => 0) ?_).cast (fun φ => by simp)
  rintro ((t|t)|q)
  · obtain ⟨p, hp⟩ := L.hgap t
    obtain ⟨A, k, hA⟩ := exists_base_one (L.tv t) (L.hunit t) p hp
    have h := W.line A k (L.lab0 (lX t)) (L.h0x t).1 (by rw [hA]; exact (L.h0x t).2)
    rw [hA] at h
    exact h
  · exact W.zero _ (L.h0y t)
  · show GReach (if q ∈ S.gauged then W.X.Φ (L.ent q) else W.Xm.Φ 0) (W.X.Φ (L.lab0 (lS q))) _
    rw [L.h0s' q]
    by_cases hq : q ∈ S.gauged
    · rw [if_pos hq]; exact GReach.refl _
    · rw [if_neg hq]; exact W.zero _ (L.hent q hq)

theorem GLabZ.adF :
    GRoute (fun r => W.X.Φ (L.labF r)) (fun r : L3 T Sl => W.Xm.Φ (labE L.tv r)) id
      (fun _ => 0) := by
  refine (GRoute.reach_all _ _ (fun _ _ => 0) ?_).cast (fun φ => by simp)
  rintro ((t|t)|q)
  · exact W.full _ (L.hFx t)
  · obtain ⟨p, hp⟩ := L.hgap t
    obtain ⟨A, k, hA⟩ := exists_base_one (L.tv t) (L.hunit t) p hp
    have h := W.perp A k (L.labF (lY t)) (L.hFy t).1 (by rw [hA]; exact (L.hFy t).2)
    rw [hA] at h
    exact h
  · exact W.full _ (L.hFs q)

theorem GLabZ.post :
    GRoute (fun r : L3 T Sl => W.Xm.Φ (labE L.tv r)) (gen2 W L.tv L.ent) id (fun _ => 0) := by
  refine (GRoute.reach_all _ _ (fun _ _ => 0) ?_).cast (fun φ => by simp)
  rintro ((t|t)|q)
  · exact GReach.refl _
  · exact GReach.refl _
  · exact W.topN _ (SLbl.top_dim _)

lemma GLabZ.cfacts (hb : Fintype.card H < Fintype.card α) (k : C) :
    (∀ x, L.z0.Mem x → (L.cen k).Mem x) ∧ L.z0.dim < (L.cen k).dim ∧
    (L.cen k).dim - L.z0.dim < Fintype.card α := by
  refine ⟨fun x hx => ?_, by rw [L.hz0]; exact L.hcen k, ?_⟩
  · rw [mem_dim_zero L.z0 L.hz0 hx]
    exact insub_zero _ _
  · have h1 : (L.cen k).dim ≤ Fintype.card H := Finset.card_le_univ _
    omega

lemma GLabZ.cdn (k : C) :
    GReach (W.X.Φ (L.cen k)) (W.X.Φ L.z0) (fun φ => bcost φ (L.cen k).dim) := by
  obtain ⟨h0, h1, h2⟩ := L.cfacts W.big k
  have h := (W.X.fmv_sup L.z0 (L.cen k) h0 h1 h2).reach
  have h' : GReach (W.X.Φ (L.cen k)) (W.X.Φ L.z0)
      (fun φ => rcost φ [(L.cen k).dim - L.z0.dim]) := GReach.ofReach h
  exact h'.cast (fun φ => by rw [rcost_single, L.hz0, Nat.sub_zero])

lemma GLabZ.cupk (k : C) :
    GReach (W.X.Φ L.z0) (W.X.Φ (L.cen k)) (fun φ => bcost φ (L.cen k).dim) := by
  obtain ⟨h0, h1, h2⟩ := L.cfacts W.big k
  have h := (W.X.fmv_sub L.z0 (L.cen k) h0 h1 h2).reach
  have h' : GReach (W.X.Φ L.z0) (W.X.Φ (L.cen k))
      (fun φ => rcost φ [(L.cen k).dim - L.z0.dim]) := GReach.ofReach h
  exact h'.cast (fun φ => by rw [rcost_single, L.hz0, Nat.sub_zero])

theorem GLabZ.cdown (Sf : L3 T Sl → CMat α) :
    GRoute (glue Sf (fun k => W.X.Φ (L.cen k))) (glue Sf (fun _ : C => W.X.Φ L.z0)) id
      (fun φ => ∑ k, bcost φ (L.cen k).dim) :=
  greach_copies Sf _ _ _ (fun k => L.cdn W k)

theorem GLabZ.cup (Sf : L3 T Sl → CMat α) :
    GRoute (glue Sf (fun _ : C => W.X.Φ L.z0)) (glue Sf (fun k => W.X.Φ (L.cen k))) id
      (fun φ => ∑ k, bcost φ (L.cen k).dim) :=
  greach_copies Sf _ _ _ (fun k => L.cupk W k)

/-- the gate `y -= K s` (or its transposed inverse) reads only ungauged slots. -/
lemma GLabZ.gate_ok {i j : L3 T Sl} (t : T) (q : Sl) (hi : i = lY t) (hj : j = lS q)
    (hM : S.toScalZ.Mt (lY t) (lS q) ≠ 0) :
    gst1 W L.tv L.ent S.gauged i = gst1 W L.tv L.ent S.gauged j := by
  subst hi; subst hj
  have hq : q ∉ S.gauged := fun hq => hM (S.hE q hq t)
  show W.Xm.Φ 0 = if q ∈ S.gauged then W.X.Φ (L.ent q) else W.Xm.Φ 0
  rw [if_neg hq]

/-- **Forward word on clean frames, gauged, sources mutable.** -/
theorem GLabZ.fwd_erasedE (c0 c1 : C → CMat α) :
    GRoute (iin W.Xm L.tv (fun q => W.X.Φ (L.ent q)) c0)
      (iout W.Xm L.tv (fun q => W.X.Φ (L.ent q).top) c1) shearF L.cost := by
  have a1 : GRoute (gst1 W L.tv L.ent S.gauged) (gst1 W L.tv L.ent S.gauged)
      (actPoint (g1 S.toScalZ.Mt)) (fun _ => 0) :=
    GRoute.gate (g1 S.toScalZ.Mt) _ _ (fun i j hij => by
      by_cases e : i = j
      · subst e; rfl
      · obtain ⟨t, q, hi, hj, hM⟩ := g1_off hij e
        exact L.gate_ok W t q hi hj hM)
  have b3 : GRoute (fun r : L3 T Sl => W.Xm.Φ (labE L.tv r))
      (fun r : L3 T Sl => W.Xm.Φ (labE L.tv r))
      (actPoint (XMut.g8z S.toScalZ.Mti)) (fun _ => 0) :=
    GRoute.gate (XMut.g8z S.toScalZ.Mti) _ _ (fun i j hij => by
      by_cases e : i = j
      · subst e; rfl
      · have f := XMut.g8z_facts S.toScalZ.hiy
        obtain ⟨h1, h2⟩ := XMut.offdiag_z f.1 f.2 hij e
        exact endZ W.Xm L.tv i j h1 h2)
  have seg1 : GRoute (gst0 W L.tv L.ent) (fun r => W.X.Φ (L.labA r))
      S.toScalZ.f1 (fun φ => rcost φ L.rA) :=
    ((((L.pre W).trans a1).trans (L.ad0 W)).trans (L.pA W.X W.big).1).cast (fun φ => by simp)
  have seg2 : GRoute (fun r => W.X.Φ (L.labA r)) (gen2 W L.tv L.ent)
      S.toScalZ.f2 (fun φ => rcost φ L.rB + rcost φ L.rF) :=
    (((((L.pB W.X W.big).1.trans (L.pF W.X W.big)).trans (L.adF W)).trans b3).trans (L.post W)).cast
      (fun φ => by simp)
  have e0 : GRoute (iin W.Xm L.tv (fun q => W.X.Φ (L.ent q)) c0)
      (glue (gst0 W L.tv L.ent) (fun k => W.X.Φ (L.cen k))) eraseC (fun _ => 0) :=
    GRoute.ofX (xerase (fun t => rfl) (fun t => rfl) (fun q => rfl))
  have L1 := glift_live (fun k => W.X.Φ (L.cen k)) seg1
  have m1 : GRoute (glue (fun r => W.X.Φ (L.labA r)) (fun k => W.X.Φ (L.cen k)))
      (glue (fun r => W.X.Φ (L.labA r)) (fun k => W.X.Φ (L.cen k)))
      (fun z => addBlk bC (ap S.Cc (z ∘ bS)) z) (fun _ => 0) :=
    GRoute.ofX (xgate_shear bC bS bS_inj S.Cc (fun k q h => congrArg W.X.Φ (L.hAc k q h)))
  have m2 := L.cdown W (fun r => W.X.Φ (L.labA r))
  have m3 : GRoute (glue (fun r => W.X.Φ (L.labA r)) (fun _ : C => W.X.Φ L.z0))
      (glue (fun r => W.X.Φ (L.labA r)) (fun _ : C => W.X.Φ L.z0))
      (fun z => addBlk bY (ap S.Jr (z ∘ bC)) z) (fun _ => 0) :=
    GRoute.ofX (xgate_shear bY bC bC_inj S.Jr (fun t k _ => congrArg W.X.Φ (L.hAy t)))
  have L2 := glift_live (fun _ : C => W.X.Φ L.z0) seg2
  have e1 : GRoute (glue (gen2 W L.tv L.ent) (fun _ : C => W.X.Φ L.z0))
      (iout W.Xm L.tv (fun q => W.X.Φ (L.ent q).top) c1) eraseC (fun _ => 0) :=
    GRoute.ofX (xerase (fun t => rfl) (fun t => rfl) (fun q => rfl))
  have tot := (((((e0.trans L1).trans m1).trans m2).trans m3).trans L2).trans e1
  refine (tot.castg (g' := shearF) (funext fun z => S.toScalZ.fwd_scalar z)).cast (fun φ => ?_)
  simp only [GLabZ.cost, zero_add, add_zero]
  ring

/-- **Backward word on clean frames, gauged, sources mutable.** -/
theorem GLabZ.bwd_erasedE (c0 c1 : C → CMat α) :
    GRoute (iin W.Xm L.tv (fun q => W.X.Φ (L.ent q)) c0)
      (iout W.Xm L.tv (fun q => W.X.Φ (L.ent q).top) c1) shearB L.cost := by
  have a1 : GRoute (gst1 W L.tv L.ent S.gauged) (gst1 W L.tv L.ent S.gauged)
      (actPoint (g1i S.toScalZ.Mt)ᵀ) (fun _ => 0) :=
    GRoute.gate (g1i S.toScalZ.Mt)ᵀ _ _ (fun i j hij => by
      by_cases e : i = j
      · subst e; rfl
      · have hij' : g1i S.toScalZ.Mt j i ≠ 0 := hij
        obtain ⟨t, q, hj, hi, hM⟩ := g1i_off hij' (Ne.symm e)
        exact (L.gate_ok W t q hj hi hM).symm)
  have b3 : GRoute (fun r : L3 T Sl => W.Xm.Φ (labE L.tv r))
      (fun r : L3 T Sl => W.Xm.Φ (labE L.tv r))
      (actPoint (XMut.h8z S.toScalZ.Mt)ᵀ) (fun _ => 0) :=
    GRoute.gate (XMut.h8z S.toScalZ.Mt)ᵀ _ _ (fun i j hij => by
      by_cases e : i = j
      · subst e; rfl
      · have hij' : XMut.h8z S.toScalZ.Mt j i ≠ 0 := hij
        have f := XMut.h8z_facts S.toScalZ.hmy
        obtain ⟨h1, h2⟩ := XMut.offdiag_z f.1 f.2 hij' (Ne.symm e)
        exact (endZ W.Xm L.tv j i h1 h2).symm)
  have seg1 : GRoute (gst0 W L.tv L.ent) (fun r => W.X.Φ (L.labA r))
      S.toScalZ.b1 (fun φ => rcost φ L.rA) :=
    ((((L.pre W).trans a1).trans (L.ad0 W)).trans (L.pA W.X W.big).2).cast (fun φ => by simp)
  have seg2 : GRoute (fun r => W.X.Φ (L.labA r)) (gen2 W L.tv L.ent)
      S.toScalZ.b2 (fun φ => rcost φ L.rB + rcost φ L.rF) :=
    (((((L.pB W.X W.big).2.trans (L.pF W.X W.big)).trans (L.adF W)).trans b3).trans (L.post W)).cast
      (fun φ => by simp)
  have e0 : GRoute (iin W.Xm L.tv (fun q => W.X.Φ (L.ent q)) c0)
      (glue (gst0 W L.tv L.ent) (fun _ : C => W.X.Φ L.z0)) eraseC (fun _ => 0) :=
    GRoute.ofX (xerase (fun t => rfl) (fun t => rfl) (fun q => rfl))
  have L1 := glift_live (fun _ : C => W.X.Φ L.z0) seg1
  have m1 : GRoute (glue (fun r => W.X.Φ (L.labA r)) (fun _ : C => W.X.Φ L.z0))
      (glue (fun r => W.X.Φ (L.labA r)) (fun _ : C => W.X.Φ L.z0))
      (fun z => addBlk bC (ap S.Jrᵀ (z ∘ bY)) z) (fun _ => 0) :=
    GRoute.ofX (xgate_shear bC bY bY_inj S.Jrᵀ (fun k t _ => (congrArg W.X.Φ (L.hAy t)).symm))
  have m2 := L.cup W (fun r => W.X.Φ (L.labA r))
  have m3 : GRoute (glue (fun r => W.X.Φ (L.labA r)) (fun k => W.X.Φ (L.cen k)))
      (glue (fun r => W.X.Φ (L.labA r)) (fun k => W.X.Φ (L.cen k)))
      (fun z => addBlk bS (ap (-S.Ccᵀ) (z ∘ bC)) z) (fun _ => 0) :=
    GRoute.ofX (xgate_shear bS bC bC_inj (-S.Ccᵀ)
      (fun q k h => (congrArg W.X.Φ (L.hAc k q (ne_of_negT S.Cc h))).symm))
  have L2 := glift_live (fun k => W.X.Φ (L.cen k)) seg2
  have e1 : GRoute (glue (gen2 W L.tv L.ent) (fun k => W.X.Φ (L.cen k)))
      (iout W.Xm L.tv (fun q => W.X.Φ (L.ent q).top) c1) eraseC (fun _ => 0) :=
    GRoute.ofX (xerase (fun t => rfl) (fun t => rfl) (fun q => rfl))
  have tot := (((((e0.trans L1).trans m1).trans m2).trans m3).trans L2).trans e1
  refine (tot.castg (g' := shearB) (funext fun z => S.toScalZ.bwd_scalar z)).cast (fun φ => ?_)
  simp only [GLabZ.cost, zero_add, add_zero]
  ring

/-- **THE GAUGED INVOCATION THEOREM, SOURCES MUTABLE.**  Same package `GInvE` as `GLab.inv`. -/
def GLabZ.inv : GInvE H T Sl C where
  tv := L.tv
  hunit := L.hunit
  hgap := L.hgap
  cost := L.cost
  ent := L.ent
  fwd := fun W N Wn hN hW Z c0 c1 =>
    gdirtyE W.Xm L.tv _ _ N Wn hN hW Z c0 c1 (L.fwd_erasedE W c0 c1)
      (shearF_comm (fun _ => 1) (fun q => N q * Z q) (fun _ => 1))
  bwd := fun W N Wn hN hW Z c0 c1 =>
    gdirtyE W.Xm L.tv _ _ N Wn hN hW Z c0 c1 (L.bwd_erasedE W c0 c1)
      (shearB_comm (fun _ => 1) (fun q => N q * Z q) (fun _ => 1))

@[simp] lemma GLabZ.inv_cost (φ : ℕ → ℝ) : L.inv.cost φ = L.cost φ := rfl
@[simp] lemma GLabZ.inv_tv : L.inv.tv = L.tv := rfl
@[simp] lemma GLabZ.inv_ent : L.inv.ent = L.ent := rfl

end Thm
end
end GX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.GX.GLabZ.fwd_erasedE
#print axioms OAI.PowerSaving.GX.GLabZ.bwd_erasedE
#print axioms OAI.PowerSaving.GX.GLabZ.inv
