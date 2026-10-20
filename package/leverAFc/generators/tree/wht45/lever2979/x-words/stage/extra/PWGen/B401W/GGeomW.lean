import Work.Bridge401.GGeom
import Work.PWGen.B401W.GNetCertW
/-! (key: p-wgen) WINDOWED copy, made by p-wgen's tools/mk_b401w.py from Work/Bridge401 (round wht45/lever, 19 Oct
2026): the package is p-warpchk's `KLW.GInvEW H T Sl C Ok` (routes only at windows with `Ok`), every declared name
carries a W, and the window facts `Ok` are threaded to the two places that call the package's routes.  The original
text follows. -/


set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace GX
open Binary Matrix Finset RAM SS CB RF BR GF
noncomputable section
open BC
section Frames401W
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
variable {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}

/-- the bank frames of the schedule `S` (`sched_frames`): before the free adapter of each sub-stage. -/
def zb401W (I : KLW.GInvEW H T Sl C Ok) (l : H) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : ∀ b, S.load (fun q => Fintype.card H - (I.ent q).dim) b + S.pad b = 73)
    (hpad : ∀ b, S.pad b < 73) : ℕ → Orth (BC.LC (Fin 8) H l) → Fin S.nb → CMat (BC.LC (Fin 8) H l) :=
  Classical.choose (sched_frames S I.ent (win401 l hH) (by rw [card_LC401 l hH]; exact hres)
    (by rw [card_LC401 l hH]; exact hpad))

/-- the bank frames of the schedule `S`: after the free adapter of each sub-stage. -/
def zp401W (I : KLW.GInvEW H T Sl C Ok) (l : H) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : ∀ b, S.load (fun q => Fintype.card H - (I.ent q).dim) b + S.pad b = 73)
    (hpad : ∀ b, S.pad b < 73) : ℕ → Orth (BC.LC (Fin 8) H l) → Fin S.nb → CMat (BC.LC (Fin 8) H l) :=
  Classical.choose (Classical.choose_spec (sched_frames S I.ent (win401 l hH) (by rw [card_LC401 l hH]; exact hres)
    (by rw [card_LC401 l hH]; exact hpad)))

lemma zfr401W_spec (I : KLW.GInvEW H T Sl C Ok) (l : H) (hH : Fintype.card H = 9) (S : BankSched Sl)
    (hres : ∀ b, S.load (fun q => Fintype.card H - (I.ent q).dim) b + S.pad b = 73)
    (hpad : ∀ b, S.pad b < 73) :
    (∀ g b, zb401W I l hH S hres hpad 0 g b = 1) ∧
      (∀ (ρ : Fin 8) (k : Fin 9) g b, GReach (zb401W I l hH S hres hpad (9 * ρ.val + k.val) g b)
        (zp401W I l hH S hres hpad (9 * ρ.val + k.val) g b) (fun _ => 0)) ∧
      (∀ (ρ : Fin 8) (k : Fin 9) g q, entW (win401 l hH k g).X (I.ent q)
          * zp401W I l hH S hres hpad (9 * ρ.val + k.val) g (S.use ρ k q)
        = zb401W I l hH S hres hpad (9 * ρ.val + k.val + 1) g (S.use ρ k q)) ∧
      (∀ (ρ : Fin 8) (k : Fin 9) g b, (∀ q, S.use ρ k q ≠ b) →
        zb401W I l hH S hres hpad (9 * ρ.val + k.val + 1) g b = zp401W I l hH S hres hpad (9 * ρ.val + k.val) g b) ∧
      (∀ g b, GReach (zb401W I l hH S hres hpad 72 g b) (kernel (BC.LC (Fin 8) H l)) (fun φ => bcost φ (S.pad b))) :=
  Classical.choose_spec (Classical.choose_spec (sched_frames S I.ent (win401 l hH)
    (by rw [card_LC401 l hH]; exact hres) (by rw [card_LC401 l hH]; exact hpad)))

end Frames401W

section Net
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
variable {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}
variable (I : KLW.GInvEW H T Sl C Ok) (base : T → OBase H) (l : H)
  (hbase : ∀ t, I.tv t = (base t).v l) (hH : Fintype.card H = 9)
include hbase

lemma g401Wgeom_x01 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    PhiC g ∅ 1 = PhiC (RrC (base t) l (0 : Fin 8) g) (qC bs401_1) (tt (I.tv t)) := by
  rw [hbase t]
  exact bankT_x0 (base t) l (0 : Fin 8) g

lemma g401Wgeom_y01 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    PhiC g ∅ (1 + tt (I.tv t)) = PhiC (RrC (base t) l (0 : Fin 8) g) (qC bs401_1) 0 := by
  rw [hbase t]
  exact bankT_y0 (base t) l (0 : Fin 8) g

lemma g401Wgeom_xx (bs : Finset (Fin 8)) (c c' : Fin 8) (hc : c ∈ bs) (hc' : c' ∉ bs)
    (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    PhiC (RrC (base t) l c g) (qC bs) 1
      = PhiC (RrC (base t) l c' g) (qC (insert c' bs)) (tt (I.tv t)) := by
  rw [hbase t]
  exact bankT_xx (base t) l bs c c' hc hc' g

lemma g401Wgeom_yy (bs : Finset (Fin 8)) (c c' : Fin 8) (hc : c ∈ bs) (hc' : c' ∉ bs)
    (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    PhiC (RrC (base t) l c g) (qC bs) (1 + tt (I.tv t))
      = PhiC (RrC (base t) l c' g) (qC (insert c' bs)) 0 := by
  rw [hbase t]
  exact bankT_yy (base t) l bs c c' hc hc' g

lemma g401Wgeom_term (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    kernel (BC.LC (Fin 8) H l) * PhiC g ∅ (tt (I.tv t))
      = shift ((BBC (base t) l g).v (Sum.inl l)) * frame (BBC (base t) l g) (univ \ {Sum.inl l}) := by
  rw [hbase t]
  exact class_term (base t) l g

include hH

/-- pair 2 before its first gate layer, `X`: ONE block of rank 16. -/
lemma g401Wgeom_aX2 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC g ∅ (tt (I.tv t))) (PhiC (RrC (base t) l (0 : Fin 8) g) (qC bs401_1) 1)
      (fun φ => bcost φ 16) := by
  rw [hbase t, classT_X0 (base t) l g, classT_X (base t) l bs401_1 (0 : Fin 8) (Finset.mem_singleton_self _) g]
  exact climb_X0 hB401 (by omega) l (BBC (base t) l g) bs401_1 _ (by rw [card_bs401_1]; omega)

/-- pair 2 before its first gate layer, `Y`: ONE block of rank 16. -/
lemma g401Wgeom_aY2 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC g ∅ 0) (PhiC (RrC (base t) l (0 : Fin 8) g) (qC bs401_1) (1 + tt (I.tv t)))
      (fun φ => bcost φ 16) := by
  rw [hbase t, classT_Y0 (base t) l g, classT_Y (base t) l bs401_1 (0 : Fin 8) (Finset.mem_singleton_self _) g]
  exact climb_Y0 hB401 (by omega) l (BBC (base t) l g) bs401_1 _ (by rw [card_bs401_1]; omega)

/-- pair 3 before its first gate layer, `X`: ONE block of rank 32. -/
lemma g401Wgeom_aX3 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC g ∅ (tt (I.tv t))) (PhiC (RrC (base t) l (2 : Fin 8) g) (qC bs401_3) 1)
      (fun φ => bcost φ 32) := by
  rw [hbase t, classT_X0 (base t) l g, classT_X (base t) l bs401_3 (2 : Fin 8) (Finset.mem_insert_self _ _) g]
  exact climb_X0 hB401 (by omega) l (BBC (base t) l g) bs401_3 _ (by rw [card_bs401_3]; omega)

/-- pair 3 before its first gate layer, `Y`: ONE block of rank 32. -/
lemma g401Wgeom_aY3 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC g ∅ 0) (PhiC (RrC (base t) l (2 : Fin 8) g) (qC bs401_3) (1 + tt (I.tv t)))
      (fun φ => bcost φ 32) := by
  rw [hbase t, classT_Y0 (base t) l g, classT_Y (base t) l bs401_3 (2 : Fin 8) (Finset.mem_insert_self _ _) g]
  exact climb_Y0 hB401 (by omega) l (BBC (base t) l g) bs401_3 _ (by rw [card_bs401_3]; omega)

/-- pair 4 before its first gate layer, `X`: ONE block of rank 48. -/
lemma g401Wgeom_aX4 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC g ∅ (tt (I.tv t))) (PhiC (RrC (base t) l (4 : Fin 8) g) (qC bs401_5) 1)
      (fun φ => bcost φ 48) := by
  rw [hbase t, classT_X0 (base t) l g, classT_X (base t) l bs401_5 (4 : Fin 8) (Finset.mem_insert_self _ _) g]
  exact climb_X0 hB401 (by omega) l (BBC (base t) l g) bs401_5 _ (by rw [card_bs401_5]; omega)

/-- pair 4 before its first gate layer, `Y`: ONE block of rank 48. -/
lemma g401Wgeom_aY4 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC g ∅ 0) (PhiC (RrC (base t) l (4 : Fin 8) g) (qC bs401_5) (1 + tt (I.tv t)))
      (fun φ => bcost φ 48) := by
  rw [hbase t, classT_Y0 (base t) l g, classT_Y (base t) l bs401_5 (4 : Fin 8) (Finset.mem_insert_self _ _) g]
  exact climb_Y0 hB401 (by omega) l (BBC (base t) l g) bs401_5 _ (by rw [card_bs401_5]; omega)

/-- pair 1 after stage 2, `X`: ONE final block of rank 48. -/
lemma g401Wgeom_fX1 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (1 : Fin 8) g) (qC bs401_2) 1) (kernel (BC.LC (Fin 8) H l))
      (fun φ => bcost φ 48) :=
  idle_Xend (base t) l hB401 (by omega) bs401_2 (1 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_2, Fintype.card_fin]; omega)

/-- pair 1 after stage 2, `Y`: ONE final block of rank 48. -/
lemma g401Wgeom_fY1 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (1 : Fin 8) g) (qC bs401_2) (1 + tt (I.tv t)))
      (frame (BBC (base t) l g) (univ \ {Sum.inl l}))
      (fun φ => bcost φ 48) := by
  rw [hbase t]
  exact idle_Yend (base t) l hB401 (by omega) bs401_2 (1 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_2, Fintype.card_fin]; omega)

/-- pair 2 after stage 4, `X`: ONE final block of rank 32. -/
lemma g401Wgeom_fX2 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (3 : Fin 8) g) (qC bs401_4) 1) (kernel (BC.LC (Fin 8) H l))
      (fun φ => bcost φ 32) :=
  idle_Xend (base t) l hB401 (by omega) bs401_4 (3 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_4, Fintype.card_fin]; omega)

/-- pair 2 after stage 4, `Y`: ONE final block of rank 32. -/
lemma g401Wgeom_fY2 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (3 : Fin 8) g) (qC bs401_4) (1 + tt (I.tv t)))
      (frame (BBC (base t) l g) (univ \ {Sum.inl l}))
      (fun φ => bcost φ 32) := by
  rw [hbase t]
  exact idle_Yend (base t) l hB401 (by omega) bs401_4 (3 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_4, Fintype.card_fin]; omega)

/-- pair 3 after stage 6, `X`: ONE final block of rank 16. -/
lemma g401Wgeom_fX3 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (5 : Fin 8) g) (qC bs401_6) 1) (kernel (BC.LC (Fin 8) H l))
      (fun φ => bcost φ 16) :=
  idle_Xend (base t) l hB401 (by omega) bs401_6 (5 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_6, Fintype.card_fin]; omega)

/-- pair 3 after stage 6, `Y`: ONE final block of rank 16. -/
lemma g401Wgeom_fY3 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (5 : Fin 8) g) (qC bs401_6) (1 + tt (I.tv t)))
      (frame (BBC (base t) l g) (univ \ {Sum.inl l}))
      (fun φ => bcost φ 16) := by
  rw [hbase t]
  exact idle_Yend (base t) l hB401 (by omega) bs401_6 (5 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_6, Fintype.card_fin]; omega)

/-- pair 4 after stage 8, `X`: ONE final block of rank 0. -/
lemma g401Wgeom_fX4 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (7 : Fin 8) g) (qC bs401_8) 1) (kernel (BC.LC (Fin 8) H l))
      (fun φ => bcost φ 0) :=
  idle_Xend (base t) l hB401 (by omega) bs401_8 (7 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_8, Fintype.card_fin]; omega)

/-- pair 4 after stage 8, `Y`: ONE final block of rank 0. -/
lemma g401Wgeom_fY4 (g : Orth (BC.LC (Fin 8) H l)) (t : T) :
    XReach (PhiC (RrC (base t) l (7 : Fin 8) g) (qC bs401_8) (1 + tt (I.tv t)))
      (frame (BBC (base t) l g) (univ \ {Sum.inl l}))
      (fun φ => bcost φ 0) := by
  rw [hbase t]
  exact idle_Yend (base t) l hB401 (by omega) bs401_8 (7 : Fin 8) (Finset.mem_insert_self _ _) g _
    (by rw [card_bs401_8, Fintype.card_fin]; omega)

/-- **The nine-stage word with entrance gauges and residual banks as a term of `G401W`** on the label space
`LC (Fin 8) H l`, for any gauged invocation package whose lines are the vectors of index `l` of orthonormal bases and
any bank schedule `S` with `load + pad = 73`, `pad < 73`. -/
def ggeom401W (S : BankSched Sl)
    (hres : ∀ b, S.load (fun q => Fintype.card H - (I.ent q).dim) b + S.pad b = 73)
    (hpad : ∀ b, S.pad b < 73)
    (hOk : ∀ (d : Orth (BC.LC (Fin 8) H l)) (q : Finset (Fin 8 × BC.Hp H l)), Ok (gwinC hB401 (by omega) d q)) :
    G401W H T Sl C (BC.LC (Fin 8) H l) (Orth (BC.LC (Fin 8) H l)) S.nb Ok where
  I := I
  ok0 := fun d => hOk d ∅
  ok1 := fun d => hOk d (qC bs401_1)
  ok2 := fun d => hOk d (qC bs401_2)
  ok3 := fun d => hOk d (qC bs401_3)
  ok4 := fun d => hOk d (qC bs401_4)
  ok5 := fun d => hOk d (qC bs401_5)
  ok6 := fun d => hOk d (qC bs401_6)
  ok7 := fun d => hOk d (qC bs401_7)
  ok8 := fun d => hOk d (qC bs401_8)
  G0 := fun d => gwinC hB401 (by omega) d ∅
  G1 := fun d => gwinC hB401 (by omega) d (qC bs401_1)
  G2 := fun d => gwinC hB401 (by omega) d (qC bs401_2)
  G3 := fun d => gwinC hB401 (by omega) d (qC bs401_3)
  G4 := fun d => gwinC hB401 (by omega) d (qC bs401_4)
  G5 := fun d => gwinC hB401 (by omega) d (qC bs401_5)
  G6 := fun d => gwinC hB401 (by omega) d (qC bs401_6)
  G7 := fun d => gwinC hB401 (by omega) d (qC bs401_7)
  G8 := fun d => gwinC hB401 (by omega) d (qC bs401_8)
  r1 := fun t => RrC (base t) l (0 : Fin 8)
  r2 := fun t => RrC (base t) l (1 : Fin 8)
  r3 := fun t => RrC (base t) l (2 : Fin 8)
  r4 := fun t => RrC (base t) l (3 : Fin 8)
  r5 := fun t => RrC (base t) l (4 : Fin 8)
  r6 := fun t => RrC (base t) l (5 : Fin 8)
  r7 := fun t => RrC (base t) l (6 : Fin 8)
  r8 := fun t => RrC (base t) l (7 : Fin 8)
  N0 := fun d q => entN (gwinC hB401 (by omega) d ∅).X (I.ent q)
  Wn0 := fun d q => entW (gwinC hB401 (by omega) d ∅).X (I.ent q)
  N1 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_1)).X (I.ent q)
  Wn1 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_1)).X (I.ent q)
  N2 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_2)).X (I.ent q)
  Wn2 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_2)).X (I.ent q)
  N3 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_3)).X (I.ent q)
  Wn3 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_3)).X (I.ent q)
  N4 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_4)).X (I.ent q)
  Wn4 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_4)).X (I.ent q)
  N5 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_5)).X (I.ent q)
  Wn5 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_5)).X (I.ent q)
  N6 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_6)).X (I.ent q)
  Wn6 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_6)).X (I.ent q)
  N7 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_7)).X (I.ent q)
  Wn7 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_7)).X (I.ent q)
  N8 := fun d q => entN (gwinC hB401 (by omega) d (qC bs401_8)).X (I.ent q)
  Wn8 := fun d q => entW (gwinC hB401 (by omega) d (qC bs401_8)).X (I.ent q)
  use := S.use
  hinj := S.hinj
  Zb := zb401W I l hH S hres hpad
  Zp := zp401W I l hH S hres hpad
  pad := S.pad
  u := fun g t => (BBC (base t) l g).v (Sum.inl l)
  YF := fun g t => frame (BBC (base t) l g) (univ \ {Sum.inl l})
  hm := two_le_card_LC hB401 (by omega)
  caX2 := fun φ => bcost φ 16
  caY2 := fun φ => bcost φ 16
  caX3 := fun φ => bcost φ 32
  caY3 := fun φ => bcost φ 32
  caX4 := fun φ => bcost φ 48
  caY4 := fun φ => bcost φ 48
  cX1 := fun φ => bcost φ 48
  cY1 := fun φ => bcost φ 48
  cX2 := fun φ => bcost φ 32
  cY2 := fun φ => bcost φ 32
  cX3 := fun φ => bcost φ 16
  cY3 := fun φ => bcost φ 16
  cX4 := fun φ => bcost φ 0
  cY4 := fun φ => bcost φ 0
  inv0 := fun g t => ⟨unframeM (g.1 * Matrix.fromBlocks (tt (I.tv t)) 0 0 (diagI ∅) * g.1ᵀ),
    frameM_unframeM _⟩
  z0 := fun g => PhiC_zero_empty g
  hN0 := fun d q => entN_spec _ _
  hW0 := fun d q => entW_spec _ _
  hN1 := fun d q => entN_spec _ _
  hW1 := fun d q => entW_spec _ _
  hN2 := fun d q => entN_spec _ _
  hW2 := fun d q => entW_spec _ _
  hN3 := fun d q => entN_spec _ _
  hW3 := fun d q => entW_spec _ _
  hN4 := fun d q => entN_spec _ _
  hW4 := fun d q => entW_spec _ _
  hN5 := fun d q => entN_spec _ _
  hW5 := fun d q => entW_spec _ _
  hN6 := fun d q => entN_spec _ _
  hW6 := fun d q => entW_spec _ _
  hN7 := fun d q => entN_spec _ _
  hW7 := fun d q => entW_spec _ _
  hN8 := fun d q => entN_spec _ _
  hW8 := fun d q => entW_spec _ _
  x01 := fun g t => g401Wgeom_x01 I base l hbase g t
  y01 := fun g t => g401Wgeom_y01 I base l hbase g t
  x12 := fun g t => g401Wgeom_xx I base l hbase bs401_1 0 1 (Finset.mem_singleton_self _) mem401_1 g t
  y12 := fun g t => g401Wgeom_yy I base l hbase bs401_1 0 1 (Finset.mem_singleton_self _) mem401_1 g t
  x23 := fun g t => g401Wgeom_xx I base l hbase bs401_2 1 2 (Finset.mem_insert_self _ _) mem401_2 g t
  y23 := fun g t => g401Wgeom_yy I base l hbase bs401_2 1 2 (Finset.mem_insert_self _ _) mem401_2 g t
  x34 := fun g t => g401Wgeom_xx I base l hbase bs401_3 2 3 (Finset.mem_insert_self _ _) mem401_3 g t
  y34 := fun g t => g401Wgeom_yy I base l hbase bs401_3 2 3 (Finset.mem_insert_self _ _) mem401_3 g t
  x45 := fun g t => g401Wgeom_xx I base l hbase bs401_4 3 4 (Finset.mem_insert_self _ _) mem401_4 g t
  y45 := fun g t => g401Wgeom_yy I base l hbase bs401_4 3 4 (Finset.mem_insert_self _ _) mem401_4 g t
  x56 := fun g t => g401Wgeom_xx I base l hbase bs401_5 4 5 (Finset.mem_insert_self _ _) mem401_5 g t
  y56 := fun g t => g401Wgeom_yy I base l hbase bs401_5 4 5 (Finset.mem_insert_self _ _) mem401_5 g t
  x67 := fun g t => g401Wgeom_xx I base l hbase bs401_6 5 6 (Finset.mem_insert_self _ _) mem401_6 g t
  y67 := fun g t => g401Wgeom_yy I base l hbase bs401_6 5 6 (Finset.mem_insert_self _ _) mem401_6 g t
  x78 := fun g t => g401Wgeom_xx I base l hbase bs401_7 6 7 (Finset.mem_insert_self _ _) mem401_7 g t
  y78 := fun g t => g401Wgeom_yy I base l hbase bs401_7 6 7 (Finset.mem_insert_self _ _) mem401_7 g t
  z0b := (zfr401W_spec I l hH S hres hpad).1
  hP0 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 0 g b
  hZ0 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 0 g q
  hZn0 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 0 g b hb
  hP1 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 1 g b
  hZ1 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 1 g q
  hZn1 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 1 g b hb
  hP2 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 2 g b
  hZ2 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 2 g q
  hZn2 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 2 g b hb
  hP3 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 3 g b
  hZ3 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 3 g q
  hZn3 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 3 g b hb
  hP4 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 4 g b
  hZ4 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 4 g q
  hZn4 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 4 g b hb
  hP5 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 5 g b
  hZ5 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 5 g q
  hZn5 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 5 g b hb
  hP6 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 6 g b
  hZ6 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 6 g q
  hZn6 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 6 g b hb
  hP7 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 7 g b
  hZ7 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 7 g q
  hZn7 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 7 g b hb
  hP8 := fun ρ g b => (zfr401W_spec I l hH S hres hpad).2.1 ρ 8 g b
  hZ8 := fun ρ g q => (zfr401W_spec I l hH S hres hpad).2.2.1 ρ 8 g q
  hZn8 := fun ρ g b hb => (zfr401W_spec I l hH S hres hpad).2.2.2.1 ρ 8 g b hb
  bpad := (zfr401W_spec I l hH S hres hpad).2.2.2.2
  aX2 := fun g t => g401Wgeom_aX2 I base l hbase hH g t
  aY2 := fun g t => g401Wgeom_aY2 I base l hbase hH g t
  aX3 := fun g t => g401Wgeom_aX3 I base l hbase hH g t
  aY3 := fun g t => g401Wgeom_aY3 I base l hbase hH g t
  aX4 := fun g t => g401Wgeom_aX4 I base l hbase hH g t
  aY4 := fun g t => g401Wgeom_aY4 I base l hbase hH g t
  fX1 := fun g t => g401Wgeom_fX1 I base l hbase hH g t
  fY1 := fun g t => g401Wgeom_fY1 I base l hbase hH g t
  fX2 := fun g t => g401Wgeom_fX2 I base l hbase hH g t
  fY2 := fun g t => g401Wgeom_fY2 I base l hbase hH g t
  fX3 := fun g t => g401Wgeom_fX3 I base l hbase hH g t
  fY3 := fun g t => g401Wgeom_fY3 I base l hbase hH g t
  fX4 := fun g t => g401Wgeom_fX4 I base l hbase hH g t
  fY4 := fun g t => g401Wgeom_fY4 I base l hbase hH g t
  term := fun g t => g401Wgeom_term I base l hbase g t

end Net

section End
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
variable {Ok : ∀ {β : Type} [Fintype β] [DecidableEq β], GWin H β → Prop}

/-- the unit price of the nine-stage word with entrance gauges and residual banks. -/
def unitPrice401W {nb : ℕ} (I : KLW.GInvEW H T Sl C Ok) (pad : Fin nb → ℕ) (φ : ℕ → ℝ) : ℝ :=
  72 * I.cost φ + (Fintype.card T : ℝ) * bankPrice377 φ + ∑ b, bcost φ (pad b)

/-- **Generalised certificate of the nine-stage word with entrance gauges and residual banks (PR #401)** for a
gauged invocation package and a bank schedule. -/
theorem bridge401W_gcert (I : KLW.GInvEW H T Sl C Ok) (hH : Fintype.card H = 9) (l : H) (S : BankSched Sl)
    (hres : S.fits I.ent)
    (hOk : ∀ (d : Orth (BC.LC (Fin 8) H l)) (q : Finset (Fin 8 × BC.Hp H l)), Ok (gwinC hB401 (by omega) d q)) :
    GCert (BC.LC (Fin 8) H l)
      (Sum.inl : BLive401 (Orth (BC.LC (Fin 8) H l)) T S.nb → BRole401 (Orth (BC.LC (Fin 8) H l)) T S.nb C)
      (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * unitPrice401W I S.pad φ) := by
  have hb : ∀ t, ∃ A : OBase H, I.tv t = A.v l :=
    fun t => base_ex (I.tv t) (I.hunit t) (I.hgap t) l
  choose base hbase using hb
  have hc := g401W_certificate (ggeom401W I base l hbase hH S hres.1 hres.2 hOk)
  have e : (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ)
        * (ggeom401W I base l hbase hH S hres.1 hres.2 hOk).unit φ)
      = (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * unitPrice401W I S.pad φ) := by
    funext φ
    simp only [G401W.unit, G401W.runUnit, ggeom401W, unitPrice401W, bankPrice377, hH, bcost_zero]
    first
      | ring
      | (norm_num; ring)
      | norm_num
  rw [e] at hc
  exact hc

end End
end
end GX
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.GX.ggeom401W
#print axioms OAI.PowerSaving.GX.bridge401W_gcert
