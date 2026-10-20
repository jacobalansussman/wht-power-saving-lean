"""Generator of Work/XHost/HGeom.lean (the hosted network on the label space LC (Fin 8) H l of Shea's
Bridge401/GGeom, with bank frames and host chains from XHost.host_frames) and Work/XHost/HEnd.lean (the end theorems).
usage: python3 -I gen_hgeom.py OUT_DIR
"""
import sys, os

OUT = sys.argv[1]
# name, k0, K, rank, X0 (X, Y), X1 (X, Y), geometry kind data
# geometry: ('A', bs index of end, c of end) | ('B', prev bs, prev c, bs, c) | ('F', bs, c)
KINDS = [
    ('A2', 0, 2, 16, ('A', 1, 0)),
    ('B2', 2, 1, 8, ('B', 1, 0, 2, 1)),
    ('A3', 0, 4, 32, ('A', 3, 2)),
    ('B3', 4, 1, 8, ('B', 3, 2, 4, 3)),
    ('A4', 0, 6, 48, ('A', 5, 4)),
    ('B4', 6, 1, 8, ('B', 5, 4, 6, 5)),
    ('F1', 3, 6, 48, ('F', 2, 1)),
    ('F2', 5, 4, 32, ('F', 4, 3)),
    ('F3', 7, 2, 16, ('F', 6, 5)),
]


def mem(bs, c):
    return '(Finset.mem_singleton_self _)' if bs == 1 else '(Finset.mem_insert_self _ _)'


def frames(geo):
    """(X0x, X0y, X1x, X1y) as Lean terms in g t"""
    if geo[0] == 'A':
        _, bs, c = geo
        return ('PhiC g ∅ (tt (I.tv t))', 'PhiC g ∅ 0',
                f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) 1',
                f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) (1 + tt (I.tv t))')
    if geo[0] == 'B':
        _, pb, pc, bs, c = geo
        return (f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) (tt (I.tv t))',
                f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) 0',
                f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) 1',
                f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) (1 + tt (I.tv t))')
    _, bs, c = geo
    return (f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) 1',
            f'PhiC (RrC (base t) l ({c} : Fin 8) g) (qC bs401_{bs}) (1 + tt (I.tv t))',
            'kernel (BC.LC (Fin 8) H l)',
            'frame (BBC (base t) l g) (univ \\ {Sum.inl l})')


def hform_lemma(name, rank, geo):
    x0x, x0y, x1x, x1y = frames(geo)
    L = []
    A = L.append
    A(f'/-- host form of the idle climb {name} (rank {rank}), both sides. -/')
    A(f'lemma hf_{name} (g : Orth (BC.LC (Fin 8) H l)) (j : Bool × T) :')
    A(f'    HostForm (x0_{name} I base l g j) (x1_{name} I base l g j) {rank} := by')
    A('  obtain ⟨s, t⟩ := j')
    A('  cases s')
    # X side
    A(f'  · show HostForm ({x0x}) ({x1x}) {rank}')
    if geo[0] == 'A':
        _, bs, c = geo
        A(f'    have e0 : {x0x} = frame (BBC (base t) l g) (JC {{l}} ∅) := by')
        A('      rw [hbase t]; exact classT_X0 (base t) l g')
        A(f'    have e1 : {x1x} = frame (BBC (base t) l g) (JC univ (qC bs401_{bs})) :=')
        A(f'      classT_X (base t) l bs401_{bs} ({c} : Fin 8) {mem(bs, c)} g')
        A('    rw [e0, e1]')
        A(f'    refine hform_of _ (Finset.disjSum_mono (Finset.subset_univ _) (Finset.empty_subset _)) {rank} ?_')
        A(f'    rw [card_JCX, card_JC, Finset.card_singleton, Finset.card_empty, card_bs401_{bs}, hH]')
    elif geo[0] == 'B':
        _, pb, pc, bs, c = geo
        A(f'    have e0 : {x0x} = frame (BBC (base t) l g) (JC univ (qC bs401_{pb})) :=')
        A(f'      (g401geom_xx I base l hbase bs401_{pb} ({pc} : Fin 8) ({c} : Fin 8) {mem(pb, pc)} mem401_{pb} g t).symm.trans')
        A(f'        (classT_X (base t) l bs401_{pb} ({pc} : Fin 8) {mem(pb, pc)} g)')
        A(f'    have e1 : {x1x} = frame (BBC (base t) l g) (JC univ (qC bs401_{bs})) :=')
        A(f'      classT_X (base t) l bs401_{bs} ({c} : Fin 8) {mem(bs, c)} g')
        A('    rw [e0, e1]')
        A(f'    refine hform_of _ (Finset.disjSum_mono (Finset.Subset.refl _) (qC_mono (Finset.subset_insert _ _))) {rank} ?_')
        A(f'    rw [card_JCX, card_JCX, card_bs401_{pb}, card_bs401_{bs}, hH]')
    else:
        _, bs, c = geo
        A(f'    have e0 : {x0x} = frame (BBC (base t) l g) (JC univ (qC bs401_{bs})) :=')
        A(f'      classT_X (base t) l bs401_{bs} ({c} : Fin 8) {mem(bs, c)} g')
        A('    rw [e0, ← frame_univ (BBC (base t) l g)]')
        A(f'    refine hform_of _ (Finset.subset_univ _) {rank} ?_')
        A(f'    rw [card_JCX, card_univ_LC, card_bs401_{bs}, hH, Fintype.card_fin]')
    # Y side
    A(f'  · show HostForm ({x0y}) ({x1y}) {rank}')
    if geo[0] == 'A':
        _, bs, c = geo
        A(f'    have e0 : {x0y} = frame (BBC (base t) l g) (JC ∅ ∅) := classT_Y0 (base t) l g')
        A(f'    have e1 : {x1y} = frame (BBC (base t) l g) (JC (univ \\ {{l}}) (qC bs401_{bs})) := by')
        A(f'      rw [hbase t]; exact classT_Y (base t) l bs401_{bs} ({c} : Fin 8) {mem(bs, c)} g')
        A('    rw [e0, e1]')
        A(f'    refine hform_of _ (Finset.disjSum_mono (Finset.empty_subset _) (Finset.empty_subset _)) {rank} ?_')
        A(f'    rw [card_JCY, card_JC, Finset.card_empty, Finset.card_empty, card_bs401_{bs}, hH]')
    elif geo[0] == 'B':
        _, pb, pc, bs, c = geo
        A(f'    have e0 : {x0y} = frame (BBC (base t) l g) (JC (univ \\ {{l}}) (qC bs401_{pb})) :=')
        A(f'      (g401geom_yy I base l hbase bs401_{pb} ({pc} : Fin 8) ({c} : Fin 8) {mem(pb, pc)} mem401_{pb} g t).symm.trans')
        A(f'        (by rw [hbase t]; exact classT_Y (base t) l bs401_{pb} ({pc} : Fin 8) {mem(pb, pc)} g)')
        A(f'    have e1 : {x1y} = frame (BBC (base t) l g) (JC (univ \\ {{l}}) (qC bs401_{bs})) := by')
        A(f'      rw [hbase t]; exact classT_Y (base t) l bs401_{bs} ({c} : Fin 8) {mem(bs, c)} g')
        A('    rw [e0, e1]')
        A(f'    refine hform_of _ (Finset.disjSum_mono (Finset.Subset.refl _) (qC_mono (Finset.subset_insert _ _))) {rank} ?_')
        A(f'    rw [card_JCY, card_JCY, card_bs401_{pb}, card_bs401_{bs}, hH]')
    else:
        _, bs, c = geo
        A(f'    have e0 : {x0y} = frame (BBC (base t) l g) (JC (univ \\ {{l}}) (qC bs401_{bs})) := by')
        A(f'      rw [hbase t]; exact classT_Y (base t) l bs401_{bs} ({c} : Fin 8) {mem(bs, c)} g')
        A('    rw [e0]')
        A(f'    refine hform_of _ (JCY_sub l _) {rank} ?_')
        A(f'    rw [card_JCY, Finset.card_sdiff_of_subset (Finset.subset_univ _), card_univ_LC, Finset.card_singleton,')
        A(f'      card_bs401_{bs}, hH, Fintype.card_fin]')
    return L


def emit():
    L = []
    A = L.append
    A('''import Work.XHost.HCert
import Work.XHost.HSched
import Work.Bridge401.End

/-!
# (key: bridge-host) The HOSTED nine-stage word on the label space `LC (Fin 8) H l` (geometry of Shea's GGeom)

Copy of the instance part of Work/Bridge401/GGeom.lean (Ryan Shea's complex588 tree; written by Nil Mamano's gen401.py
on jacklightChen's #377 geometry, Dugongue's #401 banks, Jacob Sussman's engine and stage geometry) for `GHost401`:

* `HostSched`: a host table `hostAt` (the same in every run; DaysSky's idle-axis hosting, PR #442), a bank map `bk` for
  the slots that are not hosted, the pads of the banks and the remainders `hrem` of the hosts; `HostSched.fits`: every
  bank's load plus pad is 73 (pad < 73), every host's load plus remainder is the rank of its idle climb.
* the bank frames and the host chains are ONE theorem, `XHost.host_frames` (h-L3-design's Frames.lean, after Shea's
  `bank_frames`): a bank is a host whose climb is `1 → kernel` (`hostForm_bank`), an idle climb inside one orthonormal
  basis has host form (`hostForm_frame`, from Sussman's `frame_core_free` / `core_nested`), and every gauged residual is
  `C * pk R * C'` (`entW_form`, Mamano / Shea).
* `ghostGeom`: the hosted network as a term of `GHost401`; `bridgeHost_gcert`: its `GCert`.
Written by h-L3-build/gen/gen_hgeom.py (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option maxHeartbeats 4000000
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX BC B401
noncomputable section

/-- **The residual conditions** of a host schedule for the entrance labels `ent`. -/
def HostSched.fits {Sl T H : Type} [Fintype Sl] [DecidableEq Sl] [DecidableEq T] [Fintype H] (S : HostSched Sl T)
    (ent : Sl → SLbl H) : Prop :=
  S.fitsR (fun q => Fintype.card H - (ent q).dim)

/-- HostForm of a nested pair of index sets of one orthonormal basis, with the rank of the difference. -/
lemma hform_of {α : Type} [Fintype α] [DecidableEq α] (A : OBase α) {s t : Finset α} (h : s ⊆ t) (k : ℕ)
    (hk : k + s.card = t.card) : HostForm (frame A s) (frame A t) k := by
  have e : (t \\ s).card = k := by
    have := Finset.card_sdiff_add_card_eq_card h
    omega
  rw [← e]
  exact hostForm_frame A h

/-! ## The idle climbs of the geometry have host form -/

section Geom
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
''')
    # frame definitions
    for name, k0, K, rank, geo in KINDS:
        x0x, x0y, x1x, x1y = frames(geo)
        A(f'/-- start and end frames of the idle climb {name}. -/')
        A(f'def x0_{name} (I : GInvE H T Sl C) (base : T → OBase H) (l : H) (g : Orth (BC.LC (Fin 8) H l)) (j : Bool × T) :\n    CMat (BC.LC (Fin 8) H l) :=')
        A(f'  match j with')
        A(f'  | (false, t) => {x0x}')
        A(f'  | (true, t) => {x0y}')
        A(f'def x1_{name} (I : GInvE H T Sl C) (base : T → OBase H) (l : H) (g : Orth (BC.LC (Fin 8) H l)) (j : Bool × T) :\n    CMat (BC.LC (Fin 8) H l) :=')
        A(f'  match j with')
        A(f'  | (false, t) => {x1x}')
        A(f'  | (true, t) => {x1y}')
        A('')
    A('variable (I : GInvE H T Sl C) (base : T → OBase H) (l : H)')
    A('variable (hbase : ∀ t, I.tv t = (base t).v l) (hH : Fintype.card H = 9)')
    A('include hbase hH')
    A('')
    for name, k0, K, rank, geo in KINDS:
        L.extend(hform_lemma(name, rank, geo))
        A('')
    A('end Geom')
    A('')
    # frames + instance
    A('''/-! ## The hosted network on `LC (Fin 8) H l` -/

section Net
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
variable (I : GInvE H T Sl C) (base : T → OBase H) (l : H)
  (hbase : ∀ t, I.tv t = (base t).v l) (hH : Fintype.card H = 9) (S : HostSched Sl T) (hfit : S.fits I.ent)

lemma hWnK (k0 : ℕ) : ∀ k : ℕ, ∀ d q, ∃ (C C' : CMat (BC.LC (Fin 8) H l)) (R : Finset (BC.LC (Fin 8) H l)),
    Free C ∧ Free C' ∧ R.card = Fintype.card H - (I.ent q).dim ∧ entW (win401N l hH (k0 + k) d).X (I.ent q) = C * pk R * C' :=
  fun k d q => entW_form _ (I.ent q)

lemma hWnB : ∀ k : ℕ, ∀ d q, ∃ (C C' : CMat (BC.LC (Fin 8) H l)) (R : Finset (BC.LC (Fin 8) H l)),
    Free C ∧ Free C' ∧ R.card = Fintype.card H - (I.ent q).dim ∧ entW (win401N l hH (k % 9) d).X (I.ent q) = C * pk R * C' :=
  fun k d q => entW_form _ (I.ent q)

include hbase hfit in
/-- the bank frames (`host_frames` with the banks as hosts of the climb `1 → kernel`). -/
lemma bank_ex : ∃ (Zb Zp : ℕ → Orth (BC.LC (Fin 8) H l) → Fin S.nb → CMat (BC.LC (Fin 8) H l)),
      (∀ d j, Zb 0 d j = 1) ∧
      (∀ k < 72, ∀ d j, GReach (Zb k d j) (Zp k d j) (fun _ => 0)) ∧
      (∀ k < 72, ∀ d q j, bUse S.hostAt S.bk k q = some j →
        entW (win401N l hH (k % 9) d).X (I.ent q) * Zp k d j = Zb (k + 1) d j) ∧
      (∀ k < 72, ∀ d j, (∀ q, bUse S.hostAt S.bk k q ≠ some j) → Zb (k + 1) d j = Zp k d j) ∧
      (∀ d j, GReach (Zb 72 d j) (kernel (BC.LC (Fin 8) H l)) (fun φ => bcost φ (S.pad j))) := by
  have hc : Fintype.card (BC.LC (Fin 8) H l) = 73 := card_LC401 l hH
  obtain ⟨Zb, Zp, h0, hP, hU, hN, hE⟩ := host_frames (bUse S.hostAt S.bk) (fun q => Fintype.card H - (I.ent q).dim)
    (fun k d q => entW (win401N l hH (k % 9) d).X (I.ent q)) 72
    (fun k hk => bUse_one S.hostAt S.bk S.hinj k hk) (fun k _ => hWnB I l hH k)
    (fun _ _ => (1 : CMat (BC.LC (Fin 8) H l))) (fun _ _ => kernel (BC.LC (Fin 8) H l)) (fun _ => 73)
    (fun _ _ => by rw [← hc]; exact hostForm_bank)
    (fun j => by show LH _ _ 72 j ≤ 73; have := hfit.1 j; omega)
    (fun j => by show 73 - LH _ _ 72 j < _; rw [hc]; have := hfit.1 j; have := hfit.2.1 j; omega)
  refine ⟨Zb, Zp, h0, hP, hU, hN, fun d j => ?_⟩
  have e : 73 - LH (bUse S.hostAt S.bk) (fun q => Fintype.card H - (I.ent q).dim) 72 j = S.pad j := by have := hfit.1 j; omega
  have h : GReach (Zb 72 d j) (kernel (BC.LC (Fin 8) H l))
      (fun φ => bcost φ (73 - LH (bUse S.hostAt S.bk) (fun q => Fintype.card H - (I.ent q).dim) 72 j)) := hE d j
  rw [e] at h
  exact h

def zbH := Classical.choose (bank_ex I base l hbase hH S hfit)
def zpH := Classical.choose (Classical.choose_spec (bank_ex I base l hbase hH S hfit))
def zH_spec := Classical.choose_spec (Classical.choose_spec (bank_ex I base l hbase hH S hfit))
''')
    for c, (name, k0, K, rank, geo) in enumerate(KINDS):
        A(f'''include hbase hfit in
/-- the host chains of the idle climb {name} (window: stages {k0}..{k0 + K - 1}, rank {rank}). -/
lemma chain_ex_{name} : ∃ (Zb Zp : ℕ → Orth (BC.LC (Fin 8) H l) → Bool × T → CMat (BC.LC (Fin 8) H l)),
      (∀ d j, Zb 0 d j = x0_{name} I base l d j) ∧
      (∀ k < {K}, ∀ d j, GReach (Zb k d j) (Zp k d j) (fun _ => 0)) ∧
      (∀ k < {K}, ∀ d q j, hUse S.hostAt {c} {k0} k q = some j →
        entW (win401N l hH ({k0} + k) d).X (I.ent q) * Zp k d j = Zb (k + 1) d j) ∧
      (∀ k < {K}, ∀ d j, (∀ q, hUse S.hostAt {c} {k0} k q ≠ some j) → Zb (k + 1) d j = Zp k d j) ∧
      (∀ d j, GReach (Zb {K} d j) (x1_{name} I base l d j) (fun φ => bcost φ (S.hrem {c} j))) := by
  have hc : Fintype.card (BC.LC (Fin 8) H l) = 73 := card_LC401 l hH
  have hf : ∀ j, LH (hUse S.hostAt {c} {k0}) (fun q => Fintype.card H - (I.ent q).dim) {K} j + S.hrem {c} j = {rank} := fun j => hfit.2.2 {c} j
  obtain ⟨Zb, Zp, h0, hP, hU, hN, hE⟩ := host_frames (hUse S.hostAt {c} {k0}) (fun q => Fintype.card H - (I.ent q).dim)
    (fun k d q => entW (win401N l hH ({k0} + k) d).X (I.ent q)) {K}
    (fun k _ => hUse_one S.hostAt S.hostInj {c} {k0} k) (fun k _ => hWnK I l hH {k0} k)
    (x0_{name} I base l) (x1_{name} I base l) (fun _ => {rank}) (fun d j => hf_{name} I base l hbase hH d j)
    (fun j => by show LH _ _ {K} j ≤ {rank}; have := hf j; omega)
    (fun j => by show {rank} - LH _ _ {K} j < _; rw [hc]; have := hf j; omega)
  refine ⟨Zb, Zp, h0, hP, hU, hN, fun d j => ?_⟩
  have e : {rank} - LH (hUse S.hostAt {c} {k0}) (fun q => Fintype.card H - (I.ent q).dim) {K} j = S.hrem {c} j := by have := hf j; omega
  have h : GReach (Zb {K} d j) (x1_{name} I base l d j)
      (fun φ => bcost φ ({rank} - LH (hUse S.hostAt {c} {k0}) (fun q => Fintype.card H - (I.ent q).dim) {K} j)) := hE d j
  rw [e] at h
  exact h

def zb_{name} := Classical.choose (chain_ex_{name} I base l hbase hH S hfit)
def zp_{name} := Classical.choose (Classical.choose_spec (chain_ex_{name} I base l hbase hH S hfit))
def z_{name}_spec := Classical.choose_spec (Classical.choose_spec (chain_ex_{name} I base l hbase hH S hfit))
''')
    # Hb / Hp
    A('/-- the host chains of all kinds, indexed by stage (`k - kz c` steps into the window). -/')
    for HZ, zz in (('hbH', 'zb'), ('hpH', 'zp')):
        A(f'def {HZ} (c : Fin 9) (k : ℕ) : Orth (BC.LC (Fin 8) H l) → Bool × T → CMat (BC.LC (Fin 8) H l) :=')
        A('  match c.val with')
        for c, (name, k0, K, rank, geo) in enumerate(KINDS):
            pat = str(c) if c < 8 else '_'
            A(f'  | {pat} => {zz}_{name} I base l hbase hH S hfit (k - {k0})')
        A('')
    # the instance
    A('''/-- **The hosted nine-stage word as a term of `GHost401`** on the label space `LC (Fin 8) H l`. -/
def ghostGeom : GHost401 H T Sl C (BC.LC (Fin 8) H l) (Orth (BC.LC (Fin 8) H l)) S.nb where
  I := I
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
  hostAt := S.hostAt
  hostInj := S.hostInj
  hostWin := S.hostWin
  bk := S.bk
  hinj := S.hinj
  Zb := zbH I base l hbase hH S hfit
  Zp := zpH I base l hbase hH S hfit
  pad := S.pad
  Hb := hbH I base l hbase hH S hfit
  Hp := hpH I base l hbase hH S hfit
  hrem := S.hrem
  u := fun g t => (BBC (base t) l g).v (Sum.inl l)
  YF := fun g t => frame (BBC (base t) l g) (univ \\ {Sum.inl l})
  hm := two_le_card_LC hB401 (by omega)
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
  x01 := fun g t => g401geom_x01 I base l hbase g t
  y01 := fun g t => g401geom_y01 I base l hbase g t
  x12 := fun g t => g401geom_xx I base l hbase bs401_1 0 1 (Finset.mem_singleton_self _) mem401_1 g t
  y12 := fun g t => g401geom_yy I base l hbase bs401_1 0 1 (Finset.mem_singleton_self _) mem401_1 g t
  x23 := fun g t => g401geom_xx I base l hbase bs401_2 1 2 (Finset.mem_insert_self _ _) mem401_2 g t
  y23 := fun g t => g401geom_yy I base l hbase bs401_2 1 2 (Finset.mem_insert_self _ _) mem401_2 g t
  x34 := fun g t => g401geom_xx I base l hbase bs401_3 2 3 (Finset.mem_insert_self _ _) mem401_3 g t
  y34 := fun g t => g401geom_yy I base l hbase bs401_3 2 3 (Finset.mem_insert_self _ _) mem401_3 g t
  x45 := fun g t => g401geom_xx I base l hbase bs401_4 3 4 (Finset.mem_insert_self _ _) mem401_4 g t
  y45 := fun g t => g401geom_yy I base l hbase bs401_4 3 4 (Finset.mem_insert_self _ _) mem401_4 g t
  x56 := fun g t => g401geom_xx I base l hbase bs401_5 4 5 (Finset.mem_insert_self _ _) mem401_5 g t
  y56 := fun g t => g401geom_yy I base l hbase bs401_5 4 5 (Finset.mem_insert_self _ _) mem401_5 g t
  x67 := fun g t => g401geom_xx I base l hbase bs401_6 5 6 (Finset.mem_insert_self _ _) mem401_6 g t
  y67 := fun g t => g401geom_yy I base l hbase bs401_6 5 6 (Finset.mem_insert_self _ _) mem401_6 g t
  x78 := fun g t => g401geom_xx I base l hbase bs401_7 6 7 (Finset.mem_insert_self _ _) mem401_7 g t
  y78 := fun g t => g401geom_yy I base l hbase bs401_7 6 7 (Finset.mem_insert_self _ _) mem401_7 g t
  z0b := (zH_spec I base l hbase hH S hfit).1''')
    for l in range(9):
        A(f'  hP{l} := fun ρ g b => (zH_spec I base l hbase hH S hfit).2.1 (9 * ρ.val + {l}) (by have := ρ.isLt; omega) g b')
        A(f'  hZ{l} := fun ρ g q hq => by')
        A(f'    have h := (zH_spec I base l hbase hH S hfit).2.2.1 (9 * ρ.val + {l}) (by have := ρ.isLt; omega) g q (S.bk ρ {l} q)')
        A(f'      (bUse_bank ρ {l} q hq)')
        A(f'    rw [show (9 * ρ.val + {l}) % 9 = {l} by omega] at h')
        A('    exact h')
        A(f'  hZn{l} := fun ρ g b hb => (zH_spec I base l hbase hH S hfit).2.2.2.1 (9 * ρ.val + {l}) (by have := ρ.isLt; omega) g b')
        A(f'    (fun q hq => by obtain ⟨h1, h2⟩ := bUse_to ρ {l} q b hq; exact hb q h1 h2)')
    A('  bpad := (zH_spec I base l hbase hH S hfit).2.2.2.2')
    for c, (name, k0, K, rank, geo) in enumerate(KINDS):
        sp = f'(z_{name}_spec I base l hbase hH S hfit)'
        A(f'  hH0_{name} := fun g t => ⟨{sp}.1 g (false, t), {sp}.1 g (true, t)⟩')
        for l in range(k0, k0 + K):
            k = l - k0
            A(f'  hHP_{name}_{l} := fun g j => {sp}.2.1 {k} (by decide) g j')
            A(f'  hHZ_{name}_{l} := fun g q j h => {sp}.2.2.1 {k} (by decide) g q j (hUse_of (c := {c}) (k0 := {k0}) (k := {k}) (by decide) h)')
            A(f'  hHZn_{name}_{l} := fun g j h => {sp}.2.2.2.1 {k} (by decide) g j')
            A(f'    (fun q hq => by obtain ⟨hk, e⟩ := hUse_to hq; exact h q e)')
        A(f'  hHE_{name} := fun g t => ⟨{sp}.2.2.2.2 g (false, t), {sp}.2.2.2.2 g (true, t)⟩')
    A('''  fX4 := fun g t => g401geom_fX4 I base l hbase hH g t
  fY4 := fun g t => g401geom_fY4 I base l hbase hH g t
  term := fun g t => g401geom_term I base l hbase g t

end Net

section End
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

/-- the price of the remaining climbs of kind `c` (both sides, all classes). -/
def hcostS {Sl T : Type} [Fintype Sl] [Fintype T] (S : HostSched Sl T) (c : Fin 9) (φ : ℕ → ℝ) : ℝ :=
  ∑ t : T, (bcost φ (S.hrem c (false, t)) + bcost φ (S.hrem c (true, t)))

/-- **the unit price of the hosted word**: 72 invocations, 8 runs of remaining climbs, one pad per bank. -/
def unitPriceHost (I : GInvE H T Sl C) (S : HostSched Sl T) (φ : ℕ → ℝ) : ℝ :=
  72 * I.cost φ + 8 * (hcostS S 0 φ + hcostS S 1 φ + hcostS S 2 φ + hcostS S 3 φ + hcostS S 4 φ + hcostS S 5 φ
    + hcostS S 6 φ + hcostS S 7 φ + hcostS S 8 φ) + ∑ b, bcost φ (S.pad b)

/-- **Generalised certificate of the HOSTED nine-stage word** for a gauged invocation package and a host schedule. -/
theorem bridgeHost_gcert (I : GInvE H T Sl C) (hH : Fintype.card H = 9) (l : H) (S : HostSched Sl T)
    (hfit : S.fits I.ent) :
    GCert (BC.LC (Fin 8) H l)
      (Sum.inl : BLive401 (Orth (BC.LC (Fin 8) H l)) T S.nb → BRole401 (Orth (BC.LC (Fin 8) H l)) T S.nb C)
      (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * unitPriceHost I S φ) := by
  have hb : ∀ t, ∃ A : OBase H, I.tv t = A.v l :=
    fun t => base_ex (I.tv t) (I.hunit t) (I.hgap t) l
  choose base hbase using hb
  have hc := ghost_certificate (ghostGeom I base l hbase hH S hfit)
  have e : (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ)
        * (ghostGeom I base l hbase hH S hfit).unit φ)
      = (fun φ => (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * unitPriceHost I S φ) := by
    funext φ
    simp only [GHost401.unit, GHost401.runUnit, GHost401.hcost, ghostGeom, unitPriceHost, hcostS, bcost_zero]
    ring
  rw [e] at hc
  exact hc

end End
end
end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.ghostGeom
#print axioms OAI.PowerSaving.XHost.bridgeHost_gcert
''')
    return '\n'.join(L) + '\n'


def emit_end():
    return '''import Work.XHost.HGeom

/-!
# (key: bridge-host) THE END THEOREMS for the HOSTED nine-stage word: gauged invocation package + host schedule + rate
fact => Walsh-Hadamard / kernel program (the form OpenAI's Fourier proof uses)

Copy of Work/Bridge401/End.lean (Ryan Shea's complex588 tree, written by DFT Lean (#401 L4) after Jacob Sussman's
Work/Bridge377/End.lean) with `bridgeHost_gcert` in the place of `bridge401_gcert` and the unit price
`unitPriceHost I S φ = 72 * I.cost φ + 8 * Σ_c Σ_t (bcost (hrem c (X,t)) + bcost (hrem c (Y,t))) + Σ_b bcost (pad b)`.
Hosting: DaysSky (CrocSwap/integer-mult-bounds PR #442).  Written by h-L3-build (wht43, 13 Oct 2026), prepared with
assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false

namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM RAM.Ty SS CB RF BR BG GF FoldRate GX
noncomputable section

universe U

section End
variable {H T Sl C : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]

/-- **Walsh-Hadamard transform from a gauged invocation package, a host schedule and a rate fact.** -/
theorem wht_of_ginv_host (I : GInvE H T Sl C) (hH : Fintype.card H = 9) (S : HostSched Sl T)
    (hfit : S.fits I.ent)
    (u W : ℕ) (hu : 73 = u)
    (hW : 64 * Fintype.card T + S.nb = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ, unitPriceHost I S φ = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z < 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ)) :
    ∃ solve W', WHT.WHTProgram solve W' ∧ BlockWHT.WHTTimeBoundsAt z W' := by
  obtain ⟨l⟩ : Nonempty H := Fintype.card_pos_iff.mp (by omega)
  exact BlockWHT.wht_main_of_ggroup_list u (by rw [BC.card_LC, hH, Fintype.card_fin]; omega) (by omega)
    (Sum.inl : BLive401 (Orth (BC.LC (Fin 8) H l)) T S.nb → BRole401 (Orth (BC.LC (Fin 8) H l)) T S.nb C)
    Sum.inl_injective (Fintype.card (Orth (BC.LC (Fin 8) H l))) W Orth.card_pos hW1 (by rw [card_BLive401, hW])
    _ (bridgeHost_gcert I hH l S hfit) U (fun φ => by rw [hU φ]) s z hz hz1 hnum

/-- **Kernel program from a gauged invocation package, a host schedule and a rate fact**: the form that OpenAI's
Fourier proof uses (`hills_program`), as `GX.engine_of_ginv_401`. -/
theorem engine_of_ginv_host (I : GInvE H T Sl C) (hH : Fintype.card H = 9) (S : HostSched Sl T)
    (hfit : S.fits I.ent)
    (u W : ℕ) (hu : 73 = u)
    (hW : 64 * Fintype.card T + S.nb = W) (hW1 : 1 ≤ W)
    (U : List (ℕ × ℕ))
    (hU : ∀ φ : ℕ → ℝ, unitPriceHost I S φ = listCost φ U)
    (s : ℕ) (z : ℝ) (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (hnum : ((2:ℝ)^s - 1) * BlockAccounting.moment u z U
        + (W:ℝ) * (BlockAccounting.T u z (u - 1) 1 + BlockAccounting.T u z 1 1)
        < (2:ℝ)^s * (W:ℝ))
    (cl : Paint) (m : Bool) {ι : Type U} (k : ι→ℕ) (v : ∀ i, Bits (k i) → ℂ) :
    Can m (fun i => 2^k i) (simIn cl) (bundle cl)
      (fun i => givenBits cl (k i) (v i)) (fun i => (bits (k i)).tape (ripple (k i) (v i)))
      (fun i => 2^k i*envelope z (k i)) := by
  obtain ⟨l⟩ : Nonempty H := Fintype.card_pos_iff.mp (by omega)
  exact engine_program_ggroup_list u (by rw [BC.card_LC, hH, Fintype.card_fin]; omega) (by omega)
    (Sum.inl : BLive401 (Orth (BC.LC (Fin 8) H l)) T S.nb → BRole401 (Orth (BC.LC (Fin 8) H l)) T S.nb C)
    Sum.inl_injective (Fintype.card (Orth (BC.LC (Fin 8) H l))) W Orth.card_pos hW1 (by rw [card_BLive401, hW])
    _ (bridgeHost_gcert I hH l S hfit) U (fun φ => by rw [hU φ]) s z hz hz1 hnum cl m k v

end End
end
end XHost
end PowerSaving
end OAI

#print axioms OAI.PowerSaving.XHost.wht_of_ginv_host
#print axioms OAI.PowerSaving.XHost.engine_of_ginv_host
'''


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'Work/XHost'), exist_ok=True)
    for k, v in (('HGeom', emit()), ('HEnd', emit_end())):
        p = os.path.join(OUT, 'Work/XHost', k + '.lean')
        open(p, 'w').write(v)
        print(p, v.count('\n'), 'lines')
