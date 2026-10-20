"""Generator of the HOSTED nine-stage network (same-run hosting, the idea of DaysSky's pull request #442) on Ryan Shea's complex588
tree (Work/Bridge401, written there by Nil Mamano's unpublished gen401.py; the 11 generated files are our templates).

Emits, under OUT/Work/XHost/:
  HRoles.lean  role maps of a hosted stage (slot -> bank or idle data role), injectivity, the lift of a stage route
  HDef.lean    structure GHost401 (G401 with hosted slots and host chains), prices, the 26 frame states of a run
  HStage.lean  adapters, stages, idle climbs (now: remaining climbs of the host chains), gates, finals of one run
  HCert.lean   run, 8 runs, pads, exchange, certificate (copy of GNetCert with the bank map `bk`)

Window kinds c (order of XHost.hostK / inWin): A2 B2 A3 B3 A4 B4 F1 F2 F3.
usage: python3 -I gen_hnet.py TEMPLATE_DIR OUT_DIR
"""
import sys, os, re

TPL, OUT = sys.argv[1], sys.argv[2]
PAIR = [1, 1, 1, 2, 2, 3, 3, 4, 4]
FWD = [True, False, True, False, True, False, True, False, True]
# name, pair, stages, start natural (X,Y) [N-terms], end stage e (A/B: frames of G_e (r_e t g) after it; F: kernel/YF)
KINDS = [
    ('A2', 2, [0, 1], ('N.xi0', 'N.yi0'), 1),
    ('B2', 2, [2], ('N.xi2', 'N.yi2'), 2),
    ('A3', 3, [0, 1, 2, 3], ('N.xi0', 'N.yi0'), 3),
    ('B3', 3, [4], ('N.xi4', 'N.yi4'), 4),
    ('A4', 4, [0, 1, 2, 3, 4, 5], ('N.xi0', 'N.yi0'), 5),
    ('B4', 4, [6], ('N.xi6', 'N.yi6'), 6),
    ('F1', 1, [3, 4, 5, 6, 7, 8], ('N.xo2', 'N.yo2'), None),
    ('F2', 2, [5, 6, 7, 8], ('N.xo4', 'N.yo4'), None),
    ('F3', 3, [7, 8], ('N.xo6', 'N.yo6'), None),
]
KN = [k[0] for k in KINDS]


def natexp(name, side, inside=False):
    """explicit Lean term of a natural frame N.xi2 etc. applied to g t (inside the structure: no `N.`)"""
    pre = '' if inside else 'N.'
    m = re.match(r'N\.(x|y)(i|o)(\d)', name)
    xy, io, l = m.group(1), m.group(2), int(m.group(3))
    G = f'({pre}G0 g)' if l == 0 else f'({pre}G{l} ({pre}r{l} t g))'
    if io == 'i':
        return f'{G}.Xm.Φ (tt ({pre}I.tv t))' if xy == 'x' else f'{G}.Xm.Φ 0'
    return f'{G}.Xm.Φ 1' if xy == 'x' else f'{G}.Xm.Φ (1 + tt ({pre}I.tv t))'


def window_of(p, l):
    for c, (name, pp, st, *_ ) in enumerate(KINDS):
        if pp == p and l in st:
            return c
    return None


# ---------------------------------------------------------------- states (Shea's natural states, hosting overlaid)
def shea_states():
    txt = open(os.path.join(TPL, 'GNetDef.lean')).read()
    S = {}
    for m in re.finditer(r'def G401\.S(\d+) \(ρ : Fin 8\) : BRun4 Γ T nb C → CMat α := mkR4 (.*)\n', txt):
        i = int(m.group(1)); rest = m.group(2)
        toks = []
        s = rest
        while s:
            s = s.strip()
            if s.startswith('('):
                d = 0
                for j, ch in enumerate(s):
                    d += ch == '('; d -= ch == ')'
                    if d == 0:
                        break
                toks.append(s[:j + 1]); s = s[j + 1:]
            else:
                sp = s.find(' ')
                if sp < 0:
                    toks.append(s); s = ''
                else:
                    toks.append(s[:sp]); s = s[sp:]
        assert len(toks) == 10, (i, toks)
        S[i] = toks
    assert len(S) == 26
    return S


STEPS = ['adapt0', 'stage0', 'adapt1', 'stage1', 'climbA2', 'gateA2', 'adapt2', 'stage2', 'climbB2', 'gateB2',
         'adapt3', 'stage3', 'climbA3', 'gateA3', 'adapt4', 'stage4', 'climbB3', 'gateB3', 'adapt5', 'stage5',
         'climbA4', 'gateA4', 'adapt6', 'stage6', 'climbB4', 'gateB4', 'adapt7', 'stage7', 'adapt8', 'stage8',
         'fin', 'gateC4', 'gateC3', 'gateC2']


def hosted_states():
    """state index before/after each step, frames (8 data entries) per state"""
    base = shea_states()
    cur = 0
    st = {0: list(base[0])}
    trans = []                       # (step, i_from, i_to)
    hb = lambda c, k, s: f'(fun g t => N.Hb {c} {k} g ({s}, t))'
    hp = lambda c, k, s: f'(fun g t => N.Hp {c} {k} g ({s}, t))'
    frames = list(base[0])
    for s in STEPS:
        if s.startswith('gate'):
            trans.append((s, cur, cur)); continue
        nxt = cur + 1
        fr = list(base[nxt])
        if s.startswith('adapt'):
            l = int(s[-1])
            for p in (1, 2, 3, 4):
                c = window_of(p, l)
                if c is not None and p != PAIR[l]:
                    fr[2 * (p - 1)] = hp(c, l, 'false'); fr[2 * (p - 1) + 1] = hp(c, l, 'true')
        elif s.startswith('stage'):
            l = int(s[-1])
            for p in (1, 2, 3, 4):
                c = window_of(p, l)
                if c is not None and p != PAIR[l]:
                    fr[2 * (p - 1)] = hb(c, l + 1, 'false'); fr[2 * (p - 1) + 1] = hb(c, l + 1, 'true')
        elif s.startswith('climb'):
            # the climbing pair goes to its natural frame (as Shea's); the others keep the frames of the stage before
            name = s[5:]; c = KN.index(name); p = KINDS[c][1]
            for pp in (1, 2, 3, 4):
                if pp != p:
                    fr[2 * (pp - 1)] = st[cur][2 * (pp - 1)]; fr[2 * (pp - 1) + 1] = st[cur][2 * (pp - 1) + 1]
        elif s == 'fin':
            pass
        st[nxt] = fr
        trans.append((s, cur, nxt))
        cur = nxt
    assert cur == 25
    return st, trans


# ---------------------------------------------------------------- HRoles.lean
def emit_roles():
    L = []
    A = L.append
    A('''import Work.Bridge401.GNetCert
import Work.XHost.HostK

/-!
# (key: bridge-host) Role maps of a HOSTED stage: a slot goes to a bank or to an idle data role of its run

Same-run hosting (DaysSky's idle-axis hosting, CrocSwap/integer-mult-bounds PR #442) on the nine-stage network of
Work/Bridge401 (Ryan Shea's complex588 tree; G401 by Nil Mamano on Dugongue #401 / jacklightChen #377 / Lydakis #275 /
DaysSky #363; engine, stage geometry and E8 unit by Jacob Sussman).  In stage `l` the slot `q` of invocation `g` lives in
`slotR hs u g q`: the bank `(g, u q)` when `hs q = none`, else the data role `(g, t)` of the pair `pairOf c`, side `s`
(`hs q = some (c, s, t)`, `c` the window kind of `XHost.inWin`).  `posH p hs u` is injective when the bank map is
injective on the non-hosted slots, no host is used twice, no host belongs to the acting pair `p`, and a stage sees at most
one window kind per pair (`HostOK`).  `ghliftQ<l>` is `GX.gbliftQ<l>` (Work/Bridge401/GStage.lean) for this embedding,
with the SAME scalar map `liftQ<l> u h` (the scalar maps of a stage never read the slot roles).
Written by h-L3-build/gen/gen_hnet.py (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX

/-- the pair of the window kind `c` (A2 B2 A3 B3 A4 B4 F1 F2 F3). -/
def pairOf (c : Fin 9) : ℕ :=
  match c.val with
  | 0 => 2
  | 1 => 2
  | 2 => 3
  | 3 => 3
  | 4 => 4
  | 5 => 4
  | 6 => 1
  | 7 => 2
  | _ => 3

/-- the pair acting in stage `l`. -/
def pairL (l : ℕ) : ℕ :=
  match l with
  | 0 => 1
  | 1 => 1
  | 2 => 1
  | 3 => 2
  | 4 => 2
  | 5 => 3
  | 6 => 3
  | 7 => 4
  | _ => 4

lemma win_pair : ∀ c l : Fin 9, inWin c.val l.val = true → pairOf c ≠ pairL l.val := by decide

lemma win_uniq : ∀ c c' l : Fin 9, inWin c.val l.val = true → inWin c'.val l.val = true →
    pairOf c = pairOf c' → c = c' := by decide

lemma pairOf_range : ∀ c : Fin 9, 1 ≤ pairOf c ∧ pairOf c ≤ 4 := by decide

lemma prng1 : 1 ≤ 1 ∧ 1 ≤ 4 := by decide
lemma prng2 : 1 ≤ 2 ∧ 2 ≤ 4 := by decide
lemma prng3 : 1 ≤ 3 ∧ 3 ≤ 4 := by decide
lemma prng4 : 1 ≤ 4 ∧ 4 ≤ 4 := by decide

noncomputable section

section Roles
variable {Γ T Sl C : Type} {nb : ℕ}

/-- the data role of pair `p` (1..4), side `s` (`false` = X, `true` = Y). -/
def dRole (p : ℕ) (s : Bool) (k : Γ × T) : BRun4 Γ T nb C :=
  match p, s with
  | 1, false => rqX1 k
  | 1, true => rqY1 k
  | 2, false => rqX2 k
  | 2, true => rqY2 k
  | 3, false => rqX3 k
  | 3, true => rqY3 k
  | _, false => rqX4 k
  | _, true => rqY4 k

/-- decoding of the data roles. -/
def dDec : BRun4 Γ T nb C → Option (ℕ × Bool × (Γ × T))
  | .inl (.inl (.inl (.inl (.inl (.inl k))))) => some (1, false, k)
  | .inl (.inl (.inl (.inl (.inl (.inr k))))) => some (1, true, k)
  | .inl (.inl (.inl (.inl (.inr (.inl k))))) => some (2, false, k)
  | .inl (.inl (.inl (.inl (.inr (.inr k))))) => some (2, true, k)
  | .inl (.inl (.inl (.inr (.inl k)))) => some (3, false, k)
  | .inl (.inl (.inl (.inr (.inr k)))) => some (3, true, k)
  | .inl (.inl (.inr (.inl k))) => some (4, false, k)
  | .inl (.inl (.inr (.inr k))) => some (4, true, k)
  | .inl (.inr _) => none
  | .inr _ => none

lemma dDec_dRole (p : ℕ) (hp : 1 ≤ p ∧ p ≤ 4) (s : Bool) (k : Γ × T) :
    dDec (dRole (nb := nb) (C := C) p s k) = some (p, s, k) := by
  obtain ⟨h1, h4⟩ := hp
  interval_cases p <;> cases s <;> rfl

lemma dRole_inj {p p' : ℕ} (hp : 1 ≤ p ∧ p ≤ 4) (hp' : 1 ≤ p' ∧ p' ≤ 4) {s s' : Bool} {k k' : Γ × T}
    (h : dRole (nb := nb) (C := C) p s k = dRole p' s' k') : p = p' ∧ s = s' ∧ k = k' := by
  have e := congrArg dDec h
  rw [dDec_dRole p hp, dDec_dRole p' hp'] at e
  simp only [Option.some.injEq, Prod.mk.injEq] at e
  exact e

lemma dRole_ne_S {p : ℕ} (hp : 1 ≤ p ∧ p ≤ 4) (s : Bool) (k : Γ × T) (b : Γ × Fin nb) :
    dRole (C := C) p s k ≠ rqS b := by
  intro h
  have e := congrArg dDec h
  rw [dDec_dRole p hp] at e
  exact absurd e (by simp [dDec, rqS])

lemma dRole_ne_C {p : ℕ} (hp : 1 ≤ p ∧ p ≤ 4) (s : Bool) (k : Γ × T) (b : Γ × C) :
    dRole (nb := nb) p s k ≠ rqC b := by
  intro h
  have e := congrArg dDec h
  rw [dDec_dRole p hp] at e
  exact absurd e (by simp [dDec, rqC])

/-- the role of the slot `q` of invocation `g` in a hosted stage. -/
def slotR (hs : Sl → Option (Fin 9 × (Bool × T))) (u : Sl → Fin nb) (g : Γ) (q : Sl) : BRun4 Γ T nb C :=
  match hs q with
  | some (c, s, t) => dRole (pairOf c) s (g, t)
  | none => rqS (g, u q)

lemma slotR_host {hs : Sl → Option (Fin 9 × (Bool × T))} {u : Sl → Fin nb} {g : Γ} {q : Sl} {c : Fin 9} {s : Bool}
    {t : T} (h : hs q = some (c, s, t)) : slotR (C := C) hs u g q = dRole (pairOf c) s (g, t) := by
  simp only [slotR, h]

lemma slotR_bank {hs : Sl → Option (Fin 9 × (Bool × T))} {u : Sl → Fin nb} {g : Γ} {q : Sl}
    (h : hs q = none) : slotR (C := C) hs u g q = rqS (g, u q) := by
  simp only [slotR, h]

/-- the roles of a hosted stage of pair `p` as roles of a run. -/
def posH (p : ℕ) (hs : Sl → Option (Fin 9 × (Bool × T))) (u : Sl → Fin nb) : SRole Γ T Sl C → BRun4 Γ T nb C :=
  mkS (fun g t => dRole p false (g, t)) (fun g t => dRole p true (g, t)) (fun g q => slotR hs u g q)
    (fun g c => rqC (g, c))

/-- the conditions on the slot table of one stage of pair `p`. -/
def HostOK (p : ℕ) (hs : Sl → Option (Fin 9 × (Bool × T))) (u : Sl → Fin nb) : Prop :=
  (∀ q q', hs q = none → hs q' = none → u q = u q' → q = q') ∧
  (∀ q q' j, hs q = some j → hs q' = some j → q = q') ∧
  (∀ q c j, hs q = some (c, j) → pairOf c ≠ p) ∧
  (∀ q q' c c' j j', hs q = some (c, j) → hs q' = some (c', j') → pairOf c = pairOf c' → c = c')

theorem posH_inj {p : ℕ} (hp : 1 ≤ p ∧ p ≤ 4) {hs : Sl → Option (Fin 9 × (Bool × T))} {u : Sl → Fin nb}
    (hok : HostOK p hs u) : Function.Injective (posH (Γ := Γ) (C := C) p hs u) := by
  obtain ⟨hu, hh, hP, hU⟩ := hok
  have hS : ∀ (g : Γ) (q : Sl), (hs q = none ∧ slotR (C := C) hs u g q = rqS (g, u q)) ∨
      (∃ c s t, hs q = some (c, s, t) ∧ slotR (C := C) hs u g q = dRole (pairOf c) s (g, t)) := by
    intro g q
    rcases hq : hs q with _ | ⟨c, s, t⟩
    · exact Or.inl ⟨rfl, slotR_bank hq⟩
    · exact Or.inr ⟨c, s, t, rfl, slotR_host hq⟩
  rintro (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩) (((⟨g',t'⟩|⟨g',t'⟩)|⟨g',q'⟩)|⟨g',c'⟩) h
  · change dRole p false (g, t) = dRole p false (g', t') at h
    obtain ⟨-, -, e⟩ := dRole_inj hp hp h
    rw [e]
  · change dRole p false (g, t) = dRole p true (g', t') at h
    exact absurd (dRole_inj hp hp h).2.1 (by decide)
  · change dRole p false (g, t) = slotR hs u g' q' at h
    rcases hS g' q' with ⟨-, e⟩ | ⟨c, s, t'', hq, e⟩
    · rw [e] at h
      exact absurd h (dRole_ne_S hp _ _ _)
    · rw [e] at h
      exact absurd (dRole_inj hp (pairOf_range c) h).1.symm (hP q' c (s, t'') hq)
  · change dRole p false (g, t) = rqC (g', c') at h
    exact absurd h (dRole_ne_C hp _ _ _)
  · change dRole p true (g, t) = dRole p false (g', t') at h
    exact absurd (dRole_inj hp hp h).2.1 (by decide)
  · change dRole p true (g, t) = dRole p true (g', t') at h
    obtain ⟨-, -, e⟩ := dRole_inj hp hp h
    rw [e]
  · change dRole p true (g, t) = slotR hs u g' q' at h
    rcases hS g' q' with ⟨-, e⟩ | ⟨c, s, t'', hq, e⟩
    · rw [e] at h
      exact absurd h (dRole_ne_S hp _ _ _)
    · rw [e] at h
      exact absurd (dRole_inj hp (pairOf_range c) h).1.symm (hP q' c (s, t'') hq)
  · change dRole p true (g, t) = rqC (g', c') at h
    exact absurd h (dRole_ne_C hp _ _ _)
  · change slotR hs u g q = dRole p false (g', t') at h
    rcases hS g q with ⟨-, e⟩ | ⟨c, s, t'', hq, e⟩
    · rw [e] at h
      exact absurd h.symm (dRole_ne_S hp _ _ _)
    · rw [e] at h
      exact absurd (dRole_inj (pairOf_range c) hp h).1 (hP q c (s, t'') hq)
  · change slotR hs u g q = dRole p true (g', t') at h
    rcases hS g q with ⟨-, e⟩ | ⟨c, s, t'', hq, e⟩
    · rw [e] at h
      exact absurd h.symm (dRole_ne_S hp _ _ _)
    · rw [e] at h
      exact absurd (dRole_inj (pairOf_range c) hp h).1 (hP q c (s, t'') hq)
  · change slotR hs u g q = slotR hs u g' q' at h
    rcases hS g q with ⟨hq, e⟩ | ⟨c, s, t, hq, e⟩ <;> rcases hS g' q' with ⟨hq', e'⟩ | ⟨c', s', t', hq', e'⟩ <;>
      rw [e, e'] at h
    · simp only [rqS, Sum.inl.injEq, Sum.inr.injEq, Prod.mk.injEq] at h
      obtain ⟨rfl, e2⟩ := h
      rw [hu q q' hq hq' e2]
    · exact absurd h.symm (dRole_ne_S (pairOf_range c') _ _ _)
    · exact absurd h (dRole_ne_S (pairOf_range c) _ _ _)
    · obtain ⟨e1, e2, e3⟩ := dRole_inj (pairOf_range c) (pairOf_range c') h
      simp only [Prod.mk.injEq] at e3
      obtain ⟨rfl, rfl⟩ := e3
      subst e2
      have ec := hU q q' c c' _ _ hq hq' e1
      subst ec
      rw [hh q q' _ hq hq']
  · change slotR hs u g q = rqC (g', c') at h
    rcases hS g q with ⟨-, e⟩ | ⟨c, s, t, -, e⟩
    · rw [e] at h
      exact absurd h (by simp [rqS, rqC])
    · rw [e] at h
      exact absurd h (dRole_ne_C (pairOf_range c) _ _ _)
  · change rqC (g, c) = dRole p false (g', t') at h
    exact absurd h.symm (dRole_ne_C hp _ _ _)
  · change rqC (g, c) = dRole p true (g', t') at h
    exact absurd h.symm (dRole_ne_C hp _ _ _)
  · change rqC (g, c) = slotR hs u g' q' at h
    rcases hS g' q' with ⟨-, e⟩ | ⟨c'', s, t, -, e⟩
    · rw [e] at h
      exact absurd h (by simp [rqS, rqC])
    · rw [e] at h
      exact absurd h.symm (dRole_ne_C (pairOf_range c'') _ _ _)
  · change rqC (g, c) = rqC (g', c') at h
    simp only [rqC, Sum.inr.injEq] at h
    rw [h]

def emH (p : ℕ) (hp : 1 ≤ p ∧ p ≤ 4) (hs : Sl → Option (Fin 9 × (Bool × T))) (u : Sl → Fin nb)
    (hok : HostOK p hs u) : SRole Γ T Sl C ↪ BRun4 Γ T nb C :=
  ⟨posH p hs u, posH_inj hp hok⟩

/-- the scalar maps of a stage never read the slot roles. -/
lemma GF_loc : ∀ y y' : SRole Γ T Sl C → ℂ, (∀ k, y (jX k) = y' (jX k)) → (∀ k, y (jY k) = y' (jY k)) →
    (∀ k, GF y (jX k) = GF y' (jX k)) ∧ (∀ k, GF y (jY k) = GF y' (jY k)) ∧ (∀ k, GF y (jC k) = GF y' (jC k)) := by
  intro y y' hX hY
  refine ⟨fun ⟨g, t⟩ => ?_, fun ⟨g, t⟩ => ?_, fun ⟨g, c⟩ => ?_⟩
  · exact hX (g, t)
  · show y (jY (g, t)) + y (jX (g, t)) = y' (jY (g, t)) + y' (jX (g, t))
    rw [hX, hY]
  · rfl

lemma GB_loc : ∀ y y' : SRole Γ T Sl C → ℂ, (∀ k, y (jX k) = y' (jX k)) → (∀ k, y (jY k) = y' (jY k)) →
    (∀ k, GB y (jX k) = GB y' (jX k)) ∧ (∀ k, GB y (jY k) = GB y' (jY k)) ∧ (∀ k, GB y (jC k) = GB y' (jC k)) := by
  intro y y' hX hY
  refine ⟨fun ⟨g, t⟩ => ?_, fun ⟨g, t⟩ => ?_, fun ⟨g, c⟩ => ?_⟩
  · show y (jX (g, t)) - y (jY (g, t)) = y' (jX (g, t)) - y' (jY (g, t))
    rw [hX, hY]
  · exact hY (g, t)
  · rfl
''')
    for l in range(9):
        p = PAIR[l]
        A(f'''/-- `liftQ{l}` keeps every data role outside pair {p} and every bank. -/
lemma liftQ{l}_fix (u : Sl → Fin nb) (h : (SRole Γ T Sl C → ℂ) → (SRole Γ T Sl C → ℂ)) (x : BRun4 Γ T nb C → ℂ)
    {{p' : ℕ}} (hp' : 1 ≤ p' ∧ p' ≤ 4) (hne : p' ≠ {p}) (s : Bool) (k : Γ × T) :
    liftQ{l} u h x (dRole p' s k) = x (dRole p' s k) := by
  obtain ⟨g, t⟩ := k
  obtain ⟨h1, h4⟩ := hp'
  interval_cases p' <;> cases s <;> first | rfl | exact absurd rfl hne
''')
    A('''end Roles

section Lift
variable {Γ T Sl C α : Type} {nb : ℕ} [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C]
''')
    for l in range(9):
        p = PAIR[l]
        A(f'''/-- **A hosted route of stage {l} of a run** (pair {p}): slots on banks or on idle data roles. -/
theorem ghliftQ{l} (hs : Sl → Option (Fin 9 × (Bool × T))) (u : Sl → Fin nb) (hok : HostOK {p} hs u)
    {{S0 T0 : BRun4 Γ T nb C → CMat α}}
    {{SX SY UX UY : Γ → T → CMat α}} {{SS US : Γ → Sl → CMat α}} {{SC UC : Γ → C → CMat α}}
    {{h : (SRole Γ T Sl C → ℂ) → (SRole Γ T Sl C → ℂ)}} {{c : (ℕ → ℝ) → ℝ}}
    (q : GRoute (mkS SX SY SS SC) (mkS UX UY US UC) h c)
    (hSX : ∀ g t, S0 (rqX{p} (g,t)) = SX g t) (hSY : ∀ g t, S0 (rqY{p} (g,t)) = SY g t)
    (hUX : ∀ g t, T0 (rqX{p} (g,t)) = UX g t) (hUY : ∀ g t, T0 (rqY{p} (g,t)) = UY g t)
    (hSS : ∀ g q, S0 (slotR hs u g q) = SS g q) (hUS : ∀ g q, T0 (slotR hs u g q) = US g q)
    (hSC : ∀ g c, S0 (rqC (g,c)) = SC g c) (hUC : ∀ g c, T0 (rqC (g,c)) = UC g c)
    (hfix : ∀ r, r ∉ covered (emH (Γ := Γ) (C := C) {p} prng{p} hs u hok) → S0 r = T0 r)
    (hh : ∀ y g q, h y (jS (g,q)) = y (jS (g,q)))
    (hloc : ∀ y y' : SRole Γ T Sl C → ℂ, (∀ k, y (jX k) = y' (jX k)) → (∀ k, y (jY k) = y' (jY k)) →
      (∀ k, h y (jX k) = h y' (jX k)) ∧ (∀ k, h y (jY k) = h y' (jY k)) ∧ (∀ k, h y (jC k) = h y' (jC k))) :
    GRoute S0 T0 (liftQ{l} u h) c := by
  refine GRoute.lift (emH (Γ := Γ) (C := C) {p} prng{p} hs u hok) (h := h) ?_
    (fun x => ?_) (fun x i hi => ?_) hfix
  · have e1 : S0 ∘ (emH (Γ := Γ) (C := C) {p} prng{p} hs u hok) = mkS SX SY SS SC := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · exact hSX g t
      · exact hSY g t
      · exact hSS g q
      · exact hSC g c
    have e2 : T0 ∘ (emH (Γ := Γ) (C := C) {p} prng{p} hs u hok) = mkS UX UY US UC := by
      funext r
      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
      · exact hUX g t
      · exact hUY g t
      · exact hUS g q
      · exact hUC g c
    rw [e1, e2]
    exact q
  · have hl := hloc (x ∘ posQ{l} u) (x ∘ posH {p} hs u) (fun ⟨g, t⟩ => rfl) (fun ⟨g, t⟩ => rfl)
    funext r
    rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)
    · exact hl.1 (g, t)
    · exact hl.2.1 (g, t)
    · show liftQ{l} u h x (slotR hs u g q) = h (x ∘ posH {p} hs u) (jS (g, q))
      rw [hh]
      show liftQ{l} u h x (slotR hs u g q) = x (slotR hs u g q)
      rcases hq : hs q with _ | ⟨c, s, t⟩
      · rw [slotR_bank hq]
        rfl
      · rw [slotR_host hq]
        exact liftQ{l}_fix u h x (pairOf_range c) (hok.2.2.1 q c (s, t) hq) s (g, t)
    · exact hl.2.2 (g, c)
  · rcases i with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)
''')
        for pp in (1, 2, 3, 4):
            for side in ('X', 'Y'):
                if pp == p:
                    j = 'jX' if side == 'X' else 'jY'
                    A(f'    · exact absurd ((cover_iff _ _).mpr ⟨{j} (g,t), rfl⟩) hi')
                else:
                    A('    · rfl')
        A('    · rfl')
        A('    · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hi')
        A('')
    A('''end Lift
end
end XHost
end PowerSaving
end OAI
''')
    for l in range(9):
        A(f'#print axioms OAI.PowerSaving.XHost.ghliftQ{l}')
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------- HDef.lean
def fields():
    out = []
    A = out.append
    A('  /-- the gauged invocation (`GInvE`). -/')
    A('  I : GInvE H T Sl C')
    for l in range(9):
        A(f'  G{l} : Γ → GWin H α')
    for l in range(1, 9):
        A(f'  r{l} : T → Γ ≃ Γ')
    for l in range(9):
        A(f'  N{l} : Γ → Sl → CMat α')
        A(f'  Wn{l} : Γ → Sl → CMat α')
    A('  /-- the host of slot `q` in stage `l`, the SAME in every run: (window kind, side, class); `none`: a bank. -/')
    A('  hostAt : Fin 9 → Sl → Option (Fin 9 × (Bool × T))')
    A('  hostInj : ∀ l q q\' j, hostAt l q = some j → hostAt l q\' = some j → q = q\'')
    A('  hostWin : ∀ l q c j, hostAt l q = some (c, j) → inWin c.val l.val = true')
    A('  /-- the bank of a NON-hosted slot `q` in stage `l` of run `ρ` (junk on hosted slots). -/')
    A('  bk : Fin 8 → Fin 9 → Sl → Fin nb')
    A('  hinj : ∀ ρ l q q\', hostAt l q = none → hostAt l q\' = none → bk ρ l q = bk ρ l q\' → q = q\'')
    A('  Zb : ℕ → Γ → Fin nb → CMat α')
    A('  Zp : ℕ → Γ → Fin nb → CMat α')
    A('  pad : Fin nb → ℕ')
    A('  /-- the host chains: frame of host `j` of kind `c` before (`Hb`) and after (`Hp`) the adapter of stage `k`. -/')
    A('  Hb : Fin 9 → ℕ → Γ → Bool × T → CMat α')
    A('  Hp : Fin 9 → ℕ → Γ → Bool × T → CMat α')
    A('  hrem : Fin 9 → Bool × T → ℕ')
    A('  u : Γ → T → Space α')
    A('  YF : Γ → T → CMat α')
    A('  hm : 2 ≤ Fintype.card α')
    A('  cX4 : (ℕ → ℝ) → ℝ')
    A('  cY4 : (ℕ → ℝ) → ℝ')
    A('  inv0 : ∀ g t, ∃ N : CMat α, (G0 g).Xm.Φ (tt (I.tv t)) * N = 1')
    A('  z0 : ∀ g, (G0 g).Xm.Φ 0 = 1')
    for l in range(9):
        A(f'  hN{l} : ∀ d q, (G{l} d).X.Φ (I.ent q) * N{l} d q = 1')
        A(f'  hW{l} : ∀ d q, (G{l} d).X.Φ (I.ent q).top = Wn{l} d q * (G{l} d).X.Φ (I.ent q)')
    prev = lambda l: 'g' if l == 0 else f'(r{l} t g)'
    for l in range(8):
        a = f'(G{l} {prev(l)})'
        b = f'(G{l+1} (r{l+1} t g))'
        A(f'  x{l}{l+1} : ∀ g t, {a}.Xm.Φ 1 = {b}.Xm.Φ (tt (I.tv t))')
        A(f'  y{l}{l+1} : ∀ g t, {a}.Xm.Φ (1 + tt (I.tv t)) = {b}.Xm.Φ 0')
    A('  /-- the banks (non-hosted slots only). -/')
    A('  z0b : ∀ g b, Zb 0 g b = 1')
    for l in range(9):
        A(f'  hP{l} : ∀ (ρ : Fin 8) g b, GReach (Zb (9 * ρ.val + {l}) g b) (Zp (9 * ρ.val + {l}) g b) (fun _ => 0)')
        A(f'  hZ{l} : ∀ (ρ : Fin 8) g q, hostAt {l} q = none → Wn{l} g q * Zp (9 * ρ.val + {l}) g (bk ρ {l} q) = Zb (9 * ρ.val + {l + 1}) g (bk ρ {l} q)')
        A(f'  hZn{l} : ∀ (ρ : Fin 8) g b, (∀ q, hostAt {l} q = none → bk ρ {l} q ≠ b) → Zb (9 * ρ.val + {l + 1}) g b = Zp (9 * ρ.val + {l}) g b')
    A('  bpad : ∀ g b, GReach (Zb 72 g b) (kernel α) (fun φ => bcost φ (pad b))')
    A('  /-- the host chains (XHost.host_frames): start, free adapters, hosted residuals, ONE remaining block. -/')
    for c, (name, p, st, (sx, sy), e) in enumerate(KINDS):
        k0, k1 = st[0], st[-1] + 1
        A(f'  hH0_{name} : ∀ g t, Hb {c} {k0} g (false, t) = {natexp(sx, "x", True)} ∧ Hb {c} {k0} g (true, t) = {natexp(sy, "y", True)}')
        for l in st:
            A(f'  hHP_{name}_{l} : ∀ g j, GReach (Hb {c} {l} g j) (Hp {c} {l} g j) (fun _ => 0)')
            A(f'  hHZ_{name}_{l} : ∀ g q j, hostAt {l} q = some ({c}, j) → Wn{l} g q * Hp {c} {l} g j = Hb {c} {l + 1} g j')
            A(f'  hHZn_{name}_{l} : ∀ g j, (∀ q, hostAt {l} q ≠ some ({c}, j)) → Hb {c} {l + 1} g j = Hp {c} {l} g j')
        if e is not None:
            ex, ey = f'(G{e} (r{e} t g)).Xm.Φ 1', f'(G{e} (r{e} t g)).Xm.Φ (1 + tt (I.tv t))'
        else:
            ex, ey = 'kernel α', 'YF g t'
        A(f'  hHE_{name} : ∀ g t, GReach (Hb {c} {k1} g (false, t)) ({ex}) (fun φ => bcost φ (hrem {c} (false, t))) ∧')
        A(f'    GReach (Hb {c} {k1} g (true, t)) ({ey}) (fun φ => bcost φ (hrem {c} (true, t)))')
    A('  /-- pair 4 after stage 8: ONE final block (not hosted). -/')
    A('  fX4 : ∀ g t, XReach ((G8 (r8 t g)).Xm.Φ 1) (kernel α) cX4')
    A('  fY4 : ∀ g t, XReach ((G8 (r8 t g)).Xm.Φ (1 + tt (I.tv t))) (YF g t) cY4')
    A('  term : ∀ g t, kernel α * (G0 g).Xm.Φ (tt (I.tv t)) = shift (u g t) * YF g t')
    return out


def emit_def():
    st, trans = hosted_states()
    L = []
    A = L.append
    A('''import Work.XHost.HRoles

/-!
# (key: bridge-host) The nine-stage word of PR #401 with SAME-RUN HOSTING: data and frames

Copy of Work/Bridge401/GNetDef.lean (Ryan Shea's complex588 tree; `G401` by Nil Mamano, written by his gen401.py, on
Dugongue #401 / jacklightChen #377 / Lydakis #275 / DaysSky #363; engine and E8 unit by Jacob Sussman) with idle-axis
hosting (DaysSky, CrocSwap/integer-mult-bounds PR #442): in stage `l` the slot `q` goes to the idle data role
`hostAt l q = some (c, s, t)` (window kind `c`, side `s`, class `(g, t)` of the run's replica; the same table in every
run) or, when `hostAt l q = none`, to the bank `bk ρ l q`.  The idle climbs A2..A4, B2..B4 and the finals F1..F3 become
host chains `Hb`/`Hp`: free adapters, the hosted residuals, ONE remaining block of rank `hrem c j` at the old place.
Written by h-L3-build/gen/gen_hnet.py (wht43, 13 Oct 2026), prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX
noncomputable section

structure GHost401 (H T Sl C α Γ : Type) (nb : ℕ) [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
    [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α] where''')
    L.extend(fields())
    A('''
section Net
variable {H T Sl C α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {nb : ℕ} (N : GHost401 H T Sl C α Γ nb)

/-- price of the remaining climbs of the host chains of kind `c` (both sides, all classes), per element of `Γ`. -/
def GHost401.hcost (c : Fin 9) (φ : ℕ → ℝ) : ℝ :=
  ∑ t : T, (bcost φ (N.hrem c (false, t)) + bcost φ (N.hrem c (true, t)))

/-- price of ONE run per element of `Γ`: 9 invocations, the remaining climbs, pair 4's final. -/
def GHost401.runUnit (φ : ℕ → ℝ) : ℝ :=
  9 * N.I.cost φ + (N.hcost 0 φ + N.hcost 1 φ + N.hcost 2 φ + N.hcost 3 φ + N.hcost 4 φ + N.hcost 5 φ
    + N.hcost 6 φ + N.hcost 7 φ + N.hcost 8 φ) + (Fintype.card T : ℝ) * (N.cX4 φ + N.cY4 φ)

/-- price of one unit (one element of `Γ`): 8 runs and the padding blocks of the banks. -/
def GHost401.unit (φ : ℕ → ℝ) : ℝ :=
  8 * N.runUnit φ + ∑ b, bcost φ (N.pad b)

/-- the conditions of the role map of stage `l`. -/
lemma GHost401.hok (ρ : Fin 8) (l : Fin 9) : HostOK (pairL l.val) (N.hostAt l) (N.bk ρ l) :=
  ⟨N.hinj ρ l, N.hostInj l, fun q c j h => win_pair c l (N.hostWin l q c j h),
    fun q q' c c' j j' h h' e => win_uniq c c' l (N.hostWin l q c j h) (N.hostWin l q' c' j' h') e⟩
''')
    for l in range(9):
        r = 'g' if l == 0 else f'(N.r{l} t g)'
        A(f'def GHost401.xi{l} (g : Γ) (t : T) : CMat α := (N.G{l} {r}).Xm.Φ (tt (N.I.tv t))')
        A(f'def GHost401.yi{l} (g : Γ) (t : T) : CMat α := (N.G{l} {r}).Xm.Φ 0')
    for l in (2, 4, 6, 8):
        A(f'def GHost401.xo{l} (g : Γ) (t : T) : CMat α := (N.G{l} (N.r{l} t g)).Xm.Φ 1')
        A(f'def GHost401.yo{l} (g : Γ) (t : T) : CMat α := (N.G{l} (N.r{l} t g)).Xm.Φ (1 + tt (N.I.tv t))')
    A('''
/-- the data of a replica before its run and after it. -/
def GHost401.InitD : BDat4 Γ T → CMat α := mkD4 N.xi0 N.yi0 N.xi0 N.yi0 N.xi0 N.yi0 N.xi0 N.yi0
def GHost401.FinD : BDat4 Γ T → CMat α := mkD4 (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF (fun _ _ => kernel α) N.YF
''')
    for i in range(26):
        fr = st[i]
        A(f'def GHost401.S{i} (ρ : Fin 8) : BRun4 Γ T nb C → CMat α := mkR4 {" ".join(fr[:8])} {fr[8]} {fr[9]}')
    A('''
/-- frames of all roles at the start of run `k` (`k = 8`: after the last run). -/
def GHost401.Sfull (k : Fin 9) : BRole401 Γ T nb C → CMat α :=
  mk401 (fun ρ' d => if ρ'.val < k.val then N.FinD d else N.InitD d) (fun g b => N.Zb (9 * k.val) g b) (cOne Γ C α)

/-- frames of all roles after the padding climbs. -/
def GHost401.Sfin : BRole401 Γ T nb C → CMat α :=
  mk401 (fun _ d => N.FinD d) (fun _ _ => kernel α) (cOne Γ C α)

end Net
end
end XHost
end PowerSaving
end OAI
''')
    return '\n'.join(L) + '\n', st, trans


# ---------------------------------------------------------------- HStage.lean
RCASES_R = 'rintro ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)'
RCASES_r = 'rcases r with ((((((⟨g,t⟩|⟨g,t⟩)|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|(⟨g,t⟩|⟨g,t⟩))|⟨g,b⟩)|⟨g,c⟩)'
ROLES = [(1, 'false'), (1, 'true'), (2, 'false'), (2, 'true'), (3, 'false'), (3, 'true'), (4, 'false'), (4, 'true')]


def emit_stage(st, trans):
    L = []
    A = L.append
    A('''import Work.XHost.HDef

/-!
# (key: bridge-host) The steps of one run of the HOSTED nine-stage word

Copy of Work/Bridge401/GNetStage.lean (Shea's complex588 tree; written by Nil Mamano's gen401.py) with same-run hosting:
the adapters also move every host of a window (`hHP`), a stage puts its slots on banks or hosts (`ghliftQ<l>`), the
idle climbs are the remaining climbs of the host chains (`hHE`), and the finals of pairs 1-3 are those of F1-F3.
The gate layers are Shea's, word for word.  Written by h-L3-build/gen/gen_hnet.py, prepared with assistance from
Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
set_option linter.unusedSimpArgs false
set_option maxHeartbeats 4000000
namespace OAI
namespace PowerSaving
namespace XHost
open Binary Matrix Finset RAM SS CB RF BR GF GX
noncomputable section

section Net
variable {H T Sl C α Γ : Type} [Fintype H] [DecidableEq H] [Fintype T] [DecidableEq T]
  [Fintype Sl] [DecidableEq Sl] [Fintype C] [DecidableEq C] [Fintype α] [DecidableEq α]
  [Fintype Γ] [DecidableEq Γ]
variable {nb : ℕ} (N : GHost401 H T Sl C α Γ nb)
''')
    tpl = open(os.path.join(TPL, 'GNetStage.lean')).read()
    nat_after = {0: ('N.xi1', 'N.yi1'), 1: ('N.xi2', 'N.yi2'), 2: ('N.xo2', 'N.yo2'), 3: ('N.xi4', 'N.yi4'), 4: ('N.xo4', 'N.yo4'),
                 5: ('N.xi6', 'N.yi6'), 6: ('N.xo6', 'N.yo6'), 7: ('N.xi8', 'N.yi8'), 8: ('N.xo8', 'N.yo8')}
    has_next = {0: True, 1: True, 2: False, 3: True, 4: False, 5: True, 6: False, 7: True, 8: False}
    for s, a, b in trans:
        if s.startswith('adapt'):
            l = int(s[-1])
            A(f'/-- the free adapters before stage {l} of run `ρ`: every bank from `Zb` to `Zp`, every host of a window. -/')
            A(f'lemma GHost401.adapt{l} (ρ : Fin 8) :')
            A(f'    GRoute (N.S{a} ρ) (N.S{b} ρ) id (fun _ => 0) := by')
            A(f'  have H1 := GRoute.reach_all (N.S{a} ρ) (N.S{b} ρ) (fun _ _ => 0) (by')
            A(f'      {RCASES_R}')
            for p, side in ROLES:
                c = window_of(p, l)
                if c is None or p == PAIR[l]:
                    A('      · exact GReach.refl _')
                    continue
                name, _, stw, (sx, sy), _ = KINDS[c]
                if l == stw[0]:
                    nat = natexp(sx if side == 'false' else sy, 'x' if side == 'false' else 'y')
                    comp = '1' if side == 'false' else '2'
                    A(f'      · show GReach ({nat}) (N.Hp {c} {l} g ({side}, t)) (fun _ => 0)')
                    A(f'        rw [← (N.hH0_{name} g t).{comp}]')
                    A(f'        exact N.hHP_{name}_{l} g ({side}, t)')
                else:
                    A(f'      · exact N.hHP_{name}_{l} g ({side}, t)')
            A(f'      · exact N.hP{l} ρ g b')
            A('      · exact GReach.refl _)')
            A('  exact H1.cast (fun φ => by simp)')
            A('')
        elif s.startswith('stage'):
            l = int(s[-1]); p = PAIR[l]
            fb = 'GF' if FWD[l] else 'GB'
            fwd = 'fwd' if FWD[l] else 'bwd'
            eps = '(bid Γ T)' if l == 0 else f'N.r{l}'
            nx, ny = nat_after[l]
            sl = f'slotR (N.hostAt {l}) (N.bk ρ {l}) g q'
            A(f'/-- stage {l} of run `ρ` (pair {p}): the slot `q` of `g` on its bank or on its host. -/')
            A(f'lemma GHost401.stage{l} (ρ : Fin 8) :')
            A(f'    GRoute (N.S{a} ρ) (N.S{b} ρ) (liftQ{l} (Γ := Γ) (T := T) (C := C) (N.bk ρ {l}) {fb})')
            A(f'      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by')
            A(f'  have hUS : ∀ g q, N.Wn{l} g q * N.S{a} ρ ({sl}) = N.S{b} ρ ({sl}) := by')
            A('    intro g q')
            A(f'    rcases hq : N.hostAt {l} q with _ | ⟨c, s, t⟩')
            A('    · rw [slotR_bank hq]')
            A(f'      exact N.hZ{l} ρ g q hq')
            A('    · rw [slotR_host hq]')
            A(f'      have hw := N.hostWin {l} q c (s, t) hq')
            okc = {}
            for p2 in (1, 2, 3, 4):
                c2 = window_of(p2, l)
                if c2 is not None and p2 != p:
                    okc[c2] = KINDS[c2][0]
            A('      fin_cases c')
            for c2 in range(9):
                if c2 in okc:
                    A(f'      · cases s <;> exact N.hHZ_{okc[c2]}_{l} g q _ hq')
                else:
                    A('      · exact absurd hw (by decide)')
            A(f'  have h := gistageE_{fwd} N.I {eps} N.G{l} N.N{l} N.Wn{l} N.hN{l} N.hW{l}')
            A(f'    (fun g q => N.S{a} ρ ({sl})) (cOne Γ C α) (cOne Γ C α)')
            A(f'  have h\' : GRoute (bIn N.I.tv {eps} (fun d => (N.G{l} d).Xm)')
            A(f'        (fun g q => N.S{a} ρ ({sl})) (cOne Γ C α))')
            A(f'      (mkS {nx} {ny} (fun g q => N.S{b} ρ ({sl})) (cOne Γ C α)) {fb}')
            A('      (fun φ => (Fintype.card Γ : ℝ) * N.I.cost φ) := by')
            A(f'    have e : bOut N.I.tv {eps} (fun d => (N.G{l} d).Xm)')
            A(f'        (fun g q => N.Wn{l} g q * N.S{a} ρ ({sl})) (cOne Γ C α)')
            A(f'        = mkS {nx} {ny} (fun g q => N.S{b} ρ ({sl})) (cOne Γ C α) := by')
            A('      funext r')
            A('      rcases r with (((⟨g,t⟩|⟨g,t⟩)|⟨g,q⟩)|⟨g,c⟩)')
            if has_next[l]:
                A(f'      · exact N.x{l}{l+1} g t')
                A(f'      · exact N.y{l}{l+1} g t')
            else:
                A('      · rfl')
                A('      · rfl')
            A('      · exact hUS g q')
            A('      · rfl')
            A('    rw [← e]')
            A('    exact h')
            A(f'  refine ghliftQ{l} (N.hostAt {l}) (N.bk ρ {l}) (N.hok ρ {l}) h\' (fun g t => rfl) (fun g t => rfl)')
            A('    (fun g t => rfl) (fun g t => rfl) (fun g q => rfl) (fun g q => rfl) (fun g c => rfl) (fun g c => rfl) ?_')
            A(f'    (fun _ _ _ => rfl) {fb}_loc')
            A('  intro r hr')
            A(f'  {RCASES_r}')
            for p2, side in ROLES:
                if p2 == p:
                    j = 'jX' if side == 'false' else 'jY'
                    A(f'  · exact absurd ((cover_iff _ _).mpr ⟨{j} (g,t), rfl⟩) hr')
                    continue
                c2 = window_of(p2, l)
                assert c2 is not None
                name = KINDS[c2][0]
                A(f'  · exact (N.hHZn_{name}_{l} g ({side}, t) (fun q hq => hr ((cover_iff _ _).mpr')
                A(f'      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_host hq]; rfl⟩))).symm')
            A(f'  · exact (N.hZn{l} ρ g b (fun q hq hb => hr ((cover_iff _ _).mpr')
            A(f'      ⟨jS (g, q), by show slotR _ _ g q = _; rw [slotR_bank hq, hb]; rfl⟩))).symm')
            A('  · exact absurd ((cover_iff _ _).mpr ⟨jC (g,c), rfl⟩) hr')
            A('')
        elif s.startswith('climb'):
            name = s[5:]; c = KN.index(name); _, p, stw, _, e = KINDS[c]
            k1 = stw[-1] + 1
            costs = []
            for p2, side in ROLES:
                if p2 == p:
                    costs.append(f'(fun _ t φ => bcost φ (N.hrem {c} ({side}, t)))')
                else:
                    costs.append('(fun _ _ _ => 0)')
            A(f'/-- the remaining climb of the host chains {name} (pair {p}): ONE block of rank `hrem {c} j` per role. -/')
            A(f'lemma GHost401.climb{name} (ρ : Fin 8) :')
            A(f'    GRoute (N.S{a} ρ) (N.S{b} ρ) id (fun φ => (Fintype.card Γ : ℝ) * N.hcost {c} φ) := by')
            A(f'  have H1 := GRoute.reach_all (N.S{a} ρ) (N.S{b} ρ)')
            A(f'    (mkR4 {" ".join(costs)} (fun _ _ _ => 0) (fun _ _ _ => 0)) (by')
            A(f'      {RCASES_R}')
            for p2, side in ROLES:
                if p2 != p:
                    A('      · exact GReach.refl _')
                    continue
                if side == 'false':
                    tgt = f'(N.G{e+1} (N.r{e+1} t g)).Xm.Φ (tt (N.I.tv t))'
                    A(f'      · show GReach (N.Hb {c} {k1} g (false, t)) ({tgt}) (fun φ => bcost φ (N.hrem {c} (false, t)))')
                    A(f'        rw [← N.x{e}{e+1} g t]')
                    A(f'        exact (N.hHE_{name} g t).1')
                else:
                    tgt = f'(N.G{e+1} (N.r{e+1} t g)).Xm.Φ 0'
                    A(f'      · show GReach (N.Hb {c} {k1} g (true, t)) ({tgt}) (fun φ => bcost φ (N.hrem {c} (true, t)))')
                    A(f'        rw [← N.y{e}{e+1} g t]')
                    A(f'        exact (N.hHE_{name} g t).2')
            A('      · exact GReach.refl _')
            A('      · exact GReach.refl _)')
            A('  refine H1.cast (fun φ => ?_)')
            A('  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]')
            A('  rw [sum_mkR4]')
            A('  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,')
            A('    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]')
            A('  ring')
            A('')
        elif s.startswith('gate'):
            # copy Shea's gate lemma, with the state index of this position
            m = re.search(r'/-- gate layer [^\n]*\n(lemma G401\.' + s + r' \(ρ : Fin 8\) : XRoute \(N\.S(\d+) ρ\) \(N\.S\d+ ρ\)[\s\S]*?)\n\n', tpl)
            assert m, s
            body = m.group(1)
            body = body.replace('G401.', 'GHost401.').replace(f'N.S{m.group(2)} ρ', f'N.S{a} ρ')
            A(f'/-- gate layer {s[4:]} of run `ρ` (Shea\'s, word for word): twins have equal frames. -/')
            A(body)
            A('')
        elif s == 'fin':
            costs = []
            for p2, side in ROLES:
                if p2 == 4:
                    costs.append('(fun _ _ => N.cX4)' if side == 'false' else '(fun _ _ => N.cY4)')
                else:
                    c = {1: 6, 2: 7, 3: 8}[p2]
                    costs.append(f'(fun _ t φ => bcost φ (N.hrem {c} ({side}, t)))')
            A('/-- the finals of the replica of run `ρ`: the remaining climbs of F1-F3 and pair 4\'s final block. -/')
            A('lemma GHost401.fin (ρ : Fin 8) :')
            A(f'    GRoute (N.S{a} ρ) (N.S{b} ρ) id (fun φ => (Fintype.card Γ : ℝ) *')
            A('      (N.hcost 6 φ + N.hcost 7 φ + N.hcost 8 φ + (Fintype.card T : ℝ) * (N.cX4 φ + N.cY4 φ))) := by')
            A(f'  have H1 := GRoute.reach_all (N.S{a} ρ) (N.S{b} ρ)')
            A(f'    (mkR4 {" ".join(costs)} (fun _ _ _ => 0) (fun _ _ _ => 0)) (by')
            A(f'      {RCASES_R}')
            for p2, side in ROLES:
                comp = '1' if side == 'false' else '2'
                if p2 == 4:
                    A(f'      · exact GReach.ofX (N.f{"X" if side == "false" else "Y"}4 g t)')
                else:
                    A(f'      · exact (N.hHE_F{p2} g t).{comp}')
            A('      · exact GReach.refl _')
            A('      · exact GReach.refl _)')
            A('  refine H1.cast (fun φ => ?_)')
            A('  simp_rw [mkR4_map (fun c : (ℕ → ℝ) → ℝ => c φ)]')
            A('  rw [sum_mkR4]')
            A('  simp only [GHost401.hcost, sum_const, card_univ, nsmul_eq_mul, mul_zero, add_zero, zero_add,')
            A('    Finset.sum_add_distrib, Finset.sum_const_zero, smul_zero]')
            A('  ring')
            A('')
    A('''end Net
end
end XHost
end PowerSaving
end OAI
''')
    for l in range(9):
        A(f'#print axioms OAI.PowerSaving.XHost.GHost401.stage{l}')
    A('#print axioms OAI.PowerSaving.XHost.GHost401.fin')
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------- HCert.lean
def emit_cert(trans):
    tpl = open(os.path.join(TPL, 'GNetCert.lean')).read()
    body = tpl[tpl.index('set_option linter.unusedSectionVars false'):]
    body = body.replace('namespace GX\n', 'namespace XHost\n')
    body = body.replace('open Binary Matrix Finset RAM SS CB RF BR GF\n', 'open Binary Matrix Finset RAM SS CB RF BR GF GX\n')
    body = body.replace('end GX\n', 'end XHost\n')
    body = body.replace('(N : G401 H T Sl C α Γ nb)', '(N : GHost401 H T Sl C α Γ nb)')
    body = body.replace('G401.', 'GHost401.')
    body = body.replace('N.use', 'N.bk')
    body = body.replace('OAI.PowerSaving.GX.', 'OAI.PowerSaving.XHost.')
    for name in ('A2', 'B2', 'A3', 'B3', 'A4', 'B4'):
        body = body.replace(f'GRoute.ofX (N.climb{name} ρ)', f'(N.climb{name} ρ)')
    body = body.replace('GRoute.ofX (N.fin ρ)', '(N.fin ρ)')
    body = body.replace('theorem g401_certificate', 'theorem ghost_certificate')
    body = body.replace('GHost401.g401_certificate', 'ghost_certificate').replace('XHost.g401_certificate', 'XHost.ghost_certificate')
    head = '''import Work.XHost.HStage

/-!
# (key: bridge-host) The certificate of the HOSTED nine-stage word

Copy of Work/Bridge401/GNetCert.lean (Ryan Shea's complex588 tree; written by Nil Mamano's gen401.py) for `GHost401`:
the scalar map of a run is Shea's `chainQ` of the bank map `bk` (a hosted stage has the scalar map of a banked one), so
the sweep, the exchange and the scalar identity are Shea's, word for word.  Written by h-L3-build/gen/gen_hnet.py,
prepared with assistance from Claude.  No `sorry`.
-/

'''
    return head + body


if __name__ == '__main__':
    os.makedirs(os.path.join(OUT, 'Work/XHost'), exist_ok=True)
    out = {}
    out['HRoles'] = emit_roles()
    d, st, trans = emit_def()
    out['HDef'] = d
    out['HStage'] = emit_stage(st, trans)
    out['HCert'] = emit_cert(trans)
    for k, v in out.items():
        p = os.path.join(OUT, 'Work/XHost', k + '.lean')
        open(p, 'w').write(v)
        print(p, v.count('\n'), 'lines')
