#!/usr/bin/env python3
"""v-cx-mirror / v5perx.py (own code, stdlib only; reads the U8X lowering as data, runs nothing of anyone else's).

Question: is cx-mirror's "first-order cost of building and cleaning" (462 natural / 272 best twist per call, called a
deliberately optimistic model, and "cleaning + gathering >= 272 per call" in its STATUS) a lower bound inside its OWN
pricing?  Its pricing: climbs free, a drop F -> F' costs dim F - dim F', extra (hosted) scratch registers free,
gates only between arrays at one equal frame, data arrays never drop, x_b's frame stays between <u_b> and v_b^perp.

I simulate, at the level of frames AND exact values (Fractions over the 112 x's), two multi-register schemes for the
x side of one forward call (the totals' gather, the U-drop, and the cleanup; y visits of the scatter are left out, as
in the claim's own cost table):

 S1 "per-x copies": r_b copies x_b at <u_b> (free); all copies climb to the window W; R_k = sum_b c_kb r_b at W;
    U-drop of a copy D_k of R_k from W to 0 (charged dim W, the centre); D_k back to W, D_k -= R_k; R_k -= sum c r;
    each r_b drops W -> W & v_b^perp (1 level), x_b climbs there, r_b -= x_b.
 S2 "run registers": x's with proportional coefficient vectors across the 8 totals are grouped; each group is split
    into runs (no N-edge, N[b][a] = g_b.u_a) minimising sum of rank{g_b : b in run}; one register per run gathers at
    span(u_run), the rest as S1 with one drop of rank{g_run} per run register.

Every gate is checked for frame equality, every climb for containment, every x frame for <u_b> <= F <= v_b^perp, and
at the end all scratch values are exactly zero and each R_k held exactly sigma_k at its U-drop.  Bits 0..8 window,
bit 9 = s."""
import gzip, json, os, sys, time
from fractions import Fraction as Fr
from itertools import combinations
T0 = time.time()
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(gzip.open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))))
v, R0 = d['v'], d['R']; n0 = 2 * v + R0
U = list(d['ports']); SB = 1 << 9; J = (1 << 9) - 1
def dot(a, b): return bin(a & b).count('1') & 1
# ---- phase A replay (my own, as in v1.py) to get the 8 totals with exact coefficients
reg = [dict() for _ in range(n0)]
for t in range(v): reg[t] = {t: Fr(1)}
def axpy0(t, src, c):
    dd = reg[t]
    for k, val in src.items():
        w = dd.get(k, 0) + c * val
        if w: dd[k] = w
        else: dd.pop(k, None)
for g in d['A']:
    if g[0] == 'out':
        src = dict(reg[g[2]])
        for t, a, b in g[3]: axpy0(t, src, Fr(a, b))
    else:
        for s_, a, b in g[3]: axpy0(g[2], dict(reg[s_]), Fr(a, b))
tot = [dict(reg[s_]) for _, s_, _ in sorted(d['ret'])]
K = len(tot)
assert K == 8 and all(len(T) == 56 for T in tot)
# ---- subspaces of F_2^10 as canonical reduced echelon tuples
def canon(vecs):
    """reduced row echelon (pivot = highest bit), canonical"""
    rows = []
    for x in vecs:
        for r in rows:
            if (x >> (r.bit_length() - 1)) & 1: x ^= r
        if x:
            p = x.bit_length() - 1
            rows = [r ^ x if (r >> p) & 1 else r for r in rows]
            rows.append(x)
    return tuple(sorted(rows, reverse=True))
def inside(Fa, Fb):  # Fa <= Fb
    for x in Fa:
        y = x
        for b in Fb:
            if (y >> (b.bit_length() - 1)) & 1: y ^= b
        if y: return False
    return True
def meet_perp(F, gs):  # F & {w : w.g = 0 for g in gs}
    vec = list(F)
    for g in gs:
        bad = [f for f in vec if dot(f, g)]
        if bad:
            p = bad[0]
            vec = [f ^ p if dot(f, g) else f for f in vec if f != p]
    return canon(vec)
def perp(gs, ambient):
    return meet_perp(ambient, gs)
FULL = canon([1 << i for i in range(10)]); WIN = canon([1 << i for i in range(9)]); ZERO = ()
class Bad(Exception): pass
class Sim:
    def __init__(s, Vl):
        s.V = Vl; s.cost = 0; s.udrop = 0; s.frame = {}; s.val = {}; s.kind = {}; s.ngates = 0
        for b in range(v):
            s.frame[('x', b)] = canon([U[b]]); s.val[('x', b)] = {b: Fr(1)}; s.kind[('x', b)] = 'x'
        s.ceil = {('x', b): perp([Vl[b]], FULL) for b in range(v)}
    def new(s, r):
        s.frame[r] = ZERO; s.val[r] = {}; s.kind[r] = 'scratch'
    def climb(s, r, F):
        F = canon(F)
        if not inside(s.frame[r], F): raise Bad('climb not nested %r' % (r,))
        if s.kind[r] == 'x' and not inside(F, s.ceil[r]): raise Bad('x above its ceiling %r' % (r,))
        s.frame[r] = F
    def drop(s, r, F, centre=False):
        F = canon(F)
        if s.kind[r] == 'x': raise Bad('data drop')
        if not inside(F, s.frame[r]): raise Bad('drop not nested')
        c = len(s.frame[r]) - len(F)
        if centre: s.udrop += c
        else: s.cost += c
        s.frame[r] = F
    def gate(s, t, srcs):
        F = s.frame[t]
        for r, _ in srcs:
            if s.frame[r] != F: raise Bad('gate frames differ: %r %r' % (t, r))
        dd = s.val[t]
        for r, c in srcs:
            for k, val in list(s.val[r].items()):
                w = dd.get(k, 0) + c * val
                if w: dd[k] = w
                else: dd.pop(k, None)
        s.ngates += 1
    def finish(s):
        for r in list(s.frame):
            if s.kind[r] == 'x':
                s.climb(r, s.ceil[r])
                if s.val[r] != {r[1]: 1}: raise Bad('x changed')
            else:
                if s.val[r]: raise Bad('dirty scratch %r' % (r,))
                s.climb(r, FULL)
def lines(A):
    def ap(x):
        y = 0
        for i in range(9):
            if x >> i & 1: y ^= A[i]
        return y
    return [SB | (J ^ ap(u)) for u in U]
tw = json.load(open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../second-line-mirrored/all-twists/twistall.json'))))
M = tw['best']['M']
cases = {'natural': lines([1 << i for i in range(9)]), 'best_twist': lines([(1 << i) ^ M[i] for i in range(9)])}
def check_lines(Vl):
    ok = all(dot(U[b], Vl[b]) == 0 and dot(Vl[b], Vl[b]) == 1 for b in range(v))
    return ok
def run_groups(Vl, groups):
    """groups: list of (x list, coefficient pattern per x).  Generic engine for S1 (singletons) and S2."""
    S = Sim(Vl)
    Gw = [x & J for x in Vl]                    # window part of v_b (= g_b)
    # gather: one register per group at span(u_group)
    for gi, (xs, a) in enumerate(groups):
        r = ('r', gi); S.new(r)
        F = canon([U[b] for b in xs]); S.climb(r, F)
        for b in xs:
            S.climb(('x', b), F)
            S.gate(r, [(('x', b), a[b])])
        S.climb(r, WIN)
    # totals at W: R_k = sum lambda * r_g ; lambda found from the pattern (group content = sum a_b x_b)
    for k in range(K):
        Rk = ('R', k); S.new(Rk); S.climb(Rk, WIN)
        srcs = []
        for gi, (xs, a) in enumerate(groups):
            lam = None
            for b in xs:
                c = tot[k].get(b, 0)
                l = c / a[b]
                if lam is None: lam = l
                elif l != lam: raise Bad('group not proportional for total %d' % k)
            if lam: srcs.append((('r', gi), lam))
        S.gate(Rk, srcs)
        if S.val[Rk] != tot[k]: raise Bad('R_k != sigma_k')
        # U-drop: copy D_k at W, drop to 0 (centre charge); the scatter would read D_k at frame 0
        Dk = ('D', k); S.new(Dk); S.climb(Dk, WIN); S.gate(Dk, [(Rk, Fr(1))]); S.drop(Dk, ZERO, centre=True)
        S.climb(Dk, WIN); S.gate(Dk, [(Rk, Fr(-1))])
        S.gate(Rk, [(r, -c) for r, c in srcs])
    # cleanup: each group register drops to W & g_group^perp, the x's climb there, register -= its content
    for gi, (xs, a) in enumerate(groups):
        r = ('r', gi)
        F = meet_perp(WIN, [Gw[b] for b in xs])
        S.drop(r, F)
        for b in xs:
            S.climb(('x', b), F)
        S.gate(r, [(('x', b), -a[b]) for b in xs])
    S.finish()
    return S
def rank(vecs): return len(canon(vecs))
res = {}
for name, Vl in cases.items():
    assert check_lines(Vl)
    N = [[dot(Vl[b] & J, U[a]) for a in range(v)] for b in range(v)]
    sym = all(N[a][b] == N[b][a] for a in range(v) for b in range(v))
    # S1
    S = run_groups(Vl, [([b], {b: Fr(1)}) for b in range(v)])
    s1 = dict(cleanup_drops=S.cost, udrop=S.udrop, gates=S.ngates)
    # S2: proportional classes of coefficient vectors
    vecs = {b: tuple(tot[k].get(b, Fr(0)) for k in range(K)) for b in range(v)}
    classes = []
    for b in range(v):
        for cl in classes:
            c0 = cl[0]
            same_supp = all((vecs[b][k] == 0) == (vecs[c0][k] == 0) for k in range(K))
            ratios = {vecs[b][k] / vecs[c0][k] for k in range(K) if vecs[c0][k] != 0}
            if same_supp and len(ratios) == 1:
                cl.append(b); break
        else:
            classes.append([b])
    # best partition of each class into runs minimising sum of rank{g_b}
    Gw = [x & J for x in Vl]
    def indep(T): return all(N[a][b] == 0 for a in T for b in T if a != b)
    def best_part(xs):
        from functools import lru_cache
        xs = tuple(xs)
        @lru_cache(None)
        def f(mask):
            if mask == 0: return (0, ())
            items = [xs[i] for i in range(len(xs)) if mask >> i & 1]
            i0 = (mask & -mask).bit_length() - 1
            best = None
            rest = mask & ~(1 << i0)
            sub = rest
            while True:
                T = [xs[i0]] + [xs[i] for i in range(len(xs)) if sub >> i & 1]
                if indep(T):
                    c, p = f(rest & ~sub)
                    c += rank([Gw[b] for b in T])
                    if best is None or c < best[0]: best = (c, (tuple(T),) + p)
                if sub == 0: break
                sub = (sub - 1) & rest
            return best
        return f((1 << len(xs)) - 1)
    groups = []; part_cost = 0
    for cl in classes:
        c, p = best_part(cl); part_cost += c
        for T in p:
            groups.append((list(T), {b: vecs[b][next(k for k in range(K) if vecs[T[0]][k])] for b in T}))
    S2 = run_groups(Vl, groups)
    res[name] = dict(N_symmetric=sym, S1=s1, classes=sorted(len(c) for c in classes),
                     S2=dict(runs=len(groups), run_sizes=sorted(len(g[0]) for g in groups), cleanup_drops=S2.cost,
                             predicted=part_cost, udrop=S2.udrop, gates=S2.ngates))
    print(name, '| N symmetric', sym, '| S1 per-x copies: cleanup drops', s1['cleanup_drops'], 'U-drops', s1['udrop'],
          'gates', s1['gates'], '| proportional classes', len(classes), 'sizes', sorted({len(c) for c in classes}),
          '| S2 run registers', len(groups), 'cleanup drops', S2.cost, 'U-drops', S2.udrop)
print('every gate frame-equal, every climb nested, x frames inside [<u_b>, v_b^perp], all scratch exactly 0 at the end,'
      ' R_k = sigma_k exactly at each U-drop: True (any failure raises)')
res['secs'] = round(time.time() - T0, 1)
json.dump(res, open(os.path.join(here, 'v5perx.json'), 'w'), indent=1, default=str)
print('secs', res['secs'])
