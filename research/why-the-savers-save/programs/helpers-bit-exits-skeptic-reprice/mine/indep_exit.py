#!/usr/bin/env python3
"""v-bit-exits-reprice: independent re-derivation of the exit obstruction (no d-bit-exits code imported).
Replays D, A, scatter, B of the word with frames as subspaces (integer root lines, rank over Q via fractions),
values as F2 bitmasks over all 661 columns. For each helper h: holders of e_h after B (ALL registers, targets too),
dim of the Q-span of the union of the holders' frames (J_h). For J_h of dim 7: U = J_h; inside = every register
(targets included) whose end-of-B frame lies in U; test e_h in F2-span of their post-B values.
Also: for every dim-7 frame U in the catalogue or any hyperplane at all, the dirt test only needs J_h <= U."""
import gzip, json, sys
from fractions import Fraction as Q
from collections import Counter
w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2*V + R
L = [tuple(r) for r in w['lines']]
def basis(vecs):
    B = []  # list of (pivot, row) in RREF-ish
    for v in vecs:
        r = [Q(x) for x in v]
        for p, b in B:
            if r[p]: c = r[p]; r = [x - c*y for x, y in zip(r, b)]
        p = next((i for i, x in enumerate(r) if x), None)
        if p is not None:
            c = r[p]; r = [x / c for x in r]
            B = [(pp, [x - bb[p]*y for x, y in zip(bb, r)]) for pp, bb in B]
            B.append((p, r))
    return B
def inspan(B, v):
    r = [Q(x) for x in v]
    for p, b in B:
        if r[p]: c = r[p]; r = [x - c*y for x, y in zip(r, b)]
    return not any(r)
FR = [tuple(f) for f in w['frames']]
fdim = [len(basis([L[i] for i in f])) for f in FR]
fset = [frozenset(f) for f in FR]
cur = list(w['start']); val = [1 << i for i in range(N)]
def gate(f, t, ss):
    for r in [t, *ss]:
        if cur[r] != f:
            assert fset[cur[r]] <= fset[f] and fdim[cur[r]] < fdim[f]
            cur[r] = f
    acc = 0
    for s in ss: acc ^= val[s]
    val[t] ^= acc
for g in w['D']: gate(*g)
for g in w['A']: gate(*g)
cp = [val[q] for q in w['ret']]
for t, row in enumerate(w['scatter']):
    for k in row: val[V+t] ^= cp[k]
for g in w['B']: gate(*g)
endB = list(cur); postB = list(val)
def f2span(vecs, target):
    bs = {}
    for v in vecs:
        x = v
        while x:
            p = x.bit_length()-1
            if p in bs: x ^= bs[p]
            else: bs[p] = x; break
    x = target
    while x:
        p = x.bit_length()-1
        if p in bs: x ^= bs[p]
        else: return False
    return True
jd = Counter(); res = Counter(); detail = []
for h in range(2*V, N):
    holders = [r for r in range(N) if postB[r] >> h & 1]
    lab = set()
    for r in holders: lab |= fset[endB[r]]
    B = basis([L[i] for i in sorted(lab)])
    d = len(B); jd[d] += 1
    if d < 8:
        inside = [r for r in range(N) if all(inspan(B, L[i]) for i in fset[endB[r]])]
        inside_z = [r for r in inside if not (V <= r < 2*V)]
        ok_all = f2span([postB[r] for r in inside], 1 << h)
        ok_z = f2span([postB[r] for r in inside_z], 1 << h)
        res[(d, ok_all, ok_z)] += 1
        detail.append(dict(h=h, J_dim=d, holders=len(holders), holder_kinds=Counter('x' if r < V else 'y' if r < 2*V else 'h' for r in holders),
                           own_endB_dim=fdim[endB[h]], inside=len(inside), inside_targets=len(inside)-len(inside_z), span_all=ok_all, span_z=ok_z))
print('J_h dims:', dict(sorted(jd.items())))
print('below-FULL J_h: (dim, span incl targets, span z only) ->', {str(k): v for k, v in res.items()})
for x in detail: print(x)
# also: helpers whose own end-of-B frame dim, and holders count 1
print('holders==1:', sum(1 for h in range(2*V, N) if sum(1 for r in range(N) if postB[r] >> h & 1) == 1),
      'their own endB dims', Counter(fdim[endB[h]] for h in range(2*V, N) if sum(1 for r in range(N) if postB[r] >> h & 1) == 1))
