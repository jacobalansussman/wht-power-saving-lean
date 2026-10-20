# c13_overlap.py -- B'-types (AGL(3,2) irreducibles, identified by character) in V56|S8 and in the coset module C^30
# (permutation module on the 30 singular frames, which carries the depth bound). Shared types = coupling in the subgame.
# Same question for the full game: P0-types of V56 vs the 135-frame module.
import numpy as np, itertools
from lib514 import *
import c8_s4 as C8
from libcut import isotypic
R = C8.R; Bp = C8.Bp
rng = np.random.default_rng(1)
sing = np.load('singular_lag.npy'); frames = [frozenset(int(x) for x in row) for row in sing]
fidx = {F: i for i, F in enumerate(frames)}
def perm30(g):
    M = np.zeros((30, 30))
    for i, F in enumerate(frames): M[fidx[frozenset(int(g[v]) for v in F)], i] = 1
    return M
tests = [Bp[j] for j in rng.choice(len(Bp), 40, replace=False)]
def types_of(mats_all, rep_of):
    comps = isotypic(mats_all, rng)
    out = []
    for c in comps:
        Q = c['basis'][:, :c['d']] if False else None
        # character of one irreducible piece: trace over the component / multiplicity
        B = c['basis']; chi = tuple(np.round([np.trace(B.T @ rep_of(t) @ B) / c['m'] for t in tests], 4))
        out.append((c['d'], c['m'], chi))
    return out
t56 = types_of([R.r56(b) for b in Bp], R.r56)
t30 = types_of([perm30(b) for b in Bp], perm30)
print('V56 B\'-types (dim, mult):', [(d, m) for d, m, _ in t56])
print('C^30 B\'-types (dim, mult):', [(d, m) for d, m, _ in t30])
shared = [(d, m1, m2) for d, m1, c1 in t56 for e, m2, c2 in t30 if d == e and c1 == c2]
print('shared B\'-types (dim, mult in V56, mult in C^30):', shared)
