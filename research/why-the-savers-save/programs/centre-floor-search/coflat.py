#!/usr/bin/env python3
"""ctr-search / coflat.py: for Problem D's necessary condition NC (for every set K of directions, the ports that can use
only directions in K have root rank <= |K|).  With |K| = 7 the other m - 7 directions' port sets S_W must cover a
CO-FLAT (complement of a hyperplane flat of the 112 lines).  Co-flats of size <= 69 are exactly
Gamma_a = {a} u junk(a) (53 ports, 112 of them) and N_i (56 ports, 8 of them) [E7-type hyperplanes; every other
hyperplane holds <= 42 lines (D7), checked below by enumerating the hyperplanes spanned by roots].
Computes m_d(Gamma) = max over d-dim label subspaces W of |S_W cap Gamma|, d = 1..8, for both co-flat types."""
import sys, os, itertools, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
from nd import rref_subspaces, mask
# co-flats
G4 = G
Gam = {}
for a in range(N):
    Gam[('a', a)] = (1 << a) | sum(1 << b for b in range(N) if JUNK[a][b])
for i in range(8):
    Gam[('i', i)] = sum(1 << b for b in range(N) if LINES[b][i] != 0)
# verify the hyperplane-flat claim on normals = roots/frame (E7 type) and sample others
sizes = {}
for key, m in Gam.items(): sizes[key[0]] = bin(m).count('1')
print('co-flat sizes', sizes)
# max hyperplane flat size over hyperplanes spanned by 7 independent lines: sample + all normals from E8 lattice short vectors
import random
rng = random.Random(1); R = np.array(LINES, dtype=float); best_other = 0
for _ in range(20000):
    S = rng.sample(range(N), 7)
    A = R[S]
    if np.linalg.matrix_rank(A) < 7: continue
    u, s, vt = np.linalg.svd(A); v = vt[-1]
    cnt = int(np.sum(np.abs(R @ v) < 1e-9))
    vv = v / np.max(np.abs(v))
    # is v proportional to a root or frame vector?
    isroot = any(abs(abs(np.dot(v, r)) / np.linalg.norm(r) - 1) < 1e-9 for r in list(R) + [np.eye(8)[i] for i in range(8)])
    if not isroot: best_other = max(best_other, cnt)
print('largest hyperplane flat with a non-root normal (20000 random spanning 7-sets):', best_other)
out = {}
for typ, key in (('Gamma_a', ('a', 0)), ('N_i', ('i', 0))):
    g = Gam[key]; res = {}
    for d in range(1, 9):
        cd = 9 - d; best = 0
        for rows in rref_subspaces(9, cd):
            m = (1 << N) - 1
            for r in rows: m &= mask[r]
            c = bin(m & g).count('1')
            if c > best: best = c
        res[d] = best
    out[typ] = res
    print(typ, 'size', bin(g).count('1'), 'max |S_W cap Gamma| by dim W:', res, flush=True)
json.dump(out, open(os.path.join(HERE, 'coflat.json'), 'w'), indent=1)
