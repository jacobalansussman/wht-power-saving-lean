#!/usr/bin/env python3
"""ctr-search / cover3.py: exact covering checks behind Problem D at m = 10 (see coflat.py).
(B) can three label subspaces of dim <= 6 have port sets covering a co-flat (Gamma_a, 53 ports, or N_i, 56)?
(A) can one port + a dim-6 + a dim-7 subspace cover Gamma_a (the only co-flat reachable by the 1+20+32 bound)?
Method: list every subspace W of the given dim with |S_W cap Gamma| large enough; for each pair, the uncovered rest R
must lie in ONE more subspace of the allowed dim, i.e. dim span U(R) <= that dim (exact F2 rank)."""
import sys, os, itertools, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
from nd import rref_subspaces, mask

def rank_of_mask(m):
    vs = []
    while m:
        b = (m & -m).bit_length() - 1; m &= m - 1; vs.append(LABS[b])
    return rank2_rows(vs)

def subspace_hits(dim, g, need):
    out = {}
    for rows in rref_subspaces(9, 9 - dim):
        m = (1 << N) - 1
        for r in rows: m &= mask[r]
        h = m & g
        if bin(h).count('1') >= need: out[h] = m
    return out          # distinct intersections -> one representative full mask

res = {}
t0 = time.time()
Gam_a = (1 << 0) | sum(1 << b for b in range(N) if JUNK[0][b])
N_0 = sum(1 << b for b in range(N) if LINES[b][0] != 0)
for name, g in (('Gamma_0', Gam_a), ('N_0', N_0)):
    size = bin(g).count('1')
    need6 = size - 40                     # the other two dim-6 sets cover <= 40
    H = subspace_hits(6, g, need6)
    keys = list(H.keys())
    found = 0; pairs = 0
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            u = keys[i] | keys[j]
            if bin(u).count('1') < size - 20: continue
            pairs += 1
            rest = g & ~u
            if rank_of_mask(rest) <= 6: found += 1
    res[f'B {name}'] = dict(size=size, dim6_sets_with_ge_need=len(keys), pairs_checked=pairs, covering_triples=found)
    print(f'(B) {name} ({size} ports): {len(keys)} distinct dim-6 intersections with >= {need6}; pairs checked {pairs}; '
          f'covering triples found {found} ({time.time()-t0:.0f} s)', flush=True)
# (A) one port + dim 6 + dim 7 covering Gamma_0 exactly: dim-7 must hit 32 (max), dim-6 must hit 20 (max), port 1
g = Gam_a
H7 = subspace_hits(7, g, 53 - 20 - 1)
H6 = subspace_hits(6, g, 53 - 32 - 1)
found = 0
for h7 in H7:
    rest = g & ~h7
    # rest must be covered by a dim-6 set plus one port: try removing each port of rest
    m = rest
    while m:
        b = (m & -m).bit_length() - 1; m &= m - 1
        if rank_of_mask(rest & ~(1 << b)) <= 6: found += 1; break
res['A Gamma_0'] = dict(dim7_sets=len(H7), dim6_sets=len(H6), coverings=found)
print(f'(A) Gamma_0: dim-7 intersections hitting >= 32: {len(H7)}; dim-6 hitting >= 20: {len(H6)}; coverings by port+6+7: {found}', flush=True)
json.dump(res, open(os.path.join(HERE, 'cover3.json'), 'w'), indent=1)
