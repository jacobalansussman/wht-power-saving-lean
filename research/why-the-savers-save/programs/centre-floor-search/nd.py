#!/usr/bin/env python3
"""ctr-search / nd.py: n_d = max number of ports whose labels lie in a d-dim subspace of F_2^9, for d = 5..8,
by enumerating the (9-d)-dim constraint subspaces C (ports in C^perp).  Also the root rank of the best sets."""
import sys, os, itertools, collections, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
mask = [0] * 512
for c in range(512):
    m = 0
    for b in range(N):
        if dot(LABS[b], c) == 0: m |= 1 << b
    mask[c] = m
def rref_subspaces(k, d):
    for piv in itertools.combinations(range(k - 1, -1, -1), d):
        pivset = set(piv); free = []
        for r, p in enumerate(piv): free += [(r, b) for b in range(p) if b not in pivset]
        for mk in range(1 << len(free)):
            rows = [1 << p for p in piv]
            for i, (r, b) in enumerate(free):
                if mk >> i & 1: rows[r] |= 1 << b
            yield rows
out = {}
out = {}
for cd in (tuple(int(x) for x in sys.argv[1:]) if __name__ == "__main__" else ()):          # constraint dim -> label dim 9 - cd
    hist = collections.Counter(); best = (0, None)
    for rows in rref_subspaces(9, cd):
        m = (1 << N) - 1
        for r in rows: m &= mask[r]
        n = bin(m).count('1'); hist[n] += 1
        if n > best[0]: best = (n, m)
    S = [b for b in range(N) if best[1] >> b & 1]
    rr = int(np.linalg.matrix_rank(np.array([LINES[b] for b in S], dtype=float)))
    out[9 - cd] = dict(max=best[0], hist=sorted(hist.items()), root_rank_of_a_best=rr)
    print(f'label dim {9-cd}: max ports {best[0]}; histogram {sorted(hist.items())}; root rank of one best set {rr}', flush=True)
if __name__ == "__main__": json.dump(out, open(os.path.join(HERE, 'nd_' + '_'.join(sys.argv[1:] or ['1','2','3']) + '.json'), 'w'), indent=1)
