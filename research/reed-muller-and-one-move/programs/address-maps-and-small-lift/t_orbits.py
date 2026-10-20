"""t_orbits.py: completeness check of t_small.py's enumeration: the GL(m,2)-orbits of M-multisets of nonzero directions
must have sizes summing to C(2^m - 1 + M - 1, M). Same orbit code as t_small.py."""
import sys, itertools, math
import numpy as np
from t_free import gl_elements
m, M = int(sys.argv[1]), int(sys.argv[2])
nn = 1 << m
G = gl_elements(m)
A = np.zeros((len(G), nn), dtype=np.int64)
for gi, h in enumerate(G):
    for z in range(nn):
        r = 0; i = 0; zz = z
        while zz:
            if zz & 1: r ^= h[i]
            zz >>= 1; i += 1
        A[gi, z] = r
pw = (nn ** np.arange(M)[::-1]).astype(np.int64)
seen = set(); tot = 0; ncl = 0
for ms in itertools.combinations_with_replacement(range(1, nn), M):
    key = 0
    for z in ms: key = key * nn + z
    if key in seen: continue
    keys = set((np.sort(A[:, list(ms)], axis=1) @ pw).tolist())
    seen |= keys; tot += len(keys); ncl += 1
print(f"m={m} M={M}: {ncl} classes, orbit sizes sum to {tot}; C({nn-1+M-1},{M}) = {math.comb(nn - 1 + M - 1, M)}; "
      f"match: {tot == math.comb(nn - 1 + M - 1, M)}")
