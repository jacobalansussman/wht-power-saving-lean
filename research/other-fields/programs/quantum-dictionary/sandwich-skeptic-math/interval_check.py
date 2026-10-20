#!/usr/bin/env python3
"""Independent interval-lemma check: |I(X,Y)| = number of subspaces of F_2^d; all ordered pairs at m=2,3; 41 sources at m=4.
Needs numpy."""
import time, numpy as np
def popc(a): return bin(a).count("1")
def om(u, v, m):
    M = (1 << m) - 1
    return (popc((u & M) & (v >> m)) + popc((v & M) & (u >> m))) & 1
def lag(m):
    level = {frozenset([0])}
    for k in range(m):
        level = {frozenset(S | {x ^ v for x in S}) for S in level for v in range(1, 1 << 2*m) if v not in S and all(om(v, w, m) == 0 for w in S)}
    return sorted(level, key=lambda s: sorted(s))
gal = {0: 1, 1: 2, 2: 5, 3: 16, 4: 67}
for m in (2, 3, 4):
    t = time.time(); Ls = lag(m); n = len(Ls)
    masks = [sum(1 << v for v in L) for L in Ls]
    D = np.array([[m - (popc(a & b).bit_length() - 1) for b in masks] for a in masks], dtype=np.int8)
    srcs = range(n) if m <= 3 else list(range(0, n, 57))
    bad = cnt = 0
    for a in srcs:
        for b in range(n):
            k = int(np.count_nonzero(D[a] + D[:, b] == D[a, b])); cnt += 1
            if k != gal[int(D[a, b])]: bad += 1
    print(f"m={m}: {n} frames, interval size = Galois number on {cnt} ordered pairs: {bad} failures ({time.time()-t:.1f}s)")
