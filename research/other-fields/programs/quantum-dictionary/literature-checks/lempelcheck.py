#!/usr/bin/env python3
"""Minimal N with Gamma = sum_{j<=N} z_j z_j^T over F_2 (symmetric m x m), vs rank(Gamma).
This is the number of COMMUTING directional kernels C_z (pi/4 X-type Pauli rotations) whose product is
H^m D_Gamma H^m (D_Gamma = i^{x^T Gamma x} up to Paulis), against its Clifford rank = rank(Gamma).
Expected (Lempel 1975 / Seroussi-Lempel 1980, not read first-hand): N = rank if Gamma has a nonzero diagonal entry,
rank+1 if Gamma is alternating and nonzero."""
from collections import deque
import itertools
def run(m):
    idx = [(i, j) for i in range(m) for j in range(i, m)]
    pos = {p: k for k, p in enumerate(idx)}
    def outer(z):
        v = 0
        for i in range(m):
            for j in range(i, m):
                if (z >> i) & 1 and (z >> j) & 1:
                    v |= 1 << pos[(i, j)]
        return v
    gens = [outer(z) for z in range(1, 1 << m)]
    dist = {0: 0}; dq = deque([0])
    while dq:
        g = dq.popleft()
        for h in gens:
            n = g ^ h
            if n not in dist:
                dist[n] = dist[g] + 1; dq.append(n)
    def rank_alt(g):
        rows = []
        for i in range(m):
            r = 0
            for j in range(m):
                a, b = min(i, j), max(i, j)
                if (g >> pos[(a, b)]) & 1: r |= 1 << j
            rows.append(r)
        basis = []
        for v in rows:
            for b in basis: v = min(v, v ^ b)
            if v: basis.append(v)
        alt = all(not ((g >> pos[(i, i)]) & 1) for i in range(m))
        return len(basis), alt
    bad = 0; hist = {}
    for g in range(1 << len(idx)):
        r, alt = rank_alt(g)
        pred = r + (1 if (alt and g) else 0)
        key = (r, alt and g != 0)
        hist[key] = hist.get(key, 0) + 1
        if dist[g] != pred: bad += 1
    print(f"m={m}: {1<<len(idx)} symmetric matrices; mismatches vs (rank, +1 if alternating nonzero): {bad}; classes {dict(sorted(hist.items()))}")
for m in (1, 2, 3, 4): run(m)
