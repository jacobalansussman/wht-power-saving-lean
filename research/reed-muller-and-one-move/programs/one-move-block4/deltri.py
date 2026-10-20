#!/usr/bin/env python3
"""deltri.py (t2-block4): the exact DELETION minimum for multi-target slice programs at block 4.
For a target multiset {T_1..T_W}, deletion says: for every form Q, (moves with Q(z) = 1) >= (arrays with <Q,T> = 1).
Minimise the number of moves sum_x c_x over integer counts c_x >= 0 (an ILP, HiGHS).  Every solution is re-checked
exactly in integers.  Runs over all classes {I, T2, ..., TW} up to the orthogonal group (W = 2 or 3).
Usage: deltri.py W"""
import sys, itertools, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
W = int(sys.argv[1]); m = 4
idx = [(i, j) for i in range(m) for j in range(i, m)]; NF = 1 << len(idx)
def Qval(Q, x):
    v = 0
    for b, (i, j) in enumerate(idx):
        if (Q >> b) & 1 and (x >> i) & 1 and (x >> j) & 1: v ^= 1
    return v
def mm(A, B):
    out = []
    for i in range(m):
        r = 0
        for k in range(m):
            if (A[i] >> k) & 1: r ^= B[k]
        out.append(r)
    return tuple(out)
def tr(A): return tuple(sum(((A[r] >> c) & 1) << r for r in range(m)) for c in range(m))
def rank(A):
    rows = list(A); r = 0
    for c in range(m):
        p = next((i for i in range(r, m) if (rows[i] >> c) & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(m):
            if i != r and (rows[i] >> c) & 1: rows[i] ^= rows[r]
        r += 1
    return r
def pairQT(Q, T):
    v = 0
    for b, (i, j) in enumerate(idx):
        if (Q >> b) & 1 and (T[i] >> j) & 1: v ^= 1
    return v
I = tuple(1 << i for i in range(m))
GL = [g for g in itertools.product(range(1, 1 << m), repeat=m) if rank(g) == m]
ALLT = sorted({mm(g, tr(g)) for g in GL}); O = [g for g in GL if mm(g, tr(g)) == I]
QV = np.array([[Qval(Q, x) for x in range(1, 16)] for Q in range(1, NF)])
PT = {T: np.array([pairQT(Q, T) for Q in range(1, NF)]) for T in ALLT}
seen = set(); classes = []
for rest in itertools.combinations_with_replacement(ALLT, W - 1):
    key = tuple(sorted(rest))
    if key in seen: continue
    orb = {tuple(sorted(mm(mm(g, t), tr(g)) for t in rest)) for g in O}
    seen |= orb; classes.append((key, len(orb)))
print('W=%d: %d classes of target sets {I, ...} up to the orthogonal group' % (W, len(classes)), flush=True)
from collections import Counter
dist = Counter(); low = []
for key, sz in classes:
    odd = PT[I] + sum(PT[t] for t in key)
    res = milp(c=np.ones(15), constraints=[LinearConstraint(QV, lb=odd, ub=np.inf)], integrality=np.ones(15), bounds=Bounds(0, np.inf))
    c = np.round(res.x).astype(int); M = int(c.sum())
    assert np.all(QV @ c >= odd)                       # exact integer re-check of the solution
    dist[M] += 1
    if M <= 4 * W - 1: low.append((M, key, sz, c.tolist()))
print('deletion minimum M over classes:', sorted(dist.items()))
for M, key, sz, c in sorted(low)[:40]:
    print('M_del %d  T = I,%s  class size %d  counts(dir 1..15) %s' % (M, ','.join('.'.join(map(str, t)) for t in key), sz, c))
print('classes with deletion minimum <= 4W - 1 = %d: %d of %d' % (4 * W - 1, len(low), len(classes)))
