# piece_ilp.py: fewest rank-cover pieces for one even support F at block m (uniform family), by an exact-cover ILP.
# Pieces: rectangles R x C (R, C nonempty sets of directions of F) with R u C inside F ∩ G for an odd minimum-weight
# support G (slack s for the uniform family). If the minimum is k, the uniform family needs |F| r <= k s, i.e.
# W <= r 2^(m-2) - |F| r / k.  Flats F give k = 10 (Corollary 15); this checks whether other even supports do better.
# usage: piece_ilp.py m weight howmany
import sys, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from refl import quad_forms
m, wF, howmany = map(int, sys.argv[1:4])
vals, tr = quad_forms(m)
V = vals[:, 1:]; w = V.sum(1)
oddmin = [frozenset(np.nonzero(V[q])[0].tolist()) for q in np.nonzero((tr == 1) & (w == 2 ** (m - 2)))[0]]
evs = [frozenset(np.nonzero(V[q])[0].tolist()) for q in np.nonzero((tr == 0) & (w == wF))[0]][:howmany]
for F in evs:
    P = sorted(F); idx = {x: i for i, x in enumerate(P)}; nP = len(P)
    inter = {F & G for G in oddmin if len(F & G) >= 1}
    rects = set()
    for I in inter:
        L = sorted(I)
        subs = [frozenset(c) for k in range(1, len(L) + 1) for c in itertools.combinations(L, k)]
        for R in subs:
            for C in subs:
                rects.add((R, C))
    rects = list(rects)
    A = np.zeros((nP * nP, len(rects)))
    for j, (R, C) in enumerate(rects):
        for a in R:
            for b in C:
                A[idx[a] * nP + idx[b], j] = 1
    res = milp(c=np.ones(len(rects)), constraints=LinearConstraint(A, 1, 1), integrality=np.ones(len(rects)),
               bounds=Bounds(0, 1), options={"time_limit": 120})
    k = int(round(res.fun)) if res.status in (0, 1) and res.fun is not None else None
    print(f"m={m} even support of weight {wF}: intersections with odd min-weight supports {sorted(len(I) for I in inter)[-3:]}; "
          f"candidate rectangles {len(rects)}; fewest pieces {k} (status {res.status}); uniform bound W <= "
          f"{(2 ** (m - 2) - wF / k) if k else None} r (flats: {2 ** (m - 2) - 2 ** (m - 2) / 10} r)", flush=True)
