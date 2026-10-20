# EXACT certificates for the pi classes excluded by the cut LP (LP > 11): integer dual vectors Y (= y * 1e9, rounded
# down) with  sum_X Y_X [i in X, j notin X] <= 1e9 for every arc (i,j)   and   sum_X Y_X r_pi(X) > 11e9.
# Then every arc multiset n meeting the cut constraints has  1e9 * |n| >= sum_X Y_X r(X) > 11e9, so |n| >= 12.
import sys, itertools
sys.path.insert(0, '.')
import numpy as np, cutip
from scipy.optimize import linprog
N = 8
bitperms = list(itertools.permutations(range(3)))
def act(g, S): return sum(((S >> k) & 1) << g[k] for k in range(3))
seen = set(); reps = []
for pi in itertools.permutations(range(N)):
    if pi in seen: continue
    orb = set()
    for g in bitperms:
        gi = [act(g, S) for S in range(N)]; ginv = [0]*N
        for S in range(N): ginv[gi[S]] = S
        orb.add(tuple(gi[pi[ginv[c]]] for c in range(N)))
    seen |= orb; reps.append(pi)
Ai = cutip.A.astype(np.int64)
certified = notcert = keep = 0
for pi in reps:
    r = cutip.rvec(pi).astype(np.int64)
    res = linprog(np.ones(len(cutip.arcs)), A_ub=-cutip.A, b_ub=-r, bounds=(0, None), method="highs")
    if res.fun <= 11 + 1e-7:
        keep += 1; continue
    y = np.maximum(-res.ineqlin.marginals, 0)
    Y = np.floor(y * 1e9).astype(np.int64)
    col = Y @ Ai
    if col.max() > 10**9:
        Y = (Y * 10**9) // col.max()
        col = Y @ Ai
    ok = col.max() <= 10**9 and int(Y @ r) > 11 * 10**9
    if ok: certified += 1
    else: notcert += 1; print("NOT CERTIFIED", pi, res.fun, int(Y @ r))
print(f"classes {len(reps)}; kept for search (LP <= 11) {keep}; excluded with an exact integer certificate {certified}; "
      f"excluded but certificate failed {notcert}")
