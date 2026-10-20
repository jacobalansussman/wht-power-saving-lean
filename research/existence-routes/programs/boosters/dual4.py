# e3-boost dual4.py: for registers of m = 3, 4 rows, write phi_theta = (d^theta)_d as a NONNEGATIVE combination of
# rigorous certificate rows plus the plain balance phi_1, for theta on a fine grid (dual of the booster LP).
# Reports which rows are used, so an exact proof can be built from a few fixed generators.
import sys
import numpy as np
from scipy.optimize import linprog
import certlib
from lp import rows

for m in (3, 4):
    R, names = rows(m, certlib.RIGOROUS)
    G = np.vstack([R, np.arange(1, m + 1, dtype=float)]); gn = names + ["plain"]
    used = {}
    for t in [i / 100 for i in range(1, 100)]:
        phi = np.arange(1, m + 1, dtype=float) ** t
        # min sum lambda s.t. G^T lambda = phi, lambda >= 0  (feasible iff phi in the cone)
        r = linprog(np.ones(len(G)), A_eq=G.T, b_eq=phi, bounds=[(0, None)] * len(G), method="highs")
        if r.status != 0: print(f"m={m} theta={t}: NOT in the cone"); continue
        key = tuple(gn[i] for i in np.nonzero(r.x > 1e-10)[0])
        used.setdefault(key, []).append(t)
    for k, ts in used.items():
        print(f"m={m}: generators {k}: theta in [{min(ts)}, {max(ts)}] ({len(ts)} grid points)")
