# the one-distribution LP of the deletion test in the commutative slice: min sum p_z s.t. sum_z p_z Q(z) >= 1 for every
# odd Q.  Value (moves per array) and whether the optimum is unique (range of each p_z over the optimal face).
import sys, numpy as np
from scipy.optimize import linprog
from refl import quad_forms
for m in map(int, sys.argv[1:]):
    vals, tr = quad_forms(m)
    A = vals[tr == 1][:, 1:].astype(float)
    n = A.shape[1]
    res = linprog(np.ones(n), A_ub=-A, b_ub=-np.ones(len(A)), bounds=[(0, None)] * n, method="highs")
    v = res.fun
    lo, hi = [], []
    for z in range(n):
        c = np.zeros(n); c[z] = 1
        Aub = np.vstack([-A, np.ones((1, n))]); bub = np.concatenate([-np.ones(len(A)), [v + 1e-9]])
        r1 = linprog(c, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
        r2 = linprog(-c, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
        lo.append(r1.fun); hi.append(-r2.fun)
    print(f"m={m}: LP value {v:.6f} (uniform would give {(2**m-1)/2**(m-2):.6f}); p_z range over the optimal face: "
          f"min of lows {min(lo):.6f}, max of highs {max(hi):.6f}, uniform value {1/2**(m-2):.6f} -> "
          f"{'UNIQUE (uniform)' if max(hi) - min(lo) < 1e-6 else 'not unique'}", flush=True)
