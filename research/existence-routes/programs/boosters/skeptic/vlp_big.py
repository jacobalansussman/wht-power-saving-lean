# v-e3-boost vlp_big.py: boosters with ALL jobs of rank <= rmax on a LARGE register M, against my families at block M
# (A, C_k for all k, distance sets with S = intervals, parity sets and all S of size <= 3, lines with all T).
import sys
import numpy as np
from scipy.optimize import linprog
from itertools import combinations
import vlines, vbig
def G(n, k):
    if k < 0 or k > n: return 0
    num = den = 1
    for i in range(k): num *= 2 ** (n - i) - 1; den *= 2 ** (i + 1) - 1
    return num // den
def rows(M, rmax):
    ds = list(range(1, rmax + 1))
    R = [[1] * rmax, [2 ** M - 2 ** (M - d) for d in ds]]
    for k in range(1, M + 1): R.append([G(M, k) - G(M - d, k) for d in ds])
    Ss = set()
    for i in range(M + 1):
        for j in range(i, M + 1): Ss.add(tuple(range(i, j + 1)))
    Ss.add(tuple(range(0, M + 1, 2))); Ss.add(tuple(range(1, M + 1, 2)))
    for r in (1, 2, 3):
        for S in combinations(range(M + 1), r): Ss.add(S)
    for S in Ss:
        P = vbig.dprof(M, S, ds); R.append([P[d] for d in ds])
    J = [vlines.counts(M, d) for d in ds]
    for r in range(1, 6):
        for T in combinations(vlines.TYPES, r):
            R.append([sum(c for (tu, tv), c in Jd.items() if tv in T) - sum(c for (tu, tv), c in Jd.items() if tv in T and tu in T) for Jd in J])
    R = [r for r in R if max(abs(x) for x in r) > 0]
    R = np.array([[float(x) / max(abs(y) for y in r) for x in r] for r in R])
    return R
def boost_min(m, t, R):
    phi = np.array([d ** t for d in range(1, m + 1)]); one = np.arange(1, m + 1, dtype=float)
    c = np.concatenate([phi, -phi])
    A_ub = np.vstack([np.hstack([-R, R]), np.concatenate([-one, one])]); b_ub = np.zeros(len(R) + 1)
    A_eq = np.concatenate([np.zeros(m), np.ones(m)])[None, :]
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=[(0, None)] * (2 * m), method="highs")
    return r.fun, r.x[:m] - r.x[m:]
if __name__ == "__main__":
    grid = [i / 50 for i in range(1, 50)] + [0.9983]
    for a in sys.argv[1:]:
        M, rmax = map(int, a.split(":"))
        R = rows(M, rmax)
        neg = []; worst = (0, None, None)
        for t in grid:
            f, v = boost_min(rmax, t, R)
            if f < -1e-9:
                neg.append(t)
                if f < worst[0]: worst = (f, t, v)
        print(f"my families at block {M}, ranks <= {rmax}: {len(R)} rows; booster-negative at {len(neg)}/{len(grid)}"
              + (f" ({min(neg)}..{max(neg)}); worst {worst[0]:.3g} at {worst[1]} v={np.round(worst[2]/np.abs(worst[2]).max(),4).tolist()}" if neg else ""), flush=True)
