# v-e3-boost vlp.py: my own booster LP on registers m = 3..6 with MY certificate rows only (frame count, A, C_k,
# distance sets with all S, lines with all T, on M = m..m+3 and b = 2), optionally plus the hook values quoted by certlib
# (block-4 hook (4,8), block-5 hooks, block-6 hook k=3; 'whp' mod-p values from wave16/b56-compute/pricecap.py).
import sys
import numpy as np
from scipy.optimize import linprog
from itertools import combinations
import vsurv as S, vlines
HOOKS = {4: [[27456, 48576, 68320, 86592]], 5: [[6960, 10888, 12896, 13888, 14384], [64960, 103488, 124000, 133920, 138880]],
         6: [[60512, 92688, 108864, 116928, 120960, 122976]]}
def rows(m, hooks=False, Mextra=3):
    R = [[1.0] * m]
    for M in range(m, m + Mextra + 1):
        R.append(S.Aprof(M)[:m])
        for k in range(1, M + 1): R.append(S.Cprof(M, k)[:m])
        for r in range(1, M + 2):
            for Sx in combinations(range(M + 1), r): R.append(S.Dprof(M, Sx)[:m])
        for r in range(1, 6):
            for T in combinations(vlines.TYPES, r): R.append(S.Lprof(M, set(T))[:m])
        if hooks:
            for M2, Ps in HOOKS.items():
                if M2 == M:
                    for P in Ps: R.append(P[:m])
    M = 2 * m
    R.append([S.Aprof(M)[2 * d - 1] for d in range(1, m + 1)])
    for k in range(1, M + 1): R.append([S.Cprof(M, k)[2 * d - 1] for d in range(1, m + 1)])
    for i in range(M + 1):
        for j in range(i, M + 1): R.append([S.Dprof(M, tuple(range(i, j + 1)))[2 * d - 1] for d in range(1, m + 1)])
    R = np.array([r for r in R if max(abs(x) for x in r) > 0], dtype=float)
    return R / np.abs(R).max(axis=1, keepdims=True)
def boost_min(m, t, R):
    phi = np.array([d ** t for d in range(1, m + 1)]); one = np.arange(1, m + 1, dtype=float)
    # variables v_plus (calls), v_minus (tasks) >= 0, sum tasks = 1
    c = np.concatenate([phi, -phi])
    A_ub = np.vstack([np.hstack([-R, R]), np.concatenate([-one, one])]); b_ub = np.zeros(len(R) + 1)
    A_eq = np.concatenate([np.zeros(m), np.ones(m)])[None, :]
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=[(0, None)] * (2 * m), method="highs")
    return r.fun, r.x[:m] - r.x[m:]
if __name__ == "__main__":
    grid = [i / 100 for i in range(1, 100)] + [0.9983]
    for m in map(int, sys.argv[1:] or ["3", "4", "5"]):
        for hooks in (False, True):
            R = rows(m, hooks)
            neg = []; worst = (0, None, None)
            for t in grid:
                f, v = boost_min(m, t, R)
                if f < -1e-9:
                    neg.append(t)
                    if f < worst[0]: worst = (f, t, v)
            print(f"m={m} hooks={hooks}: {len(R)} rows; booster-negative at {len(neg)}/{len(grid)} grid points"
                  + (f" ({min(neg)}..{max(neg)}); worst {worst[0]:.3g} at {worst[1]}, v = {np.round(worst[2] / np.abs(worst[2]).max(), 4).tolist()}" if neg else ""), flush=True)
