# grouped_lp.py (meta/ex-laser): does the exact tau-transfer survive GROUPED paid blocks?  Fluid LP as in
# meta/literature/tau_schedule_lp.py, but the unit at column s consumes W[k] demands of size k*s and returns P[r]
# demands of size r*s (r = 1 and r >= 2).  If min cost / n^tau_c stays bounded, the exact transfer holds in the
# fluid model for that rule [computed evidence only].  Also reports the Phi lower bound (cost >= n^tau_c).
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from math import log

def solve(n, Wt, P):
    kmax = max(Wt)
    var = []
    for k in range(2, n + 1):
        for a in range(1, k // 2 + 1): var.append(('s', k, a))
    for s in range(1, n // kmax + 1): var.append(('u', s))
    idx = {v: i for i, v in enumerate(var)}
    A = lil_matrix((n + 1, len(var))); b = np.zeros(n + 1); c = np.zeros(len(var))
    for v, i in idx.items():
        if v[0] == 's':
            _, k, a = v
            A[k, i] += 1
            for part in (a, k - a):
                if part >= 2: A[part, i] -= 1
                else: c[i] += 1
        else:
            s = v[1]
            for k, w in Wt.items(): A[k * s, i] += w
            for r, p in P.items():
                if r * s >= 2: A[r * s, i] -= p
                else: c[i] += p
    b[n] = 1
    res = linprog(c, A_eq=A[2:].tocsr(), b_eq=b[2:], bounds=(0, None), method='highs')
    assert res.status == 0, res.message
    return res.fun

def tau_c(Wt, P):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if sum(w * k ** mid for k, w in Wt.items()) > sum(p * r ** mid for r, p in P.items()): hi = mid
        else: lo = mid
    return lo

if __name__ == "__main__":
    cases = [("targets {3:1, 2:4}, paid {1:10} (control, = 1+4 in 10)", {3: 1, 2: 4}, {1: 10}),
             ("targets {3:1, 2:4}, paid {1:6, 2:2}", {3: 1, 2: 4}, {1: 6, 2: 2}),
             ("targets {3:2, 2:3}, paid {1:4, 2:4}", {3: 2, 2: 3}, {1: 4, 2: 4}),
             ("targets {4:1, 3:2}, paid {1:3, 2:2, 3:1}", {4: 1, 3: 2}, {1: 3, 2: 2, 3: 1})]
    Ns = [60, 120, 240, 480]
    for name, Wt, P in cases:
        tc = tau_c(Wt, P)
        costs = [solve(N, Wt, P) for N in Ns]
        print(f"{name}: tau_c {tc:.4f}; cost/n^tau_c {[round(cv / N ** tc, 4) for N, cv in zip(Ns, costs)]}; "
              f"slopes {[round(log(costs[i+1]/costs[i])/log(2), 4) for i in range(3)]}", flush=True)
