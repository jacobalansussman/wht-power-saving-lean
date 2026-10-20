# ex2-duality lp_combo.py: fluid LP (ex-laser's job model: splits, rule applications at every column, leaves of size 1
# cost 1) with SEVERAL rules at once and jobs larger than the top job allowed (up to cap = 3n), to test the duality's
# prediction that the exponent of a rule set is max of the intersection of the rules' Theta-sets
# (Theta = {theta in [0,1] : sum_paid p r^theta >= sum_targets w k^theta}).
# Rules are HYPOTHETICAL relations (we do not claim they exist); this tests the duality machinery, not our game.
import sys
from math import log
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

def lp(n, rules, cap_factor=3):
    cap = cap_factor * n
    rows, cols, vals, cost = [], [], [], []
    v = 0
    def add(col_entries, c):
        nonlocal v
        for size, coef in col_entries:
            if size >= 2: rows.append(size); cols.append(v); vals.append(coef)
        cost.append(c); v += 1
    for k in range(2, cap + 1):                       # split k -> a + (k - a)
        for a in range(1, k // 2 + 1):
            ent = [(k, 1.0)]; c = 0.0
            for part in (a, k - a):
                if part >= 2: ent.append((part, -1.0))
                else: c += 1.0
            add(ent, c)
    for Wt, P in rules:                                # rule at column s
        for s in range(1, cap + 1):
            if max(list(Wt) + list(P)) * s > cap: break
            ent = []; c = 0.0
            for k, w in Wt.items():
                if k * s >= 2: ent.append((k * s, w))
                else: c -= 0.0  # a size-1 target cannot be consumed by a rule in this LP (would be a free leaf); skip rule
            if any(k * s < 2 for k in Wt): continue
            for r, p in P.items():
                if r * s >= 2: ent.append((r * s, -p))
                else: c += p
            add(ent, c)
    A = coo_matrix((vals, ([r - 2 for r in rows], cols)), shape=(cap - 1, v)).tocsr()
    b = np.zeros(cap - 1); b[n - 2] = 1.0
    res = linprog(np.array(cost), A_eq=A, b_eq=b, bounds=(0, None), method='highs')
    assert res.status == 0, res.message
    return res.fun

def theta_max(rules, grid=10000):
    best = 0.0
    for i in range(grid + 1):
        t = i / grid
        if all(sum(p * r ** t for r, p in P.items()) >= sum(w * k ** t for k, w in Wt.items()) - 1e-12 for Wt, P in rules):
            best = t
    return best

R1 = ({3: 2}, {1: 5})                 # hypothetical 2 in 5 at block 3 (proved impossible in the engine game; a test rule)
R2 = ({2: 2}, {1: 1.05, 3: 1})        # 2[2] <= 1.05[1] + [3]: no saving at plain prices, D < 0 on a middle window
A5 = ({1: 1}, {2: 2})                 # [1] <= 2[2] (lower boundedness) -- needs a size-1 target, skipped by this LP
Ns = [int(x) for x in sys.argv[1:]] or [30, 60, 120, 240]
for name, rules in [("R1 alone (2 in 5 at block 3)", [R1]), ("R2 alone (middle window)", [R2]), ("R1 + R2", [R1, R2])]:
    tm = theta_max(rules)
    costs = [lp(N, rules) for N in Ns]
    print(f"{name}: duality exponent max Theta = {tm:.4f}; LP cost {[round(c, 4) for c in costs]}; local slopes "
          f"{[round(log(costs[i + 1] / costs[i]) / log(Ns[i + 1] / Ns[i]), 4) for i in range(len(Ns) - 1)]}", flush=True)
