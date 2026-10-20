# v-ex2-duality mylp.py: my own fluid LP of the job model (written from ex-laser 1.1 and ex2-duality's description,
# not from their code).  Sizes 2..cap are rows; a size-1 job is a leaf of cost 1.  Columns: splits k -> a + (k-a);
# each rule (targets T, paid P) at every column s with all sizes <= cap.  Constraint: supply - demand >= [size == n]
# (surplus may be discarded); option --eq for equality (as in ex2-duality / ex-laser).  Objective: total leaf cost.
import sys
from math import log
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

def lp(n, rules, capf=3, eq=False):
    cap = capf * n
    R, C, V, cost = [], [], [], []
    col = 0
    def push(entries, c):
        nonlocal col
        acc = {}
        for size, coef in entries: acc[size] = acc.get(size, 0.0) + coef
        for size, coef in acc.items():
            if coef != 0: R.append(size - 2); C.append(col); V.append(coef)
        cost.append(c); col += 1
    for k in range(2, cap + 1):
        for a in range(1, k // 2 + 1):
            ent, c = [(k, 1.0)], 0.0
            for part in (a, k - a):
                if part >= 2: ent.append((part, -1.0))
                else: c += 1.0
            push(ent, c)
    for T, P in rules:
        top = max(list(T) + list(P))
        for s in range(1, cap // top + 1):
            if any(k * s < 2 for k in T): continue
            ent, c = [(k * s, float(w)) for k, w in T.items()], 0.0
            for r, p in P.items():
                if r * s >= 2: ent.append((r * s, -float(p)))
                else: c += float(p)
            push(ent, c)
    A = csr_matrix((V, (R, C)), shape=(cap - 1, col))
    b = np.zeros(cap - 1); b[n - 2] = 1.0
    if eq: res = linprog(np.array(cost), A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    else: res = linprog(np.array(cost), A_ub=-A, b_ub=-b, bounds=(0, None), method="highs")
    assert res.status == 0, res.message
    return res.fun

def thmax(rules, grid=20000):
    best = None
    for i in range(grid + 1):
        t = i / grid
        if all(sum(p * r ** t for r, p in P.items()) >= sum(w * k ** t for k, w in T.items()) - 1e-12 for T, P in rules):
            best = t
    return best

R1 = ({3: 2}, {1: 5}); R2 = ({2: 2}, {1: 1.05, 3: 1}); R2b = ({4: 2}, {1: 1, 7: 1}); Rs = ({2: 2}, {1: 1, 3: 1})
GT = ({3: 1, 2: 4}, {1: 6, 2: 2}); GR = ({3: 1, 2: 2}, {1: 6}); GM = ({4: 1, 3: 2}, {1: 5, 2: 2})
cases = {"R1": [R1], "R1+R2": [R1, R2], "R1+R2b": [R1, R2b], "R1+Rs": [R1, Rs], "Rs": [Rs],
         "GT": [GT], "GR": [GR], "GM": [GM]}
which = sys.argv[1]; capf = int(sys.argv[2]); eq = sys.argv[3] == "eq"; Ns = [int(x) for x in sys.argv[4:]]
rules = cases[which]
costs = [lp(N, rules, capf, eq) for N in Ns]
sl = [round(log(costs[i + 1] / costs[i]) / log(Ns[i + 1] / Ns[i]), 4) for i in range(len(Ns) - 1)]
print(f"{which} cap={capf}n {'eq' if eq else 'ge'}: theta_max {thmax(rules):.4f}; n {Ns}; cost {[round(c, 4) for c in costs]}; slopes {sl}", flush=True)
