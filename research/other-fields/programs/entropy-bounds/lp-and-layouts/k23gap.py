# k23gap.py -- which small demand sets on K_{2,3} have cut LP < routing (min-cost form)? (calibration choice)
import itertools
import numpy as np
from scipy.optimize import linprog
L, R = ['a', 'b'], ['x', 'y', 'z']
nodes = L + R
arcs = []
for u in L:
    for v in R: arcs += [(u, v), (v, u)]
def dist(s, t):
    if s == t: return 0
    return 1 if ((s in L) != (t in L)) else 2
pairs = [(s, t) for s in nodes for t in nodes if s != t and dist(s, t) == 2]
def cutlp(dem):
    rows = []; rhs = []
    for q in range(1, 1 << 5):
        Q = {nodes[t] for t in range(5) if q >> t & 1}
        need = sum(1 for (s, t) in dem if s not in Q and t in Q)
        if need: rows.append([-1.0 if (u not in Q and v in Q) else 0.0 for (u, v) in arcs]); rhs.append(-need)
    return linprog(np.ones(len(arcs)), A_ub=np.array(rows), b_ub=np.array(rhs), bounds=(0, None), method='highs').fun
for k in (2, 3, 4):
    found = []
    for dem in itertools.combinations(pairs, k):
        c = cutlp(dem); r = 2 * k
        if c < r - 1e-9: found.append((c, dem))
    print(f"{k} sessions (all at distance 2): {len(found)} sets with cut LP < routing; examples:", found[:3])
