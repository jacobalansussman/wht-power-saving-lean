# k23two.py -- the distance-2 skeleton of the Gamma_2 star: K_{2,3} with sessions A->B (8), B->A (8) and X_i->X_j
# (4 each, all 6 ordered pairs). Routing = 2*(8+8+24) = 80; does full Shannon certify 80? (20 ground elements,
# symmetry S3 x (A<->B) of order 12, closure-reduced, PDLP with a rigorous dual bound.)
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import sys, json, time, itertools
import numpy as np
import nclp
L, R = ['A', 'B'], ['X1', 'X2', 'X3']
arcs = []
for u in L:
    for v in R: arcs += [(u, v), (v, u)]
sessions = [('A', ['B'], 8), ('B', ['A'], 8)] + [(a, [b], 4) for a in R for b in R if a != b]
net = nclp.Net(L + R, arcs, sessions)
def perm_of(mp):
    p = [None] * net.n
    for i, (s, t, r) in enumerate(sessions):
        p[i] = next(j for j, (s2, t2, r2) in enumerate(sessions) if s2 == mp[s] and t2 == [mp[u] for u in t])
    for a, (u, v) in enumerate(arcs): p[net.ns + a] = net.ns + arcs.index((mp[u], mp[v]))
    return tuple(p)
gens = [perm_of({'A': 'A', 'B': 'B', 'X1': 'X2', 'X2': 'X3', 'X3': 'X1'}),
        perm_of({'A': 'A', 'B': 'B', 'X1': 'X2', 'X2': 'X1', 'X3': 'X3'}),
        perm_of({'A': 'B', 'B': 'A', 'X1': 'X1', 'X2': 'X2', 'X3': 'X3'})]
G = nclp.close_group(gens, net.n)
print("group order", len(G)); sys.stdout.flush()
# cut LP on this network
from scipy.optimize import linprog
rows = []; rhs = []
nodes = L + R
for q in range(1, 1 << 5):
    Q = {nodes[t] for t in range(5) if q >> t & 1}
    need = sum(r for (s, t, r) in sessions if s not in Q and t[0] in Q)
    if need: rows.append([-1.0 if (u not in Q and v in Q) else 0.0 for (u, v) in arcs]); rhs.append(-need)
cut = linprog(np.ones(len(arcs)), A_ub=np.array(rows), b_ub=np.array(rhs), bounds=(0, None), method='highs').fun
print("cut LP", cut, " routing 80"); sys.stdout.flush()
t0 = time.time()
lp = nclp.build_closed_lp(net, perms=G)
primal, LB, res = nclp.solve_pdlp(lp, sum(r for (s, t, r) in sessions), time_limit=7200)
json.dump(dict(cut=cut, routing=80, shannon_primal=primal, shannon_rigorous_lb=LB), open('k23two.json', 'w'))
print("done", time.time() - t0)
