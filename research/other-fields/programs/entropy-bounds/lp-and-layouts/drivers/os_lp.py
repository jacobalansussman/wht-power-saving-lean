# os_lp.py -- calibration on the Okamura-Seymour network K_{2,3} (4 commodities, cut condition holds, routing
# needs 8 edge uses; Harvey-Kleinberg-Lehman 2006 showed coding cannot beat routing there).
# Here: min-cost form, sessions x->y, y->z, z->x, a->b at rate 1; full Shannon LP with Z3 orbit variables.
# Also: the cut-only LP on the same network, and Ingleton / DFZ cutting planes at the optimum.
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import itertools, json, sys, time
import numpy as np
import nclp

L, R = ['a', 'b'], ['x', 'y', 'z']
nodes = L + R
edges = [(u, v) for u in L for v in R]
arcs = []
for (u, v) in edges: arcs += [(u, v), (v, u)]
sessions = [('x', ['y'], 1), ('y', ['z'], 1), ('z', ['x'], 1), ('a', ['b'], 1)]
net = nclp.Net(nodes, arcs, sessions)
rot = {'a': 'a', 'b': 'b', 'x': 'y', 'y': 'z', 'z': 'x'}
perm = [None] * net.n
for i, (s, t, r) in enumerate(sessions):
    img = (rot[s], [rot[u] for u in t])
    perm[i] = next(j for j, (s2, t2, r2) in enumerate(sessions) if s2 == img[0] and t2 == img[1])
for a, (u, v) in enumerate(arcs):
    perm[net.ns + a] = net.ns + arcs.index((rot[u], rot[v]))
G = nclp.close_group([tuple(perm)], net.n)
print("group order", len(G)); sys.stdout.flush()

# cut-only LP: min sum c_a s.t. c(in(Q)) >= r(sessions with src outside Q, sink in Q), c(out(Q)) >= r(src in Q, sink out)
from scipy.optimize import linprog
rows = []; rhs = []
for q in range(1, 1 << len(nodes)):
    Q = {nodes[t] for t in range(len(nodes)) if q >> t & 1}
    need = sum(r for (s, t, r) in sessions if s not in Q and any(u in Q for u in t))
    if need > 0:
        rows.append([-1.0 if (u not in Q and v in Q) else 0.0 for (u, v) in arcs]); rhs.append(-need)
res = linprog(np.ones(len(arcs)), A_ub=np.array(rows), b_ub=np.array(rhs), bounds=(0, None), method='highs')
print("cut-only LP value", res.fun, "(routing value 8)"); sys.stdout.flush()

lp = nclp.build_lp(net, perms=G)
res = nclp.solve(lp)
out = {'cut_lp': None, 'shannon': res.fun}
h = res.x
# Ingleton / DFZ cutting planes, letters = all unions of <= 2 ground elements and the node-based sets
letters = sorted({1 << i for i in range(net.n)} | {(1 << i) | (1 << j) for i in range(net.n) for j in range(i + 1, net.n)})
cuts = []
val = res.fun
for rnd in range(4):
    new_i, vi = nclp.find_cuts(h, lp['canon'], net.n, letters, 'ingleton', tries=300000)
    new_d, vd = nclp.find_cuts(h, lp['canon'], net.n, letters, 'dfz', tries=300000)
    print(f"round {rnd}: violated Ingleton {len(new_i)} (max {vi:.2e}), DFZ {len(new_d)} (max {vd:.2e})"); sys.stdout.flush()
    if not new_i and not new_d: break
    cuts += new_i + new_d
    res = nclp.solve(lp, extra_ub=cuts); h = res.x; val = res.fun
out['with_linear_rank_cuts'] = val
# dual: which closure equalities carry weight
try:
    duals = res.eqlin.marginals
    out['n_eq_dual_nonzero'] = int(np.sum(np.abs(duals) > 1e-9))
except Exception as e:
    out['dual_err'] = str(e)
json.dump(out, open('os_result.json', 'w'), indent=1)
print(out)
