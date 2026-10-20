#!/usr/bin/env python3
"""fc-reopt: integer global pattern mix nearest (L1) to a Hall-LP class solution (build_fc.py SOL-*.json), with the exact
bin counts per capacity, the exact role inventory per width and at least 224 T empty rank-48 bins (the forced heads);
written in fill17's --mix-file format {"widths": [...], "target": [[cap, counts, n], ...]}.
usage: round_mix.py SOL.json INV.json OUT.json [--extra-rc]"""
import sys, json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
sol = json.load(open(sys.argv[1])); inv = json.load(open(sys.argv[2])); out = sys.argv[3]
W = sol['widths']; T = inv['T']
need = {int(k): v for k, v in inv['need'].items()}
bins = {int(k): v for k, v in inv['bins_by_cap'].items()}
lp = {}
for key, v in sol['global_mix'].items():
    c, q, h = key.split(':'); q = tuple(int(x) for x in q.split(','))
    lp[(int(c), q)] = lp.get((int(c), q), 0.0) + v
keys = sorted(lp); n = len(keys)
empty48 = (48, tuple(0 for _ in W))
assert empty48 in lp
# variables: n_p (n), d+_p (n), d-_p (n)
nv = 3 * n
cost = np.r_[np.zeros(n), np.ones(n), np.ones(n)]
rows, lo, hi = [], [], []
def row(coef, l, u): rows.append(coef); lo.append(l); hi.append(u)
for c in sorted(bins):
    r = np.zeros(nv)
    for i, (cc, q) in enumerate(keys):
        if cc == c: r[i] = 1
    row(r, bins[c], bins[c])
for j, w in enumerate(W):
    if need.get(w, 0) == 0 and not any(q[j] for c, q in keys): continue
    r = np.zeros(nv)
    for i, (c, q) in enumerate(keys): r[i] = q[j]
    row(r, need.get(w, 0), need.get(w, 0))
for i in range(n):           # n_p - lp_p = d+ - d-
    r = np.zeros(nv); r[i] = 1; r[n + i] = -1; r[2 * n + i] = 1
    row(r, lp[keys[i]], lp[keys[i]])
lb = np.zeros(nv); ub = np.full(nv, np.inf); lb[keys.index(empty48)] = 224 * T
integ = np.r_[np.ones(n), np.zeros(2 * n)]
res = milp(cost, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=integ, bounds=Bounds(lb, ub),
           options={'time_limit': 300})
assert res.status == 0, res.message
x = np.rint(res.x[:n]).astype(int)
target = [[c, list(q), int(v)] for (c, q), v in zip(keys, x) if v]
json.dump(dict(widths=W, target=target, source=sys.argv[1], l1_distance=float(res.fun)), open(out, 'w'), indent=1)
print('rounded mix: %d patterns, L1 distance to the LP %.3f' % (len(target), res.fun))
for t in target: print('  ', t)
