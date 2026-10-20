#!/usr/bin/env python3
"""fc-reopt: integer global pattern mixes on the line between two mixes (design A's published-style mix and the
full-cut Hall-LP mix of #2979's words): for each lambda, the integer mix nearest (L1) to (1 - lambda) A + lambda B with
the exact bin counts per capacity, the exact role inventory per width and at least 224 T empty rank-48 bins; also its
full-cut price. Written in fill17's --mix-file format.  usage: mixline.py A_MIX.json B_MIX.json INV.json OUTPREFIX L1 L2 ..."""
import sys, json, os
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
A = json.load(open(sys.argv[1])); B = json.load(open(sys.argv[2])); inv = json.load(open(sys.argv[3])); pre = sys.argv[4]
lams = [float(x) for x in sys.argv[5:]]
W = A['widths']; assert W == B['widths']; T = inv['T']
need = {int(k): v for k, v in inv['need'].items()}; bins = {int(k): v for k, v in inv['bins_by_cap'].items()}
local = {int(k): v for k, v in inv['local'].items()}
ma = {(c, tuple(q)): n for c, q, n in A['target']}; mb = {(c, tuple(q)): n for c, q, n in B['target']}
keys = sorted(set(ma) | set(mb)); n = len(keys); empty48 = (48, tuple(0 for _ in W))


def nearest(target):
    nv = 3 * n; cost = np.r_[np.zeros(n), np.ones(2 * n)]
    rows, lo, hi = [], [], []
    for c in sorted(bins):
        r = np.zeros(nv); r[[i for i, (cc, q) in enumerate(keys) if cc == c]] = 1; rows.append(r); lo.append(bins[c]); hi.append(bins[c])
    for j, w in enumerate(W):
        r = np.zeros(nv)
        for i, (c, q) in enumerate(keys): r[i] = q[j]
        rows.append(r); lo.append(need.get(w, 0)); hi.append(need.get(w, 0))
    for i in range(n):
        r = np.zeros(nv); r[i] = 1; r[n + i] = -1; r[2 * n + i] = 1; rows.append(r); lo.append(target[i]); hi.append(target[i])
    lb = np.zeros(nv); lb[keys.index(empty48)] = 224 * T
    res = milp(cost, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=np.r_[np.ones(n), np.zeros(2 * n)],
               bounds=Bounds(lb, np.full(nv, np.inf)), options={'time_limit': 120})
    assert res.status == 0, res.message
    return np.rint(res.x[:n]).astype(int), res.fun


from price_designs import design
for lam in lams:
    tgt = [(1 - lam) * ma.get(k, 0) + lam * mb.get(k, 0) for k in keys]
    x, d = nearest(tgt)
    idle = {}
    for (c, q), v in zip(keys, x):
        l = c - sum(a * w for a, w in zip(q, W))
        if l and v: idle[l] = idle.get(l, 0) + int(v)
    prof = dict(local)                      # inv local counts are totals over the 9 T calls
    for r, v in idle.items(): prof[r] = prof.get(r, 0) + v
    o = design('lambda %.3f' % lam, prof, T, 224 * T)
    out = '%s-%.3f.json' % (pre, lam)
    json.dump(dict(widths=W, target=[[c, list(q), int(v)] for (c, q), v in zip(keys, x) if v], lam=lam, idle=idle,
                   full_cut=o['full_cut']['b'], uniform=o['uniform']['b']), open(out, 'w'), indent=1)
    print('lambda %.3f  L1 %.1f  idle %s  full cut %s  hub b %s  -> %s' % (lam, d, dict(sorted(idle.items())), o['full_cut']['b'], o['uniform']['b'], out), flush=True)
