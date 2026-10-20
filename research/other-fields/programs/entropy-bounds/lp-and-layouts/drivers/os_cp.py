# os_cp.py -- Okamura-Seymour K_{2,3} calibration, full Shannon LP by cutting planes (exact same optimum as the
# full LP), then Ingleton and DFZ cutting planes; also the worst case over all rate vectors of these 4 sessions.
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import json, sys, numpy as np
import nclp
L, R = ['a', 'b'], ['x', 'y', 'z']
nodes = L + R
arcs = []
for u in L:
    for v in R: arcs += [(u, v), (v, u)]
sessions = [('x', ['y'], 1), ('y', ['z'], 1), ('z', ['x'], 1), ('a', ['b'], 1)]
net = nclp.Net(nodes, arcs, sessions)
rot = {'a': 'a', 'b': 'b', 'x': 'y', 'y': 'z', 'z': 'x'}
perm = [None] * net.n
for i, (s, t, r) in enumerate(sessions):
    perm[i] = next(j for j, (s2, t2, r2) in enumerate(sessions) if s2 == rot[s] and t2 == [rot[u] for u in t])
for a, (u, v) in enumerate(arcs): perm[net.ns + a] = net.ns + arcs.index((rot[u], rot[v]))
G = nclp.close_group([tuple(perm)], net.n)
out = {}
res, canon, key, active = nclp.solve_cp(net, perms=G, batch=60000)
out['shannon_rate1'] = res.fun
h = res.x
print("Shannon LP (rates 1):", res.fun); sys.stdout.flush()
letters = sorted({1 << i for i in range(net.n)} | {(1 << i) | (1 << j) for i in range(net.n) for j in range(i + 1, net.n)})
cuts = []
for rnd in range(3):
    ni, vi = nclp.find_cuts(h, canon, net.n, letters, 'ingleton', tries=400000)
    nd, vd = nclp.find_cuts(h, canon, net.n, letters, 'dfz', tries=400000)
    print(f"  linear-rank round {rnd}: Ingleton violated {len(ni)} (max {vi:.2e}), DFZ violated {len(nd)} (max {vd:.2e})"); sys.stdout.flush()
    if not ni and not nd: break
    cuts += ni + nd
    res, canon, key, active = nclp.solve_cp(net, perms=G, batch=60000, extra_ub=cuts)
    h = res.x; print("  value with linear-rank cuts:", res.fun); sys.stdout.flush()
out['with_linear_rank_cuts'] = res.fun
# worst case over rates: normalise routing cost (each session distance 2) to 1
res2, *_ = nclp.solve_cp(net, perms=G, batch=60000, rates_free=True, routing=[2, 2, 2, 2])
out['rates_free_min_cost_per_unit_routing'] = res2.fun
print("rates free: min LP cost / routing cost =", res2.fun)
json.dump(out, open('os_cp.json', 'w'), indent=1); print(out)
