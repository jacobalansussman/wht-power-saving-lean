# os_closed.py -- Okamura-Seymour K_{2,3} calibration with the closure-reduced full Shannon LP (same optimum as the
# full LP: h(S) = h(cl S) follows from the closure rows by submodularity), then Ingleton / DFZ cutting planes,
# then the worst case over rates. Prints the dual weights of the binding rows for the rate-1 instance.
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import json, sys, time, numpy as np
import nclp
L, R = ['a', 'b'], ['x', 'y', 'z']
arcs = []
for u in L:
    for v in R: arcs += [(u, v), (v, u)]
sessions = [('x', ['y'], 1), ('y', ['z'], 1), ('z', ['x'], 1), ('a', ['b'], 1)]
net = nclp.Net(L + R, arcs, sessions)
rot = {'a': 'a', 'b': 'b', 'x': 'y', 'y': 'z', 'z': 'x'}
perm = [None] * net.n
for i, (s, t, r) in enumerate(sessions):
    perm[i] = next(j for j, (s2, t2, r2) in enumerate(sessions) if s2 == rot[s] and t2 == [rot[u] for u in t])
for a, (u, v) in enumerate(arcs): perm[net.ns + a] = net.ns + arcs.index((rot[u], rot[v]))
G = nclp.close_group([tuple(perm)], net.n)
out = {}
lp = nclp.build_closed_lp(net, perms=G)
primal, LB, res = nclp.solve_pdlp(lp, 4)
out['shannon'] = primal; out['shannon_rigorous_lb'] = LB
json.dump(out, open('os_closed.json', 'w'), indent=1)
h = res.x; canon = lp['canon']
letters = sorted({1 << i for i in range(net.n)} | {(1 << i) | (1 << j) for i in range(net.n) for j in range(i + 1, net.n)})
cuts = []
for rnd in range(3):
    ni, vi = nclp.find_cuts(h, canon, net.n, letters, 'ingleton', tries=300000)
    nd, vd = nclp.find_cuts(h, canon, net.n, letters, 'dfz', tries=300000)
    print(f"  linear-rank round {rnd}: Ingleton violated {len(ni)} (max {vi:.2e}), DFZ violated {len(nd)} (max {vd:.2e})"); sys.stdout.flush()
    if not ni and not nd: break
    cuts += ni + nd
    primal, LB, res = nclp.solve_pdlp(lp, 4, extra_ub=cuts); h = res.x
out['with_linear_rank_cuts'] = LB
json.dump(out, open('os_closed.json', 'w'), indent=1)
lp2 = nclp.build_closed_lp(net, perms=G, rates_free=True, routing=[2, 2, 2, 2])
p2, LB2, res2 = nclp.solve_pdlp(lp2, 1.0)
out['rates_free_cost_per_unit_routing'] = LB2
json.dump(out, open('os_closed.json', 'w'), indent=1); print(out)
