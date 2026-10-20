import nclp, time
arcs = [('s1', 't1'), ('s1', 'a'), ('s2', 'a'), ('s2', 't2'), ('a', 'b'), ('b', 't1'), ('b', 't2')]
net = nclp.Net(['s1', 's2', 'a', 'b', 't1', 't2'], arcs, [('s1', ['t1', 't2'], 1), ('s2', ['t1', 't2'], 1)])
print("butterfly full:", nclp.solve(nclp.build_lp(net), verbose=False).fun, " closed:", nclp.solve(nclp.build_closed_lp(net)).fun)
net = nclp.Net([0, 1, 2, 3], [(0,1),(1,0),(1,2),(2,1),(2,3),(3,2)], [(0, [3], 1)])
print("path closed:", nclp.solve(nclp.build_closed_lp(net)).fun)
# OS: just report the size of the closed LP
L, R = ['a', 'b'], ['x', 'y', 'z']
arcs = []
for u in L:
    for v in R: arcs += [(u, v), (v, u)]
net = nclp.Net(L + R, arcs, [('x', ['y'], 1), ('y', ['z'], 1), ('z', ['x'], 1), ('a', ['b'], 1)])
rot = {'a': 'a', 'b': 'b', 'x': 'y', 'y': 'z', 'z': 'x'}
perm = [None] * net.n
for i, (s, t, r) in enumerate(net.sessions):
    perm[i] = next(j for j, (s2, t2, r2) in enumerate(net.sessions) if s2 == rot[s] and t2 == [rot[u] for u in t])
for a, (u, v) in enumerate(arcs): perm[net.ns + a] = net.ns + arcs.index((rot[u], rot[v]))
G = nclp.close_group([tuple(perm)], net.n)
t0 = time.time()
cl = nclp.closure_map(net)
import numpy as np
print("OS: distinct closed sets", len(np.unique(cl)), "of", 1 << net.n, f"({time.time()-t0:.1f}s)")
