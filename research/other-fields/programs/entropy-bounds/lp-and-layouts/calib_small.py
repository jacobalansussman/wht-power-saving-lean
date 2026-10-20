# calib_small.py -- calibration of nclp.py on networks whose answers are known.
import nclp

def undirected(edges):
    arcs = []
    for (u, v) in edges: arcs += [(u, v), (v, u)]
    return arcs

print("1. path 0-1-2-3, one session 0->3 (routing = 3)")
net = nclp.Net([0, 1, 2, 3], undirected([(0, 1), (1, 2), (2, 3)]), [(0, [3], 1)])
lp = nclp.build_lp(net); r = nclp.solve(lp)

print("2. directed butterfly, two sources each multicast to both sinks (Steiner routing = 8, coding = 7)")
arcs = [('s1', 't1'), ('s1', 'a'), ('s2', 'a'), ('s2', 't2'), ('a', 'b'), ('b', 't1'), ('b', 't2')]
net = nclp.Net(['s1', 's2', 'a', 'b', 't1', 't2'], arcs, [('s1', ['t1', 't2'], 1), ('s2', ['t1', 't2'], 1)])
lp = nclp.build_lp(net); r = nclp.solve(lp)
print("   gain detected = 8 /", r.fun, "=", 8 / r.fun)

print("3. same butterfly UNDIRECTED (both directions of each edge available; routing = 6)")
net = nclp.Net(['s1', 's2', 'a', 'b', 't1', 't2'], undirected([(u, v) for (u, v) in arcs]),
               [('s1', ['t1', 't2'], 1), ('s2', ['t1', 't2'], 1)])
lp = nclp.build_lp(net); r = nclp.solve(lp)

print("4. two-node exchange u<->v, sessions u->v and v->u (XOR on one edge is NOT free: cost 2)")
net = nclp.Net(['u', 'v'], undirected([('u', 'v')]), [('u', ['v'], 1), ('v', ['u'], 1)])
lp = nclp.build_lp(net); r = nclp.solve(lp)

print("5. directed butterfly as TWO UNICASTS s1->t2, s2->t1 with side arcs s1->t1, s2->t2 (classic coding-for-capacity;")
print("   in min-cost form routing = 6 via the middle twice, coding = 5? no: coding uses s1a,s2a,ab,bt1,bt2,s1t1,s2t2 = 7)")
net = nclp.Net(['s1', 's2', 'a', 'b', 't1', 't2'], arcs, [('s1', ['t2'], 1), ('s2', ['t1'], 1)])
lp = nclp.build_lp(net); r = nclp.solve(lp)
