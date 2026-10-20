# layout_cut_q.py -- fast certified lower bound on the directed cut LP of a fixed layout D: contract Gamma_m by the
# vector of distances to the endpoints of D (endpoints stay singletons), and solve the cut LP of the quotient
# (bundled arcs, cost 1 per bundle use). Every cut of the quotient is a cut of Gamma_m with the same arcs and the
# same requirement, so quotient value <= cut LP of D <= moves of any code for D. A quotient value > M certifies
# that no code for D has M moves.
import sys, json, random, time, collections
import numpy as np
import dpg
from layout_cut2 import cut_lp

def profile_quotient(n, adj, dist, D, extra_refs=()):
    ends = sorted({x for p in D for x in p}) + list(extra_refs)
    key = [tuple(dist[v][e] for e in ends) for v in range(n)]
    cls = {}; part = [cls.setdefault(k, len(cls)) for k in key]
    k = len(cls)
    qarcs = sorted({(part[u], part[v]) for u in range(n) for v in adj[u] if part[u] != part[v]})
    qD = [(part[s], part[t]) for (s, t) in D]
    return k, qarcs, qD, part

def certified(n, adj, dist, D, M, extra=0, rng=None):
    refs = []
    rng = rng or random.Random(0)
    for rnd in range(extra + 1):
        k, qarcs, qD, part = profile_quotient(n, adj, dist, D, refs)
        # quotient distances for seeds
        val, npat, rounds = cut_lp(k, qarcs, qD, stop_above=M)
        if val is not None and val > M + 1e-7:
            return val, k, True
        refs.append(rng.randrange(n))
    return val, k, False

if __name__ == '__main__':
    m = int(sys.argv[1])
    verts, idx, adj, dist = dpg.build(m); n = len(verts)
    tests = json.loads(sys.argv[2])
    for D, M in tests:
        t0 = time.time()
        k, qarcs, qD, part = profile_quotient(n, adj, dist, [tuple(p) for p in D])
        val, npat, rounds = cut_lp(k, qarcs, [tuple(p) for p in qD])
        print(f"layout {D}: profile quotient {k} cells, {len(qarcs)} bundles: cut LP {val:.6f} (> {M}? {val > M + 1e-7})  {time.time()-t0:.1f}s")
