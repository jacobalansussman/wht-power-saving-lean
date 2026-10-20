import time, random, sys, dpg
from layout_cut2 import cut_lp
m = int(sys.argv[1]); W = int(sys.argv[2])
verts, idx, adj, dist = dpg.build(m); n = len(verts)
arcs = [(u, v) for u in range(n) for v in adj[u]]
opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
rng = random.Random(7)
D = [(0, opp[0][0])]
while len(D) < W:
    s = rng.randrange(n); D.append((s, rng.choice(opp[s])))
t0 = time.time(); v = cut_lp(n, arcs, D, dist=dist, verbose=True); print(W, v, f"{time.time()-t0:.1f}s")
