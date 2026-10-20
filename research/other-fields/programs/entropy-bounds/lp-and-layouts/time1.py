import time, random, dpg
from layout_cut2 import cut_lp
verts, idx, adj, dist = dpg.build(3); n = len(verts)
arcs = [(u, v) for u in range(n) for v in adj[u]]
opp = [[t for t in range(n) if dist[s][t] == 3] for s in range(n)]
rng = random.Random(7)
for W in (3, 4):
    D = [(0, opp[0][0])]
    while len(D) < W:
        s = rng.randrange(n); D.append((s, rng.choice(opp[s])))
    t0 = time.time(); v = cut_lp(n, arcs, D); print(W, v, f"{time.time()-t0:.1f}s")
