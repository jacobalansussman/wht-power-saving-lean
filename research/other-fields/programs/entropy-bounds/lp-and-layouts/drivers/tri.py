# tri.py -- structured W = 5 layouts at block 3 modelled on the block-2 layout where cuts fail
# ([(0,7),(3,1),(4,10),(0,10),(10,7)]: a triangle of mutually opposite frames 0, 7, 10 with sessions 0->7, 0->10,
# 10->7, and two pendant sessions, one into the triangle). Is the directed cut LP <= 14 (cuts cannot exclude a
# coherent 5 in 14) or > 14?
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import random, sys, time, json
import dpg
from layout_cut2 import cut_lp
m = 3
verts, idx, adj, dist = dpg.build(m); n = len(verts)
arcs = [(u, v) for u in range(n) for v in adj[u]]
opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
rng = random.Random(11); t0 = time.time(); out = []
for trial in range(int(sys.argv[1])):
    a = 0; b = opp[a][rng.randrange(len(opp[a]))]
    cs = [c for c in opp[a] if dist[b][c] == m]
    c = rng.choice(cs)
    D = [(a, b), (a, c), (c, b)]
    # pendant 1: a session into c from a frame adjacent to a (as 4 -> 10 with 4 ~ 0 in the block-2 layout)
    s1 = rng.choice([v for v in adj[a] if dist[v][c] == m] or [v for v in range(n) if dist[v][c] == m])
    D.append((s1, c))
    # pendant 2: a session between two frames adjacent to the triangle
    s2 = rng.choice(adj[rng.choice([a, b, c])])
    t2 = rng.choice([t for t in opp[s2] if any(dist[t][x] == 1 for x in (a, b, c))] or opp[s2])
    D.append((s2, t2))
    v, npat, r = cut_lp(n, arcs, D, dist=dist, stop_above=14)
    out.append(dict(D=D, cut=v))
    print(f"trial {trial}: layout {D}: cut LP {'> 14 (certified ' + format(v, '.4f') + ')' if v > 14 + 1e-7 else format(v, '.6f')}  ({time.time()-t0:.0f}s)"); sys.stdout.flush()
json.dump(out, open('tri.json', 'w'))
