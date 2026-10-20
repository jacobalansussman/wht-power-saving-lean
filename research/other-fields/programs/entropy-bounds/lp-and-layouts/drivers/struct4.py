# struct4.py -- structured 3-matchings at block 4 (the only 3-session shape where cuts can fail): pairs whose
# endpoints are close to each other, 'parallel' (s2 ~ s1, t2 ~ t1) and 'crossed' (s2 ~ t1, t2 ~ s1), solved by the
# PDLP cut LP with a rigorous bound. Cut value > 11 means no coherent 3 in 11 (or 3 in 10) with that layout.
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import random, sys, time, json
import dpg, layout_cut_pdlp as L
m, W = 4, 3
verts, idx, adj, dist = dpg.build(m); n = len(verts)
arcs = [(u, v) for u in range(n) for v in adj[u]]
rng = random.Random(int(sys.argv[2])); t0 = time.time(); out = []
def near_pair(a, b, used):
    # s adjacent to a, t adjacent to b, s and t opposite, both new
    for _ in range(2000):
        s = rng.choice(adj[a]); t = rng.choice(adj[b])
        if s not in used and t not in used and s != t and dist[s][t] == m: return (s, t)
    return None
for tr in range(int(sys.argv[1])):
    s1 = 0; t1 = next(t for t in range(n) if dist[0][t] == m and rng.random() < 0.01) if False else rng.choice([t for t in range(n) if dist[0][t] == m])
    D = [(s1, t1)]
    mode = ['parallel', 'crossed', 'mixed'][tr % 3]
    used = {s1, t1}
    p2 = near_pair(s1, t1, used) if mode in ('parallel', 'mixed') else near_pair(t1, s1, used)
    if p2 is None: continue
    D.append(p2); used |= set(p2)
    p3 = near_pair(p2[0], p2[1], used) if mode == 'parallel' else near_pair(p2[1], p2[0], used)
    if p3 is None: continue
    D.append(p3)
    A, c, lo, hi, nvar, npat = L.build(n, arcs, D, W * m)
    primal, LB, st = L.solve(A, c, lo, hi, nvar, W * m)
    out.append(dict(mode=mode, D=D, primal=primal, lb=LB))
    print(f"{mode} {D}: PDLP {st}, primal {primal:.5f}, rigorous LB {LB:.5f} (routing 12)  ({time.time()-t0:.0f}s)"); sys.stdout.flush()
json.dump(out, open(f"struct4_s{sys.argv[2]}.json", "w"))
