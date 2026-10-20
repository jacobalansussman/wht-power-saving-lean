# hard4.py -- block 4: hill-climb layouts of three disjoint opposite pairs to large cross deficiency e(D) (the layouts
# Lemma W+ handles worst), then run the pull-back screen on them (exact dual; integrality: value > 11 => M >= 12).
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import sys, random, time, json
import w3cut, edef, w3screen
m = 4; g = w3cut.G(m); rng = random.Random(int(sys.argv[2])); t0 = time.time(); out = []
seen = set()
for r in range(int(sys.argv[1])):
    P = w3cut.rand_3k2(g, rng); e = edef.edef(g.D, m, P)
    for it in range(600):
        k = rng.randrange(3); end = rng.randrange(2); s, t = P[k]; x = (s, t)[end]; y = (t, s)[end]
        nb = rng.choice(g.adj[x])
        if g.D[nb, y] != m or nb in {v for p in P for v in p}: continue
        Q = list(P); Q[k] = (nb, y) if end == 0 else (y, nb); e2 = edef.edef(g.D, m, Q)
        if e2 >= e: P, e = Q, e2
    if e < 2: continue
    key = tuple(sorted(tuple(sorted(p)) for p in P))
    if key in seen: continue
    seen.add(key)
    w3screen.THR_FRAC = False
    v, lev, o = w3screen.screen_one(g, P, levels=('image3',))
    out.append(dict(P=P, e=e, value=float(v), level=lev))
    print(f"e={e} {P}: best pull-back exact value {float(v):.5f} level {lev} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(f"hard4_s{sys.argv[2]}.json", "w"))
