# wscreen.py -- session 2. Exact screen for W-array coherent units at block m (general W; w3screen.py is W = 3):
#  1. two-class stretchable layouts (Lemma W, PROOF.md section 8): cut value W*m, proved.
#  2. otherwise pull back along psi(v) = (d(u_1,v), ..., d(u_W,v)), u_i one end of pair i (up to `maxor` choices),
#     then along the profile to all 2W terminals, then (optionally) the full Gamma_m: an exact dual certificate of value
#     > W*m - 1 proves M >= W*m for every coherent unit with that layout (M is an integer).
#  Usage: wscreen.py m W N seed [maxor] [full]
import sys, json, random, time, itertools
import numpy as np
import w3cut, w3image, stretch2
from w3screen import exact_value

def screen(g, P, maxor=8, full=False):
    m = g.m; W = len(P); thr = W * m - 1; best = 0
    ors = list(itertools.product((0, 1), repeat=W))[:maxor]
    for o in ors:
        us = [P[i][o[i]] for i in range(W)]
        nH, EH, lab, keys = w3image.image_graph(g, us)
        v = exact_value(w3image.cut_value(nH, EH, [(lab[s], lab[t]) for s, t in P]))
        best = max(best, v)
        if v > thr: return v, 'imageW', nH
    nH, EH, lab, keys = w3image.image_graph(g, sorted({x for p in P for x in p}))
    v = exact_value(w3image.cut_value(nH, EH, [(lab[s], lab[t]) for s, t in P]))
    best = max(best, v)
    if v > thr: return v, 'image2W', nH
    if full:
        v = exact_value(w3image.cut_value(g.n, g.edges, P))
        best = max(best, v)
        if v > thr: return v, 'full', g.n
    return best, None, None

if __name__ == '__main__':
    m, W, N, seed = map(int, sys.argv[1:5])
    maxor = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    full = len(sys.argv) > 6 and sys.argv[6] == 'full'
    t0 = time.time(); g = w3cut.G(m); rng = random.Random(seed)
    print(f"Gamma_{m}: {g.n} frames ({time.time()-t0:.0f}s); W = {W}", flush=True)
    nst = 0; cert = []; fails = []
    for tr in range(N):
        P = stretch2.rand_layout(g, rng, W)
        if stretch2.two_class(g.D, m, P): nst += 1; continue
        v, lev, nH = screen(g, P, maxor, full)
        if lev:
            cert.append(dict(P=P, level=lev, H=nH, value=str(v), value_float=float(v)))
            print(f"  certified at level {lev} ({nH} points): {P} exact value {float(v):.5f} ({time.time()-t0:.0f}s)", flush=True)
        else:
            fails.append(dict(P=P, best=str(v))); print(f"  OPEN {P}: best exact value {float(v):.5f}", flush=True)
        print(f"  {tr+1}: stretchable {nst}, certified by pull-back {len(cert)}, open {len(fails)}", flush=True)
        json.dump(dict(m=m, W=W, N=tr + 1, seed=seed, stretchable=nst, cert=cert, fails=fails), open(f"wscreen_m{m}_W{W}_s{seed}.json", "w"))
    print(f"m={m} W={W} seed={seed}: {N} random layouts (distinct frames): two-class stretchable {nst}, certified by "
          f"pull-back {len(cert)}, open {len(fails)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(m=m, W=W, N=N, seed=seed, stretchable=nst, cert=cert, fails=fails), open(f"wscreen_m{m}_W{W}_s{seed}.json", "w"))
