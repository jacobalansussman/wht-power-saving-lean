"""lastmulti.py -- mw2-absorb4 (resumed run): does a MULTI-SUMMAND last move escape the dilemma (NOTE.md s4 covers
single-summand last moves)?  Near-miss limit + extra summands in the last orbit move: sources = the shared generator
(other letters) or the filters / layer-1 generators, letters = progress letters of their terms (any distance).
Full exact-structure test each.  Usage: lastmulti.py recipe.pkl ntries starts [cap]"""
import sys, pickle, time
from odd_twist import *

d = pickle.load(open(sys.argv[1], 'rb'))
ntries, starts = int(sys.argv[2]), int(sys.argv[3]); cap = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
m, W, H = d['m'], d['W'], d['H']
rec = [(k, list(sm), dd) for (k, sm, dd) in d['recipe']]
base = build_T(m, W, H, rec, d['skseed'], {})
G, R = base.G, base.R
tl = len(rec) - 1
kind, summ, Dm = rec[tl]
Dmat = G.D1 if Dm is None else Dm
rng = random.Random(47)
srcs = list(range(len(base.gens) - (1 + len(base.moves[-1]['images']))))
cands = []
for src in srcs:
    for (l, g) in base.gens[src]:
        if G.dist(g) < 1:
            continue
        for b in G.progress_letters(g, rng, want=8, Dmat=Dmat, steps=1):
            cands.append((src, b))
print("candidate extra summands", len(cands), "from generators", sorted(set(s for s, b in cands)), flush=True)
rng.shuffle(cands)
t0 = time.time(); res = []
for k in range(ntries):
    extra = [cands[k % len(cands)]] if k < len(cands) else [rng.choice(cands), rng.choice(cands)]
    rec2 = [(kk, list(sm), dd) for (kk, sm, dd) in rec]
    rec2[tl] = (kind, summ + extra, Dm)
    try:
        sk = build_T(m, W, H, rec2, d['skseed'], {}, R)
    except AssertionError:
        continue
    if [1 + len(mv['images']) for mv in sk.moves] != [1 + len(mv['images']) for mv in base.moves]:
        continue
    if not sk.setup_final(cap):
        print(f"try {k}: closure cap"); continue
    prep_fast(sk)
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    best = np.inf
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        x, c = vp_lm(sk, N, free_idx, x0, iters=100)
        best = min(best, c)
    dims = [(nr, nc) for (nr, nc, E) in sk.final]
    res.append((best, k, [s for s, b in extra]))
    print(f"try {k}: extra summands from {[s for s, b in extra]}: closure {dims} full res^2 {best:.6f} ({time.time()-t0:.0f}s)", flush=True)
    if best < 1e-12:
        print("*** FULL HIT ***", flush=True)
        pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "skseed": d['skseed']}, open(f"HIT_lastmulti_{k}.pkl", "wb"))
res.sort(); print("best", res[:8])
