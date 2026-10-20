"""pathscan.py -- mw2-absorb4: odd-order filter letters that keep the K pair.  In a limit skeleton, the two-term filter
b1 + a b2 of move t_f reaches the K frame through the path P = h_K h_src^-1 (h_K the K-frame term of the last generator,
h_src the filter's source term).  b2 = w b1 keeps both terms in the K frame iff u = P w P^-1 is free, P = h_K h_src^-1 b1^-1.
Here: sample many free letters w (engine + Wd monomials), keep those with u monomial, record ord(w), and run the full
exact test for each kept w (odd orders first).  Usage: pathscan.py lim.pkl t_f nsamples starts"""
import sys, pickle, time
from twist import *
d = pickle.load(open(sys.argv[1], 'rb'))
t_f, nsamp, starts = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
m, W, H = d['m'], d['W'], d['H']
rec = d['recipe']; tw = d.get('twist', {})
base = build_T(m, W, H, rec, d['skseed'], tw)
G, R = base.G, base.R; n = G.n
kind, summ, Dm = rec[t_f]
src0, b1 = summ[0]
# path: source term h_s of generator src0 (array 0), K-frame terms of the generators
srcterms = [(l, g) for (l, g) in base.gens[src0]]
Kterms = [(i, l, g) for i, T in enumerate(base.gens) for (l, g) in T if G.dist(g) == 0]
print("source terms", [(l, G.dist(g)) for (l, g) in srcterms], "K-frame terms", [(i, l) for (i, l, g) in Kterms])
Ps = []
for (i, lK, gK) in Kterms:
    for (ls, gs) in srcterms:
        if ls != lK: continue
        P = G.mul(G.MATS[gK], np.conj(G.MATS[gs]).T, np.conj(b1).T)    # u = P w P^-1 for b2 = b1 w
        Ps.append(P)
rng = random.Random(23)
keep = {}
t0 = time.time()
for it in range(nsamp):
    if it % 2 == 0:
        w = rand_free_b(G, rng, m)
    else:
        sig = list(range(n)); rng.shuffle(sig); w = G.mono(sig, [rng.randrange(4) for _ in range(n)])
    if is_scalar(w): continue
    if any(G.is_monomial(G.mul(P, w, np.conj(P).T)) for P in Ps):
        k = G.key(w)
        if k not in keep:
            keep[k] = (order_of(G, w), w)
print(f"kept {len(keep)} of {nsamp} letters; orders {sorted(set(o for o, w in keep.values()))} ({time.time()-t0:.0f}s)", flush=True)
items = sorted(keep.values(), key=lambda t: (t[0] % 2 == 0, t[0]))
res = []
for idx, (o, w) in enumerate(items):
    b2 = G.mul(w, b1)            # left ratio: b2 b1^-1 = w, and u = P w P^-1 is the K-pair letter
    rec2 = [(k, list(sm), dd) for (k, sm, dd) in rec]
    sm = rec2[t_f][1]
    if len(sm) > 1: sm[1] = (sm[1][0], b2)
    else: sm.append((src0, b2))
    try:
        sk = build_T(m, W, H, rec2, d['skseed'], tw, R)
    except Exception as e:
        continue
    if [1 + len(mv['images']) for mv in sk.moves] != [1 + len(mv['images']) for mv in base.moves]:
        continue
    if not sk.setup_final(3000):
        print(f"w {idx} order {o}: closure cap"); continue
    prep_fast(sk)
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    bb = np.inf
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        x, c = vp_lm(sk, N, free_idx, x0, iters=100)
        bb = min(bb, c)
    dims = [(nr, nc) for (nr, nc, E) in sk.final]
    res.append((bb, o, idx))
    print(f"w {idx} order {o}: closure {dims} full res^2 {bb:.4f}", flush=True)
    if bb < 1e-12:
        print("*** FULL HIT ***"); pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": tw, "skseed": d['skseed']}, open(f"HIT_path_{t_f}_{idx}.pkl", "wb"))
res.sort(); print("best", res[:8])
