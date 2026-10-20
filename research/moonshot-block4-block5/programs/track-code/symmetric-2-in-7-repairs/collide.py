"""collide.py -- mw2-absorb4: ODD-ORDER COLLISIONS at the shared move.
In the seed-66 limit the two terms of move 1's filter (b1, a b2) sit in distinct frames at distance 2, are MERGED into one
frame by the shared move (ratio of order 2: a collision / diamond), and stay merged at the K frame: K-frame part 1 + a u,
u^2 = 1; the junk forces a^2 = 1 (zero divisor).  Here: filters b_k = w^k b1 with w in B of ODD order such that
 (B) v = D1 w D1^-1 is NOT free (distinct frames at distance 2),
 (A) u1 = Dd b' v b'^-1 Dd^-1 IS free (collision at the shared move),
 (C) u = D1 d u1 d^-1 D1^-1 IS free (still one frame at K).
For every w found: 2-term filter {b1, w b1} and (order 3) 3-term filter {b1, w b1, w^2 b1}; full exact test.
Usage: collide.py lim.pkl nsamples starts"""
import sys, pickle, time
from twist import *
d = pickle.load(open(sys.argv[1], 'rb')); ns = int(sys.argv[2]); starts = int(sys.argv[3])
m, W, H = d['m'], d['W'], d['H']
rec = d['recipe']
base = build_T(m, W, H, rec, d['skseed'], {})
G, R = base.G, base.R; n = G.n
(k1, sm1, D1o), (k2, sm2, Dd), (k3, sm3, D3o) = rec[1], rec[2], rec[3]
D1 = G.D1 if D1o is None else D1o; D3 = G.D1 if D3o is None else D3o
b1 = sm1[0][1]; bp = sm2[0][1]; dl = sm3[0][1]
inv = lambda M: np.conj(M).T
rng = random.Random(31)
found = {}
t0 = time.time()
def test(w):
    v = G.mul(D1, w, inv(D1))
    if G.is_monomial(v): return None
    u1 = G.mul(Dd, bp, v, inv(bp), inv(Dd))
    if not G.is_monomial(u1): return None
    u = G.mul(D3, dl, u1, inv(dl), inv(D3))
    if not G.is_monomial(u): return None
    return order_of(G, w)
# structured candidates: D1^-1 x D1 for x free gives w with v = x free (excluded); so sample w directly
for it in range(ns):
    r = it % 3
    if r == 0: w = rand_free_b(G, rng, m)
    elif r == 1:
        sig = list(range(n)); rng.shuffle(sig); w = G.mono(sig, [rng.randrange(4) for _ in range(n)])
    else:   # products of two engine letters' commutator-like words (more odd orders)
        a_ = rand_free_b(G, rng, m); b_ = rand_free_b(G, rng, m); w = G.mul(a_, b_, inv(a_), inv(b_))
    if is_scalar(w): continue
    o = test(w)
    if o is None: continue
    k = G.key(w)
    if k not in found: found[k] = (o, w)
print(f"collision letters: {len(found)} from {ns} samples; orders {sorted(set(o for o, w in found.values()))} ({time.time()-t0:.0f}s)", flush=True)
items = sorted(found.values(), key=lambda t: (t[0] % 2 == 0, t[0]))
res = []
for idx, (o, w) in enumerate(items[:60]):
    for nterm in ([2, 3] if o % 2 == 1 else [2]):
        rec2 = [(k, list(sm), dd) for (k, sm, dd) in rec]
        src = sm1[0][0]
        terms = [(src, b1)]
        wk = np.eye(n, dtype=complex)
        for t in range(1, nterm):
            wk = G.mul(w, wk); terms.append((src, G.mul(wk, b1)))
        rec2[1] = (k1, terms, D1o)
        try:
            sk = build_T(m, W, H, rec2, d['skseed'], {}, R)
        except Exception as e:
            print("build error", e); continue
        if [1 + len(mv['images']) for mv in sk.moves] != [1 + len(mv['images']) for mv in base.moves]:
            print(f"w {idx}: orbit structure changed"); continue
        if not sk.setup_final(4000):
            print(f"w {idx} order {o} terms {nterm}: closure cap"); continue
        prep_fast(sk)
        N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
        bb = np.inf
        for st in range(starts):
            x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
            x, c = vp_lm(sk, N, free_idx, x0, iters=120)
            bb = min(bb, c)
        dims = [(nr, nc) for (nr, nc, E) in sk.final]
        res.append((bb, o, nterm, idx))
        print(f"w {idx} order {o} filter terms {nterm}: closure {dims} full res^2 {bb:.4f}", flush=True)
        if bb < 1e-12:
            print("*** FULL HIT ***", flush=True)
            pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": {}, "skseed": d['skseed']}, open(f"HIT_collide_{idx}_{nterm}.pkl", "wb"))
res.sort(); print("best", res[:10])
