"""census1.py -- w5-twoshared: frame census of the top two layers of the 2+2+1+2 shape (single letters, F diagonal:
top 0 kernelises array 0 through merge 1, top 1 kernelises array 1 through merge 2).  Routes are built first
(array 0: b0, a10, g1, mu0; array 1: b1, a21, g2, mu1), cross letters a11, a20 random progress letters.
For every top e and every non-route product y_e q_kj: its layer, kernel-frame membership, and whether it lies in a
frame of V (B'q_k'j' for the same array j) -- the only place a layer-3 junk can be absorbed (Lemma C, NOTE).
Usage: census1.py seed N"""
from ts import *
seed, N = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed)
stats = {}; t0 = time.time(); examples = []
for it in range(N):
    r0 = route(rng); r1 = route(rng)
    if r0 is None or r1 is None: continue
    b0, a10, p10, P0 = r0; b1, a21, p21, P1 = r1
    g1, mu0 = P0[rng.randrange(len(P0))]; g2, mu1 = P1[rng.randrange(len(P1))]
    h0 = G.mul(G.D1, b0); h1 = G.mul(G.D1, b1)
    a11 = prog_letter(rng, [G.intern(h1)]); a20 = prog_letter(rng, [G.intern(h0)])
    p = {(1, 0): p10, (2, 1): p21, (1, 1): G.mul(G.D1, a11, h1), (2, 0): G.mul(G.D1, a20, h0)}
    g = {1: g1, 2: g2}
    q = {kj: G.mul(G.D1, g[kj[0]], p[kj]) for kj in p}
    y = {0: G.mul(G.D1, mu0), 1: G.mul(G.D1, mu1)}
    assert is_kernel(G.mul(y[0], q[(1, 0)])) and is_kernel(G.mul(y[1], q[(2, 1)]))
    key = []
    for e, rt in ((0, (1, 0)), (1, (2, 1))):
        for kj in [(1, 0), (2, 0), (1, 1), (2, 1)]:
            if kj == rt: continue
            J = G.mul(y[e], q[kj]); L = lam(J)
            tag = str(L)
            if L == 4 and is_kernel(J): tag = 'K'
            if L == 3:
                inV = [k2 for k2 in (1, 2) if same_frame(J, q[(k2, kj[1])])]
                if inV: tag = '3V'
            key.append(tag)
    lq = tuple(lam(q[kj]) for kj in [(1, 1), (2, 0)])
    k = (lq, tuple(key))
    stats[k] = stats.get(k, 0) + 1
print(f"seed {seed} N {N} ({time.time()-t0:.0f}s); key = (lam q11, lam q20), tags for top0: q20 q11 q21, top1: q10 q20 q11")
tot = sum(stats.values())
ok = 0
for k, c in sorted(stats.items(), key=lambda t: -t[1]):
    good = all(t in ('2', '3V', 'K', '1') for t in k[1])
    ok += c if good else 0
    print(k, c, 'ALLOWED' if good else '')
print("total", tot, "allowed (no free-standing layer-3/4 junk)", ok)
