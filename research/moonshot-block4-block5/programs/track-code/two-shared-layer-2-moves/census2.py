"""census2.py -- w5-twoshared: frame census of the top two layers, 2+2+1+2 shape (engine game, single letters).
Routes first: array 0 = x0 -> s0 -> u1 -> v -> t0 (letters b0, a10, g1, mu0), array 1 = x1 -> s1 -> u2 -> v -> t1
(b1, a21, g2, mu1); cross letters a11 (s1 into u1), a20 (s0 into u2) random (any layer); variant 'priv': no cross terms.
For each top e and each non-route term q_kj of v (layer 3 only; layer <= 2 terms of v are noted), the product
J = y_e q_kj is tagged:  '<=2' (drops), 'K' (kernel frame: F gets another entry), '4' (non-kernel top-layer junk:
FORBIDDEN, Lemma A), '3V' (layer 3 in a frame of V of the same array), '3' (layer 3 elsewhere: FORBIDDEN, Lemma C),
'3c' (layer 3 where q_kj itself has layer 2: producer (C) possible).  Usage: census2.py seed N variant"""
import sys, random, time
from sp import *
seed, N, variant = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rng = random.Random(seed)
def route(rng):
    while True:
        b = letter(rng); s = mv(b, I8)
        a = prog_letter(rng, [s]); 
        if a is None: continue
        p = mv(a, s)
        for _ in range(400):
            g = letter(rng); q = mv(g, p)
            if lam(q) == 3 and dist(lag(q), LK) == 1:
                mu = to_K(rng, q)
                if mu is not None: return dict(b=b, s=s, a=a, p=p, g=g, q=q, mu=mu, y=comp(D, mu))
stats = {}; t0 = time.time(); allowed = 0; tot = 0
for it in range(N):
    R0 = route(rng); R1 = route(rng)
    Q = {(1, 0): R0['q'], (2, 1): R1['q']}
    if variant == 'cross':
        a11 = letter(rng); a20 = letter(rng)
        Q[(1, 1)] = mv(R0['g'], mv(a11, R1['s']))
        Q[(2, 0)] = mv(R1['g'], mv(a20, R0['s']))
    Y = {0: R0['y'], 1: R1['y']}
    tags = []
    for e, rt in ((0, (1, 0)), (1, (2, 1))):
        for kj in sorted(Q):
            if kj == rt: continue
            q = Q[kj]; J = comp(Y[e], q); L = lam(J); LJ = lag(J)
            if L <= 2: tag = '<=2'
            elif L == 4: tag = 'K' if LJ == LK else '4'
            else:
                if lam(q) <= 2: tag = '3c'
                elif any(lag(Q[(k2, j2)]) == LJ for (k2, j2) in Q if j2 == kj[1]): tag = '3V'
                else: tag = '3'
            tags.append((e, kj, tag))
    key = tuple(t[2] for t in tags) + (('lq', tuple(lam(Q[kj]) for kj in sorted(Q))),)
    stats[key] = stats.get(key, 0) + 1; tot += 1
    if not any(t[2] in ('4', '3') for t in tags): allowed += 1
print(f"census2 seed {seed} variant {variant}: {tot} configurations ({time.time()-t0:.0f}s)")
print("tag order: for top 0 then top 1, the non-route q_kj in sorted (k,j) order; last = layers of all q_kj")
for k, c in sorted(stats.items(), key=lambda t: -t[1])[:40]:
    print(c, k, 'ALLOWED' if not any(x in ('4', '3') for x in k[:-1]) else '')
print("ALLOWED (no forbidden junk):", allowed, "of", tot)
