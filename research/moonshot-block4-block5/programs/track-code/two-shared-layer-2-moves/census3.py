"""census3.py -- census2 with periodic output and layer-2 tags.  Layer-2 junk J of top e on array j is tagged
'2u' if its frame is a frame of some u_k' main term (B'p_k'j: absorbable by a C[B']-multiple of u_k', frame-wise),
'2' otherwise (needs an internal producer of exactly zeta J inside t_e + eps_e v).  Usage: census3.py seed N variant"""
import sys, random, time
from sp import *
seed, N, variant = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rng = random.Random(seed)
def route(rng):
    while True:
        b = letter(rng); s = mv(b, I8)
        a = prog_letter(rng, [s])
        if a is None: continue
        p = mv(a, s)
        for _ in range(400):
            g = letter(rng); q = mv(g, p)
            if lam(q) == 3 and dist(lag(q), LK) == 1:
                mu = to_K(rng, q)
                if mu is not None: return dict(b=b, s=s, a=a, p=p, g=g, q=q, mu=mu, y=comp(D, mu))
stats = {}; t0 = time.time(); allowed = 0; tot = 0
def report():
    print(f"--- census3 seed {seed} variant {variant}: {tot} configurations ({time.time()-t0:.0f}s); ALLOWED at layers 3-4: {allowed}", flush=True)
    for k, c in sorted(stats.items(), key=lambda t: -t[1])[:30]:
        print(c, k, 'ALLOWED' if not any(x in ('4', '3') for x in k[:-1]) else '', flush=True)
for it in range(N):
    R0 = route(rng); R1 = route(rng)
    Q = {(1, 0): R0['q'], (2, 1): R1['q']}; Pm = {(1, 0): R0['p'], (2, 1): R1['p']}
    if variant == 'cross':
        a11 = letter(rng); a20 = letter(rng)
        Pm[(1, 1)] = mv(a11, R1['s']); Pm[(2, 0)] = mv(a20, R0['s'])
        Q[(1, 1)] = mv(R0['g'], Pm[(1, 1)]); Q[(2, 0)] = mv(R1['g'], Pm[(2, 0)])
    Y = {0: R0['y'], 1: R1['y']}
    tags = []
    for e, rt in ((0, (1, 0)), (1, (2, 1))):
        for kj in sorted(Q):
            if kj == rt: continue
            q = Q[kj]; J = comp(Y[e], q); L = lam(J); LJ = lag(J)
            if L <= 1: tag = '<=1'
            elif L == 2: tag = '2u' if any(lag(Pm[(k2, j2)]) == LJ for (k2, j2) in Pm if j2 == kj[1]) else '2'
            elif L == 4: tag = 'K' if LJ == LK else '4'
            else:
                if lam(q) <= 2: tag = '3c'
                elif any(lag(Q[(k2, j2)]) == LJ for (k2, j2) in Q if j2 == kj[1]): tag = '3V'
                else: tag = '3'
            tags.append(tag)
    key = tuple(tags) + (('lq', tuple(lam(Q[kj]) for kj in sorted(Q))),)
    stats[key] = stats.get(key, 0) + 1; tot += 1
    if not any(t in ('4', '3') for t in tags): allowed += 1
    if tot % 40 == 0: report()
report()
