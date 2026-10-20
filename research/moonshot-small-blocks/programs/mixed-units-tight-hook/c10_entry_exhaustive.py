# c10_entry_exhaustive.py -- ms-mixed: EXHAUSTIVE single-term entry rule for a rank-2 slot, two primes.
# Classes of D b (b free) up to left multiplication by free letters (canonical form of the 8x8 matrix: rows scaled to
# leading 1 and sorted -- monomial left action), as in wave16/t514-tight wide_g2.py.  Engine: b over AGL x phases x
# shifts (random until 14 classes stop growing); widened: b over S8 x Z4^8 (random until 1680 classes stop growing).
# Test: rho(D b) E^Y_R in Y for Y = N, Nperp (rank-2 target R = D1 D2) and the same for the kernel K.
import sys, time, random, collections
from mixlib import *
t0 = time.time()
pi = int(sys.argv[1]) if len(sys.argv) > 1 else 0
p = primes_1mod4_below(1 << 20, 2)[pi]
H = Hook(p); F = H.F
Pots = {(lab, yn): Pot(H, [T], yn) for lab, T in (('K', KK), ('R', R12)) for yn in ('N', 'Nperp')}
def canon(g):
    M = g.modp(F)
    rows = []
    for r in M % p:
        nz = np.nonzero(r)[0]; inv = pow(int(r[nz[0]]), p - 2, p)
        rows.append(tuple(int(v) for v in (r * inv) % p))
    return tuple(sorted(rows))
def t2single(P, g):
    R = H.rho(g)[:, P.Yi]
    return not P.F.mm(R[P.Oi], P.E[0]).any()
rng = random.Random(5 + pi)
for game, gen, NS in (('engine', rand_engine, 4000), ('widened', rand_widened, 60000)):
    reps = {}
    for t in range(NS):
        b = gen(rng)
        k = canon(D1 @ b)
        if k not in reps: reps[k] = b
    tab = collections.Counter()
    for k, b in reps.items():
        g = D1 @ b
        tab[tuple(t2single(Pots[key], g) for key in (('K', 'N'), ('K', 'Nperp'), ('R', 'N'), ('R', 'Nperp')))] += 1
    print(f"p = {p}, {game}: {len(reps)} classes of D b from {NS} samples; (K:N, K:Nperp, R:N, R:Nperp) -> count: {dict(tab)}  [{time.time()-t0:.0f}s]", flush=True)
