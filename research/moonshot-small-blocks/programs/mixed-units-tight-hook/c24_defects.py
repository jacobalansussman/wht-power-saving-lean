# c24_defects.py -- ms-mixed: single-term T2 DEFECTS for a two-row slot: c_b = rank of pi_O rho(D b) on E_R (f = 5),
# over all 14 engine and 1680 widened classes of D b, for Y = N and Nperp.  (A waste budget for near-tight units: a unit
# with hook slack s x 80 can afford total T2 defect at most s x 80 at f = 5; this table is the input such an argument
# would need.)
import random, collections
from mixlib import *
p = primes_1mod4_below(1 << 20, 1)[0]
H = Hook(p); F = H.F
Ps = {(lab, yn): Pot(H, [T], yn) for lab, T in (('R', R12), ('K', KK)) for yn in ('N', 'Nperp')}
def canon(g):
    M = g.modp(F); rows = []
    for r in M % p:
        nz = np.nonzero(r)[0]; inv = pow(int(r[nz[0]]), p - 2, p)
        rows.append(tuple(int(v) for v in (r * inv) % p))
    return tuple(sorted(rows))
rng = random.Random(5)
for game, gen, NS in (('engine', rand_engine, 4000), ('widened', rand_widened, 60000)):
    reps = {}
    for t in range(NS):
        b = gen(rng); k = canon(D1 @ b)
        if k not in reps: reps[k] = b
    tab = collections.Counter()
    for k, b in reps.items():
        g = D1 @ b
        row = []
        for key in (('R', 'N'), ('R', 'Nperp'), ('K', 'N')):
            P = Ps[key]
            Rm = H.rho(g)[:, P.Yi]
            row.append(F.rank(F.mm(Rm[P.Oi], P.E[0])))
        tab[tuple(row)] += 1
    print(f"{game}: {len(reps)} classes; (defect R:N, defect R:Nperp, defect K:N) -> count:")
    for k, v in sorted(tab.items()): print("   ", k, v)
