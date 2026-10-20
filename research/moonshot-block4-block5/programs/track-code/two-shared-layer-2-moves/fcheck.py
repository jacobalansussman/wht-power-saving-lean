"""fcheck.py -- w5-twoshared (resumed): numerical check of the facts (F1)-(F3) behind Lemma D4 / Theorem P on private
cores (numpy-free replica of sweep.core(rng, 'priv'), from uquad.py).  Usage: fcheck.py ncores seed ..."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sp
from uquad import core_priv
D = sp.D; comp = sp.comp; lag = sp.lag; lam = sp.lam
H = sp.L0 & sp.span([sp.app(D, 1 << (4 + j)) for j in range(4)])
LD = sp.span([sp.app(D, 1 << (4 + j)) for j in range(4)])
assert H == sp.L0 & LD and len(H) == 8
bad = {}; n = 0
ncores = int(sys.argv[1])
for seed in [int(a) for a in sys.argv[2:]]:
    rng = random.Random(seed)
    for c in range(ncores):
        R0, R1, tg = core_priv(rng); R = [R0, R1]; n += 1
        U = []
        for x in (0, 1):
            m = R[x]['mu'][0]; y = comp(D, m)
            checks = {
                'F1 lam(D mu D) = 2': lam(comp(D, comp(m, D))) == 2,
                'F1 d(Lag y, L_D) = 2': sp.dist(lag(y), LD) == 2,
                'F3 d(L_D, I(q)) = 2': sp.dist(LD, sp.span([sp.app(R[x]['q'], 1 << (4 + j)) for j in range(4)])) == 2,
            }
            Ux = sp.L0 & lag(comp(D, comp(m, D)))
            checks['F2 dim U = 2'] = len(Ux) == 4
            checks['F2 U in H'] = Ux <= H
            checks['F2 U in Lag(y)'] = Ux <= lag(y)
            U.append(Ux)
            for k_, v_ in checks.items():
                if not v_: bad[k_] = bad.get(k_, 0) + 1
        if U[0] == U[1]: bad['U_0 = U_1'] = bad.get('U_0 = U_1', 0) + 1
        print(f'seed {seed} core {c}: cumulative failed checks {bad if bad else "none"}', flush=True)
print(f"{n} cores; failed checks: {bad if bad else 'none'}", flush=True)
