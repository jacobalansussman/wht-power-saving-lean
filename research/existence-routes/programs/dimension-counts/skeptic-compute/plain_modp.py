# Jacobian rank AT the plain 1-in-2 schedule (an exact point, coefficients are group elements), mod p, two primes.
# u1 = D E, u2 = D (pi u1), O = pi u2 = K E   (pi = swap of the two address bits; pi D pi = dir e2).
import sys, os, json, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import letters, rep, perm, mm
from jac import Prog
for repn in ('VV', 'VL2'):
    for p in (2000029, 1500041):
        Pr = Prog(2, repn, 1, 2, p, 3)
        pi = rep(repn, perm(4, lambda x: ((x & 1) << 1) | (x >> 1)), p); d = Pr.d
        I = np.eye(d, dtype=np.int64); Z = np.zeros((d, d), np.int64)
        Ms = [[I], [Z, pi], [Z, Z, pi]]
        _, outs = Pr.forward(Ms)
        ok = np.array_equal(outs[0] % p, Pr.K % p)
        # pi must lie in A (it is rho of a free element): check membership
        from common import Echelon
        E = Echelon(d * d, p); E.add(Pr.A.reshape(Pr.s, -1)); inA = not E.reduce(pi.reshape(1, -1)).any()
        r = Pr.jac_rank(Ms, extra=24)
        rg = Pr.jac_rank(Pr.random_point(), extra=24)
        print(json.dumps({'rep': repn, 'p': p, 'O==K exactly': bool(ok), 'pi in A': bool(inA), 'Eeff': Pr.Eeff,
                          'rank_at_plain_modp': r, 'generic_rank_modp': rg}), flush=True)
