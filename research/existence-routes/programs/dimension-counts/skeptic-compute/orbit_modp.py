# Orbit dimension of the target under free pre/post-mixing, per entry: dim(A K + K A) (U = V), mod p, two primes.
import sys, os, json, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jac import Prog
from common import rank_mod
for repn in ('VV', 'VL2'):
    for p in (2000029, 1500041):
        Pr = Prog(2, repn, 1, 1, p, 0); s = Pr.s; d = Pr.d
        AK = np.matmul(Pr.A, Pr.K) % p; KA = np.matmul(Pr.K, Pr.A) % p
        o1 = rank_mod(np.concatenate([AK.reshape(s, -1), KA.reshape(s, -1)]), p)
        print(json.dumps({'rep': repn, 'p': p, 'dimA': s, 'dim(AK+KA)': o1, 'o for W=1,2,3': [o1 * w * w for w in (1, 2, 3)]}), flush=True)
