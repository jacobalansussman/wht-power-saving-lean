# Hook certificate (V (x) Lambda^k V, Nh = {x in T}) for the 3x3 Fourier kernel F_3 (and a dense random 3x3),
# computed mod primes p = 1 mod 3 (omega a primitive cube root of unity mod p). Monomial-invariant, so it bounds the
# widest uniform F_3 game (all permutations of the 3^m column states, all diagonals).
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge
def cube_root(p):
    for g in range(2, 500):
        w = pow(g, (p - 1) // 3, p)
        if w != 1: return w
def run(m, k, p, which):
    q = 3; n = q ** m
    w = cube_root(p)
    if which == 'F3':
        F = np.array([[1, 1, 1], [1, w, w * w % p], [1, w * w % p, w]], dtype=np.int64)
    elif which == 'F3tri':  # lower-unitriangular 'zeta of the chain 0<1<2' kernel
        F = np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1]], dtype=np.int64)
    else:
        rng = np.random.default_rng(5); F = rng.integers(2, 1000, size=(3, 3)).astype(np.int64)
    I = np.eye(3, dtype=np.int64)
    def kr(ms):
        out = np.array([[1]], dtype=np.int64)
        for M in ms: out = np.kron(out, M) % p
        return out
    K1 = kr([F] + [I] * (m - 1)); T = kr([F] * m)
    subsets = list(itertools.combinations(range(n), k))
    basis = [(x, S) for x in range(n) for S in subsets]
    inN = np.array([x in S for (x, S) in basis])
    res = {}
    for name, g in (('cap', K1), ('dem', T)):
        Wg = np.kron(g, wedge(g, k, p, subsets)) % p
        res[name] = rank_mod(Wg[np.ix_(~inN, inN)], p)
    return res, int(inN.sum()), len(basis)
P1, P2 = 2147483629, 2147483587   # check p % 3 == 1 below
for p in (P1, P2): assert p % 3 == 1, p
m = int(sys.argv[1])
for k in [int(a) for a in sys.argv[2:]]:
    for which in ('F3', 'F3tri', 'rand'):
        r1, dN, dW = run(m, k, P1, which); r2, _, _ = run(m, k, P2, which)
        print(f"q=3 m={m} n={3**m} k={k} f={k+1} dimW={dW} dimN={dN} {which}: cap={r1['cap']} dem={r1['dem']} ratio={r1['dem']/r1['cap']:.4f} agree={r1==r2}", flush=True)
