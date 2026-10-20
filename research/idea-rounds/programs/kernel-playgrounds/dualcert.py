# playgrounds/dualcert.py -- monomial-invariant Venn certificates in Lambda^a V (x) Lambda^b V*  (g acts on V* by g^{-T}).
# Y_S = span{e_A (x) e*_B : |A n B| in S} is kept by every monomial matrix (a monomial's inverse-transpose is monomial
# with the same permutation). Question: does mixing V and V* give a certificate that reaches the proved value 2 at
# block 2 for the triangular kernel (where the hook gives only 1.75)? If so, its block-3 member is calibrated.
# Usage: python3 -I dualcert.py m
import sys, itertools
from fractions import Fraction
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge
P = 2147483647

def inv_T_mod(g, p):
    # inverse transpose of an integer matrix mod p (Gauss-Jordan)
    n = g.shape[0]; A = np.concatenate([g % p, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        piv = next(i for i in range(c, n) if A[i, c] % p)
        A[[c, piv]] = A[[piv, c]]
        A[c] = A[c] * pow(int(A[c, c]), p - 2, p) % p
        for i in range(n):
            if i != c and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[c]) % p
    return A[:, n:].T.copy() % p

def kr(ms, p):
    out = np.array([[1]], dtype=np.int64)
    for M in ms: out = np.kron(out, np.array(M, dtype=np.int64)) % p
    return out

def run(K, m, a, b, p=P):
    n = 2 ** m; I = [[1, 0], [0, 1]]
    K1 = kr([K] + [I] * (m - 1), p); T = kr([K] * m, p)
    SA = list(itertools.combinations(range(n), a)); SB = list(itertools.combinations(range(n), b))
    inter = np.array([len(set(A) & set(B)) for A in SA for B in SB])
    vals = sorted(set(inter.tolist()))
    mats = {}
    for name, g in (('cap', K1), ('dem', T)):
        gi = inv_T_mod(g, p)
        mats[name] = np.kron(wedge(g, a, p, SA), wedge(gi, b, p, SB)) % p
    best = (0,)
    for r in range(1, len(vals)):
        for S in itertools.combinations(vals, r):
            inY = np.isin(inter, S)
            cap = rank_mod(mats['cap'][np.ix_(~inY, inY)], p)
            if cap == 0: continue
            dem = rank_mod(mats['dem'][np.ix_(~inY, inY)], p)
            if dem / cap > best[0]: best = (dem / cap, S, cap, dem, int(inY.sum()))
    return best

if __name__ == '__main__':
    m = int(sys.argv[1]); n = 2 ** m
    for name, K in (('H', [[1, 1], [1, -1]]), ('L', [[1, 0], [1, 1]])):
        for a in range(1, n):
            for b in range(1, n):
                if len(list(itertools.combinations(range(n), a))) * len(list(itertools.combinations(range(n), b))) > 2500: continue
                best = run(K, m, a, b)
                if best[0] > 0:
                    print(f"{name} m={m} Lambda^{a} V (x) Lambda^{b} V*: best {best[0]:.4f} S={best[1]} cap={best[2]} dem={best[3]} dimY={best[4]}", flush=True)
