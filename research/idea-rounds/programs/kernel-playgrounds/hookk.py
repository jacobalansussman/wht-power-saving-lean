# playgrounds/hookk.py -- the monomial-invariant hook certificate for ANY 2x2 (or qxq) kernel.
#
# Module  W = V (x) Lambda^k V,  V = C^n, n = q^m  (a GL(V)-submodule of V^(x)(k+1), i.e. the run at f = k+1 columns).
# Kept subspace  Nh = span{ e_x (x) e_T : x in T }  -- kept by EVERY monomial matrix of V (all permutations of the
# column states, all diagonals), so it is a legal Theorem-R subspace for every game whose free steps are
# column-uniform monomials (the widest uniform game), for every kernel.
#   cap(K)  = dim(Nh + theta(K_1) Nh) - dim Nh      (K_1 = kernel on direction 1; all directions conjugate)
#   dem(K)  = dim(Nh + theta(K^(x)m) Nh) - dim Nh
#   moves per array >= dem/cap.
# Ranks are computed mod a large prime p (rank_p <= rank over the true field): dem_p is a valid LOWER bound on dem;
# cap_p equals the true cap for all but finitely many p (we use two primes to screen; exact claims say so).
# Usage: python3 -I hookk.py   (runs the table; edit KERNELS below)
import sys, itertools, random
import numpy as np

def rank_mod(A, p):
    A = np.array(A, dtype=np.int64) % p
    r = 0; rows, cols = A.shape
    for c in range(cols):
        if r == rows: break
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        i = r + piv[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        nz = np.nonzero(A[:, c])[0]
        nz = nz[nz != r]
        if len(nz):
            A[nz] = (A[nz] - np.outer(A[nz, c], A[r]) % p) % p
        r += 1
    return r

def det_mod(M, p):
    M = [list(row) for row in M]; n = len(M); d = 1
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c] % p), None)
        if piv is None: return 0
        if piv != c: M[c], M[piv] = M[piv], M[c]; d = -d
        d = d * M[c][c] % p
        inv = pow(M[c][c] % p, p - 2, p)
        for i in range(c + 1, n):
            fac = M[i][c] * inv % p
            if fac:
                for j in range(c, n): M[i][j] = (M[i][j] - fac * M[c][j]) % p
    return d % p

def wedge(g, k, p, subsets):
    """Lambda^k(g) mod p as a dense matrix indexed by sorted k-subsets (rows = image)."""
    n = g.shape[0]; S = len(subsets)
    out = np.zeros((S, S), dtype=np.int64)
    for a, T in enumerate(subsets):        # row index T
        for b, U in enumerate(subsets):    # column index U
            sub = g[np.ix_(T, U)]
            if not sub.any(): continue
            out[a, b] = det_mod(sub.tolist(), p)
    return out

def kron_list(mats):
    out = np.array([[1]], dtype=object)
    for M in mats: out = np.kron(out, np.array(M, dtype=object))
    return out

def hook(K, m, k, p, q=2):
    n = q ** m
    K = np.array(K, dtype=object)
    I = np.eye(q, dtype=object)
    K1 = kron_list([K] + [I] * (m - 1))
    T = kron_list([K] * m)
    subsets = list(itertools.combinations(range(n), k))
    basis = [(x, Tt) for x in range(n) for Tt in subsets]
    inN = np.array([x in Tt for (x, Tt) in basis])
    res = {}
    for name, g in (('cap', K1), ('dem', T)):
        gp = np.array([[int(v) % p for v in row] for row in g], dtype=np.int64)
        Wg = np.kron(gp, wedge(gp, k, p, subsets)) % p     # theta(g) on V (x) Lambda^k V
        B = Wg[np.ix_(~inN, inN)]                            # image of Nh, read in W/Nh
        res[name] = rank_mod(B, p)
    res['dimN'] = int(inN.sum()); res['dimW'] = len(basis)
    return res

if __name__ == '__main__':
    P1, P2 = 2147483647, 2147483629
    # sqrt(-1) mod a prime p = 1 mod 4 for Gaussian kernels
    def sqrtm1(p):
        for g in range(2, 200):
            if pow(g, (p - 1) // 2, p) == p - 1: return pow(g, (p - 1) // 4, p)
    random.seed(1)
    KERNELS = {
        'H [[1,1],[1,-1]] (chi=-1)': [[1, 1], [1, -1]],
        'L shear [[1,0],[1,1]] (triangular)': [[1, 0], [1, 1]],
        'Fib [[1,1],[1,0]] (triangular)': [[1, 1], [1, 0]],
        'chi=2 [[1,1],[1,2]]': [[1, 1], [1, 2]],
        'chi=-2 [[1,1],[1,-2]]': [[1, 1], [1, -2]],
        'chi=3 [[1,1],[1,3]]': [[1, 1], [1, 3]],
        'random dense': [[random.randint(2, 99) for _ in range(2)] for _ in range(2)],
    }
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    ks = [int(a) for a in sys.argv[2:]] if len(sys.argv) > 2 else [4]
    for k in ks:
        for name, K in KERNELS.items():
            r1 = hook(K, m, k, P1); r2 = hook(K, m, k, P2)
            same = (r1 == r2)
            print(f"m={m} f={k+1} k={k} dimW={r1['dimW']} dimN={r1['dimN']} | {name:38s} cap={r1['cap']} dem={r1['dem']} "
                  f"ratio={r1['dem']/r1['cap'] if r1['cap'] else float('inf'):.4f}  two-primes-agree={same}", flush=True)
