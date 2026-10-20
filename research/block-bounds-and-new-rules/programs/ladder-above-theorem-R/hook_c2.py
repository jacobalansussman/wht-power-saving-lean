# pf-ladder: price of a rank-2 block in the hook certificates at block 3 (test of the guess "char 0 forces c2 <= 1.8 c1").
# Module V (x) Lambda^k V (V = C^8), letter g acts as g (x) Lambda^k g; N = span{ e_x (x) e_T : x in T } (kept by every
# monomial matrix).  c_r = dim(N + B_r N) - dim N for B_1 = dir e1, B_2 = dir e1 dir e2, B_3 = K.  Ranks: floating SVD
# AND exact rank modulo the prime 998244353 after scaling the Gaussian-integer matrices (i -> a square root of -1 mod p).
import sys  # needs numpy
import itertools, numpy as np
P = 998244353
I_MOD = pow(3, (P - 1) // 4, P); assert I_MOD * I_MOD % P == P - 1
def dirm(n, z):
    M = np.zeros((n, n), complex)
    for x in range(n): M[x, x] += (1 + 1j); M[x ^ z, x] += (1 - 1j)     # 2 * dir z
    return M
def wedge(g, k):
    n = g.shape[0]; subs = list(itertools.combinations(range(n), k)); idx = {s: i for i, s in enumerate(subs)}
    W = np.zeros((len(subs), len(subs)), complex)
    for j, T in enumerate(subs):
        for i, S in enumerate(subs):
            W[i, j] = np.linalg.det(g[np.ix_(S, T)]) if k > 0 else 1
    return W, subs
def to_mod(A):
    re = np.rint(A.real).astype(np.int64); im = np.rint(A.imag).astype(np.int64)
    assert np.allclose(A, re + 1j * im, atol=1e-6)
    return (re % P + (im % P) * I_MOD) % P
def rank_mod(A):
    A = A.copy() % P; r = 0; rows, cols = A.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), P - 2, P); A[r] = (A[r] * inv) % P
        nz = np.nonzero(A[:, c])[0]
        for i in nz:
            if i != r: A[i] = (A[i] - A[i, c] * A[r]) % P
        r += 1
        if r == rows: break
    return r
n = 8
for k in range(1, 7):
    D1 = dirm(n, 1); D2 = dirm(n, 2); D3 = dirm(n, 4)
    mats = {1: D1, 2: D1 @ D2, 3: D1 @ D2 @ D3}
    subs = list(itertools.combinations(range(n), k))
    cols = [x * len(subs) + i for x in range(n) for i, T in enumerate(subs) if x in T]
    dimN = len(cols)
    out = {}
    for r, g in mats.items():
        Wg, _ = wedge(g, k)
        G = np.kron(g, Wg)                                   # integer Gaussian entries (scaled by 2^(r(k+1)))
        Nb = np.zeros((n * len(subs), dimN), complex)
        for j, c in enumerate(cols): Nb[c, j] = 1
        A = np.hstack([Nb, G[:, cols]])
        rf = np.linalg.matrix_rank(A, tol=1e-6 * np.abs(A).max())
        rm = rank_mod(to_mod(A).T.copy())                    # rank of the transpose = rank
        out[r] = (rf - dimN, rm - dimN)
    c1, c2, c3 = (out[r][1] for r in (1, 2, 3))
    agree = all(out[r][0] == out[r][1] for r in out)
    print(f"k={k}: dim N = {dimN}; c1, c2, c3 = {c1}, {c2}, {c3} (float and mod-p agree: {agree}); "
          f"c2/c1 = {c2 / c1:.4f}, c3/c1 = {c3 / c1:.4f}, c3 vs c1 + c2: {c3} vs {c1 + c2}")
