"""modp.py -- w3-asym: exact (mod p) version of the f = 2 widened row-space test.  Matrices of G are dyadic Gaussian
rationals: scale by 2^k (k = 12) to Z[i], map i -> sqrt(-1) mod p (p = 1 mod 4).  Blocks: D rows e_xx; S rows e_xy+e_yx;
A rows e_xy-e_yx (x<y).  member_modp(G, gens, l0, p): gens = list of term lists [(l, gid)] (coefficient 1 each,
i.e. term-split) or [(l, gid, coef)] with Gaussian-integer coef; True iff for every block the target rows lie in the
row space mod p."""
import numpy as np
def sqrtm1(p):
    for g in range(2, 1000):
        if pow(g, (p - 1) // 2, p) == p - 1:
            return pow(g, (p - 1) // 4, p)
def to_p(M, p, I, SC=4096):
    re = np.round(M.real * SC).astype(np.int64); im = np.round(M.imag * SC).astype(np.int64)
    assert np.allclose(re / SC, M.real) and np.allclose(im / SC, M.imag)
    return (re % p + (im % p) * I) % p
def rank_p(A, p):
    A = A.copy() % p; R, C = A.shape; r = 0
    A = A.astype(object) if p > 3037000499 else A.astype(np.int64)
    for c in range(C):
        if r == R: break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0: continue
        k = r + nz[0]
        if k != r: A[[r, k]] = A[[k, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        rows = np.nonzero(A[:, c])[0]
        rows = rows[rows != r]
        if rows.size:
            A[rows] = (A[rows] - np.outer(A[rows, c], A[r]) % p) % p
        r += 1
    return r
def block_rows(Mp, p):
    """Mp: 256x256 matrix mod p (rho2 of something); returns dict block -> rows mod p"""
    D = Mp[[x * 16 + x for x in range(16)]]
    xs = [(x, y) for x in range(16) for y in range(x + 1, 16)]
    S = (Mp[[x * 16 + y for x, y in xs]] + Mp[[y * 16 + x for x, y in xs]]) % p
    A = (Mp[[x * 16 + y for x, y in xs]] - Mp[[y * 16 + x for x, y in xs]]) % p
    return {'D': D, 'S': S, 'A': A}
def member_modp(G, gens, l0, p, W=2, verbose=True):
    I = sqrtm1(p)
    cache = {}
    def r2(gid):
        if gid not in cache:
            Mp = to_p(G.MATS[gid], p, I)
            cache[gid] = (np.kron(Mp, Mp)) % p
        return cache[gid]
    rows = {b: [] for b in 'DSA'}
    for T in gens:
        X = [np.zeros((256, 256), dtype=np.int64) for _ in range(W)]
        for t in T:
            l, g = t[0], t[1]; c = t[2] if len(t) > 2 else 1
            X[l] = (X[l] + c * r2(g)) % p
        for b in 'DSA':
            rows[b].append(np.hstack([block_rows(X[l], p)[b] for l in range(W)]))
    Kp = r2(G.KID)
    allok = True
    for b in 'DSA':
        S = np.vstack(rows[b]) % p
        Tm = np.zeros((block_rows(Kp, p)[b].shape[0], 256 * W), dtype=np.int64)
        Tm[:, 256 * l0:256 * (l0 + 1)] = block_rows(Kp, p)[b]
        rs = rank_p(S, p); rst = rank_p(np.vstack([S, Tm]), p)
        if verbose: print(f"   block {b}: rank S {rs}, rank [S;T] {rst}", flush=True)
        allok &= (rs == rst)
    return allok
