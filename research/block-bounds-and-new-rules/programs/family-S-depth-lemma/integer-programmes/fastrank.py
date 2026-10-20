# b4-depth4: blocked Gaussian elimination mod a prime p < 2^21, exact.
# Per panel of B columns: Gauss-Jordan on the panel (int64) finds k pivot rows; every other row r is reduced by
# C_r @ (pivot rows) with C = P[:, pivcols] @ inv(P[pivrows, pivcols]); the products are float64 BLAS calls in
# chunks of 1024 (entries < 1024 p^2 < 2^53, exact).  Returns the rank.
# (needs numpy; the original loaded a private copy of numpy here)
import numpy as np

def rank_simple(A, p, B=512):
    """Blocked elimination, simpler bookkeeping: per panel, compute an echelon basis of the panel restricted to
    the remaining rows by an explicit transform E (k x m) with E @ panel = [I-ish], then apply to all columns."""
    assert p < 2 ** 21
    A = np.array(A, dtype=np.int64) % p
    rows, cols = A.shape
    r = 0
    for c0 in range(0, cols, B):
        if r == rows: break
        c1 = min(c0 + B, cols)
        P = A[r:, c0:c1].copy()                     # (rows-r) x w
        mrows = P.shape[0]
        # Gaussian elimination on P with full row operation tracking via the identity augmentation is too big;
        # instead: find pivot rows greedily and eliminate, keeping the list of operations as dense matrices.
        piv = []                                     # (row index in P, column index in P)
        Pw = P.copy()
        used = np.zeros(mrows, dtype=bool)
        for c in range(Pw.shape[1]):
            cand = np.nonzero((Pw[:, c] != 0) & ~used)[0]
            if len(cand) == 0: continue
            i0 = int(cand[0]); used[i0] = True
            inv = pow(int(Pw[i0, c]), p - 2, p)
            Pw[i0] = (Pw[i0] * inv) % p
            f = Pw[:, c].copy(); f[i0] = 0
            nzr = np.nonzero(f)[0]
            if len(nzr): Pw[nzr] = (Pw[nzr] - (f[nzr, None] * Pw[i0][None, :]) % p) % p
            piv.append((i0, c))
        k = len(piv)
        if k == 0: continue
        # The pivot rows of P span the row space of P; rows of P are combinations P = C @ P[pivrows] restricted
        # to the panel. Solve C: P[:, pivcols] = C @ P[pivrows][:, pivcols]  (k x k invertible block).
        pr = [i for i, _ in piv]; pcs = [c for _, c in piv]
        Bk = P[np.ix_(pr, pcs)]                      # k x k, invertible mod p
        Binv = _inv_mod(Bk, p)
        C = _mm(P[:, pcs], Binv, p)                  # (rows-r) x k
        # new rows: subtract C @ (pivot rows over ALL remaining columns) from every non-pivot row
        full_piv = A[r + np.array(pr), :]            # k x cols
        rest = np.nonzero(~used)[0].astype(np.int64)
        if len(rest):
            A[r + rest, :] = (A[r + rest, :] - _mm(C[rest], full_piv, p)) % p
        # move pivot rows to the top of the active block
        order = np.concatenate([np.array(pr, dtype=np.int64), rest]) if len(rest) else np.array(pr, dtype=np.int64)
        A[r:r + mrows] = A[r + order]
        r += k
    return r

def _mm(X, Y, p, chunk=1024):
    """exact (X @ Y) mod p for int64 inputs in [0,p), p < 2^21, via float64 in chunks of the inner dimension"""
    X = np.asarray(X) % p; Y = np.asarray(Y) % p
    out = np.zeros((X.shape[0], Y.shape[1]), dtype=np.int64)
    for s in range(0, X.shape[1], chunk):
        e = min(s + chunk, X.shape[1])
        prod = X[:, s:e].astype(np.float64) @ Y[s:e].astype(np.float64)
        out = (out + (np.rint(np.fmod(prod, p))).astype(np.int64)) % p
    return out

def _inv_mod(Bk, p):
    k = Bk.shape[0]
    M = np.concatenate([Bk % p, np.eye(k, dtype=np.int64)], axis=1)
    for c in range(k):
        nz = np.nonzero(M[c:, c])[0]
        pv = c + nz[0]
        if pv != c: M[[c, pv]] = M[[pv, c]]
        M[c] = (M[c] * pow(int(M[c, c]), p - 2, p)) % p
        f = M[:, c].copy(); f[c] = 0
        nzr = np.nonzero(f)[0]
        if len(nzr): M[nzr] = (M[nzr] - (f[nzr, None] * M[c][None, :]) % p) % p
    return M[:, k:]

if __name__ == '__main__':
    import time, sys
    p = 1048573  # prime, = 1 mod 4 ? checked below
    rng = np.random.default_rng(1)
    for (rr, cc, rk) in [(50, 40, 30), (300, 200, 150), (700, 600, 600), (1200, 900, 777)]:
        X = rng.integers(0, p, (rr, rk)); Y = rng.integers(0, p, (rk, cc))
        A = _mm(X, Y, p)
        t = time.time(); r1 = rank_simple(A, p, B=128); t1 = time.time() - t
        print(rr, cc, 'expected', min(rk, rr, cc), 'rank_simple', r1, f'{t1:.2f}s')
