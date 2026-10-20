"""als.py (d-R7) -- masked alternating least squares for decoders on P112 (or a sub-unit).
Unknowns: X (n x K, serve; optional serve mask SM), R (K x n, read; read mask RM).  B = X R must have B[S,S] = 1 and
B[S,T] = 0 on orthogonal pairs (ORTH); known pairs are free.  Returns best residual and (X, R).
Residual = sqrt(sum_{orth} B^2 + sum (diag - 1)^2)."""
import numpy as np

def residual(X, R, ORTH):
    B = X @ R
    return float(np.sqrt((B[ORTH] ** 2).sum() + ((np.diag(B) - 1) ** 2).sum()))

def als(RM, ORTH, SM=None, iters=600, rng=None, w=3.0, ridge=1e-10, X0=None, R0=None, tol=1e-13, verbose=False):
    """RM: K x n bool read mask.  SM: n x K bool serve mask (None = free).  ORTH: n x n bool."""
    K, n = RM.shape
    rng = np.random.default_rng() if rng is None else rng
    O = ORTH.astype(float)
    R = rng.standard_normal((K, n)) * RM if R0 is None else R0 * RM
    X = rng.standard_normal((n, K)) if X0 is None else X0.copy()
    if SM is not None: X = X * SM
    eyeK = np.eye(K)
    RMf = RM.astype(float)
    SMf = None if SM is None else SM.astype(float)
    hist = []
    for it in range(iters):
        # X-step: for each S, min sum_{T orth S} (x.r_T)^2 + w^2 (x.r_S - 1)^2, x restricted to SM[S]
        M = np.einsum('st,kt,lt->skl', O, R, R) + w * w * np.einsum('ks,ls->skl', R, R)
        b = w * w * R.T.copy()                                  # n x K
        if SMf is not None:
            m = SMf
            M = M * m[:, :, None] * m[:, None, :] + (1 - m)[:, :, None] * eyeK[None]
            b = b * m
        M = M + ridge * eyeK[None]
        X = np.linalg.solve(M, b[:, :, None])[:, :, 0]
        # R-step: for each T, min sum_{S orth T} (x_S.r)^2 + w^2 (x_T.r - 1)^2, r restricted to RM[:, T]
        M = np.einsum('st,sk,sl->tkl', O, X, X) + w * w * np.einsum('tk,tl->tkl', X, X)
        b = w * w * X.copy()
        m = RMf.T                                               # n x K
        M = M * m[:, :, None] * m[:, None, :] + (1 - m)[:, :, None] * eyeK[None]
        b = b * m
        M = M + ridge * eyeK[None]
        R = np.linalg.solve(M, b[:, :, None])[:, :, 0].T
        if it % 20 == 19 or it == iters - 1:
            res = residual(X, R, ORTH)
            hist.append(res)
            if verbose: print(it, res, flush=True)
            if res < tol: break
            # rescale to avoid drift: balance norms per total
            nr = np.linalg.norm(R, axis=1) + 1e-300; nx = np.linalg.norm(X, axis=0) + 1e-300
            s = np.sqrt(nx / nr); R = R * s[:, None]; X = X / s[None, :]
    return residual(X, R, ORTH), X, R, hist

def als_cont(RM, ORTH, iters=1500, rng=None, w=3.0, ridge=1e-12, lam0=1e-4, lam1=1e6, hard_iters=300):
    """continuation: unmasked start, off-mask entries of R penalised with weight lam rising geometrically from lam0
    to lam1 over `iters` iterations, then `hard_iters` iterations with the hard mask.  Returns (residual, X, R, hist)."""
    K, n = RM.shape
    rng = np.random.default_rng() if rng is None else rng
    O = ORTH.astype(float)
    R = rng.standard_normal((K, n)); X = rng.standard_normal((n, K))
    off = (~RM).astype(float).T                                  # n x K
    eyeK = np.eye(K); hist = []
    for it in range(iters):
        lam = lam0 * (lam1 / lam0) ** (it / max(1, iters - 1))
        M = np.einsum('st,kt,lt->skl', O, R, R) + w * w * np.einsum('ks,ls->skl', R, R) + ridge * eyeK[None]
        X = np.linalg.solve(M, (w * w * R.T)[:, :, None])[:, :, 0]
        M = np.einsum('st,sk,sl->tkl', O, X, X) + w * w * np.einsum('tk,tl->tkl', X, X) + ridge * eyeK[None]
        M = M + lam * off[:, :, None] * eyeK[None]
        R = np.linalg.solve(M, (w * w * X)[:, :, None])[:, :, 0].T
        if it % 50 == 49:
            nr = np.linalg.norm(R, axis=1) + 1e-300; nx = np.linalg.norm(X, axis=0) + 1e-300
            s = np.sqrt(nx / nr); R = R * s[:, None]; X = X / s[None, :]
            hist.append((round(lam, 6), residual(X, R * RM, ORTH)))
    res, X, R, h2 = als(RM, ORTH, iters=hard_iters, rng=rng, w=w, ridge=ridge, X0=X, R0=R)
    return res, X, R, hist + [('hard', res)]
