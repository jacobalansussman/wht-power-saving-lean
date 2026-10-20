# solve_h.py -- does a call sequence compute H^(x)4 up to monomials for SOME gate parameters?
# Gate k = diag(alpha_k, beta_k) . H on its pair (all scalings between calls are free, so this is
# general).  W ~ H16 (real Sylvester up to monomials) iff every normalised cross ratio
# X_kl = W_kl W_00 / (W_0l W_k0) is +-1 and the +-1 matrix X is a Hadamard matrix.  We minimise
# sum |X_kl^2 - 1|^2 over complex (alpha, beta) by Levenberg-Marquardt from random starts.
# A residual ~1e-20 with X X^T = 16 I is a numerical certificate to be made exact separately;
# failure from many starts is evidence only.  Memory: tiny.  Usage: solve_h.py words.txt [starts]
import sys, numpy as np
from scipy.optimize import least_squares
H = np.array([[1, 1], [1, -1]], complex)
def W_of(calls, ab):
    W = np.eye(16, dtype=complex)
    for k, (u, v) in enumerate(calls):
        ru, rv = W[u].copy(), W[v].copy()
        W[u] = ab[2*k] * (ru + rv); W[v] = ab[2*k+1] * (ru - rv)
    return W
def resid(x, calls):
    n = len(calls); ab = x[:2*n] + 1j * x[2*n:]
    W = W_of(calls, ab)
    if np.min(np.abs(W)) < 1e-300: return np.full(450, 1e3)
    X = (W * W[0, 0]) / np.outer(W[:, 0], W[0, :])
    r = (X[1:, 1:] ** 2 - 1).ravel()
    return np.concatenate([r.real, r.imag])
def solve(calls, starts=20, seed=1):
    rng = np.random.default_rng(seed); best = None
    for s in range(starts):
        n = len(calls)
        x0 = np.concatenate([rng.normal(1, 0.7, 2*n), rng.normal(0, 0.7, 2*n)])
        try:
            res = least_squares(resid, x0, args=(calls,), method='lm', max_nfev=4000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        except Exception as e:
            continue
        c = 0.5 * np.sum(res.fun ** 2)
        if best is None or c < best[0]: best = (c, res.x)
        if c < 1e-20: break
    c, x = best
    ab = x[:2*len(calls)] + 1j * x[2*len(calls):]
    W = W_of(calls, ab); X = (W * W[0, 0]) / np.outer(W[:, 0], W[0, :])
    Xr = np.round(X.real); had = np.allclose(Xr @ Xr.T, 16 * np.eye(16)) and np.allclose(X, Xr, atol=1e-6)
    return c, had
if __name__ == '__main__':
    starts = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    for line in open(sys.argv[1]):
        line = line.strip()
        if not line or line.startswith('#'): continue
        name, seq = line.split(':', 1)
        calls = [tuple(map(int, p.split('-'))) for p in seq.split()]
        c, had = solve(calls, starts)
        print(f"{name.strip()} (g={len(calls)}): best residual {c:.3e}; +-1 Hadamard reached: {had}", flush=True)
