"""n1.py: the H-free necessary condition N1 for U_m^(r) in the commutative slice (PROOF.md Lemma 5).

For every odd-trace form Q of minimum weight (support = an affine (m-2)-flat P not through 0, odd number of basis
vectors), F_Q = A_P (I + H_P)^{-1} B_P^T = I_W with W = r 2^(m-2) = |P| r  =>  det A_P det B_P = 1.
For every even-trace minimum-weight form: F_Q = 0  =>  det A_P det B_P = 0.
(The product of the two determinants does not depend on the order of the moves.)
So: find A, B in C^{W x M} (M = r(2^m-1)) with det(A_P) det(B_P) = [P odd] for all minimum-weight flats P.
"""
import itertools, sys, time, warnings
import numpy as np
warnings.filterwarnings("ignore")
from refl import quad_forms


def min_flats(m):
    vals, tr = quad_forms(m)
    w = vals[:, 1:].sum(axis=1)
    wmin = 1 << (m - 2)
    out = []
    for f in np.nonzero(w == wmin)[0]:
        out.append((np.nonzero(vals[f, 1:])[0] + 1, int(tr[f])))
    return out


def adjT(Mx):
    """cofactor matrix (d det / d M) of a square matrix (batched over first axis)."""
    n = Mx.shape[-1]
    C = np.zeros_like(Mx)
    for i in range(n):
        for j in range(n):
            minor = np.delete(np.delete(Mx, i, axis=-2), j, axis=-1)
            C[..., i, j] = (-1) ** (i + j) * np.linalg.det(minor)
    return C


def build(m, r):
    flats = min_flats(m)
    pts = list(range(1, 1 << m))
    col = {}
    k = 0
    for z in pts:
        col[z] = list(range(k, k + r)); k += r
    M = k; W = r * (1 << (m - 2))
    idx = np.array([sum((col[z] for z in P), []) for (P, _) in flats])
    phi = np.array([p for (_, p) in flats], dtype=float)
    return M, W, idx, phi


def resid(x, M, W, idx, phi, jac=False):
    A = x[:W * M].reshape(W, M); B = x[W * M:].reshape(W, M)
    AD = np.transpose(A[:, idx], (1, 0, 2)); BD = np.transpose(B[:, idx], (1, 0, 2))
    da = np.linalg.det(AD); db = np.linalg.det(BD)
    r = da * db - phi
    if not jac:
        return r
    nf = len(phi)
    J = np.zeros((nf, 2 * W * M), dtype=complex)
    CA = adjT(AD) * db[:, None, None]; CB = adjT(BD) * da[:, None, None]
    for q in range(nf):
        cols = (np.arange(W)[:, None] * M + idx[q][None, :])
        J[q, cols.ravel()] = CA[q].ravel()
        J[q, W * M + cols.ravel()] = CB[q].ravel()
    return r, J


def solve(x, M, W, idx, phi, iters=3000):
    mu = 1e-3
    r, J = resid(x, M, W, idx, phi, True); f = np.vdot(r, r).real
    for it in range(iters):
        JH = J.conj().T; JJ = JH @ J; g = JH @ r; dg = np.real(np.diag(JJ))
        while True:
            xn = x + np.linalg.solve(JJ + mu * np.diag(dg + 1e-12), -g)
            rn = resid(xn, M, W, idx, phi); fn = np.vdot(rn, rn).real
            if fn < f:
                x, f = xn, fn; mu = max(mu / 3, 1e-12); break
            mu *= 4
            if mu > 1e12:
                return x, f
        if f < 1e-24:
            return x, f
        r, J = resid(x, M, W, idx, phi, True)
    return x, f


if __name__ == "__main__":
    m = int(sys.argv[1]); r = int(sys.argv[2]); R = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    M, W, idx, phi = build(m, r)
    print(f"m={m} r={r}: W={W} M={M} min-weight flats={len(phi)} (odd {int(phi.sum())}, even {int(len(phi)-phi.sum())})", flush=True)
    rng = np.random.default_rng(7)
    for k in range(R):
        x0 = (rng.standard_normal(2 * W * M) + 1j * rng.standard_normal(2 * W * M)) / np.sqrt(2)
        t = time.time(); x, f = solve(x0, M, W, idx, phi)
        rr = resid(x, M, W, idx, phi)
        print(f"  restart {k}: |res|^2 = {f:.3e}  max|res| = {np.abs(rr).max():.2e}  t = {time.time()-t:.1f}s", flush=True)
        if f < 1e-20:
            np.save(f"n1_sol_m{m}_r{r}_{k}.npy", x)
