"""refl.py (ex2-families, 14 Oct 2026): the commutative slice as a reflection problem.

Commutative slice at block m: moves t_z (z in F_2^m \ 0) generate an abelian group A; shifts t_z^2 are free;
no address maps, no phase letters.  Restricted to the characters of A that are trivial on the shifts
(binary quadratic forms Q with Q(0) = 0, lambda_z = (-1)^Q(z)), a program with ordered moves z_1..z_M and any
scratch is EXACTLY described by (PROOF.md, Lemma 1):
    A, B in C^{W x M}, H strictly upper triangular M x M (any values),
    F_Q = A_D (I + H_D)^{-1} B_D^T   with D = {i : Q(z_i) = 1} (in program order),
and the program is a unit (on these characters) iff F_Q = [tr Q odd] I_W for every Q, tr Q = sum_j Q(e_j).
This is a NECESSARY condition for a commutative-slice unit (the mu_4 characters add more equations).

Solver: complex Levenberg-Marquardt on all forms (holomorphic residuals, analytic Jacobian).
"""
import itertools, sys, time
import numpy as np


def quad_forms(m):
    """All binary quadratic forms Q(z) = sum_{j<=k} q_jk z_j z_k on F_2^m (Q(0) = 0).
    Returns (vals, tr): vals[f, z] = Q_f(z) for z = 0..2^m-1, tr[f] = sum_j Q_f(e_j) mod 2."""
    pairs = [(j, k) for j in range(m) for k in range(j, m)]
    n = 1 << m
    zb = np.array([[(z >> j) & 1 for j in range(m)] for z in range(n)], dtype=np.int64)
    mon = np.stack([zb[:, j] * zb[:, k] for (j, k) in pairs], axis=1)  # n x P
    P = len(pairs)
    coeffs = np.array(list(itertools.product([0, 1], repeat=P)), dtype=np.int64)  # F x P
    vals = (coeffs @ mon.T) % 2  # F x n
    diag_idx = [pairs.index((j, j)) for j in range(m)]
    tr = coeffs[:, diag_idx].sum(axis=1) % 2
    return vals, tr


def deletion_bound(m, S):
    """two-point / deletion test: max W allowed = min over odd-trace Q of #{i : Q(S_i) = 1}."""
    vals, tr = quad_forms(m)
    w = vals[:, S].sum(axis=1)
    return int(w[tr == 1].min())


class System:
    def __init__(self, m, S, W, forms=None):
        self.m, self.S, self.W = m, list(S), W
        M = self.M = len(S)
        vals, tr = quad_forms(m)
        sub = vals[:, S]
        Ds, phis = [], []
        for f in range(vals.shape[0]):
            D = np.nonzero(sub[f])[0]
            if len(D) == 0:
                if tr[f] == 1:
                    raise ValueError("odd form vanishing on S: infeasible by the deletion test")
                continue
            Ds.append(D); phis.append(int(tr[f]))
        if forms is not None:
            Ds = [Ds[i] for i in forms]; phis = [phis[i] for i in forms]
        self.Ds, self.phis = Ds, phis
        self.nA = W * M
        self.pairs = [(i, j) for i in range(M) for j in range(i + 1, M)]
        self.pidx = {p: k for k, p in enumerate(self.pairs)}
        self.nvar = 2 * self.nA + len(self.pairs)
        self.nres = len(Ds) * W * W
        # precompute H index arrays per form
        self.hidx = []
        for D in Ds:
            aa, bb, kk = [], [], []
            for a in range(len(D)):
                for b in range(a + 1, len(D)):
                    aa.append(a); bb.append(b); kk.append(2 * self.nA + self.pidx[(D[a], D[b])])
            self.hidx.append((np.array(aa, int), np.array(bb, int), np.array(kk, int)))

    def unpack(self, x):
        W, M = self.W, self.M
        A = x[:self.nA].reshape(W, M)
        B = x[self.nA:2 * self.nA].reshape(W, M)
        H = np.zeros((M, M), dtype=complex)
        iu = np.triu_indices(M, 1)
        H[iu] = x[2 * self.nA:]
        return A, B, H

    def residual(self, x, jac=False):
        W, M = self.W, self.M
        A, B, H = self.unpack(x)
        r = np.zeros(self.nres, dtype=complex)
        J = np.zeros((self.nres, self.nvar), dtype=complex) if jac else None
        I = np.eye(W)
        for q, D in enumerate(self.Ds):
            d = len(D)
            X = np.linalg.inv(np.eye(d) + H[np.ix_(D, D)])
            AD, BD = A[:, D], B[:, D]
            XB = X @ BD.T          # d x W
            AX = AD @ X            # W x d
            F = AX @ BD.T
            R = F - self.phis[q] * I
            base = q * W * W
            r[base:base + W * W] = R.ravel()
            if jac:
                # d F[p, c] / d A[p, D[a]] = XB[a, c]
                for p in range(W):
                    rows = base + p * W + np.arange(W)          # (p, c) for c
                    cols = p * M + D                              # A[p, D[a]]
                    J[np.ix_(rows, cols)] += XB.T                 # [c, a]
                # d F[rr, p] / d B[p, D[a]] = AX[rr, a]
                for p in range(W):
                    rows = base + np.arange(W) * W + p            # (rr, p)
                    cols = self.nA + p * M + D
                    J[np.ix_(rows, cols)] += AX                   # [rr, a]
                aa, bb, kk = self.hidx[q]
                if len(kk):
                    # d F[rr, c] / d H[D[a], D[b]] = - AX[rr, a] XB[b, c]
                    blk = -np.einsum('ra,ac->rca', AX[:, aa], XB[bb, :])  # W x W x npairs
                    J[base:base + W * W, kk] += blk.reshape(W * W, len(kk))
        return (r, J) if jac else r

    def solve(self, x0, iters=300, tol=1e-13, verbose=False):
        x = x0.copy(); mu = 1e-3
        r, J = self.residual(x, jac=True); f = np.vdot(r, r).real
        for it in range(iters):
            JH = J.conj().T
            g = JH @ r
            JJ = JH @ J
            dg = np.real(np.diag(JJ)).copy()
            while True:
                step = np.linalg.solve(JJ + mu * np.diag(dg + 1e-12), -g)
                xn = x + step
                rn = self.residual(xn)
                fn = np.vdot(rn, rn).real
                if fn < f:
                    x, f = xn, fn; mu = max(mu / 3, 1e-12); break
                mu *= 4
                if mu > 1e12:
                    return x, f, it
            if verbose and it % 20 == 0:
                print(it, f, mu, flush=True)
            if f < tol:
                return x, f, it
            r, J = self.residual(x, jac=True)
        return x, f, iters


def random_x(sys_, rng, scale=1.0, hscale=1.0):
    n = sys_.nvar
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    x[:2 * sys_.nA] *= scale
    x[2 * sys_.nA:] *= hscale
    return x


def run(m, S, W, restarts=10, seed=0, iters=300, hscale=1.0, label=""):
    rng = np.random.default_rng(seed)
    sy = System(m, S, W)
    best = None
    t0 = time.time()
    fs = []
    for k in range(restarts):
        x0 = random_x(sy, rng, 1.0, hscale)
        x, f, it = sy.solve(x0, iters=iters)
        fs.append(f)
        if best is None or f < best[1]:
            best = (x, f)
    print(f"{label} m={m} W={W} M={len(S)} S={S} forms={len(sy.Ds)} vars={sy.nvar} "
          f"best={best[1]:.3e} all={[float('%.2e' % v) for v in sorted(fs)]} t={time.time()-t0:.1f}s", flush=True)
    return best, sy
