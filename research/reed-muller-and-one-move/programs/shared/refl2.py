"""refl2.py: batched version of refl.py (same system, forms grouped by support size; same residual ordering).
Optional subset of forms (screening: a subset is still a NECESSARY condition)."""
import time
import numpy as np
from refl import quad_forms


class System2:
    def __init__(self, m, S, W, subset=None, rng=None):
        self.m, self.S, self.W = m, list(S), W
        M = self.M = len(S)
        vals, tr = quad_forms(m)
        sub = vals[:, S]
        Ds, phis = [], []
        for f in range(vals.shape[0]):
            D = np.nonzero(sub[f])[0]
            if len(D) == 0:
                if tr[f] == 1:
                    raise ValueError("odd form vanishing on S")
                continue
            Ds.append(D); phis.append(int(tr[f]))
        if subset is not None and subset < len(Ds):
            rng = rng or np.random.default_rng(0)
            keep = rng.choice(len(Ds), size=subset, replace=False)
            keep.sort()
            Ds = [Ds[i] for i in keep]; phis = [phis[i] for i in keep]
        self.Ds, self.phis = Ds, np.array(phis)
        self.nA = W * M
        iu = np.triu_indices(M, 1)
        self.iu = iu
        pidx = -np.ones((M, M), dtype=int)
        pidx[iu] = np.arange(len(iu[0]))
        self.nvar = 2 * self.nA + len(iu[0])
        self.nres = len(Ds) * W * W
        # group by size
        self.groups = []
        sizes = np.array([len(D) for D in Ds])
        for d in sorted(set(sizes.tolist())):
            q = np.nonzero(sizes == d)[0]
            Dg = np.array([Ds[i] for i in q])            # nf x d
            base = q * W * W
            aa, bb = np.triu_indices(d, 1)
            hcols = 2 * self.nA + pidx[Dg[:, aa], Dg[:, bb]]   # nf x npair
            self.groups.append((d, q, Dg, base, aa, bb, hcols))

    def unpack(self, x):
        W, M = self.W, self.M
        A = x[:self.nA].reshape(W, M)
        B = x[self.nA:2 * self.nA].reshape(W, M)
        H = np.zeros((M, M), dtype=complex)
        H[self.iu] = x[2 * self.nA:]
        return A, B, H

    def residual(self, x, jac=False):
        W, M = self.W, self.M
        A, B, H = self.unpack(x)
        r = np.zeros(self.nres, dtype=complex)
        J = np.zeros((self.nres, self.nvar), dtype=complex) if jac else None
        I = np.eye(W)
        pW = np.arange(W)
        for (d, q, Dg, base, aa, bb, hcols) in self.groups:
            nf = len(q)
            HD = H[Dg[:, :, None], Dg[:, None, :]]                     # nf d d
            X = np.linalg.inv(np.eye(d)[None] + HD)
            AD = np.transpose(A[:, Dg], (1, 0, 2))                       # nf W d
            BD = np.transpose(B[:, Dg], (1, 0, 2))                       # nf W d
            XB = X @ np.transpose(BD, (0, 2, 1))                         # nf d W
            AX = AD @ X                                                  # nf W d
            F = AX @ np.transpose(BD, (0, 2, 1))                         # nf W W
            R = F - self.phis[q][:, None, None] * I[None]
            rows_all = base[:, None, None] + pW[None, :, None] * W + pW[None, None, :]   # nf W W
            r[rows_all.ravel()] = R.ravel()
            if jac:
                # A part: rows base + p W + c, cols p M + D[a], val XB[a, c]
                rows = base[:, None, None, None] + pW[None, :, None, None] * W + pW[None, None, None, :]  # nf p a c
                rows = np.broadcast_to(rows, (nf, W, d, W))
                cols = pW[None, :, None, None] * M + Dg[:, None, :, None]
                cols = np.broadcast_to(cols, (nf, W, d, W))
                vals = np.broadcast_to(XB[:, None, :, :], (nf, W, d, W))
                J[rows.ravel(), cols.ravel()] = vals.ravel()
                # B part: rows base + rr W + p, cols nA + p M + D[a], val AX[rr, a]   (index order nf rr p a)
                rows = base[:, None, None, None] + pW[None, :, None, None] * W + pW[None, None, :, None]
                rows = np.broadcast_to(rows, (nf, W, W, d))
                cols = self.nA + pW[None, None, :, None] * M + Dg[:, None, None, :]
                cols = np.broadcast_to(cols, (nf, W, W, d))
                vals = np.broadcast_to(AX[:, :, None, :], (nf, W, W, d))
                J[rows.ravel(), cols.ravel()] = vals.ravel()
                if len(aa):
                    # H part: rows base + rr W + c, cols hcols[k], val -AX[rr, a_k] XB[b_k, c]  (nf rr c k)
                    vals = -AX[:, :, None, aa] * np.transpose(XB[:, bb, :], (0, 2, 1))[:, None, :, :]
                    rows = np.broadcast_to(rows_all[:, :, :, None], vals.shape)
                    cols = np.broadcast_to(hcols[:, None, None, :], vals.shape)
                    J[rows.ravel(), cols.ravel()] = vals.ravel()
        return (r, J) if jac else r

    def solve(self, x0, iters=300, tol=1e-13, verbose=False, tmax=None):
        t0 = time.time()
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
            if f < tol or (tmax and time.time() - t0 > tmax):
                return x, f, it
            r, J = self.residual(x, jac=True)
        return x, f, iters


def random_x(sy, rng, scale=1.0, hscale=1.0):
    n = sy.nvar
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    x[:2 * sy.nA] *= scale
    x[2 * sy.nA:] *= hscale
    return x
