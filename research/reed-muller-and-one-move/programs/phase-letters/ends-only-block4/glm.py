"""glm.py (rm2-solve): Gauss-Newton / Levenberg-Marquardt for the multi-target slice, with J^H J built from per-form
Gram matrices (no Jacobian is ever stored).

F(Q) = A T (I + H T)^-1 B^T = sum over time-increasing chains s_1 < .. < s_k inside D_Q of
       (-1)^(k-1) a_(s1) H_(s1 s2) .. H_(s(k-1) sk) b_(sk)^T.
Prefix vectors L_i = t_i (a_i - sum_(j<i) H_ji L_j), suffix vectors R_i = t_i (b_i - sum_(k>i) H_ik R_k); then
F = L B^T = A R^T, and
    dF/dA_(c,i) = e_c R_i^T,   dF/dB_(c,i) = L_i e_c^T,   dF/dH_(ij) = -L_i R_j^T.
So every block of J^H J is a sum over forms of products of the Gram matrices GL = L^H L, GR = R^H R and of L, R
themselves: one GEMM per block with the forms as the inner dimension.
"""
import time, sys
import numpy as np


class Model:
    def __init__(self, t, p, W, M, dtype=np.complex128):
        self.t = np.ascontiguousarray(t, dtype=float); self.p = np.ascontiguousarray(p, dtype=float)
        self.W, self.M = W, M; self.dtype = dtype
        self.iu = np.triu_indices(M, 1)
        self.nA = W * M; self.nH = len(self.iu[0]); self.nvar = 2 * self.nA + self.nH
        self.w = np.ones(W)
        self.pen = None

    # ---- packing
    def unpack(self, x):
        W, M = self.W, self.M
        A = x[:self.nA].reshape(W, M); B = x[self.nA:2 * self.nA].reshape(W, M)
        H = np.zeros((M, M), dtype=x.dtype); H[self.iu] = x[2 * self.nA:]
        return A, B, H

    def pack(self, A, B, H):
        return np.concatenate([A.ravel(), B.ravel(), H[self.iu]])

    # ---- forward
    def LR(self, A, B, H, t):
        n = t.shape[0]; W, M = self.W, self.M
        L = np.zeros((n, W, M), dtype=self.dtype); R = np.zeros((n, W, M), dtype=self.dtype)
        for i in range(M):
            v = A[:, i][None, :]
            if i: v = v - L[:, :, :i] @ H[:i, i]
            L[:, :, i] = t[:, i, None] * v
        for i in range(M - 1, -1, -1):
            v = B[:, i][None, :]
            if i < M - 1: v = v - R[:, :, i + 1:] @ H[i, i + 1:]
            R[:, :, i] = t[:, i, None] * v
        return L, R

    def Fval(self, x, idx=None, chunk=4096):
        A, B, H = self.unpack(x)
        t = self.t if idx is None else self.t[idx]
        out = []
        for s in range(0, t.shape[0], chunk):
            L, _ = self.LR_fwd(A, H, t[s:s + chunk])
            out.append(np.einsum('nai,bi->nab', L, B))
        return np.concatenate(out)

    def LR_fwd(self, A, H, t):
        n = t.shape[0]; W, M = self.W, self.M
        L = np.zeros((n, W, M), dtype=self.dtype)
        for i in range(M):
            v = A[:, i][None, :]
            if i: v = v - L[:, :, :i] @ H[:i, i]
            L[:, :, i] = t[:, i, None] * v
        return L, None

    def cost(self, x, idx=None, chunk=4096, per_array=False):
        A, B, H = self.unpack(x)
        t = self.t if idx is None else self.t[idx]; p = self.p if idx is None else self.p[idx]
        c = 0.0; pa = np.zeros(self.W)
        ar = np.arange(self.W)
        for s in range(0, t.shape[0], chunk):
            L, _ = self.LR_fwd(A, H, t[s:s + chunk])
            F = np.einsum('nai,bi->nab', L, B)
            F[:, ar, ar] -= p[s:s + chunk]
            a2 = np.abs(F) ** 2
            rows = a2.sum(axis=(0, 2))
            if per_array: pa += rows; c += rows.sum()
            else: c += (rows * self.w).sum()
        if self.pen is not None and not per_array: c += float((self.pen * np.abs(x) ** 2).sum())
        return (c, pa) if per_array else c

    # ---- normal equations
    def normal(self, x, idx=None, chunk=1024):
        W, M = self.W, self.M; nA = self.nA; iu = self.iu
        A, B, H = self.unpack(x)
        t = self.t if idx is None else self.t[idx]; p = self.p if idx is None else self.p[idx]
        ar = np.arange(W)
        HH = np.zeros((M * M, M * M), dtype=self.dtype)
        AH = np.zeros((W * M, M * M), dtype=self.dtype)     # [(c,k),(i,l)] = sum L[c,k] GR[i,l]
        BH = np.zeros((M * M, W * M), dtype=self.dtype)     # [(i,k),(c,l)] = sum GL[i,k] R[c,l]
        AB = np.zeros((W * M, W * M), dtype=self.dtype)     # [(d,i),(c,k)] = sum conj R[d,i] L[c,k]
        sGL = np.zeros((M, M), dtype=self.dtype); sGR = np.zeros((M, M), dtype=self.dtype)
        gA = np.zeros((W, M), dtype=self.dtype); gB = np.zeros((W, M), dtype=self.dtype)
        gH = np.zeros((M, M), dtype=self.dtype)
        cost = 0.0
        for s in range(0, t.shape[0], chunk):
            tt = t[s:s + chunk]; n = tt.shape[0]
            L, R = self.LR(A, B, H, tt)
            F = np.einsum('nai,bi->nab', L, B)
            F[:, ar, ar] -= p[s:s + chunk]
            r = F
            w = self.w
            cost += ((np.abs(r) ** 2).sum(axis=(0, 2)) * w).sum()
            Lc = L.conj(); Rc = R.conj()
            wL = w[None, :, None] * L
            rw = w[None, :, None] * r
            GL = np.einsum('nai,nak->nik', Lc, wL); GR = np.einsum('nai,nak->nik', Rc, R)
            GLf = GL.reshape(n, M * M); GRf = GR.reshape(n, M * M)
            Lf = wL.reshape(n, W * M); Rf = R.reshape(n, W * M)
            HH += GLf.T @ GRf
            AH += Lf.T @ GRf
            BH += GLf.T @ Rf
            AB += Rc.reshape(n, W * M).T @ Lf
            sGL += GL.sum(0); sGR += GR.sum(0)
            gA += np.einsum('nbi,ncb->ci', Rc, rw)
            gB += np.einsum('nai,nac->ci', Lc, rw)
            gH -= np.einsum('nai,nab,nbj->ij', Lc, rw, Rc, optimize=True)
        # assemble N = J^H J
        nv = self.nvar; N = np.zeros((nv, nv), dtype=self.dtype)
        I_W = np.eye(W)
        N[:nA, :nA] = np.kron(np.diag(self.w), sGR)                                  # (c,i),(d,k): d_cd GR[i,k]
        N[nA:2 * nA, nA:2 * nA] = np.kron(I_W, sGL)
        # AB: N[(c,i),(d,k)] = sum conj R[d,i] L[c,k]  = AB[(d,i),(c,k)]
        ABr = AB.reshape(W, M, W, M)                                     # d,i,c,k
        N[:nA, nA:2 * nA] = ABr.transpose(2, 1, 0, 3).reshape(nA, nA)    # c,i,d,k
        # AH: N[(c,i),(k,l)] = - sum L[c,k] GR[i,l] = -AH[(c,k),(i,l)]
        AHr = AH.reshape(W, M, M, M)                                     # c,k,i,l
        AHt = -AHr.transpose(0, 2, 1, 3)                                 # c,i,k,l
        N[:nA, 2 * nA:] = AHt.reshape(nA, M, M)[:, iu[0], iu[1]]
        # BH: N[(c,i),(k,l)] = - sum GL[i,k] R[c,l] = -BH[(i,k),(c,l)]
        BHr = BH.reshape(M, M, W, M)                                     # i,k,c,l
        BHt = -BHr.transpose(2, 0, 1, 3)                                 # c,i,k,l
        N[nA:2 * nA, 2 * nA:] = BHt.reshape(nA, M, M)[:, iu[0], iu[1]]
        # HH: N[(i,j),(k,l)] = sum GL[i,k] GR[j,l] = HH[(i,k),(j,l)]
        HHr = HH.reshape(M, M, M, M)                                     # i,k,j,l
        HHt = HHr.transpose(0, 2, 1, 3)                                  # i,j,k,l
        N[2 * nA:, 2 * nA:] = HHt[iu[0], iu[1]][:, iu[0], iu[1]]
        # Hermitian completion
        lo = np.tril_indices(nv, -1)
        N[lo] = N.T.conj()[lo]
        g = np.concatenate([gA.ravel(), gB.ravel(), gH[iu]])
        if self.pen is not None:
            N[np.diag_indices(nv)] += self.pen; g = g + self.pen * x; cost += float((self.pen * np.abs(x) ** 2).sum())
        return N, g, cost


def lm(model, x0, iters=200, idx=None, lam=1e-2, tol=1e-26, tlimit=1e9, verbose=1, log=sys.stdout, tag="",
       eval_idx=None, extrap=True):
    x = x0.copy()
    t0 = time.time()
    N, g, c = model.normal(x, idx)
    hist = [c]
    for it in range(iters):
        d = np.real(np.diag(N)).copy()
        accepted = False
        for _ in range(12):
            Ad = N + np.diag(lam * (d + 1e-9 * (d.max() + 1e-30)))
            try:
                dx = -np.linalg.solve(Ad, g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            xn = x + dx
            cn = model.cost(xn, idx)
            if np.isfinite(cn) and cn < c:
                # extrapolation along the accepted step (cheap: cost only)
                if extrap:
                    s_ = 1.0
                    while s_ < 64:
                        xe = x + 2 * s_ * dx; ce = model.cost(xe, idx)
                        if np.isfinite(ce) and ce < cn: xn, cn, s_ = xe, ce, 2 * s_
                        else: break
                x = xn; lam = max(lam / 3, 1e-12); accepted = True
                break
            lam *= 4
        if not accepted:
            if verbose: print(f"{tag} it {it} no step accepted, lam {lam:.2e}, cost {c:.6e}", file=log, flush=True)
            if lam > 1e12: break
            continue
        N, g, c = model.normal(x, idx)
        hist.append(c)
        if verbose and (it % 10 == 0 or c < 1e-10):
            print(f"{tag} it {it} cost {c:.6e} lam {lam:.1e} t {time.time() - t0:.0f}s", file=log, flush=True)
        if c < tol or time.time() - t0 > tlimit: break
    return x, c, hist


def random_start(model, rng, sa=1.0, sh=0.3):
    W, M = model.W, model.M
    def cn(*s): return (rng.standard_normal(s) + 1j * rng.standard_normal(s)) / np.sqrt(2)
    A = sa * cn(W, M); B = sa * cn(W, M); H = np.triu(sh * cn(M, M), 1)
    return model.pack(A, B, H)


def lm_weighted(model, x0, epochs=10, per=30, gamma=1.0, idx=None, tlimit=1e9, verbose=1, log=sys.stdout, tag="",
                wclip=(0.2, 5.0)):
    """anti-split heuristic: LM in epochs; after each, re-weight arrays by their unweighted error share."""
    x = x0.copy(); t0 = time.time(); hist = []
    for ep in range(epochs):
        x, c, h = lm(model, x, iters=per, idx=idx, verbose=0, tlimit=tlimit - (time.time() - t0))
        cu, pa = model.cost(x, idx, per_array=True)
        hist.append(cu)
        if verbose: print(f"{tag} epoch {ep} unweighted {cu:.6e} w {np.array2string(model.w, precision=2)} "
                          f"t {time.time() - t0:.0f}s", file=log, flush=True)
        if cu < 1e-24 or time.time() - t0 > tlimit: break
        sh = (pa + 1e-30) / (pa.mean() + 1e-30)
        w = np.clip(sh ** gamma, *wclip); model.w = w / w.mean()
    model.w = np.ones(model.W)
    x, c, h = lm(model, x, iters=per * 3, idx=idx, verbose=0, tlimit=max(60, tlimit - (time.time() - t0)))
    return x, c, hist + h
