#!/usr/bin/env python3
"""toylab.py (p6-refute, prop6 round, 15 Oct 2026): ABSTRACT FRAME GAMES and a two-array unit solver.

A frame game is (G, P, H, T): a finite group G, a free subgroup P, paid letters H = a union of double cosets P h P
(closed under inverses), a target T.  Frames are the right cosets P g; a paid move takes P g to P h g.  The paid
length l(g) is the frame-graph distance from P to P g.  The engine at block m is the case G = Sp(2m, 2), P = Siegel
parabolic, H = P T_z P (sl-combinatorics Lemma 1.1(d); fp-build's frame game).

A W-array unit: A(end) = F_M Y_M F_(M-1) ... Y_1 F_0 = T * I_W in M_W(C[G]) (a formal identity), where every F_k is
in M_W(C[P]) (free mixing; any final free monomial is absorbed into F_M) and Y_k = diag(h_k, 1, .., 1) with h_k a
chosen double-coset representative (other members and other slots are absorbed into the neighbouring F's).
The identity is checked irrep by irrep (C[G] = sum_rho End V_rho), i.e. as a formal identity in C[G].

Here: groups from permutation generators or from F_2 matrices, irreps by the right-commutant eigen method,
double cosets, frame distances, and an ALS solver (each F_k enters linearly) with random restarts.
"""
import itertools, sys, time, json
import numpy as np

# --------------------------------------------------------------------------------------------- groups
class Group:
    def __init__(self, gens, mul, ident, name=""):
        """gens: list of hashable elements; mul(a, b) = a*b; ident: identity."""
        self.name = name
        elems = [ident]; idx = {ident: 0}; fr = [ident]
        while fr:
            nf = []
            for a in fr:
                for g in gens:
                    c = mul(g, a)
                    if c not in idx:
                        idx[c] = len(elems); elems.append(c); nf.append(c)
            fr = nf
        self.elems, self.idx, self.n = elems, idx, len(elems)
        n = self.n
        self.mt = np.empty((n, n), dtype=np.int32)
        for i, a in enumerate(elems):
            for j, b in enumerate(elems):
                self.mt[i, j] = idx[mul(a, b)]
        self.e = 0
        self.inv = np.empty(n, dtype=np.int32)
        for i in range(n):
            self.inv[i] = int(np.where(self.mt[i] == 0)[0][0])

    def m(self, i, j): return int(self.mt[i, j])

    def subgroup(self, gen_idx):
        S = {0}; fr = [0]
        while fr:
            nf = []
            for a in fr:
                for g in gen_idx:
                    c = self.m(g, a)
                    if c not in S: S.add(c); nf.append(c)
            fr = nf
        return sorted(S)

def perm_mul(a, b):            # (a*b)(i) = a(b(i))
    return tuple(a[b[i]] for i in range(len(a)))

def perm_group(gens, name=""):
    n = len(gens[0])
    return Group([tuple(g) for g in gens], perm_mul, tuple(range(n)), name)

def cyc(n, *cycles):
    p = list(range(n))
    for c in cycles:
        for k in range(len(c)):
            p[c[k]] = c[(k + 1) % len(c)]
    return tuple(p)

def f2mat_mul(A, B):           # matrices over F_2 as tuples of row bitmasks
    out = []
    for row in A:
        acc = 0; j = 0; r = row
        while r:
            if r & 1: acc ^= B[j]
            r >>= 1; j += 1
        out.append(acc)
    return tuple(out)

def f2_group(gens, name=""):
    n = len(gens[0])
    return Group([tuple(g) for g in gens], f2mat_mul, tuple(1 << i for i in range(n)), name)

# --------------------------------------------------------------------------------------------- frames
class FrameGame:
    def __init__(self, G, P, hreps):
        """P: sorted list of element indices; hreps: list of paid representatives (indices)."""
        self.G, self.P = G, P
        Pset = set(P)
        # right cosets P g: label each element by its coset
        self.coset = -np.ones(G.n, dtype=np.int32); reps = []
        for g in range(G.n):
            if self.coset[g] < 0:
                c = len(reps); reps.append(g)
                for p in P: self.coset[G.m(p, g)] = c
        self.reps, self.ncos = reps, len(reps)
        # paid letters: union of double cosets of the hreps, closed under inverse
        H = set()
        for h in hreps:
            for p in P:
                for q in P:
                    H.add(G.m(G.m(p, h), q))
        Hinv = {int(G.inv[x]) for x in H}
        if not Hinv <= H:
            H |= Hinv
        assert not (H & Pset), "a paid letter lies in P"
        self.H = sorted(H); self.hreps = hreps
        # frame graph: P g -> P h g ; distances from the coset P (index coset[0])
        adj = [set() for _ in range(self.ncos)]
        for c, g in enumerate(reps):
            for h in self.H:
                adj[c].add(int(self.coset[G.m(h, g)]))
        self.adj = adj
        d = [-1] * self.ncos; c0 = int(self.coset[0]); d[c0] = 0; fr = [c0]
        while fr:
            nf = []
            for c in fr:
                for c2 in adj[c]:
                    if d[c2] < 0: d[c2] = d[c] + 1; nf.append(c2)
            fr = nf
        self.dist = d
        self.diam = max(d)
        self.ell = np.array([d[self.coset[g]] for g in range(G.n)])

    def double_cosets(self):
        G, P = self.G, self.P
        seen = -np.ones(G.n, dtype=np.int32); out = []
        for g in range(G.n):
            if seen[g] < 0:
                k = len(out); out.append(g)
                for p in P:
                    for q in P:
                        seen[G.m(G.m(p, g), q)] = k
        return out, seen

# --------------------------------------------------------------------------------------------- irreps
def irreps(G, seed=1, tol=1e-6):
    """All irreps (one per class) as arrays rho[g] (n x d x d), via eigenspaces of a random Hermitian element of
    the right-regular action (which commutes with the left-regular action)."""
    n = G.n; rng = np.random.default_rng(seed)
    b = rng.normal(size=n) + 1j * rng.normal(size=n)
    # Hermitian: c_g = b_g + conj(b_{g^-1})
    c = b + np.conj(b[G.inv])
    # right action R(g) e_x = e_{x g^-1}; H = sum_g c_g R(g); H[x g^-1, x] += c_g
    Hm = np.zeros((n, n), dtype=complex)
    for g in range(n):
        gi = int(G.inv[g])
        rows = G.mt[:, gi]          # x * g^-1
        Hm[rows, np.arange(n)] += c[g]
    assert np.allclose(Hm, Hm.conj().T)
    w, V = np.linalg.eigh(Hm)
    # cluster eigenvalues
    clusters = []; i = 0
    while i < n:
        j = i + 1
        while j < n and abs(w[j] - w[i]) < tol * max(1.0, abs(w[i])): j += 1
        clusters.append((i, j)); i = j
    reps = []; chars = []
    for (i, j) in clusters:
        Q = V[:, i:j]; d = j - i
        # rho(g) = Q^H L(g) Q, with L(g) e_x = e_{g x}: (L(g) Q)[g x] = Q[x]
        rho = np.empty((n, d, d), dtype=complex)
        for g in range(n):
            LQ = np.empty_like(Q); LQ[G.mt[g, :]] = Q
            rho[g] = Q.conj().T @ LQ
        ch = np.trace(rho, axis1=1, axis2=2)
        new = True
        for c2 in chars:
            if len(c2) == len(ch) and np.allclose(c2, ch, atol=1e-6):
                new = False; break
        if new:
            chars.append(ch); reps.append(rho)
    tot = sum(r.shape[1] ** 2 for r in reps)
    assert tot == n, ("irreps incomplete", tot, n, [r.shape[1] for r in reps])
    # homomorphism check on random pairs
    for _ in range(20):
        a, bb = rng.integers(n), rng.integers(n)
        for r in reps:
            assert np.allclose(r[G.m(a, bb)], r[a] @ r[bb], atol=1e-8)
    return reps

# --------------------------------------------------------------------------------------------- solver
class Solver:
    """Two-or-more-array units F_M Y_M ... Y_1 F_0 = T I_W, all moves on slot 0 (wlog, swaps are free)."""
    def __init__(self, fg, reps, W=2):
        self.fg, self.G, self.W = fg, fg.G, W
        self.reps = reps
        self.P = fg.P
        self.rP = [r[self.P] for r in reps]            # rho(p) for p in P: (|P|, d, d)
        self.wts = [np.sqrt(r.shape[1]) for r in reps]  # Plancherel weights

    def Fimg(self, Fc, k):
        # Fc: (W, W, |P|) -> block matrix (W d, W d) for irrep k
        rP = self.rP[k]; d = rP.shape[1]; W = self.W
        X = np.einsum('ijp,pab->iajb', Fc, rP)
        return X.reshape(W * d, W * d)

    def Yimg(self, h, k):
        r = self.reps[k]; d = r.shape[1]; W = self.W
        Y = np.eye(W * d, dtype=complex); Y[:d, :d] = r[h]
        return Y

    def target(self, T, k):
        r = self.reps[k]
        return np.kron(np.eye(self.W), r[T])

    def product(self, Fs, hs, k, lo=0, hi=None):
        """product F_hi Y_hi ... Y_(lo+1) F_lo   (hs[j] is the letter of move j+1, between F_j and F_(j+1))"""
        if hi is None: hi = len(Fs) - 1
        X = self.Fimg(Fs[lo], k)
        for j in range(lo + 1, hi + 1):
            X = self.Fimg(Fs[j], k) @ self.Yimg(hs[j - 1], k) @ X
        return X

    def residual(self, Fs, hs, T):
        tot = 0.0; nrm = 0.0
        for k in range(len(self.reps)):
            E = self.product(Fs, hs, k) - self.target(T, k)
            tot += self.wts[k] ** 2 * np.sum(np.abs(E) ** 2)
            nrm += self.wts[k] ** 2 * self.W * self.reps[k].shape[1]
        return np.sqrt(tot / nrm)

    def als(self, hs, T, iters=300, restarts=10, seed=0, lam=1e-9, tol=1e-12, verbose=False, init=None):
        W = self.W; nP = len(self.P); M = len(hs)
        rng = np.random.default_rng(seed)
        best = (np.inf, None)
        for rs in range(restarts):
            if init is not None and rs == 0:
                Fs = [f.copy() for f in init]
            else:
                Fs = [(rng.normal(size=(W, W, nP)) + 1j * rng.normal(size=(W, W, nP))) / np.sqrt(nP) for _ in range(M + 1)]
            res = np.inf; hist = []
            for it in range(iters):
                for kk in list(range(M + 1)) + list(range(M - 1, 0, -1)):
                    Fs[kk] = self._update(Fs, hs, T, kk, lam)
                res = self.residual(Fs, hs, T); hist.append(res)
                if res < tol: break
                if it > 40 and hist[-40] - res < 1e-6 * hist[-40]: break
            if verbose: print(f"   restart {rs}: residual {res:.3e} after {it+1} sweeps", flush=True)
            if res < best[0]: best = (res, [f.copy() for f in Fs])
            if res < tol: break
        return best

    def _update(self, Fs, hs, T, kk, lam):
        W = self.W; nP = len(self.P); M = len(hs)
        rows = []; rhs = []
        for k in range(len(self.reps)):
            d = self.reps[k].shape[1]
            # L = F_M Y ... F_(kk+1) Y_(kk+1)   ;   R = F_(kk-1) ... (identity if kk == 0)
            if kk < M:
                L = self.product(Fs, hs, k, kk + 1, M) @ self.Yimg(hs[kk], k)
            else:
                L = np.eye(W * d, dtype=complex)
            if kk > 0:
                R = self.Yimg(hs[kk - 1], k) @ self.product(Fs, hs, k, 0, kk - 1)
            else:
                R = np.eye(W * d, dtype=complex)
            rP = self.rP[k]
            Lb = L.reshape(W * d, W, d)      # L[:, i-block]
            Rb = R.reshape(W, d, W * d)      # R[j-block, :]
            # col(i,j,p) = Lb[:, i, :] @ rP[p] @ Rb[j, :, :]
            C = np.einsum('xia,pab,jby->ijpxy', Lb, rP, Rb, optimize=True)
            rows.append(self.wts[k] * C.reshape(W * W * nP, -1).T)
            rhs.append(self.wts[k] * self.target(T, k).reshape(-1))
        A = np.vstack(rows); b = np.concatenate(rhs)
        if lam > 0:
            A = np.vstack([A, np.sqrt(lam) * np.eye(A.shape[1])]); b = np.concatenate([b, np.zeros(A.shape[1])])
        f, *_ = np.linalg.lstsq(A, b, rcond=None)
        return f.reshape(W, W, nP)

# --------------------------------------------------------------------------------------------- game builders
def dihedral(n):
    r = cyc(n, list(range(n))); s = tuple((-i) % n for i in range(n))
    G = perm_group([r, s], f"D{n}")
    P = G.subgroup([G.idx[s]])
    return G, P, [G.idx[r]]

def johnson(n, k):
    """S_n with P = S_k x S_(n-k); paid = the transposition (k-1, k) double coset (one element in, one out)."""
    gens = [cyc(n, [0, 1]), cyc(n, list(range(n)))]
    G = perm_group(gens, f"J({n},{k})")
    pg = []
    if k >= 2: pg += [G.idx[cyc(n, [0, 1])], G.idx[cyc(n, list(range(k)))]]
    if n - k >= 2: pg += [G.idx[cyc(n, [k, k + 1])], G.idx[cyc(n, list(range(k, n)))]]
    P = G.subgroup(pg)
    return G, P, [G.idx[cyc(n, [k - 1, k])]]

def matching(k):
    """S_2k with P = S_2 wr S_k (stabiliser of the matching {01, 23, ...}); paid = transposition (1 2)."""
    n = 2 * k
    G = perm_group([cyc(n, [0, 1]), cyc(n, list(range(n)))], f"Match({k})")
    pg = [G.idx[cyc(n, [0, 1])]]
    if k >= 2:
        pg.append(G.idx[cyc(n, [0, 2], [1, 3])])
        pg.append(G.idx[tuple(((i + 2) % n) for i in range(n))])
    P = G.subgroup(pg)
    return G, P, [G.idx[cyc(n, [1, 2])]]

def sp4():
    """Sp(4,2) directly: vectors (a|b), bits 0,1 = a, bits 2,3 = b; T_z: (a|b) -> (a + (z.b) z | b)."""
    m = 2
    def T(z):
        rows = []
        for i in range(m): rows.append((1 << i) | ((z << m) if (z >> i) & 1 else 0))
        for i in range(m): rows.append(1 << (m + i))
        return tuple(rows)
    def ph(S):   # (a|b) -> (a | b + S a), S symmetric given as list of row masks
        rows = [1 << i for i in range(m)] + [1 << (m + i) for i in range(m)]
        # x-part rows unchanged; z-part row i gets + sum_j S[i][j] a_j
        rows = [ (1 << i) for i in range(m)] + [ (1 << (m + i)) | S[i] for i in range(m)]
        return tuple(rows)
    def addr(Mr, MiT):
        return tuple(Mr + [r << m for r in MiT])
    gens = [T(1), T(2), T(3), ph([1, 0]), ph([0, 2]), ph([2, 1])]
    G = f2_group(gens, "Sp(4,2)")
    # P = stabiliser of L0 = {(0|b)}: elements fixing span(bits 2,3) as a set: columns... use the action on vectors
    def apply(A, v):
        return sum((bin(A[i] & v).count("1") & 1) << i for i in range(len(A)))
    L0 = {0, 4, 8, 12}
    P = [i for i, A in enumerate(G.elems) if {apply(A, v) for v in L0} == L0]
    return G, P, [G.idx[T(1)]]

if __name__ == "__main__":
    pass

# --------------------------------------------------------------------------------------------- Levenberg-Marquardt
def _design(S, Fs, hs, T, kk):
    """design matrix A_kk (rows over irreps, weighted) and target vector b, so that residual = A_kk f_kk - b."""
    W = S.W; nP = len(S.P); M = len(hs)
    rows = []; rhs = []
    for k in range(len(S.reps)):
        d = S.reps[k].shape[1]
        L = S.product(Fs, hs, k, kk + 1, M) @ S.Yimg(hs[kk], k) if kk < M else np.eye(W * d, dtype=complex)
        R = S.Yimg(hs[kk - 1], k) @ S.product(Fs, hs, k, 0, kk - 1) if kk > 0 else np.eye(W * d, dtype=complex)
        Lb = L.reshape(W * d, W, d); Rb = R.reshape(W, d, W * d)
        C = np.einsum('xia,pab,jby->ijpxy', Lb, S.rP[k], Rb, optimize=True)
        rows.append(S.wts[k] * C.reshape(W * W * nP, -1).T)
        rhs.append(S.wts[k] * S.target(T, k).reshape(-1))
    return np.vstack(rows), np.concatenate(rhs)

def lm(S, hs, T, Fs, iters=200, tol=1e-13, verbose=False):
    M = len(hs); W = S.W; nP = len(S.P); blk = W * W * nP
    mu = 1e-3
    def pack(Fs): return np.concatenate([f.reshape(-1) for f in Fs])
    def unpack(x): return [x[k * blk:(k + 1) * blk].reshape(W, W, nP) for k in range(M + 1)]
    x = pack(Fs); hist = []
    A0, b = _design(S, unpack(x), hs, T, 0)
    r = A0 @ x[:blk] - b; cost = np.vdot(r, r).real
    for it in range(iters):
        Fs = unpack(x)
        J = np.hstack([_design(S, Fs, hs, T, k)[0] for k in range(M + 1)])
        g = J.conj().T @ r; Hm = J.conj().T @ J
        while True:
            dx = np.linalg.solve(Hm + mu * np.diag(np.diag(Hm).real + 1e-12), -g)
            xn = x + dx; Fn = unpack(xn)
            A0n, _ = _design(S, Fn, hs, T, 0)
            rn = A0n @ xn[:blk] - b; cn = np.vdot(rn, rn).real
            if cn < cost:
                x, r, cost = xn, rn, cn; mu = max(mu / 3, 1e-12); break
            mu *= 4
            if mu > 1e12: break
        nrm = np.sqrt(cost / np.vdot(b, b).real); hist.append(nrm)
        if verbose and it % 10 == 0: print(f"    LM it {it}: rel residual {nrm:.3e} mu {mu:.1e}", flush=True)
        if nrm < tol or mu > 1e12: break
        if it >= 30 and nrm > 1e-3 and hist[-21] - nrm < 1e-3 * hist[-21]: break   # stagnation
    return np.sqrt(cost / np.vdot(b, b).real), unpack(x)

def search(S, hs, T, restarts=10, seed=0, als_iters=60, lm_iters=150, verbose=False):
    best = (np.inf, None)
    for rs in range(restarts):
        res, Fs = S.als(hs, T, iters=als_iters, restarts=1, seed=seed * 1000 + rs)
        res2, Fs2 = lm(S, hs, T, Fs, iters=lm_iters)
        if verbose: print(f"   restart {rs}: ALS {res:.2e} -> LM {res2:.2e}", flush=True)
        if res2 < best[0]: best = (res2, Fs2)
        if res2 < 1e-12: break
    return best

def search_lm(S, hs, T, restarts=10, seed=0, lm_iters=200, verbose=False, scale=1.0):
    W = S.W; nP = len(S.P); M = len(hs); rng = np.random.default_rng(seed)
    best = (np.inf, None); allres = []
    for rs in range(restarts):
        Fs = [scale * (rng.normal(size=(W, W, nP)) + 1j * rng.normal(size=(W, W, nP))) / np.sqrt(nP) for _ in range(M + 1)]
        res, Fs2 = lm(S, hs, T, Fs, iters=lm_iters)
        allres.append(res)
        if verbose: print(f"   restart {rs}: LM {res:.2e}", flush=True)
        if res < best[0]: best = (res, Fs2)
        if res < 1e-12: break
    return best[0], best[1], allres

# --------------------------------------------------------------------------------------------- with scratch arrays
class SolverS:
    """Wt arrays (D data arrays holding the inputs at the start, Wt - D scratch arrays starting empty).
    Unit: rows 0..D-1 of F_M Y_M ... Y_1 F_0 [I_D; 0] equal T I_D (scratch rows end as anything).
    Every F_k in M_Wt(C[P]) (free mixing; copies and overwrites allowed); moves on slot 0 (wlog)."""
    def __init__(self, fg, reps, Wt=3, D=2):
        self.fg, self.G, self.Wt, self.D = fg, fg.G, Wt, D
        self.reps = reps; self.P = fg.P
        self.rP = [r[self.P] for r in reps]
        self.wts = [np.sqrt(r.shape[1]) for r in reps]
    def Fimg(self, Fc, k):
        rP = self.rP[k]; d = rP.shape[1]; W = self.Wt
        return np.einsum('ijp,pab->iajb', Fc, rP).reshape(W * d, W * d)
    def Yimg(self, h, k):
        r = self.reps[k]; d = r.shape[1]
        Y = np.eye(self.Wt * d, dtype=complex); Y[:d, :d] = r[h]; return Y
    def A0(self, k):
        d = self.reps[k].shape[1]
        A = np.zeros((self.Wt * d, self.D * d), dtype=complex); A[:self.D * d, :] = np.eye(self.D * d); return A
    def target(self, T, k):
        return np.kron(np.eye(self.D), self.reps[k][T])
    def design(self, Fs, hs, T, kk):
        Wt, D, nP, M = self.Wt, self.D, len(self.P), len(hs)
        rows = []; rhs = []
        for k in range(len(self.reps)):
            d = self.reps[k].shape[1]
            L = np.eye(Wt * d, dtype=complex)[:D * d]
            for j in range(M, kk, -1):
                L = L @ self.Fimg(Fs[j], k) @ self.Yimg(hs[j - 1], k)
            R = self.A0(k)
            for j in range(0, kk):
                R = self.Yimg(hs[j], k) @ self.Fimg(Fs[j], k) @ R
            Lb = L.reshape(D * d, Wt, d); Rb = R.reshape(Wt, d, D * d)
            C = np.einsum('xia,pab,jby->ijpxy', Lb, self.rP[k], Rb, optimize=True)
            rows.append(self.wts[k] * C.reshape(Wt * Wt * nP, -1).T)
            rhs.append(self.wts[k] * self.target(T, k).reshape(-1))
        return np.vstack(rows), np.concatenate(rhs)
    def lm(self, hs, T, Fs, iters=200, tol=1e-13):
        M = len(hs); blk = self.Wt * self.Wt * len(self.P); mu = 1e-3
        pack = lambda Fs: np.concatenate([f.reshape(-1) for f in Fs])
        unpack = lambda x: [x[k * blk:(k + 1) * blk].reshape(self.Wt, self.Wt, len(self.P)) for k in range(M + 1)]
        x = pack(Fs); A0m, b = self.design(unpack(x), hs, T, 0); hist = []
        r = A0m @ x[:blk] - b; cost = np.vdot(r, r).real; bn = np.vdot(b, b).real
        for it in range(iters):
            Fs = unpack(x)
            J = np.hstack([self.design(Fs, hs, T, k)[0] for k in range(M + 1)])
            g = J.conj().T @ r; Hm = J.conj().T @ J
            while True:
                dx = np.linalg.solve(Hm + mu * np.diag(np.diag(Hm).real + 1e-12), -g)
                xn = x + dx; A0n, _ = self.design(unpack(xn), hs, T, 0)
                rn = A0n @ xn[:blk] - b; cn = np.vdot(rn, rn).real
                if cn < cost: x, r, cost = xn, rn, cn; mu = max(mu / 3, 1e-12); break
                mu *= 4
                if mu > 1e12: break
            nrm = np.sqrt(cost / bn); hist.append(nrm)
            if nrm < tol or mu > 1e12: break
            if it >= 30 and nrm > 1e-3 and hist[-21] - nrm < 1e-3 * hist[-21]: break   # stagnation
        return np.sqrt(cost / bn), unpack(x)
    def search(self, hs, T, restarts=4, seed=0, iters=200):
        rng = np.random.default_rng(seed); best = (np.inf, None); allr = []
        for rs in range(restarts):
            Fs = [(rng.normal(size=(self.Wt, self.Wt, len(self.P))) + 1j * rng.normal(size=(self.Wt, self.Wt, len(self.P))))
                  / np.sqrt(len(self.P)) for _ in range(len(hs) + 1)]
            res, Fs2 = self.lm(hs, T, Fs, iters=iters); allr.append(res)
            if res < best[0]: best = (res, Fs2)
            if res < 1e-12: break
        return best[0], best[1], allr

def init_F(rng, W, nP, kind, eps=0.3):
    if kind == "rand":
        return (rng.normal(size=(W, W, nP)) + 1j * rng.normal(size=(W, W, nP))) / np.sqrt(nP)
    F = np.zeros((W, W, nP), dtype=complex)
    if kind == "near_id":
        for i in range(W): F[i, i, 0] = 1.0          # P[0] is the identity (index 0 in the group)
    elif kind == "mono":
        perm = rng.permutation(W)
        for i in range(W): F[i, perm[i], rng.integers(nP)] = np.exp(2j * np.pi * rng.random())
    return F + eps * (rng.normal(size=(W, W, nP)) + 1j * rng.normal(size=(W, W, nP))) / np.sqrt(nP)

def search_lm2(S, hs, T, restarts=10, seed=0, lm_iters=300, kinds=("rand", "near_id", "mono"), eps=0.3):
    W = S.W; nP = len(S.P); M = len(hs); rng = np.random.default_rng(seed)
    best = (np.inf, None); allres = []
    for rs in range(restarts):
        kind = kinds[rs % len(kinds)]
        Fs = [init_F(rng, W, nP, kind, eps) for _ in range(M + 1)]
        res, Fs2 = lm(S, hs, T, Fs, iters=lm_iters)
        allres.append((kind, res))
        if res < best[0]: best = (res, Fs2)
        if res < 1e-12: break
    return best[0], best[1], allres

def idempotent_pool(G, P, rng, n=12):
    """idempotents e_Q = mean over a random small subgroup Q of P, and 1 - e_Q, as coefficient vectors on P"""
    pos = {p: i for i, p in enumerate(P)}
    pool = []
    for _ in range(n):
        gens = [P[rng.integers(len(P))] for _ in range(rng.integers(1, 3))]
        Q = G.subgroup(gens)
        e = np.zeros(len(P), dtype=complex)
        for q in Q: e[pos[q]] += 1.0 / len(Q)
        one = np.zeros(len(P), dtype=complex); one[pos[0]] = 1.0
        pool += [e, one - e]
    return [v for v in pool if np.linalg.norm(v) > 1e-9]

def init_idem(rng, W, nP, pool, P, G, eps=0.2):
    F = np.zeros((W, W, nP), dtype=complex)
    pos = {p: i for i, p in enumerate(P)}
    for i in range(W):
        for j in range(W):
            if rng.random() < 0.7:
                e = pool[rng.integers(len(pool))]
                # e * p for a random p: shift coefficients by right multiplication
                p = P[rng.integers(nP)]
                v = np.zeros(nP, dtype=complex)
                for a in range(nP):
                    if abs(e[a]) > 0: v[pos[G.m(P[a], p)]] += e[a]
                F[i, j] = (rng.normal() + 1j * rng.normal()) * v
    return F + eps * (rng.normal(size=(W, W, nP)) + 1j * rng.normal(size=(W, W, nP))) / np.sqrt(nP)

def search_idem(S, hs, T, restarts=10, seed=0, lm_iters=300, eps=0.2):
    W = S.W; nP = len(S.P); M = len(hs); rng = np.random.default_rng(seed)
    pool = idempotent_pool(S.G, S.P, rng)
    best = (np.inf, None); allres = []
    for rs in range(restarts):
        Fs = [init_idem(rng, W, nP, pool, S.P, S.G, eps) for _ in range(M + 1)]
        res, Fs2 = lm(S, hs, T, Fs, iters=lm_iters)
        allres.append(res)
        if res < best[0]: best = (res, Fs2)
        if res < 1e-12: break
    return best[0], best[1], allres

def solverS_search2(S, hs, T, restarts=8, seed=0, iters=300, kinds=("idem", "mono", "rand", "near_id"), eps=0.25):
    """SolverS search with structured starts."""
    rng = np.random.default_rng(seed); W = S.Wt; nP = len(S.P)
    pool = idempotent_pool(S.G, S.P, rng)
    best = (np.inf, None); allr = []
    for rs in range(restarts):
        kind = kinds[rs % len(kinds)]
        if kind == "idem":
            Fs = [init_idem(rng, W, nP, pool, S.P, S.G, eps) for _ in range(len(hs) + 1)]
        else:
            Fs = [init_F(rng, W, nP, kind, eps) for _ in range(len(hs) + 1)]
        res, Fs2 = S.lm(hs, T, Fs, iters=iters); allr.append((kind, res))
        if res < best[0]: best = (res, Fs2)
        if res < 1e-12: break
    return best[0], best[1], allr
