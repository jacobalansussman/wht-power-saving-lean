# sxcore.py -- the V(x)V stand-in of ex-dimension (block 2), W arrays in s shared paid moves, rebuilt with an EXACT basis
# of A = span rho(B): 88 group elements rho(b) = b (x) b (entries 0, +-1, +-i).  R = rho(D0), Kr = rho(D1 D0) have
# entries in Z[i]/4.  Unknowns: coefficients of each mixing matrix in that basis.  Equations: outputs = targets.
import itertools, numpy as np
m = 2; n = 4
def perm(p):
    M = np.zeros((n, n), complex)
    for x in range(n): M[p(x), x] = 1
    return M
X0 = perm(lambda x: x ^ 1); X1 = perm(lambda x: x ^ 2)
CN01 = perm(lambda x: x ^ ((x & 1) << 1)); CN10 = perm(lambda x: x ^ ((x >> 1) & 1))
S0 = np.diag([1j if (x & 1) else 1 for x in range(n)]); S1 = np.diag([1j if (x & 2) else 1 for x in range(n)])
gens = [X0, X1, CN01, CN10, S0, S1, 1j * np.eye(n)]
def key(M): return tuple(np.round(M.flatten() * 2).astype(complex).tolist())
seen = {key(np.eye(n)): np.eye(n, dtype=complex)}; frontier = [np.eye(n, dtype=complex)]
while frontier:
    nf = []
    for g in frontier:
        for h in gens:
            M = h @ g; k = key(M)
            if k not in seen: seen[k] = M; nf.append(M)
    frontier = nf
Bgrp = list(seen.values())
ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
D0 = ca * np.eye(n) + cb * X0; D1 = ca * np.eye(n) + cb * X1
K = D1 @ D0
rho = lambda g: np.kron(g, g)
d = 16
R = rho(D0); Kr = rho(K)
# exact basis of A from group elements (greedy Gram-Schmidt selection, the matrices themselves are exact)
flat = np.array([rho(b).flatten() for b in Bgrp]); chosen = []; Q = np.zeros((0, d * d), complex)
for k in range(len(Bgrp)):
    v = flat[k] - (Q.conj() @ flat[k]) @ Q if len(Q) else flat[k].copy()
    if np.linalg.norm(v) > 1e-6:
        chosen.append(k); Q = np.vstack([Q, v / np.linalg.norm(v)])
dimA = len(chosen)
Pa = np.array([rho(Bgrp[k]) for k in chosen])           # (dimA, d, d), exact
Bchosen = [Bgrp[k] for k in chosen]

class Sys:
    def __init__(self, W, s):
        self.W, self.s = W, s
        self.slots = [W + j for j in range(s)] + [W + s] * W
        self.P = sum(self.slots) * dimA
        self.E = []
        for l in range(W):
            X = np.zeros((d, d * W), complex); X[:, l * d:(l + 1) * d] = np.eye(d); self.E.append(X)
        self.T = []
        for t in range(W):
            X = np.zeros((d, d * W), complex); X[:, t * d:(t + 1) * d] = Kr; self.T.append(X)
        self.Tv = np.concatenate([x.reshape(-1) for x in self.T])
        # index map: variable k -> (slot, source, basis)
        self.idx = [(ks, i, a) for ks in range(len(self.slots)) for i in range(self.slots[ks]) for a in range(dimA)]
    def unpack(self, p):
        out = []; k = 0
        for n_ in self.slots: out.append(p[k:k + n_ * dimA].reshape(n_, dimA)); k += n_ * dimA
        return out
    def forward(self, p):
        W, s = self.W, self.s
        coef = self.unpack(p); mats = [np.tensordot(c, Pa, axes=(1, 0)) for c in coef]
        srcs = list(self.E)
        for j in range(s):
            c = sum(mats[j][i] @ srcs[i] for i in range(len(srcs))); srcs.append(R @ c)
        outs = [sum(mats[s + t][i] @ srcs[i] for i in range(len(srcs))) for t in range(W)]
        return mats, srcs, np.concatenate([o.reshape(-1) for o in outs])
    def res(self, p): return self.forward(p)[2] - self.Tv
    def jac(self, p):
        W, s = self.W, self.s
        mats, srcs, _ = self.forward(p)
        cols = []
        for kslot in range(len(self.slots)):
            for i in range(self.slots[kslot]):
                dsrc = [np.zeros((dimA, d, d * W), complex) for _ in range(W)]
                for j in range(s):
                    dc = np.zeros((dimA, d, d * W), complex)
                    if kslot < j or True:
                        for ii in range(len(dsrc)):
                            if np.any(dsrc[ii]): dc += np.matmul(mats[j][ii], dsrc[ii])
                    if j == kslot: dc += np.matmul(Pa, srcs[i])
                    dsrc.append(np.matmul(R, dc))
                outs = []
                for t in range(W):
                    kk = s + t; do = np.zeros((dimA, d, d * W), complex)
                    for ii in range(len(dsrc)):
                        if np.any(dsrc[ii]): do += np.matmul(mats[kk][ii], dsrc[ii])
                    if kk == kslot: do += np.matmul(Pa, srcs[i])
                    outs.append(do.reshape(dimA, -1))
                cols.append(np.concatenate(outs, 1))
        return np.concatenate(cols, 0).T

def gn(S, p, free, iters=60, tol=1e-13, lam=1e-6, verbose=False):
    """damped Gauss-Newton (Levenberg-Marquardt in dual form, dp = J^H (J J^H + lam I)^-1 (-r)) on the free
    coordinates (boolean mask); the others stay fixed at their values in p"""
    p = p.copy(); r = S.res(p); nr = np.linalg.norm(r)
    fi = np.where(free)[0]
    for it in range(iters):
        J = S.jac(p)[:, fi]; G = J @ J.conj().T; sc = np.trace(G).real / len(G)
        while True:
            y = np.linalg.solve(G + lam * sc * np.eye(len(G)), -r); dp = J.conj().T @ y
            q = p.copy(); q[fi] += dp; r2 = S.res(q); n2 = np.linalg.norm(r2)
            if n2 < nr: p, r, nr = q, r2, n2; lam = max(lam / 10, 1e-16); break
            lam *= 10
            if lam > 1e6: break
        if verbose: print('   it', it, '%.3e' % nr, 'lam %.0e' % lam, flush=True)
        if nr < tol or lam > 1e6: break
    return p, nr
