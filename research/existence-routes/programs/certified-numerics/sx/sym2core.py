# sym2core.py -- the V(x)V stand-in reduced to Sym^2 V = U6 (+) U4 (decomp.py: A = End(L2) (+) End(U6) (+) End(U4), and
# on Lambda^2 the target already lies in A, so the Lambda^2 part needs no move).  Everything here is EXACT data
# (Gaussian rationals as Fractions), with float copies for the numerics.
import numpy as np
from fractions import Fraction as Fr
exec(open('sxcore.py').read())
# exact Sym^2 basis (columns of Bs, 16 x 10): e_a e_b + e_b e_a (a < b) and e_a e_a
pairs = [(a, a) for a in range(4)] + [(a, b) for a in range(4) for b in range(a + 1, 4)]
Bs = np.zeros((16, 10))
for c, (a, b) in enumerate(pairs):
    Bs[a * 4 + b, c] += 1
    if a != b: Bs[b * 4 + a, c] += 1
Bpinv = np.linalg.solve(Bs.T @ Bs, Bs.T)      # entries 1, 1/2: exact in binary
def restr(M): return Bpinv @ M @ Bs
assert all(np.allclose(M @ Bs, Bs @ restr(M)) for M in [R, Kr] + list(Pa))
As = np.array([restr(M) for M in Pa])         # A|Sym2 (spans 52 dims)
Rs, Ks = restr(R), restr(Kr)
# central idempotent onto U6: commutant of A|Sym2 is 2-dim (I and the projector)
rows = [np.kron(M, np.eye(10)) - np.kron(np.eye(10), M.T) for M in As]
C = np.concatenate(rows, 0); _, sv, vh = np.linalg.svd(C); null = vh[(sv > 1e-9 * sv[0]).sum():]
Z = null[0].reshape(10, 10) + 0.37 * null[1].reshape(10, 10)
ev = np.linalg.eigvals(Z); ev = ev[np.argsort(ev.real)]
g1, g2 = ev[:4], ev[4:]                     # cluster by sorting; the 6-fold and the 4-fold eigenvalue
if np.ptp(g1.real) + np.ptp(g1.imag) > 1e-6: g1, g2 = ev[:6], ev[6:]
l6, l4 = (g1.mean(), g2.mean()) if len(g1) == 6 else (g2.mean(), g1.mean())
E6 = (Z - l4 * np.eye(10)) / (l6 - l4)
def recq(x, den=64):
    r = Fr(float(x.real)).limit_denominator(den); i = Fr(float(x.imag)).limit_denominator(den)
    assert abs(float(r) + 1j * float(i) - x) < 1e-9, x
    return r, i
E6q = [[recq(E6[a, b]) for b in range(10)] for a in range(10)]
E6 = np.array([[float(r) + 1j * float(i) for (r, i) in row] for row in E6q])
assert np.allclose(E6 @ E6, E6) and abs(np.trace(E6) - 6) < 1e-9 and all(np.allclose(E6 @ M, M @ E6) for M in As)
# exact bases of U6 = range E6 and U4 = range (I - E6): greedy independent columns
def cols(Pm, k):
    ch = []
    for j in range(10):
        if np.linalg.matrix_rank(Pm[:, ch + [j]], tol=1e-9) == len(ch) + 1: ch.append(j)
        if len(ch) == k: break
    return Pm[:, ch]
Ch = np.concatenate([cols(E6, 6), cols(np.eye(10) - E6, 4)], 1)
Chi = np.linalg.inv(Ch)
Rp, Kp = Chi @ Rs @ Ch, Chi @ Ks @ Ch
Ap = np.array([Chi @ M @ Ch for M in As])
assert np.allclose(Ap[:, :6, 6:], 0) and np.allclose(Ap[:, 6:, :6], 0)
BLK = [(i, j) for i in range(6) for j in range(6)] + [(6 + i, 6 + j) for i in range(4) for j in range(4)]   # 52 coords
nA = len(BLK)
def Emat(k):
    M = np.zeros((10, 10), complex); M[BLK[k]] = 1; return M
Pb = np.array([Emat(k) for k in range(nA)])
class Sys2:
    def __init__(self, W, s):
        self.W, self.s = W, s; dd = 10
        self.slots = [W + j for j in range(s)] + [W + s] * W
        self.P = sum(self.slots) * nA
        self.E = []
        for l in range(W):
            X = np.zeros((dd, dd * W), complex); X[:, l * dd:(l + 1) * dd] = np.eye(dd); self.E.append(X)
        self.T = []
        for t in range(W):
            X = np.zeros((dd, dd * W), complex); X[:, t * dd:(t + 1) * dd] = Kp; self.T.append(X)
        self.Tv = np.concatenate([x.reshape(-1) for x in self.T])
        self.idx = [(ks, i, a) for ks in range(len(self.slots)) for i in range(self.slots[ks]) for a in range(nA)]
    def unpack(self, p):
        out = []; k = 0
        for n_ in self.slots: out.append(p[k:k + n_ * nA].reshape(n_, nA)); k += n_ * nA
        return out
    def forward(self, p):
        W, s = self.W, self.s
        coef = self.unpack(p); mats = [np.tensordot(c, Pb, axes=(1, 0)) for c in coef]
        srcs = list(self.E)
        for j in range(s):
            c = sum(mats[j][i] @ srcs[i] for i in range(len(srcs))); srcs.append(Rp @ c)
        outs = [sum(mats[s + t][i] @ srcs[i] for i in range(len(srcs))) for t in range(W)]
        return mats, srcs, np.concatenate([o.reshape(-1) for o in outs])
    def res(self, p): return self.forward(p)[2] - self.Tv
    def jac(self, p):
        W, s = self.W, self.s; dd = 10
        mats, srcs, _ = self.forward(p)
        cols_ = []
        for kslot in range(len(self.slots)):
            for i in range(self.slots[kslot]):
                dsrc = [np.zeros((nA, dd, dd * W), complex) for _ in range(W)]
                for j in range(s):
                    dc = np.zeros((nA, dd, dd * W), complex)
                    for ii in range(len(dsrc)):
                        if np.any(dsrc[ii]): dc += np.matmul(mats[j][ii], dsrc[ii])
                    if j == kslot: dc += np.matmul(Pb, srcs[i])
                    dsrc.append(np.matmul(Rp, dc))
                outs = []
                for t in range(W):
                    kk = s + t; do = np.zeros((nA, dd, dd * W), complex)
                    for ii in range(len(dsrc)):
                        if np.any(dsrc[ii]): do += np.matmul(mats[kk][ii], dsrc[ii])
                    if kk == kslot: do += np.matmul(Pb, srcs[i])
                    outs.append(do.reshape(nA, -1))
                cols_.append(np.concatenate(outs, 1))
        return np.concatenate(cols_, 0).T
