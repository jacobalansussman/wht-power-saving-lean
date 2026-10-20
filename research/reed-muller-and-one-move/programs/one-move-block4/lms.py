#!/usr/bin/env python3
"""lms.py (t2-block4): the MULTI-TARGET commutative slice at block m, exact formulation, complex LM (analytic Jacobian).
Algebra C[Sym_m(F_2)] (e3-rmbuild Corollary R1/R2): moves multiply a register by [z z^T]; free = scalar mixing;
array l must receive [T_l] y_l, T_l = h_l h_l^T a non-alternating invertible symmetric matrix.
Why it matters [argued here, t2-block4]: [h h^T] = [h][I][h^-1] in C[Sym x| GL], so a multi-target slice program is a
genuine no-phase engine unit (address maps h_l used only before and after: feed [h_l] y_l, finish with [h_l^-1]).
Characters: for every binary quadratic form Q (all 2^(m(m+1)/2) - 1 nonzero ones), with D = {t : Q(z_t) = 1},
    A_D (I - G_DD)^-1 B_D = diag_l <Q, T_l>                                                        (*)
(the scalar case of lmh.py; same derivation).  Unknowns: A (W x M), B (M x W), G strictly lower (M x M), complex.
Usage: lms.py --m 4 --targets I,h... --dirs z1,...  (targets as comma list of 'I' or 'r0.r1.r2.r3' row bitmasks of T)
"""
import argparse, time, json, sys
import numpy as np
np.seterr(all='ignore')

class Slice:
    def __init__(s, m, targets, dirs):
        s.m = m; s.W = W = len(targets); s.M = M = len(dirs); s.dirs = dirs
        idx = [(i, j) for i in range(m) for j in range(i, m)]
        nf = 1 << len(idx)
        def Qval(Q, x):
            v = 0
            for b, (i, j) in enumerate(idx):
                if (Q >> b) & 1 and (x >> i) & 1 and (x >> j) & 1: v ^= 1
            return v
        def pairQT(Q, T):   # <Q, T> = sum_{i<=j} q_ij T_ij
            v = 0
            for b, (i, j) in enumerate(idx):
                if (Q >> b) & 1 and (T[i] >> j) & 1: v ^= 1
            return v
        s.Qval = Qval; s.pairQT = pairQT
        s.Pi = np.array([[Qval(Q, z) for z in dirs] for Q in range(1, nf)], dtype=float)       # (F, M)
        s.O = np.array([[pairQT(Q, T) for T in targets] for Q in range(1, nf)], dtype=float)   # (F, W)
        s.F = nf - 1
        s.tril = [(t, u) for t in range(M) for u in range(t)]
        s.nA = W * M; s.nB = M * W; s.nG = len(s.tril); s.nu = s.nA + s.nB + s.nG
        s.ti = np.array(s.tril, dtype=int).reshape(-1, 2)
    def deletion_ok(s):
        return bool(np.all(s.Pi.sum(1) >= s.O.sum(1)))
    def unpack(s, x):
        W, M = s.W, s.M
        A = x[:s.nA].reshape(W, M); B = x[s.nA:s.nA + s.nB].reshape(M, W)
        G = np.zeros((M, M), dtype=complex)
        if s.nG: G[s.ti[:, 0], s.ti[:, 1]] = x[s.nA + s.nB:]
        return A, B, G
    def eval(s, x, jac=True):
        W, M, F = s.W, s.M, s.F
        A, B, G = s.unpack(x)
        L = np.eye(M)[None] - s.Pi[:, :, None] * G[None]
        Y = np.linalg.solve(L, s.Pi[:, :, None] * B[None])          # (F, M, W)
        R = np.einsum('lt,ftw->flw', A, Y) - np.einsum('fl,lw->flw', s.O, np.eye(W))
        if not jac: return R
        Z = np.linalg.solve(L.transpose(0, 2, 1), np.broadcast_to(A.T, (F, M, W))).transpose(0, 2, 1)   # (F, W, M)
        V = Z * s.Pi[:, None, :]
        J = np.zeros((F, W, W, s.nu), dtype=complex)
        # A[l, t]: R[l, :] += Y[t, :]
        for l in range(W):
            J[:, l, :, l * M:(l + 1) * M] = Y.transpose(0, 2, 1)
        # B[t, w]: R[:, w] += V[:, t]
        for t in range(M):
            for w in range(W):
                J[:, :, w, s.nA + t * W + w] = V[:, :, t]
        if s.nG:
            # G[t, u]: R += V[:, t] Y[u, :]
            J[..., s.nA + s.nB:] = np.einsum('fap,fpb->fabp', V[:, :, s.ti[:, 0]], Y[:, s.ti[:, 1], :])
        return R, J.reshape(F * W * W, s.nu)

def lm(S, x0, iters=400, tol=1e-26, maxtime=600, stall=1e-4, lam0=1e-2):
    x = x0.copy(); R, J = S.eval(x); r = R.ravel(); cost = float(np.vdot(r, r).real)
    lam = lam0; t0 = time.time(); hist = []
    it = 0
    for it in range(iters):
        JH = J.conj().T; JHJ = JH @ J; g = JH @ r
        d = np.real(np.diag(JHJ)).copy(); d[d < 1e-12] = 1e-12
        ok = False
        while lam < 1e12:
            try: dx = np.linalg.solve(JHJ + lam * np.diag(d), -g)
            except np.linalg.LinAlgError: lam *= 10; continue
            xn = x + dx; rn = S.eval(xn, jac=False).ravel(); cn = float(np.vdot(rn, rn).real)
            if cn < cost: lam = max(lam / 3, 1e-15); ok = True; break
            lam *= 4
        if not ok: break
        x = xn; cost = cn; R, J = S.eval(x); r = R.ravel(); hist.append(cost)
        if cost < tol or time.time() - t0 > maxtime: break
        if it > 60 and hist[-50] - cost < stall * hist[-50] and cost > 1e-8: break
    return x, cost, it

def parse_T(spec, m):
    if spec == 'I': return tuple(1 << i for i in range(m))
    return tuple(int(v) for v in spec.split('.'))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, default=4); ap.add_argument('--targets', default='I')
    ap.add_argument('--dirs', required=True); ap.add_argument('--restarts', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--iters', type=int, default=400)
    ap.add_argument('--check', action='store_true'); ap.add_argument('--scale', type=float, default=0.5)
    a = ap.parse_args()
    T = [parse_T(t, a.m) for t in a.targets.split(',')]
    dirs = [int(v) for v in a.dirs.split(',')]
    S = Slice(a.m, T, dirs)
    print('m=%d targets=%s dirs=%s deletion_ok=%s unknowns=%d' % (a.m, T, dirs, S.deletion_ok(), S.nu), flush=True)
    rng = np.random.default_rng(a.seed)
    if a.check:
        x = (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)) * 0.4
        R, J = S.eval(x); dx = rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu); eps = 1e-6
        fd = (S.eval(x + eps * dx, jac=False) - R).ravel() / eps
        print('FD rel err %.2e' % (np.linalg.norm(fd - J @ dx) / np.linalg.norm(fd))); return
    for rs in range(a.restarts):
        x0 = (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)) * a.scale
        x, cost, it = lm(S, x0, iters=a.iters)
        print('restart %d cost %.6e iters %d' % (rs, cost, it), flush=True)

if __name__ == '__main__':
    main()
