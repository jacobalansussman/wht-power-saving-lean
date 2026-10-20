#!/usr/bin/env python3
"""vlm.py (v-t2-block4): my own solver for the multi-target commutative slice at block 4 (exact formulation).
Derivation (mine): a program's registers are, as formal elements, the inputs plus scalar combinations of the move
increments d_t = ([z_t z_t^T] - 1) v_t, with Q-independent scalars. At the character Q ([S] -> (-1)^<Q,S>), d_t =
-2 pi_t v_t with pi_t = Q(z_t). Writing v_t = (B e)_t + sum_{u<t} G'_tu d_u and output_l = e_l + sum_t A'_lt d_t and
rescaling (G = -2 G', A = -2 A'):  u = (I - Pi G)^-1 Pi B,  A u = O_Q  (O_Q = diag of <Q, T_l>), for all 1023 Q.
Solved as real nonlinear least squares (scipy least_squares, method 'lm', analytic Jacobian); cost = sum |residual|^2
(same normalisation as lms.py).  Usage:
  vlm.py --targets I,r.r.r.r,... --dirs z1,... [--restarts k] [--seed s] [--save file.npy]
  vlm.py --eval x.npy --targets ... --dirs ...   (cost of a saved complex vector in lms.py's packing: A, B, G-tril)"""
import sys, argparse, time
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from vcore import setup, pairing
np.seterr(all='ignore')
D = setup(4); P = D['P']; Qs = D['Qs']

def parse(sp):
    if sp == 'I': return np.eye(4, dtype=np.int64)
    rs = [int(v) for v in sp.split('.')]
    return np.array([[(rs[i] >> j) & 1 for j in range(4)] for i in range(4)], dtype=np.int64)

class Prob:
    def __init__(s, targets, dirs):
        s.W = len(targets); s.M = len(dirs)
        s.Pi = P[:, [z - 1 for z in dirs]].astype(float)                 # (1023, M)
        s.O = np.stack([pairing(Qs, S) for S in targets], 1).astype(float)  # (1023, W)
        s.tl = np.array([(t, u) for t in range(s.M) for u in range(t)], dtype=int).reshape(-1, 2)
        s.n = 2 * s.W * s.M + len(s.tl)
    def split(s, z):
        W, M = s.W, s.M
        A = z[:W * M].reshape(W, M); B = z[W * M:2 * W * M].reshape(M, W)
        G = np.zeros((M, M), complex)
        if len(s.tl): G[s.tl[:, 0], s.tl[:, 1]] = z[2 * W * M:]
        return A, B, G
    def res(s, z, jac=False):
        W, M = s.W, s.M; A, B, G = s.split(z); Pi = s.Pi
        L = np.eye(M)[None] - Pi[:, :, None] * G[None]                     # (F, M, M)
        U = np.linalg.solve(L, Pi[:, :, None] * B[None])                   # (F, M, W)
        R = np.einsum('lt,ftw->flw', A, U) - s.O[:, :, None] * np.eye(W)[None]
        if not jac: return R
        # V = A L^-1 Pi : solve L^T X = A^T  -> X^T = A L^-1
        X = np.linalg.solve(L.transpose(0, 2, 1), np.broadcast_to(A.T.astype(complex), (len(Pi), M, W)))
        V = X.transpose(0, 2, 1) * Pi[:, None, :]                          # (F, W, M)
        F = len(Pi); J = np.zeros((F, W, W, s.n), complex)
        for l in range(W): J[:, l, :, l * M:(l + 1) * M] = U.transpose(0, 2, 1)
        for t in range(M):
            for w in range(W): J[:, :, w, W * M + t * W + w] = V[:, :, t]
        if len(s.tl):
            J[..., 2 * W * M:] = V[:, :, s.tl[:, 0]][:, :, None, :] * U[:, s.tl[:, 1], :].transpose(0, 2, 1)[:, None, :, :]
        return R, J.reshape(F * W * W, s.n)

def solve(pr, z0, maxfev=3000):
    def f(x):
        z = x[:pr.n] + 1j * x[pr.n:]; r = pr.res(z).ravel(); return np.concatenate([r.real, r.imag])
    def jf(x):
        z = x[:pr.n] + 1j * x[pr.n:]; _, J = pr.res(z, True)
        return np.block([[J.real, -J.imag], [J.imag, J.real]])
    x0 = np.concatenate([z0.real, z0.imag])
    sol = least_squares(f, x0, jac=jf, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=maxfev)
    z = sol.x[:pr.n] + 1j * sol.x[pr.n:]
    return z, float(np.sum(np.abs(pr.res(z)) ** 2)), sol.nfev

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--targets', required=True); ap.add_argument('--dirs', required=True)
    ap.add_argument('--restarts', type=int, default=3); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--scale', type=float, default=0.5); ap.add_argument('--eval', default='')
    ap.add_argument('--save', default=''); ap.add_argument('--check', action='store_true')
    ap.add_argument('--maxfev', type=int, default=3000)
    a = ap.parse_args()
    tg = [parse(t) for t in a.targets.split(',')]; dirs = [int(v) for v in a.dirs.split(',')]
    pr = Prob(tg, dirs)
    dele = bool(np.all(pr.Pi.sum(1) >= pr.O.sum(1)))
    print('W=%d M=%d dirs=%s deletion_ok=%s unknowns=%d' % (pr.W, pr.M, dirs, dele, pr.n), flush=True)
    rng = np.random.default_rng(a.seed)
    if a.check:
        z = (rng.standard_normal(pr.n) + 1j * rng.standard_normal(pr.n)) * 0.4
        R, J = pr.res(z, True); dz = rng.standard_normal(pr.n) + 1j * rng.standard_normal(pr.n); e = 1e-6
        fd = (pr.res(z + e * dz) - R).ravel() / e
        print('FD rel err %.2e' % (np.linalg.norm(fd - J @ dz) / np.linalg.norm(fd))); sys.exit()
    if a.eval:
        z = np.load(a.eval)
        print('eval %s: cost %.8f |x| %.4f' % (a.eval, float(np.sum(np.abs(pr.res(z)) ** 2)), np.linalg.norm(z))); sys.exit()
    best = None
    for r in range(a.restarts):
        t0 = time.time()
        z0 = (rng.standard_normal(pr.n) + 1j * rng.standard_normal(pr.n)) * a.scale
        z, c, nf = solve(pr, z0, a.maxfev)
        print('restart %d cost %.6e |x| %.3e nfev %d  %.1f s' % (r, c, np.linalg.norm(z), nf, time.time() - t0), flush=True)
        if best is None or c < best[1]: best = (z, c)
    if a.save: np.save(a.save, best[0])
