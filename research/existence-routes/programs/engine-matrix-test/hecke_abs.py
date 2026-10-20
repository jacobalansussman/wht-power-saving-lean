#!/usr/bin/env python3
"""hecke_abs.py: the trivially pinned component of the engine identity in spectral form, with general inner links.
For a word with M moves (all = t up to free letters), W arrays:
   psi_j(e Z e) = tau + A R_j B = kappa_j I_W  (j = 0..m),  R_j = mu_j I + L_j,
   (L_j)_{i,i'} = psi_j(h_{i i'}),  h_{i i'} in H_{<= i-i'+1} (distance filtration), psi_0(h) = 0   (i > i'),
mu_j = omega_j(1) - 1, kappa_j = omega_j(m).  j = 0 gives tau = I.  Proved lower bound: M >= max(2W, W+m-1).
This searches for solutions ("abstract fakes"; realisability of the h's by actual free steps is checked elsewhere).
usage: hecke_abs.py m W M [trials]"""
import sys, time
import numpy as np
from scipy.optimize import least_squares
from sp import spherical

def omega_all(m):
    th, om = spherical(m)
    return om  # om[j, d]

def setup(m, W, M):
    om = omega_all(m)
    mu = om[1:, 1] - 1.0; kap = om[1:, m]
    pairs = [(i, ip) for i in range(M) for ip in range(i)]
    # for each pair, basis of H_{<=D} cap ker psi_0: u_d - u_0, d = 1..min(D, m), D = i - ip + 1
    pb = []
    for (i, ip) in pairs:
        D = min(i - ip + 1, m)
        pb.append([om[1:, d] - om[1:, 0] for d in range(1, D + 1)])  # spectral values at j = 1..m
    return om, mu, kap, pairs, pb

def solve(m, W, M, trials=10, seed=0, cplx=False):
    om, mu, kap, pairs, pb = setup(m, W, M)
    nA = W * M; nB = M * W; nh = sum(len(b) for b in pb)
    n = nA + nB + nh
    rng = np.random.default_rng(seed)
    def unpack(z):
        if cplx: z = z[:n] + 1j * z[n:]
        A = z[:nA].reshape(W, M); B = z[nA:nA + nB].reshape(M, W)
        L = np.zeros((m, M, M), dtype=z.dtype)
        k = nA + nB
        for (i, ip), b in zip(pairs, pb):
            x = z[k:k + len(b)]; k += len(b)
            L[:, i, ip] = sum(xx * bb for xx, bb in zip(x, b))
        return A, B, L
    def res(z):
        A, B, L = unpack(z)
        out = []
        for j in range(m):
            R = mu[j] * np.eye(M) + L[j]
            out.append((A @ R @ B - (kap[j] - 1) * np.eye(W)).ravel())
        r = np.concatenate(out)
        return np.concatenate([r.real, r.imag]) if cplx else r
    best = None
    for t in range(trials):
        z0 = rng.normal(size=2 * n if cplx else n)
        r = least_squares(res, z0, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        nr = np.linalg.norm(r.fun)
        if best is None or nr < best[0]: best = (nr, r.x)
        if nr < 1e-12: break
    # fresh re-check of the returned point (residual recomputed from the stored numbers)
    return best[0], np.linalg.norm(res(best[1]))

if __name__ == "__main__":
    m, W, M = map(int, sys.argv[1:4]); trials = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    t0 = time.time()
    nr, nr2 = solve(m, W, M, trials)
    print("m=%d W=%d M=%d (bound max(2W, W+m-1) = %d): best residual %.2e (recheck %.2e) %s  [%.1fs]"
          % (m, W, M, max(2 * W, W + m - 1), nr, nr2, "ABSTRACT FAKE" if nr2 < 1e-10 else "none found", time.time() - t0), flush=True)
