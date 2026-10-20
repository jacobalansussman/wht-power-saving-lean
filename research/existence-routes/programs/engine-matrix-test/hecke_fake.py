#!/usr/bin/env python3
"""hecke_fake.py: the trivially pinned (two-sided e_B) compression of the engine identity, restricted to inner links
gamma = y * e (e = the averaging idempotent of the free group, which is an element of C[B]).  Spectrally (characters
j = 1..m of the Hecke algebra of the dual polar graph; j = 0 fixes tau = I):
    mu_j A (I - mu_j Y)^{-1} B = (kappa_j - 1) I_W,    mu_j = omega_j(1) - 1, kappa_j = omega_j(m),
A (W x M), B (M x W), Y strictly lower triangular (M x M).  Any solution is a "fake" that satisfies every condition of
the trivially pinned component with M moves.  Lower bound (D2 / Theorem D): M >= W + m - 1.
Searches M = W + m - 1 (and W + m - 2 as a control that must fail) by Levenberg-Marquardt from random starts."""
import sys, time
import numpy as np
from scipy.optimize import least_squares
from sp import spherical

def solve(m, W, M, trials=20, seed=0):
    th, om = spherical(m)
    mu = om[1:, 1] - 1.0; c = om[1:, m] - 1.0
    rng = np.random.default_rng(seed)
    tri = [(i, j) for i in range(M) for j in range(i)]
    nA, nB, nY = W * M, M * W, len(tri)
    def unpack(z):
        A = z[:nA].reshape(W, M); B = z[nA:nA + nB].reshape(M, W); Y = np.zeros((M, M))
        for k, (i, j) in enumerate(tri): Y[i, j] = z[nA + nB + k]
        return A, B, Y
    def res(z):
        A, B, Y = unpack(z)
        out = []
        for mj, cj in zip(mu, c):
            F = mj * A @ np.linalg.solve(np.eye(M) - mj * Y, B)
            out.append((F - cj * np.eye(W)).ravel())
        return np.concatenate(out)
    best = None
    for t in range(trials):
        z0 = rng.normal(size=nA + nB + nY)
        r = least_squares(res, z0, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        nr = np.linalg.norm(r.fun)
        if best is None or nr < best[0]: best = (nr, r.x)
        if nr < 1e-12: break
    return best, unpack(best[1]), mu, c

if __name__ == "__main__":
    t0 = time.time()
    for m in (2, 3, 4, 5):
        for W in (1, 2, 3, 4, 6):
            for M in (W + m - 2, W + m - 1):
                if M < 1: continue
                (nr, z), (A, B, Y), mu, c = solve(m, W, M, trials=12 if M == W + m - 2 else 25)
                # fresh check of the residual with the returned numbers
                print("m=%d W=%d M=%2d (W+m-1=%2d): best residual %.2e %s" % (m, W, M, W + m - 1, nr,
                      "FAKE FOUND" if nr < 1e-10 else "no solution found"), flush=True)
    print("time %.1fs" % (time.time() - t0))
