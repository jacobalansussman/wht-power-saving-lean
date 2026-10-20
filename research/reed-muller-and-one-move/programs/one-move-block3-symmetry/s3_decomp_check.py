#!/usr/bin/env python3
"""s3_decomp_check.py (t2-block3): numerical calibration of the per-character decomposition (PROOF.md Lemma S1).

A random Z7-symmetric program (Singer lift c acting through the Levi, c'_h = K^-1 c_h^-1 K) is run on random inputs in
the f-fold action (f = 1, 2), together with its claimed per-character programs:
    Y_{r,chi} = sum_h chi(h)^-1 c_h^-1 x_{r,h};   move  Y_r := m Y_r;   mix  Y_r' := sum_{s,k} mu_{s,k} chi(k) b_{s,k} c^k Y_s;
    inputs: regular array orbit l -> yhat_{l,chi} = sum_h chi(h)^-1 c'_h y_{l,h};  fixed array f -> 7 P_chi y_f.
Check: the Fourier combinations of the full run equal the per-character runs, register by register, at every step.
Also checks the intertwining Q'_chi K = K P_chi (the reason fixed arrays end correctly).
"""
import numpy as np
from lib3 import *

rng = np.random.default_rng(7)
c = np.load("singer_c.npy")
cp = Kinv @ np.linalg.inv(c) @ K
n = 7
z7 = np.exp(2j * np.pi / n)
# a few free letters
FREE = [X(3), Z(5), phase_letter(6), phase_letter(1, 1), perm_mat(lambda x: x ^ ((x & 1) << 1)), I8]
MOVES = [dirz(z) for z in range(1, 8)]

def kron_pow(g, f):
    out = np.array([[1.0 + 0j]])
    for _ in range(f):
        out = np.kron(out, g)
    return out

def run(f, nreg=6, nsteps=25):
    d = 8 ** f
    L = lambda g: kron_pow(g, f)
    ch = [np.linalg.matrix_power(c, h) for h in range(n)]
    cph = [np.linalg.matrix_power(cp, h) for h in range(n)]
    Lc = [L(g) for g in ch]; Lci = [L(np.linalg.inv(g)) for g in ch]; Lcp = [L(g) for g in cph]
    # arrays: one regular orbit (7 arrays y_{0,h}) and one fixed array yf
    yreg = [rng.normal(size=d) + 1j * rng.normal(size=d) for _ in range(n)]
    yf = rng.normal(size=d) + 1j * rng.normal(size=d)
    # registers: orbit 0 = regular input, orbit 1 = fixed-array input copies, others start at 0
    x = np.zeros((nreg, n, d), dtype=complex)
    for h in range(n):
        x[0, h] = Lc[h] @ Lcp[h] @ yreg[h]
        x[1, h] = Lc[h] @ Lcp[h] @ yf
    Y = np.zeros((n, nreg, d), dtype=complex)        # per character
    for k in range(n):
        Y[k, 0] = sum(z7 ** (-k * h) * Lcp[h] @ yreg[h] for h in range(n))
        Y[k, 1] = sum(z7 ** (-k * h) * Lcp[h] @ yf for h in range(n))
    def fourier(xr, k):
        return sum(z7 ** (-k * h) * Lci[h] @ xr[h] for h in range(n))
    maxerr = 0.0
    for step in range(nsteps):
        if rng.random() < 0.4:
            r = rng.integers(nreg); m = MOVES[rng.integers(7)]
            Lm = L(m)
            for h in range(n):
                x[r, h] = Lc[h] @ Lm @ Lci[h] @ x[r, h]
            for k in range(n):
                Y[k, r] = Lm @ Y[k, r]
        else:
            r = rng.integers(nreg)
            terms = [(rng.integers(nreg), rng.integers(n), rng.normal() + 1j * rng.normal(), FREE[rng.integers(len(FREE))])
                     for _ in range(rng.integers(1, 4))]
            newx = np.zeros((n, d), dtype=complex)
            for h in range(n):
                for (s, kk, mu, b) in terms:
                    newx[h] += mu * (Lc[h] @ L(b) @ Lci[h]) @ x[s, (h + kk) % n]
            x[r] = newx
            newY = np.zeros((n, d), dtype=complex)
            for k in range(n):
                for (s, kk, mu, b) in terms:
                    newY[k] += mu * z7 ** (k * kk) * (L(b) @ Lc[kk]) @ Y[k, s]
            Y[:, r] = newY
        for rr in range(nreg):
            for k in range(n):
                e = np.abs(fourier(x[rr], k) - Y[k, rr]).max() / (1 + np.abs(Y[k, rr]).max())
                maxerr = max(maxerr, e)
    # intertwining Q'_chi K = K P_chi in the f-fold action
    LK = L(K); ierr = 0.0
    for k in range(n):
        Q = sum(z7 ** (-k * h) * Lci[h] for h in range(n)) / n
        P = sum(z7 ** (-k * h) * Lcp[h] for h in range(n)) / n
        ierr = max(ierr, np.abs(Q @ LK - LK @ P).max())
    return maxerr, ierr

for f in (1, 2):
    e, ie = run(f)
    print("f = %d: max relative mismatch full-run Fourier vs per-character runs %.2e; |Q'K - KP| %.2e" % (f, e, ie))
