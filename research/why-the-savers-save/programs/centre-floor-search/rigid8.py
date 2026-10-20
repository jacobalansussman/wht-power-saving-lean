#!/usr/bin/env python3
"""ctr-search / rigid8.py: numerical sanity check of the cited rigidity (hub #1212, #1815): rank-8 realisations of the
P112 pattern (M_aa = 1, M_ab = 0 on orthogonal pairs, junk free) are Gram up to diagonal gauge.  Alternating least
squares from random starts for M = Y^T X (8 x 112 each); at convergence test the gauge invariant M_ab M_ba = 1/4 on
every junk pair (true for D1 (G/4) D2 with M_aa = 1).  Evidence only, not a proof."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
from p112 import *
rng = np.random.default_rng(11)
Eb = E.copy(); target = np.eye(N)
rows = [np.nonzero(Eb[a])[0] for a in range(N)]
for trial in range(6):
    X = rng.normal(size=(8, N)); Y = rng.normal(size=(8, N))
    for it in range(3000):
        for a in range(N):                        # y_a: <y_a, x_b> = delta_ab on E[a]
            S = rows[a]; Y[:, a] = np.linalg.lstsq(X[:, S].T, target[a, S], rcond=None)[0]
        for b in range(N):
            S = rows[b]; X[:, b] = np.linalg.lstsq(Y[:, S].T, target[S, b], rcond=None)[0]
        if it % 100 == 0:
            M = Y.T @ X; err = np.max(np.abs((M - target)[Eb]))
            if err < 1e-12: break
    M = Y.T @ X; err = np.max(np.abs((M - target)[Eb]))
    inv = (M * M.T)[~Eb]
    print(f'trial {trial}: pattern error {err:.1e}; junk M_ab*M_ba in [{inv.min():.6f}, {inv.max():.6f}] (Gram: 0.25)', flush=True)
