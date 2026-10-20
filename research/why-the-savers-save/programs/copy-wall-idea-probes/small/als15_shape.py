"""als15_shape.py: shape of the real rank-15 decoders of the E8 bit pattern found by ALS (H2).
Is the support a union of octad blocks (the spread decoder: cost 15 x 8 = 120), or something else?"""
import sys; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from common import *
import numpy as np
Z = np.array(NONORTH); rows = [np.where(Z[S])[0] for S in range(V)]; cols = [np.where(Z[:, T])[0] for T in range(V)]
Lq = np.array(L, dtype=float)
for seed in (1, 2):
    rng = np.random.default_rng(seed); r = 15
    X = rng.standard_normal((V, r)); Y = rng.standard_normal((V, r))
    for it in range(400):
        for S in range(V):
            A = np.vstack([Y[rows[S]], 3 * Y[S]]); b = np.zeros(len(rows[S]) + 1); b[-1] = 3
            X[S] = np.linalg.lstsq(A, b, rcond=None)[0]
        for T in range(V):
            A = np.vstack([X[cols[T]], 3 * X[T]]); b = np.zeros(len(cols[T]) + 1); b[-1] = 3
            Y[T] = np.linalg.lstsq(A, b, rcond=None)[0]
    B = X @ Y.T
    res = np.sqrt((B[Z] ** 2).sum() + ((np.diag(B) - 1) ** 2).sum())
    nz = np.abs(B) > 1e-6
    offdiag_orth = int(nz[np.array(ORTH)].sum())
    # best possible star cost with this row space: greedy over sparse vectors is hard; report the read-span of
    # each row of B (a term serving one target) and of a random basis of the row space
    sp = [np.linalg.matrix_rank(Lq[np.where(nz[S])[0]]) for S in range(V)]
    print('seed %d residual %.1e  nonzero orth entries %d (spread decoder: 840)  row read-span dims: min %d max %d'
          % (seed, res, offdiag_orth, min(sp), max(sp)))
    print('   symmetric up to scaling?', np.allclose(nz, nz.T), '  row support sizes', sorted(set(nz.sum(1)))[:10])
