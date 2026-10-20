"""als_minrank.py (g-small, H2): is there a REAL matrix of rank r with the E8 bit pattern
(B[S,S] = 1, B[S,T] = 0 whenever lines S, T are not orthogonal)?  F_2 minrank is 9 (hub SAT); clique bound 8;
symmetric PSD needs rank >= theta = 15.  Alternating least squares, random restarts; prints the best residual.
usage: venv/python -I als_minrank.py r [restarts] [iters] [pattern: bit|cx]"""
import sys; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from common import *
import numpy as np
r = int(sys.argv[1]); R = int(sys.argv[2]) if len(sys.argv) > 2 else 5; IT = int(sys.argv[3]) if len(sys.argv) > 3 else 400
pat = sys.argv[4] if len(sys.argv) > 4 else 'bit'
Z = np.array(NONORTH if pat == 'bit' else ORTH)           # entries forced to zero
rows = [np.where(Z[S])[0] for S in range(V)]
cols = [np.where(Z[:, T])[0] for T in range(V)]
rng = np.random.default_rng(1)
best = 1e9
for rs in range(R):
    X = rng.standard_normal((V, r)); Y = rng.standard_normal((V, r))
    for it in range(IT):
        for S in range(V):                 # x_S: Y[rows]x = 0, y_S.x = 1  (least squares, weight on diag)
            A = np.vstack([Y[rows[S]], 3 * Y[S]]); b = np.zeros(len(rows[S]) + 1); b[-1] = 3
            X[S] = np.linalg.lstsq(A, b, rcond=None)[0]
        for T in range(V):
            A = np.vstack([X[cols[T]], 3 * X[T]]); b = np.zeros(len(cols[T]) + 1); b[-1] = 3
            Y[T] = np.linalg.lstsq(A, b, rcond=None)[0]
        B = X @ Y.T
        res = np.sqrt((B[Z] ** 2).sum() + ((np.diag(B) - 1) ** 2).sum())
    sv = np.linalg.svd(B, compute_uv=False)
    best = min(best, res)
    print('rank %d restart %d residual %.3e  min|diag| %.3f  sv tail %s' % (r, rs, res, np.abs(np.diag(B)).min(), np.round(sv[r-2:r+1], 4)))
print('BEST', r, pat, best)
