"""numerical minrank of the cx pattern restricted to the ROWS T = C_h (all 112 columns), h a port label (anisotropic,
omega(C_h) = 4).  If rank 4..7 fits, the counting step K_h >= omega(T_h) is far from tight for such classes."""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-isotropic-cuts'))
from base import *
import numpy as np
which = sys.argv[1]; r = int(sys.argv[2]); R = int(sys.argv[3]); IT = int(sys.argv[4])
h = lab[0] if which == 'port' else [x for x in range(1, 512) if not dot9(x, x) and len(C_of([x])) == 60][0]
rowsT = C_of([h]); m = len(rowsT); print(which, 'h', h, '|C_h|', m, 'omega', omega(rowsT), 'rank', r)
Z = np.zeros((m, N), bool); Dg = np.zeros((m, N), bool)
for a, S in enumerate(rowsT):
    for k in orth[S]: Z[a, k] = True
    Dg[a, S] = True
rng = np.random.default_rng(3); W = 3.0
best = 1e9
for rs in range(R):
    X = rng.standard_normal((m, r)); Y = rng.standard_normal((N, r))
    for it in range(IT):
        for a in range(m):
            cols = np.where(Z[a])[0]; S = rowsT[a]
            A = np.vstack([Y[cols], W * Y[S][None]]); b = np.zeros(len(cols) + 1); b[-1] = W
            X[a] = np.linalg.lstsq(A, b, rcond=None)[0]
        for T in range(N):
            rr = np.where(Z[:, T])[0]; d = np.where(Dg[:, T])[0]
            A = np.vstack([X[rr], W * X[d]]) if len(d) else X[rr]
            b = np.zeros(A.shape[0]);
            if len(d): b[-1] = W
            if A.shape[0] == 0: continue
            Y[T] = np.linalg.lstsq(A, b, rcond=None)[0]
        B = X @ Y.T
        res = np.sqrt((B[Z] ** 2).sum() + ((B[Dg] - 1) ** 2).sum())
        if res < 1e-11: break
    best = min(best, res)
    print('  restart %d iters %d residual %.3e' % (rs, it + 1, res), flush=True)
print('BEST', which, r, '%.3e' % best)
