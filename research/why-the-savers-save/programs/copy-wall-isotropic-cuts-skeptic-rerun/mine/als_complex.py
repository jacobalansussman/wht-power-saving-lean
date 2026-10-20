"""complex-valued masked ALS (the agent's s2_als.py, but X, Y over C with complex random starts)."""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-isotropic-cuts'))
from base import *
import numpy as np
design = sys.argv[1]; R = int(sys.argv[2]); IT = int(sys.argv[3]); seed = int(sys.argv[4])
Z = np.zeros((N, N), bool)
for i in range(N):
    for k in orth[i]: Z[i, k] = True
iso = [h for h in range(1, 512) if not dot9(h, h)]
def Cmask(h): return np.array([not dot9(lab[i], h) for i in range(N)])
special = [147, 169, 202, 240, 260, 318, 349, 359]
c60 = [h for h in iso if len(C_of([h])) == 60]
c56 = [h for h in iso if len(C_of([h])) == 56 and h not in special]
c52 = [h for h in iso if len(C_of([h])) == 52]
allm = np.ones(N, bool)
D = {'free8': [allm] * 8,
     'k8special': [allm] * 7 + [Cmask(special[0])],
     'k8c60': [allm] * 7 + [Cmask(c60[0])],
     'k8c60b': [allm] * 7 + [Cmask(c60[57])],
     'k8two': [allm] * 6 + [Cmask(special[0]), Cmask(special[1])],
     'k8c60x4': [allm] * 4 + [Cmask(c60[i]) for i in (0, 20, 40, 60)],
     'k8mixed': [allm] * 4 + [Cmask(special[2]), Cmask(c60[5]), Cmask(c56[3]), Cmask(c52[7])],
     }
rm = D[design]; K = len(rm)
M = np.array(rm).T
rows = [np.where(Z[S])[0] for S in range(N)]
rng = np.random.default_rng(seed); W = 3.0
Gm = np.array([[G[P[i]][P[k]] for k in range(N)] for i in range(N)], float)
for rs in range(R):
    X = (rng.standard_normal((N, K)) + 1j * rng.standard_normal((N, K))) * M
    Y = (rng.standard_normal((N, K)) + 1j * rng.standard_normal((N, K))) * M
    for it in range(IT):
        for S in range(N):
            a = np.where(M[S])[0]
            A = np.vstack([Y[rows[S]][:, a], W * Y[S, a][None]]); b = np.zeros(len(rows[S]) + 1, complex); b[-1] = W
            X[S] = 0; X[S, a] = np.linalg.lstsq(A, b, rcond=None)[0]
        for T in range(N):
            a = np.where(M[T])[0]
            A = np.vstack([X[rows[T]][:, a], W * X[T, a][None]]); b = np.zeros(len(rows[T]) + 1, complex); b[-1] = W
            Y[T] = 0; Y[T, a] = np.linalg.lstsq(A, b, rcond=None)[0]
        B = X @ Y.T
        res = np.sqrt((np.abs(B[Z]) ** 2).sum() + (np.abs(np.diag(B) - 1) ** 2).sum())
        if res < 1e-11: break
    dg = np.diag(B)
    nat = np.abs(B * B.T - np.outer(dg, dg) * Gm ** 2 / np.outer(np.diag(Gm), np.diag(Gm))).max()
    print('%s (complex) restart %d iters %d residual %.3e rank %d naturality defect %.2e' % (design, rs, it + 1, res,
          np.linalg.matrix_rank(B, tol=1e-6), nat), flush=True)
