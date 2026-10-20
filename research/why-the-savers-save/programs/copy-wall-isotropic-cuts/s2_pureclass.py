"""s2_pureclass.py -- is a chain-compatible PURE line-class decoder realisable at all?  Pick isotropic h's whose C_h
cover P112, give each port one class, K_h = omega(T_h) + extra totals per class (read C_h, serve T_h), masked ALS.
usage: VENV -I s2_pureclass.py [extra] [restarts] [iters]"""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from base import *
import numpy as np
from itertools import combinations
extra = int(sys.argv[1]) if len(sys.argv) > 1 else 0
R = int(sys.argv[2]) if len(sys.argv) > 2 else 3; IT = int(sys.argv[3]) if len(sys.argv) > 3 else 600
iso = [h for h in range(1, 512) if not dot9(h, h)]
Cs = {h: set(C_of([h])) for h in iso}
best = None
for tri in combinations(sorted(iso, key=lambda h: -len(Cs[h]))[:60], 3):
    if len(Cs[tri[0]] | Cs[tri[1]] | Cs[tri[2]]) == N:
        best = tri; break
print('covering triple', best, [len(Cs[h]) for h in best])
T = {}; seen = set()
for h in best:
    T[h] = sorted(Cs[h] - seen); seen |= Cs[h]
rm, sm = [], []
for h in best:
    w = omega(T[h]); k = w + extra
    print('class h %d: |T_h| %d omega %d totals %d (full-span cost %d)' % (h, len(T[h]), w, k, 7 * k))
    for _ in range(k):
        rm.append(np.array([i in Cs[h] for i in range(N)])); sm.append(np.array([i in set(T[h]) for i in range(N)]))
K = len(rm); RM = np.array(rm).T; SM = np.array(sm).T
Z = np.zeros((N, N), bool)
for i in range(N):
    for k in orth[i]: Z[i, k] = True
rows = [np.where(Z[S])[0] for S in range(N)]
rng = np.random.default_rng(7); W = 3.0; bestres = 1e9
for rs in range(R):
    X = rng.standard_normal((N, K)) * SM; Y = rng.standard_normal((N, K)) * RM
    for it in range(IT):
        for S in range(N):
            a = np.where(SM[S])[0]
            A = np.vstack([Y[rows[S]][:, a], W * Y[S, a]]); b = np.zeros(len(rows[S]) + 1); b[-1] = W
            X[S] = 0; X[S, a] = np.linalg.lstsq(A, b, rcond=None)[0]
        for Tt in range(N):
            a = np.where(RM[Tt])[0]
            A = np.vstack([X[rows[Tt]][:, a], W * X[Tt, a]]); b = np.zeros(len(rows[Tt]) + 1); b[-1] = W
            Y[Tt] = 0; Y[Tt, a] = np.linalg.lstsq(A, b, rcond=None)[0]
        B = X @ Y.T
        res = np.sqrt((B[Z] ** 2).sum() + ((np.diag(B) - 1) ** 2).sum())
        if res < 1e-11: break
    bestres = min(bestres, res)
    print('K %d restart %d iters %d residual %.3e' % (K, rs, it + 1, res), flush=True)
print('BEST K %d residual %.3e' % (K, bestres))
