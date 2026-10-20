# small5.py -- is there a contraction of Gamma_2 with rho > 1.5 whose entropy LP is small enough to solve
# (<= 20 ground elements = ordered adjacent cell pairs + ordered cell pairs carrying sessions)?
import random, math, sys, collections
import numpy as np
import dpg, quot3
n, Dm, OPP, ADJ = quot3.setup(2)
npairs = OPP.sum(); ne = ADJ.sum() // 2
def evalp(part, k):
    M = np.zeros((k, n), dtype=np.int64); M[part, np.arange(n)] = 1
    P = M @ OPP @ M.T; C = M @ ADJ @ M.T
    used = sorted(set(part.tolist()))
    if len(used) < k: return 0, 99
    off = ~np.eye(k, dtype=bool)
    D = quot3.qmetric(np.where(off, C, 0))
    if (D >= 10 ** 6).any(): return 0, 99
    r = quot3.rho_of(P, C, D, npairs, ne)
    elems = int(((C > 0) & off).sum() + ((P > 0) & off).sum())
    return r, elems
best = {}
rng = random.Random(1)
for k in (5, 6):
    for restart in range(60):
        part = np.array([rng.randrange(k) for _ in range(n)])
        r, e = evalp(part, k); cur = r - 0.05 * max(0, e - 18)
        T = 0.05
        for it in range(4000):
            v = rng.randrange(n); old = part[v]; part[v] = rng.randrange(k)
            r2, e2 = evalp(part, k); val = r2 - 0.05 * max(0, e2 - 18)
            if val >= cur or rng.random() < math.exp((val - cur) / T):
                cur = val
                if r2 > 1.5 + 1e-9:
                    key = (k, e2)
                    if key not in best or r2 > best[key][0]: best[key] = (r2, part.copy())
            else: part[v] = old
            T = max(1e-4, T * 0.999)
for key in sorted(best): print(f"cells {key[0]}, LP elements {key[1]}: best rho {best[key][0]:.5f}")
