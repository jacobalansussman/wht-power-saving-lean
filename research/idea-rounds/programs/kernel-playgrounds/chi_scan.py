# Block-3 f=5 hook for dense 2x2 kernels K = [[1,1],[1,c]] over F_p for several c (cross-ratio chi = c), including
# the special chi with K^2 in N K N (c^2+c+2=0) and with K D K in N (always), and random c. Memory-light block build.
import sys, itertools, random
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge
def sqrt_mod(a, p):
    # Tonelli-Shanks
    a %= p
    if pow(a, (p - 1) // 2, p) != 1: return None
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
    M, c, t, R = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (M - i - 1), p); M, c, t, R = i, b * b % p, t * b * b % p, R * b % p
    return R
def kr(ms, p):
    out = np.array([[1]], dtype=np.int64)
    for M in ms: out = np.kron(out, M) % p
    return out
def hook(K, p, m=3, k=4):
    n = 2 ** m; I = np.eye(2, dtype=np.int64)
    subsets = list(itertools.combinations(range(n), k))
    res = []
    for g in (kr([K, I, I], p), kr([K, K, K], p)):
        Wk = wedge(g, k, p, subsets)
        rowsN = [(x, i) for x in range(n) for i, T in enumerate(subsets) if x not in T]
        colsN = [(y, j) for y in range(n) for j, S in enumerate(subsets) if y in S]
        rx = np.array([r[0] for r in rowsN]); ri = np.array([r[1] for r in rowsN])
        cy = np.array([c[0] for c in colsN]); cj = np.array([c[1] for c in colsN])
        res.append(rank_mod((g[np.ix_(rx, cy)] * Wk[np.ix_(ri, cj)]) % p, p))
    return res
for p in (2147483647, 2147483629, 2147483587):
    r7 = sqrt_mod(-7, p)
    cs = {}
    if r7 is not None:
        inv2 = pow(2, p - 2, p)
        cs['c^2+c+2=0 (K^2 ~ K)'] = (-1 + r7) * inv2 % p
    ri = sqrt_mod(-1, p)
    random.seed(p)
    cs['random c'] = random.randrange(2, p - 1)
    cs['c = -1 (H)'] = p - 1
    cs['c = 2'] = 2
    for name, c in cs.items():
        K = np.array([[1, 1], [1, c]], dtype=np.int64)
        cap, dem = hook(K, p)
        print(f"p={p} {name:24s} cap={cap} dem={dem} ratio={dem/cap:.4f}", flush=True)
