# Hook block built without the full Kronecker product (memory-light): B[(x,T),(y,S)] = g[x,y] * Lambda^k(g)[T,S]
# for (x,T) not in Nh, (y,S) in Nh.  Used for n = 27 (F_3 at block 3) with small k.
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge
P = 2147483629  # = 1 mod 3
def cube_root(p):
    for g in range(2, 500):
        w = pow(g, (p - 1) // 3, p)
        if w != 1: return w
def kr(ms, p):
    out = np.array([[1]], dtype=np.int64)
    for M in ms: out = np.kron(out, M) % p
    return out
def block_rank(g, k, p):
    n = g.shape[0]
    subsets = list(itertools.combinations(range(n), k))
    Wk = wedge(g, k, p, subsets)
    rowsN = [(x, i) for x in range(n) for i, T in enumerate(subsets) if x not in T]
    colsN = [(y, j) for y in range(n) for j, S in enumerate(subsets) if y in S]
    rx = np.array([r[0] for r in rowsN]); ri = np.array([r[1] for r in rowsN])
    cy = np.array([c[0] for c in colsN]); cj = np.array([c[1] for c in colsN])
    B = (g[np.ix_(rx, cy)] * Wk[np.ix_(ri, cj)]) % p
    return rank_mod(B, p), len(colsN)
which, m, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
w = cube_root(P)
if which == 'F3':
    q = 3; K = np.array([[1, 1, 1], [1, w, w * w % P], [1, w * w % P, w]], dtype=np.int64)
elif which == 'H':
    q = 2; K = np.array([[1, 1], [1, P - 1]], dtype=np.int64)
elif which == 'L':
    q = 2; K = np.array([[1, 0], [1, 1]], dtype=np.int64)
I = np.eye(q, dtype=np.int64)
K1 = kr([K] + [I] * (m - 1), P); T = kr([K] * m, P)
cap, dN = block_rank(K1, k, P); dem, _ = block_rank(T, k, P)
print(f"{which} q={q} m={m} n={q**m} k={k} f={k+1} dimN={dN}: cap={cap} dem={dem} ratio={dem/cap:.4f}", flush=True)
