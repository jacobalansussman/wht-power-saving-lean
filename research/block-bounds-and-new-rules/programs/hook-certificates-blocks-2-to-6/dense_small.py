# Dense floating-point demand of g_L = prod_{b in L} dir(e_b) on the hook N in V (x) Lambda^k V, for small cases only
# (shares no code with the C programs).  Usage: python3 dense_small.py m k lab   (needs numpy; m = 3, or m = 4 with k = 2)
# Prints the rank of the block of G = g_L (x) Lambda^k g_L from N to its complement, with the singular-value gap.
import sys, itertools, numpy as np
m, k, lab = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); n = 1 << m
g = np.eye(n, dtype=complex)
# numpy 2.0 with Apple's Accelerate raises spurious divide/overflow warnings in small complex matmuls on M-series
# Macs; the products are right (the unitarity check below), so those warnings are silenced here only.
with np.errstate(all='ignore'):
    for b in range(m):
        if lab >> b & 1:
            X = np.zeros((n, n))
            for x in range(n): X[x ^ (1 << b), x] = 1
            g = (((1 + 1j) / 2) * np.eye(n) + ((1 - 1j) / 2) * X) @ g
    unitary = np.allclose(g.conj().T @ g, np.eye(n))
assert unitary and np.isfinite(g).all()
subs = list(itertools.combinations(range(n), k))
Lk = np.array([[np.linalg.det(g[np.ix_(S, T)]) for T in subs] for S in subs])
G = np.kron(g, Lk)                      # index (x, T) -> x * len(subs) + T
inN = np.array([x in T for x in range(n) for T in subs])
B = G[np.ix_(~inN, inN)]
sv = np.linalg.svd(B, compute_uv=False)
kept, dropped = sv[sv > 1e-8], sv[sv <= 1e-8]
print(f"m={m} k={k} lab={lab}: dim N {inN.sum()}, demand (dense rank) {len(kept)}, smallest kept sv {kept.min():.3e}, "
      f"largest dropped {dropped.max() if len(dropped) else 0:.1e}")
