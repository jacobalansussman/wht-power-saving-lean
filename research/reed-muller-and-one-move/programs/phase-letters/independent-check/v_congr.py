import numpy as np, itertools, random
from vb import *
res = []
for m in (2, 3, 4):
    S = sym_nonsing(m); keys = {key(s) for s in S}
    alt = [s for s in S if not s.diagonal().any()]
    G = gl(m)
    orbI = {key((A @ A.T) % 2) for A in G}
    orbs = []
    rem = set(keys)
    while rem:
        k0 = next(iter(rem)); S0 = np.array(k0).reshape(m, m)
        o = {key((A @ S0 @ A.T) % 2) for A in G}
        orbs.append((len(o), int(not S0.diagonal().any()))); rem -= o
    res.append(f"m={m}: nonsingular symmetric {len(S)}, alternating {len(alt)}, orbit of I {len(orbI)}, "
               f"orbits (size, alternating) {sorted(orbs)}, |GL|={len(G)}")
# m = 5: count nonsingular symmetric and the orbit of I by random sampling
m = 5
S = sym_nonsing(m); alt = [s for s in S if not s.diagonal().any()]
rng = np.random.default_rng(1); orbI = set(); n = 0
while len(orbI) < len(S) and n < 400000:
    A = rng.integers(0, 2, (m, m))
    if rank2(A) == m:
        orbI.add(key((A @ A.T) % 2))
    n += 1
res.append(f"m=5: nonsingular symmetric {len(S)}, alternating {len(alt)}, orbit of I reached {len(orbI)} after {n} draws")
# m = 6: nonsingular alternating
m = 6; P = [(i, j) for i in range(m) for j in range(i + 1, m)]; cnt = 0
for bits in itertools.product((0, 1), repeat=len(P)):
    A = np.zeros((m, m), dtype=np.int64)
    for b, (i, j) in zip(bits, P): A[i, j] = A[j, i] = b
    if rank2(A) == m: cnt += 1
res.append(f"m=6: nonsingular alternating {cnt}; |GL(6,2)|/|Sp(6,2)| = {20158709760 // 1451520}")
print("\n".join(res))
