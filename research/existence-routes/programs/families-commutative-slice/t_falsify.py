# falsification test for Theorem 1: units that certainly exist (disjoint unions of W representations of I_m, one
# register per array, plain-type) must never be EXCLUDED by the clique test.  Random representations of sizes m..m+4.
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from clique_test import clique_test
m, ntests, seed = map(int, sys.argv[1:4])
rng = np.random.default_rng(seed)
def cls(z):
    v = [(z >> j) & 1 for j in range(m)]
    return tuple(v[j] * v[k] for j in range(m) for k in range(j, m))
I = cls((1 << m) - 1)
I = tuple(1 if j == k else 0 for j in range(m) for k in range(j, m))
pool = []
while len(pool) < 300:
    k = int(rng.integers(m, m + 5))
    zs = [int(z) for z in rng.integers(1, 1 << m, size=k)]
    s = np.zeros(len(I), dtype=int)
    for z in zs: s = (s + np.array(cls(z))) % 2
    if tuple(s.tolist()) == I: pool.append(zs)
bad = 0; tight_hits = 0
for t in range(ntests):
    W = int(rng.integers(1, 5))
    S = sum((pool[i] for i in rng.integers(0, len(pool), size=W)), [])
    rng.shuffle(S)
    v = clique_test(m, S, W, verbose=False)
    if v[0] == "EXCLUDED": bad += 1; print("FALSIFIED?", m, W, S, v)
    if v[0] == "passes" and v[1] == W: tight_hits += 1
print(f"m={m}: {ntests} existing unions tested, EXCLUDED {bad}, clique exactly W in {tight_hits}")
