#!/usr/bin/env python3
"""cx-medium / e7pool.py: the cost spectrum of the E7 candidate totals of decoders.py's pool (support size, hub cost =
dim span of gathered labels), and which pool directions have support inside 1 + V for a 4-dim even subspace V (the
15-port, cost-5 totals that cliquefloor.py's bound points to)."""
import sys, os
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels, ip, rank2
from decoders import pool
L = design('E7'); labs, k, B = labels(L); n = len(L)
A = np.array(L, float); r = np.linalg.matrix_rank(A)
c = pool(8, A, r)
spec = Counter()
cheap = []
for f in c:
    S = [a for a in range(n) if ip(f, L[a]) != 0]
    if not S: continue
    h = rank2([labs[a] for a in S])
    spec[(len(S), h)] += 1
    if h <= 5: cheap.append((f, len(S), h))
print('pool', len(c), 'spectrum (support, hub cost) -> count:', dict(sorted(spec.items())))
print('cost <= 5 candidates:', len(cheap), cheap[:5])
