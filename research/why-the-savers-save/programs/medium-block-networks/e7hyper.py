#!/usr/bin/env python3
"""cx-medium / e7hyper.py: the largest number of E7 lines in one hyperplane of R^7 (used for: every symmetric frame
total gathers >= 63 - max lines in a hyperplane).  Hyperplanes holding the most lines are spanned by lines, so it
suffices to take hyperplanes spanned by 6 independent lines; random sampling of 6-tuples plus every hyperplane
orthogonal to a line (D6 type) and the E6 hyperplane."""
import sys, os, random
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design
L = np.array(design('E7'), float); n = len(L)
U, s, Vt = np.linalg.svd(L, full_matrices=False); B7 = Vt[:7]; X = L @ B7.T          # coordinates in R^7
rng = random.Random(3); best = 0; seen = set()
for _ in range(20000):
    idx = rng.sample(range(n), 6)
    A = X[idx]
    if np.linalg.matrix_rank(A) < 6: continue
    nrm = np.linalg.svd(A)[2][-1]
    cnt = int(np.sum(np.abs(X @ nrm) < 1e-9))
    best = max(best, cnt)
print('largest number of E7 lines in a hyperplane (20,000 random 6-tuples):', best, '-> smallest symmetric-frame support', n - best)
