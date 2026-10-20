# test_rank.py: rank_modp_fast (blocked, float64 BLAS) against a plain Python Gaussian elimination mod p, on random
# low-rank integer matrices and on Hadamard-minor matrices of the kind venn_mask.py uses.  python3 -I test_rank.py
import sys, random, warnings
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from venn_mask import rank_modp_fast, P

def rank_slow(A, p=P):
    M = [[int(x) % p for x in row] for row in A]; r = 0; rows = len(M); cols = len(M[0]) if rows else 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(r + 1, rows):
            if M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
        if r == rows: break
    return r

rng = np.random.default_rng(3); bad = 0; tot = 0
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for trial in range(60):
        R = int(rng.integers(1, 160)); C = int(rng.integers(1, 160)); k = int(rng.integers(0, min(R, C) + 1))
        A = rng.integers(-3, 4, size=(R, k)) @ rng.integers(-3, 4, size=(k, C)) if k else np.zeros((R, C), dtype=np.int64)
        if trial % 3 == 0:  # sparse +-1 entries, like products of minors
            A = A * (rng.random((R, C)) < 0.3)
        a, b = rank_modp_fast(A), rank_slow(A)
        tot += 1; bad += (a != b)
print(f"rank_modp_fast vs plain elimination mod {P}: {tot - bad}/{tot} agree")
