# self-test of bareiss_rank (hook_frames_check.py) against mod-p ranks on random low-rank integer matrices
import random, numpy as np, sys
sys.path.insert(0, '.')
from hook_frames_check import bareiss_rank, rank_mod
rng = random.Random(7); P = 2147483629; bad = 0
for trial in range(60):
    n = rng.choice([5, 9, 20, 40]); r = rng.randint(0, n)
    A = [[rng.randint(-3, 3) for _ in range(r)] for _ in range(n)]
    B = [[rng.randint(-3, 3) for _ in range(n)] for _ in range(r)]
    M = [[sum(A[i][k]*B[k][j] for k in range(r)) for j in range(n)] for i in range(n)]
    # zero some columns / duplicate rows to force skipped pivots
    if n > 5:
        for i in range(n): M[i][1] = 0
    rb = bareiss_rank(M); rp = rank_mod(np.array(M, dtype=np.int64) % P, P)
    if rb != rp: bad += 1; print('MISMATCH', n, r, rb, rp)
print('bareiss self-test: 60 random matrices, mismatches vs mod-p rank:', bad)
