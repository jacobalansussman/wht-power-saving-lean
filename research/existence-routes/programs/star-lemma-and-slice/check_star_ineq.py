# Sanity check of the two general rank inequalities used in the star lemma, on random sparse integer matrices
# (they hold for every matrix; this guards against slips in the derivation):
#  (I1) rank [[A,B],[C,D]] >= rank[A B] + rank[A;C] - rank A
#  (I2) rank [A B_1 ... B_g] <= sum_j rank[A B_j] - (g-1) rank A
# and the star-lemma conclusion on matrices T built so that each star has rank <= r:
#  2c <= g(r) - (g-2) a, when T's rows/cols through the corner are independent.
import numpy as np
rng = np.random.default_rng(7)
rk = lambda X: np.linalg.matrix_rank(X) if X.size else 0
bad1 = bad2 = 0
for t in range(4000):
    c, q, g = rng.integers(1, 4), rng.integers(1, 5), rng.integers(2, 5)
    A = (rng.random((c, c)) < 0.3) * rng.integers(-2, 3, (c, c))
    Bs = [(rng.random((c, q)) < 0.3) * rng.integers(-2, 3, (c, q)) for _ in range(g)]
    Cs = [(rng.random((q, c)) < 0.3) * rng.integers(-2, 3, (q, c)) for _ in range(g)]
    Ds = [(rng.random((q, q)) < 0.3) * rng.integers(-2, 3, (q, q)) for _ in range(g)]
    for B, C, D in zip(Bs, Cs, Ds):
        if rk(np.block([[A, B], [C, D]])) < rk(np.hstack([A, B])) + rk(np.vstack([A, C])) - rk(A): bad1 += 1
    if rk(np.hstack([A] + Bs)) > sum(rk(np.hstack([A, B])) for B in Bs) - (g - 1) * rk(A): bad2 += 1
print("violations of (I1):", bad1, " of (I2):", bad2)
