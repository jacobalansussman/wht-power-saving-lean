# pf-ladder: coordinate certificates on the FRAME permutation module (engine game).  V = Q[frames], Y_A = span of the
# frames at distance in A from L0 (A a subset of {0..m}); kept by every free letter of the engine game (they fix L0 and
# act on frames by isometries of the dual polar graph DSp(2m,2)).  f(g) = |S_A(L0) minus S_A(gL0)| depends only on
# i = dist(L0, gL0):  f_i = sum_{a in A} k_a - sum_{a,b in A} p^i_{ab}.  Direction: i = 1; kernel: i = m; rank-2 block: 2.
# Intersection array of DSp(2m,2): b_i = 2^(i+1) (2^(m-i) - 1), c_i = 2^i - 1 (checked against m = 3: 14,12,8 / 1,3,7).
import itertools, sys
from fractions import Fraction as Fr

def pmats(m):
    b = [2 ** (i + 1) * (2 ** (m - i) - 1) for i in range(m + 1)]
    c = [2 ** i - 1 for i in range(m + 1)]
    k = b[0]; a = [k - b[i] - c[i] for i in range(m + 1)]
    n = m + 1
    M = [[[Fr(0)] * n for _ in range(n)] for _ in range(n)]
    for j in range(n): M[0][j][j] = Fr(1)
    for j in range(n):
        if j >= 1: M[1][j - 1][j] = Fr(b[j - 1])
        M[1][j][j] = Fr(a[j])
        if j + 1 < n: M[1][j + 1][j] = Fr(c[j + 1])
    def mul(X, Y): return [[sum(X[r][t] * Y[t][s] for t in range(n)) for s in range(n)] for r in range(n)]
    for i in range(1, m):
        P = mul(M[1], M[i])
        M[i + 1] = [[(P[r][s] - b[i - 1] * M[i - 1][r][s] - a[i] * M[i][r][s]) / c[i + 1] for s in range(n)] for r in range(n)]
    kk = [M[i][0][i] for i in range(n)]          # k_i = p^0_{ii} = M[i][0][i]
    return M, kk, b, c, a

for m in range(2, 9):
    M, kk, b, c, a = pmats(m)
    n = m + 1
    assert sum(kk) == __import__('math').prod(2 ** i + 1 for i in range(1, m + 1)), (m, kk)
    # careful: M[x] is the matrix of multiplication by A_x: (A_x A_y) = sum_k M[x][k][y] A_k, so p^k_{xy} = M[x][k][y]
    best = []
    for size in range(1, n):
        for A in itertools.combinations(range(n), size):
            tot = sum(kk[x] for x in A)
            f = [tot - sum(M[x][i][y] for x in A for y in A) for i in range(n)]
            if f[1] > 0: best.append((Fr(f[m], f[1]), A, f))
    best.sort(reverse=True)
    print(f"m={m}: frames {sum(kk)}, sphere sizes {[int(x) for x in kk]}")
    for r, A, f in best[:4]:
        print(f"   A={A}: f_1..f_m = {[int(x) for x in f[1:]]}, ratio f_m/f_1 = {float(r):.4f} ({r}); "
              f"f_2/f_1 = {float(Fr(f[2], f[1])) if m >= 2 else 0:.3f}")
