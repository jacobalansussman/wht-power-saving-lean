# v-ex2-duality a5_exact.py: EXACT (Gaussian-integer, int64 matmul, no BLAS) check of spectrum A5 and of
# certificate A's prices.  Matrices are pairs (Re, Im) of integer arrays.
import sys
import numpy as np
from functools import reduce
def gm(A, B): return (A[0] @ B[0] - A[1] @ B[1], A[0] @ B[1] + A[1] @ B[0])
def gk(A, B): return (np.kron(A[0], B[0]) - np.kron(A[1], B[1]), np.kron(A[0], B[1]) + np.kron(A[1], B[0]))
def geq(A, B): return np.array_equal(A[0], B[0]) and np.array_equal(A[1], B[1])
def gs(c, A): return (c * A[0], c * A[1])
I2 = (np.eye(2, dtype=np.int64), np.zeros((2, 2), np.int64))
X = (np.array([[0, 1], [1, 0]], np.int64), np.zeros((2, 2), np.int64))
Zm = (np.array([[1, 0], [0, -1]], np.int64), np.zeros((2, 2), np.int64))
D2 = (np.array([[1, 1], [1, 1]], np.int64), np.array([[1, -1], [-1, 1]], np.int64))   # 2*dir = (1+i)I + (1-i)X
S = (np.array([[1, 0], [0, 0]], np.int64), np.array([[0, 0], [0, 1]], np.int64))      # diag(1, i)
def on(m, j, A): return reduce(gk, [A if t == j else I2 for t in range(m)])
for m in range(1, 7):
    K2 = reduce(gm, [on(m, j, D2) for j in range(m)])          # 2^m K
    Xall = reduce(gk, [X] * m)
    S1 = on(m, 0, S); S13 = gm(gm(S1, S1), S1)
    Id = (np.eye(2 ** m, dtype=np.int64), np.zeros((2 ** m, 2 ** m), np.int64))
    inv_ok = geq(gm(K2, gm(K2, Xall)), gs(4 ** m, Id))        # K * (K X_all) = I
    lhs = gm(gm(gm(gm(S1, K2), S1), gm(K2, Xall)), S13)       # 4^m * S1 K S1 K^-1 S1^3
    rhs = gs(4 ** m // 2, on(m, 0, D2))                       # 4^m * dir_e1
    print(f"A5 exact m={m}: K K X_all = I: {inv_ok}; S1 K S1 K^-1 S1^3 = dir_e1: {geq(lhs, rhs)}")
# certificate A: conjugation action on Pauli labels; rank-r block = dir_e1...dir_er; count labels of L0 = {Z^z}
def pauli(m, x, z): return reduce(gk, [gm(X if (x >> t) & 1 else I2, Zm if (z >> t) & 1 else I2) for t in range(m)])
for w in (3, 4, 5):
    n = 2 ** w
    for r in range(1, w + 1):
        g2 = reduce(gm, [on(w, j, D2) for j in range(r)])                 # 2^r g
        # g^-1 = (dir_e1..dir_er)^-1 = X_{e1..er} * g  (dir^2 = X): g^-1 * 2^r = Xr @ g2
        Xr = reduce(gk, [X if t < r else I2 for t in range(w)])
        g2i = gm(Xr, g2)
        assert geq(gm(g2, g2i), gs(4 ** r, (np.eye(n, dtype=np.int64), np.zeros((n, n), np.int64))))
        moved = 0
        for z in range(1, n):
            P = gm(gm(g2, pauli(w, 0, z)), g2i)                            # 4^r g Z^z g^-1
            # is it a multiple of a diagonal (Z-type) Pauli?  off-diagonal part zero <=> stays in L0
            off = (P[0] - np.diag(np.diag(P[0]))), (P[1] - np.diag(np.diag(P[1])))
            if np.any(off[0]) or np.any(off[1]): moved += 1
        print(f"cert A exact, width {w}: rank-{r} block moves {moved} of the {n - 1} labels of L0; formula 2^w - 2^(w-r) = {n - 2 ** (w - r)}")
