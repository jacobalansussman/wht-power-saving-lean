"""feedback_exact.py -- exact (Fraction) re-check of feedback_state.py Part B at f = 2.

Run: python3 feedback_exact.py   (pure Python fractions, a few seconds)
Closure H(1) = M + C (I - D)^-1 B of (i) the simultaneous square L(x)L (data 16, state 48) and (ii) the per-column
cascade (data 16, state 32), compared with (I + J)(x)(I + J).  lam = 3.
"""
from fractions import Fraction as Fr
from itertools import product

lam = 3


def mat(r, c, f=lambda i, j: 0):
    return [[Fr(f(i, j)) for j in range(c)] for i in range(r)]


def mul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    Bt = list(zip(*B))
    return [[sum((A[i][t] * Bt[j][t] for t in range(k) if A[i][t]), Fr(0)) for j in range(m)] for i in range(n)]


def kron(A, B):
    return [[A[i // len(B)][j // len(B[0])] * B[i % len(B)][j % len(B[0])] for j in range(len(A[0]) * len(B[0]))]
            for i in range(len(A) * len(B))]


def solve(A, B):
    n = len(A)
    M = [A[i][:] + B[i][:] for i in range(n)]
    m = len(B[0])
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        inv = 1 / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


def closure(K, nd):
    n = len(K)
    M = [row[:nd] for row in K[:nd]]
    C = [row[nd:] for row in K[:nd]]
    B = [row[:nd] for row in K[nd:]]
    D = [row[nd:] for row in K[nd:]]
    r = n - nd
    ImD = [[(Fr(1) if i == j else Fr(0)) - D[i][j] for j in range(r)] for i in range(r)]
    X = solve(ImD, B)
    CX = mul(C, X)
    return [[M[i][j] + CX[i][j] for j in range(nd)] for i in range(nd)]


I4 = mat(4, 4, lambda i, j: int(i == j))
J = mat(4, 4, lambda i, j: (1 if j == i + 1 else 0) + (lam if (i == 3 and j == 0) else 0))
IJ = [[I4[i][j] + J[i][j] for j in range(4)] for i in range(4)]
target = kron(IJ, IJ)

# (i) simultaneous square of L = [[I, J],[I, 0]]
L = [[(I4[i][j] if i < 4 and j < 4 else J[i][j - 4] if i < 4 else I4[i - 4][j] if j < 4 else Fr(0))
      for j in range(8)] for i in range(8)]
L2 = kron(L, L)
idx = list(product(range(8), repeat=2))
data = [k for k, t in enumerate(idx) if t[0] < 4 and t[1] < 4]
state = [k for k, t in enumerate(idx) if not (t[0] < 4 and t[1] < 4)]
perm = data + state
Lp = [[L2[a][b] for b in perm] for a in perm]
H_sim = closure(Lp, 16)
diff = max(abs(H_sim[i][j] - target[i][j]) for i in range(16) for j in range(16))
print(f"simultaneous L(x)L: state {len(state)}; closure at z=1 equals (I+J)(x)(I+J): {diff == 0}; max |diff| = {float(diff):.4g}")

# (ii) per-column cascade
n = 16
J1 = kron(J, I4)
J2 = kron(I4, J)
Z = mat(n, n)
Id = mat(n, n, lambda i, j: int(i == j))


def block(rows):
    return [sum((blk[i] for blk in row), []) for row in rows for i in range(n)]


K1 = block([[Id, J1, Z], [Id, Z, Z], [Z, Z, Id]])
K2 = block([[Id, Z, J2], [Z, Id, Z], [Id, Z, Z]])
tick = mul(K2, K1)
H_cas = closure(tick, 16)
diff2 = max(abs(H_cas[i][j] - target[i][j]) for i in range(16) for j in range(16))
print(f"per-column cascade: state 32; closure at z=1 equals (I+J)(x)(I+J): {diff2 == 0}")
