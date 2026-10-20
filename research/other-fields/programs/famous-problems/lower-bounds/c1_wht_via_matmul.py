"""c1: the reduction chain  WHT lower bound  =>  matrix-multiplication lower bound  (iso-lower, 14 Oct 2026).

Claim checked on small cases (exact integer arithmetic):
 (a) Alman 2021, Lemma 8.3 with k = 2: H_{n^2} x = vec(H_n X H_n^T), i.e. the WHT of size N = n^2 is two
     n x n x n matrix products with one factor fixed (H_n), plus a reshape (free in a circuit).
 (b) Linearisation: any arithmetic circuit (+, -, *, scalar) that, with one input matrix fixed, computes a
     linear map can be turned into a linear circuit with at most 3x as many gates, by keeping for every gate
     its constant term (known in advance) and its degree-1 part:
        (u+v)_1 = u_1 + v_1,   (lam u)_1 = lam u_1,   (u*v)_1 = u0*v_1 + v0*u_1   (<= 3 linear gates).
     Tested on recursive Strassen (7 products per 2x2 level) computing A X for fixed A = H_4 and H_8,
     and the resulting linear circuit is checked against H_n X exactly.
So L_lin(WHT_{n^2}) <= 2 * 3 * MM_arith(n,n,n), and a WHT lower bound N (log N)^eps transfers to
MM_arith(n,n,n) >= n^2 (2 log n)^eps / 6: superlinear in the input size.
"""
import itertools, random
from fractions import Fraction

def hadamard(k):
    H = [[1]]
    for _ in range(k):
        H = [[*r, *r] for r in H] + [[*r, *[-v for v in r]] for r in H]
    return H

def matmul(A, B):
    return [[sum(A[i][t] * B[t][j] for t in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

# (a) ------------------------------------------------------------------------------------------
for a in (1, 2, 3):
    n = 2 ** a
    Hn, HN = hadamard(a), hadamard(2 * a)
    x = [random.randint(-9, 9) for _ in range(n * n)]
    X = [x[i * n:(i + 1) * n] for i in range(n)]            # row-major reshape
    Y = matmul(matmul(Hn, X), [list(r) for r in zip(*Hn)])  # H X H^T
    y_direct = [sum(HN[i][j] * x[j] for j in range(n * n)) for i in range(n * n)]
    y_mm = [Y[i][j] for i in range(n) for j in range(n)]
    print(f"(a) n={n}: WHT_{n*n} == vec(H X H^T): {y_direct == y_mm}")

# (b) ------------------------------------------------------------------------------------------
# A tiny straight-line-program builder that records gates and evaluates symbolically as
# (constant, linear-part) pairs, the linear part being a dict var -> coefficient.
class Prog:
    def __init__(self):
        self.arith_gates = 0
        self.lin_gates = 0
    def const(self, c):
        return (Fraction(c), {})
    def var(self, name):
        return (Fraction(0), {name: Fraction(1)})
    def _lin_add(self, p, q, sp=1, sq=1):
        out = dict(p)
        for k, v in q.items():
            out[k] = out.get(k, 0) * sp + v * sq if k in out else v * sq
        return {k: v for k, v in out.items() if v != 0}
    def add(self, u, v):
        self.arith_gates += 1
        # linear circuit: u_1 + v_1 costs 1 gate if both nonzero
        if u[1] and v[1]:
            self.lin_gates += 1
        return (u[0] + v[0], self._lin_add(u[1], v[1]))
    def sub(self, u, v):
        self.arith_gates += 1
        if u[1] and v[1]:
            self.lin_gates += 1
        return (u[0] - v[0], self._lin_add(u[1], {k: -c for k, c in v[1].items()}))
    def mul(self, u, v):
        self.arith_gates += 1
        u0, u1 = u
        v0, v1 = v
        # (u*v)_1 = u0 v1 + v0 u1 ; u1*v1 is degree 2 and must cancel in the outputs of a linear map
        terms = []
        if u0 != 0 and v1:
            terms.append({k: u0 * c for k, c in v1.items()})
        if v0 != 0 and u1:
            terms.append({k: v0 * c for k, c in u1.items()})
        self.lin_gates += len(terms) + max(0, len(terms) - 1)   # scalings + one addition
        lin = {}
        for t in terms:
            lin = self._lin_add(lin, t)
        quad_free = (not u1) or (not v1)
        if not quad_free:
            self.has_quadratic = True
        return (u0 * v0, lin)

def strassen(P, A, B):
    n = len(A)
    if n == 1:
        return [[P.mul(A[0][0], B[0][0])]]
    h = n // 2
    def blk(M, i, j):
        return [row[j * h:(j + 1) * h] for row in M[i * h:(i + 1) * h]]
    def madd(X, Y):
        return [[P.add(X[i][j], Y[i][j]) for j in range(h)] for i in range(h)]
    def msub(X, Y):
        return [[P.sub(X[i][j], Y[i][j]) for j in range(h)] for i in range(h)]
    A11, A12, A21, A22 = blk(A, 0, 0), blk(A, 0, 1), blk(A, 1, 0), blk(A, 1, 1)
    B11, B12, B21, B22 = blk(B, 0, 0), blk(B, 0, 1), blk(B, 1, 0), blk(B, 1, 1)
    M1 = strassen(P, madd(A11, A22), madd(B11, B22))
    M2 = strassen(P, madd(A21, A22), B11)
    M3 = strassen(P, A11, msub(B12, B22))
    M4 = strassen(P, A22, msub(B21, B11))
    M5 = strassen(P, madd(A11, A12), B22)
    M6 = strassen(P, msub(A21, A11), madd(B11, B12))
    M7 = strassen(P, msub(A12, A22), madd(B21, B22))
    C11 = madd(msub(madd(M1, M4), M5), M7)
    C12 = madd(M3, M5)
    C21 = madd(M2, M4)
    C22 = madd(madd(msub(M1, M2), M3), M6)
    C = [[None] * n for _ in range(n)]
    for i in range(h):
        for j in range(h):
            C[i][j], C[i][j + h], C[i + h][j], C[i + h][j + h] = C11[i][j], C12[i][j], C21[i][j], C22[i][j]
    return C

for a in (2, 3):
    n = 2 ** a
    Hn = hadamard(a)
    P = Prog()
    A = [[P.const(Hn[i][j]) for j in range(n)] for i in range(n)]          # fixed factor: constants
    X = [[P.var((i, j)) for j in range(n)] for i in range(n)]              # variable factor
    C = strassen(P, A, X)
    ok = True
    for i in range(n):
        for j in range(n):
            c0, c1 = C[i][j]
            want = {(t, j): Fraction(Hn[i][t]) for t in range(n)}
            if c0 != 0 or c1 != want:
                ok = False
    print(f"(b) n={n}: Strassen arithmetic gates (incl. gates on constants) = {P.arith_gates}; "
          f"linear-circuit gates after linearisation = {P.lin_gates}; "
          f"ratio = {P.lin_gates / P.arith_gates:.3f} (claim <= 3); outputs == H_n X exactly: {ok}")
