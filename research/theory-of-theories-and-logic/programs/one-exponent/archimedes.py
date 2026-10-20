# meta-spectrum archimedes.py: the boundedness identities of the composition semiring, checked exactly.
# Letters (engine, w-bounds section 1): dir_z = ca*I + cb*X_z, ca = (1+i)/2, cb = (1-i)/2; kernel K = dir_e1 ... dir_em;
# free: Paulis, the phase letter S_j = diag(1, i) on bit j, address maps.
# Claims checked (exact arithmetic in Q(i), matrices of size 2^m, m = 1..3):
#  (a) dir_z^2 = X_z, so K^2 = X_all (a free Pauli) and K^-1 = K X_all.
#  (b) dir_e1 = S_1 . K . S_1^-1 . K^-1 . S_1^-1 up to a scalar in {1,i,-1,-i}   (exponents found by search below)
#      i.e. one direction = 2 kernel calls + free letters:  [1] <= 2 [m]  ("l x >= 1" of Strassen's Definition 2.4).
#  (c) K = dir_e1 ... dir_em:  [m] <= m [1]  ("l >= x").
from fractions import Fraction as F
from itertools import product

class G:  # Gaussian rational a + b i
    __slots__ = ('a', 'b')
    def __init__(s, a, b=0): s.a, s.b = F(a), F(b)
    def __add__(s, o): return G(s.a + o.a, s.b + o.b)
    def __sub__(s, o): return G(s.a - o.a, s.b - o.b)
    def __mul__(s, o): return G(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a)
    def __eq__(s, o): return s.a == o.a and s.b == o.b
    def iszero(s): return s.a == 0 and s.b == 0

Z0, ONE, I = G(0), G(1), G(0, 1)
def mat(n): return [[Z0] * n for _ in range(n)]
def mul(A, B):
    n = len(A); C = mat(n)
    for r in range(n):
        for k in range(n):
            if A[r][k].iszero(): continue
            a = A[r][k]
            for c in range(n):
                if not B[k][c].iszero(): C[r][c] = C[r][c] + a * B[k][c]
    return C
def eye(n):
    E = mat(n)
    for r in range(n): E[r][r] = ONE
    return E
def scal(s, A): return [[s * x for x in row] for row in A]
def eq(A, B): return all(A[r][c] == B[r][c] for r in range(len(A)) for c in range(len(A)))

ca, cb = G(F(1, 2), F(1, 2)), G(F(1, 2), F(-1, 2))
def Xz(m, z):
    n = 2 ** m; A = mat(n)
    for x in range(n): A[x ^ z][x] = ONE
    return A
def dirz(m, z):
    n = 2 ** m; A = scal(ca, eye(n)); X = Xz(m, z)
    return [[A[r][c] + cb * X[r][c] for c in range(n)] for r in range(n)]
def S(m, j, e):  # diag(1, i)^e on bit j
    n = 2 ** m; A = mat(n); ie = [ONE, I, G(-1), G(0, -1)][e % 4]
    for x in range(n): A[x][x] = ie if (x >> j) & 1 else ONE
    return A

for m in (1, 2, 3):
    n = 2 ** m
    K = eye(n)
    for j in range(m): K = mul(K, dirz(m, 1 << j))
    Xall = Xz(m, n - 1)
    ok_a = eq(mul(dirz(m, 1), dirz(m, 1)), Xz(m, 1)) and eq(mul(K, K), Xall)
    Kinv = mul(K, Xall); ok_inv = eq(mul(K, Kinv), eye(n))
    D1 = dirz(m, 1)
    found = []
    for a, b, c in product(range(4), repeat=3):
        W = mul(mul(mul(mul(S(m, 0, a), K), S(m, 0, b)), Kinv), S(m, 0, c))
        for s in (ONE, I, G(-1), G(0, -1)):
            if eq(scal(s, W), D1): found.append((a, b, c, (s.a, s.b)))
    print(f"m={m}: dir^2 = X and K^2 = X_all: {ok_a}; K*K X_all = I: {ok_inv}; "
          f"dir_e1 = s * S^a K S^b K^-1 S^c for (a,b,c,s) in {found[:4]}{' ...' if len(found) > 4 else ''} ({len(found)} solutions)")
    # (c) is the definition of K.
