"""lawN_cex.py (mw2-singer4): explicit counterexample to the empirical regularity 'top-augmentation rank <= nK1'.
Block m = 3, W = 2, 7 moves:  p1 = D1 x1, p2 = D2 x2, u = D3 (p1 + p2) [the ONLY generator with terms one move from K],
cleanups q1 = D3 x1, q2 = D3 x2, tops t1 = D2 u, t2 = D1 u;  out1 = t1 - X^2 q2,  out2 = t2 - X^1 q1.
dir z = a I + b X^z with a^2 + b^2 = 0, 2ab = 1 (so dir z ^2 = X^z), K = dir e1 dir e2 dir e3 (8 x 8 matrices).
Checks the formal identity out_l = K x_l coefficient by coefficient (pure Python complex arithmetic)."""
import cmath
n = 8
a = cmath.exp(-1j * cmath.pi / 4) / 2 ** 0.5
b = 1j * a
def X(z):
    return [[1 if r == (c ^ z) else 0 for c in range(n)] for r in range(n)]
def add(A, B, s=1):
    return [[A[r][c] + s * B[r][c] for c in range(n)] for r in range(n)]
def mul(*Ms):
    R = Ms[0]
    for M in Ms[1:]:
        R = [[sum(R[r][k] * M[k][c] for k in range(n)) for c in range(n)] for r in range(n)]
    return R
def scal(s, A):
    return [[s * x for x in row] for row in A]
I = X(0)
D = {z: add(scal(a, I), scal(b, X(z))) for z in (1, 2, 4)}
assert max(abs(x - y) for r1, r2 in zip(mul(D[1], D[1]), X(1)) for x, y in zip(r1, r2)) < 1e-12
K = mul(D[1], D[2], D[4])
Z = [[0] * n for _ in range(n)]
# coefficients (on x1, on x2) of each generator
p1 = (D[1], Z); p2 = (Z, D[2])
u = (mul(D[4], p1[0]), mul(D[4], p2[1]))
q1 = (D[4], Z); q2 = (Z, D[4])
t1 = (mul(D[2], u[0]), mul(D[2], u[1]))
t2 = (mul(D[1], u[0]), mul(D[1], u[1]))
out1 = (add(t1[0], mul(X(2), q2[0]), -1), add(t1[1], mul(X(2), q2[1]), -1))
out2 = (add(t2[0], mul(X(1), q1[0]), -1), add(t2[1], mul(X(1), q1[1]), -1))
def err(A, B):
    return max(abs(x - y) for r1, r2 in zip(A, B) for x, y in zip(r1, r2))
print('out1 = K x1 + 0 x2 :', err(out1[0], K) < 1e-12 and err(out1[1], Z) < 1e-12)
print('out2 = 0 x1 + K x2 :', err(out2[0], Z) < 1e-12 and err(out2[1], K) < 1e-12)
print('top augmentation vectors: t1 -> (1, 0), t2 -> (0, 1): rank 2; generators with terms one move from K: only u (nK1 = 1)')
