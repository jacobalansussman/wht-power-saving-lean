"""sk_fredkin.py -- skeptic: independent check of tb4-nonlinear's block-3 Fredkin route (exact Gaussian integers x 2^k)."""
import itertools
n = 8
def mat_mul(P, Q):
    return [[sum(P[i][k] * Q[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
# scaled by 2: a = 1+i, b = 1-i  (true a = (1+i)/2, b = (1-i)/2)
a, b = complex(1, 1), complex(1, -1)
def perm(f):
    return [[1 if f(j) == i else 0 for j in range(n)] for i in range(n)]
def step(f):   # 2 * (a I + b P_f)
    P = perm(f)
    return [[(a if i == j else 0) + b * P[i][j] for j in range(n)] for i in range(n)]
X = lambda z: (lambda x: x ^ z)
t1 = lambda x: x ^ (1 if not (x & 4) else 2)
t2 = lambda x: x ^ (2 if not (x & 4) else 1)
# kernel x 8 = (2C)^{(x)3}, C = [[a,b],[b,a]]
def K8():
    return [[ (lambda x, y: __import__('functools').reduce(lambda u, v: u * v, [(a if ((x >> k) & 1) == ((y >> k) & 1) else b) for k in range(3)]))(i, j) for j in range(n)] for i in range(n)]
K = K8()
P = mat_mul(step(X(4)), mat_mul(step(t2), step(t1)))
print("D(e3)(a+b t2)(a+b t1) == K exactly:", P == K)
Q = mat_mul(step(t2), step(t1)); R = mat_mul(step(X(2)), step(X(1)))
print("(a+b t2)(a+b t1) == D(e2)D(e1) exactly (so the extra route only changes the intermediate frame):", Q == R)
# frame of a+b t1: rows up to Z4 phase, as a sorted multiset
def frame(M):
    rows = []
    for r in M:
        first = next(v for v in r if v != 0)
        for u in (1, 1j, -1, -1j):
            w = first * u
            if w.real > 0 and w.imag >= 0:
                break
        rows.append(tuple((v * u) for v in r))
    return tuple(sorted(rows, key=lambda t: [(z.real, z.imag) for z in t]))
f1 = frame(step(t1))
aff = [frame(step(X(z))) for z in range(1, 8)]
print("depth-1 frame B(a+b t1) equals some affine frame B D(z)?", f1 in aff)
# is t1 affine?
print("t1 affine?", all(t1(x ^ y) ^ t1(0) == t1(x) ^ t1(y) for x in range(n) for y in range(n)))
