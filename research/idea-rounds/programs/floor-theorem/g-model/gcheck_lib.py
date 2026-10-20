# g-model (skeptic of f-model), 10 Oct 2026.  Pure Python, own code (not f-model's).
# Entries are Python complex numbers with dyadic parts (every sum and product here is exact).
# Rank: EXACT over Q(i) (Fractions, real 2N x 2N embedding) for N <= 16; for N = 64 modulo two primes
# (a rank modulo p is a lower bound for the true rank; the second prime differs from f-model's).
import itertools, sys
from fractions import Fraction as Fr
I = 1j
P1, P2 = 998244353, 469762049
def sq(p):
    for g in range(2, 200):
        r = pow(g, (p - 1) // 4, p)
        if r * r % p == p - 1: return r
SQ = {P1: sq(P1), P2: sq(P2)}
def mm(A, B):
    Bt = list(zip(*B)); return [[sum(a * b for a, b in zip(r, c) if a and b) for c in Bt] for r in A]
def kron(A, B): return [[a * b for a in ra for b in rb] for ra in A for rb in B]
def eye(n): return [[1 + 0j if i == j else 0j for j in range(n)] for i in range(n)]
def kp(g, f):
    r = [[1 + 0j]]
    for _ in range(f): r = kron(r, g)
    return r
def rank_exact(M):
    N, C = len(M), len(M[0])
    R = [[Fr(z.real) for z in row] + [Fr(-z.imag) for z in row] for row in M] + \
        [[Fr(z.imag) for z in row] + [Fr(z.real) for z in row] for row in M]
    r = 0
    for c in range(2 * C):
        k = next((j for j in range(r, 2 * N) if R[j][c]), None)
        if k is None: continue
        R[r], R[k] = R[k], R[r]; pv = R[r]; inv = 1 / pv[c]; nz = [t for t in range(c, 2 * C) if pv[t]]
        for j in range(2 * N):
            if j != r and R[j][c]:
                fct = R[j][c] * inv; row = R[j]
                for t in nz: row[t] -= fct * pv[t]
        r += 1
    assert r % 2 == 0
    return r // 2
def rank_p(M, p):
    s = SQ[p]
    def cm(z):
        a = Fr(z.real); b = Fr(z.imag)
        return (a.numerator * pow(a.denominator, p - 2, p) + b.numerator * pow(b.denominator, p - 2, p) * s) % p
    R = [[cm(z) for z in row] for row in M]; r = 0; N = len(R); C = len(R[0])
    for c in range(C):
        k = next((j for j in range(r, N) if R[j][c]), None)
        if k is None: continue
        R[r], R[k] = R[k], R[r]; inv = pow(R[r][c], p - 2, p); pv = [x * inv % p for x in R[r]]; R[r] = pv
        nz = [t for t in range(c, C) if pv[t]]
        for j in range(N):
            if j != r and R[j][c]:
                fct = R[j][c]; row = R[j]
                for t in nz: row[t] = (row[t] - fct * pv[t]) % p
        r += 1
    return r
def rk(M):
    if len(M) <= 16: return "%d exact" % rank_exact(M)
    return "%d,%d mod 2 primes" % (rank_p(M, P1), rank_p(M, P2))
def comm(A, Pd):                      # [A, diag(Pd)]
    n = len(A); return [[A[r][c] * (Pd[c] - Pd[r]) for c in range(n)] for r in range(n)]
def commM(A, B):
    X, Y = mm(A, B), mm(B, A); return [[x - y for x, y in zip(rx, ry)] for rx, ry in zip(X, Y)]
def iszero(M): return all(z == 0 for row in M for z in row)
def setup(m):
    q = 2 ** m
    shift = lambda z: [[1 + 0j if y == x ^ z else 0j for x in range(q)] for y in range(q)]
    def dirz(z):
        S = shift(z); return [[((1 + I) / 2 if y == x else 0j) + (1 - I) / 2 * S[y][x] for x in range(q)] for y in range(q)]
    phase = lambda z, c: [[(I ** ((bin(z & x).count("1") + c) % 2) if x == y else 0j) for x in range(q)] for y in range(q)]
    def perm(G):
        def img(x):
            y = 0
            for j in range(m):
                if x >> j & 1: y ^= G[j]
            return y
        return [[1 + 0j if y == img(x) else 0j for x in range(q)] for y in range(q)]
    free = [shift(1 << j) for j in range(m)] + [phase(z, c) for z in range(q) for c in (0, 1) if (z, c) != (0, 0)]
    for a in range(m):
        for b in range(m):
            if a != b:
                G = [1 << j for j in range(m)]; G[a] ^= 1 << b; free.append(perm(G))
    return q, shift, dirz, phase, free
def blk(dirz, zs, q):
    g = eye(q)
    for z in zs: g = mm(dirz(z), g)
    return g
def digits(a, q, f): return [(a // q ** (f - 1 - c)) % q for c in range(f)]
