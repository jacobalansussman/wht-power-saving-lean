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
TUP = {1: [(1,)], 2: [(1,), (3,), (1, 2), (3, 1)], 3: [(1,), (5,), (1, 2), (3, 5), (1, 2, 4), (1, 3, 5)]}

print("A. Weights rank[g^(x)f, Pi12] (witness: first two digits equal).  Predicted 2(q - 2^(m-k)) q^(f-2).")
for m, f in ((1, 2), (2, 2), (3, 2), (1, 3), (2, 3), (1, 4)):
    q, shift, dirz, phase, free = setup(m); n = q ** f
    Pd = [1 if digits(a, q, f)[0] == digits(a, q, f)[1] else 0 for a in range(n)]
    out = []
    for zs in TUP[m]:
        k = len(zs); out.append("k=%d z=%s: %s (pred %d)" % (k, zs, rk(comm(kp(blk(dirz, zs, q), f), Pd)), 2 * (q - 2 ** (m - k)) * q ** (f - 2)))
    fr = all(iszero(comm(kp(b, f), Pd)) for b in free)
    print(" m=%d f=%d n=%d: free generators blind: %s; %s" % (m, f, n, fr, "; ".join(out))); sys.stdout.flush()

print()
print("B. A block applied ONLY on the fibres that meet the witness set (digits 1,2 in the same coset of V), f = 2.")
for m in (2, 3):
    q, shift, dirz, phase, free = setup(m); n = q * q
    Pd = [1 if a // q == a % q else 0 for a in range(n)]
    for zs in TUP[m]:
        k = len(zs)
        if k == m: continue
        V = {0}
        for z in zs: V |= {v ^ z for v in V}
        S = [((a // q) ^ (a % q)) in V for a in range(n)]
        g2 = kp(blk(dirz, zs, q), 2)
        PS = [[(g2[r][c] if (S[r] and S[c]) else (1 + 0j if (r == c and not S[r]) else 0j)) for c in range(n)] for r in range(n)]
        leak = any(g2[r][c] != 0 for r in range(n) for c in range(n) if S[r] != S[c])
        print(" m=%d k=%d z=%s: entries touched %d of %d (fraction %s); block crosses the set boundary: %s; weight of the partial block: %s; weight of the whole-array block: %d"
              % (m, k, zs, sum(S), n, Fr(sum(S), n), leak, rk(comm(PS, Pd)), 2 * (q - 2 ** (m - k)))); sys.stdout.flush()

print()
print("C. Hatch H2 for more sizes: one-column move = G (a I + b * shift on ALL columns) G, G the address map x_d <- x_d + x_c.")
for m, f in ((1, 2), (1, 3), (2, 2), (1, 4)):
    q, shift, dirz, phase, free = setup(m); n = q ** f
    okall = True; tot = eye(n); nadd = 0
    for c in range(f):
        def Gimg(a):
            d = digits(a, q, f); d = [d[t] ^ (d[c] if t != c else 0) for t in range(f)]
            return sum(d[t] * q ** (f - 1 - t) for t in range(f))
        G = [[1 + 0j if y == Gimg(x) else 0j for x in range(n)] for y in range(n)]
        for j in range(m):
            z = 1 << j; A = kp(shift(z), f)
            mid = [[((1 + I) / 2 if r == cc else 0j) + (1 - I) / 2 * A[r][cc] for cc in range(n)] for r in range(n)]
            step = mm(G, mm(mid, G))
            one = [[1 + 0j]]
            for t in range(f): one = kron(one, dirz(z) if t == c else eye(q))
            okall = okall and step == one; tot = mm(step, tot); nadd += 1
    K = blk(dirz, [1 << j for j in range(m)], q)
    print(" m=%d f=%d n=%d: every one-column move equals (address map)(one whole-array addition)(address map): %s; product = kernel on all columns: %s; whole-array additions used: %d = m*f; blocks used: 0"
          % (m, f, n, okall, tot == kp(K, f), nadd)); sys.stdout.flush()

print()
print("D. Masks from the engine's own letters: (1 + Ph + Ph^2 + Ph^3)/4 with Ph = phase(z,0) on every column (rational coefficients).")
for m, f in ((1, 2), (2, 2), (1, 3)):
    q, shift, dirz, phase, free = setup(m); n = q ** f
    Ph = kp(phase(1, 0), f); acc = eye(n); pw = eye(n)
    for _ in range(3):
        pw = mm(Ph, pw); acc = [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(acc, pw)]
    Mk = [[z / 4 for z in row] for row in acc]
    diag01 = all((Mk[r][c] in (0, 1)) if r == c else Mk[r][c] == 0 for r in range(n) for c in range(n))
    print(" m=%d f=%d: is a 0/1 diagonal matrix: %s; keeps %d of %d addresses" % (m, f, diag01, sum(1 for r in range(n) if Mk[r][r] == 1), n))

print()
print("E. Column permutations as free steps (f = 3).  Witness 'all three digits equal' and witness 'number of equal pairs'.")
for m in (1, 2):
    q, shift, dirz, phase, free = setup(m); f = 3; n = q ** f
    dg = [digits(a, q, f) for a in range(n)]
    Pall = [1 if d[0] == d[1] == d[2] else 0 for d in dg]
    Pnum = [(d[0] == d[1]) + (d[0] == d[2]) + (d[1] == d[2]) for d in dg]
    P12 = [1 if d[0] == d[1] else 0 for d in dg]
    def colperm(s):
        img = lambda a: sum(dg[a][s[t]] * q ** (f - 1 - t) for t in range(f))
        return [[1 + 0j if y == img(x) else 0j for x in range(n)] for y in range(n)]
    sw12, sw23 = colperm((1, 0, 2)), colperm((0, 2, 1))
    b1 = all(iszero(comm(kp(b, f), Pall)) and iszero(comm(kp(b, f), Pnum)) for b in free)
    print(" m=%d: swap(2,3) blind to 'first two equal': %s; swaps blind to 'all equal': %s; to 'number of equal pairs': %s; free generators blind to both: %s"
          % (m, iszero(comm(sw23, P12)), iszero(comm(sw12, Pall)) and iszero(comm(sw23, Pall)), iszero(comm(sw12, Pnum)) and iszero(comm(sw23, Pnum)), b1))
    for zs in TUP[m]:
        g3 = kp(blk(dirz, zs, q), f)
        print("   k=%d z=%s: weight under 'all equal': %s; under 'number of equal pairs': %s" % (len(zs), zs, rk(comm(g3, Pall)), rk(comm(g3, Pnum)))); sys.stdout.flush()
