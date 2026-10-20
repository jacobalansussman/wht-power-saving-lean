# f-model checks, 10 Oct 2026.  Pure Python.  Matrices are built with Python complex numbers whose parts are
# dyadic rationals (exact in binary floating point: every product and sum here is exact), equality tests are
# exact.  Ranks are computed modulo two primes p = 1 mod 4 (i -> a square root of -1, 1/2 -> inverse of 2);
# a rank modulo p is a rigorous LOWER bound for the rank over Q(i); the predicted value is a derived formula.
import itertools, math, random, sys
from fractions import Fraction
I = 1j
PR = [998244353, 1004535809]
def sqrtm1(p):
    for g in range(2, 50):
        r = pow(g, (p - 1) // 4, p)
        if r * r % p == p - 1: return r
SQ = {p: sqrtm1(p) for p in PR}
def cmod(z, p):
    a = Fraction(z.real); b = Fraction(z.imag)
    return (a.numerator * pow(a.denominator, p - 2, p) + b.numerator * pow(b.denominator, p - 2, p) * SQ[p]) % p
def matmod(M, p): return [[cmod(z, p) for z in row] for row in M]
def rank_mod(M, p):
    M = [row[:] for row in M]; r = 0; rows = len(M); cols = len(M[0])
    for c in range(cols):
        if r == rows: break
        k = next((j for j in range(r, rows) if M[j][c]), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]
        inv = pow(M[r][c], p - 2, p); pr = [x * inv % p for x in M[r]]; M[r] = pr
        nzc = [j for j in range(c, cols) if pr[j]]
        for j in range(rows):
            if j != r and M[j][c]:
                f = M[j][c]; row = M[j]
                for t in nzc: row[t] = (row[t] - f * pr[t]) % p
        r += 1
    return r
def mm(A, B):
    Bt = list(zip(*B)); return [[sum(a * b for a, b in zip(row, col) if a and b) for col in Bt] for row in A]
def kron(A, B): return [[a * b for a in ra for b in rb] for ra in A for rb in B]
def eye(n): return [[1 + 0j if i == j else 0j for j in range(n)] for i in range(n)]
def sub(A, B): return [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
def setup(m):
    q = 2 ** m
    def shift(z): return [[1 + 0j if y == x ^ z else 0j for x in range(q)] for y in range(q)]
    def dirz(z):
        S = shift(z); return [[((1 + I) / 2 if y == x else 0j) + (1 - I) / 2 * S[y][x] for x in range(q)] for y in range(q)]
    def phase(z, c): return [[(I ** ((bin(z & x).count("1") + c) % 2) if x == y else 0j) for x in range(q)] for y in range(q)]
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
def prod(ms, q):
    g = eye(q)
    for x in ms: g = mm(x, g)
    return g
def indep(zs):
    span = {0}
    for z in zs:
        if z in span: return False
        span |= {s ^ z for s in span}
    return True

print("PART 1. Witness = 'first two digits equal' (two columns, f = 2).  comm(g) = rank of [g(x)g, Pi].")
random.seed(11)
for m in (1, 2, 3):
    q, shift, dirz, phase, free = setup(m)
    Pi = [1 if (a // q) == (a % q) else 0 for a in range(q * q)]          # diagonal 0/1 projector on V(x)V
    def commutator(g2): return [[g2[r][c] * (Pi[c] - Pi[r]) for c in range(q * q)] for r in range(q * q)]
    def comm_rank(g):
        C = commutator(kron(g, g)); return [rank_mod(matmod(C, p), p) for p in PR]
    def commutes(g): return all(z == 0 for row in commutator(kron(g, g)) for z in row)
    ok = all(commutes(b) for b in free)
    pm = list(range(q)); random.shuffle(pm)
    Pm = [[1 + 0j if y == pm[x] else 0j for x in range(q)] for y in range(q)]
    Dg = [[(complex(random.randint(1, 8), random.randint(1, 8)) if x == y else 0j) for x in range(q)] for y in range(q)]
    ok2 = commutes(Pm) and commutes(Dg) and commutes(mm(Pm, Dg))
    print(f" m={m} q={q}: all {len(free)} free generators commute with the witness: {ok}; a random permutation of the digit and a random diagonal also: {ok2}")
    for k in range(1, m + 1):
        tup = [zs for zs in itertools.permutations(range(1, q), k) if indep(zs)]
        if len(tup) > 60: tup = random.sample(tup, 60)
        outs = {}
        for zs in tup:
            rs = tuple(comm_rank(prod([dirz(z) for z in zs], q))); outs[rs] = outs.get(rs, 0) + 1
        print(f"   blocks of rank k={k}: {len(tup)} independent tuples tried; (rank mod p1, rank mod p2): count = {outs}; predicted 2(q - 2^(m-k)) = {2*(q - 2**(m-k))}")
    K = prod([dirz(1 << j) for j in range(m)], q)
    print(f"   kernel K: {comm_rank(K)}; predicted 2(q-1) = {2*(q-1)}")
    b1 = random.choice(free); b2 = random.choice(free); g = mm(b1, mm(dirz(1), b2))
    print(f"   free * move * free: {comm_rank(g)} (same as one move: {2*(q - 2**(m-1))})")
    if m >= 2:
        d1 = dirz(1); top = q >> 1
        Pp = [[(d1[y][x] if (x & top) else (1 + 0j if x == y else 0j)) if (x & top) == (y & top) else 0j for x in range(q)] for y in range(q)]
        print(f"   NOT column-uniform free step, for contrast: a single-column shift S(x)1 commutes with the witness: {all(z == 0 for row in commutator(kron(shift(1), eye(q))) for z in row)}")
        print(f"   partial call (one move, applied only to the half of the addresses with top bit 1), as P(x)P: {comm_rank(Pp)}")
    sys.stdout.flush()

print()
print("PART 2. Where does the game start to bind?  Is g^{(x)f} a combination of b^{(x)f}, b free (modulo scalars)?")
def closure(gens, q):
    def norm(M):
        c = next(z for row in M for z in row if z != 0); return tuple(tuple(z / c for z in row) for row in M)
    start = norm(eye(q)); seen = {start}; todo = [start]
    while todo:
        g = todo.pop()
        for b in gens:
            h = norm(mm(b, [list(r) for r in g]))
            if h not in seen: seen.add(h); todo.append(h)
    return [[list(r) for r in g] for g in seen]
def kp(g, f):
    r = [[1 + 0j]]
    for _ in range(f): r = kron(r, g)
    return r
def span_test(vecs, targets, p):
    basis = {}                                    # pivot -> reduced row
    def reduce(v):
        v = v[:]
        for c in range(len(v)):
            if v[c]:
                if c in basis:
                    f = v[c]; b = basis[c]
                    for t in range(c, len(v)):
                        if b[t]: v[t] = (v[t] - f * b[t]) % p
                else: return c, v
        return None, v
    for v in vecs:
        c, w = reduce(v)
        if c is not None:
            inv = pow(w[c], p - 2, p); basis[c] = [x * inv % p for x in w]
    return len(basis), [reduce(t)[0] is None for t in targets]
for m, fs, primes in ((1, (1, 2, 3, 4), PR), (2, (1, 2), PR[:1])):
    q, shift, dirz, phase, free = setup(m)
    B = closure(free, q); K = prod([dirz(1 << j) for j in range(m)], q)
    print(f" m={m}: free maps modulo scalars: {len(B)} elements")
    for f in fs:
        out = []
        for p in primes:
            vecs = [[cmod(z, p) for row in kp(b, f) for z in row] for b in B]
            tg = [[cmod(z, p) for row in kp(g, f) for z in row] for g in (dirz(1), K)]
            out.append(span_test(vecs, tg, p))
        print(f"   f={f} (array length {q**f}): dim span = {out[0][0]} of {q**(2*f)}; one move in span: {out[0][1][0]}; kernel in span: {out[0][1][1]}; primes agree: {all(o == out[0] for o in out)} ({len(primes)} prime(s))")
        sys.stdout.flush()

print()
print("PART 3. What the floor  sum_j (1 - 2^-k_j) >= W (1 - 2^-m)  gives for the engine's exponent z (price: sum_j (k_j/m)^z = W)")
for m in (3, 4, 5, 6, 8, 10, 12, 16, 20, 30, 45):
    zs = [math.log((1 - 2.0 ** -m) / (1 - 2.0 ** -k)) / math.log(m / k) for k in range(1, m)]
    lam = (1 - 2.0 ** -m) / (1 - 2.0 ** -(m - 1))
    print(f"   m={m:2d}: blocks per array >= {lam:.12f};  z >= {min(zs):.3e} (attained by rank k={1 + zs.index(min(zs))});  saving exponent 1 - z <= {1 - min(zs):.6f}")
