"""hensel_demo.py -- v-ml-fooling: a ONE-PRIME existence proof for a SQUARE slice (K3_coord, 12 equations, 12 unknowns).
Certificate: a point x mod p with f(x) = 0 mod p and det Jf(x) != 0 mod p (integer-coefficient equations). Multivariate
Hensel lifts it to a zero in Z_p^12, a field of characteristic 0 containing Q, so 1 is not in the ideal over Q and (weak
Nullstellensatz) a complex zero exists. The point is FOUND here from ml-fooling's matrices (how a certificate is found
does not matter); it is then CHECKED without them, and lifted to p^16 as a sanity check."""
import json, sys
from fractions import Fraction
from math import lcm
C = json.load(open('../../data/two-certificate-theorem/certs/K3_coord.json'))
n = len(C['unknowns'])
eqs = []
for f in C['equations']:
    den = 1
    for e, c in f: den = lcm(den, Fraction(c).denominator)
    eqs.append([(tuple(e), int(Fraction(c) * den)) for e, c in f])
def ev(f, x, m): 
    s = 0
    for e, c in f:
        t = c
        for i, k in enumerate(e): t = t * pow(x[i], k, m) % m
        s += t
    return s % m
def jac(x, m):
    J = []
    for f in eqs:
        row = []
        for j in range(n):
            s = 0
            for e, c in f:
                if e[j] == 0: continue
                t = c * e[j]
                for i, k in enumerate(e):
                    kk = k - 1 if i == j else k
                    t = t * pow(x[i], kk, m) % m
                s += t
            row.append(s % m)
        J.append(row)
    return J
def det_mod(A, p):
    A = [r[:] for r in A]; N = len(A); d = 1
    for c in range(N):
        piv = next((r for r in range(c, N) if A[r][c] % p), None)
        if piv is None: return 0
        if piv != c: A[c], A[piv] = A[piv], A[c]; d = -d
        d = d * A[c][c] % p; inv = pow(A[c][c], -1, p)
        for r in range(c + 1, N):
            f = A[r][c] * inv % p
            A[r] = [(A[r][j] - f * A[c][j]) % p for j in range(N)]
    return d % p
def solve_mod(A, b, m_p, p):   # solve A y = b mod p (A invertible mod p)
    N = len(A); M = [A[i][:] + [b[i]] for i in range(N)]
    for c in range(N):
        piv = next(r for r in range(c, N) if M[r][c] % p)
        M[c], M[piv] = M[piv], M[c]; inv = pow(M[c][c], -1, p)
        M[c] = [v * inv % p for v in M[c]]
        for r in range(N):
            if r != c and M[r][c] % p:
                f = M[r][c]; M[r] = [(M[r][j] - f * M[c][j]) % p for j in range(N + 1)]
    return [M[i][N] for i in range(N)]
Ms = [[[Fraction(x) for x in row] for row in M] for M in C['matrices']]; r = len(Ms[0])
def isprime(q): return q > 1 and all(q % d for d in range(2, int(q ** .5) + 1))
for p in range(10007, 40000, 2):
    if not isprime(p): continue
    try:
        Mp = [[[x.numerator * pow(x.denominator, -1, p) % p for x in row] for row in M] for M in Ms]
    except ValueError: continue
    S = [[sum((t + 1) * Mp[t][i][j] for t in range(n)) % p for j in range(r)] for i in range(r)]
    found = None
    for lam in range(p):
        A = [[(S[i][j] - (lam if i == j else 0)) % p for j in range(r)] for i in range(r)]
        if det_mod(A, p) == 0: found = lam; break
    if found is None: continue
    # eigenvector: null vector of A mod p
    A = [[(S[i][j] - (found if i == j else 0)) % p for j in range(r)] for i in range(r)]
    # brute small null-space solve: fix last coordinate free
    M = [row[:] for row in A]; piv = []; rr = 0
    for c in range(r):
        pr = next((i for i in range(rr, r) if M[i][c] % p), None)
        if pr is None: continue
        M[rr], M[pr] = M[pr], M[rr]; inv = pow(M[rr][c], -1, p); M[rr] = [v * inv % p for v in M[rr]]
        for i in range(r):
            if i != rr and M[i][c]: f = M[i][c]; M[i] = [(M[i][j] - f * M[rr][j]) % p for j in range(r)]
        piv.append(c); rr += 1
    free = [c for c in range(r) if c not in piv][0]; v = [0] * r; v[free] = 1
    for i, c in enumerate(piv): v[c] = (-M[i][free]) % p
    k = next(i for i in range(r) if v[i])
    x = [sum(Mp[t][k][j] * v[j] for j in range(r)) * pow(v[k], -1, p) % p for t in range(n)]
    ok = all(ev(f, x, p) == 0 for f in eqs); dj = det_mod(jac(x, p), p)
    print(f"p = {p}: point mod p {x}; f(x) = 0 mod p: {ok}; det J(x) mod p = {dj}")
    if ok and dj:
        # Hensel/Newton lift to p^16
        m = p; X = x[:]
        for step in range(4):
            m2 = m * m
            F = [ev(f, X, m2) for f in eqs]   # divisible by m
            J = jac(X, m)
            y = solve_mod([[a % m for a in row] for row in J], [(-Fi) % m2 for Fi in F], None, m) if False else None
            # solve J * d = -F/m  mod m  (m = p^k, J invertible mod p => mod p^k via Newton on inverse) -- do mod p^k with adjugate trick
            # simple: solve mod m using Fraction-free Gauss with unit pivots mod p
            Nn = n; Aug = [J[i][:] + [(-(F[i] // m)) % m] for i in range(Nn)]
            for c in range(Nn):
                pr = next(i for i in range(c, Nn) if Aug[i][c] % p)
                Aug[c], Aug[pr] = Aug[pr], Aug[c]; inv = pow(Aug[c][c], -1, m); Aug[c] = [vv * inv % m for vv in Aug[c]]
                for i in range(Nn):
                    if i != c and Aug[i][c] % m: f_ = Aug[i][c]; Aug[i] = [(Aug[i][j] - f_ * Aug[c][j]) % m for j in range(Nn + 1)]
            d = [Aug[i][Nn] for i in range(Nn)]
            X = [(X[i] + m * d[i]) % m2 for i in range(Nn)]; m = m2
        print(f"   lifted to p^16 ({m.bit_length()} bits): f(X) = 0 mod p^16 for all 12 equations: {all(ev(f, X, m) == 0 for f in eqs)}")
        break
