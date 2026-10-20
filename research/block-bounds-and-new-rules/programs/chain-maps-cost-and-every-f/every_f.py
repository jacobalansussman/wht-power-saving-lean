# xc-legit (wave 17): "correct at EVERY f" for words with cross-column letters, as a FINITE exact check.
#
# A word whose letters are finite-state across the digit columns (column-uniform letters, chain maps row ^= row << e,
# their inverses, carries, column-0 boundary letters, neighbour phases) acts at f columns as a matrix product
# operator (MPO):  <out| Word_f |in> = alpha^T  T[out_0,in_0] ... T[out_{f-1},in_{f-1}]  beta,   T[o,i] n x n.
# Two such families agree for EVERY f (and every out/in) iff they agree on all column strings of length < n1 + n2
# (weighted-automaton equivalence: Schuetzenberger 1961, Tzeng 1992); the test below computes the forward-reachable
# span of the difference automaton exactly over Q(i) and needs no enumeration of f.
#  Part A: the g-hypotheses 19:16 identity  P D_z^(x)f P^-1 = X^(z x 1) . D_{z e_0}^-1 . CZ^W_z  (P: row ^= row << 1 on
#          every row; CZ^W_z = Walsh-basis CZ chain on y_j = z.w_j) -> proved for every f by the automaton test,
#          and cross-checked by brute force at f = 0..4 against matrices built independently.
#  Part B: why small f alone is not a certificate: K^(x)f . E_e (E_e: row0 ^= row1 << e) equals K^(x)f for f <= e and
#          differs at f = e+1; the automaton test finds it.
# Usage: python3 -I every_f.py        (m = 2; largest automaton n = 20; brute force up to 256 x 256; < 50 MB)
from fractions import Fraction as Fr
import itertools

class G:   # exact Gaussian rationals
    __slots__ = ('re', 'im')
    def __init__(s, re=0, im=0): s.re = Fr(re); s.im = Fr(im)
    def __add__(s, o): return G(s.re + o.re, s.im + o.im)
    def __sub__(s, o): return G(s.re - o.re, s.im - o.im)
    def __mul__(s, o): return G(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    def __neg__(s): return G(-s.re, -s.im)
    def inv(s):
        d = s.re * s.re + s.im * s.im; return G(s.re / d, -s.im / d)
    def iszero(s): return s.re == 0 and s.im == 0
    def __repr__(s): return f"({s.re}+{s.im}i)"
Z0, ONE = G(0), G(1)
A_, B_ = G(Fr(1, 2), Fr(1, 2)), G(Fr(1, 2), Fr(-1, 2))      # a = (1+i)/2, b = (1-i)/2

m = 2; d = 1 << m

class MPO:
    """T: dict (o,i) -> n x n list-of-lists of G; alpha, beta: length-n lists."""
    def __init__(s, n, T, alpha, beta): s.n, s.T, s.alpha, s.beta = n, T, alpha, beta
    def __matmul__(s, o):      # operator product s . o  (s applied after o)
        n = s.n * o.n
        T = {}
        for oo in range(d):
            for ii in range(d):
                M = [[Z0] * n for _ in range(n)]
                for mid in range(d):
                    P, Q = s.T[(oo, mid)], o.T[(mid, ii)]
                    for a in range(s.n):
                        for b in range(s.n):
                            if P[a][b].iszero(): continue
                            for c in range(o.n):
                                for e in range(o.n):
                                    if Q[c][e].iszero(): continue
                                    M[a * o.n + c][b * o.n + e] = M[a * o.n + c][b * o.n + e] + P[a][b] * Q[c][e]
                T[(oo, ii)] = M
        al = [x * y for x in s.alpha for y in o.alpha]; be = [x * y for x in s.beta for y in o.beta]
        return MPO(n, T, al, be)

def local(Mfun):              # bond-1 MPO of a column-uniform local matrix
    return MPO(1, {(o, i): [[Mfun(o, i)]] for o in range(d) for i in range(d)}, [ONE], [ONE])

def Cz(z):    return lambda o, i: A_ if o == i else (B_ if o == i ^ z else Z0)
def Czinv(z): return lambda o, i: B_ if o == i else (A_ if o == i ^ z else Z0)
def Xz(z):    return lambda o, i: ONE if o == i ^ z else Z0
def Kfull():  # C^(x)m on one column: product of C along each bit
    def k(o, i):
        v = ONE
        for b in range(m): v = v * (A_ if ((o >> b) & 1) == ((i >> b) & 1) else B_)
        return v
    return k
def Hu(scale):
    return lambda o, i: G(Fr((-1) ** bin(o & i).count('1'), scale))

def chain(transpose=False):   # P[y,x] = prod_j [y_j = x_j ^ x_{j-1}], state = x_{j-1}
    T = {}
    for y in range(d):
        for x in range(d):
            M = [[Z0] * d for _ in range(d)]
            for s in range(d):
                if y == x ^ s: M[s][x] = ONE
            T[(x, y) if transpose else (y, x)] = M
    return MPO(d, T, [ONE] + [Z0] * (d - 1), [ONE] * d)

def col0(Mfun):               # Mfun on column 0 only (state 0 = column 0 not yet read)
    T = {}
    for o in range(d):
        for i in range(d):
            T[(o, i)] = [[Z0, Mfun(o, i)], [Z0, ONE if o == i else Z0]]
    return MPO(2, T, [ONE, Z0], [ONE, ONE])

def czdiag(z):                # diag (-1)^(sum_j y_j y_{j+1}),  y_j = parity(z & w_j); state = y_{j-1}
    T = {}
    for o in range(d):
        for i in range(d):
            M = [[Z0, Z0], [Z0, Z0]]
            if o == i:
                y = bin(z & i).count('1') & 1
                for s in (0, 1): M[s][y] = G((-1) ** (s * y))
            T[(o, i)] = M
    return MPO(2, T, [ONE, Z0], [ONE, ONE])

def shiftreg(e):              # E_e: bit 0 of column j += bit 1 of column j-e (row0 ^= row1 << e); state = last e row-1 bits
    ns = 1 << e; T = {}
    for o in range(d):
        for i in range(d):
            M = [[Z0] * ns for _ in range(ns)]
            for s in range(ns):
                oldest = (s >> (e - 1)) & 1
                if o == i ^ oldest:
                    M[s][((s << 1) | ((i >> 1) & 1)) & (ns - 1)] = ONE
            T[(o, i)] = M
    return MPO(ns, T, [ONE] + [Z0] * (ns - 1), [ONE] * ns)

def equivalent(A, B):
    """Exact weighted-automaton equivalence of two MPO families over all lengths f >= 0.
    Returns (True, depth) or (False, witness length)."""
    n = A.n + B.n
    alpha = A.alpha + B.alpha; beta = A.beta + [-x for x in B.beta]
    def step(v, key):
        TA, TB = A.T[key], B.T[key]; w = [Z0] * n
        for a in range(A.n):
            if v[a].iszero(): continue
            for b in range(A.n):
                if not TA[a][b].iszero(): w[b] = w[b] + v[a] * TA[a][b]
        for a in range(B.n):
            if v[A.n + a].iszero(): continue
            for b in range(B.n):
                if not TB[a][b].iszero(): w[A.n + b] = w[A.n + b] + v[A.n + a] * TB[a][b]
        return w
    basis = []   # list of (pivot, vector) in reduced form
    def reduce(v):
        v = v[:]
        for p, bv in basis:
            if not v[p].iszero():
                c = v[p]; v = [x - c * y for x, y in zip(v, bv)]
        return v
    def insert(v):
        v = reduce(v)
        p = next((k for k, x in enumerate(v) if not x.iszero()), None)
        if p is None: return False
        c = v[p].inv(); v = [x * c for x in v]
        for idx, (q, bv) in enumerate(basis):
            if not bv[p].iszero():
                cc = bv[p]; basis[idx] = (q, [x - cc * y for x, y in zip(bv, v)])
        basis.append((p, v)); return True
    level = [alpha]; insert(alpha); depth = 0
    def dot(v):
        s = Z0
        for x, y in zip(v, beta): s = s + x * y
        return s
    if not dot(alpha).iszero(): return (False, 0)
    while level:
        depth += 1; nxt = []
        for v in level:
            for key in A.T:
                w = step(v, key)
                if not dot(w).iszero(): return (False, depth)
                if insert(w): nxt.append(w)
        level = nxt
    return (True, depth - 1)

def expand(Mp, f):            # full d^f x d^f matrix from the MPO (column 0 = least significant digit)
    N = d ** f; out = [[Z0] * N for _ in range(N)]
    for oi in range(N):
        os_ = [(oi // d ** j) % d for j in range(f)]
        for ii in range(N):
            is_ = [(ii // d ** j) % d for j in range(f)]
            v = Mp.alpha
            for j in range(f):
                M = Mp.T[(os_[j], is_[j])]
                v = [sum((v[a] * M[a][b] for a in range(Mp.n)), Z0) for b in range(Mp.n)]
            out[oi][ii] = sum((x * y for x, y in zip(v, Mp.beta)), Z0)
    return out

def direct_lhs(z, f):         # independent brute force of P D_z^(x)f P^-1 via permutation and butterflies
    N = d ** f
    def L(x):  cols = [(x // d ** j) % d for j in range(f)]; return sum((cols[j] ^ (cols[j - 1] if j else 0)) * d ** j for j in range(f))
    Pm = [[ONE if y == L(x) else Z0 for x in range(N)] for y in range(N)]
    Pinv = [[Pm[x][y] for x in range(N)] for y in range(N)]
    Dm = [[Z0] * N for _ in range(N)]
    for o in range(N):
        for i in range(N):
            v = ONE
            for j in range(f):
                oj, ij = (o // d ** j) % d, (i // d ** j) % d
                v = v * (A_ if oj == ij else (B_ if oj == ij ^ z else Z0))
                if v.iszero(): break
            Dm[o][i] = v
    def mm(X, Y):
        return [[sum((X[r][t] * Y[t][c] for t in range(N) if not X[r][t].iszero()), Z0) for c in range(N)] for r in range(N)]
    return mm(mm(Pm, Dm), Pinv)

def same(X, Y): return all((x - y).iszero() for rx, ry in zip(X, Y) for x, y in zip(rx, ry))

if __name__ == '__main__':
    print(f"m = {m}: Part A, P D_z^(x)f P^-1 = X^(z x 1) . D_(z e_0)^-1 . CZ^W_z for every f")
    for z in (1, 2, 3):
        lhs = chain() @ local(Cz(z)) @ chain(transpose=True)
        czw = local(Hu(d)) @ czdiag(z) @ local(Hu(1))
        rhs = local(Xz(z)) @ col0(Czinv(z)) @ czw
        res = equivalent(lhs, rhs)
        bf = all(same(direct_lhs(z, f), expand(rhs, f)) for f in range(0, 4))
        bf2 = all(same(direct_lhs(z, f), expand(lhs, f)) for f in range(0, 4))
        print(f"  z={z}: automata n={lhs.n}+{rhs.n}; every-f equivalence: {res}; brute force f=0..3 vs RHS: {bf}, vs LHS MPO: {bf2}",
              flush=True)
        # the paid part alone is not column-uniform: is the conjugated move equal to (free column-uniform) x D^(x)f ?
        alt = local(Xz(z)) @ local(Czinv(z))      # a column-uniform guess: X . D^-1 on every column
        print(f"        control: same identity with D^-1 on EVERY column instead of column 0: {equivalent(lhs, alt @ czw)}")
    print("Part B: K^(x)f . E_e versus K^(x)f (E_e = row0 ^= row1 << e)")
    for e in (1, 2, 3, 4):
        res = equivalent(local(Kfull()) @ shiftreg(e), local(Kfull()))
        print(f"  e={e}: automata n={1 << e}+1; result {res}  (expected: first difference at f = {e + 1})", flush=True)
