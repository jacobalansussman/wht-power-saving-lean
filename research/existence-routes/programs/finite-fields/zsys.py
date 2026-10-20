"""zsys.py -- ex2-modp: small exact toolkit (pure Python 3, no numpy) for polynomial systems with integer/rational
coefficients, their reductions mod p, p-adic (Hensel) lifting, rational reconstruction and EXACT verification.

Polynomials: dict {exponent tuple (len n): coefficient}.  Coefficients are int (mod m) or Fraction (exact).
A 'system' is (names, [poly, ...]).

Main entry points
  parse_ms(path)                -> names, p, polys mod p          (msolve input format, as nag-rerun exports it)
  to_Q(sys_p1, sys_p2)          -> names, polys over Q            (rational reconstruction mod p1, checked mod p2)
  ev(poly, x, m=None)           -> value (exact if m is None, else mod m)
  jac(polys, n)                 -> Jacobian as polys
  rank_mod(M, p), solve_mod(M, r, p)
  hensel(polys, xbar, p, K, fixed)   digit-by-digit lifting with the Jacobian frozen at xbar (see PROOF.md Prop. 4)
  ratrec(a, m)                  -> Fraction or None
"""
from fractions import Fraction
import re, itertools, random

# ---------------------------------------------------------------- parsing
def parse_ms(path):
    txt = open(path).read().strip().split('\n')
    names = [v.strip() for v in txt[0].split(',')]
    p = int(txt[1].strip())
    body = '\n'.join(txt[2:])
    polys = []
    for chunk in body.split(','):
        chunk = chunk.strip()
        if not chunk:
            continue
        polys.append(parse_poly(chunk, names, p))
    return names, p, polys

def parse_poly(s, names, p=None):
    idx = {v: i for i, v in enumerate(names)}
    n = len(names); P = {}
    s = s.replace(' ', '').replace('-', '+-')
    for term in s.split('+'):
        if not term:
            continue
        fs = term.split('*'); c = 1; e = [0] * n
        for f in fs:
            if re.fullmatch(r'-?\d+', f):
                c *= int(f)
            elif f == '-':
                c = -c
            else:
                neg = f.startswith('-')
                if neg:
                    f = f[1:]; c = -c
                if '^' in f:
                    v, k = f.split('^'); e[idx[v]] += int(k)
                else:
                    e[idx[f]] += 1
        k = tuple(e)
        P[k] = P.get(k, 0) + c
        if p is not None:
            P[k] %= p
    return {k: v for k, v in P.items() if v}

# ---------------------------------------------------------------- reconstruction
def ratrec(a, m):
    """rational reconstruction: Fraction r/s with |r|,|s| <= sqrt(m/2) and r = a s mod m, or None"""
    a %= m
    B = int((m // 2) ** 0.5)
    r0, r1, s0, s1 = m, a, 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B:
        return None
    return Fraction(r1, s1)

def to_Q(sysA, sysB):
    """sysA, sysB = (names, p, polys) for the same system mod two primes; returns (names, polys over Q)."""
    (nA, pA, PA), (nB, pB, PB) = sysA, sysB
    assert nA == nB and len(PA) == len(PB)
    out = []
    for X, Y in zip(PA, PB):
        Q = {}
        for k, v in X.items():
            r = ratrec(v, pA)
            assert r is not None, ('no rational reconstruction', k, v)
            assert (r.numerator - r.denominator * Y.get(k, 0)) % pB == 0, ('primes disagree', k, r)
            Q[k] = r
        assert set(Y) <= set(X), 'monomial present mod p2 only'
        out.append(Q)
    return nA, out

# ---------------------------------------------------------------- arithmetic
def ev(P, x, m=None):
    tot = 0
    for k, c in P.items():
        t = c
        for xi, e in zip(x, k):
            if e:
                t = t * (pow(xi, e, m) if m else xi ** e)
                if m: t %= m
        tot += t
        if m: tot %= m
    if m:
        if isinstance(tot, Fraction):
            tot = tot.numerator * pow(tot.denominator, -1, m) % m
        return tot % m
    return tot

def coef_mod(c, m):
    if isinstance(c, Fraction):
        return c.numerator * pow(c.denominator, -1, m) % m
    return c % m

def red(P, m):
    out = {}
    for k, c in P.items():
        v = coef_mod(c, m)
        if v: out[k] = v
    return out

def deriv(P, i):
    out = {}
    for k, c in P.items():
        if k[i]:
            k2 = list(k); k2[i] -= 1; k2 = tuple(k2)
            out[k2] = out.get(k2, 0) + c * k[i]
    return {k: v for k, v in out.items() if v}

def jac(polys, n):
    return [[deriv(P, i) for i in range(n)] for P in polys]

def subst(P, i, val):
    """substitute x_i = val (exact) -> polynomial in the same variables (x_i absent)"""
    out = {}
    for k, c in P.items():
        k2 = list(k); e = k2[i]; k2[i] = 0; k2 = tuple(k2)
        out[k2] = out.get(k2, 0) + c * (val ** e)
    return {k: v for k, v in out.items() if v}

def degree(P):
    return max((sum(k) for k in P), default=-1)

# ---------------------------------------------------------------- linear algebra mod p and over Q
def rref_mod(M, p, ncols=None):
    """row-reduce a copy of M (list of lists) mod p; return (R, pivots)"""
    R = [[x % p for x in row] for row in M]
    ncols = len(R[0]) if R else 0 if ncols is None else ncols
    piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(R)) if R[i][c]), None)
        if pr is None:
            continue
        R[r], R[pr] = R[pr], R[r]
        inv = pow(R[r][c], -1, p)
        R[r] = [x * inv % p for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c]:
                f = R[i][c]
                R[i] = [(a - f * b) % p for a, b in zip(R[i], R[r])]
        piv.append(c); r += 1
        if r == len(R):
            break
    return R, piv

def rank_mod(M, p):
    if not M: return 0
    return len(rref_mod(M, p)[1])

def solve_mod(M, rhs, p):
    """one solution x of M x = rhs mod p (free variables 0), or None if inconsistent"""
    if not M:
        return []
    n = len(M[0])
    A = [list(row) + [b] for row, b in zip(M, rhs)]
    R, piv = rref_mod(A, p, n + 1)
    if n in piv:
        return None
    x = [0] * n
    for i, c in enumerate(piv):
        x[c] = R[i][n]
    return x

def solve_Q(M, rhs):
    """exact solution over Q (free vars 0) or None"""
    n = len(M[0]) if M else 0
    A = [[Fraction(x) for x in row] + [Fraction(b)] for row, b in zip(M, rhs)]
    piv = []; r = 0
    for c in range(n + 1):
        pr = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if pr is None:
            continue
        if c == n:
            return None
        A[r], A[pr] = A[pr], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c); r += 1
    x = [Fraction(0)] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    return x

# ---------------------------------------------------------------- Hensel lifting (PROOF.md Proposition 4)
def hensel(polys, n, xbar, p, K, fixed=()):
    """Lift xbar (a zero mod p of every poly, integer polys) to a zero mod p^K, digit by digit, with the Jacobian
    frozen at xbar.  Coordinates in `fixed` are kept equal to their given integer values (not lifted).
    Returns ('lifted', x mod p^K) or ('obstructed', k, residual) where k is the first precision p^(k+1) that cannot be
    reached along this path.  If J(xbar) restricted to the free coordinates is injective mod p, the lift is unique, so
    'obstructed' PROVES that xbar has no p-adic lift with these fixed coordinates (Prop. 4)."""
    free = [i for i in range(n) if i not in fixed]
    J = jac(polys, n)
    Jb = [[ev(J[r][i], xbar, p) for i in free] for r in range(len(polys))]
    x = list(xbar)
    pk = p
    for k in range(1, K):
        vals = [ev(P, x, pk * p) for P in polys]
        if any(v % pk for v in vals):
            return ('bad-start', k, None)
        r = [(-(v // pk)) % p for v in vals]
        d = solve_mod(Jb, r, p)
        if d is None:
            return ('obstructed', k, r)
        for j, i in enumerate(free):
            x[i] = x[i] + pk * d[j]
        pk *= p
    return ('lifted', x, pk)

def injective_on(polys, n, xbar, p, fixed=()):
    free = [i for i in range(n) if i not in fixed]
    J = jac(polys, n)
    Jb = [[ev(J[r][i], xbar, p) for i in free] for r in range(len(polys))]
    return rank_mod(Jb, p) == len(free), rank_mod(Jb, p), len(free)

# ---------------------------------------------------------------- linear-in-beta structure
def split_linear(names, polys, beta_prefixes=('b', 'c')):
    """indices of alpha (others) and beta variables; checks every poly has degree <= 1 in the betas jointly"""
    bi = [i for i, v in enumerate(names) if v[0] in beta_prefixes]
    ai = [i for i in range(len(names)) if i not in bi]
    for P in polys:
        for k in P:
            assert sum(k[i] for i in bi) <= 1, 'not linear in beta'
    return ai, bi

def beta_system(names, polys, ai, bi, alpha, m=None):
    """for a numeric alpha (dict index->value), return (A, rhs) with A beta = rhs, exact or mod m"""
    A = [[0] * len(bi) for _ in polys]; rhs = [0] * len(polys)
    pos = {b: j for j, b in enumerate(bi)}
    for r, P in enumerate(polys):
        for k, c in P.items():
            t = c if m is None else coef_mod(c, m)
            for i in ai:
                if k[i]:
                    t = t * (alpha[i] ** k[i]) if m is None else t * pow(alpha[i], k[i], m) % m
            bj = next((b for b in bi if k[b]), None)
            if bj is None:
                rhs[r] = rhs[r] - t
                if m is not None: rhs[r] %= m
            else:
                A[r][pos[bj]] = A[r][pos[bj]] + t
                if m is not None: A[r][pos[bj]] %= m
    return A, rhs

def load_Q(prefix):
    """prefix.p1.ms + prefix.p2.ms -> (names, polys over Q)"""
    return to_Q(parse_ms(prefix + '.p1.ms'), parse_ms(prefix + '.p2.ms'))

def joint(systems):
    """systems = [(names, polys)] for the arrays of one skeleton: the alphas (a*) are shared, the betas get suffixes.
    Returns (names, polys)."""
    alln = []
    for l, (names, _) in enumerate(systems):
        for v in names:
            w = v if v[0] == 'a' else f"{v}_{l}"
            if w not in alln:
                alln.append(w)
    out = []
    for l, (names, polys) in enumerate(systems):
        mp = [alln.index(v if v[0] == 'a' else f"{v}_{l}") for v in names]
        for P in polys:
            Q = {}
            for k, c in P.items():
                e = [0] * len(alln)
                for i, x in enumerate(k):
                    e[mp[i]] += x
                Q[tuple(e)] = c
            out.append(Q)
    return alln, out
