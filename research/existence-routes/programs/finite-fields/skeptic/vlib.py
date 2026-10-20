"""vlib.py -- v-ex2-modp (skeptic): an independent small exact toolkit.  Written from scratch; does NOT import
ex2-modp's zsys.py or lift_cert.py.  Polynomials are dicts {exponent tuple: coefficient}."""
from fractions import Fraction
import re

P1, P2 = 1073741789, 1073741741

def tokenize_term(term):
    """split a product term into factors, handling a leading sign"""
    return [f for f in term.split('*') if f != '']

def parse_poly_str(s, names, mod=None):
    """parse '3*a1*a1*b2+-5*a2^3+7' style (msolve or Singular).  Repeated factors ACCUMULATE exponents.
    Coefficients: integers or a/b fractions.  Returns dict exps->Fraction (or int mod `mod`)."""
    idx = {v: i for i, v in enumerate(names)}
    n = len(names)
    s = s.replace(' ', '').replace('\n', '').replace('\t', '')
    if s in ('', '0'):
        return {}
    # split on + and - at top level (no parentheses expected except Singular '(..)' coefficients)
    s = s.replace('(', '').replace(')', '')
    terms = []
    cur = ''
    for i, ch in enumerate(s):
        if ch in '+-' and i > 0 and s[i - 1] not in '^*/+-':
            terms.append(cur); cur = ch
        else:
            cur += ch
    terms.append(cur)
    P = {}
    for t in terms:
        if t in ('', '+'):
            continue
        sign = 1
        while t and t[0] in '+-':
            if t[0] == '-':
                sign = -sign
            t = t[1:]
        c = Fraction(sign)
        e = [0] * n
        for f in tokenize_term(t):
            if re.fullmatch(r'\d+', f):
                c *= int(f)
            elif re.fullmatch(r'\d+/\d+', f):
                a, b = f.split('/'); c *= Fraction(int(a), int(b))
            else:
                if '^' in f:
                    v, k = f.split('^')
                    e[idx[v]] += int(k)
                else:
                    e[idx[f]] += 1
        k = tuple(e)
        P[k] = P.get(k, 0) + c
    out = {}
    for k, v in P.items():
        if mod is not None:
            v = v.numerator * pow(v.denominator, -1, mod) % mod
        if v:
            out[k] = v
    return out

def read_ms(path):
    """msolve input file -> (names, p, [poly mod p]).  Polynomials are separated by ',' at line ends."""
    L = open(path).read().split('\n')
    names = [v.strip() for v in L[0].split(',') if v.strip()]
    p = int(L[1].strip())
    body = '\n'.join(L[2:])
    polys = []
    for chunk in body.split(','):
        if chunk.strip():
            polys.append(parse_poly_str(chunk, names, p))
    return names, p, polys

def ratrec(a, m):
    """my own rational reconstruction: r/s with |r|, |s| <= floor(sqrt(m/2)), r = a s (mod m)"""
    a %= m
    import math
    B = math.isqrt(m // 2)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B:
        return None
    f = Fraction(r1, s1)
    if (f.numerator - a * f.denominator) % m:
        return None
    return f

def to_Q(sysA, sysB):
    nA, pA, PA = sysA; nB, pB, PB = sysB
    assert nA == nB, 'variable lists differ'
    assert len(PA) == len(PB), 'row counts differ'
    out = []; maxnum = 0; maxden = 1
    for X, Y in zip(PA, PB):
        assert set(X) == set(Y), 'monomial supports differ between primes'
        Q = {}
        for k, v in X.items():
            r = ratrec(v, pA)
            assert r is not None, 'no reconstruction'
            assert (r.numerator - r.denominator * Y[k]) % pB == 0, 'primes disagree'
            Q[k] = r
            maxnum = max(maxnum, abs(r.numerator)); maxden = max(maxden, r.denominator)
        out.append(Q)
    return nA, out, (maxnum, maxden)

def mul(P, Q):
    out = {}
    for k1, c1 in P.items():
        for k2, c2 in Q.items():
            k = tuple(a + b for a, b in zip(k1, k2))
            out[k] = out.get(k, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}

def add_into(T, P):
    for k, v in P.items():
        w = T.get(k, 0) + v
        if w:
            T[k] = w
        else:
            T.pop(k, None)

def ev(P, x):
    tot = 0
    for k, c in P.items():
        t = c
        for xi, e in zip(x, k):
            if e:
                t = t * xi ** e
        tot = tot + t
    return tot

def ev_mod(P, x, p):
    tot = 0
    for k, c in P.items():
        t = c % p if not isinstance(c, Fraction) else c.numerator * pow(c.denominator, -1, p) % p
        for xi, e in zip(x, k):
            if e:
                t = t * pow(xi, e, p) % p
        tot = (tot + t) % p
    return tot
