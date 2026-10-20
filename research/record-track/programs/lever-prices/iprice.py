#!/usr/bin/env python3
"""lg-g1 (gate G1): independent rigorous re-pricing of the mixed-target moment root.
Own code: exact rational interval arithmetic (Fractions, outward rounding to 2^-320), logs from atanh series
with tail bounds, exp from Taylor with tail bound. No Decimal, no float in any decision.
Price: b = 1 - tau, tau the root of  sum_children n_r r^tau = sum_arrays k_w^tau  (hub moment root when all k_w = m).
F(b) = sum_r n_r r^(1-b) - sum_k W_k k^(1-b); b* = sup{b : F(b) < 0}; certified by F(lo) < 0 < F(hi).
Profile read from the hub's own SPEC.json (expected/global_profile, expected/live) given as argv[1]."""
import sys, json
from fractions import Fraction as Q
P = 320; S = Q(1, 2**P)
def dn(x): n = x.numerator * 2**P // x.denominator; return Q(n, 2**P)
def up(x): n = -((-x.numerator * 2**P) // x.denominator); return Q(n, 2**P)
class I:
    __slots__ = ('lo', 'hi')
    def __init__(s, lo, hi=None):
        if hi is None: hi = lo
        assert lo <= hi; s.lo, s.hi = dn(Q(lo)), up(Q(hi))
    def __add__(s, o): o = o if isinstance(o, I) else I(o); return I(s.lo + o.lo, s.hi + o.hi)
    __radd__ = __add__
    def __neg__(s): return I(-s.hi, -s.lo)
    def __sub__(s, o): o = o if isinstance(o, I) else I(o); return I(s.lo - o.hi, s.hi - o.lo)
    def __rsub__(s, o): return I(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, I) else I(o); c = [s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi]; return I(min(c), max(c))
    __rmul__ = __mul__
def atanh(x):  # x rational, |x| <= 1/3; sum x^(2k+1)/(2k+1) for k < 80, tail <= |x|^161/(1-x^2)
    x = Q(x); s = Q(0); p = x
    for k in range(80): s += p / (2*k + 1); p *= x * x
    t = abs(x)**161 / (1 - x*x)
    return I(s - t, s + t)
LN2 = 2 * atanh(Q(1, 3)); LN3 = LN2 + 2 * atanh(Q(1, 5)); LN5 = 2 * LN2 + 2 * atanh(Q(1, 9))
LN7 = 3 * LN2 + 2 * atanh(Q(-1, 15)); LN73 = 3 * LN2 + 2 * LN3 + 2 * atanh(Q(1, 145))
def ln(n):
    out = I(0); m = n
    for p, L in ((2, LN2), (3, LN3), (5, LN5), (7, LN7), (73, LN73)):
        while m % p == 0: out = out + L; m //= p
    assert m == 1, n
    return out
def expi(x):  # x interval with |x| < 1/2: Taylor to 40 terms on a 2^-400 grid; error <= 41*2*2^-400 + 2|x|^41/41!
    import math
    G = 2**400
    def ser(v):
        s = Q(0); t = Q(1)
        for k in range(1, 41):
            s += t; n = t * v / k; t = Q((n.numerator * G) // n.denominator, G)
        return s, abs(v)**41 / Q(math.factorial(41)) * 2 + Q(82, G)
    a, ea = ser(x.lo); b, eb = ser(x.hi)
    return I(a - ea, b + eb)  # exp increasing
LNC = {}
def pw(k, tau):  # k^tau for tau in [0,1]: k * exp(-(1-tau) ln k)
    if k == 1: return I(1)
    if k not in LNC: LNC[k] = ln(k)
    return k * expi(-(1 - tau) * LNC[k])
def F(children, targets, b):
    tau = I(1 - Q(b))
    s = I(0)
    for r, n in children.items():
        if n: s = s + n * pw(r, tau)
    for k, n in targets.items():
        if n: s = s - n * pw(k, tau)
    return s
def root(children, targets, lo=Q(0), hi=Q(1, 100), width=Q(1, 10**16)):
    assert F(children, targets, lo).hi < 0, 'F(lo) not < 0'
    assert F(children, targets, hi).lo > 0, 'F(hi) not > 0'
    while hi - lo > width:
        mid = (lo + hi) / 2; v = F(children, targets, mid)
        if v.hi < 0: lo = mid
        elif v.lo > 0: hi = mid
        else: raise RuntimeError('undecided at %s' % mid)
    return lo, hi
def sep(children, targets):
    return max(r for r, n in children.items() if n) < min(k for k, n in targets.items() if n)
def add(d, r, n): d = dict(d); d[r] = d.get(r, 0) + n; return d
def fmt(x, d=13): return ('%.' + str(d) + 'e') % (x.numerator / x.denominator) if False else '{:.13e}'.format(float(x))
def dec(x, digits=16):  # exact decimal string of a Fraction truncated
    n = x.numerator * 10**(digits) // x.denominator; s = str(n).rjust(digits + 1, '0'); return '0.' + s[1:] if s[0] == '0' else s
