# fc-gate route 2 (own code, independent of mixed.py's Decimal root and of x-pack's gen_rates/gen_sharp enclosures):
# exact rational interval arithmetic (outward dyadic rounding at 2^-PREC), logs by atanh series with remainder,
# exps by Taylor with remainder + squaring.  Mixed two-target rate fact of Theorem 1G in e-r45m1's normal form:
#     F(z) = sum_r n_r (r/m)^z - (W - H) - H (kc/m)^z  <= 0
# root tau (F(tau) = 0), saving b = 1 - tau to 60+ digits, largest grid N with the rate fact PROVED at
# z1 = 1 - (100 N + 2)/10^12 (e-r45m1's convention), and N + 1 certified FALSE (at z1(N+1) and at 1 - (N+1)/10^10).
# usage: python3 -I -B route2.py  (profiles built below from the counts; no input files)
from fractions import Fraction as Q
import json, sys
PREC = 330
S = 1 << PREC
def dn(x): return Q((x.numerator * S) // x.denominator, S)
def up(x): return -dn(-x)
class I:
    __slots__ = ('lo', 'hi')
    def __init__(s, lo, hi=None):
        s.lo = lo; s.hi = lo if hi is None else hi
        assert s.lo <= s.hi
    def __add__(a, b):
        b = b if isinstance(b, I) else I(Q(b)); return I(dn(a.lo + b.lo), up(a.hi + b.hi))
    def __sub__(a, b):
        b = b if isinstance(b, I) else I(Q(b)); return I(dn(a.lo - b.hi), up(a.hi - b.lo))
    def __mul__(a, b):
        b = b if isinstance(b, I) else I(Q(b))
        p = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi]; return I(dn(min(p)), up(max(p)))
    __rmul__ = __mul__
_LN = {}
def atanh_enc(t):          # 0 <= t <= 1/3 rational; returns interval for atanh(t)
    assert 0 <= t <= Q(1, 3)
    s = Q(0); p = t; j = 0
    while True:
        s += p / (2 * j + 1); j += 1; p *= t * t
        if p < Q(1, 1 << (PREC + 8)): break
    rem = p / ((2 * j + 1) * (1 - t * t))      # tail sum_{i>=j} t^(2i+1)/(2i+1) <= p/((2j+1)(1-t^2))
    return I(dn(s), up(s + rem))
def ln_int(n):             # ln n for integer n >= 1
    if n in _LN: return _LN[n]
    if n == 1: v = I(Q(0))
    else:
        ln2 = I(Q(2)) * atanh_enc(Q(1, 3))
        k = n.bit_length() - 1; y = Q(n, 1 << k)            # y in [1, 2)
        if y > Q(4, 3):                                    # keep t <= 1/7
            k += 1; y = y / 2                              # y in (2/3, 1]: ln y = -2 atanh((1-y)/(1+y))
            t = (1 - y) / (1 + y); v = ln2 * k - I(Q(2)) * atanh_enc(t)
        else:
            t = (y - 1) / (y + 1); v = ln2 * k + I(Q(2)) * atanh_enc(t)
    _LN[n] = v; return v
def exp_small(x):          # |x| <= 1/16 rational, enclosure of e^x
    assert abs(x) <= Q(1, 16)
    s = Q(0); t = Q(1); j = 0
    while True:
        s += t; j += 1; t = t * x / j
        if abs(t) < Q(1, 1 << (PREC + 8)): break
    rem = abs(t) * Q(16, 15)                              # |tail| <= |t| / (1 - |x|)
    return I(dn(s - rem), up(s + rem))
def exp_iv(v):             # e^v for an interval v (monotone)
    def one(x, lower):
        k = 0
        while abs(x) > Q(1, 16) * (1 << k): k += 1
        e = exp_small(x / (1 << k)); r = e.lo if lower else e.hi
        for _ in range(k): r = dn(r * r) if lower else up(r * r)
        return r
    return I(one(v.lo, True), one(v.hi, False))
def pow_ratio(r, m, z):    # (r/m)^z, z rational
    return exp_iv((ln_int(r) - ln_int(m)) * z)
def F(p, z):
    m = p['m']; acc = I(Q(0))
    for r, n in p['children'].items():
        acc = acc + pow_ratio(r, m, z) * n
    return acc - (p['W'] - p['H']) - pow_ratio(p['kc'], m, z) * p['H']
def holds(p, z): return F(p, z).hi <= 0          # rate fact proved at z
def fails(p, z): return F(p, z).lo > 0           # rate fact certified false at z
def root(p, digits=64):
    lo, hi = Q(97, 100), Q(1)                    # F(lo) > 0 > F(hi)
    assert fails(p, lo) and holds(p, hi)
    while hi - lo > Q(1, 10 ** digits):
        mid = (lo + hi) / 2; mid = Q(round(mid * 10 ** (digits + 4)), 10 ** (digits + 4))
        if fails(p, mid): lo = mid
        elif holds(p, mid): hi = mid
        else: raise SystemExit('interval too wide at %s' % float(mid))
    return lo, hi
def z1(N): return 1 - Q(100 * N + 2, 10 ** 12)
def scan(p, guess):
    N = guess
    while not holds(p, z1(N)): N -= 1
    while holds(p, z1(N + 1)): N += 1
    return N
def sdec(x, d=66):
    s = str(int(x * 10 ** d)); s = s.rjust(d + 1, '0'); return s[:-d] + '.' + s[-d:]
def mk(base, extra, W, H, kc, m=73):
    ch = dict(base)
    for r, n in extra.items(): ch[r] = ch.get(r, 0) + n
    return dict(children=dict(sorted(ch.items())), W=W, H=H, kc=kc, m=m)
base = {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920}   # #3945 PROFILE without the 137,952 idle 48s
W = 322560
P = {
  'full_B':  mk(base, {24: 2 * 57312}, W, 80640, 25),                 # design B full cut
  'half_B':  mk(base, {24: 80640, 48: 57312}, W, 80640, 49),          # design B half cut (control 1.7581520e-3)
  'full_pub': mk(base, {24: 2 * 80640}, W, 57312, 25),                # #3945 published fill, full cut (control 1.7307680e-3)
  'unif':    dict(children=dict(sorted({**base, 48: 137952}.items())), W=W, H=0, kc=25, m=73),   # control 1.7306632e-3
}
out = {}
for k, p in P.items():
    if k == 'unif':
        lo, hi = Q(97, 100), Q(1)
        while hi - lo > Q(1, 10 ** 64):
            mid = (lo + hi) / 2; mid = Q(round(mid * 10 ** 68), 10 ** 68)
            if F(p, mid).lo > 0: lo = mid
            elif F(p, mid).hi <= 0: hi = mid
            else: raise SystemExit('wide')
        out[k] = dict(b_lo=sdec(1 - hi), b_hi=sdec(1 - lo)); print(k, out[k], flush=True); continue
    lo, hi = root(p)
    b_lo, b_hi = 1 - hi, 1 - lo
    g = int(b_lo * 10 ** 10)
    N = scan(p, g)
    o = dict(profile=p, b_lo=sdec(b_lo), b_hi=sdec(b_hi), N=N, N_proved=holds(p, z1(N)),
             N1_refused_at_z1=fails(p, z1(N + 1)), N1_false_at_grid=fails(p, 1 - Q(N + 1, 10 ** 10)),
             root_minus_N_grid_steps=float(b_lo * 10 ** 10 - N),
             separated=max(p['children']) < p['kc'], rank_mass=sum(r * n for r, n in p['children'].items()),
             count=sum(p['children'].values()))
    out[k] = o; print(k, json.dumps({x: y for x, y in o.items() if x != 'profile'}), flush=True)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else '/dev/null', 'w'), indent=1)
