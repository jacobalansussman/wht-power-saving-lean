#!/usr/bin/env python3
"""hx-price (wht45/hubx, 20 Oct 2026): full-cut price of a counted design, by three routes, with controls.

Input: the JSON of tools/hx_heads.py (own head recount: paid_profile, idle_profile, literal_heads, W).
The full cut is built HERE from the counts (fc-reopt's pricing code is not used):
  children = paid profile, minus the literal heads (rank 48, removed: their arrays get target 73 - 48 = 25),
             with every other idle piece of rank > 24 split into pieces of 24 (48 -> 24 + 24, 42 -> 24 + 18, ...);
  targets  = {73: W - heads, 25: heads}.
Rule: tau is the root of sum_children n_r r^tau = sum_targets W_k k^tau; b = 1 - tau (complex side).
Routes:
  A  the moment-engine re-implementation tools/lever/pricer.py (its exact rational
     log_bounds / exp_bounds: Mercator logs, 24-term exp, 224-bit dyadic outward rounding), applied to the mixed form
     sum_r n_r (r/73)(73/r)^b  vs  W73 + W25 (25/73)(73/25)^b;  certified on the 1e-24 grid (b24 and b24 + 1e-24).
  B  tools/lever mixed.py (Decimal, 60 digits) on the children/targets built here, and lever_price.variants()'s
     'full_cut' row on (profile, W, heads, non) as that tool builds it.
  C  own interval arithmetic in exact rationals (atanh-series logs, argument-reduced Taylor exp, 2^-640 dyadic outward
     rounding) on the unnormalised form sum n_r r^tau - sum W_k k^tau: a 70-digit root by interval bisection, and the
     sign at every grid point that matters.
Grid facts certified by A and C: the largest N with N/1e10 < b (and N + 1 fails), and the Lean stack's rate point
convention (fc-rate NOTE section 2): rate point (100N + 2)/1e12 below the root, and (100(N+1) + 2)/1e12 above it.
usage: hx_price.py HEADS.json [--label L] [--controls]
"""
import sys, json, io, contextlib
from fractions import Fraction as F
from decimal import Decimal, getcontext

import os as _os
LEVER = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', '..', '..', 'tools', 'lever')   # the repository's tools/lever
LR = LEVER

# ---------------------------------------------------------------- route C: own exact interval arithmetic
BITS = 640
DY = 1 << BITS
def dn(x): return F(x.numerator * DY // x.denominator, DY)
def up(x): return -dn(-x)

_lncache = {}
def ln_iv(q):
    """[lo, hi] enclosing ln q for rational q >= 1: ln q = k ln 2 + 2 atanh((z-1)/(z+1)), z = q / 2^k in [1, 2)."""
    q = F(q)
    if q in _lncache: return _lncache[q]
    assert q >= 1
    # exact partial sums (no rounding inside, so lo is a true lower bound): recompute with exact powers
    def atanh2x(y, n):
        s = F(0); t = y; y2 = y * y
        for j in range(n):
            s += t / (2 * j + 1); t = t * y2
        return 2 * s, 2 * s + 2 * t / ((2 * n + 1) * (1 - y2))
    if 2 not in _lncache:
        lo2, hi2 = atanh2x(F(1, 3), 150)   # ln 2 = 2 atanh(1/3); (1/9)^150 ~ 1e-143
        _lncache[F(2)] = (dn(lo2), up(hi2))
    l2, h2 = _lncache[F(2)]
    k = 0; z = q
    while z >= 2: z /= 2; k += 1
    y = (z - 1) / (z + 1)
    lo, hi = atanh2x(y, 150) if y else (F(0), F(0))
    r = (dn(k * l2 + lo), up(k * h2 + hi))
    _lncache[q] = r
    return r

def exp_iv(lo, hi):
    """[l, u] enclosing exp([lo, hi]) for 0 <= lo <= hi <= 8: exp(x) = exp(x / 2^10)^(2^10), Taylor with remainder."""
    assert 0 <= lo <= hi <= 8
    R = 10
    def ser(x, upper):
        x = (up if upper else dn)(x / (1 << R))
        s = F(1); t = F(1)
        for n in range(1, 40):
            t = (up if upper else dn)(t * x / n); s += t
        if upper:                       # remainder <= t x/40 / (1 - x/41)
            s += up(t * x / 40 / (1 - x / 41))
        v = (up if upper else dn)(s)
        for _ in range(R): v = (up if upper else dn)(v * v)
        return v
    return ser(lo, False), ser(hi, True)

def f_iv(children, targets, tau):
    """interval of sum n_r r^tau - sum W_k k^tau (tau rational in (0.9, 1))."""
    lo = hi = F(0)
    for r, n in children.items():
        if r == 1: lo += n; hi += n; continue
        a, b = ln_iv(r); e, g = exp_iv(dn(tau * a), up(tau * b)); lo += n * e; hi += n * g
    for k, n in targets.items():
        a, b = ln_iv(k); e, g = exp_iv(dn(tau * a), up(tau * b)); lo -= n * g; hi -= n * e
    return lo, hi

def sign_C(children, targets, b):
    """-1 if certainly b < b* (f(1-b) < 0), +1 if certainly b > b*, 0 if undecided."""
    lo, hi = f_iv(children, targets, 1 - F(b))
    return -1 if hi < 0 else (1 if lo > 0 else 0)

def root_C(children, targets, lo=F(1, 1000), hi=F(3, 1000), iters=240):
    assert sign_C(children, targets, lo) == -1 and sign_C(children, targets, hi) == 1
    for _ in range(iters):
        mid = (lo + hi) / 2; s = sign_C(children, targets, mid)
        if s == 0: break
        if s < 0: lo = mid
        else: hi = mid
    return lo, hi

# ---------------------------------------------------------------- route A: pricer.py's exact enclosures, mixed form
sys.path.insert(0, LEVER)
import pricer                                     # m-anatomy's exact re-implementation of the hub engine
def mixed_moment_A(children, targets, b, m=73):
    """exact enclosures of  CH(b) = sum_r n_r (r/m)(m/r)^b  and  TG(b) = sum_k W_k (k/m)(m/k)^b."""
    def side(d):
        lo = hi = F(0)
        for r, n in d.items():
            if r == m: lo += n; hi += n; continue
            assert 0 < r < m
            l, u = pricer.log_bounds(F(m, r)); e, f = pricer.exp_bounds(b * l, b * u)
            p = F(r * n, m); lo += p * e; hi += p * f
        return pricer._floor(lo), pricer._ceil(hi)
    return side(children), side(targets)

def sign_A(children, targets, b):
    (cl, ch), (tl, th) = mixed_moment_A(children, targets, F(b))
    return -1 if ch < tl else (1 if cl > th else 0)

def grid24_A(children, targets, approx):
    """certify b24 = floor(b* 1e24)/1e24 from a Decimal approximation (re-centred if one step off)."""
    q = int(Decimal(approx) * Decimal(10) ** 24)
    for d in (0, -1, 1, -2, 2):
        b = F(q + d, 10 ** 24)
        if sign_A(children, targets, b) == -1 and sign_A(children, targets, b + F(1, 10 ** 24)) == 1: return b
    raise SystemExit('route A: 1e-24 grid point not certified')

# ---------------------------------------------------------------- route B: tools/lever mixed.py / lever_price.py
sys.path.insert(0, LR)
import mixed as MX
with contextlib.redirect_stdout(io.StringIO()):
    import lever_price as LP                      # its module body prints its controls; silenced here

# ---------------------------------------------------------------- build the full cut from counts
def full_cut(paid, idle, heads, W, cap=24):
    paid = {int(r): int(n) for r, n in paid.items()}; idle = {int(r): int(n) for r, n in idle.items()}
    local = dict(paid)
    for r, n in idle.items():
        local[r] -= n; assert local[r] >= 0
        if not local[r]: del local[r]
    assert max(local) <= cap, 'a non-idle paid block above the cap'
    ch = dict(local)
    idle2 = dict(idle); idle2[48] = idle2.get(48, 0) - heads; assert idle2[48] >= 0, 'heads exceed idle 48s'
    split = {}
    for r, n in idle2.items():
        if not n: continue
        rest = r
        while rest > 0:
            p = min(cap, rest); ch[p] = ch.get(p, 0) + n; rest -= p
        if r > cap: split[r] = n
    tg = {73: W - heads, 25: heads} if heads else {73: W}
    mass_before = sum(r * n for r, n in paid.items()); mass_after = sum(r * n for r, n in ch.items())
    assert mass_before - mass_after == 48 * heads, 'rank mass drops by exactly 48 per head'
    return {r: n for r, n in sorted(ch.items()) if n}, tg, dict(local=local, idle_after_heads=idle2, split_idle_pieces=split,
                                                                mass_before=mass_before, mass_after=mass_after)

def price_all(label, ch, tg):
    getcontext().prec = 60
    bB = MX.root(ch, tg)                                            # route B (Decimal 60, as the tool sets it)
    lo, hi = root_C(ch, tg)                                         # route C (own intervals)
    bC = (lo + hi) / 2
    b24 = grid24_A(ch, tg, str(bB))                                 # route A (1e-24 grid, certified)
    # grid facts
    N = int(lo * 10 ** 10)
    facts = {}
    for name, b in (('N/1e10', F(N, 10 ** 10)), ('(N+1)/1e10', F(N + 1, 10 ** 10)),
                    ('rate point (100N+2)/1e12', F(100 * N + 2, 10 ** 12)), ('rate point (100(N+1)+2)/1e12', F(100 * (N + 1) + 2, 10 ** 12)),
                    ('seam (100N+1)/1e12', F(100 * N + 1, 10 ** 12))):
        facts[name] = dict(A=sign_A(ch, tg, b), C=sign_C(ch, tg, b))
    ok = (facts['N/1e10'] == dict(A=-1, C=-1) and facts['(N+1)/1e10'] == dict(A=1, C=1)
          and facts['rate point (100N+2)/1e12'] == dict(A=-1, C=-1) and facts['rate point (100(N+1)+2)/1e12'] == dict(A=1, C=1))
    def d(x, p=70):
        getcontext().prec = p; return str(Decimal(x.numerator) / Decimal(x.denominator))
    rC = (d(lo), d(hi))
    return dict(label=label, children=ch, targets=tg, max_paid_rank=max(ch), min_target=min(tg),
                separated=max(ch) < min(tg),
                b_routeB_mixed60=str(bB), b_routeC_interval_lo=rC[0], b_routeC_interval_hi=rC[1],
                b_routeA_grid24=d(b24, 30), b_routeA_grid24_certified=True,
                B_vs_C_abs_diff=str(abs(Decimal(str(bB)) - Decimal(rC[0]))),
                grid_numerator_1e10=N, steps_above_N=str(round((Decimal(rC[0]) * Decimal(10) ** 10 - N), 6)),
                grid_signs_minus_below_plus_above=facts, grid_facts_ok=ok)

def design(label, paid, idle, heads, W):
    ch, tg, info = full_cut(paid, idle, heads, W)
    out = dict(label=label, W=W, heads=heads, construction=info)
    out['full_cut'] = price_all(label + ' full cut', ch, tg)
    prof = {int(r): int(n) for r, n in paid.items()}
    non = int(idle.get(48, idle.get('48', 0))) - heads
    getcontext().prec = 60
    with contextlib.redirect_stdout(io.StringIO()):
        v = LP.variants(label, prof, W, heads, non, quiet=True)
    getcontext().prec = 60
    out['lever_price_variants'] = {k: '{:.40e}'.format(v[k]['b']) for k in ('uniform', 'half_cut', 'full_cut')}
    out['lever_price_variants']['breakeven_full_heads'] = v['breakeven_full']
    out['uniform_hub_engine_grid24'] = pricer.dec(pricer.price_b(prof, 73, W), 30)   # pricer.price_b, 1e-24 certified
    return out

if __name__ == '__main__':
    res = []
    if '--controls' in sys.argv:
        P3945 = {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 48: 137952}
        res.append(design('design B (#3945 profile, 80,640 heads, 57,312 other 48s)', P3945, {48: 137952}, 80640, 322560))
        P2979 = {1: 1645680, 2: 470196, 3: 159468, 4: 71928, 5: 13716, 8: 5328, 48: 28190}
        res.append(design('design A (#2979 profile, 16,576 heads, 11,614 other 48s)', P2979, {48: 28190}, 16576, 66304))
    else:
        h = json.load(open(sys.argv[1])); label = sys.argv[sys.argv.index('--label') + 1] if '--label' in sys.argv else sys.argv[1]
        res.append(design(label, h['paid_profile'], h['idle_profile'], h['literal_heads'], h['W']))
        res[-1]['input'] = sys.argv[1]
    print(json.dumps(res, indent=1, default=str))
