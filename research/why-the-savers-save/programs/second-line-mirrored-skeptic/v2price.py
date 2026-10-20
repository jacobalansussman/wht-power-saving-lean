#!/usr/bin/env python3
"""v-cx-mirror / v2price.py (own code, stdlib Decimal).  Independent re-implementation of the Lean fold hypothesis
pricing used by an-ours / an2-word / v-an2-word / cx-mirror:
   (2^40 - 1) * sum_r n_r (r/m)^(1-a) + W ((m-1)/m)^(1-a) + W (1/m)^(1-a) < 2^40 W ,  a = sup by bisection.
Profile: the U8 local profile per call (words.json, read as data) x 72 calls, extra pieces, rest of m W - 8 D at
rank 48 (+ one remainder piece).  I solve with a secant/bisection hybrid on my own function, not the claim's code.
Checks: calibration 1.6647327e-3 (m 73, D 320), ceiling at C = 128 (m 74, D = 1792 - 9 C), the break-even C against
the hub 1.6797550e-3, and two sensitivity cases the claim did not state: (i) the extra charge of C - 128 as rank-1
pieces (drops of one level each, the cheapest kind of drop), (ii) the room if the V-centre itself cost less (C_V < 64)."""
import json, os
from decimal import Decimal as X, getcontext
getcontext().prec = 60
BC = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs') + '/')
prof = {int(k): int(n) for k, n in json.load(open(BC + 'words.json'))['869/870/890_U8_exit_frames(E1)']['profile'].items()}
NT, WR = 8, 896; W = NT * WR; HUB = X('1.6797550e-3')
TWO40 = X(2) ** 40
def pw(r, m, e): return (X(r) / X(m)).ln() * e
def F(H, m, a):
    e = 1 - a
    s = sum(X(c) * pw(r, m, e).exp() for r, c in H.items())
    return (TWO40 - 1) * s + X(W) * (pw(m - 1, m, e).exp() + pw(1, m, e).exp()) - TWO40 * X(W)
def solve(H, m):
    lo, hi = X(0), X('0.02')
    assert F(H, m, lo) < 0 < F(H, m, hi)
    for _ in range(80):
        mid = (lo + hi) / 2
        (lo, hi) = (mid, hi) if F(H, m, mid) < 0 else (lo, mid)
    return lo
def hist(m, D, extra):
    H = {}
    for r, c in prof.items(): H[r] = H.get(r, 0) + 72 * c
    for r, c in extra.items(): H[r] = H.get(r, 0) + 72 * c
    rest = m * W - NT * D - sum(r * c for r, c in H.items())
    assert rest >= 0
    q, rr = divmod(rest, 48)
    H[48] = H.get(48, 0) + q
    if rr: H[rr] = H.get(rr, 0) + 1
    return H
res = {}
a = solve(hist(73, 320, {}), 73); res['calib'] = float(a); print('calibration m73 D320: %.7e' % a)
a = solve(hist(74, 1792 - 9 * 128, {8: 8}), 74); res['ceiling'] = float(a); print('ceiling C=128 (8 extra rank-8 pieces per call), D %d: %.7e  x%.3f hub' % (1792 - 9 * 128, a, a / HUB))
def breakeven(kind):
    best = None
    for C in range(128, 200):
        if kind == 'rank8':
            ex = {8: 8 + (C - 128) // 8}
            if (C - 128) % 8: ex[(C - 128) % 8] = 1
        elif kind == 'rank1':
            ex = {8: 8, 1: C - 128}
        else:
            ex = {8: 8}
        p = solve(hist(74, 1792 - 9 * C, ex), 74)
        if p > HUB: best = (C, float(p))
        else: break
    return best
for kind in ('rank8', 'rank1', 'absorbed'):
    b = breakeven(kind); res['breakeven_' + kind] = b
    print('largest C beating the hub, extra as %s: C = %d (room %d per call beyond two P112 centres), price there %.7e' % (kind, b[0], b[0] - 128, b[1]))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v2price.json'), 'w'), indent=1)
