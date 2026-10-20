#!/usr/bin/env python3
"""cx-mirror / price.py (own code, stdlib Decimal).  What a mirrored word may cost and still beat the record.

Model: the Lean fold hypothesis (Work/BlockAccounting/Moment.lean, as used by an-ours, an2-word, v-an2-word):
    (2^40 - 1) * sum_r n_r (r/m)^(1-a) + W ((m-1)/m)^(1-a) + W (1/m)^(1-a) < 2^40 W ,   root a by bisection;
profile = the U8 word's per-call local profile (fsp-price/data/words.json, U8_exit_frames(E1)) x 72 calls, plus the
extra centre pieces of a two-sided word, the rest of m W - T D filled with rank-48 idle pieces (an-ours' exact fill).
Two-sided at m = 74 (one new local coordinate s): each data array pays 72 of 74, so D per replica = 1792 - 9 C,
C = the whole centre charge per call (64 for the U-centre + everything the second line adds).
Cases: calibration (one-sided, m 73, D 320 -> 1.6647327e-3); two-sided with the V-centre as a second P112 code
(8 pieces of rank 8 per call, C = 128: the first-order ceiling); break-even C against our exact fill and the hub,
with the extra charge as rank-8 pieces (like a centre) and, more generously, absorbed into the fill."""
import json, os, math
from decimal import Decimal as Dd, getcontext
getcontext().prec = 50
BC = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs') + '/')
wj = json.load(open(BC + 'words.json'))
prof = {int(k): n for k, n in wj['869/870/890_U8_exit_frames(E1)']['profile'].items()}
T, Wrep = 8, 896; W = T * Wrep; HUB = 1.6797550e-3; OURS = 1.6647327e-3

def g(H, m, a, s=40):
    e = 1 - a; lm = Dd(m)
    tot = sum(Dd(n) * ((Dd(r) / lm).ln() * e).exp() for r, n in H.items())
    fill = ((Dd(m - 1) / lm).ln() * e).exp() + ((Dd(1) / lm).ln() * e).exp()
    two = Dd(2) ** s
    return (two - 1) * tot + Dd(W) * fill - two * Dd(W)

def root(H, m):
    lo, hi = Dd(0), Dd('0.05')
    assert g(H, m, lo) < 0 < g(H, m, hi)
    for _ in range(90):
        mid = (lo + hi) / 2
        if g(H, m, mid) < 0: lo = mid
        else: hi = mid
    return float(lo)

def build(m, Drep, extra):
    """extra: {rank: pieces per call} added for each of the 72 calls"""
    H = {r: 72 * n for r, n in prof.items()}
    for r, n in extra.items(): H[r] = H.get(r, 0) + 72 * n
    rem = m * W - T * Drep - sum(r * n for r, n in H.items())
    assert rem >= 0, rem
    q, rr = divmod(rem, 48); H[48] = H.get(48, 0) + q
    if rr: H[rr] = H.get(rr, 0) + 1
    assert m * W - sum(r * n for r, n in H.items()) == T * Drep
    return H

out = {}
a0 = root(build(73, 320, {}), 73); out['calibration'] = a0
print('calibration (one-sided, m 73, D 320): %.7e' % a0); assert abs(a0 - OURS) < 1e-10
def two(C, as_centre=True):
    Drep = 1792 - 9 * C
    extra = {8: (C - 64) // 8} if as_centre and C - 64 >= 8 else {}
    if as_centre and (C - 64) % 8: extra[(C - 64) % 8] = extra.get((C - 64) % 8, 0) + 1
    return root(build(74, Drep, extra), 74), Drep
a, Dr = two(128); out['ceiling_C128'] = dict(D_per_replica=Dr, a=a)
print('two-sided, V-centre = a second P112 code (C = 128 per call), m 74, D %d: %.7e  (x%.3f of the hub)' % (Dr, a, a / HUB))
for mode in (True, False):
    for target, name in ((OURS, 'our exact fill 1.6647327e-3'), (HUB, 'the hub 1.6797550e-3')):
        lo, hi = 64, 199
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if two(mid, mode)[0] > target: lo = mid
            else: hi = mid
        key = 'breakeven_vs_%s_%s' % ('ours' if target == OURS else 'hub', 'centre_pieces' if mode else 'absorbed')
        out[key] = dict(C_max=lo, D_per_replica=1792 - 9 * lo, extra_over_V_centre_per_call=lo - 128)
        print('largest C per call that still beats %s (%s): %d  -> D %d per replica; room beyond a P112-sized V-centre: %d per call'
              % (name, 'extra as rank-8 centre pieces' if mode else 'extra absorbed into the fill', lo, 1792 - 9 * lo, lo - 128))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'price.json'), 'w'), indent=1)
