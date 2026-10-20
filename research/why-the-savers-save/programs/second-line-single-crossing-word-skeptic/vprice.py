#!/usr/bin/env python3
"""v-an2-word: prices of the two-sided cases, own code (stdlib Decimal), same pricing model as an2-word/quotient.py
(the Lean fold hypothesis of Work/BlockAccounting/Moment.lean, as printed in RateU8X0HQx.lean and re-derived by
v-an-build/calc/vcalib.py):
    (2^40 - 1) * sum_r n_r (r/m)^(1-a) + W ((m-1)/m)^(1-a) + W (1/m)^(1-a) < 2^40 W ,
root a by bisection.  Profile: the U8 word's local per-call profile (words.json, U8_exit_frames(E1)) x 72 calls, plus
optional rank-1 crossing pieces, the rest of m*W - T*D filled with rank-48 idle pieces (an-ours' exact fill).
Cases:
  calibration     one-sided, m 73, D 320 / replica                      (must give 1.6647327e-3)
  lemma floor     two-sided, m 74, D 208; crossings as rank-1 pieces    (an2-word: 1.0545261e-3)
  lemma floor'    two-sided, m 74, D 208; crossings absorbed into the fill (no rank-1 pieces; most optimistic)
  trivial lift    two-sided, m 74, D 320; each X array climbs s and drops s as two rank-1 pieces (896 / replica)
  trivial lift'   two-sided, m 74, D 320; the 896 crossings absorbed into the fill (most optimistic)
  first order     a ~ T D / N, N = sum n r ln(m/r) (narrow mass) for each case, as a cross-check.
usage: python3 -I vprice.py OUT.json"""
import json, math, os, sys
from decimal import Decimal as Dd, getcontext
getcontext().prec = 60


def g(H, W, m, a, s=40):
    e = 1 - a; lm = Dd(m)
    tot = sum(Dd(n) * ((Dd(r) / lm).ln() * e).exp() for r, n in H.items())
    fill = ((Dd(m - 1) / lm).ln() * e).exp() + ((Dd(1) / lm).ln() * e).exp()
    two = Dd(2) ** s
    return (two - 1) * tot + Dd(W) * fill - two * Dd(W)


def root(H, W, m):
    lo, hi = Dd(0), Dd('0.05')
    assert g(H, W, m, lo) < 0 < g(H, W, m, hi)
    for _ in range(100):
        mid = (lo + hi) / 2
        if g(H, W, m, mid) < 0: lo = mid
        else: hi = mid
    return float(lo)


BC = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs') + '/')
wj = json.load(open(BC + 'words.json'))
prof = {int(k): n for k, n in wj['869/870/890_U8_exit_frames(E1)']['profile'].items()}
T, Wrep = 8, 896
W = T * Wrep


def build(m, Drep, rank1):
    H = {r: 72 * n for r, n in prof.items()}
    if rank1: H[1] = H.get(1, 0) + rank1
    rem = m * W - T * Drep - sum(r * n for r, n in H.items())
    assert rem >= 0
    q, rr = divmod(rem, 48)
    H[48] = H.get(48, 0) + q
    if rr: H[rr] = H.get(rr, 0) + 1
    assert m * W - sum(r * n for r, n in H.items()) == T * Drep
    return H


cases = [('calibration (one-sided, m 73, D 320)', 73, 320, 0),
         ('lemma floor (m 74, D 208, 9x112 rank-1 per replica)', 74, 208, T * 9 * 112),
         ("lemma floor' (m 74, D 208, crossings absorbed)", 74, 208, 0),
         ('trivial lift (m 74, D 320, 2 rank-1 per X array)', 74, 320, T * 448 * 2),
         ("trivial lift' (m 74, D 320, crossings absorbed)", 74, 320, 0)]
out = {}
for name, m, Drep, r1 in cases:
    H = build(m, Drep, r1)
    a = root(H, W, m)
    N = sum(n * r * math.log(m / r) for r, n in H.items())
    out[name] = dict(m=m, D_per_replica=Drep, rank1_pieces=r1, root=a, first_order_TD_over_N=T * Drep / N)
    print('%-58s root %.7e   first order T*D/N %.4e' % (name, a, T * Drep / N))
assert abs(out['calibration (one-sided, m 73, D 320)']['root'] - 1.6647327e-3) < 1e-10
hub = 1.6797550e-3
for k, v in out.items():
    v['below_hub_1.6797550e-3'] = v['root'] < hub
    v['below_one_sided_exact_fill'] = v['root'] < 1.6647327e-3 or k.startswith('calibration')
json.dump(out, open(sys.argv[1], 'w'), indent=1)
