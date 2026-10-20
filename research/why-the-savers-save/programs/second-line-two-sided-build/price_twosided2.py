#!/usr/bin/env python3
"""an-build: price of two-sided slack on U8X WITH the GB crossing copies (mintrans2_U8X.json), our fold (s = 40), same
idealisation as an-ours' exact-fill ceiling (idle remainder in rank-48 pieces; banks break even).  Per word: 72 calls,
of which 32 are GB (4 per run: stages 1, 3, 5, 7; Gates.lean liftQ1/3/5/7 GB); each GB call pays 112 rank-1 crossing
copies (<s> -> 0).  Y-class residuals only on Y-side children (half the idle cover); S-class (switching) residuals only on
banks (no data array can switch s); the S slots are then clean in GF calls (banks restore), so GF needs no crossing.
Bank models for the switching residuals: strict (one per bank), per-run (one per bank and run), relative (plain width)."""
import sys, os, json
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
from anatomy_ours import price, rank, narrow, add, LOCAL
HERE = os.path.dirname(os.path.abspath(__file__))
WD, CALLS, GB, CHILD = 7168, 72, 32, 387072
mt = json.load(open(os.path.join(HERE, 'mintrans2_U8X.json')))
loc = {int(k): v for k, v in mt['profile_GF'].items()}
Yw, Sw, Sn, Xc = mt['width_sum'].get('Y', 0), mt['width_sum'].get('S', 0), mt['switching'], mt['crossing_copies_per_GB_call']
Xw = mt['width_sum'].get('X', 0)
m = 74
D_expect = m * WD - CALLS * mt['rank_minus_widths'] - CHILD - GB * Xc
def fill(rem, piece=48):
    q, r = divmod(rem, piece); return {piece: q, r: 1} if r else ({piece: q} if q else {})
res = {}
for bm in ('strict', 'per-run', 'relative'):
    side = CHILD // 2
    hostY = min(side, CALLS * Yw); hostX = min(side, CALLS * Xw)
    overflow = CALLS * (Yw + Xw) - hostY - hostX
    idle = CHILD - hostY - hostX
    S_n, S_tot = CALLS * Sn, CALLS * Sw; w = Sw // Sn
    pads = {}; nb = 0
    if bm == 'strict':
        nb = S_n; pads[m - w] = S_n
        q, r = divmod(overflow, m); nb += q + (1 if r else 0)
        if r: pads[m - r] = pads.get(m - r, 0) + 1
    elif bm == 'per-run':
        nbs = -(-S_n // 8); free = nbs * (m - 8 * w); put = min(free, overflow); rest = overflow - put
        base, extra = divmod(free - put, nbs)
        if base: pads[base] = nbs - extra
        if extra: pads[base + 1] = pads.get(base + 1, 0) + extra
        q, r = divmod(rest, m); nb = nbs + q + (1 if r else 0)
        if r: pads[m - r] = pads.get(m - r, 0) + 1
    else:
        q, r = divmod(overflow + S_tot, m); nb = q + (1 if r else 0)
        if r: pads[m - r] = 1
    H = add({kk: CALLS * vv for kk, vv in loc.items()}, {1: GB * Xc}, fill(idle), {kk: c for kk, c in pads.items() if kk > 0})
    W = WD + nb; D = m * W - rank(H)
    assert D == D_expect, (D, D_expect)
    a = float(price(H, W, m=m))
    res[bm] = dict(W=W, banks=nb, D=D, N=round(narrow(H, m), 1), a=a)
    print('two-sided U8X with GB crossing copies, bank model %-9s: W %6d banks %5d D %5d N %10.1f  a = %.7e  (x %.3f over 1.6647327e-3, x %.3f over the hub record 1.6797550e-3)'
          % (bm, W, nb, D, narrow(H, m), a, a / 1.6647327e-3, a / 1.6797550e-3))
json.dump(dict(D_word=D_expect, models=res, inputs=dict(Y=Yw, S=Sw, switching=Sn, crossing_per_GB=Xc)), open(os.path.join(HERE, 'price_twosided2.json'), 'w'), indent=1)
