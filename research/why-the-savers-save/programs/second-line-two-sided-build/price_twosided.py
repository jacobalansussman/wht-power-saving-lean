#!/usr/bin/env python3
"""an-build: price of two-sided slack on the U8X local word (with source copies), our fold (an-ours/anatomy_ours.price,
s = 40), in the same idealisation as an-ours' exact-fill ceiling (idle remainder in rank-48 pieces; banks break even).
Per word (T = 8 replicas, 72 calls): data arrays 7168; children (idle climbs) 387,072 rank, half on the X side, half on Y.
Two-sided (m = 74): X arrays never hold s, Y arrays hold s from time 0, so
  - Y-class residuals (hold s from the start; mintrans_U8X.json) can be hosted on Y-side children only;
  - X-class residuals on X-side children only (there are none on U8X);
  - S-class (switching) residuals on no data array; banks only.
Bank models for S residuals: strict (one switch per bank: a bank's frame gains s once), per-run (one switch per bank and
run, if the replicas' s-lines differ), relative (a switch is just one more unit of width).
Everything else as in host_price.py: hist = 72 x local profile + idle remainders + bank pads; W = 7168 + banks."""
import sys, os, json, math
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
from anatomy_ours import price, rank, narrow, add, LOCAL
HERE = os.path.dirname(os.path.abspath(__file__))
WD, CALLS, CHILD = 7168, 72, 387072
def fill_pieces(rem, piece=48):
    q, r = divmod(rem, piece); return {piece: q, r: 1} if r else ({piece: q} if q else {})
def model(local, m, Yw, Xw, Sw, Sn, bankmode, D_expect):
    """Yw, Xw, Sw: residual width per call by class; Sn: switching residuals per call."""
    side = CHILD // 2
    hostY = min(side, CALLS * Yw); hostX = min(side, CALLS * Xw)
    overflow = CALLS * (Yw + Xw) - hostY - hostX
    idle = CHILD - hostY - hostX
    pads = {}; nb = 0
    S_tot = CALLS * Sw; S_n = CALLS * Sn
    if bankmode == 'strict':          # one switching residual per bank; pad = m - w; overflow in full banks
        w = Sw // Sn; nb += S_n; pads[m - w] = pads.get(m - w, 0) + S_n
        q, r = divmod(overflow, m); nb += q + (1 if r else 0)
        if r: pads[m - r] = pads.get(m - r, 0) + 1
    elif bankmode == 'per-run':       # one switch per bank and run: 8 switches per bank, rest filled with overflow
        w = Sw // Sn; per = 8; nbs = -(-S_n // per); nb += nbs
        free = nbs * (m - per * w)                    # room left in the switch banks
        put = min(free, overflow); rest = overflow - put
        padtot = free - put
        base, extra = divmod(padtot, nbs)
        if base: pads[base] = pads.get(base, 0) + nbs - extra
        if extra: pads[base + 1] = pads.get(base + 1, 0) + extra
        q, r = divmod(rest, m); nb += q + (1 if r else 0)
        if r: pads[m - r] = pads.get(m - r, 0) + 1
    elif bankmode == 'relative':      # a switch is plain width: everything in full banks
        tot = overflow + S_tot; q, r = divmod(tot, m); nb += q + (1 if r else 0)
        if r: pads[m - r] = pads.get(m - r, 0) + 1
    H = {k: CALLS * n for k, n in local.items()}
    H = add(H, fill_pieces(idle), {k: c for k, c in pads.items() if k > 0})
    W = WD + nb
    D = m * W - rank(H)
    assert D == D_expect, (D, D_expect)
    a = float(price(H, W, m=m))
    return dict(m=m, W=W, banks=nb, D=D, N=round(narrow(H, m), 1), a=a, overflow=overflow, idle_left=idle)
out = {}
# calibration: one-sided U8 exact fill through the same function (all residuals hostable, both sides)
U8 = LOCAL['U8E1']
cal = model(U8, 73, 3778 // 2, 3778 - 3778 // 2, 0, 1, 'relative', 2560)
print('calibration, one-sided U8 exact fill (residuals 3778/call, either side):  a = %.7e  (an-ours 1.6647327e-3)' % cal['a'], cal)
assert abs(cal['a'] - 1.6647327e-3) < 2e-10
out['calibration'] = cal
# one-sided U8X (with the 112 copies) for reference: the copies' price
U8X = {1: 2660, 2: 338, 3: 446, 4: 162, 5: 72, 8: 120}
ref = model(U8X, 73, (3778 + 1008) // 2, 3778 + 1008 - (3778 + 1008) // 2, 0, 1, 'relative', 2560)
print('one-sided U8X (with copies), exact fill:                                  a = %.7e' % ref['a'])
out['one_sided_U8X'] = ref
mt = json.load(open(os.path.join(HERE, 'mintrans_U8X.json')))
loc2 = {int(k): v for k, v in mt['profile'].items()}
Yw, Sw, Sn = mt['width_sum'].get('Y', 0), mt['width_sum'].get('S', 0), mt['switching']
Xw = mt['width_sum'].get('X', 0)
D2 = 74 * WD - CALLS * mt['rank_minus_widths'] - CHILD
print('two-sided on U8X: Y width/call', Yw, ' X', Xw, ' S', Sw, '(', Sn, 'switching)  D per word', D2)
for bm in ('strict', 'per-run', 'relative'):
    r = model(loc2, 74, Yw, Xw, Sw, Sn, bm, D2)
    out['two_sided_' + bm] = r
    print('  bank model %-9s: W %6d banks %5d  D %5d  N %10.1f  a = %.7e  (x %.3f over the one-sided U8 ceiling)' %
          (bm, r['W'], r['banks'], r['D'], r['N'], r['a'], r['a'] / cal['a']))
json.dump(out, open(os.path.join(HERE, 'price_twosided.json'), 'w'), indent=1)
