#!/usr/bin/env python3
"""g-break: what-if prices on #942 for the breaks found (m-anatomy pricer; hosting models as u-helpers-hub).
 pair14: pair-row decomposition, 13 helpers/block (two copies) -> 182 plane helpers; 210 mean plane helpers removed;
         penalised: every line-visiting chain gets one extra run (21 line units/block: one extra rank-1 child whose width comes out of the idle stock).
 mscratch56: the 56 cold 0->4 scratch helpers (runs 4 x rank 1) retired; 56 idle M copies take their plane role,
         penalised: one extra run split in each M-copy chain (one extra rank-1 child whose width comes out of the idle stock).
WHAT-IF prices, not constructions.  Usage: python3 -I price_breaks.py"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('PRICER_DIR', '.'))
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as F
from collections import Counter
cen = json.load(open(os.path.join(here, '..', '..', 'helpers-census', 'u8_census.json')))
P942 = {1: 848490, 2: 112554, 3: 148518, 4: 53946, 5: 23976, 8: 2664, 48: 11086}; M, W = 73, 33152
b0 = price_b(P942, M, W)
PR = {int(r): x for r, x in cen['plane_runs'].items()}
def build(removed_runs_global, holes, splits2):
    out = {}
    freed = sum(r*x for r, x in removed_runs_global.items()) - splits2
    for mode in ('optimistic', 'hole'):
        G = Counter(P942)
        for r, x in removed_runs_global.items(): G[r] -= x
        G[1] += splits2   # extra rank-1 runs, their width taken from the idle stock
        if mode == 'optimistic':
            G[48] -= 11086; G[1] -= 6; hosted = 37*9*3778 - freed; resd = Counter()
            for cap, nb in ((8, 672*37), (16, 448*37), (32, 448*37), (48, 448*37)):
                full = min(nb, hosted//cap); hosted -= full*cap; left = nb-full
                if left and hosted: resd[cap-hosted] += 1; left -= 1; hosted = 0
                resd[cap] += left
            for r, x in resd.items(): G[r] += x
        else:
            hs = Counter(holes); tot = sum(w*c for w, c in hs.items()); nh = sum(hs.values()); f2 = tot - splits2
            lo = f2//nh; hi = f2 - lo*nh; G[lo] += nh-hi; G[lo+1] += hi
        G = {r: x for r, x in G.items() if x}
        assert sum(r*x for r, x in G.items()) == M*W-11840 and min(G.values()) > 0, G
        b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
        out[mode] = dict(b=dec(b, 13), b_pct=round(float((b/b0-1)*100), 3), kappa=dec(k, 16))
    return out
res = {}
N = 210
d = {r: round(F(333*N*x, 392)) for r, x in PR.items()}; fr = sum(r*x for r, x in d.items())
lo = fr//(333*N); hi = fr - lo*333*N
holes = Counter({lo: 333*N-hi}); holes[lo+1] += hi
res['pair14_gross'] = build(d, holes, 0)
res['pair14_penalised_294_splits_per_call'] = build(d, holes, 294*333)
d56 = {1: 56*4*333}
res['mscratch56_gross'] = build(d56, {4: 56*333}, 0)
res['mscratch56_penalised_56_splits_per_call'] = build(d56, {4: 56*333}, 56*333)
for k, x in res.items(): print(k, json.dumps(x))
json.dump(res, open(os.path.join(here, 'price_breaks.json'), 'w'), indent=1)
