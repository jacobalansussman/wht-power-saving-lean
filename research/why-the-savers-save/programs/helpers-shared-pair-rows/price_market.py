#!/usr/bin/env python3
"""d-pair-market: WHAT-IF price of the model count 448 (56 fewer helpers per call than U8's 504) on #942.
(a) 56 typical helpers (mean U8 helper runs x 56/504), as u-helpers-hub/price_removals.py;
(b) 56 of the 112 non-plane (M-copy) helpers, at their class-mean runs (half of {1:384,2:64,3:32,4:24}), the class
    whose role (own-block second halves) the model shares; holes by width at half of {8:56,5:32,4:24}.
Hosting: optimistic (re-pack relaxation) and hole, as in price_removals.py.  Control reproduced first.
Usage: python3 -I price_market.py -> price_market.json"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('PRICER_DIR', '.'))
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as F
from collections import Counter
cen = json.load(open(os.path.join(here, '..', 'helpers-census', 'u8_census.json')))
P942 = {1: 848490, 2: 112554, 3: 148518, 4: 53946, 5: 23976, 8: 2664, 48: 11086}; M, W = 73, 33152
b0 = price_b(P942, M, W); k0 = kappa_from(b0, C_PR589)[0]
assert b0 == F(1664727258740678855772, 10**24) and k0 == F(1661960547713947, 10**18)
def go(runs_global, holes):
    res = {}
    for mode in ('optimistic', 'hole'):
        G = Counter(P942)
        for r, x in runs_global.items(): G[r] -= x
        freed = sum(r*x for r, x in runs_global.items())
        if mode == 'optimistic':
            G[48] -= 11086; G[1] -= 6; hosted = 37*9*3778 - freed; resd = Counter()
            for cap, nb in ((8, 672*37), (16, 448*37), (32, 448*37), (48, 448*37)):
                full = min(nb, hosted//cap); hosted -= full*cap; left = nb-full
                if left and hosted: resd[cap-hosted] += 1; left -= 1; hosted = 0
                resd[cap] += left
            for r, x in resd.items(): G[r] += x
        else:
            for w, k in holes.items(): G[w] += k
        G = {r: x for r, x in G.items() if x}
        assert sum(r*x for r, x in G.items()) == M*W-11840 and min(G.values()) > 0
        b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
        res[mode] = dict(b=dec(b, 13), b_gain_pct=round(float((b/b0-1)*100), 3), kappa=dec(k, 16),
                         kappa_gain_pct=round(float((k/k0-1)*100), 3), binding='bit' if b > C_PR589 else 'complex')
    return res
N = 56
HR = {int(k): x for k, x in cen['helper_runs_total'].items()}
typ = {r: round(F(333*N*x, 504)) for r, x in HR.items()}; freed = sum(r*n for r, n in typ.items())
lo = freed // (333*N); hi = freed - lo*333*N
out = dict(control=dict(b=dec(b0, 25), kappa=dec(k0, 19)),
           typical56=go(typ, {lo: 333*N - hi, lo+1: hi}),
           mcopy56=go({int(r): 333*x//2 for r, x in cen['nonplane_runs'].items()}, {8: 333*28, 5: 333*16, 4: 333*12}))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(here, 'price_market.json'), 'w'), indent=1)
