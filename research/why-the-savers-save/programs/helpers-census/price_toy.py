#!/usr/bin/env python3
"""u-helpers-hub: price each variant of the packet-plane visit toy (toy_visits.json) on #942.
A variant's least helper count H replaces U8's 392 plane-visiting helpers; the 392-H removed helpers are priced as
mean plane helpers (runs = U8 plane-helper runs/392 per call, from u8_census.json) in all 333 calls, with the two
hosting models of price_removals.py. Gross WHAT-IF prices: they do not charge whatever the variant needs instead
(descents, output conjugation, a different packet decomposition).  Usage: python3 -I price_toy.py"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('PRICER_DIR', '.'))
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as F
from collections import Counter
cen = json.load(open(os.path.join(here, 'u8_census.json'))); toy = json.load(open(os.path.join(here, 'toy_visits.json')))
P942 = {1: 848490, 2: 112554, 3: 148518, 4: 53946, 5: 23976, 8: 2664, 48: 11086}; M, W = 73, 33152
b0 = price_b(P942, M, W)
PR = {int(r): x for r, x in cen['plane_runs'].items()}
def price(N):
    d = {r: round(F(333*N*x, 392)) for r, x in PR.items()}; freed = sum(r*x for r, x in d.items())
    res = {}
    for mode in ('optimistic', 'hole'):
        G = Counter(P942)
        for r, x in d.items(): G[r] -= x
        if mode == 'optimistic':
            G[48] -= 11086; G[1] -= 6; hosted = 37*9*3778 - freed; resd = Counter()
            for cap, nb in ((8, 672*37), (16, 448*37), (32, 448*37), (48, 448*37)):
                full = min(nb, hosted//cap); hosted -= full*cap; left = nb-full
                if left and hosted: resd[cap-hosted] += 1; left -= 1; hosted = 0
                resd[cap] += left
            for r, x in resd.items(): G[r] += x
        else:
            holes = 333*N; lo = freed//holes; hi = freed - lo*holes; G[lo] += holes-hi; G[lo+1] += hi
        G = {r: x for r, x in G.items() if x}
        assert sum(r*x for r, x in G.items()) == M*W-11840 and min(G.values()) > 0
        b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
        res[mode] = dict(b=dec(b, 13), b_pct=round(float((b/b0-1)*100), 3), kappa=dec(k, 16))
    return res
out = {}
for name, cnt in toy['least_helpers'].items():
    H = cnt['all14']; N = 392 - H
    out[name] = dict(least_plane_helpers=H, removed_vs_U8=N, **(price(N) if N > 0 else {}))
    print(name, json.dumps(out[name]))
out['hub_stronger_cap_matching_28'] = dict(least_plane_helpers=364, removed_vs_U8=28, **price(28))
print('hub 28-helper output-borrowing matching', json.dumps(out['hub_stronger_cap_matching_28']))
json.dump(out, open(os.path.join(here, 'price_toy.json'), 'w'), indent=1)
