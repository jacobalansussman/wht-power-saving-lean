#!/usr/bin/env python3
"""u-helpers-hub: price removing the 112 U8 helpers that never visit a packet plane (all entering helpers; widths
56x8, 32x5, 24x4; actual runs {1:384,2:64,3:32,4:24} per call from u8_census.json) on #942, in all 333 calls.
Hosting models as in price_removals.py: optimistic (re-pack relaxation, bins filled smallest capacity first) and
hole (each removed helper leaves a hole of its own width, paid as one child). WHAT-IF prices, not a construction.
Usage: python3 -I price_classes.py   -> price_classes.json"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('PRICER_DIR', '.'))
from pricer import price_b, kappa_from, C_PR589, dec
from collections import Counter
cen = json.load(open(os.path.join(here, 'u8_census.json')))
P942 = {1: 848490, 2: 112554, 3: 148518, 4: 53946, 5: 23976, 8: 2664, 48: 11086}; M, W = 73, 33152
b0 = price_b(P942, M, W)
def go(runs, holes_by_width):
    res = {}
    for mode in ('optimistic', 'hole'):
        G = Counter(P942)
        for r, x in runs.items(): G[int(r)] -= 333*x
        freed = 333*sum(int(r)*x for r, x in runs.items())
        if mode == 'optimistic':
            G[48] -= 11086; G[1] -= 6; hosted = 37*9*3778 - freed; resd = Counter()
            for cap, nb in ((8, 672*37), (16, 448*37), (32, 448*37), (48, 448*37)):
                full = min(nb, hosted//cap); hosted -= full*cap; left = nb-full
                if left and hosted: resd[cap-hosted] += 1; left -= 1; hosted = 0
                resd[cap] += left
            for r, x in resd.items(): G[r] += x
        else:
            for w, k in holes_by_width.items(): G[w] += 333*k
        G = {r: x for r, x in G.items() if x}
        assert sum(r*x for r, x in G.items()) == M*W-11840 and min(G.values()) > 0
        b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
        res[mode] = dict(b=dec(b, 13), pct=round(float((b/b0-1)*100), 3), kappa=dec(k, 16))
    return res
out = dict(nonplane112=go(cen['nonplane_runs'], {8: 56, 5: 32, 4: 24}))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(here, 'price_classes.json'), 'w'), indent=1)
