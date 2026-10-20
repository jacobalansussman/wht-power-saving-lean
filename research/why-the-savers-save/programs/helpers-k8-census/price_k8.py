#!/usr/bin/env python3
"""k8-pairs: WHAT-IF prices of N fewer helpers per call on the #1042 control (k8 word, T16 owner map, 144 calls).
#1042 = 144 x local {1:2412,2:418,3:432,4:142,5:82,8:8} + idle {1:2249,2:1,3:400,7:3,16:47,48:4850}, m 73, W 14336,
deficit 5120 (hub research/community-frontier1042-audit/INDEPENDENT-T16.json).  Control reproduced first (exact).
Removed helpers: (a) 'typical' = N x the mean k8 helper runs (census: {1:1675,2:272,3:365,4:54,5:40} over 480);
(b) 'mcopy' = N/112 of the 112 non-plane helpers' runs {1:384,2:64,3:32,4:24} (widths {8:56,5:32,4:24}), the class
d-pair-market's own-block saving is set against.  Hosting: 'hole' = #1042's owner map kept, each removed helper
leaves one idle child of its width; 'bins48' = the freed width re-packed into full rank-48 children (remainder one
child) - an optimistic re-pack, not a constructed map.  Same method as d-pair-market/price_market.py.
Usage: python3 -I price_k8.py -> price_k8.json"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as F
from collections import Counter
P1042 = {1: 349577, 2: 60193, 3: 62608, 4: 20448, 5: 11808, 7: 3, 8: 1152, 16: 47, 48: 4850}; M, W, D = 73, 14336, 5120
CALLS = 144
b0 = price_b(P1042, M, W); k0 = kappa_from(b0, C_PR589)[0]
assert b0 == F(418747253984418082117, 250000000000000000000000) and k0 == F(1672188119205421, 10**18)
HR = {1: 1675, 2: 272, 3: 365, 4: 54, 5: 40}; assert sum(r*x for r, x in HR.items()) == 3730
NP = {1: 384, 2: 64, 3: 32, 4: 24}; NPW = {8: 56, 5: 32, 4: 24}
def go(runs_global, holes):
    res = {}
    freed = sum(r*x for r, x in runs_global.items())
    for mode in ('hole', 'bins48'):
        G = Counter(P1042)
        for r, x in runs_global.items(): G[r] -= x
        if mode == 'hole':
            for w, k in holes.items(): G[w] += k
        else:
            q, rem = divmod(freed, 48); G[48] += q
            if rem: G[rem] += 1
        G = {r: x for r, x in G.items() if x}
        assert sum(r*x for r, x in G.items()) == M*W - D and min(G.values()) > 0
        b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
        res[mode] = dict(b=dec(b, 13), b_gain_pct=round(float((b/b0-1)*100), 3), kappa=dec(k, 16),
                         kappa_gain=float(k - k0), kappa_gain_pct=round(float((k/k0-1)*100), 3),
                         binding='bit' if b > C_PR589 else 'complex')
    return res
out = dict(control=dict(b=dec(b0, 25), kappa=dec(k0, 19)))
for N in (32, 36, 40, 56):
    typ = {r: round(F(CALLS*N*x, 480)) for r, x in HR.items()}; fr = sum(r*n for r, n in typ.items())
    lo = fr // (CALLS*N); hi = fr - lo*CALLS*N
    mc = {r: round(F(CALLS*N*x, 112)) for r, x in NP.items()}; frm = sum(r*n for r, n in mc.items())
    # holes of the M-copy class: N/112 of the widths {8:56,5:32,4:24}, made to carry exactly the freed width
    hw = {w: round(F(CALLS*N*x, 112)) for w, x in NPW.items()}
    diff = frm - sum(w*x for w, x in hw.items())
    if diff > 0: hw[diff] = hw.get(diff, 0) + 1              # rounding remainder as one extra small hole
    elif diff < 0: hw[5] -= 1; hw[5 + diff] = hw.get(5 + diff, 0) + 1
    out[f'N{N}'] = dict(typical=go(typ, {lo: CALLS*N - hi, lo+1: hi}), mcopy=go(mc, {w: x for w, x in hw.items() if x}))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(here, 'price_k8.json'), 'w'), indent=1)
