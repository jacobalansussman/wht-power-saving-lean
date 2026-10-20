#!/usr/bin/env python3
"""k8-pairs: WHAT-IF of the pair market on the 8 non-K blocks (N = 32 or 36 fewer helpers per call) STACKED with
d-hpe-multigate's twin pairs g13 (574/575, PR #1090's pair) + g14 (580/581), which sit in K-blocks 8 and 7 (disjoint
from the non-K blocks).  Twin local {1:2400,2:446,3:416,4:142,5:82,8:8} (width 3726, from v-hpe-multigate-rules).
hole = #1042 map kept (removed helpers -> one idle child of their width; twins' 4 freed units per call -> rank-1 holes);
bins48 = all freed width re-packed into rank-48 children (optimistic).  Removed runs: M-copy class (N/112 of
{1:384,2:64,3:32,4:24}) and typical (N/480 of the k8 helper runs).
Usage: python3 -I price_k8_twins.py -> price_k8_twins.json"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as F
from collections import Counter
M, W, D, CALLS = 73, 14336, 5120, 144
K8 = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}; TW = {1: 2400, 2: 446, 3: 416, 4: 142, 5: 82, 8: 8}
IDLE = {1: 2249, 2: 1, 3: 400, 7: 3, 16: 47, 48: 4850}
P1042 = Counter(IDLE)
for r, x in K8.items(): P1042[r] += CALLS*x
b0 = price_b(dict(P1042), M, W); k0 = kappa_from(b0, C_PR589)[0]; assert k0 == F(1672188119205421, 10**18)
HR = {1: 1675, 2: 272, 3: 365, 4: 54, 5: 40}; NP = {1: 384, 2: 64, 3: 32, 4: 24}; NPW = {8: 56, 5: 32, 4: 24}
out = dict(control=dict(b=dec(b0, 25), kappa=dec(k0, 19)))
for N in (0, 32, 36):
    for kind in ('mcopy', 'typical'):
        if N == 0 and kind == 'typical': continue
        if kind == 'mcopy':
            runs = {r: round(F(CALLS*N*x, 112)) for r, x in NP.items()}
            holes = {w: round(F(CALLS*N*x, 112)) for w, x in NPW.items()}
            diff = sum(r*n for r, n in runs.items()) - sum(w*x for w, x in holes.items())
            if diff > 0: holes[diff] = holes.get(diff, 0) + 1
            elif diff < 0: holes[5] -= 1; holes[5 + diff] = holes.get(5 + diff, 0) + 1
        else:
            runs = {r: round(F(CALLS*N*x, 480)) for r, x in HR.items()}; fr = sum(r*n for r, n in runs.items())
            lo = fr // (CALLS*N); hi = fr - lo*CALLS*N; holes = {lo: CALLS*N - hi, lo + 1: hi}
        for mode in ('hole', 'bins48'):
            G = Counter(IDLE)
            for r, x in TW.items(): G[r] += CALLS*x
            for r, x in runs.items(): G[r] -= x
            freed = sum(r*x for r, x in runs.items()) + CALLS*(3730 - 3726)
            if mode == 'hole':
                for w, k in holes.items(): G[w] += k
                G[1] += CALLS*(3730 - 3726)
            else:
                q, rem = divmod(freed, 48); G[48] += q
                if rem: G[rem] += 1
            G = {r: x for r, x in G.items() if x}
            assert sum(r*x for r, x in G.items()) == M*W - D and min(G.values()) > 0
            b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
            out[f'twins+N{N}/{kind}/{mode}'] = dict(b=dec(b, 13), b_gain_pct=round(float((b/b0-1)*100), 3), kappa=dec(k, 16),
                                                     kappa_gain=float(k - k0), binding='bit' if b > C_PR589 else 'complex')
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(here, 'price_k8_twins.json'), 'w'), indent=1)
