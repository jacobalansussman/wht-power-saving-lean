#!/usr/bin/env python3
"""cx-medium / ceilings.py (own code, stdlib; python3 -I).  The most generous ceilings of the one-pair template:
junk cleaning FREE (no helpers, no banks), every data array's 3 sigma levels in the two largest pieces its gates allow
(one cut must fall in the middle call's range, so {sigma, 2 sigma}), the centre's copies as pieces of their rank.
P112 / E8 at block 25 with the fold (an-ours anatomy_ours.price); P112 with U8X0's own data pieces; E7 at block 19 as
a function of the centre charge c per call (uncapped moment root: the fold's bisection stops at 0.01)."""
import sys, os, json
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
from anatomy_ours import price, rank
def root(H, W, m):
    f = lambda a: sum(n * (r / m) ** (1 - a) for r, n in H.items()) - W
    lo, hi = 0.0, 0.9
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if f(mid) > 0 else (mid, hi)
    return lo
out = {}
H = {8: 224 + 24, 16: 224}; out['P112_25_free_junk_best_pieces'] = float(price(H, 224, m=25)); assert 25 * 224 - rank(H) == 32
x = {1: 442, 2: 12, 3: 30, 4: 50, 5: 28}; y = {1: 352, 2: 56, 3: 32, 4: 64, 5: 16}
H2 = {}
for h in (x, y):
    for r, n in h.items(): H2[r] = H2.get(r, 0) + 3 * n
H2[8] = H2.get(8, 0) + 24
out['P112_25_free_junk_U8X0_pieces'] = float(price(H2, 224, m=25))
H3 = {8: 240, 16: 240, 9: 24}; out['E8_120_25_free_junk_best_pieces'] = float(price(H3, 240, m=25))
for c in (42, 41, 40, 39, 38, 37, 36, 35, 30):
    T = 6; Hc = {6: 126 * T + 3 * c * T / 6, 12: 126 * T}
    out[f'E7_19_free_junk_best_pieces_c{c}'] = root(Hc, 126 * T, 19) if c < 42 else 0.0
for k, v in out.items(): print(f"{k:42s} {v:.4e}")
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ceilings.json'), 'w'), indent=1)

# ---- the 112-port design at block 25 as a function of a hypothetical centre charge c per call (8 copies of rank c/8),
#      junk cleaning free, data pieces {8, 16}: where it crosses the hub and 1e-2; the clique floor is 52.36 (=> 53)
for c in (64, 63, 62, 61, 60, 58, 56, 53):
    H = {8: 224, 16: 224}
    H[c / 8] = H.get(c / 8, 0) + 24          # (c = 64: rank-8 copies join the rank-8 data pieces)
    out[f'P112_25_free_junk_best_pieces_c{c}'] = root(H, 224, 25)
    print(f"P112 block 25, centre {c}/call, junk free, pieces {{8,16}}: D = {224 - 3 * c}, root {out[f'P112_25_free_junk_best_pieces_c{c}']:.4e}")
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ceilings.json'), 'w'), indent=1)
