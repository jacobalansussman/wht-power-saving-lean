#!/usr/bin/env python3
"""v-slack: arithmetic on the known units (from the release README section 6 table and an-ours identity.out), parity of
the move deficit, and recomputation of every ceiling / threshold number quoted by the two claims."""
from math import log
import json

units = {
    # name: (h, v, m, R, cst)  -- release README sec. 6-7: W = 4v + R, D = 4v - 5 cst
    'fourth result (E8)': (9, 120, 45, 783, 72),
    'second/third result (m=110)': (22, 1320, 110, 9412, 440),
    'first result (m=80)': (16, 560, 80, 7122, 241),
}
readme = {  # the README's own table: W, D, moves of W m
    'fourth result (E8)': (1263, 120, 56715, 56835),
    'second/third result (m=110)': (14692, 3080, 1613040, 1616120),
    'first result (m=80)': (9362, 1035, 747925, 748960),
}
out = {'units': {}}
for k, (h, v, m, R, cst) in units.items():
    W = 4 * v + R; D = 4 * v - 5 * cst; Wm = W * m; M = Wm - D
    rw, rd, rM, rWm = readme[k]
    smc = 1 - log(M / W) / log(m)
    out['units'][k] = {'m': m, 'W': W, 'D': D, 'Wm': Wm, 'M': M, 'matches_README_table': (W, D, M, Wm) == (rw, rd, rM, rWm),
                       'D_parity': 'odd' if D % 2 else 'even', 'share D/(Wm)': D / Wm,
                       'move-count exponent': smc,
                       'ThmR ceiling ln(m/2.9447)/ln m': log(m / 2.9447) / log(m),
                       'one-overlap common-frame move-count ceiling -ln(1-1/m)/ln m': -log(1 - 1 / m) / log(m)}
out['units']['U8 per replica (an-ours identity.out)'] = {'data arrays': 896, 'D per replica': 320, 'D_parity': 'even'}
out['units']['stand-in 2 in 3 (not engine)'] = {'W': 2, 'm': 2, 'M': 3, 'D': 1, 'D_parity': 'odd'}
# Theorem R ceilings
def s_of(per_array, m): return 1 - log(per_array) / log(m)
out['ThmR'] = {'block3 14/5': s_of(14 / 5, 3), 'block4 41/13': s_of(41 / 13, 4), 'block5 2.99172': s_of(2.99172, 5),
               'block6 2.97290': s_of(2.97290, 6), 'block8 2.9447': s_of(2.9447, 8)}
# thresholds: overlap k per array, free bridges: s <= -ln(1 - k/m)/ln m ; largest m with s >= target
def maxm(k, target):
    best = None
    for m in range(k + 1, 5000):
        if -log(1 - k / m) / log(m) >= target: best = m
    return best
out['thresholds'] = {f'overlap {k}, s={t}': maxm(k, t) for k in (1, 2, 4) for t in (0.1, 0.01)}
# needs tables
ms = [3, 4, 5, 6, 8, 16, 25, 41, 73]
out['needs_deficit_per_array'] = {str(s): {m: round(m - m ** (1 - s), 3) for m in ms} for s in (0.1, 0.01)}
out['needs_rawpacking'] = {str(s): {m: round(m ** s - 1, 4) for m in (4, 8, 16, 73)} for s in (0.1, 0.01)}
out['needs_coherent_avg_level'] = {str(s): {m: round(m / 2 * (1 - m ** (-s)), 4) for m in (4, 8, 16, 73)} for s in (0.1, 0.01)}
out['two-array targets'] = {f'2 in {M} at block {m}': 1 - log(M / 2) / log(m)
                            for (M, m) in [(7, 4), (7, 5), (8, 5), (6, 5), (8, 6), (9, 6), (10, 6), (7, 6)]}
out['block3: 0.1 needs M/W <='] = 3 ** 0.9
out['block4: 0.1 needs M/W <='] = 4 ** 0.9
out['169/57'] = 169 / 57
print(json.dumps(out, indent=1))
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
