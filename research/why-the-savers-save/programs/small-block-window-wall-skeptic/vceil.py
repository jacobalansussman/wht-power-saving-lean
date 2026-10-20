#!/usr/bin/env python3
"""v-an2-smallnet / vceil.py (own code, stdlib).  Re-solves the moment equation sum_children c (r/m)^(1-a) = W
(W = 2 per pair) by Newton/bisection, for an2-smallnet's model: data 2(m-u) levels and centre 2c levels per pair, all
in children of rank sigma.  Also a sensitivity line: the same pair with each data array's levels in ONE child
(rank m-u), to show how much the ceilings depend on the sigma-block grouping assumption."""
import json, math, sys
def solve(ch, W=2.0):
    lo, hi = -5.0, 0.999999
    for _ in range(300):
        a = (lo + hi) / 2
        if sum(c * x ** (1 - a) for x, c in ch) > W: hi = a
        else: lo = a
    return (lo + hi) / 2
out = {}
out['calibration m16 u2 r4 free bridge'] = solve([(4/16, 2*(16-2)/4)])
out['closed form 1 - ln3.5/ln4'] = 1 - math.log(3.5)/math.log(4)
for name, sg, c in (('E7', 6, 42/63), ('P112', 8, 64/112), ('E8', 8, 64/120)):
    for m in (1 + 2*sg, 1 + 3*sg):
        for u in (1, 1.5, 2):
            out[f'{name} m={m} u={u} sigma-blocks'] = solve([(sg/m, 2*(m-u)/sg), (sg/m, 2*c/sg)])
        out[f'{name} m={m} u=1 ONE data child of rank m-1 (sensitivity)'] = solve([((m-1)/m, 2.0), (sg/m, 2*c/sg)])
for k, v in out.items(): print(f'{k:60s} {v:.4f}')
json.dump(out, open(sys.argv[1], 'w'), indent=1)
