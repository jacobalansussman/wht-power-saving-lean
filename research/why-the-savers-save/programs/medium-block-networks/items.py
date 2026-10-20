#!/usr/bin/env python3
"""cx-medium / items.py (own code, stdlib; python3 -I).  Every cost of every template point in blocks 13..25, itemised
per replica (one replica = the 2Pv data arrays of one run), and the narrow-mass decomposition of the one configuration
with D > 0 (P112 at block 25, U8X0 word, helpers on banks with ideal packing).
Inputs: natural.json (coordinate decoders), the U8X0 per-call class histogram (csx-rerun/csx_checks.json, our exact
checker run), price_p112.json.  Template: P bank pairs, 2P+1 stages, window k, sigma = k-1, m = 1 + sigma(2P+1);
free directions = 2Pv (one per data array); centre = (2P+1) c."""
import json, os, math
HERE = os.path.dirname(os.path.abspath(__file__))
nat = json.load(open(os.path.join(HERE, 'natural.json')))
rows = []
def row(design, P, c, label):
    n = nat[design]['n']; k = nat[design]['k']; s = k - 1
    m = 1 + s * (2 * P + 1); free = 2 * P * n; cent = (2 * P + 1) * c
    rows.append(dict(design=design, P=P, m=m, data_arrays=2 * P * n, free=free, centre=cent, D_before_junk=free - cent, c=c, decoder=label,
                     break_even_c=free / (2 * P + 1)))
for d, P in (('E7', 1), ('D6', 1), ('A5', 1), ('D6', 2), ('A5', 2), ('P112', 1), ('E8', 1)):
    row(d, P, nat[d]['hub_cost'], 'coordinate decoder (natural.py)')
row('P112', 1, 64, "hub's decoder (U8X0, Lean-checked)")
row('E7', 1, 42, 'hypothetical t(k-1) decoder (not found)')
row('D6', 2, 24, 'hypothetical t(k-1) decoder (not found)')
print(f"{'design':6s} {'P':>2s} {'m':>3s} {'data':>5s} {'free':>5s} {'centre':>6s} {'D':>5s}  c/call (break-even)  decoder")
for r in rows:
    print(f"{r['design']:6s} {r['P']:2d} {r['m']:3d} {r['data_arrays']:5d} {r['free']:5d} {r['centre']:6d} {r['D_before_junk']:5d}  {r['c']:3d} ({r['break_even_c']:.1f})  {r['decoder']}")
# P112 at 25: narrow mass by class (U8X0 classes from our exact checker), per replica (3 calls), ideal banks
cx = json.load(open(os.path.join(HERE, 'csx-rerun', 'csx_checks.json')))['U8X0']['hist_by_class']
m = 25
N = {}; R = {}
for cl, h in cx.items():
    R[cl] = 3 * sum(int(r) * c for r, c in h.items())
    N[cl] = 3 * sum(int(r) * c * math.log(m / int(r)) for r, c in h.items())
W = 224 + 3 * 3778 / 25
D = 32
tot = sum(N.values())
print('\nP112, block 25, one replica (3 calls), U8X0 word, helpers on banks (ideal packing):')
print(f"  W = 224 data + {3*3778/25:.2f} banks = {W:.2f};  paid levels: x {R['x']}, y {R['y']}, helpers {R['s']}, copies {R['c']}; D = {D}")
for cl, lab in (('x', 'sources'), ('y', 'targets'), ('s', 'helpers (on banks)'), ('c', 'centre copies')):
    print(f"  narrow mass {lab:20s} {N[cl]:9.1f}  ({100*N[cl]/tot:5.1f} %)")
print(f"  first order a = D/N = {D/tot:.4e} (the fold gives 7.3057223e-4)")
print(f"  with helpers free (no banks, no helper mass): D/N = {D/(tot-N['s']):.4e}")
json.dump(dict(rows=rows, p112_25=dict(W=W, D=D, levels=R, narrow=N, first_order=D / tot, helpers_free=D / (tot - N['s']))),
          open(os.path.join(HERE, 'items.json'), 'w'), indent=1)

# ---- what P112 at block 25 would need to beat the hub, and the most generous ceiling
HUB = 1.6797550e-3
budget = D / HUB - (N['x'] + N['y'] + N['c'])
print(f"\nto beat the hub at block 25 (first order, D = 32 fixed, U8X0's data and copy blocks kept):"
      f" helper narrow mass per replica <= {budget:.0f}, i.e. <= {100*budget/N['s']:.1f} % of U8X0's {N['s']:.0f}")
# most generous: helpers free, every data array's 24 levels in the fewest children the stage gates allow.
# A data array has a gate in each of its 3 calls; its frames in stage s lie in [8s(+1), 8s+8(+1)], so one cut lies in
# the middle stage's range: the best split of 24 levels is {8, 16}.
best = min(sum(r * math.log(m / r) for r in parts) for parts in ([8, 16], [16, 8], [12, 12], [24]) if parts != [24])
ceil = D / (224 * best + N['c'])
print(f"most generous ceiling at block 25 (helpers free, data in children {{8, 16}}, hub's centre): a <= {ceil:.3e}")
for name, Dd, nd, cst in (('E8 (120 ports, coordinate decoder 72)', 24, 240, 216),):
    print(f"  {name}: D = {Dd}, ceiling a <= {Dd / (nd * best + cst * math.log(25 / 9)):.3e}")
json.dump(dict(rows=rows, p112_25=dict(W=W, D=D, levels=R, narrow=N, first_order=D / tot, helpers_free=D / (tot - N['s']),
               beat_hub_helper_budget=budget, beat_hub_fraction=budget / N['s'], generous_ceiling=ceil)),
          open(os.path.join(HERE, 'items.json'), 'w'), indent=1)
