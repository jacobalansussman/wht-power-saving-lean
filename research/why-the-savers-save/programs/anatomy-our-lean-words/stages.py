#!/usr/bin/env python3
"""an-ours: the stage-count knob of the P112 nine-stage family (own first-order arithmetic).
Geometry (pr717 P4; our identity.py): m = 1 + 8(2P+1), per replica 2Pv data arrays, 2P+1 calls,
D/replica = 2Pv - (2P+1) cst, each data array pays m - 1 levels, 8 per active stage in-call (one pair active per stage),
so idle cover per replica = 2v * 8(2P+1)(P-1); helpers' leftover width per call Rw must be hosted in it:
  (2P+1) Rw <= 2v * 8(2P+1)(P-1)   <=>   Rw <= 16 v (P-1).
"""
import math, json
v, cst = 112, 64
rows = []
for P in range(1, 9):
    m = 1 + 8 * (2 * P + 1)
    delta = 1 - (2 * P + 1) * cst / (2 * P * v)
    cover_per_call = 16 * v * (P - 1)
    rows.append(dict(P=P, stages=2 * P + 1, m=m, delta=delta, delta_over_m=delta / m,
                     rank1_log=math.log(m), cover_per_call=cover_per_call,
                     fits_U8=cover_per_call >= 3778, fits_atm3=cover_per_call >= 4060, cover_over_U8=cover_per_call / 3778))
    r = rows[-1]
    print(f"P={P} stages={2*P+1} m={m:3d} delta={delta:.4f} delta/m={delta/m:.5f} ln m={math.log(m):.3f} "
          f"idle cover per call={cover_per_call:5d} (U8 leftovers 3778: ratio {cover_per_call/3778:.2f})")
print('limit P->inf: delta -> 1 - cst/v =', 1 - cst / v)
json.dump(rows, open('stages.json', 'w'), indent=1)
