#!/usr/bin/env python3
"""vtables.py (v-fp-capacity): recompute every number in fp-capacity's tables (NOTE sections 4-5, PROOF Cor C1,
Prop C) from the definitions: saving s with M = W m^(1-s); floor L_m (moves per array) from meta/MAP.md
(2, 14/5, 41/13, 2.99172, 2.97290, 2.9447 for m >= 4)."""
from math import log
from fractions import Fraction as Fr
floors = {2: 2.0, 3: 14 / 5, 4: 41 / 13, 5: 2.99172, 6: 2.97290, 8: 2.9447, 16: 2.9447, 73: 2.9447}
print("block floor  s_max   kappa=m/L-1  0.1:per-move 0.1:per-array  0.01:per-move 0.01:per-array")
for m, L in floors.items():
    smax = 1 - log(L) / log(m); kap = m / L - 1
    r = []
    for s in (0.1, 0.01):
        need = m ** s - 1; ex = "excl" if need > kap + 1e-12 else ""
        r.append("%.4f %s / %.3f" % (need, ex, m - m ** (1 - s)))
    print("%3d  %.5f  %.4f  %.4f   %s   %s" % (m, L, smax, kap, r[0], r[1]))
print("kappa exact: block3 =", Fr(3) / Fr(14, 5) - 1, " block4 =", Fr(4) / Fr(41, 13) - 1)
print("width cap k=2: ln2/ln m:", {m: round(log(2) / log(m), 4) for m in (2, 3, 4, 5, 6)})
for (m, W, M) in ((4, 2, 7), (5, 2, 6), (5, 2, 7), (5, 2, 8), (6, 2, 7), (6, 2, 8), (6, 2, 9), (6, 2, 10), (3, 5, 14), (4, 3, 10), (4, 4, 13)):
    s = 1 - log(M / W) / log(m); print("target %d in %d at block %d: Delta=%d s=%.4f  floor ok: %s" % (W, M, m, W * m - M, s, M >= W * floors.get(m, 2.9447) - 1e-9))
